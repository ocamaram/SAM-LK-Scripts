#!/usr/bin/env python3
"""
shade_postprocess.py
Post-procesado de resultados del 3D Shade Calculator de SAM.

Salidas:
  summary_statistics.csv        Resumen estadístico por grupo (subarray+string)
  seasonal_hourly_curves.png    Curvas de sombra horaria promedio por estación
  heatmap_panel_geometry.png    Heatmap estacional y anual sobre geometría de paneles
"""

import argparse
import os
import sys
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle

# ── ESTACIONES (hemisferio norte) ───────────────────────────────────────────
SEASONS = {
    'Invierno': [12, 1, 2],
    'Primavera': [3, 4, 5],
    'Verano':    [6, 7, 8],
    'Otoño':     [9, 10, 11],
}
SEASON_COLORS = {
    'Invierno': '#4e9af1',
    'Primavera': '#4caf50',
    'Verano':    '#f4a261',
    'Otoño':     '#e07b39',
}


def parse_group(name):
    """
    Parsea el identificador de grupo a (subarray_1based, string_1based).
    Formatos soportados:
      '3.5'               → SAM interno 0-based  → (4, 6)
      'SA2_ST5'           → nuestro formato       → (2, 5)
      'Subarray 4, String 6' → nombre display SAM → (4, 6)
    Devuelve None si no encaja ningún patrón.
    """
    name = str(name).strip()
    # '3.5' o '0.0' → subarray.string 0-based
    m = re.match(r'^(\d+)\.(\d+)$', name)
    if m:
        return (int(m.group(1)) + 1, int(m.group(2)) + 1)
    # 'SA2_ST5'
    m = re.match(r'SA(\d+)_ST(\d+)', name)
    if m:
        return (int(m.group(1)), int(m.group(2)))
    # 'Subarray 4, String 6'
    m = re.match(r'Subarray\s+(\d+),\s+String\s+(\d+)', name)
    if m:
        return (int(m.group(1)), int(m.group(2)))
    return None


def group_label(name):
    """Convierte cualquier formato de grupo a 'SA<N>_ST<M>' legible."""
    parsed = parse_group(name)
    if parsed:
        return f'SA{parsed[0]}_ST{parsed[1]}'
    return name


def build_time_index(n):
    return pd.date_range('2001-01-01', periods=n, freq='h')


# ── ESTADÍSTICAS ────────────────────────────────────────────────────────────

def compute_statistics(ts, diff):
    stats = pd.DataFrame({
        'sombra_media_%':   (ts.mean()          * 100).round(2),
        'sombra_mediana_%': (ts.median()         * 100).round(2),
        'sombra_p90_%':     (ts.quantile(0.90)   * 100).round(2),
        'sombra_max_%':     (ts.max()             * 100).round(2),
        'h_con_sombra':     (ts > 0).sum(),
        'h_sombra_>10%':    (ts > 0.10).sum(),
        'h_sombra_>50%':    (ts > 0.50).sum(),
    })
    # Añadir sombra difusa si los grupos coinciden
    if diff is not None:
        col = diff.columns[0]
        stats['sombra_difusa_%'] = stats.index.map(
            lambda g: diff.loc[g, col] if g in diff.index else np.nan
        )
    return stats


# ── CURVAS HORARIAS POR ESTACIÓN ─────────────────────────────────────────────

def seasonal_hourly_avg(ts):
    """Dict {season: DataFrame(24h × groups)} con media de sombra."""
    result = {}
    for season, months in SEASONS.items():
        mask = ts.index.month.isin(months)
        result[season] = (
            ts[mask]
            .groupby(ts[mask].index.hour)
            .mean()
        )
    return result


def plot_seasonal_curves(seas_avg, groups, out_path):
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), sharey=True, sharex=True)
    axes = axes.flatten()

    n_groups = len(groups)
    cmap_grp = plt.get_cmap('tab20', n_groups)

    for ax, (season, df) in zip(axes, seas_avg.items()):
        # Líneas individuales por grupo
        for i, grp in enumerate(groups):
            ax.plot(df.index, df[grp] * 100,
                    color=cmap_grp(i), alpha=0.55, linewidth=1.1)

        # Media global del array (línea gruesa)
        array_mean = df[groups].mean(axis=1)
        ax.plot(array_mean.index, array_mean * 100,
                color='black', linewidth=2.2, linestyle='--', label='Media array')

        ax.set_title(season, fontsize=12, fontweight='bold',
                     color=SEASON_COLORS[season])
        ax.set_xlim(0, 23)
        ax.set_ylim(0, 100)
        ax.set_xticks(range(0, 24, 3))
        ax.set_xticklabels([f'{h:02d}:00' for h in range(0, 24, 3)])
        ax.set_xlabel('Hora del día')
        ax.set_ylabel('Sombra (%)')
        ax.grid(True, alpha=0.25)
        ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%g%%'))

    # Leyenda de grupos (solo si hay pocos)
    if n_groups <= 20:
        handles = [
            plt.Line2D([0], [0], color=cmap_grp(i), linewidth=2, label=g)
            for i, g in enumerate(groups)
        ]
        handles.append(
            plt.Line2D([0], [0], color='black', linewidth=2.2,
                       linestyle='--', label='Media array')
        )
        fig.legend(handles=handles, loc='lower center',
                   ncol=min(8, n_groups + 1), fontsize=8,
                   bbox_to_anchor=(0.5, -0.01))

    fig.suptitle('Sombra directa horaria promedio por estación',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Curvas estacionales → {out_path}')


# ── HEATMAP ──────────────────────────────────────────────────────────────────

def build_grids(seas_avg, ts, groups):
    """Dict {label: ndarray(n_sa × n_st)} con sombra media (%)."""
    coords = {g: parse_group(g) for g in groups}
    valid  = {g: c for g, c in coords.items() if c is not None}

    if not valid:
        print('AVISO: ningún grupo sigue el patrón SA_ST — heatmap omitido.')
        return None, 0, 0

    n_sa = max(c[0] for c in valid.values())
    n_st = max(c[1] for c in valid.values())

    grids = {}
    for season, df in seas_avg.items():
        grid = np.full((n_sa, n_st), np.nan)
        for g, (sa, st) in valid.items():
            grid[sa - 1, st - 1] = df[g].mean() * 100
        grids[season] = grid

    # Anual
    grid_ann = np.full((n_sa, n_st), np.nan)
    for g, (sa, st) in valid.items():
        grid_ann[sa - 1, st - 1] = ts[g].mean() * 100
    grids['Anual'] = grid_ann

    return grids, n_sa, n_st


def plot_heatmaps(grids, n_sa, n_st, panel_w, panel_l, out_path):
    labels  = list(grids.keys())   # 4 estaciones + Anual
    n_plots = len(labels)          # 5

    all_vals = np.concatenate([g[~np.isnan(g)] for g in grids.values()])
    vmin, vmax = np.percentile(all_vals, 2), np.percentile(all_vals, 98)

    cmap = LinearSegmentedColormap.from_list(
        'shade', ['#fffde7', '#ffcc02', '#e65100'], N=256
    )
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm   = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    # ── Dimensiones físicas exactas por subplot ───────────────
    cell_in  = 0.55           # pulgadas por metro de panel
    ax_w     = n_st * panel_w * cell_in   # ancho del axes en pulgadas
    ax_h     = n_sa * panel_l * cell_in   # alto del axes en pulgadas

    n_cols, n_rows = 3, 2
    cbar_w   = 0.25           # ancho colorbar en pulgadas
    col_gap  = 0.80           # espacio horizontal entre subplots
    row_gap  = 0.70           # espacio vertical entre subplots
    pad_l    = 0.55           # margen izquierdo
    pad_r    = 0.15           # margen derecho
    pad_top  = 0.55           # margen superior (título)
    pad_bot  = 0.40           # margen inferior

    slot_w   = ax_w + cbar_w + 0.15   # ancho de cada slot (axes + colorbar)
    fig_w    = pad_l + n_cols * slot_w + (n_cols - 1) * col_gap + pad_r
    fig_h    = pad_top + n_rows * ax_h + (n_rows - 1) * row_gap + pad_bot

    fig = plt.figure(figsize=(fig_w, fig_h))

    for idx, label in enumerate(labels):
        row = idx // n_cols
        col = idx % n_cols

        # Posición del axes en coordenadas normalizadas de figura
        ax_left   = (pad_l + col * (slot_w + col_gap)) / fig_w
        ax_bottom = (pad_bot + (n_rows - 1 - row) * (ax_h + row_gap)) / fig_h
        ax_rect   = [ax_left, ax_bottom, ax_w / fig_w, ax_h / fig_h]

        ax = fig.add_axes(ax_rect)
        grid = grids[label]

        for sa in range(n_sa):
            for st in range(n_st):
                val = grid[sa, st]
                if np.isnan(val):
                    continue
                rect = Rectangle(
                    xy=(st * panel_w, sa * panel_l),
                    width=panel_w, height=panel_l,
                    facecolor=cmap(norm(val)),
                    edgecolor='white', linewidth=0.5
                )
                ax.add_patch(rect)

        ax.set_xlim(0, n_st * panel_w)
        ax.set_ylim(n_sa * panel_l, 0)   # Y invertido: SA1 arriba
        # Sin aspect='equal': las proporciones correctas vienen del figsize

        season_color = SEASON_COLORS.get(label, '#333333')
        ax.set_title(label, fontsize=10, fontweight='bold', color=season_color, pad=4)
        ax.set_xlabel('String', fontsize=8, labelpad=3)
        ax.set_ylabel('Subarray', fontsize=8, labelpad=3)

        ax.set_xticks([(st + 0.5) * panel_w for st in range(n_st)])
        ax.set_xticklabels([f'ST{st+1}' for st in range(n_st)], fontsize=6)
        ax.set_yticks([(sa + 0.5) * panel_l for sa in range(n_sa)])
        ax.set_yticklabels([f'SA{sa+1}' for sa in range(n_sa)], fontsize=6)

        # Colorbar en axes propio, adyacente al subplot
        cbar_left   = ax_left + ax_w / fig_w + 0.01
        cbar_rect   = [cbar_left, ax_bottom, cbar_w / fig_w, ax_h / fig_h]
        cax = fig.add_axes(cbar_rect)
        cbar = fig.colorbar(sm, cax=cax)
        cbar.set_label('Sombra (%)', fontsize=7)
        cbar.ax.tick_params(labelsize=6)

    fig.text(0.5, 1 - pad_top / fig_h / 2, 'Heatmap de sombra — geometría de paneles',
             ha='center', va='top', fontsize=12, fontweight='bold')

    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Heatmap → {out_path}')


# ── MAIN ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Post-proceso sombras SAM')
    parser.add_argument('--ts',     required=True,       help='CSV series temporales')
    parser.add_argument('--diff',   required=True,       help='CSV sombra difusa')
    parser.add_argument('--out',    required=True,       help='Carpeta de salida')
    parser.add_argument('--width',  type=float, default=1.0, help='Ancho del panel (m)')
    parser.add_argument('--length', type=float, default=1.0, help='Largo del panel (m)')
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # ── Carga ────────────────────────────────────────────────────
    print('Cargando datos...')
    ts_raw = pd.read_csv(args.ts, index_col=0)
    ts_raw.index = build_time_index(len(ts_raw))

    # Renombrar columnas al formato legible SA<N>_ST<M>
    ts_raw.columns = [group_label(c) for c in ts_raw.columns]
    groups = list(ts_raw.columns)
    print(f'  {len(groups)} grupos: {groups[:4]}{"..." if len(groups) > 4 else ""}')

    # Los valores del time series están en % (0-100) → normalizar a fracción 0-1
    if ts_raw.max().max() > 1.5:
        ts_raw = ts_raw / 100.0
        print('  Valores time series convertidos de % a fracción (0-1).')

    try:
        diff_raw = pd.read_csv(args.diff, index_col=0)
        # Renombrar índice del diffuse al mismo formato SA<N>_ST<M>
        diff_raw.index = [group_label(i) for i in diff_raw.index]
        # Los valores del diffuse son fracción (0-1) → convertir a %
        col = diff_raw.columns[0]
        if diff_raw[col].max() <= 1.5:
            diff_raw[col] = diff_raw[col] * 100.0
            print('  Valores diffuse convertidos de fracción a %.')
    except Exception:
        diff_raw = None
        print('  AVISO: no se pudo cargar el CSV de sombra difusa.')

    # ── Estadísticas ─────────────────────────────────────────────
    print('Calculando estadísticas...')
    stats = compute_statistics(ts_raw, diff_raw)
    stats_path = os.path.join(args.out, 'summary_statistics.csv')
    stats.to_csv(stats_path)
    print(f'  Estadísticas → {stats_path}')

    # ── Curvas estacionales ───────────────────────────────────────
    print('Generando curvas horarias por estación...')
    seas_avg = seasonal_hourly_avg(ts_raw)
    plot_seasonal_curves(
        seas_avg, groups,
        os.path.join(args.out, 'seasonal_hourly_curves.png')
    )

    # ── Heatmaps ─────────────────────────────────────────────────
    print('Generando heatmaps...')
    grids, n_sa, n_st = build_grids(seas_avg, ts_raw, groups)
    if grids:
        plot_heatmaps(
            grids, n_sa, n_st,
            args.width, args.length,
            os.path.join(args.out, 'heatmap_panel_geometry.png')
        )

    print('\nCompletado.')


if __name__ == '__main__':
    main()
