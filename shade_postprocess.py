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
from matplotlib.colors import LinearSegmentedColormap

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


def plot_heatmaps(grids, n_sa, n_st, out_path):
    labels  = list(grids.keys())          # 4 estaciones + Anual
    n_plots = len(labels)                  # 5

    # Layout 2×3 (última celda vacía)
    n_cols, n_rows = 3, 2
    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(5 * n_cols, 4.5 * n_rows),
                              gridspec_kw={'wspace': 0.45, 'hspace': 0.55})
    axes_flat = axes.flatten()

    all_vals = np.concatenate([g[~np.isnan(g)] for g in grids.values()])
    vmin, vmax = np.percentile(all_vals, 2), np.percentile(all_vals, 98)

    cmap = LinearSegmentedColormap.from_list(
        'shade', ['#fffde7', '#ffcc02', '#e65100'], N=256
    )

    for idx, label in enumerate(labels):
        ax = axes_flat[idx]
        grid = grids[label]

        im = ax.imshow(grid, cmap=cmap, vmin=vmin, vmax=vmax,
                       aspect='auto', origin='upper')

        # Anotar cada celda
        for r in range(n_sa):
            for c in range(n_st):
                val = grid[r, c]
                if not np.isnan(val):
                    txt_color = 'white' if val > (vmin + 0.65 * (vmax - vmin)) else 'black'
                    ax.text(c, r, f'{val:.1f}',
                            ha='center', va='center',
                            fontsize=7.5, color=txt_color, fontweight='bold')

        season_color = SEASON_COLORS.get(label, '#333333')
        ax.set_title(label, fontsize=11, fontweight='bold', color=season_color)
        ax.set_xlabel('String', fontsize=9)
        ax.set_ylabel('Subarray', fontsize=9)
        ax.set_xticks(range(n_st))
        ax.set_xticklabels([f'ST{i+1}' for i in range(n_st)], fontsize=7)
        ax.set_yticks(range(n_sa))
        ax.set_yticklabels([f'SA{i+1}' for i in range(n_sa)], fontsize=7)

        cbar = plt.colorbar(im, ax=ax, shrink=0.85, pad=0.03)
        cbar.set_label('Sombra (%)', fontsize=8)
        cbar.ax.tick_params(labelsize=7)

    # Ocultar celda sobrante
    for idx in range(n_plots, n_rows * n_cols):
        axes_flat[idx].set_visible(False)

    fig.suptitle('Heatmap de sombra sobre geometría de paneles\n(Subarray × String)',
                 fontsize=13, fontweight='bold')
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  Heatmap → {out_path}')


# ── MAIN ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Post-proceso sombras SAM')
    parser.add_argument('--ts',   required=True, help='CSV series temporales')
    parser.add_argument('--diff', required=True, help='CSV sombra difusa')
    parser.add_argument('--out',  required=True, help='Carpeta de salida')
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
            os.path.join(args.out, 'heatmap_panel_geometry.png')
        )

    print('\nCompletado.')


if __name__ == '__main__':
    main()
