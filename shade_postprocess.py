#!/usr/bin/env python3
"""
shade_postprocess.py
Post-procesado de resultados del 3D Shade Calculator de SAM.

Modos de entrada:
  --ts <archivo.csv>           Un único CSV (array ≤32 paneles)
  --batches <b0.csv> <b1.csv>  Varios CSVs de batches que se fusionan antes del análisis

Salidas:
  summary_statistics.csv        Resumen estadístico por panel
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


MAX_SURFACES = 32   # límite SAM: 4 subarrays × 8 strings
MAX_STRINGS  = 8    # strings por subarray en SAM


def nice_tick_step(n, max_ticks=20):
    """Paso de ticks tal que n/paso ≤ max_ticks."""
    for step in [1, 2, 5, 10, 20, 25, 50, 100, 200]:
        if n <= step * max_ticks:
            return step
    return max(1, n // max_ticks)


def set_axis_ticks(ax, n_st, n_sa, st_labels, sa_labels, panel_w, panel_l, fs=7):
    step_st = nice_tick_step(n_st)
    step_sa = nice_tick_step(n_sa)
    ax.set_xticks([(st + 0.5) * panel_w for st in range(0, n_st, step_st)])
    ax.set_xticklabels([st_labels[st] for st in range(0, n_st, step_st)], fontsize=fs)
    ax.set_yticks([(sa + 0.5) * panel_l for sa in range(0, n_sa, step_sa)])
    ax.set_yticklabels([sa_labels[sa] for sa in range(0, n_sa, step_sa)], fontsize=fs)


def local_to_global_label(local_name, batch_id, cols):
    """
    SAM genera nombres de grupo a partir de Subarray/String locales del batch,
    ignorando el campo Group. Esta función convierte la etiqueta local al panel
    real en el array completo usando batch_id y el número de columnas.
    """
    parsed = parse_group(local_name)
    if parsed is None:
        return local_name
    local_sub, local_str = parsed                          # 1-based
    local_idx  = (local_sub - 1) * MAX_STRINGS + (local_str - 1)
    panel_idx  = batch_id * MAX_SURFACES + local_idx
    global_row = panel_idx // cols
    global_col = panel_idx % cols
    return f'SA{global_row + 1}_ST{global_col + 1}'


def load_timeseries(args):
    """
    Carga la serie temporal desde --ts (fichero único) o --batches (varios ficheros).
    En modo batches, renombra las columnas de etiquetas locales SAM a etiquetas globales.
    """
    if args.batches:
        dfs = []
        for path in args.batches:
            # Extraer batch_id del nombre de archivo (shade_batch<N>_...)
            m = re.search(r'shade_batch(\d+)_', os.path.basename(path))
            if not m:
                print(f'ERROR: no se puede extraer batch_id de {os.path.basename(path)}')
                print('El archivo debe llamarse shade_batch<N>_az<A>_inc<I>.csv')
                sys.exit(1)
            batch_id = int(m.group(1))

            df = pd.read_csv(path, index_col=0)
            df.columns = [local_to_global_label(c, batch_id, args.cols) for c in df.columns]
            dfs.append(df)
            print(f'  Batch {batch_id} cargado: {os.path.basename(path)}  ({len(df.columns)} paneles)')

        ts_raw = pd.concat(dfs, axis=1)
        duplicates = ts_raw.columns[ts_raw.columns.duplicated()].tolist()
        if duplicates:
            print(f'ERROR: etiquetas globales duplicadas: {duplicates}')
            print('Comprueba que rows, cols y batch_id sean coherentes en todos los batches.')
            sys.exit(1)
    else:
        ts_raw = pd.read_csv(args.ts, index_col=0)
        ts_raw.columns = [group_label(c) for c in ts_raw.columns]

    ts_raw.index = build_time_index(len(ts_raw))
    return ts_raw


# ── ESTADÍSTICAS ────────────────────────────────────────────────────────────

def compute_statistics(ts):
    return pd.DataFrame({
        'sombra_media_%':   (ts.mean()          * 100).round(2),
        'sombra_mediana_%': (ts.median()         * 100).round(2),
        'sombra_p90_%':     (ts.quantile(0.90)   * 100).round(2),
        'sombra_max_%':     (ts.max()             * 100).round(2),
        'h_con_sombra':     (ts > 0).sum(),
        'h_sombra_>10%':    (ts > 0.10).sum(),
        'h_sombra_>50%':    (ts > 0.50).sum(),
    })


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
    # sharex=False para que cada subplot muestre sus propias etiquetas de hora
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), sharey=True, sharex=False)
    axes = axes.flatten()

    n_groups = len(groups)
    cmap_grp = plt.get_cmap('tab20', n_groups)

    for ax, (season, df) in zip(axes, seas_avg.items()):
        for i, grp in enumerate(groups):
            ax.plot(df.index, df[grp] * 100,
                    color=cmap_grp(i), alpha=0.55, linewidth=1.1)

        array_mean = df[groups].mean(axis=1)
        ax.plot(array_mean.index, array_mean * 100,
                color='black', linewidth=2.2, linestyle='--')

        ax.set_title(season, fontsize=14, fontweight='bold', color='black')
        ax.set_xlim(0, 23)
        ax.set_ylim(0, 100)
        ax.set_xticks(range(0, 24, 3))
        ax.set_xticklabels([f'{h:02d}:00' for h in range(0, 24, 3)], fontsize=13)
        ax.set_xlabel('Hora', fontsize=14)
        ax.set_ylabel('Sombra (%)', fontsize=14)
        ax.tick_params(axis='y', labelsize=13)
        ax.grid(True, alpha=0.25)
        ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%g%%'))

    # Leyenda fija con dos entradas: paneles individuales y media del array
    legend_handles = [
        plt.Line2D([0], [0], color='steelblue', alpha=0.6, linewidth=1.5,
                   label=f'Panel individual (N={n_groups})'),
        plt.Line2D([0], [0], color='black', linewidth=2.2, linestyle='--',
                   label='Media del array'),
    ]
    fig.legend(handles=legend_handles, loc='lower center', ncol=2,
               fontsize=14, bbox_to_anchor=(0.5, -0.02),
               frameon=True, edgecolor='#cccccc')

    fig.suptitle('Sombra directa horaria promedio por estación',
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(out_path, dpi=600, bbox_inches='tight')
    plt.close()
    print(f'  Curvas estacionales → {out_path}')


# ── HEATMAP ──────────────────────────────────────────────────────────────────

def build_grids(seas_avg, ts, groups):
    """Dict {label: ndarray(n_sa × n_st)} con sombra media (%)."""
    coords = {g: parse_group(g) for g in groups}
    valid  = {g: c for g, c in coords.items() if c is not None}

    if not valid:
        print('AVISO: ningún grupo sigue el patrón SA_ST — heatmap omitido.')
        return None, 0, 0, [], []

    n_sa = max(c[0] for c in valid.values())
    n_st = max(c[1] for c in valid.values())

    grids = {}
    for season, df in seas_avg.items():
        grid = np.full((n_sa, n_st), np.nan)
        for g, (sa, st) in valid.items():
            grid[sa - 1, st - 1] = df[g].mean() * 100
        grids[season] = grid

    grid_ann = np.full((n_sa, n_st), np.nan)
    for g, (sa, st) in valid.items():
        grid_ann[sa - 1, st - 1] = ts[g].mean() * 100
    grids['Anual'] = grid_ann

    # Eliminar filas y columnas completamente vacías
    ref = grid_ann
    row_mask = ~np.all(np.isnan(ref), axis=1)
    col_mask = ~np.all(np.isnan(ref), axis=0)
    row_idx  = np.where(row_mask)[0]   # índices 0-based de filas con datos
    col_idx  = np.where(col_mask)[0]

    grids = {k: v[np.ix_(row_mask, col_mask)] for k, v in grids.items()}
    sa_labels = [f'SA{i+1}' for i in row_idx]
    st_labels = [f'ST{j+1}' for j in col_idx]

    return grids, len(row_idx), len(col_idx), sa_labels, st_labels


def season_slug(label):
    """Convierte etiqueta de estación/anual a slug ASCII para nombres de archivo."""
    replacements = str.maketrans('ñáéíóúÑÁÉÍÓÚ', 'naeiouNAEIOU')
    return label.lower().translate(replacements)


def plot_single_heatmap(grids, label, n_sa, n_st, sa_labels, st_labels, panel_w, panel_l, out_path):
    """Heatmap independiente para un label (estación o 'Anual'), escala común a todos los grids."""
    grid = grids[label]
    all_vals = np.concatenate([g[~np.isnan(g)] for g in grids.values()])
    vmin, vmax = np.percentile(all_vals, 2), np.percentile(all_vals, 98)

    cmap = LinearSegmentedColormap.from_list(
        'shade', ['#fffde7', '#ffcc02', '#e65100'], N=256
    )
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm   = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cell_in  = 0.80
    ax_w     = n_st * panel_w * cell_in
    ax_h     = n_sa * panel_l * cell_in

    fscale   = max(1.0, ax_w / 12.0)
    fs_title = round(10 * fscale)
    fs_label = round(9  * fscale)
    fs_tick  = round(8  * fscale)

    cbar_w   = 0.30 * fscale
    cbar_gap = 0.12 * fscale
    pad_l    = max(0.70, fs_label / 72 * 2.5)
    pad_r    = 0.20
    pad_top  = max(0.65, fs_title / 72 * 2.0)
    pad_bot  = max(0.55, fs_label / 72 * 2.0)

    fig_w = pad_l + ax_w + cbar_gap + cbar_w + pad_r
    fig_h = pad_top + ax_h + pad_bot

    fig = plt.figure(figsize=(fig_w, fig_h))

    ax_rect = [pad_l / fig_w, pad_bot / fig_h, ax_w / fig_w, ax_h / fig_h]
    ax = fig.add_axes(ax_rect)

    for sa in range(n_sa):
        for st in range(n_st):
            val = grid[sa, st]
            if np.isnan(val):
                continue
            rect = Rectangle(
                xy=(st * panel_w, sa * panel_l),
                width=panel_w, height=panel_l,
                facecolor=cmap(norm(val)),
                edgecolor='white', linewidth=0.8
            )
            ax.add_patch(rect)

    ax.set_xlim(0, n_st * panel_w)
    ax.set_ylim(0, n_sa * panel_l)
    title = 'Sombra anual media' if label == 'Anual' else f'Sombra media — {label}'
    ax.set_title(title, fontsize=fs_title, fontweight='bold', color='black',
                 pad=round(7 * fscale))
    ax.set_xlabel('String',   fontsize=fs_label, labelpad=round(4 * fscale))
    ax.set_ylabel('Subarray', fontsize=fs_label, labelpad=round(4 * fscale))
    set_axis_ticks(ax, n_st, n_sa, st_labels, sa_labels, panel_w, panel_l, fs=fs_tick)

    cbar_left = pad_l / fig_w + ax_w / fig_w + cbar_gap / fig_w
    cax = fig.add_axes([cbar_left, pad_bot / fig_h, cbar_w / fig_w, ax_h / fig_h])
    cbar = fig.colorbar(sm, cax=cax)
    cbar_label = 'Sombra media anual (%)' if label == 'Anual' else f'Sombra media {label} (%)'
    cbar.set_label(cbar_label, fontsize=fs_label)
    cbar.ax.tick_params(labelsize=fs_tick)

    plt.savefig(out_path, dpi=600, bbox_inches='tight')
    plt.close()
    print(f'  Heatmap {label} → {out_path}')


def plot_heatmaps(grids, n_sa, n_st, sa_labels, st_labels, panel_w, panel_l, out_path):
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
    cell_in  = 0.70           # pulgadas por metro de panel
    ax_w     = n_st * panel_w * cell_in
    ax_h     = n_sa * panel_l * cell_in

    fscale   = max(1.0, ax_w / 12.0)
    fs_super = round(10 * fscale)
    fs_title = round(8  * fscale)
    fs_label = round(8  * fscale)
    fs_tick  = round(7  * fscale)

    n_cols, n_rows = 3, 2
    cbar_w   = 0.25 * fscale
    cbar_gap = 0.10 * fscale
    col_gap  = 0.80 * fscale
    row_gap  = 0.70 * fscale
    pad_l    = max(0.55, fs_label / 72 * 2.5)
    pad_r    = 0.15
    pad_top  = max(0.55, fs_super / 72 * 2.5)
    pad_bot  = max(0.40, fs_label / 72 * 2.0)

    slot_w   = ax_w + cbar_w + cbar_gap
    fig_w    = pad_l + n_cols * slot_w + (n_cols - 1) * col_gap + pad_r
    fig_h    = pad_top + n_rows * ax_h + (n_rows - 1) * row_gap + pad_bot

    fig = plt.figure(figsize=(fig_w, fig_h))

    for idx, label in enumerate(labels):
        row = idx // n_cols
        col = idx % n_cols

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
        ax.set_ylim(0, n_sa * panel_l)   # SA1 abajo, coincide con eje Y de la escena 3D

        ax.set_title(label, fontsize=fs_title, fontweight='bold', color='black',
                     pad=round(4 * fscale))
        ax.set_xlabel('String',   fontsize=fs_label, labelpad=round(3 * fscale))
        ax.set_ylabel('Subarray', fontsize=fs_label, labelpad=round(3 * fscale))

        set_axis_ticks(ax, n_st, n_sa, st_labels, sa_labels, panel_w, panel_l, fs=fs_tick)

        cbar_left = ax_left + ax_w / fig_w + cbar_gap / fig_w
        cax = fig.add_axes([cbar_left, ax_bottom, cbar_w / fig_w, ax_h / fig_h])
        cbar = fig.colorbar(sm, cax=cax)
        cbar.set_label('Sombra (%)', fontsize=fs_label)
        cbar.ax.tick_params(labelsize=fs_tick)

    fig.text(0.5, 1 - pad_top / fig_h / 2, 'Heatmap de sombra — geometría de paneles',
             ha='center', va='top', fontsize=fs_super, fontweight='bold')

    plt.savefig(out_path, dpi=600, bbox_inches='tight')
    plt.close()
    print(f'  Heatmap → {out_path}')


# ── HEATMAPS DE IRRADIANCIA ──────────────────────────────────────────────────

def build_irradiance_grids(shade_grids, irradiance_ref):
    """
    Convierte grids de sombra (%) a irradiancia bloqueada (W/m²):
        irradiance_blocked = irradiance_ref × shade_fraction
    """
    return {label: grid * irradiance_ref / 100.0 for label, grid in shade_grids.items()}


def plot_irradiance_heatmaps(irr_grids, n_sa, n_st, sa_labels, st_labels,
                              panel_w, panel_l, irradiance_ref, out_path):
    """Heatmaps estacionales + anual de irradiancia bloqueada por sombras."""
    labels  = list(irr_grids.keys())
    n_plots = len(labels)

    all_vals = np.concatenate([g[~np.isnan(g)] for g in irr_grids.values()])
    vmin, vmax = np.percentile(all_vals, 2), np.percentile(all_vals, 98)

    cmap = LinearSegmentedColormap.from_list(
        'irradiance', ['#fff9c4', '#ff9800', '#b71c1c'], N=256
    )
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm   = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cell_in  = 0.70
    ax_w     = n_st * panel_w * cell_in
    ax_h     = n_sa * panel_l * cell_in

    fscale   = max(1.0, ax_w / 12.0)
    fs_super = round(10 * fscale)
    fs_title = round(8  * fscale)
    fs_label = round(8  * fscale)
    fs_tick  = round(7  * fscale)

    n_cols, n_rows = 3, 2
    cbar_w   = 0.25 * fscale
    cbar_gap = 0.10 * fscale
    col_gap  = 0.80 * fscale
    row_gap  = 0.70 * fscale
    pad_l    = max(0.55, fs_label / 72 * 2.5)
    pad_r    = 0.15
    pad_top  = max(0.65, fs_super / 72 * 2.5)
    pad_bot  = max(0.40, fs_label / 72 * 2.0)

    slot_w   = ax_w + cbar_w + cbar_gap
    fig_w    = pad_l + n_cols * slot_w + (n_cols - 1) * col_gap + pad_r
    fig_h    = pad_top + n_rows * ax_h + (n_rows - 1) * row_gap + pad_bot

    fig = plt.figure(figsize=(fig_w, fig_h))

    for idx, label in enumerate(labels):
        row = idx // n_cols
        col = idx % n_cols

        ax_left   = (pad_l + col * (slot_w + col_gap)) / fig_w
        ax_bottom = (pad_bot + (n_rows - 1 - row) * (ax_h + row_gap)) / fig_h
        ax_rect   = [ax_left, ax_bottom, ax_w / fig_w, ax_h / fig_h]

        ax = fig.add_axes(ax_rect)
        grid = irr_grids[label]

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
        ax.set_ylim(0, n_sa * panel_l)

        ax.set_title(label, fontsize=fs_title, fontweight='bold', color='black',
                     pad=round(4 * fscale))
        ax.set_xlabel('String',   fontsize=fs_label, labelpad=round(3 * fscale))
        ax.set_ylabel('Subarray', fontsize=fs_label, labelpad=round(3 * fscale))

        set_axis_ticks(ax, n_st, n_sa, st_labels, sa_labels, panel_w, panel_l, fs=fs_tick)

        cbar_left = ax_left + ax_w / fig_w + cbar_gap / fig_w
        cax = fig.add_axes([cbar_left, ax_bottom, cbar_w / fig_w, ax_h / fig_h])
        cbar = fig.colorbar(sm, cax=cax)
        cbar.set_label('Irradiancia bloqueada (W/m²)', fontsize=fs_label)
        cbar.ax.tick_params(labelsize=fs_tick)

    fig.text(
        0.5, 1 - pad_top / fig_h / 2,
        f'Irradiancia bloqueada por sombras — ref. {irradiance_ref:.0f} W/m²',
        ha='center', va='top', fontsize=fs_super, fontweight='bold'
    )

    plt.savefig(out_path, dpi=600, bbox_inches='tight')
    plt.close()
    print(f'  Heatmap irradiancia → {out_path}')


def plot_irradiance_single_heatmap(irr_grids, label, n_sa, n_st, sa_labels, st_labels,
                                    panel_w, panel_l, irradiance_ref, out_path):
    """Heatmap independiente de irradiancia bloqueada para un label (estación o 'Anual')."""
    grid = irr_grids[label]
    all_vals = np.concatenate([g[~np.isnan(g)] for g in irr_grids.values()])
    vmin, vmax = np.percentile(all_vals, 2), np.percentile(all_vals, 98)

    cmap = LinearSegmentedColormap.from_list(
        'irradiance', ['#fff9c4', '#ff9800', '#b71c1c'], N=256
    )
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm   = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cell_in  = 0.80
    ax_w     = n_st * panel_w * cell_in
    ax_h     = n_sa * panel_l * cell_in

    fscale   = max(1.0, ax_w / 12.0)
    fs_title = round(10 * fscale)
    fs_label = round(9  * fscale)
    fs_tick  = round(8  * fscale)

    cbar_w   = 0.30 * fscale
    cbar_gap = 0.12 * fscale
    pad_l    = max(0.70, fs_label / 72 * 2.5)
    pad_r    = 0.20
    pad_top  = max(0.65, fs_title / 72 * 2.0)
    pad_bot  = max(0.55, fs_label / 72 * 2.0)

    fig_w = pad_l + ax_w + cbar_gap + cbar_w + pad_r
    fig_h = pad_top + ax_h + pad_bot

    fig = plt.figure(figsize=(fig_w, fig_h))

    ax_rect = [pad_l / fig_w, pad_bot / fig_h, ax_w / fig_w, ax_h / fig_h]
    ax = fig.add_axes(ax_rect)

    for sa in range(n_sa):
        for st in range(n_st):
            val = grid[sa, st]
            if np.isnan(val):
                continue
            rect = Rectangle(
                xy=(st * panel_w, sa * panel_l),
                width=panel_w, height=panel_l,
                facecolor=cmap(norm(val)),
                edgecolor='white', linewidth=0.8
            )
            ax.add_patch(rect)

    ax.set_xlim(0, n_st * panel_w)
    ax.set_ylim(0, n_sa * panel_l)
    period = 'anual media' if label == 'Anual' else label
    ax.set_title(
        f'Irradiancia bloqueada {period} — ref. {irradiance_ref:.0f} W/m²',
        fontsize=fs_title, fontweight='bold', color='black', pad=round(7 * fscale)
    )
    ax.set_xlabel('String',   fontsize=fs_label, labelpad=round(4 * fscale))
    ax.set_ylabel('Subarray', fontsize=fs_label, labelpad=round(4 * fscale))
    set_axis_ticks(ax, n_st, n_sa, st_labels, sa_labels, panel_w, panel_l, fs=fs_tick)

    cbar_left = pad_l / fig_w + ax_w / fig_w + cbar_gap / fig_w
    cax = fig.add_axes([cbar_left, pad_bot / fig_h, cbar_w / fig_w, ax_h / fig_h])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label('Irradiancia bloqueada media (W/m²)', fontsize=fs_label)
    cbar.ax.tick_params(labelsize=fs_tick)

    plt.savefig(out_path, dpi=600, bbox_inches='tight')
    plt.close()
    print(f'  Heatmap irradiancia {label} → {out_path}')


# ── MAIN ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Post-proceso sombras SAM')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--ts',      help='CSV series temporales (modo single-batch)')
    group.add_argument('--batches', nargs='+', metavar='CSV',
                       help='CSVs de múltiples batches a combinar')
    parser.add_argument('--out',    required=True, help='Carpeta de salida')
    parser.add_argument('--width',  type=float, default=1.0, help='Ancho del panel (m)')
    parser.add_argument('--length', type=float, default=1.0, help='Largo del panel (m)')
    parser.add_argument('--cols',       type=int,   default=8,
                        help='Número de columnas del array completo (necesario con --batches)')
    parser.add_argument('--irradiance', type=float, default=1000.0,
                        help='Irradiancia de referencia en W/m² para heatmaps de irradiancia (default: 1000)')
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # ── Carga y fusión ───────────────────────────────────────
    mode = f'{len(args.batches)} batches' if args.batches else '1 archivo'
    print(f'Cargando datos ({mode})...')
    ts_raw = load_timeseries(args)
    groups = list(ts_raw.columns)
    print(f'  {len(groups)} paneles en total: {groups[:4]}{"..." if len(groups) > 4 else ""}')

    # Los valores del time series están en % (0-100) → normalizar a fracción 0-1
    if ts_raw.max().max() > 1.5:
        ts_raw = ts_raw / 100.0
        print('  Valores convertidos de % a fracción (0-1).')

    # ── Estadísticas ─────────────────────────────────────────
    print('Calculando estadísticas...')
    stats = compute_statistics(ts_raw)
    stats_path = os.path.join(args.out, 'summary_statistics.csv')
    stats.to_csv(stats_path)
    print(f'  Estadísticas → {stats_path}')

    # ── Curvas estacionales ───────────────────────────────────
    print('Generando curvas horarias por estación...')
    seas_avg = seasonal_hourly_avg(ts_raw)
    plot_seasonal_curves(
        seas_avg, groups,
        os.path.join(args.out, 'seasonal_hourly_curves.png')
    )

    # ── Heatmaps ─────────────────────────────────────────────
    print('Generando heatmaps...')
    grids, n_sa, n_st, sa_labels, st_labels = build_grids(seas_avg, ts_raw, groups)
    if grids:
        plot_heatmaps(
            grids, n_sa, n_st, sa_labels, st_labels,
            args.width, args.length,
            os.path.join(args.out, 'heatmap_panel_geometry.png')
        )
        for lbl in list(SEASONS.keys()) + ['Anual']:
            slug = season_slug(lbl)
            plot_single_heatmap(
                grids, lbl, n_sa, n_st, sa_labels, st_labels,
                args.width, args.length,
                os.path.join(args.out, f'heatmap_{slug}.png')
            )

        # ── Heatmaps de irradiancia ───────────────────────────
        print(f'Generando heatmaps de irradiancia (ref. {args.irradiance:.0f} W/m²)...')
        irr_grids = build_irradiance_grids(grids, args.irradiance)
        plot_irradiance_heatmaps(
            irr_grids, n_sa, n_st, sa_labels, st_labels,
            args.width, args.length, args.irradiance,
            os.path.join(args.out, 'heatmap_irradiance.png')
        )
        for lbl in list(SEASONS.keys()) + ['Anual']:
            slug = season_slug(lbl)
            plot_irradiance_single_heatmap(
                irr_grids, lbl, n_sa, n_st, sa_labels, st_labels,
                args.width, args.length, args.irradiance,
                os.path.join(args.out, f'heatmap_irradiance_{slug}.png')
            )

    print('\nCompletado.')


if __name__ == '__main__':
    main()
