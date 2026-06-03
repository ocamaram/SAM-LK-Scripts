#!/usr/bin/env python3
"""
shade_postprocess.py
Post-procesado de resultados del 3D Shade Calculator de SAM.

Modos de entrada:
  --ts <archivo.csv>           Un único CSV (array ≤32 paneles)
  --batches <b0.csv> <b1.csv>  Varios CSVs de batches que se fusionan antes del análisis

Opciones de irradiancia:
  --irradiance-ts <archivo.csv>  Serie temporal de irradiancia horaria (W/m²) exportada
                                  desde SAM. Si se omite, se usa --irradiance constante.
  --irradiance-col <nombre>      Columna del CSV de irradiancia (default: primera columna).
  --irradiance <W/m²>            Irradiancia constante de referencia (solo si no hay --irradiance-ts).

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


def _fmt_m(v):
    """Formatea un valor en metros: entero si es entero, un decimal si no."""
    return f'{v:.0f}' if v % 1 == 0 else f'{v:.1f}'


def set_axis_ticks(ax, n_st, n_sa, panel_w, panel_l, spacing_x=None, spacing_y=None, fs=7):
    """Dibuja ticks con distancias en metros en lugar de etiquetas SA/ST."""
    sp_y = spacing_y if spacing_y is not None else panel_w   # separación entre columnas
    sp_x = spacing_x if spacing_x is not None else panel_l   # separación entre filas

    step_st = nice_tick_step(n_st)
    step_sa = nice_tick_step(n_sa)

    x_idxs = list(range(0, n_st, step_st)) + [n_st]
    y_idxs = list(range(0, n_sa, step_sa)) + [n_sa]

    ax.set_xticks([st * panel_w for st in x_idxs])
    ax.set_xticklabels([_fmt_m(st * sp_y) for st in x_idxs], fontsize=fs)
    ax.set_yticks([sa * panel_l for sa in y_idxs])
    ax.set_yticklabels([_fmt_m(sa * sp_x) for sa in y_idxs], fontsize=fs)


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


def load_irradiance_ts(path, col_name=None):
    """
    Carga la serie temporal de irradiancia horaria desde un CSV exportado por SAM.
    La primera columna se trata como índice (SAM exporta horas como primera columna).
    Devuelve una Series con los valores en las unidades originales del CSV.
    """
    df = pd.read_csv(path, index_col=0)
    if col_name:
        if col_name not in df.columns:
            print(f'ERROR: columna "{col_name}" no encontrada en {path}.')
            print(f'  Columnas disponibles: {list(df.columns)}')
            sys.exit(1)
        irr = df[col_name]
    else:
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) == 0:
            print(f'ERROR: no se encontraron columnas numéricas en {path}.')
            sys.exit(1)
        if len(numeric_cols) > 1:
            print(f'  Irradiancia: columna auto-seleccionada "{numeric_cols[0]}" '
                  f'(usa --irradiance-col para especificar otra).')
        irr = df[numeric_cols[0]]

    irr = irr.fillna(0.0).clip(lower=0.0)
    return irr.reset_index(drop=True)


def group_batches_by_orientation(batch_paths):
    """
    Agrupa los ficheros de batch por orientación (sufijo az<A>_inc<I>).
    Devuelve un dict ordenado {az_inc: [path, ...]} con los ficheros de cada grupo.
    Si no hay sufijo de orientación, todos van a la clave ''.
    """
    orient_re = re.compile(r'shade_batch\d+_(az\d+_inc\d+)\.csv$')
    groups = {}
    for path in batch_paths:
        m = orient_re.search(os.path.basename(path))
        key = m.group(1) if m else ''
        groups.setdefault(key, []).append(path)
    return dict(sorted(groups.items()))


def load_timeseries(args, batch_paths=None):
    """
    Carga la serie temporal desde --ts (fichero único) o una lista de batch CSVs.
    batch_paths sobreescribe args.batches cuando se llama por orientación.
    """
    paths = batch_paths if batch_paths is not None else args.batches
    if paths:
        dfs = []
        for path in paths:
            m = re.search(r'shade_batch(\d+)_', os.path.basename(path))
            if not m:
                print(f'ERROR: no se puede extraer batch_id de {os.path.basename(path)}')
                sys.exit(1)
            batch_id = int(m.group(1))

            if os.path.getsize(path) == 0:
                print(f'  Batch {batch_id} omitido: {os.path.basename(path)}  (fichero vacío)')
                continue

            df = pd.read_csv(path, index_col=0)
            if df.empty or len(df.columns) == 0:
                print(f'  Batch {batch_id} omitido: {os.path.basename(path)}  (sin datos)')
                continue

            df.columns = [local_to_global_label(c, batch_id, args.cols) for c in df.columns]
            dfs.append(df)
            print(f'  Batch {batch_id} cargado: {os.path.basename(path)}  ({len(df.columns)} paneles)')

        if not dfs:
            return None

        ts_raw = pd.concat(dfs, axis=1)
        duplicates = ts_raw.columns[ts_raw.columns.duplicated()].tolist()
        if duplicates:
            print(f'ERROR: etiquetas globales duplicadas: {duplicates[:8]}...')
            print('Comprueba que --cols sea correcto.')
            sys.exit(1)
    else:
        ts_raw = pd.read_csv(args.ts, index_col=0)
        ts_raw.columns = [group_label(c) for c in ts_raw.columns]

    ts_raw.index = build_time_index(len(ts_raw))
    return ts_raw


# ── ESTADÍSTICAS ────────────────────────────────────────────────────────────

def compute_statistics(ts, irr_ts=None):
    stats = {
        'sombra_media_%':   (ts.mean()          * 100).round(2),
        'sombra_mediana_%': (ts.median()         * 100).round(2),
        'sombra_p90_%':     (ts.quantile(0.90)   * 100).round(2),
        'sombra_max_%':     (ts.max()             * 100).round(2),
        'h_con_sombra':     (ts > 0).sum(),
        'h_sombra_>10%':    (ts > 0.10).sum(),
        'h_sombra_>50%':    (ts > 0.50).sum(),
    }

    if irr_ts is not None:
        irr_arr = irr_ts.values
        irr_total = irr_arr.sum()   # Wh/m² totales del período
        if irr_total > 0:
            blocked = {col: float((irr_arr * ts[col].values).sum()) for col in ts.columns}
            stats['sombra_irr_ponderada_%'] = pd.Series(
                {col: round(v / irr_total * 100, 2) for col, v in blocked.items()}
            )
            stats['energia_bloqueada_kWh_m2'] = pd.Series(
                {col: round(v / 1000, 2) for col, v in blocked.items()}
            )
            stats['energia_incidente_kWh_m2'] = pd.Series(
                {col: round((irr_total - v) / 1000, 2) for col, v in blocked.items()}
            )

    return pd.DataFrame(stats)


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


def plot_seasonal_curves(seas_avg, groups, out_path, dpi=150):
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
    plt.savefig(out_path, dpi=dpi, bbox_inches='tight')
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


def plot_single_heatmap(grids, label, n_sa, n_st, sa_labels, st_labels, panel_w, panel_l, out_path,
                         vmin=None, vmax=None, dpi=150, spacing_x=None, spacing_y=None):
    """Heatmap independiente para un label (estación o 'Anual'), escala común a todos los grids."""
    grid = grids[label]
    if vmin is None or vmax is None:
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
    ax.set_xlabel('x (m)', fontsize=fs_label, labelpad=round(4 * fscale))
    ax.set_ylabel('y (m)', fontsize=fs_label, labelpad=round(4 * fscale))
    set_axis_ticks(ax, n_st, n_sa, panel_w, panel_l,
                   spacing_x=spacing_x, spacing_y=spacing_y, fs=fs_tick)

    cbar_left = pad_l / fig_w + ax_w / fig_w + cbar_gap / fig_w
    cax = fig.add_axes([cbar_left, pad_bot / fig_h, cbar_w / fig_w, ax_h / fig_h])
    cbar = fig.colorbar(sm, cax=cax)
    cbar_label = 'Sombra media anual (%)' if label == 'Anual' else f'Sombra media {label} (%)'
    cbar.set_label(cbar_label, fontsize=fs_label)
    cbar.ax.tick_params(labelsize=fs_tick)

    plt.savefig(out_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f'  Heatmap {label} → {out_path}')


def plot_heatmaps(grids, n_sa, n_st, sa_labels, st_labels, panel_w, panel_l, out_path, dpi=150,
                  spacing_x=None, spacing_y=None):
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
        ax.set_xlabel('x (m)', fontsize=fs_label, labelpad=round(3 * fscale))
        ax.set_ylabel('y (m)', fontsize=fs_label, labelpad=round(3 * fscale))

        set_axis_ticks(ax, n_st, n_sa, panel_w, panel_l,
                       spacing_x=spacing_x, spacing_y=spacing_y, fs=fs_tick)

        cbar_left = ax_left + ax_w / fig_w + cbar_gap / fig_w
        cax = fig.add_axes([cbar_left, ax_bottom, cbar_w / fig_w, ax_h / fig_h])
        cbar = fig.colorbar(sm, cax=cax)
        cbar.set_label('Sombra (%)', fontsize=fs_label)
        cbar.ax.tick_params(labelsize=fs_tick)

    fig.text(0.5, 1 - pad_top / fig_h / 2, 'Heatmap de sombra — geometría de paneles',
             ha='center', va='top', fontsize=fs_super, fontweight='bold')

    plt.savefig(out_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f'  Heatmap → {out_path}')


# ── HEATMAPS DE IRRADIANCIA ──────────────────────────────────────────────────

def build_irr_received_ts(ts, irr_ts):
    """
    Calcula la irradiancia recibida hora a hora para cada panel:
        irr_received(t) = irradiance(t) × (1 - shade_fraction(t))
    Devuelve un DataFrame con las mismas dimensiones que ts (en W/m²).
    """
    irr_arr = irr_ts.values.reshape(-1, 1)
    return pd.DataFrame(
        irr_arr * (1.0 - ts.values),
        index=ts.index,
        columns=ts.columns,
    )


def build_irradiance_grids_from_ts(ts, irr_ts, groups):
    """
    Construye grids de irradiancia media recibida (W/m²) ponderando la fracción de
    sombra hora a hora por la irradiancia real.  Devuelve la misma estructura que
    build_grids(): (grids_dict, n_sa, n_st, sa_labels, st_labels).
    """
    irr_received_ts = build_irr_received_ts(ts, irr_ts)

    coords = {g: parse_group(g) for g in groups}
    valid  = {g: c for g, c in coords.items() if c is not None}
    if not valid:
        print('AVISO: ningún grupo sigue el patrón SA_ST — heatmap irradiancia omitido.')
        return None, 0, 0, [], []

    n_sa = max(c[0] for c in valid.values())
    n_st = max(c[1] for c in valid.values())

    grids = {}
    for season, months in SEASONS.items():
        mask = irr_received_ts.index.month.isin(months)
        seasonal_mean = irr_received_ts[mask].mean()
        grid = np.full((n_sa, n_st), np.nan)
        for g, (sa, st) in valid.items():
            grid[sa - 1, st - 1] = seasonal_mean[g]
        grids[season] = grid

    grid_ann = np.full((n_sa, n_st), np.nan)
    for g, (sa, st) in valid.items():
        grid_ann[sa - 1, st - 1] = irr_received_ts[g].mean()
    grids['Anual'] = grid_ann

    ref = grid_ann
    row_mask = ~np.all(np.isnan(ref), axis=1)
    col_mask = ~np.all(np.isnan(ref), axis=0)
    row_idx  = np.where(row_mask)[0]
    col_idx  = np.where(col_mask)[0]

    grids = {k: v[np.ix_(row_mask, col_mask)] for k, v in grids.items()}
    sa_labels = [f'SA{i+1}' for i in row_idx]
    st_labels = [f'ST{j+1}' for j in col_idx]

    return grids, len(row_idx), len(col_idx), sa_labels, st_labels


def build_irradiance_grids(shade_grids, irradiance_ref):
    """
    Convierte grids de sombra (%) a irradiancia recibida (W/m²):
        irradiance_received = irradiance_ref × (1 - shade_fraction)
    Fallback cuando no se dispone de serie temporal de irradiancia.
    """
    return {label: (1.0 - grid / 100.0) * irradiance_ref for label, grid in shade_grids.items()}


def plot_irradiance_heatmaps(irr_grids, n_sa, n_st, sa_labels, st_labels,
                              panel_w, panel_l, irradiance_ref=None, out_path=None, dpi=150,
                              spacing_x=None, spacing_y=None):
    """Heatmaps estacionales + anual de irradiancia recibida (W/m²)."""
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
        ax.set_xlabel('x (m)', fontsize=fs_label, labelpad=round(3 * fscale))
        ax.set_ylabel('y (m)', fontsize=fs_label, labelpad=round(3 * fscale))

        set_axis_ticks(ax, n_st, n_sa, panel_w, panel_l,
                       spacing_x=spacing_x, spacing_y=spacing_y, fs=fs_tick)

        cbar_left = ax_left + ax_w / fig_w + cbar_gap / fig_w
        cax = fig.add_axes([cbar_left, ax_bottom, cbar_w / fig_w, ax_h / fig_h])
        cbar = fig.colorbar(sm, cax=cax)
        cbar.set_label('Irradiancia recibida (W/m²)', fontsize=fs_label)
        cbar.ax.tick_params(labelsize=fs_tick)

    if irradiance_ref is not None:
        suptitle = f'Irradiancia recibida — ref. {irradiance_ref:.0f} W/m²'
    else:
        suptitle = 'Irradiancia media recibida (ponderada por irradiancia horaria real)'
    fig.text(
        0.5, 1 - pad_top / fig_h / 2,
        suptitle,
        ha='center', va='top', fontsize=fs_super, fontweight='bold'
    )

    plt.savefig(out_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f'  Heatmap irradiancia → {out_path}')


def plot_irradiance_single_heatmap(irr_grids, label, n_sa, n_st, sa_labels, st_labels,
                                    panel_w, panel_l, irradiance_ref=None, out_path=None,
                                    vmin=None, vmax=None, dpi=150, spacing_x=None, spacing_y=None):
    """Heatmap independiente de irradiancia recibida para un label (estación o 'Anual')."""
    grid = irr_grids[label]
    if vmin is None or vmax is None:
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
    if irradiance_ref is not None:
        title_str = f'Irradiancia recibida {period} — ref. {irradiance_ref:.0f} W/m²'
    else:
        title_str = f'Irradiancia media recibida — {period}'
    ax.set_title(title_str, fontsize=fs_title, fontweight='bold', color='black',
                 pad=round(7 * fscale))
    ax.set_xlabel('x (m)', fontsize=fs_label, labelpad=round(4 * fscale))
    ax.set_ylabel('y (m)', fontsize=fs_label, labelpad=round(4 * fscale))
    set_axis_ticks(ax, n_st, n_sa, panel_w, panel_l,
                   spacing_x=spacing_x, spacing_y=spacing_y, fs=fs_tick)

    cbar_left = pad_l / fig_w + ax_w / fig_w + cbar_gap / fig_w
    cax = fig.add_axes([cbar_left, pad_bot / fig_h, cbar_w / fig_w, ax_h / fig_h])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label('Irradiancia recibida media (W/m²)', fontsize=fs_label)
    cbar.ax.tick_params(labelsize=fs_tick)

    plt.savefig(out_path, dpi=dpi, bbox_inches='tight')
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
    parser.add_argument('--width',     type=float, default=1.0, help='Ancho del panel (m)')
    parser.add_argument('--length',    type=float, default=1.0, help='Largo del panel (m)')
    parser.add_argument('--spacing-x', dest='spacing_x', type=float, default=None,
                        help='Separación entre filas (m). Default: igual a --length.')
    parser.add_argument('--spacing-y', dest='spacing_y', type=float, default=None,
                        help='Separación entre columnas (m). Default: igual a --width.')
    parser.add_argument('--cols',       type=int,   default=8,
                        help='Número de columnas del array completo (necesario con --batches)')
    parser.add_argument('--irradiance-ts', dest='irradiance_ts', metavar='CSV',
                        help='CSV de irradiancia horaria exportado desde SAM. '
                             'Si se especifica, los heatmaps usan irradiancia ponderada real.')
    parser.add_argument('--irradiance-col', dest='irradiance_col', metavar='COLUMNA',
                        help='Nombre de la columna de irradiancia en --irradiance-ts '
                             '(default: primera columna numérica).')
    parser.add_argument('--irradiance-area', dest='irradiance_area', type=float, metavar='M2',
                        help='Área total del array en m². Si se especifica, convierte los valores '
                             'de --irradiance-ts de kW a W/m² dividiéndolos por esta área.')
    parser.add_argument('--irradiance', type=float, default=1000.0,
                        help='Irradiancia constante de referencia en W/m² (solo si no se usa '
                             '--irradiance-ts, default: 1000)')
    parser.add_argument('--dpi', type=int, default=150,
                        help='Resolución de las imágenes de salida (default: 150). Usa 300 para calidad de impresión.')
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # ── Carga de irradiancia horaria (opcional, común a todas las orientaciones) ──
    irr_ts_raw = None
    if args.irradiance_ts:
        print(f'Cargando irradiancia horaria desde {args.irradiance_ts}...')
        irr_ts_raw = load_irradiance_ts(args.irradiance_ts, args.irradiance_col)
        if args.irradiance_area:
            if args.irradiance_area <= 0:
                print('ERROR: --irradiance-area debe ser un valor positivo.')
                sys.exit(1)
            irr_ts_raw = irr_ts_raw * 1000.0 / args.irradiance_area
            print(f'  Conversión kW → W/m²: ÷ {args.irradiance_area} m²')

    # ── Determinar grupos de orientación ─────────────────────
    if args.batches:
        orient_groups = group_batches_by_orientation(args.batches)
        multi = len(orient_groups) > 1
    else:
        orient_groups = {'': None}   # modo --ts
        multi = False

    all_stats = []
    all_ts_parts = []   # para curva estacional combinada

    for orient, batch_paths in orient_groups.items():
        prefix = f'{orient}_' if (multi and orient) else ''
        header = f'Orientación {orient}' if (multi and orient) else 'datos'
        n_files = len(batch_paths) if batch_paths else 1
        print(f'\nCargando {header} ({n_files} batches)...')

        ts = load_timeseries(args, batch_paths=batch_paths)
        if ts is None:
            print(f'  Sin datos para {orient}, omitido.')
            continue
        panels = list(ts.columns)
        print(f'  {len(panels)} paneles: {panels[:4]}{"..." if len(panels) > 4 else ""}')

        if ts.max().max() > 1.5:
            ts = ts / 100.0
            print('  Valores convertidos de % a fracción (0-1).')

        # Alinear irradiancia con este TS
        irr_ts = None
        if irr_ts_raw is not None:
            if len(irr_ts_raw) != len(ts):
                print(f'ERROR: CSV de irradiancia ({len(irr_ts_raw)} filas) ≠ '
                      f'sombras ({len(ts)} filas).')
                sys.exit(1)
            irr_ts = irr_ts_raw.copy()
            irr_ts.index = ts.index
            print(f'  Irradiancia media en horas de sol: {irr_ts[irr_ts > 0].mean():.1f} W/m²')

        all_ts_parts.append(ts)

        # ── Estadísticas ─────────────────────────────────────
        print('  Calculando estadísticas...')
        stats = compute_statistics(ts, irr_ts=irr_ts)
        all_stats.append(stats)

        # ── Curvas estacionales (por orientación si hay varias) ──
        seas_avg = seasonal_hourly_avg(ts)
        if multi:
            plot_seasonal_curves(
                seas_avg, panels,
                os.path.join(args.out, f'{prefix}seasonal_hourly_curves.png'),
                dpi=args.dpi,
            )

        # ── Heatmaps de sombra ────────────────────────────────
        print('  Generando heatmaps de sombra...')
        grids, n_sa, n_st, sa_labels, st_labels = build_grids(seas_avg, ts, panels)
        if grids:
            sp_x = args.spacing_x
            sp_y = args.spacing_y
            plot_heatmaps(
                grids, n_sa, n_st, sa_labels, st_labels,
                args.width, args.length,
                os.path.join(args.out, f'{prefix}heatmap_panel_geometry.png'),
                dpi=args.dpi, spacing_x=sp_x, spacing_y=sp_y,
            )
            all_shade_vals = np.concatenate([g[~np.isnan(g)] for g in grids.values()])
            shade_vmin = float(np.percentile(all_shade_vals, 2))
            shade_vmax = float(np.percentile(all_shade_vals, 98))
            del all_shade_vals
            for lbl in list(SEASONS.keys()) + ['Anual']:
                slug = season_slug(lbl)
                plot_single_heatmap(
                    grids, lbl, n_sa, n_st, sa_labels, st_labels,
                    args.width, args.length,
                    os.path.join(args.out, f'{prefix}heatmap_{slug}.png'),
                    vmin=shade_vmin, vmax=shade_vmax, dpi=args.dpi,
                    spacing_x=sp_x, spacing_y=sp_y,
                )

            # ── Heatmaps de irradiancia ───────────────────────
            if irr_ts is not None:
                print('  Generando heatmaps de irradiancia ponderada...')
                irr_grids, irr_n_sa, irr_n_st, irr_sa_lbl, irr_st_lbl = \
                    build_irradiance_grids_from_ts(ts, irr_ts, panels)
                irr_ref_arg = None
            else:
                print(f'  Generando heatmaps de irradiancia (ref. {args.irradiance:.0f} W/m²)...')
                irr_grids = build_irradiance_grids(grids, args.irradiance)
                irr_n_sa, irr_n_st = n_sa, n_st
                irr_sa_lbl, irr_st_lbl = sa_labels, st_labels
                irr_ref_arg = args.irradiance

            if irr_grids:
                plot_irradiance_heatmaps(
                    irr_grids, irr_n_sa, irr_n_st, irr_sa_lbl, irr_st_lbl,
                    args.width, args.length, irr_ref_arg,
                    os.path.join(args.out, f'{prefix}heatmap_irradiance.png'),
                    dpi=args.dpi, spacing_x=sp_x, spacing_y=sp_y,
                )
                all_irr_vals = np.concatenate([g[~np.isnan(g)] for g in irr_grids.values()])
                irr_vmin = float(np.percentile(all_irr_vals, 2))
                irr_vmax = float(np.percentile(all_irr_vals, 98))
                del all_irr_vals
                for lbl in list(SEASONS.keys()) + ['Anual']:
                    slug = season_slug(lbl)
                    plot_irradiance_single_heatmap(
                        irr_grids, lbl, irr_n_sa, irr_n_st, irr_sa_lbl, irr_st_lbl,
                        args.width, args.length, irr_ref_arg,
                        os.path.join(args.out, f'{prefix}heatmap_irradiance_{slug}.png'),
                        vmin=irr_vmin, vmax=irr_vmax, dpi=args.dpi,
                        spacing_x=sp_x, spacing_y=sp_y,
                    )

    # ── Estadísticas combinadas ───────────────────────────────
    print('\nGuardando estadísticas combinadas...')
    combined_stats = pd.concat(all_stats)
    stats_path = os.path.join(args.out, 'summary_statistics.csv')
    combined_stats.to_csv(stats_path)
    print(f'  Estadísticas → {stats_path}')

    # ── Curva estacional combinada (todos los paneles) ────────
    if all_ts_parts:
        ts_all = pd.concat(all_ts_parts, axis=1)
        print('Generando curvas horarias estacionales (todos los paneles)...')
        seas_all = seasonal_hourly_avg(ts_all)
        plot_seasonal_curves(
            seas_all, list(ts_all.columns),
            os.path.join(args.out, 'seasonal_hourly_curves.png'),
            dpi=args.dpi,
        )

    print('\nCompletado.')


if __name__ == '__main__':
    main()
