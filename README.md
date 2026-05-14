# SAM LK Scripts

Librería de scripts LK para [NREL System Advisor Model (SAM)](https://sam.nrel.gov/), orientada al análisis de sombreado 3D y simulación de sistemas fotovoltaicos.

---

## Instalación de dependencias Python

```bash
pip install -r requirements.txt
```

---

## Scripts disponibles

### `SAM_3D_Shading.lk` — Generador de superficies activas y análisis de sombras

Genera una cuadrícula de superficies activas en el **3D Shade Calculator** de SAM, lanza el análisis de series temporales y exporta los porcentajes de sombra por subarray y string a CSV.

#### Características

- Crea `rows × cols` superficies activas con orientación (azimut + inclinación) configurable.
- Asigna automáticamente cada superficie a su subarray y string correspondientes.
- **Preserva la geometría de obstáculos existente** (edificios, árboles, cajas…) — solo elimina las superficies activas previas.
- Ejecuta el análisis de sombra directa en series temporales (`direct_shade`) con el paso temporal deseado.
- Exporta los resultados a CSV y **lanza automáticamente** `shade_postprocess.py` al finalizar.

#### Parámetros

| Parámetro | Descripción | Valor por defecto |
|---|---|---|
| `rows` | Número de filas de superficies | `4` |
| `cols` | Número de columnas de superficies | `8` |
| `x0`, `y0`, `z0` | Origen de la cuadrícula (m) | `0.0` |
| `spacing_x` | Separación entre filas (m) | `1.0` |
| `spacing_y` | Separación entre columnas (m) | `1.0` |
| `azimut` | Orientación solar de la superficie — 0=N, 90=E, 180=S, 270=O (°) | `180` |
| `inclinacion` | Inclinación — 0=horizontal, 90=vertical/fachada (°) | `30` |
| `timestep_min` | Paso temporal del análisis de series temporales (min) | `60` |

> **Límites SAM:** máximo 4 subarrays × 8 strings = 32 superficies activas.

#### Uso

1. Abre el **3D Shade Calculator** desde la pestaña *Shading* de tu caso SAM.
2. Configura la **ubicación** en la pestaña *Location* (necesaria para el cálculo solar).
3. Añade opcionalmente los **obstáculos** de sombra (edificios, árboles…) en la pestaña *3D Scene*.
4. Abre la pestaña *Scripting*, carga `SAM_3D_Shading.lk` y ajusta los parámetros de la sección `PARÁMETROS`.
5. Ejecuta el script. Se abrirá un diálogo para elegir la carpeta de exportación.

#### Archivos de salida

| Archivo | Contenido |
|---|---|
| `shade_timeseries_az<A>_inc<I>.csv` | 8 760 filas × n_grupos columnas con el factor de sombra directa horario (0–1) por subarray+string |
| `summary_statistics.csv` | Resumen estadístico por grupo (generado por `shade_postprocess.py`) |
| `seasonal_hourly_curves.png` | Curvas horarias promedio de sombra por estación |
| `heatmap_panel_geometry.png` | Heatmap estacional y anual sobre la geometría de paneles |

El nombre de cada columna sigue el formato `SA<subarray>_ST<string>` (ej. `SA1_ST3`).

---

### `shade_postprocess.py` — Post-procesado de resultados

Script Python que consume los CSVs generados por el LK y produce visualizaciones de análisis. Se puede ejecutar manualmente o es llamado automáticamente desde el LK si `py_script` apunta a él.

```bash
python3 shade_postprocess.py \
  --ts     shade_timeseries_az180_inc30.csv \
  --out    ./resultados \
  --width  1.0 \
  --length 2.0
```

> El script se ejecuta automáticamente al final de `SAM_3D_Shading.lk`. También puede lanzarse manualmente con el comando anterior.

#### Salidas

| Archivo | Descripción |
|---|---|
| `summary_statistics.csv` | Media, mediana, P90, máximo y horas con sombra >10 % y >50 % por grupo |
| `seasonal_hourly_curves.png` | 4 subplots (Invierno/Primavera/Verano/Otoño): sombra media horaria por grupo + media del array; eje X en formato HH:MM |
| `heatmap_panel_geometry.png` | Grid Subarray × String con valor de sombra (%) por estación y anual; proporciones reales del panel |

---

## Requisitos

- SAM 2022.11.21 o superior (testado en SAM 2024.12.12).
- El script debe ejecutarse desde el **scripting del 3D Shade Calculator**, no desde el scripting general de SAM.

---

## Roadmap

- [ ] Script de análisis paramétrico (ubicación × inclinación)
- [ ] Script de generación de obstáculos de sombra desde CSV

---

## Licencia

MIT
