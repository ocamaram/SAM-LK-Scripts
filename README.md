# SAM LK Scripts

Librería de scripts LK para [NREL System Advisor Model (SAM)](https://sam.nrel.gov/), orientada al análisis de sombreado 3D y simulación de sistemas fotovoltaicos.

---

## Instalación de dependencias Python

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

---

## Scripts disponibles

### `SAM_3D_Shading.lk` — Generador de superficies activas y análisis de sombras

Genera una cuadrícula de superficies activas en el **3D Shade Calculator** de SAM, lanza el análisis de series temporales y exporta los porcentajes de sombra por subarray y string a CSV.

#### Características

- Crea `rows × cols` superficies activas con orientación (azimut + inclinación) configurable.
- Asigna automáticamente cada superficie a su subarray y string correspondientes.
- **Sistema de batches automático**: SAM limita el análisis a 32 superficies (4 subarrays × 8 strings). El script divide el array en tantos batches como sean necesarios, los analiza secuencialmente y fusiona los resultados en Python.
- **Preserva la geometría de obstáculos existente** (edificios, árboles, cajas…) usando un archivo de escena temporal; solo las superficies activas se recrean en cada batch.
- Exporta los resultados a CSV y **lanza automáticamente** `shade_postprocess.py` al finalizar.

#### Parámetros

| Parámetro | Descripción | Valor por defecto |
|---|---|---|
| `rows` | Número de filas de superficies | `10` |
| `cols` | Número de columnas de superficies | `8` |
| `x0`, `y0`, `z0` | Origen de la cuadrícula (m) | `0.0` |
| `panel_length` | Largo del panel — dirección de filas (m) | `2.0` |
| `panel_width` | Ancho del panel — dirección de columnas (m) | `1.0` |
| `spacing_x` | Separación entre filas (m, ≥ `panel_length`) | `panel_length` |
| `spacing_y` | Separación entre columnas (m, ≥ `panel_width`) | `panel_width` |
| `azimut` | Orientación solar — 0=N, 90=E, 180=S, 270=O (°) | `180` |
| `inclinacion` | Inclinación — 0=horizontal, 90=vertical (°) | `30` |
| `timestep_min` | Paso temporal del análisis (min) | `60` |
| `irradiance_ref` | Irradiancia de referencia para heatmaps (W/m²), usado si no se dispone de serie temporal | `1000.0` |

> **Límites SAM:** máximo 32 superficies activas por análisis (4 subarrays × 8 strings). Arrays más grandes se procesan en múltiples batches automáticamente.

#### Uso

1. Abre el **3D Shade Calculator** desde la pestaña *Shading* de tu caso SAM.
2. Configura la **ubicación** en la pestaña *Location* (necesaria para el cálculo solar).
3. Añade opcionalmente **obstáculos** de sombra (edificios, árboles…) en la pestaña *3D Scene* con la escena limpia (sin superficies activas).
4. Abre la pestaña *Scripting*, carga `SAM_3D_Shading.lk` y ajusta los parámetros de la sección `PARÁMETROS`.
5. Ejecuta el script. Los resultados se guardan en la subcarpeta configurada en `results_dir`.

> **Protocolo para modificar obstáculos:** pulsa *New* en SAM para limpiar la escena, añade los obstáculos nuevos y ejecuta el script. En ejecuciones normales (sin cambiar obstáculos) puedes ejecutar directamente.

#### Archivos de salida

| Archivo | Contenido |
|---|---|
| `shade_batch<N>_az<A>_inc<I>.csv` | Series temporales del factor de sombra directa (0–100) por batch |
| `summary_statistics.csv` | Resumen estadístico por panel (generado por `shade_postprocess.py`) |
| `seasonal_hourly_curves.png` | Curvas horarias promedio de sombra por estación |
| `heatmap_panel_geometry.png` | Heatmap de sombra (%) estacional y anual |
| `heatmap_<estacion>.png` | Heatmap individual por estación y anual (alta resolución) |
| `heatmap_irradiance.png` | Heatmap de irradiancia recibida (W/m²) estacional y anual |
| `heatmap_irradiance_<estacion>.png` | Heatmap de irradiancia individual por estación y anual |

---

### `shade_postprocess.py` — Post-procesado de resultados

Script Python que consume los CSVs generados por el LK y produce estadísticas y visualizaciones. Se ejecuta automáticamente al final del LK, pero también puede lanzarse manualmente.

#### Modo single-batch (array ≤ 32 paneles)

```bash
python3 shade_postprocess.py \
  --ts      Results/shade_batch0_az180_inc30.csv \
  --out     Results/ \
  --width   1.0 \
  --length  1.0
```

#### Modo multi-batch con irradiancia horaria real

```bash
python3 shade_postprocess.py \
  --batches Results/shade_batch0_az225_inc0.csv \
             Results/shade_batch1_az225_inc0.csv \
  --out      Results/ \
  --width    1.0 \
  --length   1.0 \
  --cols     205 \
  --irradiance-ts  POA_total_radiation_nominal.csv
```

Si SAM exporta la irradiancia en kW totales del array (no en W/m²), añade `--irradiance-area` con el área total en m²:

```bash
  --irradiance-ts  POA_total_radiation_nominal.csv \
  --irradiance-area 1254
```

#### Argumentos

| Argumento | Descripción | Por defecto |
|---|---|---|
| `--ts` | CSV de un único batch (modo single) | — |
| `--batches` | Uno o más CSVs de batches a combinar | — |
| `--out` | Carpeta de salida | requerido |
| `--width` | Ancho del panel (m) | `1.0` |
| `--length` | Largo del panel (m) | `1.0` |
| `--spacing-x` | Separación entre filas (m). Si difiere de `--length`, corrige las distancias en los ejes | `= --length` |
| `--spacing-y` | Separación entre columnas (m). Si difiere de `--width`, corrige las distancias en los ejes | `= --width` |
| `--cols` | Columnas del array completo (necesario con `--batches`) | `8` |
| `--irradiance-ts` | CSV con irradiancia horaria exportada desde SAM (W/m² o kW). Activa ponderación temporal real | — |
| `--irradiance-col` | Nombre de la columna en `--irradiance-ts` (auto-detecta si se omite) | primera columna numérica |
| `--irradiance-area` | Área total del array (m²). Convierte los valores de `--irradiance-ts` de kW a W/m² | — |
| `--irradiance` | Irradiancia constante de referencia (W/m²), solo si no se usa `--irradiance-ts` | `1000.0` |
| `--dpi` | Resolución de las imágenes de salida | `150` |

#### Cómo obtener la irradiancia horaria desde SAM

En los resultados de simulación de SAM, busca la variable **"Array POA front-side total radiation nominal (kW)"** y expórtala a CSV. Esta variable contiene la irradiancia total incidente sobre el array antes de aplicar sombras. Si el CSV contiene valores en kW (no W/m²), usa `--irradiance-area` con el área total del array.

> La ponderación temporal es importante porque una sombra del 50% al mediodía (irradiancia ~800 W/m²) implica cuatro veces más pérdida de energía que la misma sombra a primera hora de la mañana (~200 W/m²). La media simple de fracción de sombra subestima sistemáticamente el impacto real.

#### Salidas

| Archivo | Descripción |
|---|---|
| `summary_statistics.csv` | Estadísticas por panel: media, mediana, P90, máximo, horas con sombra >10%/>50%. Si se usa `--irradiance-ts`: añade `sombra_irr_ponderada_%`, `energia_bloqueada_kWh_m2` y `energia_incidente_kWh_m2` |
| `seasonal_hourly_curves.png` | 4 subplots estacionales: sombra media horaria por panel + media del array |
| `heatmap_panel_geometry.png` | Grid con sombra (%) por estación y anual; ejes en metros |
| `heatmap_<estacion>.png` | Heatmap individual de sombra por estación y anual |
| `heatmap_irradiance.png` | Grid con irradiancia media recibida (W/m²) por estación y anual |
| `heatmap_irradiance_<estacion>.png` | Heatmap individual de irradiancia por estación y anual |

Los ejes de todos los heatmaps muestran distancias en metros (x, y). Si la separación entre paneles difiere del tamaño del panel, especifica `--spacing-x` / `--spacing-y` para obtener coordenadas físicas correctas.

---

## Requisitos

- SAM 2022.11.21 o superior (testado en SAM 2024.12.12).
- El script `.lk` debe ejecutarse desde el **scripting del 3D Shade Calculator**, no desde el scripting general de SAM.

---

## Roadmap

- [ ] Script de análisis paramétrico (ubicación × inclinación)
- [ ] Script de generación de obstáculos de sombra desde CSV

---

## Licencia

MIT
