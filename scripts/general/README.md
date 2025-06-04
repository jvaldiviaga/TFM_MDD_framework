# 📂 scripts/general

Este directorio contiene scripts de uso general dentro del framework multiómico del TFM. Son independientes del tipo de ómica (transcriptómica, metilación, genómica, neuroimagen) y automatizan tareas clave como descarga de datos, detección del tipo de ómica, preprocesamiento, evaluación de calidad y mapeo de identificadores para análisis funcional.

Todos los scripts son reutilizables, modulares y compatibles con ejecución directa por consola o integración en flujos automatizados tipo Snakemake.

# scripts/general/README.md

Este directorio contiene scripts generales y reutilizables para tareas de preprocesamiento, control de calidad, validación cruzada, mapeo de IDs y generación de informes en el contexto del framework multiómico para el estudio del MDD.

| Script                              | Descripción                                                                   | Inputs                                              | Outputs                                                   | Ejemplo de uso                                                                                                                                                         |
| ----------------------------------- | ----------------------------------------------------------------------------- | --------------------------------------------------- | --------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `prepare_geo_dataset.sh`            | Prepara estructura de carpetas para un GSE dado                               | ID del dataset                                      | Carpetas `data/`, `results/` con subcarpetas vacías       | `bash scripts/general/prepare_geo_dataset.sh GSE98793`                                                                                                                 |
| `download_metadata_geo.R`           | Descarga metadatos de GEO y guarda `series_matrix.txt.gz`                     | ID del dataset GEO                                  | `*_metadata.csv`, `*_series_matrix.txt.gz`                | `Rscript scripts/general/download_metadata_geo.R GSE98793`                                                                                                             |
| `download_raw_geo.R`                | Descarga datos crudos (.CEL) de GEO automáticamente                           | ID del dataset                                      | Carpeta con archivos `.CEL`                               | `Rscript scripts/general/download_raw_geo.R GSE98793`                                                                                                                  |
| `preprocessing_master.R`            | Detecta el tipo de datos (CEL, matrix) y lanza el preprocesamiento adaptativo | Ruta al dataset (`data/transcriptomics/GSE[ID]`)    | Archivos preprocesados, normalizados y QC previo          | `Rscript scripts/general/preprocessing_master.R data/transcriptomics/GSE98793`                                                                                         |
| `qc_post_normalization.R`           | Evalúa calidad tras normalización y genera informe interpretativo             | ID del dataset                                      | `qc_interpretacion.txt`, visualizaciones                  | `Rscript scripts/general/qc_post_normalization.R GSE98793`                                                                                                             |
| `map_gene_ids.R`                    | Mapea probes a ENTREZID y SYMBOL (basado en plataforma)                       | ID del dataset                                      | Archivos `*_mapped.csv`                                   | `Rscript scripts/general/map_gene_ids.R GSE98793`                                                                                                                      |
| `check_ready_for_DEG.R`             | Verifica formato, log2, número de muestras y metadatos antes de DEG           | `*_expression.csv`, `*_metadata.csv`                | `*_qc_ready_for_DEG.txt`                                  | `Rscript scripts/general/check_ready_for_DEG.R data/.../expression.csv data/.../metadata.csv`                                                                          |
| `check_metadata_columns.R`          | Revisa automáticamente los nombres de columnas clínicas y su formato          | Archivo de metadatos                                | Reporte por consola                                       | `Rscript scripts/general/check_metadata_columns.R metadata.csv`                                                                                                        |
| `filter_metadata_and_expression.R`  | Filtra metadatos y matriz de expresión por grupos                             | Metadata, grupos A/B, archivo expresión             | Archivos filtrados                                        | `Rscript filter_metadata_and_expression.R metadata.csv metadata_filtered.csv --group_column=X.disease --group_A="MDD" --group_B="CTL" --expression_csv=expression.csv` |
| `filter_metadata_by_groups.R`       | Versión simplificada sin filtrado de expresión                                | Metadata, grupos A/B                                | `*_metadata_filtered.csv`                                 | `Rscript filter_metadata_by_groups.R ...`                                                                                                                              |
| `robust_biomarker_selection.R`      | Aplica LASSO (50 seeds) y RF (20 runs) para selección robusta                 | ID del dataset                                      | `lasso_selected_features.csv`, `rf_selected_features.csv` | `Rscript robust_biomarker_selection.R GSE98793`                                                                                                                        |
| `combine_panels_expandido.R`        | Une paneles robustos de múltiples datasets                                    | Archivos `panel_robusto_expandido.csv`              | `combined_panel.csv`                                      | `Rscript combine_panels_expandido.R`                                                                                                                                   |
| `validate_interomics.R`             | Compara panel transcriptómico y MEG por mecanismos fisiopatológicos           | Archivos de genes, features MEG, rutas GO           | Resultados de validación cruzada                          | `Rscript validate_interomics.R`                                                                                                                                        |
| `generate_biomarker_report.R`       | Genera informe clínico e interpretativo de los biomarcadores (genes y MEG)    | Archivos de salida de DEG, enrichment y regresiones | `biomarker_report.pdf`, `summary.csv`                     | `Rscript generate_biomarker_report.R`                                                                                                                                  |
| `interpret_biomarkers_multiomics.R` | Mapea rutas GO a mecanismos del MDD                                           | `go_enrichment_results.csv`                         | Tabla de mecanismos por ruta/gen                          | `Rscript interpret_biomarkers_multiomics.R`                                                                                                                            |

**Nota**: Todos los scripts son modulares y pueden integrarse en Snakemake. Requieren el entorno `tfm_mdd_framework` activo, salvo indicación contraria.

---

Para dudas sobre la configuración del entorno y dependencias, ver `environment.yml` o contactar con el autor del proyecto.

---

## 🌐 FASE 1 — Descarga y preparación de datos desde GEO

### `prepare_geo_dataset.sh`

Script en Bash para automatizar la descarga de datasets desde GEO:

- Detecta el tipo de ómica según el contenido del `.tar`
- Descarga metadatos y datos crudos
- Llama automáticamente a `download_raw_geo.R` si la descarga directa falla
- Registra el tipo de ómica en `dataset_type.txt`
- Organiza archivos en carpetas limpias (`final/`, `series_matrix.txt.gz`)

**Uso:**
```bash
bash scripts/general/prepare_geo_dataset.sh GSEXXXXX
```

---

### `download_metadata_geo.R`

Descarga el archivo `series_matrix.txt.gz` y extrae automáticamente los metadatos clínicos en formato `.csv` si están presentes en un `ExpressionSet`.

**Input:**

- ID de dataset GEO (ej. GSE98793)

**Output:**

- `GSEXXXXX_series_matrix.txt.gz`
- `GSEXXXXX_metadata.csv`

**Uso:**
```bash
Rscript scripts/general/download_metadata_geo.R GSEXXXXX
```

---

### `download_raw_geo.R`

Script auxiliar que descarga los archivos suplementarios (`.tar`) desde GEO usando `getGEOSuppFiles()` de `GEOquery`.

**Input:**

- ID de dataset GEO

**Output:**

- Archivos crudos del estudio

**Uso:**
```bash
Rscript scripts/general/download_raw_geo.R GSEXXXXX
```

---

## ⚙️ FASE 2 — Preprocesamiento adaptativo

### `preprocessing_master.R`

Script maestro que detecta automáticamente el tipo de datos en una carpeta GEO (microarrays, RNA-seq, metilación, imagen, etc.) y lanza el script de preprocesamiento correspondiente.

- Detecta tipo de ómica por extensión o `dataset_type.txt`
- Llama al script apropiado (`preprocessing_microarray.R`, etc.)
- Evita reprocesar si ya existe `.rds`
- Registra logs de ejecución

**Input:**

- Carpeta de datos `GSEXXXXX/` con archivos crudos + metadatos

**Output:**

- Matriz de expresión normalizada `.rds`
- Gráficos y logs generados por el script correspondiente

**Uso:**
```bash
Rscript scripts/general/preprocessing_master.R data/transcriptomics/GSEXXXXX/
```

---

## 🧠 FASE 3 — Evaluación de calidad post-normalización

### `qc_post_normalization.R`

Evalúa la calidad técnica de los datos normalizados:

- Genera boxplot, SD por muestra, carga total
- PCA coloreado por variables clínicas detectadas automáticamente
- Heatmap anotado por metadatos
- Detección de variable con mayor separación en PC1
- Estimación de clústeres jerárquicos
- Ejecución de `arrayQualityMetrics` si n ≤ 100
- Generación de informes automáticos interpretables

**Input:**

- `normalized_expression.rds`
- `*_metadata.csv` o `*_metadata_clean.csv` (opcional pero recomendado)

**Output:**

- `boxplot_normalized.png`, `heatmap_normalized.png`, `pca_<variable>.png`
- `qc_interpretacion.txt`, `interpretacion_automatizada.txt`
- `AQM_Report/` (si aplica)

**Uso:**
```bash
Rscript scripts/general/qc_post_normalization.R GSEXXXXX
```

---

## 🔍 FASE 4 — Revisión de metadatos

### `check_metadata_columns.R`

Revisa las columnas del archivo `*_metadata_clean.csv` y ayuda a identificar variables clínicas válidas.

- Interpreta nombres (diagnóstico, ansiedad, sexo, edad, etc.)
- Muestra número de niveles y ejemplos de valores
- Útil antes de lanzar el análisis de expresión diferencial

**Input:**

- `*_metadata_clean.csv`

**Uso:**
```bash
Rscript scripts/general/check_metadata_columns.R GSEXXXXX
```

---

## 🧬 FASE 5 — Mapeo de identificadores

### `map_gene_ids.R`

Mapea los genes diferenciales (`DEG_results.csv`, y variantes UP/DOWN) a `SYMBOL` y `ENTREZID`, necesarios para enriquecimiento funcional.

- Itera sobre múltiples configuraciones (`no_covariates`, `with_batch`, etc.)
- Usa `org.Hs.eg.db` y, si falla, tabla de anotación GPL de fallback
- Aplica el mapeo a DEG global y a UP/DOWN

**Input:**

- `DEG_results.csv`, `DEG_upregulated.csv`, `DEG_downregulated.csv`

**Output:**

- `DEG_results_mapped.csv`, `*_mapped.csv`
- `mapping_summary.txt`

**Uso:**
```bash
Rscript scripts/general/map_gene_ids.R GSEXXXXX
```

---

## 📌 Notas técnicas

- Todos los scripts están preparados para ejecución directa (`Rscript script.R`) o integración en Snakemake.
- Son independientes del tipo de ómica, y reutilizables en flujos de transcriptómica, metilación, GWAS o neuroimagen.
- Las rutas siguen la estructura modular del framework: `~/TFM_MDD/data/` y `~/TFM_MDD/results/`.
