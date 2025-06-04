# 📂 scripts/transcriptomics

Este directorio contiene los scripts específicos para el análisis transcriptómico dentro del framework multiómico del TFM. Cubre desde el preprocesamiento de datos crudos hasta el análisis funcional e interpretación fisiopatológica, con validación cruzada de biomarcadores. Todos los scripts son modulares y compatibles con ejecución directa o automatizada (Snakemake).

# scripts/transcriptomics/README.md

Este directorio contiene los scripts específicos para el análisis de datos transcriptómicos en el marco del estudio de biomarcadores para la depresión mayor (MDD). Incluye preprocesamiento, análisis de expresión diferencial, enriquecimiento funcional y validación cruzada entre datasets.

| Script                                      | Descripción                                                     | Inputs                               | Outputs                                        | Ejemplo de uso                                                               |
| ------------------------------------------- | --------------------------------------------------------------- | ------------------------------------ | ---------------------------------------------- | ---------------------------------------------------------------------------- |
| `preprocessing_microarray.R`                | Aplica preprocesamiento técnico adaptativo a archivos .CEL      | Archivos `.CEL` de Affymetrix        | Matriz normalizada, QC previo, informe técnico | `Rscript preprocessing_microarray.R data/transcriptomics/GSE98793/`          |
| `preprocessing_rnaseq.R`                    | Preprocesamiento adaptativo para RNA-seq ya cuantificado        | `expression.csv`, `metadata.csv`     | Datos filtrados y log-transformados            | `Rscript preprocessing_rnaseq.R expression.csv metadata.csv`                 |
| `process_series_matrix.R`                   | Extrae matriz y metadatos desde `series_matrix.txt.gz`          | Archivo `.txt.gz` descargado de GEO  | `*_expression.csv`, `*_metadata.csv`           | `Rscript process_series_matrix.R GSE39653`                                   |
| `differential_expression.R`                 | Lógica principal para DEG con y sin covariables (limma, DESeq2) | `expression.rds`, `metadata.csv`     | Resultados DEG con/ sin ajuste                 | `Rscript differential_expression.R GSE98793 diagnosis control age,gender`    |
| `run_differential_expression_interactive.R` | Script interactivo para ejecución manual del análisis DEG       | Ninguno (usa archivos predefinidos)  | Resultados en subcarpetas                      | `source("run_differential_expression_interactive.R")` en R                   |
| `enrichment_analysis.R`                     | Análisis funcional por GO, KEGG y Reactome, up/down separados   | Archivos `*_mapped.csv`              | `.csv` de rutas enriquecidas, gráficos         | `Rscript enrichment_analysis.R GSE98793`                                     |
| `run_enrichment_interactive.R`              | Versión interactiva para análisis funcional                     | Archivos DEG ya mapeados             | `.txt` interpretativo, dotplots                | `source("run_enrichment_interactive.R")` en R                                |
| `validate_logfc_direction.R`                | Compara dirección del efecto (logFC) entre datasets             | Base y dataset de validación         | `genes_validados.csv`, resumen                 | `Rscript validate_logfc_direction.R GSE98793 with_covars GSE39653 no_covars` |
| `validate_panel_overlap.R`                  | Verifica si genes validados aparecen en el panel principal      | Panel principal y genes validados    | `overlap_summary.txt`                          | `Rscript validate_panel_overlap.R`                                           |
| `validate_inter_dataset_enrichment.R`       | Consolida evidencia multigen entre datasets                     | Resultados de validación por dataset | `inter_dataset_robustness_summary.csv`         | `Rscript validate_inter_dataset_enrichment.R`                                |

**Subcarpeta `processed_data/`**: contiene scripts utilitarios para conversiones de formato u otras tareas técnicas específicas como `convert_expression_csv_to_rds.R`.

---

Para ejecutar correctamente estos scripts es necesario activar el entorno `tfm_mdd_framework` previamente:

```bash
conda activate tfm_mdd_framework
```

Todos los scripts han sido diseñados para ser reutilizables y compatibles con el flujo de Snakemake (`workflow/`).

---

## 🔁 FASE 1 — Preprocesamiento

### `preprocessing_microarray.R`

Preprocesa datos `.CEL` (Affymetrix) con control de calidad adaptativo:

- Genera visualizaciones de QC (boxplot, PCA, heatmap)
- Evalúa métricas técnicas (SD, carga, efecto batch)
- Decide entre normalización `rma()` o `vsn()`
- Aplica `ComBat()` si se detecta batch
- Exporta resultados y justificación técnica

**Input:**  
- Archivos `.CEL` en `data/transcriptomics/GSEXXXXX/`  
- `GSEXXXXX_metadata.csv`

**Output:**  
- `normalized_expression.rds`  
- Gráficos: `boxplot_normalized.png`, `heatmap_normalized.png`  
- Informes: `decision_normalizacion.txt`, `qc_interpretacion.txt`

---

## 📊 FASE 2 — Análisis de expresión diferencial

### `differential_expression.R`

Ejecuta análisis diferencial con `limma`:

- Comparación entre grupos clínicos definidos por el usuario
- Ajuste por covariables (sexo, edad, batch, etc.)
- Selección automática del coeficiente relevante
- Exporta resultados completos y subconjuntos UP/DOWN

**Input:**  
- `normalized_expression.rds`  
- `*_metadata_clean.csv`

**Output:**  
- `DEG_results.csv`, `DEG_upregulated.csv`, `DEG_downregulated.csv`  
- Gráficos: `volcano_plot.png`, `heatmap_top50.png`, `pca_top100.png`  
- Informes: `summary_DEG.txt`, `log_comando.txt`

---

### `run_differential_expression_interactive.R`

Asistente interactivo por consola:

- Lista todas las variables disponibles en el metadata
- Permite elegir variable a comparar, grupo de referencia y covariables
- Llama automáticamente a `differential_expression.R`

---

## 🧬 FASE 3 — Análisis de enriquecimiento funcional

### `enrichment_analysis.R`

Realiza enriquecimiento funcional sobre genes significativos (`adj.P.Val < 0.05`):

- ORA sobre GO (BP, MF, CC), KEGG y Reactome
- GSEA si hay estadístico `t`
- Análisis por subconjunto (ALL / UP / DOWN)
- Interpretación fisiopatológica automática (MDD)
- Análisis GO/KEGG por gen, glosario fisiopatológico, salida textual

**Input:**  
- `DEG_results_mapped.csv` generado tras `map_gene_ids.R`

**Output:**  
- Enriquecimientos `.csv`, `.rds` por base y grupo  
- Visualizaciones: `*_dotplot.png`, `*_emapplot.png`, `*_cnetplot.png`  
- Interpretaciones:  
  - `interpretacion_fisiopatologica.txt`  
  - `interpretacion_fisiopatologica_UP.txt`  
  - `ruta_por_gen.csv`, `interpretacion_por_gen.txt`  
  - `summary_enrichment.txt`

---

### `run_enrichment_interactive.R`

Selector interactivo que permite:

- Elegir una subcarpeta de resultados DEG (con/sin covariables)
- Lanzar automáticamente `enrichment_analysis.R`

---

## 🔁 FASE 4 — Validación cruzada de biomarcadores

### `validate_logfc_direction.R`

Compara la dirección del cambio (signo de logFC) entre dos análisis:

- Indica si el gen mantiene la dirección (`✅`) o cambia (`❌`)
- Calcula diferencia absoluta de logFC
- Genera tabla comparativa y visualización

**Input:**  
- Dos archivos `DEG_results_mapped.csv`

**Output:**  
- CSV de validación cruzada  
- Gráfico `logFC_comparison_*.png`  
- Estadísticas de concordancia por consola

---

