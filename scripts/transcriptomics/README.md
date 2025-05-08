# 📂 scripts/transcriptomics

Este directorio contiene los scripts específicos para el análisis transcriptómico dentro del framework multiómico del TFM. Cubre desde el preprocesamiento de datos crudos hasta el análisis funcional e interpretación fisiopatológica, con validación cruzada de biomarcadores. Todos los scripts son modulares y compatibles con ejecución directa o automatizada (Snakemake).

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

