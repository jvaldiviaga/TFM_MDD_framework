# 📂 scripts/general

Este directorio contiene scripts de uso general dentro del framework multiómico del TFM. Son independientes del tipo de ómica (transcriptómica, metilación, genómica, neuroimagen) y automatizan tareas clave como descarga de datos, detección del tipo de ómica, preprocesamiento, evaluación de calidad y mapeo de identificadores para análisis funcional.

Todos los scripts son reutilizables, modulares y compatibles con ejecución directa por consola o integración en flujos automatizados tipo Snakemake.

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
