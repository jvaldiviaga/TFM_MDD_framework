# 📂 pipeline

Esta carpeta está destinada a contener archivos de automatización como flujos de trabajo Snakemake, Makefiles o scripts maestros de ejecución por lotes.

Permite orquestar todo el procesamiento y análisis desde un punto centralizado y reproducible.

# FLUJO GENERAL DEL PIPELINE (`pipeline/README.md`)

Este directorio documenta el flujo completo del framework automatizado desarrollado para el análisis multiómico de la depresión mayor (MDD). Se incluyen tanto los pasos para transcriptómica como neuroimagen MEG, organizados por fases modulares.

---

## 🔬 1. TRANSCRIPTÓMICA — FLUJO GENERAL

### Paso 1: Descarga y preparación del dataset

```bash
bash scripts/general/prepare_geo_dataset.sh GSE[ID]  # estructura carpetas
Rscript scripts/general/download_metadata_geo.R GSE[ID]  # metadatos
```

### Paso 2: Preprocesamiento adaptativo

```bash
Rscript scripts/general/preprocessing_master.R data/transcriptomics/GSE[ID]
```

### Paso 3: Evaluación post-normalización

```bash
Rscript scripts/general/qc_post_normalization.R GSE[ID]
```

### Paso 4: Análisis de expresión diferencial

```r
setwd("scripts/transcriptomics/")
source("run_differential_expression_interactive.R")
```

### Paso 5: Mapeo de IDs

```bash
Rscript scripts/general/map_gene_ids.R GSE[ID]
```

### Paso 6: Enriquecimiento funcional

```bash
conda activate tfm_mdd_framework
R
source("run_enrichment_interactive.R")
```

### Paso 7: Validación cruzada entre datasets

```bash
Rscript scripts/transcriptomics/validate_logfc_direction.R BASE_ID BASE_CONFIG VAL_ID VAL_CONFIG
```

---

## 🧠 2. NEUROIMAGEN MEG — FLUJO GENERAL

### Paso 1: Fusión de archivos por sujeto

```bash
conda activate tfm_meg_env
python scripts/neuroimaging/load_all_meg_subjects.py --data_root data/neuroimaging/ds005356
```

### Paso 2: QC visual previo

```bash
python scripts/neuroimaging/run_qc_all_meg.py \
  --input_dir data/neuroimaging/ds005356/derivatives/MEGcombined \
  --output_dir results/ds005356/qc_pre_filtering
```

### Paso 3: Evaluación efecto batch

```bash
python scripts/neuroimaging/analyze_batch_meg_raw.py \
  --data_dir data/neuroimaging/ds005356/derivatives/MEGcombined \
  --metadata_file data/neuroimaging/Code/MEG\ MDD\ IDs\ and\ Quex.xlsx \
  --output_dir results/ds005356/analysis_batch
```

### Paso 4: Preprocesamiento completo

```bash
python scripts/neuroimaging/run_all_preprocessing.py
```

### Paso 5: QC post-preprocesamiento

```bash
python scripts/neuroimaging/run_qc_post_all.py
```

### Paso 6: Epoching y calidad de epochs

```bash
python scripts/neuroimaging/run_all_epoching.py
python scripts/neuroimaging/run_all_epoch_quality.py
```

### Paso 7: Filtros adicionales y extracción de features

```bash
python scripts/neuroimaging/run_all_filtering.py
bash scripts/neuroimaging/run_all_extract_features.sh
```

### Paso 8: Unión con datos clínicos y análisis grupal

```bash
bash scripts/neuroimaging/run_merge_with_clinical.sh
bash scripts/neuroimaging/run_group_stats_meg.sh
```

---

## 🔁 3. VALIDACIÓN Y CONSOLIDACIÓN

### Paso 1: Selección robusta de biomarcadores por dataset

```bash
Rscript scripts/general/robust_biomarker_selection.R GSE[ID]
```

### Paso 2: Combinación de paneles y comparación interdataset

```bash
Rscript scripts/general/combine_panels_expandido.R
```

### Paso 3: Validación cruzada multiómica y reporte clínico

```bash
Rscript scripts/general/validate_interomics.R
Rscript scripts/general/generate_biomarker_report.R
```

---

Este flujo puede ejecutarse de forma modular o integrarse en un `workflow/Snakefile` automatizado. Para más detalles, ver los README individuales de cada subcarpeta.


