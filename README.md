Perfecto, Jessica. Aquí tienes el contenido completo, revisado y listo para pegar directamente en la página de GitHub en formato Markdown (`README.md`). Incluye todos los apartados esenciales: objetivo, estructura del repositorio, ejecución, entorno, scripts, autoría, licencia y cita sugerida.

---

```markdown
# TFM_MDD_framework

Framework multiómico para la identificación de biomarcadores en depresión mayor (MDD), centrado en transcriptómica, con integración de análisis adaptativo y automatizado para datos de microarrays y RNA-seq.

---

## 🧠 Objetivo del proyecto

Desarrollar un framework reproducible y modular que permita el preprocesamiento, análisis diferencial, enriquecimiento funcional e interpretación de datos ómicos relacionados con la depresión, con foco inicial en datos transcriptómicos (Affymetrix y RNA-seq).

---

## 📁 Estructura del repositorio

```

TFM\_MDD/
├── data/                  # Datos brutos y metadatos (no incluidos en Git)
│   └── transcriptomics/
│
├── notebooks/             # Exploraciones interactivas (Jupyter)
│   └── README.md
│
├── pipeline/              # Automatización (Snakemake u otros)
│   └── README.md
│
├── results/               # Resultados generados (no incluidos)
│   ├── GSE98793/
│   ├── GSE44593/
│   └── log\_preprocessing.txt
│
├── scripts/               # Scripts organizados por tipo
│   ├── general/           # Descarga, mapeo, QC y master
│   └── transcriptomics/   # Expresión diferencial, enriquecimiento y validación
│
├── environment.yml        # Entorno reproducible Conda
└── .gitignore             # Exclusión de carpetas y archivos pesados

````

---

## ⚙️ Instalación del entorno

Requiere Conda:

```bash
conda env create -f environment.yml
conda activate tfm_mdd
````

---

## 🚀 Cómo ejecutar los análisis

### Preprocesamiento adaptativo

```bash
Rscript scripts/general/preprocessing_master.R
```

### Análisis de expresión diferencial

```bash
Rscript scripts/transcriptomics/differential_expression.R
```

### Enriquecimiento funcional (GO/KEGG)

```bash
Rscript scripts/transcriptomics/enrichment_analysis.R
```

> También disponibles versiones interactivas para ejecución manual en RStudio:
>
> * `run_differential_expression_interactive.R`
> * `run_enrichment_interactive.R`

---

## 🧬 Scripts incluidos

**General:**

* `prepare_geo_dataset.sh`: descarga automatizada desde GEO
* `download_raw_geo.R`, `download_metadata_geo.R`: extracción estructurada
* `map_gene_ids.R`: mapeo de identificadores a Entrez/SYMBOL
* `preprocessing_master.R`: orquestador general por dataset
* `qc_post_normalization.R`: validación técnica post-normalización

**Transcriptómicos:**

* `preprocessing_microarray.R`, `preprocessing_rnaseq.R`: flujo adaptativo según plataforma
* `differential_expression.R`: análisis DESeq2 o limma según tipo de datos
* `enrichment_analysis.R`: enriquecimiento GO/KEGG
* `validate_logfc_direction.R`: verificación de sentido de regulación

---

## 📄 Datos

Los datos transcriptómicos (microarrays y RNA-seq) y los resultados generados no se incluyen en el repositorio por motivos de tamaño y reproducibilidad. Su estructura está descrita en:

* `data/README.md`: organización por tipo y GSE ID
* `results/README.md`: estructura de salida por análisis

---

## 👩‍💻 Autoría

Este proyecto forma parte del Trabajo de Fin de Máster de:

**Jessica Valdivia**
Máster en Bioinformática y Bioestadística
Universitat Oberta de Catalunya (UOC)

---

## 📄 Licencia

MIT License — Puedes usar, adaptar y distribuir este framework con fines académicos o personales.

---

## 📚 Cita sugerida

Si utilizas o te basas en este framework, por favor cita:

```
Valdivia, J. (2025). TFM_MDD_framework: Framework multiómico para la identificación de biomarcadores en depresión mayor. Trabajo de Fin de Máster, UOC.
```

```

---
