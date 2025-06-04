# 🧬 TFM MDD Framework - Snakemake Workflow

Este repositorio contiene el flujo completo de análisis multiómico para el TFM sobre biomarcadores de la depresión (MDD), incluyendo:

- Transcriptómica (microarrays)
- Neuroimagen funcional (MEG)
- Validación cruzada inter-ómica

## 📁 Estructura del proyecto
```
TFM_MDD/
├── workflow/
│   ├── Snakefile               <- Snakefile principal con reglas incluidas
│   ├── config/config.yaml      <- Configuración del pipeline (datasets, rutas, entornos)
│   ├── rules/                  <- Reglas separadas por módulo
│   │   ├── transcriptomics.smk
│   │   ├── neuroimaging.smk
│   │   ├── interomics.smk
│   │   └── common.smk
│   └── envs/                   <- Entornos Conda por módulo
│       ├── tfm_mdd_framework.yml
│       └── tfm_meg_env.yml
```

## ⚙️ Requisitos
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) o [Anaconda](https://www.anaconda.com/)
- [Snakemake](https://snakemake.readthedocs.io/) (instalable con `conda install -c bioconda snakemake`)
- [Graphviz](https://graphviz.org/download/) (opcional, para visualizar el flujo)

## 📦 Instalación de entornos Conda
Desde el directorio `workflow/`:

```bash
conda env create -f envs/tfm_mdd_framework.yml
conda env create -f envs/tfm_meg_env.yml
```

## 🚀 Ejecución del pipeline completo
Desde `~/TFM_MDD/workflow/`:

```bash
snakemake -j 4 --use-conda --rerun-incomplete --printshellcmds
```

## 🧪 Ejecutar una regla concreta (ejemplo QC MEG)
```bash
snakemake qc_meg_raw -j 1 --use-conda
```

## 🗂️ Modificar datasets y parámetros
Edita el archivo:
```
workflow/config/config.yaml
```
Para cambiar los IDs de datasets, rutas base o parámetros.

## 📤 Resultados esperados
- `results/<dataset>/normalized/` → expresión normalizada
- `results/<dataset>/deg/` → genes diferenciales
- `results/<dataset>/enrichment/` → análisis GO/KEGG
- `results/ds005356/features/` → features MEG
- `results/interomics/` → informe PDF final con biomarcadores multiómicos

---

## ❓ Soporte y contacto
Este workflow está diseñado para ser usado por personal técnico o clínico. No requiere conocimientos de programación, pero sí familiaridad básica con archivos y terminal Linux.

Para dudas contactar con la autora del TFM o revisar los scripts en `scripts/`.

---

© 2025. Proyecto académico. No redistribuir sin permiso.
