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

# ACTUALIZACIÓN!!!

# workflow/README.md

Este directorio contiene el archivo `Snakefile` y las reglas necesarias para ejecutar el pipeline completo de análisis transcriptómico y neurofisiológico de forma automatizada con Snakemake.

---

## ⚠️ Estado actual del `Snakefile`

> **⚠️ IMPORTANTE:** A fecha de la última actualización, el `Snakefile` y los flujos declarados en este directorio **no han sido validados completamente** debido a:
>
> * Cambios recientes en la estructura de carpetas y nombres de scripts.
> * Modificaciones en los argumentos y rutas que aún no han sido integradas al 100%.
> * Falta de pruebas de ejecución completas de principio a fin.

Por tanto, **no se garantiza su funcionamiento actual**, aunque refleja la lógica general del flujo implementado.

---

## 📁 Archivos presentes

* `Snakefile`: definición de las reglas del flujo (por dataset, pasos modulares)
* `config/config.yaml`: parámetros generales, rutas a datasets y opciones
* `workflow_ready.zip`: versión provisional empaquetada del flujo

---

## 🧪 Pendientes para validación completa

* Actualizar reglas con nombres de scripts definitivos
* Sincronizar con los README actuales de cada subcarpeta (`scripts/`)
* Probar la ejecución completa con `snakemake -n` y luego en seco (`-j 1`)

---

## 💡 Ejemplo de ejecución esperada (una vez validado)

```bash
conda activate tfm_mdd_framework
snakemake -s Snakefile -j 4 --use-conda --config dataset=GSE98793
```

---

Para cualquier duda sobre la integración o ejecución del flujo, contactar con el autor del proyecto o consultar los `README.md` en `scripts/` y `pipeline/`.

