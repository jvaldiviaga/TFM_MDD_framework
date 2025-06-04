# scripts/neuroimaging/README.md

Este directorio contiene todos los scripts relacionados con el análisis de neuroimagen funcional MEG del dataset ds005356 (OpenNeuro/NEMAR). El flujo automatizado permite fusionar archivos, realizar control de calidad (QC), preprocesamiento, segmentación, extracción de características y análisis estadístico, todo ello con fines de detección de biomarcadores funcionales del MDD.

---

## 🧪 Flujo MEG completo — Scripts y utilidad

| Script                            | Descripción                                                                         | Inputs                                        | Outputs                                                | Ejemplo de uso                                                                        |
| --------------------------------- | ----------------------------------------------------------------------------------- | --------------------------------------------- | ------------------------------------------------------ | ------------------------------------------------------------------------------------- |
| `load_all_meg_subjects.py`        | Fusiona automáticamente los archivos MEG por sujeto (raw\.fif)                      | `data_root/` con subcarpetas por sujeto       | Archivos combinados en `MEGcombined/` y resumen `.csv` | `python load_all_meg_subjects.py --data_root data/neuroimaging/ds005356`              |
| `run_qc_all_meg.py`               | Ejecuta QC espectral visual masivo sobre archivos combinados                        | `MEGcombined/`                                | `qc_pre_filtering/` con imágenes PSD por sujeto        | `python run_qc_all_meg.py --input_dir ... --output_dir ...`                           |
| `analyze_batch_meg_raw.py`        | Evalúa el efecto batch técnico antes del preprocesamiento (scanner, duración, etc.) | Archivos combinados `.fif`, metadatos `.xlsx` | Gráficos PCA, dendrogramas, logs en `analysis_batch/`  | `python analyze_batch_meg_raw.py --data_dir ... --metadata_file ... --output_dir ...` |
| `run_all_preprocessing.py`        | Aplica Maxwell filter, ICA y otros pasos previos                                    | Archivos combinados `.fif`                    | Archivos `*_meg_clean.fif` por sujeto                  | `python run_all_preprocessing.py`                                                     |
| `run_qc_post_all.py`              | Evalúa calidad de los archivos limpios (`*_meg_clean.fif`)                          | Archivos preprocesados                        | `qc_post_summary.csv`, figuras PSD                     | `python run_qc_post_all.py`                                                           |
| `run_all_epoching.py`             | Realiza segmentación temporal en función de eventos                                 | Archivos limpios                              | Archivos `.fif` epocados (`*_epo.fif`)                 | `python run_all_epoching.py`                                                          |
| `run_all_epoch_quality.py`        | Evalúa calidad de los epochs (duración, picos, artefactos)                          | Epochs                                        | Métricas por sujeto en `.csv`                          | `python run_all_epoch_quality.py`                                                     |
| `run_all_filtering.py`            | Aplica filtrado adaptativo (notch, high-pass, etc.)                                 | Epochs                                        | Epochs filtrados                                       | `python run_all_filtering.py`                                                         |
| `run_all_extract_features.sh`     | Extrae features funcionales tipo RWP, GFP, P300, etc.                               | Epochs filtrados                              | Archivos `.csv` por sujeto con features                | `bash run_all_extract_features.sh`                                                    |
| `combine_features_safely.py`      | Une todos los archivos de features en una tabla general                             | Carpeta con `.csv` de features por sujeto     | `features_all_subjects.csv`                            | `python combine_features_safely.py --input_dir ... --output_file ...`                 |
| `run_merge_with_clinical.sh`      | Automatiza la unión de features MEG con datos clínicos                              | Features y metadatos                          | Tabla final para análisis grupal                       | `bash run_merge_with_clinical.sh`                                                     |
| `merge_features_with_clinical.py` | Script ejecutado dentro del `.sh` anterior                                          | Igual que anterior                            | Igual que anterior                                     | -                                                                                     |
| `run_group_stats_meg.sh`          | Lanza el análisis de comparación funcional entre grupos                             | Tabla combinada                               | Resultados estadísticos, gráficas, logs                | `bash run_group_stats_meg.sh`                                                         |

---

## 🔧 Requisitos

* Entorno Conda: `tfm_meg_env`
* Biblioteca base: `MNE-Python`
* Otros requerimientos especificados en `environment_meg.yml`

Para activar el entorno antes de cualquier ejecución:

```bash
conda activate tfm_meg_env
```

---

Todos los scripts han sido validados en el dataset `ds005356`. El flujo permite integración posterior con capas moleculares para análisis multiómico.
