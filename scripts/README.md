# 📂 scripts

Este directorio contiene todos los scripts que componen el framework de análisis multiómico.

Los scripts están organizados por función y tipo de ómica:

- `scripts/general/` → scripts reutilizables en cualquier ómica (descarga, preprocesamiento maestro, QC, etc.)
- `scripts/transcriptomics/` → scripts específicos para análisis de datos de expresión génica

---

## Convenciones generales

- Todos los scripts están preparados para ejecución desde terminal con `Rscript script.R`
- Se pueden integrar en flujos tipo Snakemake
- La estructura modular facilita su extensión a nuevas capas ómicas
