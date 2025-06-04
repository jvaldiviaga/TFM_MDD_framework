# Reglas de Snakemake para la validación cruzada inter-ómica (transcriptómica + neuroimagen)
# Estas reglas automatizan la comparación entre genes diferenciales y señales MEG, integrando biomarcadores finales

# ========================================
# 1. VALIDACIÓN ENTRE GENES Y FEATURES MEG
# ========================================
rule validate_interomics:
    input:
        script = f"{config['paths']['scripts']}/general/validate_interomics.R",
        genes = f"{config['paths']['results']}/transcriptomics/VALIDATION_INTER_ENRICHMENT/inter_dataset_robustness_summary.csv",
        features = f"{config['paths']['results']}/{config['dataset_meg']}/stats/group_stats_results.csv"
    output:
        table = f"{config['paths']['results']}/interomics/association_table.csv"
    log:
        f"{config['paths']['logs']}/interomics/validate_interomics.log"
    conda:
        config['conda_envs']['interomics']
    shell:
        "Rscript {input.script} > {log} 2>&1"

# ===============================
# 2. INFORME FINAL MULTI-ÓMICO
# ===============================
rule generate_biomarker_report:
    input:
        script = f"{config['paths']['scripts']}/general/generate_biomarker_report.R",
        genes = f"{config['paths']['results']}/transcriptomics/VALIDATION_INTER_ENRICHMENT/inter_dataset_robustness_summary.csv",
        features = f"{config['paths']['results']}/{config['dataset_meg']}/stats/group_stats_results.csv"
    output:
        report = f"{config['paths']['results']}/interomics/biomarker_report.pdf"
    log:
        f"{config['paths']['logs']}/interomics/biomarker_report.log"
    conda:
        config['conda_envs']['interomics']
    shell:
        "Rscript {input.script} > {log} 2>&1"
