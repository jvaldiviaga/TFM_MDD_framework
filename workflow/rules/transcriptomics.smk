# Reglas de Snakemake para el flujo transcriptómico completo

# Wildcards explícitos para que Snakemake los detecte correctamente
wildcard_constraints:
    dataset = "[A-Za-z0-9_]+"

# 1. Descarga de datos GEO
rule download_geo_data:
    input:
        script = config["paths"]["scripts"] + "/transcriptomics/prepare_geo_dataset.sh"
    output:
        tarball = config["paths"]["data"] + "/transcriptomics/{dataset}/{dataset}_RAW.tar"
    params:
        dataset = "{dataset}"
    shell:
        "bash {input.script} {params.dataset}"

# 2. Descarga de metadatos
rule download_metadata_geo:
    input:
        script = config["paths"]["scripts"] + "/transcriptomics/download_metadata_geo.R"
    output:
        metadata = config["paths"]["data"] + "/transcriptomics/{dataset}/{dataset}_metadata.csv"
    params:
        dataset = "{dataset}"
    shell:
        "Rscript {input.script} {params.dataset}"

# 3. Preprocesamiento y normalización
rule preprocess_expression_data:
    input:
        script = config["paths"]["scripts"] + "/transcriptomics/preprocessing_master.R",
        raw_data = config["paths"]["data"] + "/transcriptomics/{dataset}/{dataset}_RAW.tar"
    output:
        normalized = config["paths"]["results"] + "/{dataset}/normalized/normalized_expression_combat.rds"
    log:
        config["paths"]["logs"] + "/{dataset}/preprocessing.log"
    params:
        dataset = "{dataset}"
    conda:
        config["conda_envs"]["transcriptomics"]
    shell:
        "Rscript {input.script} {params.dataset} > {log} 2>&1"

# 4. QC post-normalización
rule qc_post_normalization:
    input:
        script = config["paths"]["scripts"] + "/transcriptomics/qc_post_normalization.R",
        expr = config["paths"]["results"] + "/{dataset}/normalized/normalized_expression_combat.rds"
    output:
        interpret = config["paths"]["results"] + "/{dataset}/qc_post/qc_interpretacion.txt"
    log:
        config["paths"]["logs"] + "/{dataset}/qc_post.log"
    params:
        dataset = "{dataset}"
    conda:
        config["conda_envs"]["transcriptomics"]
    shell:
        "Rscript {input.script} {params.dataset} > {log} 2>&1"

# 5. Expresión diferencial
rule differential_expression:
    input:
        script = config["paths"]["scripts"] + "/transcriptomics/differential_expression.R",
        expr = config["paths"]["results"] + "/{dataset}/normalized/normalized_expression_combat.rds"
    output:
        deg = config["paths"]["results"] + "/{dataset}/deg/DEG_results.csv"
    log:
        config["paths"]["logs"] + "/{dataset}/deg.log"
    params:
        dataset = "{dataset}"
    conda:
        config["conda_envs"]["transcriptomics"]
    shell:
        "Rscript {input.script} {params.dataset} > {log} 2>&1"


# 6. Enriquecimiento funcional (solo si hay DEGs significativos)
rule enrichment_analysis:
    input:
        script = config["paths"]["scripts"] + "/transcriptomics/enrichment_analysis.R",
        deg = config["paths"]["results"] + "/{dataset}/deg/DEG_results.csv"
    output:
        go_csv = config["paths"]["results"] + "/{dataset}/enrichment/enrichment_results_GO.csv"
    log:
        config["paths"]["logs"] + "/{dataset}/enrichment.log"
    params:
        dataset = "{dataset}"
    conda:
        config["conda_envs"]["transcriptomics"]
    shell:
        """
        # Revisa si hay DEGs significativos antes de ejecutar el script
        if [ $(awk -F',' 'NR>1 && $5 < 0.05 {{print}}' {input.deg} | wc -l) -gt 0 ]; then
            Rscript {input.script} {params.dataset} > {log} 2>&1
        else
            echo "❌ No DEGs significativos. Enriquecimiento omitido." > {log}
            mkdir -p $(dirname {output.go_csv})
            touch {output.go_csv}
        fi
        """
# 7. Validación cruzada (logFC)
rule validate_logfc_direction:
    input:
        script = config["paths"]["scripts"] + "/general/validate_logfc_direction.R",
        deg = config["paths"]["results"] + "/{dataset}/deg/DEG_results.csv"
    output:
        val_file = config["paths"]["results"] + "/{dataset}/validation/logfc_validation.csv"
    log:
        config["paths"]["logs"] + "/{dataset}/validation.log"
    params:
        dataset = "{dataset}"
    conda:
        config["conda_envs"]["transcriptomics"]
    shell:
        "Rscript {input.script} {params.dataset} > {log} 2>&1"
