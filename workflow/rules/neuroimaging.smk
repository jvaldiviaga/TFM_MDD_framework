# Reglas de Snakemake para el flujo MEG (neuroimagen funcional)
# Automatiza el control de calidad, preprocesamiento, epoching, extracción de features y estadísticas grupales
# Todas las reglas están explicadas para usuarios técnicos o clínicos sin experiencia informática

# =====================
# 1. QC PREVIO (RAW)
# =====================
rule qc_meg_raw:
    input:
        script = f"{config['paths']['scripts']}/neuroimaging/qc_meg_raw.py",
        raw_dir = f"{config['paths']['data']}/neuroimaging/{config['dataset_meg']}/derivatives/MEGcombined"
    output:
        summary = f"{config['paths']['results']}/{config['dataset_meg']}/qc/qc_psd_summary.csv"
    log:
        f"{config['paths']['logs']}/{config['dataset_meg']}/qc_raw.log"
    conda:
        config['conda_envs']['neuroimaging']
    shell:
        "python3 {input.script} --input_dir {input.raw_dir} --output_dir {config[\"paths\"][\"results\"]}/{config[\"dataset_meg\"]}/qc > {log} 2>&1"

# =======================================
# 2. PREPROCESAMIENTO INDIVIDUAL (clean)
# =======================================
rule preprocess_meg_subjects:
    input:
        script = f"{config['paths']['scripts']}/neuroimaging/preprocess_meg_subject.py",
        qc = rules.qc_meg_raw.output.summary
    output:
        done = f"{config['paths']['results']}/{config['dataset_meg']}/preprocessing/.done"
    log:
        f"{config['paths']['logs']}/{config['dataset_meg']}/preprocessing.log"
    conda:
        config['conda_envs']['neuroimaging']
    shell:
        "python3 {input.script} --input_dir {config[\"paths\"][\"data\"]}/neuroimaging/{config[\"dataset_meg\"]}/derivatives/MEGcombined --output_dir {config[\"paths\"][\"results\"]}/{config[\"dataset_meg\"]}/preprocessing > {log} 2>&1 && touch {output.done}"

# ===============================
# 3. EPOCHING DE LOS SUJETOS
# ===============================
rule epoch_meg_subjects:
    input:
        script = f"{config['paths']['scripts']}/neuroimaging/epoch_meg_subject.py",
        preproc_done = rules.preprocess_meg_subjects.output.done
    output:
        done = f"{config['paths']['results']}/{config['dataset_meg']}/epochs/.done"
    log:
        f"{config['paths']['logs']}/{config['dataset_meg']}/epoching.log"
    conda:
        config['conda_envs']['neuroimaging']
    shell:
        "python3 {input.script} --input_dir {config[\"paths\"][\"results\"]}/{config[\"dataset_meg\"]}/preprocessing --output_dir {config[\"paths\"][\"results\"]}/{config[\"dataset_meg\"]}/epochs > {log} 2>&1 && touch {output.done}"

# ================================
# 4. EXTRACCIÓN DE FEATURES
# ================================
rule extract_meg_features:
    input:
        script = f"{config['paths']['scripts']}/neuroimaging/extract_features_meg.py",
        epoch_done = rules.epoch_meg_subjects.output.done
    output:
        features = f"{config['paths']['results']}/{config['dataset_meg']}/features/features_all_subjects.csv"
    log:
        f"{config['paths']['logs']}/{config['dataset_meg']}/features.log"
    conda:
        config['conda_envs']['neuroimaging']
    shell:
        "python3 {input.script} --input_dir {config[\"paths\"][\"results\"]}/{config[\"dataset_meg\"]}/epochs --output_file {output.features} > {log} 2>&1"

# ============================================
# 5. ESTADÍSTICAS GRUPALES ENTRE MDD Y CTL
# ============================================
rule group_stats_meg:
    input:
        script = f"{config['paths']['scripts']}/neuroimaging/group_stats_meg.py",
        features = rules.extract_meg_features.output.features
    output:
        stats = f"{config['paths']['results']}/{config['dataset_meg']}/stats/group_stats_results.csv"
    log:
        f"{config['paths']['logs']}/{config['dataset_meg']}/group_stats.log"
    conda:
        config['conda_envs']['neuroimaging']
    shell:
        "python3 {input.script} --input_file {input.features} --output_file {output.stats} > {log} 2>&1"
