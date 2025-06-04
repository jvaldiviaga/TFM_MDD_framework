#!/bin/bash
# ===================================================================================================
# Script: run_all_epoching.py
# Objetivo:
# Epocar automáticamente todos los sujetos MEG del dataset (e.g., ds005356) usando el script
# 'epoch_meg_subject.py'
# ===================================================================================================

import os
import subprocess

# -------------------- PARÁMETROS --------------------
data_dir = os.path.expanduser("~/TFM_MDD/results/ds005356/preprocessing")
output_dir = os.path.expanduser("~/TFM_MDD/results/ds005356/epochs")
bids_root = os.path.expanduser("~/TFM_MDD/data/neuroimaging/ds005356")
script_path = os.path.expanduser("~/TFM_MDD/scripts/neuroimaging/epoch_meg_subject.py")

os.makedirs(output_dir, exist_ok=True)

# Precomputar lista de sujetos válidos
subjects = [f for f in sorted(os.listdir(data_dir)) if f.endswith("_meg_clean.fif")]
total = len(subjects)

# -------------------- LOOP POR SUJETOS --------------------
for idx, fname in enumerate(subjects):
    subj_id = fname.split("_")[0]  # ej: sub-M87100058
    input_file = os.path.join(data_dir, fname)
    output_file = os.path.join(output_dir, f"{subj_id}_epo.fif")
    log_file = os.path.join(output_dir, f"{subj_id}_epoching.log")

    # Saltar si ya está epocado
    if os.path.exists(output_file):
        print(f"⏩ {subj_id} ya epocado. Se omite.")
        continue

    # Detectar ruta al archivo de eventos .tsv
    events_file = os.path.join(
        bids_root, subj_id, "ses-01", "meg",
        f"{subj_id}_ses-01_task-pst_run-1_events.tsv"
    )
    
    if not os.path.exists(events_file):
        print(f"⚠️ No se encontró el archivo de eventos para {subj_id}. Se omite.")
        continue

    print(f"\n[{idx+1}/{total}] 🔁 Epocando: {subj_id}")

    cmd = [
        "python3", script_path,
        "--input_file", input_file,
        "--output_file", output_file,
        "--annotations_file", events_file,
        "--tmin", "-0.2",
        "--tmax", "0.8",
        "--export_events"
    ]

    with open(log_file, "w") as logf:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            print(line, end="")       # imprime en consola
            logf.write(line)          # guarda en log
        process.wait()
        if process.returncode != 0:
            print(f"❌ Error al epocar {subj_id}. Revisa el log: {log_file}")
            continue

    print(f"✅ Finalizado {subj_id} → Log: {log_file}")

print("\n🎯 Epocado finalizado para todos los sujetos.")


# ============================================================
# FUNCIONES DEL SCRIPT
# ============================================================
# • Automatiza el epocado por lote para datasets MEG grandes.
# • Se adapta al formato BIDS sin necesidad de renombrar ni convertir archivos.
# • Genera archivos .epo.fif compatibles con análisis funcional y multiclase.
# • Exporta eventos a .csv por sujeto y guarda logs individualizados.

# ============================================================
# LIMITACIONES Y POSIBLES MEJORAS
# ============================================================
# - No soporta aún análisis por sesiones múltiples (solo ses-01).
#   → Posible extensión si se desea procesar otras runs o tareas automáticamente.

# - No permite control fino de condiciones experimentales ni filtrado a nivel de script.
#   → El análisis condicional se recomienda hacerlo después, usando MNE-Python o scripts ERP/ERF.

# ============================================================
# ESTADO ACTUAL
# ============================================================
# El script está listo para usarse como paso intermedio en el pipeline MEG del TFM.
# Integra validación de eventos, exportación, rechazo automático y registro de logs.
# Es compatible con estructura BIDS y puede integrarse fácilmente en pipelines reproducibles.
