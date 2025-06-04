#!/usr/bin/env python3
# ============================================================
# Script: run_all_epoch_quality.py
#
# OBJETIVO:
# Ejecutar evaluate_epochs_quality.py para todos los archivos .epo.fif
# generados en el análisis previo. Detecta automáticamente los ensayos
# malos por métricas y guarda CSV + heatmap por sujeto.
# ============================================================

import os
import subprocess

# -------------------- PARÁMETROS --------------------
epochs_dir = os.path.expanduser("~/TFM_MDD/results/ds005356/epochs")
script_path = os.path.expanduser("~/TFM_MDD/scripts/neuroimaging/evaluate_epochs_quality.py")

# -------------------- LOOP POR ARCHIVOS --------------------
epo_files = sorted([f for f in os.listdir(epochs_dir) if f.endswith("_epo.fif")])
total = len(epo_files)

for idx, fname in enumerate(epo_files):
    subj_id = fname.split("_")[0]
    input_path = os.path.join(epochs_dir, fname)
    print(f"\n[{idx+1}/{total}] 🔎 Evaluando calidad: {subj_id}")

    cmd = [
        "python3", script_path,
        "--input_file", input_path,
        "--output_dir", epochs_dir,
        "--z_thresh", "3.0"  # Puedes ajustar si lo deseas
    ]

    log_file = os.path.join(epochs_dir, f"{subj_id}_quality.log")

    with open(log_file, "w") as logf:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            print(line, end="")
            logf.write(line)
        process.wait()

        if process.returncode != 0:
            print(f"❌ Error al evaluar calidad de {subj_id}. Revisa: {log_file}")
            continue

    print(f"✅ Evaluación completada para {subj_id} → log: {log_file}")

print("\n🎯 Calidad de epochs evaluada para todos los sujetos.")
