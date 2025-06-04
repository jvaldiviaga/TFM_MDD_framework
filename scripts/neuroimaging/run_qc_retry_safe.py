#!/usr/bin/env python3
# ============================================================
# Script: run_qc_retry_safe.py
#
# OBJETIVO:
# Reintentar el QC solo para sujetos fallidos previamente,
# usando el flag --safe_mode (crop + preload=False)
# ============================================================

import os
import subprocess

# Ruta a donde están los resultados
output_dir = os.path.expanduser("~/TFM_MDD/results/ds005356/qc_pre_filtering")
input_dir = os.path.expanduser("~/TFM_MDD/data/neuroimaging/ds005356/derivatives/MEGcombined")

# Leer lista de fallidos
failed_path = os.path.join(output_dir, "qc_failed_subjects.txt")
if not os.path.exists(failed_path):
    print("❌ No se encontró qc_failed_subjects.txt. No hay nada que reintentar.")
    exit(0)

with open(failed_path, "r") as f:
    lines = f.readlines()

# Procesar cada sujeto fallido
for line in lines:
    parts = line.strip().split(",")
    if not parts:
        continue
    fname = parts[0]
    if not fname.endswith("_meg.fif"):
        continue

    subject_id = fname.replace("_meg.fif", "")
    input_path = os.path.join(input_dir, fname)

    if not os.path.exists(input_path):
        print(f"⚠️ Archivo no encontrado: {input_path}. Saltando.")
        continue

    print(f"🔁 Reintentando en modo seguro: {fname}")
    cmd = [
        "python3", "qc_meg_raw.py",
        "--input_file", input_path,
        "--output_dir", output_dir,
        "--safe_mode"
    ]

    try:
        subprocess.run(cmd, check=True)
        print(f"✅ QC relanzado con éxito para {subject_id}")

        # Eliminar del archivo de fallidos
        lines = [l for l in lines if not l.startswith(fname)]
        with open(failed_path, "w") as f_out:
            f_out.writelines(lines)

    except subprocess.CalledProcessError as e:
        print(f"❌ Sigue fallando: {subject_id} → {e}")




