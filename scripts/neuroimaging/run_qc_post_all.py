#!/usr/bin/env python3
# ============================================================
# Script: run_qc_post_all.py
#
# OBJETIVO:
# Ejecutar el QC post-preprocesado para todos los sujetos de forma automatizada,
# llamando al nuevo script qc_meg_post.py (renombrado desde v2).
# Incluye generación de imágenes PSD, métricas por banda, evaluación de artefactos
# y análisis por grupo si hay metadatos disponibles.
# ============================================================

import os
import subprocess

# -------------------- PARÁMETROS --------------------

# Ruta a los archivos *_meg_clean.fif
input_dir = os.path.expanduser("~/TFM_MDD/results/ds005356/preprocessing")

# Directorio donde guardar los resultados del QC post
output_dir = os.path.expanduser("~/TFM_MDD/results/ds005356/qc_post")

# Ruta al script principal de QC post (renombrado qc_meg_post.py)
qc_script = os.path.expanduser("~/TFM_MDD/scripts/neuroimaging/qc_meg_post.py")

# -------------------- VERIFICACIONES --------------------
# Creamos el directorio de salida si no existe
os.makedirs(output_dir, exist_ok=True)

# Verificamos que el script de QC exista
if not os.path.exists(qc_script):
    raise FileNotFoundError(f"❌ No se encuentra el script de QC: {qc_script}")

# -------------------- COMANDO DE EJECUCIÓN MULTISUJETO --------------------

# Ruta al QC previo
qc_pre_path = os.path.expanduser("~/TFM_MDD/results/ds005356/qc_pre_filtering/qc_psd_summary.csv")
if not os.path.exists(qc_pre_path):
    raise FileNotFoundError(f"❌ No se encuentra el QC previo: {qc_pre_path}")

# Ejecutar qc_meg_post.py sujeto por sujeto
for fname in sorted(os.listdir(input_dir)):
    if fname.endswith("_meg_clean.fif"):
        input_file = os.path.join(input_dir, fname)
        print(f"\n🚀 Procesando sujeto: {fname}")
        cmd = [
            "python3", qc_script,
            "--input_file", input_file,
            "--output_dir", output_dir,
            "--qc_pre_path", qc_pre_path
        ]
        print("🔗 Comando:", " ".join(cmd))
        subprocess.run(cmd)
