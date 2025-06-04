#!/usr/bin/env python3
# ============================================================
# Script: run_all_preprocessing.py
#
# OBJETIVO:
# Automatizar el preprocesamiento técnico de todos los archivos
# MEG crudos fusionados (*.fif), aplicando:
#   - Detección de canales malos
#   - Maxwell filtering
#   - Filtro notch (60 Hz)
#   - Filtro pasa banda (1–40 Hz)
#
# FLUJO:
# 1. Detecta todos los archivos *_meg.fif sin filtrar
# 2. Llama al script `preprocess_meg_subject.py` para cada uno
# 3. Guarda archivos *_meg_filtered.fif y log por sujeto
#
# IMPORTANTE:
# Este paso se debe ejecutar DESPUÉS del control de calidad visual (QC)
# realizado con `qc_meg_raw.py` para validar los datos crudos.
# ============================================================

import os               # Para manejo de rutas y archivos
import subprocess       # Para ejecutar otros scripts con parámetros
import pandas as pd     # Para leer archivos QC en CSV

# ========================== DEFINIR RUTAS ==========================

# Carpeta con los archivos crudos fusionados (*.fif)
input_dir = os.path.expanduser("~/TFM_MDD/data/neuroimaging/ds005356/derivatives/MEGcombined")

# Carpeta de salida con archivos ya filtrados y limpios
output_dir = os.path.expanduser("~/TFM_MDD/results/ds005356/preprocessing")
os.makedirs(output_dir, exist_ok=True)

# Log de sujetos que fallaron durante el preprocesamiento
failed_log_path = os.path.join(output_dir, "failed_subjects.txt")

# ========================== CARGA DEL QC PREVIO ==========================

qc_ready_path = os.path.join(os.path.dirname(output_dir), "qc_pre_filtering", "qc_ready_flag.csv")
if os.path.exists(qc_ready_path):
    qc_flags = pd.read_csv(qc_ready_path)
    qc_flags["subject"] = qc_flags["subject"].astype(str).str.strip()
    qc_flags.set_index("subject", inplace=True)
    print(f"📋 QC flags cargados para {len(qc_flags)} sujetos.")
else:
    qc_flags = pd.DataFrame()
    print("⚠️ No se encontró el archivo qc_ready_flag.csv. Se procesarán todos los sujetos.")

# ========================== PROCESAMIENTO EN LOTE ==========================

# Bucle sobre todos los archivos .fif en el input_dir
for fname in sorted(os.listdir(input_dir)):

    # Selección de archivos crudos válidos (evita re-procesar)
    if fname.endswith("_meg.fif") and not fname.endswith("_meg_clean.fif"):

        subject_id = fname.split("_")[0]  # sub-XXXXXX
        input_path = os.path.join(input_dir, fname)
        output_file = os.path.join(output_dir, f"{subject_id}_meg_clean.fif")

        # Si ya fue preprocesado, se omite
        if os.path.exists(output_file):
            print(f"⏩ {subject_id} ya tiene archivo limpio. Saltando.")
            continue

        # Aviso si el QC marcó como "no apto"
        if subject_id in qc_flags.index:
            flag = str(qc_flags.loc[subject_id, "qc_ready_for_preprocessing"]).lower()
            if flag in ["false", "0", "no"]:
                print(f"⚠️ {subject_id} fue marcado como NO listo en el QC previo, pero se procesará igualmente.")

        # Construcción del comando
        cmd = [
            "python3", "preprocess_meg_subject.py",
            "--input_file", input_path,
            "--output_dir", output_dir
        ]

        # Ejecución con control de errores
        print(f"🚀 Procesando: {subject_id}")
        try:
            subprocess.run(cmd, check=True)
            print(f"✅ {subject_id} procesado correctamente.")
        except subprocess.CalledProcessError as e:
            error_msg = e.stderr.strip() if e.stderr else "❌ Error sin salida estándar."
            print(f"❌ Error procesando {subject_id}: {error_msg}")
            with open(failed_log_path, "a") as f:
                f.write(f"{subject_id},{error_msg}\n")

# ========================== FINALIZACIÓN ==========================

print("✅ Preprocesamiento completado.")
print(f"📂 Resultados guardados en: {output_dir}")
print(f"📜 Log de errores (si los hubo): {failed_log_path}")
