#!/usr/bin/env python3
# ============================================================
# Script: run_qc_all_meg.py
#
# OBJETIVO:
# Automatizar la ejecución por lotes del script `qc_meg_raw.py`,
# aplicándolo a todos los archivos MEG combinados (*.fif) que se
# encuentren en una carpeta específica del dataset (por ejemplo,
# derivados de `load_all_meg_subjects.py`).
#
# Este paso genera:
# - Curvas de densidad espectral multicanal (PSD)
# - Curvas de depuración espectral (perfil invertido)
# - CSVs interpretativos para análisis posterior (flags, etiquetas, etc.)
# - Logs detallados por sujeto
#
# MEJORAS:
# ✅ Verifica si el QC está completo (imagen PSD + CSVs) antes de ejecutar.
# ✅ Registra sujetos con errores o split incompleto.
# ✅ Genera logs de éxito/fallo.
#
# Y guarda logs de:
# - Éxitos
# - Fallos
# - Sujetos omitidos
# ============================================================
# Esta evaluación se realiza **antes del preprocesamiento** técnico (filtros, ICA),
# con el objetivo de detectar sujetos corruptos o de baja calidad y evitar
# aplicar transformaciones sobre datos inservibles.
# ============================================================

import os                # Para gestionar rutas y archivos
import argparse          # Para gestionar argumentos desde terminal
import subprocess        # Para ejecutar scripts externos (qc_meg_raw.py)
import mne               # Para verificar legibilidad de archivos .fif
import pandas as pd      # Para manejar CSVs

# ========================== PARÁMETROS DE ENTRADA ==========================
# Definimos los argumentos esperados:
parser = argparse.ArgumentParser(description="Lanzar QC por lotes sobre archivos .fif combinados.")
parser.add_argument('--input_dir', type=str, required=True,
                    help='Ruta con archivos *_meg.fif combinados (sin preprocesar)')
parser.add_argument('--output_dir', type=str, required=True,
                    help='Ruta donde se guardarán los gráficos, CSVs y logs generados por el QC')
args = parser.parse_args()

# Expandimos posibles tildes (~) en las rutas y nos aseguramos de que existan
input_dir = os.path.expanduser(args.input_dir)
output_dir = os.path.expanduser(args.output_dir)
os.makedirs(output_dir, exist_ok=True)

# ========================== FUNCIÓN DE CHEQUEO DE QC COMPLETO ==========================
# Esta función determina si el QC de un sujeto ya está completo.
# Requiere:
# - La imagen de densidad espectral (PSD)
# - Entrada en qc_psd_summary.csv (artefactos detectados)
# - Entrada en qc_bad_channels.csv (canales ruidosos detectados)

def check_qc_complete(output_dir, subject_id):
    psd_img = os.path.join(output_dir, f"{subject_id}_psd.png")                      # Imagen PSD
    summary_csv = os.path.join(output_dir, "qc_psd_summary.csv")                    # CSV general
    bads_csv = os.path.join(output_dir, "qc_bad_channels.csv")                      # CSV de canales malos

    # Si no existen los 3 archivos necesarios, el QC no está completo
    if not (os.path.exists(psd_img) and os.path.exists(summary_csv) and os.path.exists(bads_csv)):
        return False

    # Cargamos los CSVs y comprobamos si el sujeto está incluido en ambos
    try:
        summary_df = pd.read_csv(summary_csv)
        bads_df = pd.read_csv(bads_csv)
        return subject_id in summary_df["subject"].values and subject_id in bads_df["subject"].values
    except Exception:
        return False  # Error leyendo CSVs → tratamos como QC incompleto

# ========================== DETECCIÓN DE ARCHIVOS .FIF ==========================
# Mostramos mensaje de inicio y recopilamos todos los archivos *_meg.fif
print(f"📁 Buscando archivos .fif en: {input_dir}")
# ➤ Mensaje de consola para saber dónde se está buscando.

fif_files = sorted([f for f in os.listdir(input_dir) if f.endswith("_meg.fif")])
# ➤ Filtra todos los archivos en la carpeta input_dir que terminen en "_meg.fif"
# ➤ Esto evita accidentalmente seleccionar archivos preprocesados, epocados o filtrados
# ➤ Se ordenan alfabéticamente para tener un control reproducible del orden de ejecución

print(f"🔎 Se encontraron {len(fif_files)} archivos para procesar.\n")
# ➤ Muestra cuántos archivos .fif serán evaluados en este batch

# ========================== BUCLE PRINCIPAL ==========================
# Iteramos sobre cada archivo .fif encontrado
for fname in fif_files:
    input_path = os.path.join(input_dir, fname)               # Ruta completa del archivo .fif
    subject_id = fname.replace("_meg.fif", "")                # Extraemos ID del sujeto (ej. sub-XXXXX)

    # ⏩ Si el QC ya está completamente hecho, se omite
    if check_qc_complete(output_dir, subject_id):
        print(f"⏩ {subject_id} ya tiene QC completo. Se omite.")
        continue

    # ⚠️ Comprobamos que el archivo no esté corrupto ni sea un split incompleto
    try:
        raw = mne.io.read_raw_fif(input_path, preload=False, verbose="ERROR")
    except Exception as e:
        print(f"⚠️ Archivo ilegible: {fname} → {e}")
        with open(os.path.join(output_dir, "qc_skipped_split_or_corrupt.txt"), "a") as f:
            f.write(f"{fname},{str(e)}\n")
        continue  # Pasamos al siguiente sujeto

    # 🚀 Ejecutamos el QC individual para este sujeto usando qc_meg_raw.py
    cmd = [
        "python3", "qc_meg_raw.py",
        "--input_file", input_path,
        "--output_dir", output_dir
    ]

    print(f"🚀 Ejecutando QC para: {fname}")
    try:
        subprocess.run(cmd, check=True)  # check=True lanza error si el QC falla
        print(f"✅ QC finalizado con éxito para: {fname}")
        with open(os.path.join(output_dir, "qc_successful_subjects.txt"), "a") as f:
            f.write(f"{fname}\n")

    except subprocess.CalledProcessError as e:
        print(f"❌ Error durante QC de {fname}: {e}")
        with open(os.path.join(output_dir, "qc_failed_subjects.txt"), "a") as f:
            f.write(f"{fname},{str(e)}\n")
        continue  # Seguimos con el siguiente sujeto

# ========================== FINALIZACIÓN ==========================
# Mensaje de cierre al terminar todos los sujetos
print("\n🎯 QC por lotes completado.")
print(f"📄 Resultados, logs e informes en: {output_dir}")