#!/usr/bin/env python3
# ============================================================
# Script: load_all_meg_subjects.py (versión avanzada)
# Objetivo: Cargar y fusionar automáticamente los archivos .fif
#           divididos (split-*) de cada sujeto del dataset MEG
#           (formato BIDS, como ds005356), validando su
#           consistencia técnica y generando un .fif combinado
#           por sujeto, junto con un resumen en .csv y log.
# ============================================================
# FLUJO Y OBJETIVO:
# Cargar todos los sujetos MEG (formato BIDS) del dataset ds005356,
# detectar y fusionar automáticamente los archivos .fif split-*,
# validar consistencia técnica (canales, frecuencia),
# guardar un archivo combinado .fif por sujeto,
# y generar un resumen .csv + log .txt de todo el proceso.
#
# PUNTOS FUERTES:
# - Robusto, generalizable y automatizable
# - Detecta errores de entrada sin romperse
# - Compatible con Snakemake, cluster o ejecución por lotes
#
# ⚠️IMITACIONES ACTUALES:
# - Solo evalúa 1 sesión (`--session`)
# - No permite filtrar sujetos por ahora
# - No evalúa calidad de señal aún (eso se hace después)
# ============================================================

import os               # Gestión de rutas y archivos
import mne              # Procesamiento de datos MEG/EEG
import glob             # Búsqueda de archivos tipo wildcard (*.fif)
import argparse         # Argumentos por consola para mayor flexibilidad
import pandas as pd     # Resumen final en .csv
from datetime import datetime  # Marca temporal para logs

# ============================================================
#  1. PARSEO DE ARGUMENTOS DINÁMICOS (flexibilidad del script)
# ============================================================

parser = argparse.ArgumentParser(description="Fusionar archivos .fif split por sujeto (MEG, BIDS).")
parser.add_argument('--data_root', type=str, required=True,
                    help='Ruta raíz del dataset (ej. ~/TFM_MDD/data/neuroimaging/ds005356)')
parser.add_argument('--session', type=str, default='ses-01',
                    help='Sesión BIDS (por defecto: ses-01)')
parser.add_argument('--log_file', type=str, default='meg_combined_log.txt',
                    help='Nombre del archivo de log textual (opcional)')
args = parser.parse_args()

# Resolución de rutas
data_root = os.path.expanduser(args.data_root)
session = args.session
log_path = os.path.join(data_root, args.log_file)

# ============================================================
#  ESTRUCTURA DE SALIDA
# ============================================================

# Carpeta de salida para guardar los archivos combinados
deriv_dir = os.path.join(data_root, "derivatives", "MEGcombined")
os.makedirs(deriv_dir, exist_ok=True)  # ⚙️ Se crea si no existe

summary = []        # Lista con info por sujeto para el CSV
log_lines = []      # Lista de líneas para el log .txt

# ============================================================
#  2. DETECTAR CARPETAS DE SUJETOS (formato sub-XXXXX)
# ============================================================

subjects = sorted([
    d for d in os.listdir(data_root)
    if d.startswith("sub-") and os.path.isdir(os.path.join(data_root, d))
])

log_lines.append(f"[{datetime.now()}] 📁 Detectados {len(subjects)} sujetos en {data_root}\n")

# ============================================================
#  3. BUCLE PRINCIPAL — CARGA Y FUSIÓN DE CADA SUJETO
# ============================================================

for subject in subjects:
    print(f"\n🔄 Procesando sujeto: {subject}")
    subject_dir = os.path.join(data_root, subject, session, "meg")
    out_path = os.path.join(deriv_dir, f"{subject}_meg.fif")

    # Verificación anticipada: ¿ya existe un .fif combinado válido?
    if os.path.exists(out_path):
        try:
            raw_check = mne.io.read_raw_fif(out_path, preload=False, verbose='ERROR')

            # Detectamos si es un archivo split (incompleto)
            if any("split" in fn for fn in raw_check.filenames):
                raise ValueError("Archivo dividido detectado")

            duration_check = raw_check.n_times / raw_check.info["sfreq"]
            msg = f"⏩ {subject} ya combinado previamente ({round(duration_check, 2)} s). Se omite."
            print(msg)
            log_lines.append(msg)
            summary.append({
                "subject": subject,
                "status": "skipped_existing",
                "n_files": "already_done",
                "duration_sec": round(duration_check, 2),
                "n_channels": raw_check.info["nchan"],
                "error_msg": ""
            })
            continue  # Saltar a siguiente sujeto

        except Exception as e:
            print(f"⚠️ Archivo existente es inválido (split roto o corrupto). Se regenerará → {e}")
            log_lines.append(f"{subject} → archivo inválido, se rehace: {e}")

    # Verificar carpeta MEG
    if not os.path.exists(subject_dir):
        msg = f"⚠️ {subject}: carpeta MEG no encontrada en {subject_dir}"
        print(msg)
        log_lines.append(msg)
        summary.append({
            "subject": subject, "status": "no_meg_folder", "n_files": 0,
            "duration_sec": None, "n_channels": None, "error_msg": "no_meg_folder"
        })
        continue

    # Buscar archivos split del sujeto
    fif_files = sorted(glob.glob(os.path.join(subject_dir, "*_meg.fif")))
    if len(fif_files) == 0:
        msg = f"⚠️ {subject}: no se encontraron archivos .fif"
        print(msg)
        log_lines.append(msg)
        summary.append({
            "subject": subject, "status": "no_fif_files", "n_files": 0,
            "duration_sec": None, "n_channels": None, "error_msg": "no_fif_files"
        })
        continue

    # ========== COMBINACIÓN ==========
    try:
        print(f"📥 Archivos encontrados: {len(fif_files)} → combinando...")
        raws = []
        ref_info = None

        for f in fif_files:
            print(f"  ↪️ {os.path.basename(f)}")
            raw = mne.io.read_raw_fif(f, preload=True, verbose='ERROR')

            if ref_info is None:
                ref_info = (raw.info['nchan'], raw.info['sfreq'])
            else:
                if raw.info['nchan'] != ref_info[0] or raw.info['sfreq'] != ref_info[1]:
                    raise ValueError(f"Inconsistencia entre archivos: {f} → {raw.info['nchan']} ch / {raw.info['sfreq']} Hz")

            raws.append(raw)

        raw_combined = mne.concatenate_raws(raws)

        # Guardar archivo combinado
        print(f"💾 Guardando archivo en: {out_path}")
        raw_combined.save(out_path, overwrite=True)

        duration = raw_combined.n_times / raw_combined.info["sfreq"]
        print(f"✅ {subject} combinado correctamente ({round(duration, 2)} s)")

        summary.append({
            "subject": subject,
            "status": "ok",
            "n_files": len(fif_files),
            "duration_sec": round(duration, 2),
            "n_channels": raw_combined.info["nchan"],
            "error_msg": ""
        })
        log_lines.append(f"{subject} → combinación exitosa")

    except Exception as e:
        msg = f"❌ {subject} error crítico: {str(e)}"
        print(msg)
        log_lines.append(msg)
        summary.append({
            "subject": subject,
            "status": "error",
            "n_files": len(fif_files),
            "duration_sec": None,
            "n_channels": None,
            "error_msg": str(e)
        })

# ============================================================
#  4. GUARDAR SALIDAS (CSV + LOG)
# ============================================================

# Guardar el resumen como archivo tabular
df = pd.DataFrame(summary)
csv_path = os.path.join(data_root, "meg_combined_summary.csv")
df.to_csv(csv_path, index=False)

# Guardar el log en texto plano
with open(log_path, "w") as f:
    f.write("\n".join(log_lines))

print(f"📄 Resumen guardado en: {csv_path}")
print(f"🗒️  Log guardado en: {log_path}")
