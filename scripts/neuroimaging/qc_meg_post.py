#!/usr/bin/env python3
# ============================================================
# Script: qc_meg_post.py
#
# OBJETIVO:
# Evaluar la calidad técnica de los archivos MEG ya preprocesados (*.fif)
# Detectar si los filtros aplicados han eliminado efectivamente los artefactos
# detectados en el QC previo (60Hz, drift, músculo, clipping, etc.).
# ============================================================

import os
import argparse
import mne
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from scipy.signal import find_peaks
import json

# ========================== PARÁMETROS DE ENTRADA ==========================

parser = argparse.ArgumentParser(description="Evaluación post-preprocesamiento MEG.")
parser.add_argument('--input_file', type=str, required=True, help='Archivo .fif ya preprocesado')
parser.add_argument('--output_dir', type=str, required=True, help='Directorio para guardar resultados')
parser.add_argument('--qc_pre_path', type=str, required=True, help='Ruta al qc_psd_summary.csv preprocesado')
args = parser.parse_args()

input_path = os.path.expanduser(args.input_file)
output_dir = os.path.expanduser(args.output_dir)
qc_pre_path = os.path.expanduser(args.qc_pre_path)
os.makedirs(output_dir, exist_ok=True)
subject_id = os.path.basename(input_path).split("_")[0]

# ========================== CARGA DEL LOG DE PREPROCESAMIENTO ==========================
log_preproc = os.path.join(os.path.dirname(input_path), "log_preproc_meg.csv")

ica_applied = "no"
n_ica_removed = 0
has_eog = False
has_ecg = False

if os.path.exists(log_preproc):
    try:
        df_log = pd.read_csv(log_preproc)
        row_log = df_log[df_log["subject"] == subject_id]
        if not row_log.empty:
            ica_applied = row_log["ica_applied"].values[0]
            n_ica_removed = int(row_log.get("ica_n_components_removed", pd.Series([0])).values[0])
            has_eog = row_log.get("eog_channels_present", pd.Series(["no"])).values[0] == "yes"
            has_ecg = row_log.get("ecg_channels_present", pd.Series(["no"])).values[0] == "yes"
    except Exception as e:
        print(f"⚠️ No se pudo leer el log de preprocesamiento: {e}")

# ========================== CARGA DE RAW Y GENERACIÓN PSD ==========================
print(f"📂 Cargando archivo preprocesado: {input_path}")
try:
    raw = mne.io.read_raw_fif(input_path, preload=True)
except Exception as e:
    print(f"❌ No se pudo cargar el archivo: {e}")
    exit(1)

fig_psd = raw.plot_psd(fmax=65, average=False, show=False)
psd_img_path = os.path.join(output_dir, f"{subject_id}_psd_post.png")
fig_psd.savefig(psd_img_path, dpi=150)
plt.close(fig_psd)

# ========================== ANÁLISIS DE IMAGEN PSD ==========================
# Analiza imagen PSD → perfil invertido → detección de picos → artefactos

print("🔬 Analizando imagen PSD post-preprocesamiento...")
img = Image.open(psd_img_path).convert('L')
arr = np.array(img)
region = arr[150:350, 200:800]
profile = (region.max() - region.mean(axis=0))
adaptive_threshold = max(10, np.percentile(profile, 90) * 0.1)
peaks, _ = find_peaks(profile, prominence=adaptive_threshold)
freq_max = 65
peak_freqs = [round(freq_max * i / (len(profile) - 1), 1) for i in peaks]

# Guardar los picos detectados por PSD
pd.DataFrame({"frequency": peak_freqs}).to_csv(
    os.path.join(output_dir, f"{subject_id}_peaks_post.csv"), index=False)

# ========================== FASE FINAL: INTERPRETACIÓN AUTOMÁTICA ==========================
# Comparamos estado pre vs post y determinamos si el preprocesamiento ha sido efectivo

try:
    df_pre = pd.read_csv(qc_pre_path)
    row_pre = df_pre[df_pre["subject"] == subject_id]
except Exception as e:
    print(f"❌ Error al leer QC previo: {e}")
    row_pre = pd.DataFrame()

interpretacion = []
status_flags = []

# 1. Pico a 60 Hz — solo evaluamos si persiste, ya que el QC previo aplicó notch antes de medir
has_60hz_now = any(59.8 <= f <= 60.2 for f in peak_freqs)
print(f"🔎 Picos detectados post-PSD: {peak_freqs}")
if not has_60hz_now:
    interpretacion.append("✔️ sin pico 60 Hz tras preprocesamiento")
else:
    interpretacion.append("⚠️ persistencia del pico 60 Hz")
    status_flags.append(False)

# 2. Drift
drift_post = profile[:20].mean() > profile[20:40].mean() + 2
drift_pre = row_pre.get("drift_present", pd.Series([False])).values[0]
if drift_pre and not drift_post:
    interpretacion.append("✔️ drift corregido")
    status_flags.append(True)
elif drift_post:
    interpretacion.append("⚠️ drift residual")
    status_flags.append(False)

# 3. Ruido muscular
muscle_post = any(52 <= f <= 60 for f in peak_freqs)
muscle_pre = row_pre.get("muscle_present", pd.Series([False])).values[0]
if muscle_pre and not muscle_post:
    interpretacion.append("✔️ ruido muscular corregido")
    status_flags.append(True)
elif muscle_post:
    interpretacion.append("⚠️ ruido muscular persistente")
    status_flags.append(False)

# 4. Canales ruidosos
n_bad_pre = row_pre.get("n_bad_channels_pre", pd.Series([0])).values[0]
if n_bad_pre > 3 and len(peak_freqs) < 3:
    interpretacion.append("✔️ canales estabilizados tras Maxwell")
    status_flags.append(True)
elif len(peak_freqs) > 5:
    interpretacion.append("⚠️ perfil PSD inestable → revisar canales")
    status_flags.append(False)

# 5. ICA aplicada
ica_ok = ica_applied == "yes" and n_ica_removed > 0
if has_eog and not ica_ok:
    interpretacion.append("⚠️ posible EOG sin corregir")
    status_flags.append(False)
if has_ecg and not ica_ok:
    interpretacion.append("⚠️ posible ECG sin corregir")
    status_flags.append(False)
if ica_ok:
    interpretacion.append("✔️ ICA aplicada correctamente")
    status_flags.append(True)

# Evaluación global
if all(status_flags):
    qc_post_status = "🟢 corregido"
elif any(status_flags):
    qc_post_status = "🟡 parcialmente corregido"
else:
    qc_post_status = "🔴 sin corregir"

# Forzamos booleanos a tipo nativo Python (para evitar error en json)
has_60hz_now = bool(has_60hz_now)
drift_post = bool(drift_post)
muscle_post = bool(muscle_post)

# ========================== EXPORTACIÓN FINAL ==========================
summary_str = " + ".join(interpretacion)
print(f"🧠 Interpretación final para {subject_id}: {summary_str} | Estado: {qc_post_status}")

# Exportar interpretación como CSV
pd.DataFrame([{
    "subject": subject_id,
    "qc_post_interpretation": summary_str,
    "qc_post_status": qc_post_status
}]).to_csv(
    os.path.join(output_dir, "qc_post_interpretation.csv"),
    mode='a', index=False,
    header=not os.path.exists(os.path.join(output_dir, "qc_post_interpretation.csv"))
)

# Exportar log como JSON para trazabilidad futura
with open(os.path.join(output_dir, f"{subject_id}_qc_post_log.json"), "w") as f:
    json.dump({
        "subject": subject_id,
        "has_60hz_now": has_60hz_now,
        "drift_post": drift_post,
        "muscle_post": muscle_post,
        "ica_applied": ica_applied,
        "n_ica_removed": n_ica_removed,
        "qc_post_status": qc_post_status
    }, f, indent=4)
