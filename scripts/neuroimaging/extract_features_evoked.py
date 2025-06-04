#!/usr/bin/env python3
# ============================================================
# Script: extract_features_evoked.py (versión biomarcadores)
# Objetivo: Extraer métricas clínicas como GFP, P300, slope, RWP por canal y condición
# ============================================================

import os
import argparse
import mne
import numpy as np
import pandas as pd

# -------------------------------
# FUNCIONES AUXILIARES
# -------------------------------
def gfp(data):  # GFP: raíz de varianza entre canales por tiempo
    return np.sqrt((data ** 2).mean(axis=0))

def slope_25_75(sig, times):
    peak_idx = np.argmax(sig)
    peak_val = sig[peak_idx]
    lower = 0.25 * peak_val
    upper = 0.75 * peak_val

    try:
        idx_25 = np.where(sig[:peak_idx] >= lower)[0][0]
        idx_75 = np.where(sig[:peak_idx] >= upper)[0][0]
        t25 = times[idx_25]
        t75 = times[idx_75]
        return (upper - lower) / (t75 - t25) if (t75 - t25) > 0 else np.nan
    except:
        return np.nan

def rwp_duration(sig, times):
    peak_idx = np.argmax(sig)
    half_peak = 0.5 * sig[peak_idx]
    above = sig >= half_peak
    indices = np.where(above)[0]
    if len(indices) < 2:
        return np.nan
    return times[indices[-1]] - times[indices[0]]

def p300_metrics(sig, times):
    win = (times >= 0.25) & (times <= 0.5)
    peak_val = np.max(sig[win])
    latency = times[win][np.argmax(sig[win])]
    return peak_val, latency

# -------------------------------
# ARGUMENTOS
# -------------------------------
parser = argparse.ArgumentParser(description="Extrae biomarcadores evocados por canal y condición.")
parser.add_argument("--input_file", required=True, help="Archivo .epo_filtered.fif")
parser.add_argument("--output_file", required=True, help="CSV de salida")
args = parser.parse_args()

# -------------------------------
# CARGA
# -------------------------------
print(f"📂 Cargando archivo de epochs: {args.input_file}")
epochs = mne.read_epochs(args.input_file, preload=True)
subject_id = os.path.basename(args.input_file).split("_")[0].replace("sub-", "")
times = epochs.times
print(f"🧠 Sujeto: {subject_id}")
print(f"🔢 Ensayos disponibles: {len(epochs)}")

# -------------------------------
# CONDICIONES
# -------------------------------
condiciones = {
    "cond_01": "FB/loss", "cond_02": "FB/win",
    "cond_03": "cue/AB",  "cond_04": "cue/BA",
    "cond_05": "cue/CD",  "cond_06": "cue/DC",
    "cond_07": "cue/EF",  "cond_08": "cue/FE"
}

features_flat = {}

# -------------------------------
# EXTRACCIÓN POR CONDICIÓN
# -------------------------------
for cond_label, event_key in condiciones.items():
    if event_key in epochs.event_id:
        try:
            epochs_cond = epochs[event_key]
            if len(epochs_cond) < 5:
                print(f"⚠️ {cond_label} ({event_key}): pocos ensayos ({len(epochs_cond)}). Saltando.")
                continue

            evoked = epochs_cond.average()
            data = evoked.data  # (n_channels, n_times)
            ch_names = evoked.ch_names

            gfp_vec = gfp(data)

            for i, chan in enumerate(ch_names):
                sig = data[i, :]
                gfp_chan = gfp_vec
                p300_amp, p300_lat = p300_metrics(sig, times)
                slope = slope_25_75(sig, times)
                duration = rwp_duration(sig, times)

                features_flat[f"{cond_label}_P300_amp_{chan}"] = p300_amp
                features_flat[f"{cond_label}_P300_latency_{chan}"] = p300_lat
                features_flat[f"{cond_label}_slope_{chan}"] = slope
                features_flat[f"{cond_label}_rwp_duration_{chan}"] = duration
                features_flat[f"{cond_label}_GFP_max_{chan}"] = np.max(gfp_chan)

            print(f"✅ {cond_label} ({event_key}): {len(epochs_cond)} ensayos analizados.")
        except Exception as e:
            print(f"❌ Error en {cond_label} ({event_key}): {e}")
    else:
        print(f"⚠️ {cond_label} ({event_key}) no presente.")

# -------------------------------
# GUARDADO
# -------------------------------
if features_flat:
    df_out = pd.DataFrame([features_flat])
    df_out.insert(0, "subject", subject_id)
    df_out.to_csv(args.output_file, index=False)
    print(f"💾 Features guardados en: {args.output_file}")
else:
    print("❌ No se extrajeron features válidos.")
