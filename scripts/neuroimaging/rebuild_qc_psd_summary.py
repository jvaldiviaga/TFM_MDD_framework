#!/usr/bin/env python3

import os
import pandas as pd

# Ruta a la carpeta donde están los CSVs generados por qc_meg_raw.py
qc_dir = os.path.expanduser("~/TFM_MDD/results/ds005356/qc_pre_filtering")

# Carga de los archivos auxiliares que contienen la información original
binary_df = pd.read_csv(os.path.join(qc_dir, "qc_post_binary_summary.csv"))
events_df = pd.read_csv(os.path.join(qc_dir, "qc_event_check.csv"))
bads_df = pd.read_csv(os.path.join(qc_dir, "qc_bad_channels.csv"))

# Limpieza de IDs
for df in [binary_df, events_df, bads_df]:
    df["subject"] = df["subject"].astype(str).str.strip().str.lower()

# Merge robusto por sujeto
df = binary_df.merge(events_df, on="subject", how="outer")
df = df.merge(bads_df[["subject", "n_bad", "n_flat"]].rename(columns={
    "n_bad": "n_bad_channels_pre",
    "n_flat": "n_flat_channels_pre"
}), on="subject", how="left")

# Simulación de métricas espectrales faltantes (si no tienes `low_freq_power_pre`, etc.)
df["low_freq_power_pre"] = 0.0
df["high_freq_power_pre"] = 0.0
df["psd_max"] = 0.0
df["psd_median"] = 0.0

# Renombramos y ordenamos columnas como las espera el script de preproc
df_renamed = df.rename(columns={
    "muscle_noise_>50Hz": "muscle_present",
    "low_freq_drift_<1.5Hz": "drift_present",
    "clipping_detected": "has_clipping",
    "stim_channel_present": "events_ok"
})[[
    "subject", "drift_present", "muscle_present", "has_clipping",
    "low_freq_power_pre", "high_freq_power_pre", "psd_max", "psd_median",
    "n_bad_channels_pre", "n_flat_channels_pre", "events_ok"
]]

# Exporta el CSV corregido
output_path = os.path.join(qc_dir, "qc_psd_summary.csv")
df_renamed.to_csv(output_path, index=False)
print(f"✅ Archivo reconstruido: {output_path}")
