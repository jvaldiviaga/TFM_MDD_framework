#!/usr/bin/env python3
# ============================================================
# Script: evaluate_epochs_quality.py
#
# OBJETIVO:
# Evaluar la calidad de cada ensayo epocado a partir de métricas reales
# y generar un CSV con etiquetas "good"/"bad", junto con un heatmap visual.
#
# FUNCIONALIDADES:
# 1. Carga un archivo .epo.fif
# 2. Calcula por ensayo:
#    - RMS (Root Mean Square)
#    - Amplitud máxima
#    - Desviación estándar
#    - Z-score global
# 3. Marca como outliers aquellos que excedan umbral (por defecto: RMS > 3*SD)
# 4. Exporta:
#    - CSV con columnas: ensayo, métricas, calidad
#    - Heatmap .png para revisión visual rápida
# ============================================================

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mne

# ---------------------- PARÁMETROS DE ENTRADA ----------------------
parser = argparse.ArgumentParser(description="Evaluar calidad de epochs MEG")
parser.add_argument('--input_file', required=True, help='Archivo .epo.fif de entrada')
parser.add_argument('--output_dir', required=True, help='Directorio de salida para CSV y heatmap')
parser.add_argument('--z_thresh', type=float, default=3.0, help='Umbral z-score para marcar outliers (default=3.0)')
args = parser.parse_args()

# ---------------------- CARGA DE EPOCHS ----------------------
print(f"📂 Cargando epochs: {args.input_file}")
epochs = mne.read_epochs(args.input_file, preload=True)
data = epochs.get_data()  # shape: (n_epochs, n_channels, n_times)

n_epochs = data.shape[0]
print(f"🔢 Total de ensayos: {n_epochs}")

# ---------------------- CÁLCULO DE MÉTRICAS POR ENSAYO ----------------------
print("🧠 Calculando métricas por ensayo...")

rms_vals = np.sqrt(np.mean(data ** 2, axis=(1, 2)))
max_amp = np.max(np.abs(data), axis=(1, 2))
std_vals = np.std(data, axis=(1, 2))

# Z-score del RMS
rms_mean = np.mean(rms_vals)
rms_std = np.std(rms_vals)
z_scores = (rms_vals - rms_mean) / rms_std

# Etiqueta "bad" si z_score > threshold
quality = ['bad' if abs(z) > args.z_thresh else 'good' for z in z_scores]

# ---------------------- EXPORTACIÓN A CSV ----------------------
subject = os.path.basename(args.input_file).split("_")[0]
out_csv = os.path.join(args.output_dir, f"{subject}_epoch_quality.csv")

df = pd.DataFrame({
    "epoch": np.arange(n_epochs),
    "rms": rms_vals,
    "max_amplitude": max_amp,
    "std": std_vals,
    "zscore_rms": z_scores,
    "epoch_quality": quality
})
df.to_csv(out_csv, index=False)
print(f"📝 Guardado CSV: {out_csv}")

# ---------------------- GENERACIÓN DE HEATMAP ----------------------
print("📊 Generando heatmap...")
fig, ax = plt.subplots(figsize=(8, 5))
im = ax.imshow(np.stack([rms_vals, max_amp, std_vals, z_scores]).T,
               aspect='auto', cmap='viridis', interpolation='nearest')

ax.set_title(f"Calidad de ensayos: {subject}")
ax.set_xlabel("Métrica")
ax.set_ylabel("Epoch")
ax.set_xticks([0, 1, 2, 3])
ax.set_xticklabels(["RMS", "MaxAmp", "STD", "Z(RMS)"])

# Marcar líneas rojas donde hay ensayos malos
for i, q in enumerate(quality):
    if q == 'bad':
        ax.axhline(i, color='red', linewidth=0.4, alpha=0.7)

plt.colorbar(im, ax=ax, shrink=0.8, label='Valor')

out_img = os.path.join(args.output_dir, f"{subject}_epoch_quality_heatmap.png")
plt.tight_layout()
plt.savefig(out_img, dpi=150)
plt.close()
print(f"🖼️ Guardado heatmap: {out_img}")

print("✅ Evaluación completada.")
