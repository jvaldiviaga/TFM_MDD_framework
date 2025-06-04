#!/usr/bin/env python3
# ============================================================
# Script: run_all_filtering.py
#
# OBJETIVO:
# Evaluar y filtrar ensayos epocados (.epo.fif) para todos los sujetos
# usando los resultados previos de calidad generados por evaluate_epochs_quality.py
#
# Este script guía al usuario paso a paso:
#  - Carga cada archivo .epo.fif
#  - Carga su CSV de calidad (por ensayo)
#  - Muestra un resumen y una imagen opcional
#  - Permite guardar un nuevo archivo solo con los "good"
#
# ➕ Pensado para usuarios sin experiencia informática
# ============================================================

import os
from filter_epochs_by_quality import filter_epochs_by_quality

# ---------------------- CONFIGURACIÓN ----------------------
# Ruta donde están los archivos epocados y los resultados de calidad
epochs_dir = os.path.expanduser("~/TFM_MDD/results/ds005356/epochs")

# Buscar todos los archivos de tipo *_epo.fif
files = sorted(f for f in os.listdir(epochs_dir) if f.endswith("_epo.fif"))
total = len(files)

# Contadores para resumen final
n_total = 0
n_good_saved = 0
n_skipped = 0

# ---------------------- LOOP POR SUJETOS ----------------------
for i, fname in enumerate(files, start=1):
    subj_id = fname.split("_")[0]
    epo_path = os.path.join(epochs_dir, fname)
    csv_path = os.path.join(epochs_dir, f"{subj_id}_epoch_quality.csv")
    heatmap_path = os.path.join(epochs_dir, f"{subj_id}_epoch_quality_heatmap.png")

    print("="*60)
    print(f"[{i}/{total}] 🔎 Evaluando calidad para: {subj_id}")
    print("- Archivo de ensayos:      ", epo_path)
    print("- Resultados de calidad:   ", csv_path)
    print("- Imagen de calidad:       ", heatmap_path if os.path.exists(heatmap_path) else "❌ No disponible")

    if not os.path.exists(csv_path):
        print(f"⚠️  No se encontró el archivo CSV de calidad para {subj_id}. Se omite.")
        n_skipped += 1
        continue

    n_total += 1
    filter_epochs_by_quality(epo_path, csv_path, heatmap_path)

    # Verifica si se guardó el archivo *_epo_good.fif
    filtered_path = epo_path.replace("_epo.fif", "_epo_good.fif")
    if os.path.exists(filtered_path):
        n_good_saved += 1

# ---------------------- RESUMEN FINAL ----------------------
print("\n" + "="*60)
print("📊 RESUMEN FINAL:")
print(f"🧠 Sujetos procesados:     {n_total}")
print(f"✅ Archivos guardados con solo 'good': {n_good_saved}")
print(f"⏭️ Sujetos omitidos (sin CSV): {n_skipped}")
print("🎯 Filtrado de calidad finalizado.\n")
