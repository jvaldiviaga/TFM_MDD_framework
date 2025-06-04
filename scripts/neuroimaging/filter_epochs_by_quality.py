#!/usr/bin/env python3
# ============================================================
# Script: filter_epochs_by_quality.py
#
# OBJETIVO:
# Filtrar ensayos MEG marcados como "bad" tras evaluación de calidad
# Generar nuevo archivo .epo_filtered.fif solo con ensayos válidos
# Justificar todos los umbrales aplicados para análisis clínico
# ============================================================

import os
import pandas as pd
import mne
import numpy as np

def filter_epochs_by_quality(epo_path, quality_csv_path, heatmap_path=None):
    subject = os.path.basename(epo_path).split("_")[0]
    print(f"\n🔎 Evaluando: {subject}")

    # --- Cargar epochs ---
    print(f"📂 Cargando archivo de epochs: {epo_path}")
    epochs = mne.read_epochs(epo_path, preload=True, verbose=False)
    print("🔍 Total de ensayos en epochs:", len(epochs))

    # --- Cargar métricas de calidad ---
    if not os.path.exists(quality_csv_path):
        print(f"❌ No se encontró el CSV de calidad: {quality_csv_path}")
        return

    df = pd.read_csv(quality_csv_path)
    print("🔍 Total de filas en CSV:", len(df))

    # --- Verificación de alineación ---
    if len(df) != len(epochs):
        print("❌ El número de filas del CSV no coincide con los ensayos del .fif")
        print("🔎 Verifica que el archivo CSV y el .epo.fif sean del mismo sujeto y no prefiltrados.")
        return
    else:
        print("✅ Verificación superada: el CSV y el archivo .fif están alineados (1 fila por ensayo)")


    # --- Clasificación de severidad por ensayo ---
    # UMBRALES JUSTIFICADOS:
    # - max_amplitude > 4000 → señal potencialmente artefactual, pero tolerable ("warning")
    # - max_amplitude > 8000 → artefacto severo (probable movimiento, saturación) → "bad"
    # - zscore_rms > 3       → ensayo estadísticamente aberrante (outlier) → "bad"
    # Referencias:
    # • MNE-Python Documentation
    # • Cavanagh et al. (2011) — Electrophysiology of error and conflict monitoring
    # • Luck, S. J. (2014) — An Introduction to the Event-Related Potential Technique

    df["severity"] = "good"
    df.loc[(df["max_amplitude"] > 4000) & (df["max_amplitude"] <= 8000), "severity"] = "warning"
    df.loc[(df["zscore_rms"] > 3) | (df["max_amplitude"] > 8000), "severity"] = "bad"

    n_total = len(df)
    n_good = (df["severity"] == "good").sum()
    n_warn = (df["severity"] == "warning").sum()
    n_bad = (df["severity"] == "bad").sum()
    z_bad = df[(df["severity"] == "bad") & (df["zscore_rms"] > 3)]["zscore_rms"].round(2).tolist()

    print(f"🔢 Total de ensayos: {n_total}")
    print(f"✔️ Buenos (good): {n_good}")
    print(f"⚠️ Aviso  (warning): {n_warn} → amplitud moderada")
    print(f"❌ Malos  (bad): {n_bad} → zscore > 3 o max_amplitude > 8000")
    print(f"📈 Z-scores de malos: {z_bad[:5]}{' ...' if len(z_bad) > 5 else ''}")

    # --- UMBRAL DE CALIDAD GLOBAL ---
    # Justificación: si más del 20 % de los ensayos son "bad", se considera el sujeto clínicamente no fiable
    bad_ratio = n_bad / n_total
    if bad_ratio > 0.20:
        print(f"\n❌ {subject} descartado: {bad_ratio*100:.1f}% de ensayos marcados como 'bad'")
        return

    # --- Selección y guardado de ensayos válidos ---
    # --- Selección y guardado de ensayos válidos ---
    valid_rows = df[df["severity"].isin(["good", "warning"])]
    selected_epochs = valid_rows["epoch"].astype(int).tolist()

    # Verificación rápida
    if max(selected_epochs, default=-1) >= len(epochs):
       print(f"❌ Índices fuera de rango detectados. Epochs tiene {len(epochs)} ensayos.")
       return

    
    try:
        epochs_filtered = epochs[selected_epochs]
    except Exception as e:
        print(f"❌ Error al acceder a los ensayos seleccionados con indexado directo: {e}")
        return


    out_path = epo_path.replace("_epo.fif", "_epo_filtered.fif")
    epochs_filtered.save(out_path, overwrite=True)
    print(f"✅ Archivo guardado con ensayos filtrados en: {out_path}")
