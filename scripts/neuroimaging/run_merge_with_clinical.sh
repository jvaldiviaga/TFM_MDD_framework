#!/bin/bash
# ============================================================
# Script: run_merge_with_clinical.sh
# Objetivo: Fusionar features_all_subjects.csv con Excel clínico real
# ============================================================

# Ruta al archivo CSV combinado con features
features_csv=~/TFM_MDD/results/ds005356/features/features_all_subjects.csv

# Archivo Excel clínico real
clinical_excel=~/TFM_MDD/data/neuroimaging/ds005356/Code/MEG\ MDD\ IDs\ and\ Quex.xlsx

# Archivo final de salida
output_file=~/TFM_MDD/results/ds005356/features/features_combined_with_clinical.csv

# Script de fusión (ahora actualizado con el nuevo comportamiento)
script_path=~/TFM_MDD/scripts/neuroimaging/merge_features_with_clinical.py

# Ejecutar fusión
python3 "$script_path" \
  --features_csv "$features_csv" \
  --clinical_excel "$clinical_excel" \
  --output_file "$output_file" \
  --id_column URSI \
  --id_prefix M871 \
  --id_digits 5
