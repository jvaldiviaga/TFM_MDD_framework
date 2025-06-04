#!/bin/bash
# ============================================================
# Script: run_all_extract_features.sh (versión final)
# Objetivo: Extraer features por sujeto y combinarlos en un único CSV
# ============================================================

# Ruta a los archivos .epo_filtered.fif
input_dir=~/TFM_MDD/results/ds005356/epochs

# Carpeta de salida
output_dir=~/TFM_MDD/results/ds005356/features
mkdir -p "$output_dir"

# Script de extracción
script_path=~/TFM_MDD/scripts/neuroimaging/extract_features_evoked.py

# Eliminar CSVs anteriores (opcional)
rm -f "$output_dir"/*_features.csv

# Bucle por sujetos
for epo_file in "$input_dir"/sub-*_epo_filtered.fif; do
  subj=$(basename "$epo_file" | cut -d'_' -f1)
  output_csv="$output_dir/${subj}_features.csv"

  echo "🔍 Procesando $subj → $output_csv"
  python3 "$script_path" \
    --input_file "$epo_file" \
    --output_file "$output_csv"
  echo "✅ $subj completado"
done

# Fusión de todos los CSVs individuales
echo "📊 Combinando todos los archivos *_features.csv"
combined_csv="$output_dir/features_all_subjects.csv"
head -n 1 $(ls "$output_dir"/*_features.csv | head -n 1) > "$combined_csv"
tail -n +2 -q "$output_dir"/*_features.csv >> "$combined_csv"
echo "✅ Combinado guardado en: $combined_csv"
