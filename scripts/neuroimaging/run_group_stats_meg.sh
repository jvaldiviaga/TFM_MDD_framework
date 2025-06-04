#!/bin/bash
# ============================================================
# Script: run_group_stats_meg.sh
# Objetivo: Ejecutar análisis lineal y LASSO para MEG funcional
# Genera panel robusto para validación cruzada con transcriptómica
# ============================================================

INPUT=~/TFM_MDD/results/ds005356/features/features_combined_with_clinical.csv
OUTPUT_DIR=~/TFM_MDD/results/ds005356/stats

echo "🚀 Ejecutando análisis lineal multivariado (MDD vs CTL)..."
python3 ~/TFM_MDD/scripts/neuroimaging/group_stats_meg.py \
  --input_csv "$INPUT" \
  --output_dir "$OUTPUT_DIR" \
  --model linear

echo "🚀 Ejecutando selección funcional con LASSO..."
python3 ~/TFM_MDD/scripts/neuroimaging/group_stats_meg.py \
  --input_csv "$INPUT" \
  --output_dir "$OUTPUT_DIR" \
  --model lasso

echo "✅ Análisis funcional completo terminado. Revisa:"
echo "📄 $OUTPUT_DIR/regression_results_linear.csv"
echo "📄 $OUTPUT_DIR/lasso_selected_features.csv (panel funcional para validación cruzada)"
echo "📄 $OUTPUT_DIR/vars_utilizadas.txt (trazabilidad)"
