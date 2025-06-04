#!/usr/bin/env python3
# ============================================================
# Script: combine_features_safely.py
#
# Objetivo:
# Unir todos los *_features.csv por sujeto en un solo archivo seguro y robusto.
# ============================================================

import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Unir todos los CSVs de features individuales en uno solo")
parser.add_argument('--input_dir', required=True, help='Directorio con los *_features.csv')
parser.add_argument('--output_file', required=True, help='Ruta de salida del archivo combinado')
args = parser.parse_args()

combined = []
for fname in sorted(os.listdir(args.input_dir)):
    if fname.endswith("_features.csv"):
        fpath = os.path.join(args.input_dir, fname)
        df = pd.read_csv(fpath)
        combined.append(df)

if combined:
    df_all = pd.concat(combined, axis=0, ignore_index=True)
    df_all.to_csv(args.output_file, index=False)
    print(f"✅ Archivo combinado guardado en: {args.output_file}")
else:
    print("❌ No se encontraron archivos *_features.csv para combinar.")
