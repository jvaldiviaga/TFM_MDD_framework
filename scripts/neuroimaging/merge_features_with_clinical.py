#!/usr/bin/env python3
# ============================================================
# Script: merge_combined_features_with_clinical.py
#
# Objetivo:
# Fusionar un único archivo de features combinados (todos los sujetos)
# con un archivo Excel clínico, de forma robusta y sin depender de nombres fijos.
# ============================================================

import argparse
import pandas as pd

# ---------------------- ARGUMENTOS ----------------------
# Argumentos configurables por línea de comandos para mayor flexibilidad

parser = argparse.ArgumentParser(description="Fusionar features combinados con Excel clínico")

# Ruta al CSV combinado con features por sujeto
parser.add_argument('--features_csv', required=True, help='Archivo CSV con features de todos los sujetos')

# Archivo Excel con información clínica, demográfica, etc.
parser.add_argument('--clinical_excel', required=True, help='Archivo Excel con metadatos clínicos')

# Nombre del archivo CSV de salida combinado
parser.add_argument('--output_file', default="features_combined.csv", help='Ruta del archivo final de salida')

# Columna en el Excel con el identificador numérico del sujeto
parser.add_argument('--id_column', default="URSI", help='Nombre de la columna de ID en el Excel clínico')

# Prefijo común para los sujetos (ej. M871)
parser.add_argument('--id_prefix', default="M871", help='Prefijo para formar el nombre del sujeto completo')

# Número de dígitos a completar con ceros (ej. 58 → 00058)
parser.add_argument('--id_digits', type=int, default=5, help='Número de dígitos para rellenar el ID')

# Cargar argumentos del usuario
args = parser.parse_args()

# ---------------------- CARGA DE DATOS ----------------------

# Cargar CSV con features combinados de todos los sujetos
print("📂 Cargando archivo de features:", args.features_csv)
df_features = pd.read_csv(args.features_csv)

# Cargar archivo Excel clínico
print("📂 Cargando archivo clínico:", args.clinical_excel)
try:
    df_clinical = pd.read_excel(args.clinical_excel)
except Exception as e:
    raise ValueError(f"❌ Error al cargar el archivo Excel: {e}")

# ---------------------- PROCESAMIENTO DEL ID ----------------------

# Verificar que la columna de ID esté presente
if args.id_column not in df_clinical.columns:
    raise ValueError(f"❌ La columna '{args.id_column}' no existe en el Excel clínico.")

# Normalizar los IDs (rellenar con ceros + prefijo)
df_clinical[args.id_column] = df_clinical[args.id_column].astype(str).str.zfill(args.id_digits)
df_clinical['subject'] = args.id_prefix + df_clinical[args.id_column]

# Eliminar columnas innecesarias o trial-wise en el Excel
exclude_keywords = ['trial', 'reaction', 'RT', 'PE', 'correct', 'feedback', 'block', 'stim', 'condition']
cols_to_exclude = [col for col in df_clinical.columns if any(k in col.lower() for k in exclude_keywords)]
cols_to_keep = [col for col in df_clinical.columns if col not in cols_to_exclude and col != args.id_column and col != "subject"]
df_clinical = df_clinical[['subject'] + cols_to_keep]

# ---------------------- FUSIÓN ----------------------

# Combinar por columna 'subject'
df_final = df_features.merge(df_clinical, on="subject", how="left")

# Diagnóstico
n_total = len(df_final)
n_missing = df_final['Group'].isnull().sum() if 'Group' in df_final.columns else 0
print(f"✅ Fusión completada. Sujetos totales: {n_total}, con datos clínicos faltantes: {n_missing}")

# ---------------------- GUARDADO ----------------------

# Guardar archivo combinado
df_final.to_csv(args.output_file, index=False)
print(f"📄 Archivo final guardado en: {args.output_file}")


# ---------------------- JUSTIFICACIÓN TÉCNICA ----------------------
# Este script permite fusionar de forma automatizada las métricas funcionales obtenidas por sujeto (como amplitud, latencia o GFP) con las variables clínicas y demográficas del estudio.
# Está diseñado para ser reutilizable en diferentes contextos (MEG, EEG, fMRI) y distintos datasets, sin depender de nombres fijos. 
# Esto permite preparar una matriz completa de análisis donde cada fila representa un sujeto, y cada columna una métrica funcional o variable clínica, lista para análisis estadístico,
# comparaciones por grupo, regresiones clínicas o aplicación de algoritmos de machine learning.
