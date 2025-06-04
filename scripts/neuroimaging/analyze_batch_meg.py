#!/usr/bin/env python3
# ============================================================
# Script: analyze_batch_meg.py
#
# OBJETIVO:
# Evaluar el efecto batch tras el preprocesamiento MEG, utilizando variables técnicas
# como número de proyecciones, duración de grabación y posición de cabeza.
# Se integran metadatos clínicos para visualizar posibles sesgos sistemáticos 
# (e.g., por escáner) que puedan afectar a la comparación entre grupos clínicos.
#
# Este análisis es disruptivo porque tradicionalmente se realiza solo QC previo, 
# y rara vez se reevalúa si el batch persiste tras el preprocesamiento. Aquí se 
# propone como paso obligatorio antes de cualquier análisis funcional.
# ============================================================

import os
import mne
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import argparse

# =================== PARÁMETROS DE ENTRADA ===================
parser = argparse.ArgumentParser(description="Evaluar efecto batch en datos MEG preprocesados")
parser.add_argument('--data_dir', required=True, help='Directorio con archivos *_meg_clean.fif')
parser.add_argument('--metadata_file', required=True, help='Archivo Excel con URSI, Group, Scanner, etc.')
parser.add_argument('--output_dir', required=True, help='Ruta de salida para resultados del análisis batch')
args = parser.parse_args()

# Expandimos rutas, creamos carpeta de salida
data_dir = os.path.expanduser(args.data_dir)
output_dir = os.path.expanduser(args.output_dir)
os.makedirs(output_dir, exist_ok=True)

# =================== CARGA DE METADATOS ===================
# Leemos Excel y generamos columna subject compatible con nombres *_meg_clean.fif
df_meta = pd.read_excel(args.metadata_file)
df_meta.columns = df_meta.columns.str.strip().str.lower()
df_meta['subject'] = "M871" + df_meta['ursi'].astype(str).str.zfill(5)
df_meta['subject'] = df_meta['subject'].str.upper().str.strip()

# =================== EXTRACCIÓN TÉCNICA DE CADA SUJETO ===================
info_list = []

for fname in sorted(os.listdir(data_dir)):
    if fname.endswith("_meg_clean.fif"):
        path = os.path.join(data_dir, fname)
        subj = fname.split("_")[0].upper().strip()

        try:
            raw = mne.io.read_raw_fif(path, preload=False, verbose="ERROR")
            meas_date = raw.info['meas_date']
            n_projs = len(raw.info['projs'])
            duration = raw.n_times / raw.info['sfreq']
            proj_types = ','.join(set(p['desc'] for p in raw.info['projs'])) if raw.info['projs'] else 'None'

            # Extracción de posición de cabeza (head_x, head_y, head_z)
            if raw.info.get("dev_head_t") is not None:
                head_pos = raw.info["dev_head_t"]["trans"][:3, 3]
            else:
                head_pos = [np.nan, np.nan, np.nan]

            info_list.append({
                "subject": subj,
                "meas_date": str(meas_date.date()) if meas_date else "unknown",
                "n_projs": n_projs,
                "proj_types": proj_types,
                "duration_sec": round(duration, 2),
                "head_x": head_pos[0],
                "head_y": head_pos[1],
                "head_z": head_pos[2]
            })

        except Exception as e:
            print(f"❌ Error procesando {subj}: {e}")

df_info = pd.DataFrame(info_list)

# =================== UNIÓN CON METADATOS CLÍNICOS ===================
df = df_info.merge(df_meta, how='left', on='subject')

# Normalización de variables categóricas
df['scanner'] = df['scanner'].astype(str).str.strip()
df['group'] = df['group'].astype(str).str.strip()

# =================== PCA SOBRE VARIABLES TÉCNICAS ===================
# Seleccionamos las variables técnicas más relevantes
features = ['n_projs', 'duration_sec', 'head_x', 'head_y', 'head_z']
X = StandardScaler().fit_transform(df[features])  # Estandarización

# Reducción de dimensionalidad
pca = PCA(n_components=2)
components = pca.fit_transform(X)
df[['PC1', 'PC2']] = components

# =================== VISUALIZACIÓN DEL PCA ===================
sns.set(style="whitegrid", font_scale=1.1)
plt.figure(figsize=(10, 6))
sns.scatterplot(
    data=df,
    x="PC1", y="PC2",
    hue="group",         # Color por grupo clínico
    style="scanner",     # Forma por escáner
    s=100
)
plt.title("PCA de características técnicas: Grupo vs Scanner")
plt.xlabel("PC1 ({}% var)".format(round(pca.explained_variance_ratio_[0]*100, 1)))
plt.ylabel("PC2 ({}% var)".format(round(pca.explained_variance_ratio_[1]*100, 1)))
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "pca_batch_effect_group_scanner.png"), dpi=150)
plt.close()

# =================== DENDROGRAMA JERÁRQUICO ===================
from scipy.cluster.hierarchy import linkage, dendrogram

Z = linkage(X, method='ward')
labels = df['subject'].tolist()

plt.figure(figsize=(12, 6))
dendrogram(
    Z,
    labels=labels,
    leaf_rotation=90,
    leaf_font_size=9,
    color_threshold=0
)
plt.title("🔍 Dendrograma jerárquico por características técnicas")
plt.xlabel("Sujeto")
plt.ylabel("Distancia (Ward linkage)")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "dendrogram_batch_effect.png"), dpi=150)
plt.close()

# =================== TEST MANOVA POR SCANNER ===================
from statsmodels.multivariate.manova import MANOVA

try:
    df_manova = df.dropna(subset=features + ['scanner']).copy()
    formula = ' + '.join(features) + ' ~ scanner'
    manova = MANOVA.from_formula(formula, data=df_manova)
    result = manova.mv_test()
    print("📊 Resultado MANOVA por scanner:\n", result)
except Exception as e:
    print(f"⚠️ No se pudo calcular MANOVA: {e}")

# =================== EXPORTACIÓN ===================
df.to_csv(os.path.join(output_dir, "batch_analysis_summary.csv"), index=False)

# =================== RESUMEN FINAL ===================
print("✅ Análisis de batch completado.")
print("📁 Resultados guardados en:", output_dir)

# =================== LIMITACIONES Y MEJORAS FUTURAS ===================
# - Este script solo evalúa variables técnicas; no analiza señal funcional.
# - No aplica corrección del batch, solo lo detecta.
# - Requiere dev_head_t para head_pos; si falta, los valores serán NaN.
# - Se puede extender con correlaciones a ERF, STC o conectividad funcional.
# - Futuro: integrar neuroCombat para armonización multi-sujeto.
