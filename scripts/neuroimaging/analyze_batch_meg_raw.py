#!/usr/bin/env python3
# ============================================================
# Script: analyze_batch_meg_raw.py
#
# OBJETIVO:
# Detectar posibles efectos batch (scanner, duración, fecha de adquisición)
# en los archivos MEG crudos (.fif) ANTES del preprocesado.
# Esto permite anticipar problemas sistemáticos y adaptar decisiones técnicas
# como el filtrado, corrección por ICA o exclusión.
#
# === FUNCIONES IMPLEMENTADAS ===
# ✅ Extracción de metadatos técnicos desde raw (duración, fecha, n_projs)
# ✅ Fusión con metadatos clínicos (scanner, grupo)
# ✅ PCA con visualización por grupo y por scanner
# ✅ Dendrograma jerárquico de similitud técnica
# ✅ Cálculo de silhouette score para detección objetiva de batch
# ✅ Exportación a batch_raw_summary.csv + PNGs
#
# === DISRUPTIVO RESPECTO A PIPELINES COMUNES ===
# ✅ Evalúa el efecto batch antes del preprocesamiento
# ✅ Ayuda a decidir si conviene armonizar o modificar la estrategia técnica
# ✅ Reduce gasto computacional en sujetos mal registrados
#
# === LIMITACIONES (CONSCIENTES) ===
# ❌ No evalúa diferencias funcionales (solo técnicas)
# ❌ No armoniza automáticamente (p.ej., neuroComBat)
# ✅ Es un bloque exploratorio previo → guía decisiones clínicas posteriores
# ============================================================

# Evaluamos el efecto batch tanto antes como después del preprocesado:
# - Antes: para detectar sesgos técnicos estructurales (scanner, duración, n_projs) que puedan comprometer el análisis,
# - Después: para verificar si dichos sesgos persisten tras la limpieza (ICA, filtros) y decidir si aplicar corrección (e.g., neuroCombat).

# ======================== LIBRERÍAS ========================
import mne
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import linkage, dendrogram
import argparse
import os
import seaborn as sns

# ===================== PARÁMETROS DE ENTRADA =====================
parser = argparse.ArgumentParser(description="Evaluar efecto batch en MEG crudo.")
parser.add_argument('--data_dir', type=str, required=True,
                    help='Directorio con archivos *_meg.fif (crudos)')
parser.add_argument('--metadata_file', type=str, required=True,
                    help='Archivo Excel con columnas: URSI, Group, Scanner')
parser.add_argument('--output_dir', type=str, required=True,
                    help='Directorio de salida para resultados')
args = parser.parse_args()

data_dir = os.path.expanduser(args.data_dir)
output_dir = os.path.expanduser(args.output_dir)
os.makedirs(output_dir, exist_ok=True)

# ===================== CARGA Y FUSIÓN DE METADATOS CLÍNICOS =====================
df_meta = pd.read_excel(args.metadata_file)
df_meta.columns = df_meta.columns.str.strip().str.lower()  # 🔧 Convertimos todas a minúscula

# Construir ID compatible con archivos .fif (usa 'ursi' en minúscula ahora)
df_meta['subject'] = 'M871' + df_meta['ursi'].astype(str).str.zfill(5)

# Renombrar campos si fuera necesario (opcional en este caso)
# df_meta = df_meta.rename(columns={'group': 'group', 'scanner': 'scanner'})

# Mostramos columnas detectadas (útil para debug)
print("🧩 Columnas detectadas en df_meta:", df_meta.columns.tolist())

# Validación de duplicados
if df_meta['subject'].duplicated().any():
    raise ValueError("❌ Hay IDs duplicados en df_meta['subject']. Revisa el Excel.")

# Verificamos columnas obligatorias
required_cols = ['subject', 'group', 'scanner']
missing = [col for col in required_cols if col not in df_meta.columns]
if missing:
    raise ValueError(f"❌ Faltan columnas clave en el Excel: {missing}")

# ===================== NOTA SOBRE GENERALIZACIÓN =====================
# Este bloque está adaptado al dataset AIM/Reward.
# Si se utiliza otro Excel, es necesario revisar:
# - Los nombres de columnas
# - El encoding invisible (espacios, tildes, mayúsculas)
# - Que 'group' y 'scanner' no contengan NaN ni codificaciones no estándar


# ===================== EXTRACCIÓN DE INFORMACIÓN TÉCNICA =====================
info_list = []
for fname in sorted(os.listdir(data_dir)):
    if fname.startswith("sub-") and fname.endswith("_meg.fif"):
        print(f"📂 Leyendo archivo: {fname}", flush=True)  # Muestra nombre antes de intentar abrir

        path = os.path.join(data_dir, fname)
        subj = os.path.basename(fname).split("_")[0].replace("sub-", "").strip()

        try:
            print(f"🔄 Extrayendo info técnica de {subj}...", flush=True)
            raw = mne.io.read_raw_fif(path, preload=False, verbose="ERROR")

            meas_date = raw.info['meas_date']
            duration = raw.n_times / raw.info['sfreq']
            sfreq = raw.info['sfreq']
            n_projs = len(raw.info['projs'])
            proj_types = ','.join(set(p['desc'] for p in raw.info['projs'])) if raw.info['projs'] else 'None'

            info_list.append({
                "subject": subj,
                "meas_date": str(meas_date.date()) if meas_date else "unknown",
                "duration_sec": round(duration, 2),
                "sfreq": sfreq,
                "n_projs": n_projs,
                "proj_types": proj_types
            })

            print(f"✅ {subj} procesado", flush=True)

        except Exception as e:
            print(f"❌ Error leyendo {fname}: {e}", flush=True)


df_info = pd.DataFrame(info_list)
print("✅ IDs en df_info (de archivos):", df_info['subject'].unique()[:5])
print("✅ IDs en df_meta (desde Excel):", df_meta['subject'].unique()[:5])

# ===================== FUSIÓN Y LIMPIEZA =====================
# Homogeneizamos los IDs a mayúsculas sin espacios
df_info['subject'] = df_info['subject'].astype(str).str.upper().str.strip()
df_meta['subject'] = df_meta['subject'].astype(str).str.upper().str.strip()

# Hacemos el merge clínico-técnico por ID homogéneo
df = df_info.merge(df_meta, how='left', on='subject')

# Normalizamos columnas
df['group'] = df['group'].astype(str).str.strip()
df['scanner'] = pd.Categorical(df['scanner'])  # Categoría necesaria para silhouette_score

# Eliminamos sujetos sin metadatos clínicos válidos
df = df.dropna(subset=['group', 'scanner'])

if df.empty:
    raise ValueError("❌ No hay sujetos válidos tras el merge. Revisa IDs o columnas del Excel.")

# Diagnóstico en consola
print(f"👥 Sujetos crudos detectados: {len(df_info)}")
print(f"📋 Sujetos con metadatos válidos: {len(df)}")
print("📊 Grupos detectados:", df['group'].unique())
print("🧪 Scanners detectados:", df['scanner'].unique())

# ===================== PCA =====================
features = ['duration_sec', 'sfreq', 'n_projs']
X = StandardScaler().fit_transform(df[features])
pca = PCA(n_components=2)
df[['PC1', 'PC2']] = pca.fit_transform(X)

# PCA combinado
sns.set(style="whitegrid", font_scale=1.1)
plt.figure(figsize=(10, 6))
sns.scatterplot(data=df, x="PC1", y="PC2", hue="group", style="scanner", palette="Set2", s=100)
plt.title("PCA de sujetos MEG crudos: duración / scanner / proyecciones")
plt.xlabel(f"PC1 ({round(pca.explained_variance_ratio_[0]*100, 1)}%)")
plt.ylabel(f"PC2 ({round(pca.explained_variance_ratio_[1]*100, 1)}%)")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "pca_batch_effect_raw.png"), dpi=150)
plt.close()

# PCA por grupo
plt.figure(figsize=(10, 6))
sns.scatterplot(data=df, x="PC1", y="PC2", hue="group", palette="Set2", s=100)
plt.title("PCA coloreado por grupo clínico")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "pca_by_group.png"), dpi=150)
plt.close()

# PCA por scanner
plt.figure(figsize=(10, 6))
sns.scatterplot(data=df, x="PC1", y="PC2", hue=df['scanner'].astype(str), palette="Set1", s=100)
plt.title("PCA coloreado por scanner (batch)")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "pca_by_scanner.png"), dpi=150)
plt.close()

# ===================== DENDROGRAMA =====================
Z = linkage(X, method='ward')
plt.figure(figsize=(12, 6))
dendrogram(Z, labels=df['subject'].tolist(), leaf_rotation=90, leaf_font_size=9)
plt.title("Dendrograma jerárquico: agrupación técnica")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "dendrogram_batch_raw.png"), dpi=150)
plt.close()

# ===================== DETECCIÓN AUTOMÁTICA DE BATCH =====================
try:
    sil_score = silhouette_score(X, df['scanner'].cat.codes)
    print(f"🔍 Silhouette score por scanner (batch): {sil_score:.2f}")
    if sil_score > 0.3:
        print("⚠️ Posible efecto batch detectado.")
    else:
        print("✅ No se detecta efecto batch claro.")
except Exception as e:
    print(f"⚠️ No se pudo calcular silhouette score: {e}")

# ===================== EXPORTACIÓN =====================
df.to_csv(os.path.join(output_dir, "batch_raw_summary.csv"), index=False)
print("✅ Resultados guardados en:", output_dir)
