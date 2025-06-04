#!/usr/bin/env python3
# ============================================================
# Script: group_stats_meg.py
# Objetivo: Comparar MDD vs CTL con regresión lineal o LASSO.
# Biomarcadores: P300, latency, slope, GFP, RWP.
# ============================================================

import argparse
import pandas as pd
import numpy as np
import os
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.linear_model import LogisticRegressionCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer

# ---------------------- ARGUMENTOS ----------------------
parser = argparse.ArgumentParser(description="Análisis funcional MEG clínico")
parser.add_argument('--input_csv', required=True)
parser.add_argument('--output_dir', required=True)
parser.add_argument('--model', choices=['linear', 'lasso'], default='linear')
parser.add_argument('--group_column', default='Group')
parser.add_argument('--control_label', default='CTL')
parser.add_argument('--patient_label', default='MDD')
args = parser.parse_args()

# ---------------------- CARGA ----------------------
os.makedirs(args.output_dir, exist_ok=True)
df = pd.read_csv(args.input_csv)
df = df[df[args.group_column].isin([args.control_label, args.patient_label])].copy()
df['group_binary'] = df[args.group_column].map({args.control_label: 0, args.patient_label: 1})

# ---------------------- FUNCIONES ----------------------
def traducir_feature(feature_name):
    partes = feature_name.split("_")
    if len(partes) < 3:
        return feature_name
    cond = partes[0].replace("cond", "").zfill(2)
    canal = partes[-1]
    tipo = "GFP approx."
    if "slope" in feature_name.lower(): tipo = "pendiente de pico"
    elif "latency" in feature_name.lower(): tipo = "latencia del P300"
    elif "rwp" in feature_name.lower(): tipo = "duración RWP"
    elif "p300_amp" in feature_name.lower(): tipo = "amplitud del P300"
    elif "gfp" in feature_name.lower(): tipo = "GFP"
    return f"cue-{cond} – {canal} ({canal[:3]}) – {tipo}"

# ---------------------- VARIABLES ----------------------
functional_vars = [c for c in df.columns if any(k in c.lower() for k in [
    'p300_amp', 'latency', 'rwp', 'slope', 'gfp'
]) and c.startswith("cond_")]

covariates = [
    'Age', 'Sex', 'Scanner',
    'MASQ_Anxious_Arousal', 'MASQ_Anhedonic_Depression',
    'DARS_TOTAL_SCALE', 'TEPS_Total', 'PC1'
]
covariates = [c for c in covariates if c in df.columns]

with open(os.path.join(args.output_dir, 'vars_utilizadas.txt'), 'w') as f:
    f.write("VARIABLES FUNCIONALES:\n" + "\n".join(functional_vars) + "\n\n")
    f.write("COVARIABLES CLÍNICAS:\n" + "\n".join(covariates))

# ---------------------- MODO LINEAL ----------------------
if args.model == 'linear':
    results = []
    for var in functional_vars:
        if df[var].isnull().any(): continue
        formula = f"{var} ~ group_binary + " + " + ".join(covariates)
        model = smf.ols(formula, data=df).fit()
        results.append({
            'Feature': var,
            'Coef_Group_MDD': model.params.get('group_binary', np.nan),
            'P_value': model.pvalues.get('group_binary', np.nan),
            'R2': model.rsquared
        })

    df_results = pd.DataFrame(results).sort_values('P_value')
    df_results['Interpretation'] = df_results['Feature'].apply(traducir_feature)

    print("\n🧠 TOP 10 RESULTADOS LINEALES:")
    for _, row in df_results.head(10).iterrows():
        print(f"→ {row['Interpretation']}: Coef={row['Coef_Group_MDD']:.3f}, P={row['P_value']:.4f}, R²={row['R2']:.3f}")
    if len(df_results) > 10:
        print(f"... ({len(df_results)-10} resultados adicionales ocultos; revisa el CSV completo)")

    print("\n🧠 RESULTADOS LINEALES (P < 0.1):")
    for _, row in df_results[df_results['P_value'] < 0.1].head(50).iterrows():
        desc = traducir_feature(row['Feature'])
        print(f"{desc}: Coef={row['Coef_Group_MDD']:.4f}, P={row['P_value']:.4f}, R²={row['R2']:.3f}")

    df_results.to_csv(os.path.join(args.output_dir, 'regression_results_linear.csv'), index=False)
    print(f"✅ Resultados lineales guardados en: {args.output_dir}/regression_results_linear.csv")

    # ---------------------- VISUALIZACIÓN ----------------------
    if not df_results.empty:
        feature_to_plot = df_results.iloc[0]['Feature']
        plt.figure(figsize=(6, 4))
        sns.boxplot(data=df, x=args.group_column, y=feature_to_plot, palette="Set2")
        sns.stripplot(data=df, x=args.group_column, y=feature_to_plot, color='black', alpha=0.3)
        plt.title(traducir_feature(Path(feature_to_plot).stem))
        plot_name = feature_to_plot.replace(" ", "_") + "_boxplot.png"
        plot_path = os.path.join(args.output_dir, plot_name)
        plt.tight_layout()
        plt.savefig(plot_path)
        plt.close()
        print(f"📊 Visualización guardada en: {plot_path}")

# ---------------------- MODO LASSO ----------------------
elif args.model == 'lasso':
    X = df[functional_vars].copy()
    y = df['group_binary']
    if X.empty:
        print("❌ Matriz vacía. No se puede ejecutar LASSO.")
        exit(1)

    pipeline = Pipeline([
        ('imputer', SimpleImputer(strategy='mean')),
        ('scaler', StandardScaler()),
        ('lasso', LogisticRegressionCV(
            penalty='l1',
            solver='saga',
            Cs=np.logspace(-4, 2, 20),
            cv=5,
            max_iter=10000,
            scoring='accuracy',
            n_jobs=-1
        ))
    ])

    pipeline.fit(X, y)
    coefs = pipeline.named_steps['lasso'].coef_[0]
    selected = [{'Feature': f, 'LASSO_Coef': round(c, 4)} for f, c in zip(X.columns, coefs) if c != 0]

    if selected:
        df_selected = pd.DataFrame(selected)
        df_selected['Interpretation'] = df_selected['Feature'].apply(traducir_feature)
        df_selected = df_selected.sort_values('LASSO_Coef', key=lambda x: abs(x), ascending=False)
        df_selected.to_csv(os.path.join(args.output_dir, 'lasso_selected_features.csv'), index=False)
        print(f"✅ Biomarcadores seleccionados por LASSO guardados en: {args.output_dir}/lasso_selected_features.csv")
        print("\n🔎 FEATURES SELECCIONADOS POR LASSO:")
        for _, row in df_selected.iterrows():
            print(f"→ {row['Interpretation']}: Coef={row['LASSO_Coef']}")
    else:
        print("⚠️ LASSO no seleccionó ningún feature.")
