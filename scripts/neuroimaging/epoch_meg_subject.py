#!/usr/bin/env python3
# ============================================================
# Script: epoch_meg_subject.py (versión definitiva)
#
# OBJETIVO:
# Epocar señal MEG preprocesada en torno a eventos definidos,
# conservando todos los ensayos y delegando la validación de calidad
# a un análisis posterior (evaluate_epochs_quality.py).
# ============================================================

import os
import argparse
import mne
import pandas as pd
import numpy as np

mne.set_log_level('WARNING')

# ---------------------- ARGUMENTOS DE ENTRADA ----------------------
parser = argparse.ArgumentParser(description="Epocar datos MEG con soporte BIDS y multiclase")
parser.add_argument('--input_file', required=True, help='Archivo .fif de entrada (preprocesado)')
parser.add_argument('--output_file', required=True, help='Ruta del archivo de salida .epo.fif')
parser.add_argument('--annotations_file', help='(Opcional) archivo .tsv o .csv con anotaciones BIDS')
parser.add_argument('--behavior_file', help='(Opcional) archivo CSV con variables por trial')
parser.add_argument('--tmin', type=float, default=-0.2, help='Inicio del epoch (ej: -0.2)')
parser.add_argument('--tmax', type=float, default=0.8, help='Fin del epoch (ej: +0.8)')
parser.add_argument('--stim_channel', default='STI 014', help='Canal de eventos si no hay anotaciones')
parser.add_argument('--export_events', action='store_true', help='Si se activa, guarda los eventos en CSV')
args = parser.parse_args()

# ---------------------- CARGA DEL ARCHIVO .fif ----------------------
print("📂 Cargando archivo:", args.input_file)
raw = mne.io.read_raw_fif(args.input_file, preload=True)

# ---------------------- DETECCIÓN DE EVENTOS ----------------------
print("🔍 Detectando eventos para epocado...")

if args.annotations_file:
    print(f"📂 Usando anotaciones externas: {args.annotations_file}")
    if not os.path.exists(args.annotations_file):
        raise FileNotFoundError(f"❌ No se encontró: {args.annotations_file}")
    
    ext = os.path.splitext(args.annotations_file)[-1]
    if ext == '.tsv':
        df_ann = pd.read_table(args.annotations_file)
    elif ext == '.csv':
        df_ann = pd.read_csv(args.annotations_file)
    else:
        raise ValueError("❌ Anotaciones externas deben ser .csv o .tsv")

    required_cols = {'onset', 'duration'}
    if not required_cols.issubset(df_ann.columns):
        raise ValueError("❌ Faltan columnas requeridas: onset y duration")

    if 'description' in df_ann.columns:
        desc_col = 'description'
    elif 'trial_type' in df_ann.columns:
        desc_col = 'trial_type'
    elif 'value' in df_ann.columns:
        desc_col = 'value'
    else:
        raise ValueError("❌ No se encontró columna de descripción de evento")

    annotations = mne.Annotations(
        onset=df_ann['onset'].values,
        duration=df_ann['duration'].values,
        description=df_ann[desc_col].astype(str).values
    )
    raw.set_annotations(annotations)
    events, event_id = mne.events_from_annotations(raw)
    print(f"✅ Eventos cargados desde '{desc_col}':", event_id)

else:
    if 'STI 014' in raw.ch_names:
        print("✅ Canal STI 014 detectado. Extrayendo eventos...")
        events = mne.find_events(raw, stim_channel=args.stim_channel, verbose=True)
        unique_codes = sorted(set(events[:, 2]))
        event_id = {f"event_{code}": code for code in unique_codes}
        print(f"📌 Mapeo automático:", event_id)
    elif raw.annotations and len(raw.annotations) > 0:
        print("⚠️ Sin canal STI 014. Usando anotaciones internas...")
        events, event_id = mne.events_from_annotations(raw)
        print(f"✅ Eventos detectados:", event_id)
    else:
        raise RuntimeError(f"❌ No se encontraron eventos ni anotaciones en {args.input_file}")

print(f"🔢 Total de eventos detectados: {len(events)}")

# ---------------------- EXPORTACIÓN OPCIONAL DE EVENTOS ----------------------
if args.export_events:
    event_df = pd.DataFrame(events, columns=["sample", "prev", "event_code"])
    out_csv = args.output_file.replace(".epo.fif", "_events.csv")
    event_df.to_csv(out_csv, index=False)
    print(f"📝 Eventos exportados en: {out_csv}")

# ---------------------- FILTRADO DE EVENTOS EN RANGO ----------------------
sfreq = raw.info['sfreq']
start_limit = 0 - args.tmin
end_limit = raw.times[-1] - args.tmax

filtered_events = []
for ev in events:
    ev_time = ev[0] / sfreq
    if start_limit <= ev_time <= end_limit:
        filtered_events.append(ev)

filtered_events = np.array(filtered_events)
print(f"🔍 Eventos dentro de rango válido: {len(filtered_events)} / {len(events)}")

if len(filtered_events) == 0:
    raise RuntimeError("❌ No hay eventos dentro del rango permitido. Revisa tmin/tmax o duración del raw.")

# ---------------------- CREACIÓN DE EPOCHS ----------------------
epochs = mne.Epochs(
    raw, filtered_events, event_id,
    tmin=args.tmin, tmax=args.tmax,
    baseline=(args.tmin, 0),
    reject=None, preload=True, detrend=1, verbose=False
)

# ---------------------- METADATOS POR TRIAL ----------------------
if args.behavior_file:
    df_behavior = pd.read_csv(args.behavior_file)
    if "onset" in df_behavior.columns:
        print("🔍 Ordenando metadatos por onset...")
        df_behavior = df_behavior.sort_values("onset").reset_index(drop=True)
    if len(df_behavior) != len(epochs):
        print(f"⚠️ Mismatch: {len(df_behavior)} trials vs {len(epochs)} epochs. Se ignoran metadatos.")
    else:
        epochs.metadata = df_behavior
        print("✅ Metadatos por trial integrados.")

# ---------------------- GUARDADO ----------------------
epochs.save(args.output_file, overwrite=True)
print(f"✅ Epochs guardados en: {args.output_file}")

# ---------------------- LOG DE RESUMEN ----------------------
summary = {
    "subject": os.path.basename(args.input_file).split("_")[0],
    "n_epochs": len(epochs),
    "n_events": len(events),
    "n_dropped": len(events) - len(filtered_events),
    "event_types": ";".join(event_id.keys())
}

log_path = os.path.join(os.path.dirname(args.output_file), "log_epoching.csv")
df_summary = pd.DataFrame([summary])
df_summary.to_csv(log_path, mode='a', header=not os.path.exists(log_path), index=False)
print(f"📊 Resumen del epocado: {summary}")
