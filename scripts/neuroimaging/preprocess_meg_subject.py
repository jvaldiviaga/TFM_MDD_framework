
#!/usr/bin/env python3
# ============================================================
# Script: preprocess_meg_subject.py (versión QC-aware final)
#
# OBJETIVO:
# Preprocesar automáticamente un archivo MEG crudo (*.fif),
# aplicando sólo los filtros necesarios según el QC previo.
# Reutiliza la detección de canales malos ya realizada
# y ajusta los pasos en función de problemas reales detectados.
#
# JUSTIFICACIÓN DEL FLUJO:
# ✅ Evita redundancia: no recalcula bad channels
# ✅ Es adaptativo: aplica solo los filtros requeridos por el QC (salvo la señal de red a 60Hz que siempre la trata)
# ✅ Mejora trazabilidad: enlaza el preproc con el informe de QC
# ✅ Enfoque clínico: reduce sobreprocesamiento innecesario
#
# FUENTES:
# Basado en AIM_2_MEG_PREPROC_CARC.py + flujos MNE/Elekta + resultados de qc_meg_raw.py
# ============================================================
import os
import argparse
import mne
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ========================== PARÁMETROS DE ENTRADA ==========================
parser = argparse.ArgumentParser(description="Preprocesamiento MEG adaptativo basado en QC.")
parser.add_argument('--input_file', type=str, required=True, help='Archivo .fif crudo (ej. sub-XXXX_meg.fif)')
parser.add_argument('--output_dir', type=str, required=True, help='Directorio donde guardar el archivo filtrado')
parser.add_argument('--log_file', type=str, default="log_preproc_meg.csv", help='Archivo CSV para registrar resumen')
args = parser.parse_args()

# Construcción de rutas absolutas
input_path = os.path.expanduser(args.input_file)
output_dir = os.path.expanduser(args.output_dir)
os.makedirs(output_dir, exist_ok=True)
log_path = os.path.join(output_dir, args.log_file)

# Extracción del ID del sujeto
subject_id = os.path.basename(input_path).split("_")[0]

# ========================== FASE 0B: OMITIR SI YA EXISTE EL ARCHIVO LIMPIO ==========================
filtered_path = os.path.join(output_dir, f"{subject_id}_meg_clean.fif")
if os.path.exists(filtered_path):
    print(f"⏭️ Archivo ya existente: {filtered_path} → sujeto omitido.")
    exit(0)

# ========================== FASE 1: CARGA DEL ARCHIVO RAW ==========================
print(f"📂 Cargando archivo: {input_path}")
raw = mne.io.read_raw_fif(input_path, preload=False)
original_duration = raw.n_times / raw.info["sfreq"]

if original_duration > 600:
    print(f"⚠️ Duración original: {original_duration:.1f}s → recorte obligatorio a 600s (validado por bibliografía)")
    raw.crop(tmax=600.0)

raw.load_data()

# ========================== FASE 2: LECTURA DEL QC ==========================
# Esta fase consulta los resultados del QC previo para adaptar el preprocesamiento.
# No se descarta ningún sujeto: incluso si el QC sugiere que es no procesable, se intenta procesar y se registra todo.
# El QC se basa en el análisis espectral (PSD) realizado en `qc_meg_raw.py`, y guarda tanto variables binarias como métricas numéricas continuas.

print("📋 Leyendo informes de calidad técnica (QC)...")

# Ruta al resumen PSD por sujeto generado por el script qc_meg_raw.py
qc_path = os.path.join(os.path.dirname(output_dir), "qc_pre_filtering", "qc_post_binary_summary.csv")
qc_df = pd.read_csv(qc_path)

# Filtramos la fila correspondiente al sujeto actual
qc_row = qc_df[qc_df["subject"] == subject_id]

# Si no se encuentra, abortamos el flujo: es un fallo estructural del pipeline (no del sujeto)
if qc_row.empty:
    print(f"❌ No se encontró entrada de QC para {subject_id}. Abortando.")
    exit(1)

# ========================================================================
# EXTRACCIÓN DE VARIABLES BINARIAS DE ANOMALÍA
# Estas variables se generaron a partir del análisis visual automatizado
# de la PSD en escala de grises. Cada una representa una anomalía técnica clara:
# ========================================================================

drift_present = bool(qc_row["low_freq_drift_<1.5Hz"].values[0])         # Drift técnico o desacoplamiento
muscle_present = bool(qc_row["muscle_noise_>50Hz"].values[0])           # Ruido EMG
has_clipping = bool(qc_row["clipping_detected"].values[0])              # Saturación de señal
n_bad_pre = int(qc_row["n_bad_channels"].values[0])                     # Canales ruidosos
events_ok = bool(qc_row["stim_channel_present"].values[0]) if "stim_channel_present" in qc_row else True
# Nota: events_ok indica si el canal STI014 está presente y tiene eventos válidos.
# Su ausencia puede inutilizar el epocado posterior, por lo que se considera señal de alerta.

# ========================================================================
# DECISIÓN ADAPTATIVA: ¿HAY ANOMALÍA TÉCNICA DETECTADA?
anomaly_detected = any([
    drift_present,
    muscle_present,
    has_clipping,
    n_bad_pre > 3,
    not events_ok
])

# ========================================================================
# CARGA DE CANALES MALOS DETECTADOS POR MAXWELL
# Si existen, se anotan en raw.info["bads"] y se aplicará Maxwell más adelante.
# ========================================================================

bad_channels = []
bad_path = os.path.join(os.path.dirname(output_dir), "qc_pre_filtering", "qc_bad_channels.csv")
if os.path.exists(bad_path):
    bad_df = pd.read_csv(bad_path)
    match = bad_df[bad_df["subject"] == subject_id]
    if not match.empty:
        bad_str = match.iloc[0]["bad_channels"]
        if isinstance(bad_str, str) and bad_str.strip():
            bad_channels = bad_str.split(",")
            raw.info["bads"] = bad_channels
            print(f"⚠️ Canales malos anotados desde QC: {bad_channels}")
        else:
            print(f"ℹ️ {subject_id} no tiene canales marcados como malos.")
else:
    print("ℹ️ No se encontró archivo de canales malos (se asume ninguno).")


# ========================== FASE 3: APLICACIÓN DE MAXWELL FILTER ==========================
# Esta fase aplica corrección por artefactos externos mediante Maxwell filtering (SSS).
# Se aplica si hay canales ruidosos detectados o si existe dev_head_t.
# También detecta y guarda la matriz de posición de la cabeza (`dev_head_t`) si está presente,
# lo cual es útil para análisis posteriores (e.g., estimación de fuentes, coregistro anatómico).

print("🧲 Evaluando aplicación de Maxwell filtering...")

if bad_channels or raw.info.get('dev_head_t'):
    print("🔧 Se aplicará Maxwell filtering porque hay canales ruidosos o se detectó head_pos.")
    try:
        raw = mne.preprocessing.maxwell_filter(raw, cross_talk=None, calibration=None)
        maxwell_applied = "yes"
    except Exception as e:
        print(f"❌ Error al aplicar Maxwell filtering: {e}")
        maxwell_applied = "no"
else:
    print("ℹ️ No se aplica Maxwell: no hay canales malos ni head_pos.")
    maxwell_applied = "no"

# Comprobación de la matriz de posición de cabeza (head_pos), útil para trazabilidad y análisis morfológico
if raw.info.get('dev_head_t'):
    head_pos_info = raw.info['dev_head_t']['trans']
    head_pos_detected = "yes"
    print("📐 Se detectó información de posición de cabeza (head_pos).")
    print(f"🧭 Matriz de transformación (head_pos):\n{head_pos_info}")

    # Guardamos esta matriz como archivo de texto para revisión manual o reconstrucción posterior
    head_pos_path = os.path.join(output_dir, f"{subject_id}_head_pos_matrix.txt")
    np.savetxt(head_pos_path, head_pos_info, fmt="%.6f")
    print(f"📁 Matriz head_pos guardada en: {head_pos_path}")
else:
    print("⚠️ No se encontró dev_head_t → no se detectó posición de cabeza.")
    head_pos_detected = "no"


# ========================== FASE 4: FILTRADO ADAPTATIVO SEGÚN QC ==========================
# Esta fase aplica filtros para corregir artefactos de frecuencia específicos.
# - El notch de 60 Hz se aplica **siempre**, independientemente del QC, porque es una fuente común de ruido eléctrico.
# - El bandpass (1–40 Hz) se aplica **solo si se detectaron anomalías** como drift o ruido muscular.

print("⚙️ Ejecutando filtrado adaptativo...")

# ================== FILTRO NOTCH (60 Hz) ==================
# Se aplica SIEMPRE porque el ruido de red es sistémico (incluso si no es visible en PSD).
# Se omite la comprobación de `60Hz_peak` porque este ya no se usa como criterio de decisión.
print("⚡ Aplicando filtro notch a 60 Hz (forzado por política estándar)")
try:
    raw.notch_filter(freqs=60)
    notch_applied = "yes"
except Exception as e:
    print(f"❌ Error al aplicar notch: {e}")
    notch_applied = "no"

# ================== FILTRO BANDPASS (1–40 Hz) ==================
# Se aplica solo si el QC detectó artefactos de tipo muscular (>50 Hz), drift (<1.5 Hz), clipping o múltiples canales ruidosos.
# La variable `anomaly_detected` fue calculada previamente en FASE 2.

if anomaly_detected:
    print("🎚️ Anomalías técnicas detectadas → aplicando filtro pasa banda [1–40 Hz]")
    try:
        raw.filter(l_freq=1.0, h_freq=40.0, fir_design='firwin')
        bandpass_applied = "yes"
    except Exception as e:
        print(f"❌ Error al aplicar filtro bandpass: {e}")
        bandpass_applied = "no"
else:
    print("✅ QC no indica artefactos → se omite el filtro bandpass")
    bandpass_applied = "no"

# ========================== FASE 4B: ICA ADAPTATIVA PARA EOG / ECG ==========================
# Esta fase permite eliminar artefactos fisiológicos típicos de MEG:
# - Parpadeos (EOG)
# - Latidos cardíacos (ECG)
# La ICA se ejecuta SIEMPRE que existan canales EOG o ECG en el archivo.
# No depende del QC: se asume que si hay canales fisiológicos, deben corregirse por defecto.

# Detectamos si hay canales EOG/ECG disponibles
has_eog = any('EOG' in ch for ch in raw.ch_names)
has_ecg = any('ECG' in ch for ch in raw.ch_names)

# Inicializamos flags y contenedores para seguimiento del proceso
ica_applied = "no"
ica = None
exclude = []           # Componentes a eliminar
fallback_used = False  # Por si se usa la segunda estrategia (n fija)

if has_eog or has_ecg:
    print("🧠 Ejecutando ICA (se han detectado canales EOG o ECG)...")

    # Primer intento: ICA con n_components=0.99 (explicación del 99% de la varianza)
    try:
        ica = mne.preprocessing.ICA(n_components=0.99, method='fastica', random_state=42)
        ica.fit(raw)
        fallback_used = False  # No se usó fallback si esto funciona

        # Si el número de componentes útiles es muy bajo, se considera fallido (posible overfitting)
        if ica.n_components_ <= 1:
            raise ValueError(f"⚠️ Solo {ica.n_components_} componente explica >99% de varianza.")
    except Exception as e:
        print(f"⚠️ ICA con varianza 99% falló: {e}")
        print("🔁 Reintentando con n_components=20 como fallback (estrategia robusta para datos ruidosos)...")

        # Segundo intento: ICA con número fijo de componentes
        try:
            ica = mne.preprocessing.ICA(n_components=20, method='fastica', random_state=42)
            ica.fit(raw)
            fallback_used = True
        except Exception as e2:
            print(f"❌ ICA fallback también falló: {e2}")
            ica = None

    # Si ICA fue exitosa, buscamos artefactos en los componentes
    if ica is not None:
        if has_eog:
            eog_inds, _ = ica.find_bads_eog(raw)
            if eog_inds:
                exclude.extend(eog_inds)
                print(f"👁️ Componentes con artefactos EOG detectados: {eog_inds}")
            else:
                print("✅ No se detectaron artefactos EOG.")

        if has_ecg:
            ecg_inds, _ = ica.find_bads_ecg(raw)
            if ecg_inds:
                exclude.extend(ecg_inds)
                print(f"❤️ Componentes con artefactos ECG detectados: {ecg_inds}")
            else:
                print("✅ No se detectaron artefactos ECG.")

        # Si hay componentes a eliminar, se aplica la ICA y se guarda todo
        if exclude:
            ica.exclude = exclude
            raw = ica.apply(raw)
            ica_applied = "yes"
            print(f"✅ ICA aplicada. Componentes eliminados: {exclude}")

            # Alerta si se eliminaron demasiados (puede indicar señal pobre o ICA inestable)
            if len(exclude) > 10:
                print(f"⚠️ ¡Atención! Se eliminaron {len(exclude)} componentes ICA en {subject_id}")

            # Guardamos el objeto ICA para trazabilidad y reproducibilidad
            ica_path = os.path.join(output_dir, f"{subject_id}_ica.fif")
            try:
                ica.save(ica_path, overwrite=True)
                print(f"📁 Objeto ICA guardado en: {ica_path}")
            except Exception as e:
                print(f"⚠️ Error al guardar el objeto ICA: {e}")

            # Guardamos una imagen con los primeros componentes para revisión visual rápida
            ica_fig_path = os.path.join(output_dir, f"{subject_id}_ica_components.png")
            try:
                ica.plot_components(inst=raw, picks=range(min(20, ica.n_components_)), show=False)
                plt.savefig(ica_fig_path, dpi=150)
                plt.close()
                print(f"🖼️ Imagen de componentes ICA guardada en: {ica_fig_path}")
            except Exception as e:
                print(f"⚠️ No se pudo guardar la imagen ICA: {e}")
        else:
            print("ℹ️ ICA ejecutada, pero no fue necesario eliminar ningún componente.")
else:
    print("ℹ️ ICA omitida: no hay canales EOG o ECG disponibles.")

# ========================== FASE 4C: POST-ICA HEAR (opcional) ==========================
# Este paso adicional evalúa si aún persisten artefactos transitorios de alta varianza
# como pops, movimientos bruscos o deriva residual. Si se detectan, se loguea para revisión futura.

print("🔍 Evaluando presencia de artefactos de alta varianza post-ICA (HEAR)...")

window_size = int(raw.info['sfreq'] * 2)  # 2 segundos
data = raw.get_data()
n_channels, n_times = data.shape

# Deslizamiento por ventanas
high_var_windows = 0
for start in range(0, n_times - window_size, window_size):
    window = data[:, start:start + window_size]
    if np.var(window) > 1e-22:  # Umbral conservador, ajustable
        high_var_windows += 1

if high_var_windows > 0:
    print(f"⚠️ HEAR: Se detectaron {high_var_windows} ventanas con alta varianza post-ICA.")
else:
    print("✅ HEAR: No se detectaron artefactos de alta varianza post-ICA.")

# ========================== FASE 5: GUARDADO DEL ARCHIVO LIMPIO ==========================
# Guardamos el archivo ya preprocesado (tras ICA, filtros, etc.) en el directorio indicado.
# Siempre se sobrescribe cualquier archivo anterior con el mismo nombre.

filtered_path = os.path.join(output_dir, f"{subject_id}_meg_clean.fif")

print(f"💾 Guardando archivo preprocesado en: {filtered_path}")
try:
    raw.save(filtered_path, overwrite=True, buffer_size_sec=60.0)
except Exception as e:
    print(f"❌ ERROR al guardar el archivo: {e}")

# Verificamos si efectivamente se guardó el archivo (por seguridad)
if os.path.exists(filtered_path):
    print(f"✅ Archivo guardado con éxito: {filtered_path}")
else:
    print(f"❌ ERROR: El archivo no se guardó. Revisa permisos, espacio en disco o errores en MNE.")

# ========================== FASE 6: LOG CLÍNICO DETALLADO Y TRAZABILIDAD ==========================
# Esta fase genera un resumen completo del preprocesamiento, útil para revisión clínica, trazabilidad
# en interfaz gráfica, o informes automatizados por sujeto. El log incluye:
# - Artefactos detectados en el QC
# - Filtros aplicados
# - Maxwell, ICA y componentes eliminados
# - Duración, canales malos, head_pos y resumen clínico

log_data = pd.DataFrame([{
    "subject": subject_id,
    "input_file": os.path.basename(input_path),
    "n_channels": raw.info["nchan"],
    "bad_channels": ",".join(bad_channels),
    "duration_sec": round(raw.n_times / raw.info["sfreq"], 2),
    "maxwell_applied": maxwell_applied,
    "notch_applied": notch_applied,
    "bandpass_applied": bandpass_applied,
    "ica_applied": "yes" if ica_applied == "yes" else "no",
    "ica_method_fallback": "yes" if fallback_used else "no",
    "ica_n_components_removed": len(exclude) if exclude else 0,
    "head_pos_detected": head_pos_detected,
    "qc_drift_present": drift_present,
    "qc_muscle_present": muscle_present,
    "qc_clipping_present": has_clipping,
    "qc_n_bad_channels": n_bad_pre,
    "qc_events_ok": events_ok,
    "qc_summary": " + ".join([
        flag for flag, present in {
            "drift": drift_present,
            "muscle": muscle_present,
            "clipping": has_clipping,
            ">3 bad chans": n_bad_pre > 3,
            "no events": not events_ok
        }.items() if present
    ]) if anomaly_detected else "✔️ sin anomalías relevantes"
}])

# Guardamos o actualizamos el CSV de log
log_data.to_csv(log_path, mode='a', header=not os.path.exists(log_path), index=False)
print(f"📝 Registro clínico-técnico guardado en: {log_path}")


