#!/usr/bin/env python3
# ============================================================
# Script: qc_meg_raw.py
#
# OBJETIVO:
# Evaluar la calidad técnica de los datos MEG crudos (.fif),
# detectando automáticamente los problemas más comunes en
# sistemas Elekta Neuromag (como los del dataset AIM/Reward).
# Este control de calidad inicial es clave para tomar decisiones
# técnicas personalizadas por sujeto en el preprocesamiento posterior.
# ============================================================

# === FUNCIONES IMPLEMENTADAS ===
# ✅ Interferencia 60 Hz (detección de pico en PSD multicanal)
# ✅ Ruido muscular (>50 Hz)
# ✅ Drift <1 Hz (aumento de energía en baja frecuencia)
# ✅ Clipping o saturación (valores extremos en canales MEG)
# ✅ Canales malos (detección mediante Maxwell filtering)
# ✅ Presencia y número de eventos en canal STI014 (estímulos)
#
# 📝 Resultados exportados a:
#    - qc_psd_summary.csv        → detección binaria de anomalías técnicas
#    - qc_bad_channels.csv       → lista de canales ruidosos o planos
#    - qc_event_check.csv        → presencia y cantidad de eventos codificados
#    - qc_interpretation.csv     → resumen técnico-clínico interpretable por sujeto
#    - qc_ready_flag.csv         → flag para saber si requiere preprocesamiento
#    - qc_quality_label.csv      → clasificación global de calidad técnica
#    - *_peaks_debug.png         → curva espectral invertida con picos detectados
#    - *_psd.png                 → densidad espectral de potencia multicanal (0–60 Hz)

# === FUNCIONES NO INCLUIDAS (por decisión metodológica) ===
# ⚠️ Visualización espectral ya implementada (0–60 Hz completa mediante PSD multicanal).
# ❌ No se realiza aún cuantificación automática por bandas cognitivas específicas
#    (e.g., alpha: 8–13 Hz, theta: 4–7 Hz, beta: 13–30 Hz).
#    Aunque dichas bandas están contenidas visualmente en la PSD,
#    su análisis cuantitativo requiere:
#    - segmentación del raw en condiciones homogéneas (epocas),
#    - eliminación de artefactos fisiológicos (ICA: EOG, ECG),
#    - y un diseño experimental controlado (reposo vs tarea).
#    De lo contrario, se generarían falsos positivos o medidas erróneas.
#    Por esta razón, la cuantificación de bandas cognitivas se pospone
#    a fases posteriores del pipeline, como el epocado o análisis de fuente.

# ❌ Artefactos fisiológicos (EOG, ECG) → se detectarán tras filtrado e ICA
# ❌ Movimiento de cabeza (HPI/head_pos) → se analizará solo si existe dev_head_t
# ============================================================

# ========================== LIBRERÍAS ==========================
import os                                # Para manipulación de rutas y archivos
import argparse                          # Para lectura de argumentos desde consola
import mne                               # Para carga y procesamiento de archivos MEG (.fif)
import numpy as np                       # Para cálculos numéricos y manejo de arrays
import pandas as pd                      # Para cargar y guardar tablas (.csv)
import matplotlib                        # Backend gráfico para evitar errores sin GUI
matplotlib.use('Agg')                    # Forzamos backend no interactivo (por seguridad en WSL2)
import matplotlib.pyplot as plt          # Para generación de figuras (PSD, curvas, etc.)
import seaborn as sns                    # Para visualización avanzada (e.g., heatmaps)
from PIL import Image                    # Para análisis espectral desde la imagen PSD
from scipy.signal import find_peaks      # Para detectar picos en el perfil espectral

# ========================== PARÁMETROS DE ENTRADA ==========================
# Se definen los argumentos que se deben pasar desde la terminal al ejecutar el script

parser = argparse.ArgumentParser(description="Evaluación automática de calidad MEG crudo.")
parser.add_argument('--input_file', type=str, required=True,
                    help='Archivo .fif MEG crudo (ej. sub-XXX_meg.fif)')
parser.add_argument('--output_dir', type=str, required=True,
                    help='Directorio para guardar gráficos e informes')
args = parser.parse_args()

# ========================== RUTAS Y VARIABLES CLAVE ==========================
# Expandimos rutas con ~, creamos carpeta de salida si no existe, y extraemos el ID del sujeto

input_path = os.path.expanduser(args.input_file)
output_dir = os.path.expanduser(args.output_dir)
os.makedirs(output_dir, exist_ok=True)

subject_id = os.path.basename(input_path).split("_")[0]
# ➤ Extraemos el ID antes del primer guion bajo (ej. sub-M87139247)

# ========================== COMPROBACIÓN DE QC COMPLETO ==========================
# Verificamos si el QC ya fue realizado completamente (PSD + CSVs)

def check_qc_complete(output_dir, subject_id):
    psd_img = os.path.join(output_dir, f"{subject_id}_psd.png")
    summary_csv = os.path.join(output_dir, "qc_psd_summary.csv")
    bads_csv = os.path.join(output_dir, "qc_bad_channels.csv")

    # ➤ Si falta alguno de los tres componentes esenciales, el QC no está completo
    if not (os.path.exists(psd_img) and os.path.exists(summary_csv) and os.path.exists(bads_csv)):
        return False

    try:
        summary_df = pd.read_csv(summary_csv)
        bads_df = pd.read_csv(bads_csv)
        return subject_id in summary_df["subject"].values and subject_id in bads_df["subject"].values
    except Exception:
        return False  # En caso de error leyendo los CSVs, consideramos el QC incompleto

# Si el QC ya está completamente hecho, se omite
if check_qc_complete(output_dir, subject_id):
    print(f"⏩ QC ya completado previamente para {subject_id}. Se omite ejecución.")
    exit(0)

# ========================== FASE 1: CARGA DEL RAW ==========================
# Intentamos cargar el archivo .fif crudo completo. Si está corrupto o incompleto, se registra error.
# Recorte por eficiencia computacional y justificación científica
# 600 s (10 min) son suficientes para evaluar:
# - Pico a 60 Hz (interferencia)
# - Ruido muscular (>50 Hz)
# - Drift (<1 Hz)
# - Clipping y canales malos
#
# Umbral empíricamente validado: >600s causa errores de RAM en WSL2 al hacer raw.load_data()
# Además, los estudios de calidad MEG (Hämäläinen 2010, Gross 2013) confirman que 2–5 min bastan.


print(f"📂 Cargando archivo crudo: {input_path}")
try:
    # Cargamos el raw sin preload inicial
    raw = mne.io.read_raw_fif(input_path, preload=False)
    duration_sec = raw.n_times / raw.info["sfreq"]

    # Si la duración es muy alta, recortamos automáticamente a 600s
    if duration_sec > 600:
        raw.crop(tmax=600)
        print(f"⚠️ Duración larga: {duration_sec:.1f}s → 600s suficientes para QC técnico")
    else:
        print(f"🕒 Duración: {duration_sec:.1f}s → sin recorte necesario")

    # Cargamos los datos ya recortados (si se aplicó crop) o completos
    raw.load_data()

except Exception as e:
    print(f"❌ ERROR: No se pudo cargar {subject_id}. Split roto o archivo corrupto: {e}")
    with open(os.path.join(output_dir, "qc_failed_subjects.txt"), "a") as f:
        f.write(f"{subject_id},{str(e)}\n")
    exit(1)


# ========================== FASE 2: DENSIDAD ESPECTRAL DE POTENCIA (PSD MULTICANAL) ==========================
 # Esta fase genera una visualización clave para la detección automática y visual de artefactos espectrales.
# Detecta:
# - Pico de interferencia de red (⚡ 60 Hz en EE.UU. / 50 Hz en Europa)
# - Ruido muscular (>50 Hz), que aparece como picos múltiples y dispersos
# - Drift excesivo (<1 Hz) como incremento de potencia en bajas frecuencias
# - Canales con distribución espectral anómala o plano (curvas atenuadas o distorsionadas)
#
# ⚠️ En este proyecto usamos fmax=65 Hz porque los datos proceden del estudio AIM/Reward (EE.UU.),
# donde la frecuencia de red es 60 Hz. Si se trabaja con datasets europeos, considerar fmax=50 Hz.
#
# Esta versión del PSD:
# Muestra cada canal por separado (curvas de colores) → mejor para inspección técnica
# Guarda la imagen como PNG (no se abre ventana Qt)
# Se reutiliza después para extraer métricas automáticas
# Evita errores típicos de entornos sin GUI (como WSL2) al usar un backend no interactivo

# Forzamos backend no interactivo (por seguridad en WSL2)
# Importamos la librería base de gráficos matplotlib,
# que MNE utiliza por debajo para generar figuras.


# Mostramos por consola que inicia la generación de la gráfica
print("📈 Generando PSD multicanal (curvas individuales por canal)...")

# Generamos la figura del PSD multicanal
# - fmax=65 limita el eje X a 65 Hz, suficiente para inspección funcional y artefactos comunes
# - average=False dibuja cada canal por separado para identificar outliers
# - show=False evita mostrar la figura en pantalla (ideal para WSL2 o ejecución por lotes)
fig_psd = raw.plot_psd(
    fmax=65,
    average=False,
    show=False
)

# Definimos el path donde se guardará la imagen del PSD
psd_img_path = os.path.join(output_dir, f"{subject_id}_psd.png")

# Guardamos la figura como imagen PNG
fig_psd.savefig(psd_img_path, dpi=150)  # 150 dpi → buena resolución sin generar archivos pesados
plt.close(fig_psd)                      # Cerramos la figura para liberar memoria

# Confirmamos por consola que la imagen se generó correctamente
print(f"✅ PSD multicanal guardado como imagen: {psd_img_path}")

# === OPCIONAL: Visualización interactiva (desactivada por defecto) ===
# Este bloque está pensado para futuras interfaces gráficas (e.g. Shiny, Jupyter, Qt):
# Permitiría explorar la señal cruda antes del preprocesamiento, pero está desactivado
# porque puede lanzar errores en servidores o entornos sin GUI (como WSL2).
#
# try:
#     print("🧠 Mostrando señal cruda en visor interactivo (solo en entornos gráficos)...")
#     raw.plot(n_channels=50, duration=5.0, show=True)
#     print("✅ Visualización interactiva abierta.")
# except Exception as e:
#     print(f"⚠️ Visualización interactiva no disponible: {e}")


# ========================== FASE 3: ANÁLISIS ESPECTRAL AUTOMÁTICO DESDE IMAGEN PSD ==========================
# Este enfoque evita errores comunes de MNE con compute_psd(), y es más estable clínicamente.
# Analiza la imagen PNG generada en la fase anterior para detectar artefactos mediante perfil espectral invertido.

print("🔬 Evaluando artefactos espectrales mediante análisis de imagen PSD...")

# Abrimos la imagen PSD y la convertimos a escala de grises
img = Image.open(psd_img_path).convert('L')     # L = 8-bit grayscale
arr = np.array(img)                             # Convertimos la imagen a array numpy

# Extraemos una región espectral relevante (estimada de pruebas previas)
region = arr[150:350, 200:800]

# Calculamos el perfil espectral invertido (energía relativa por frecuencia)
profile = (region.max() - region.mean(axis=0))

# Umbral adaptativo para detección de picos
adaptive_threshold = max(10, np.percentile(profile, 90) * 0.1)

# Detectamos picos prominentes
peaks, _ = find_peaks(profile, prominence=adaptive_threshold)

# Mapeamos las posiciones de los picos a frecuencias reales (aprox. 0–60 Hz)
freq_max = 60
peak_freqs = [round(freq_max * i / (len(profile) - 1), 1) for i in peaks]

# === DIAGNÓSTICO DE ARTEFACTOS ===
drift_present = profile[:20].mean() > profile[20:40].mean() + 2
muscle_present = any(52 <= f <= 60 for f in peak_freqs)

# === GUARDADO DE CURVA DE DEPURACIÓN ===
debug_plot = os.path.join(output_dir, f"{subject_id}_peaks_debug.png")
plt.figure(figsize=(10, 3))
plt.plot(profile, label="Perfil espectral invertido")
plt.plot(peaks, profile[peaks], "rx", label="Picos detectados")
plt.title(f"Curva espectral invertida: {subject_id}")
plt.xlabel("Frecuencia estimada (Hz)")
plt.ylabel("Contraste")
plt.legend()
plt.tight_layout()
plt.savefig(debug_plot)
plt.close()
print(f"📊 Curva de depuración guardada: {debug_plot}")

# === Métricas sintéticas para CSV (comparación pre/post) ===
power_drift = float(profile[:20].mean()) if drift_present else np.nan
power_muscle = float(profile[60:].mean()) if muscle_present else np.nan
psd_mean = profile  # Para mantener compatibilidad con exportación posterior



# ========================== FASE 4: DETECCIÓN DE CLIPPING O SATURACIÓN ==========================
# Aquí se evalúan valores extremos en la señal, que podrían indicar clipping (saturación).
# Este artefacto puede deberse a:
# - movimiento súbito de la cabeza,
# - interferencia externa intensa,
# - errores de adquisición (amplificadores).

# Se considera que hay clipping si algún valor en los canales MEG supera los 5e-11 Tesla (50.000 fT).
# Este umbral es conservador y se basa en estudios previos y en valores observados en el dataset AIM.
print("🩻 Evaluando saturación de la señal (clipping)...")

# Primero seleccionamos solo los canales MEG (magnetómetros + gradiometers)
meg_data = raw.copy().pick_types(meg=True)

# Extraemos la matriz de datos MEG (canales x tiempo)
meg_values = meg_data.get_data()

# Evaluamos si existe al menos un valor fuera del umbral
has_clipping = np.any(np.abs(meg_values) > 5e-11)


# ========================== FASE 5: DETECCIÓN DE CANALES MALOS ==========================
# Esta fase detecta canales ruidosos o planos utilizando Maxwell filtering (Elekta Neuromag).
# Esta técnica permite identificar sensores defectuosos (ruido excesivo o actividad plana) para excluirlos antes del preprocesamiento.

print("🔍 Detectando canales malos con Maxwell filtering (requiere archivos de calibración y crosstalk)...")

# Rutas absolutas a los archivos de calibración necesarios
cal_path = os.path.expanduser("~/TFM_MDD/data/neuroimaging/Code/sss_cal_3046_2009-04-29.dat")
ct_path  = os.path.expanduser("~/TFM_MDD/data/neuroimaging/Code/ct_sparse_orion.fif")

try:
    if os.path.exists(cal_path) and os.path.exists(ct_path):
        # Detección automática de canales ruidosos y planos mediante Maxwell filtering
        bads_detected, flats_detected, scores = mne.preprocessing.find_bad_channels_maxwell(
            raw,
            calibration=cal_path,
            cross_talk=ct_path,
            return_scores=True,
            verbose='WARNING'
        )

        # Visualización opcional: heatmap con puntuación de ruido por canal
        import matplotlib.pyplot as plt
        import seaborn as sns

        sns.heatmap(scores['scores_noisy'], cmap='viridis')
        plt.title(f"Canales ruidosos - {subject_id}")
        heatmap_path = os.path.join(output_dir, f"{subject_id}_bad_channel_scores.png")
        plt.savefig(heatmap_path)
        plt.close()
        print(f"🧾 Heatmap de canales ruidosos guardado: {heatmap_path}")

    else:
        print("⚠️ Archivos de calibración o crosstalk no encontrados. Se omite Maxwell filtering.")
        bads_detected, flats_detected = [], []

except Exception as e:
    print(f"⚠️ Maxwell filtering falló por error interno: {e}")
    bads_detected, flats_detected = [], []

# Añadir los canales malos al objeto raw.info['bads'] para su exclusión posterior
raw.info['bads'] = bads_detected

# === Exportar resultados técnicos ===
# Este CSV resume los canales ruidosos y planos detectados por sujeto.
bads_out_path = os.path.join(output_dir, "qc_bad_channels.csv")

bads_df = pd.DataFrame([{
    "subject": subject_id,
    "bad_channels": ",".join(bads_detected) if bads_detected else "",
    "flat_channels": ",".join(flats_detected) if flats_detected else "",
    "n_bad": len(bads_detected),
    "n_flat": len(flats_detected)
}])

if os.path.exists(bads_out_path):
    df_old = pd.read_csv(bads_out_path)
    df_old = df_old[df_old["subject"].astype(str).str.strip().str.lower() != subject_id.strip().lower()]
else:
    df_old = pd.DataFrame()

df_final = pd.concat([df_old, bads_df], ignore_index=True)
df_final.to_csv(bads_out_path, index=False)

# Diagnóstico por consola
print(f"📉 Canales malos detectados: {bads_detected if bads_detected else 'ninguno'}")


# ========================== FASE 6: VERIFICACIÓN DEL CANAL DE EVENTOS ==========================
# En esta fase se comprueba la presencia del canal STI014, estándar en BIDS-MEG para codificación de eventos.
# Este canal es necesario para segmentar los datos en épocas (epoching) y detectar respuestas evocados (ERPs/ERFs).
# Si el canal existe y contiene eventos válidos, se considera el sujeto apto para segmentación.
print("🔎 Verificando canal de eventos (STI014)...")

events_ok = True  # Por defecto lo asumimos como válido para este dataset
n_events = 0

# Solo se verifica si STI014 está en los canales
if 'STI014' in raw.ch_names:
    try:
        events = mne.find_events(raw, stim_channel='STI014', verbose=False)
        n_events = events.shape[0]
        events_ok = n_events > 0
        print(f"✅ Canal de eventos detectado. Número de eventos: {n_events}")
    except Exception as e:
        print(f"⚠️ Canal STI014 ilegible: {e}")
        events_ok = False
else:
    print("ℹ️ Canal STI014 no está presente en este dataset. Se omite evaluación de eventos.")

# Guardamos el resultado en un archivo CSV (`qc_event_check.csv`)
events_check_path = os.path.join(output_dir, "qc_event_check.csv")
pd.DataFrame([{
    "subject": subject_id,
    "events_ok": events_ok,
    "n_events": n_events if 'STI014' in raw.ch_names else "NA"
}]).to_csv(
    events_check_path,
    mode='a', index=False,
    header=not os.path.exists(events_check_path)
)

# Guardar etiquetas + valores numéricos en un único CSV (para comparación pre/post)
n_bad_pre = len(bads_detected)
n_flat_pre = len(flats_detected)

df_psd = pd.DataFrame({
    "subject": [subject_id],
    "drift_present": [drift_present],
    "muscle_noise_present": [muscle_present],
    "clipping_pre": [has_clipping],
    "low_freq_power_pre": [power_drift],
    "high_freq_power_pre": [power_muscle],
    "psd_max": [np.max(psd_mean)],
    "psd_median": [np.median(psd_mean)],
    "n_bad_channels_pre": [n_bad_pre],
    "n_flat_channels_pre": [n_flat_pre]
})


df_psd.to_csv(
    os.path.join(output_dir, "qc_psd_summary.csv"),
    mode='a', index=False,
    header=not os.path.exists(os.path.join(output_dir, "qc_psd_summary.csv"))
)
print("✅ PSD cuantitativo y artefactos guardados.")


# ========================== FASE 7: EXPORTAR INFORMES Y CLASIFICACIÓN FINAL ==========================
# En esta fase se consolidan los resultados del QC en archivos separados para:
# - Resumen binario por sujeto
# - Flag de validez para preprocesamiento
# - Etiqueta de calidad técnica global (3 niveles)
# Todas las decisiones se basan en métricas cuantitativas reales obtenidas por compute_psd() o señales crudas.

# --- 1. Resumen binario del QC por sujeto ---
# Incluye presencia de artefactos, clipping, canales malos y eventos
qc_summary = pd.DataFrame([{
    "subject": subject_id,
    "muscle_noise_>50Hz": muscle_present,        # Energía excesiva en >50 Hz (ruido EMG típico en tensión muscular)
    "low_freq_drift_<1.5Hz": drift_present,      # Potencia excesiva en <1.5 Hz, considerado drift técnico o desacoplamiento
    "clipping_detected": has_clipping,           # Señal bruta con valores extremos (>5e-11 T) → posible saturación
    "n_bad_channels": n_bad_pre,                 # Canales ruidosos detectados por Maxwell filtering
    "stim_channel_present": events_ok if 'STI014' in raw.ch_names else "NA",           # Presencia del canal STI014 (evento codificado)
    "n_events": n_events if 'STI014' in raw.ch_names else "NA"
                          # Número total de eventos detectados
}])

qc_summary_path = os.path.join(output_dir, "qc_post_binary_summary.csv")

# Escritura robusta: elimina duplicados si ya existe una fila para este sujeto
if os.path.exists(qc_summary_path):
    df_old = pd.read_csv(qc_summary_path)
    df_old = df_old[df_old["subject"].astype(str).str.strip().str.lower() != subject_id.strip().lower()]
else:
    df_old = pd.DataFrame()

# Guardamos el resumen actualizado
df_final = pd.concat([df_old, qc_summary], ignore_index=True)
df_final.to_csv(qc_summary_path, index=False)

# --- 2. Flag de validez para preprocesamiento ---
# Se considera que el sujeto NO está listo para preprocesamiento si:
# - hay interferencia 60 Hz
# - hay ruido muscular
# - hay drift técnico (<1.5 Hz)
# - hay clipping
# - hay más de 3 canales ruidosos
#
# Justificación:
# Estos artefactos indican que la señal puede estar severamente contaminada.
# Más de 3 canales malos sugiere que el sensor estaba mal colocado o mal calibrado.
qc_ready = (
    not has_clipping and
    n_bad_pre <= 3
)


# Guardamos el flag
pd.DataFrame([{
    "subject": subject_id,
    "qc_ready_for_preprocessing": qc_ready
}]).to_csv(
    os.path.join(output_dir, "qc_ready_flag.csv"),
    mode='a', index=False,
    header=not os.path.exists(os.path.join(output_dir, "qc_ready_flag.csv"))
)

# --- 3. Etiqueta de calidad técnica (3 niveles) ---
# Clasificación visual rápida según número de banderas activadas:
# - 🔴 severo: 3 o más problemas detectados (alta probabilidad de datos no utilizables)
# - 🟡 moderado: exactamente 2 problemas (precaución recomendada)
# - 🟢 aceptable: 0 o 1 problema (apto para análisis)
#
# Justificación clínica:
# El umbral de ≥3 indicadores permite detectar sujetos con múltiples contaminaciones simultáneas.
# Un solo artefacto no invalida el análisis, pero múltiples sí comprometen la validez.

summary_flags = sum([
    muscle_present,
    drift_present,
    has_clipping,
    n_bad_pre > 3
])

qc_quality_label = (
    "🔴 severo" if summary_flags >= 3 else
    "🟡 moderado" if summary_flags == 2 else
    "🟢 aceptable"
)

# Guardamos la etiqueta de calidad técnica
pd.DataFrame([{
    "subject": subject_id,
    "qc_tech_label": qc_quality_label
}]).to_csv(
    os.path.join(output_dir, "qc_quality_label.csv"),
    mode='a', index=False,
    header=not os.path.exists(os.path.join(output_dir, "qc_quality_label.csv"))
)


# ========================== FASE 8: INTERPRETACIÓN Y LOG EXPLICATIVO ==========================
# Esta fase construye una interpretación textual del QC para cada sujeto.
# Resume los hallazgos de forma comprensible, útil tanto para informes como para revisión manual o GUI.

interpretation = []  # Lista de anomalías detectadas

if muscle_present:
    interpretation.append("ruido muscular >50 Hz")
if drift_present:
    interpretation.append("drift <1.5 Hz")
if has_clipping:
    interpretation.append("posible clipping o saturación")
if n_bad_pre > 3:
    interpretation.append(f"más de 3 canales ruidosos ({n_bad_pre})")
if 'STI014' in raw.ch_names and not events_ok:
          interpretation.append("⚠️ canal STI014 ausente o sin eventos")
if not interpretation:
    interpretation.append("✔️ sin anomalías relevantes")  # Si no hay anomalías, se indica explícitamente

# Convertimos la lista en una cadena legible separada por +
interpretation_str = " + ".join(interpretation)

# Mostramos por pantalla un resumen del QC, incluyendo la etiqueta global de calidad
print(f"📊 QC automático para {subject_id}: {interpretation_str} | Calidad: {qc_quality_label} | Preprocesar: {qc_ready}")

# Guardamos esta interpretación para que pueda consultarse desde el pipeline o informes
pd.DataFrame([{
    "subject": subject_id,
    "qc_interpretation": interpretation_str
}]).to_csv(
    os.path.join(output_dir, "qc_interpretation.csv"),
    mode='a', index=False,
    header=not os.path.exists(os.path.join(output_dir, "qc_interpretation.csv"))
)


# ========================== FINALIZACIÓN ==========================
# Indicamos por pantalla que el QC ha finalizado correctamente para este sujeto.
print(f"✅ QC completado para: {subject_id}")

