#!/usr/bin/env Rscript
# ==========================================================
# Script: run_differential_expression_interactive.R
# Objetivo: Guiar al usuario en la exploración del metadata
#           y lanzar el análisis de expresión diferencial
# ==========================================================

# ============================================================
# OBJETIVO:
# Este script interactivo permite al usuario explorar las variables clínicas
# disponibles en el metadata y lanzar el análisis de expresión diferencial
# mediante `differential_expression.R` con parámetros personalizados,
# sin necesidad de modificar código directamente.

# JUSTIFICACIONES DE DISEÑO:
# - Ofrece una descripción interpretada de cada columna (e.g. "diagnosis", "sex", "score"...)
#   usando patrones textuales conocidos para orientar al usuario.
# - Permite seleccionar:
#     • Variable clínica a comparar (e.g., diagnóstico)
#     • Grupo de referencia (por defecto el más frecuente)
#     • Covariables de ajuste (sexo, edad, batch, etc.)
# - Valida visualmente los valores únicos por variable antes de ejecutar el análisis.
# - Simplifica el uso para usuarios no familiarizados con R y evita errores por selección incorrecta.

# MOTIVACIONES:
# - Mejora la reproducibilidad al registrar los parámetros utilizados en cada análisis.
# - Sirve como frontend accesible (CLI) del script `differential_expression.R`.
# - Permite adaptar fácilmente el pipeline a distintos datasets sin necesidad de reprogramar.

# POSIBLES MEJORAS FUTURAS:
# - Detección automática de la mejor variable clínica (binaria o balanceada).
# - Sugerencia automática del grupo de referencia (el más frecuente).
# - Filtrado inteligente de covariables (solo categóricas o relevantes).
# - Soporte para archivos de configuración por dataset (e.g. `GSE12345_config.json`).
# - Registro automático de inputs y selección como log (`log_inputs.txt`).
# - Validación del diseño experimental antes de lanzar el análisis (n, desequilibrio, NA...).
# - Exploración opcional del diseño (frecuencias, balance, visualización previa).
# - Versión futura con interfaz gráfica (Shiny) para ejecución desde navegador.
# ============================================================


# Oculta los mensajes de carga de paquetes para mantener la salida limpia
suppressPackageStartupMessages({
  library(readr)     # Para leer archivos CSV
  library(stringr)   # Para manejo de strings (no usado explícitamente pero útil si se extiende)
})

# Función para interpretar el significado probable de una columna del metadata
interpretar <- function(nombre) {
  n <- tolower(nombre)  # Convierte el nombre a minúsculas para comparación más flexible
  if (grepl("group|diagnos|status", n)) return("🧠 Diagnóstico clínico (MDD vs control)")
  if (grepl("anxiety", n)) return("😰 Ansiedad comórbida (sí/no)")
  if (grepl("severity|level", n)) return("📉 Gravedad clínica (leve/moderada/severa)")
  if (grepl("sex|gender", n)) return("🚻 Sexo biológico")
  if (grepl("age", n)) return("🎂 Edad del sujeto")
  if (grepl("batch|platform|array", n)) return("🧪 Lote técnico / batch")
  if (grepl("score|scale|hamd|bdi", n)) return("📊 Escala clínica")
  return("⚠️ Descripción no disponible")  # Por defecto si no matchea ningún patrón conocido
}

# === Paso 1: Solicita el ID del dataset (ej. GSE98793) ===
cat("🧬 Bienvenida al asistente de análisis de expresión diferencial\n")
dataset_id <- readline("📂 Introduce el ID del dataset (ej. GSE98793): ")

# === Paso 2: Construye rutas a archivos relevantes ===
base_dir <- file.path("~/TFM_MDD/data/transcriptomics", dataset_id)
meta_clean <- file.path(base_dir, paste0(dataset_id, "_metadata_clean.csv"))
meta_raw   <- file.path(base_dir, paste0(dataset_id, "_metadata.csv"))
expr_file_combat <- file.path("~/TFM_MDD/results", dataset_id, "normalized_expression_combat.rds")
expr_file_raw    <- file.path("~/TFM_MDD/results", dataset_id, "normalized_expression.rds")


# Usa la versión con ComBat si existe, si no, la normal
expr_file <- if (file.exists(expr_file_combat)) expr_file_combat else expr_file_raw
cat("🧬 Matriz de expresión utilizada:", basename(expr_file), "\n")

# === Paso 3: Verifica que el archivo de expresión exista ===
if (!file.exists(expr_file)) stop("❌ No se encontró la matriz de expresión normalizada: ", expr_file)

# === Paso 4: Carga metadata limpio si existe, o el crudo si no ===
if (file.exists(meta_clean)) {
  metadata <- read_csv(meta_clean, show_col_types = FALSE)
} else if (file.exists(meta_raw)) {
  cat("⚠️ Metadata limpio no encontrado. Mostrando columnas desde el archivo original.\n")
  metadata <- read_csv(meta_raw, show_col_types = FALSE)
} else {
  stop("❌ No se encontró ningún archivo de metadata disponible.")
}


# Normaliza nombres de columnas para evitar problemas con espacios o caracteres inválidos
colnames(metadata) <- make.names(colnames(metadata))
columnas <- colnames(metadata)

# === Paso 5: Muestra e interpreta cada columna disponible ===
cat("\n📋 Variables disponibles:\n")
for (i in seq_along(columnas)) {
  cat(sprintf("%2d) %s → %s\n", i, columnas[i], interpretar(columnas[i])))
}


# === Paso 6: Muestra ejemplos de valores únicos para cada variable ===
cat("\n🔍 Ejemplos de valores por variable:\n")
for (col in columnas) {
  vals <- unique(na.omit(metadata[[col]]))
  cat("•", col, ": ", paste(head(vals, 5), collapse = ", "))  # Muestra los primeros 5 valores únicos
  if (length(vals) > 5) cat(", ...")
  cat("\n")
}

# === Paso 7: El usuario selecciona la variable clínica para comparar ===
var_idx <- as.integer(readline("\n❓ ¿Qué variable quieres comparar? (número): "))
group_col <- columnas[var_idx]

# === Paso 7.1: Muestra los niveles posibles de esa variable ===
niveles <- unique(na.omit(metadata[[group_col]]))
cat("\n📊 Valores posibles para", group_col, ":\n")
for (i in seq_along(niveles)) {
  cat(sprintf("%2d) %s\n", i, niveles[i]))
}

# === Paso 8: El usuario selecciona el grupo de referencia ===
ref_idx <- as.integer(readline("⚙️  ¿Qué valor será el grupo de referencia? (número): "))
ref_group <- niveles[ref_idx]

# === Paso 9: El usuario elige covariables para ajustar el modelo ===
cat("\n🎯 ¿Qué covariables quieres ajustar? (elige varias separadas por coma, o 0 para ninguna):\n")
for (i in seq_along(columnas)) {
  cat(sprintf("%2d) %s → %s\n", i, columnas[i], interpretar(columnas[i])))
}
cov_idx <- readline("🧩 Covariables (ej. 5,7) o 0: ")

# Procesamiento del input de covariables
covariables <- ""
if (cov_idx != "0") {
  cov_indices <- as.integer(strsplit(cov_idx, ",")[[1]])
  covariables <- paste(columnas[cov_indices], collapse = ",")
}

# === Paso 10: Ejecuta el script de análisis con los parámetros seleccionados ===
cmd <- sprintf("Rscript differential_expression.R %s %s \"%s\" \"%s\"",
               dataset_id, group_col, ref_group, covariables)

# Muestra el comando que se ejecutará
cat("\n🚀 Ejecutando análisis con:\n", cmd, "\n\n")
system(cmd)
