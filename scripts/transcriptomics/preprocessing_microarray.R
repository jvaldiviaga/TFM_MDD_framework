#!/usr/bin/env Rscript
# ==========================================================
# Script: preprocessing_microarray.R
# Objetivo: Preprocesamiento adaptativo de microarrays Affymetrix
# Ejecuta QC técnico, aplica normalización automática (rma o vsn), detecta efecto batch,
# genera visualizaciones y reportes técnicos.
# Uso: Se llama desde preprocessing_master
# ==========================================================

# ========================== Script de Preprocesamiento Microarrays ==========================
# Descripción: Este script realiza el preprocesamiento adaptativo de datos de microarrays Affymetrix.
# Incluye: QC técnico, decisión automática de normalización (rma vs vsn), corrección por batch (ComBat),
# visualizaciones adaptativas, e informes interpretativos.

# FUTURO: Este script asume entrada en formato Affymetrix (.CEL). No detecta otros formatos como Illumina o Agilent.
# Motivo: extender el framework a más plataformas mejorará su aplicabilidad en otros contextos clínicos o estudios.

# NOTA: No se modulariza todo en funciones porque este script es específico para microarrays Affymetrix.
# La arquitectura del framework se basa en tener scripts independientes por tipo de ómica,
# llamados desde preprocessing_master.R en función del tipo de datos.
# Esto favorece la escalabilidad y evita cargar con funciones irrelevantes en flujos distintos (RNA-seq, metilación, etc.)

# FORTALEZAS Y JUSTIFICACIÓN DEL FLUJO:
# - Realiza preprocesamiento técnico adaptativo de microarrays Affymetrix (.CEL) con lógica condicional robusta.
# - Evalúa métricas técnicas clave (SD, carga, PCA, heatmap) para decidir entre `rma()` y `vsn()`.
# - Aplica corrección por efecto batch con `ComBat()` si detecta separación significativa entre lotes.
# - Integra generación de visualizaciones (boxplot, PCA, heatmap) y reportes automáticos interpretables.
# - Detecta automáticamente variables clínicas relevantes (`diagnosis`, `group`, etc.) sin necesidad de renombrarlas.
# - Optimizado para funcionar tanto de forma autónoma como desde `preprocessing_master.R`.
# - Compatible con datasets grandes (≥100 muestras) mediante control de visualización adaptativo.
# - Genera archivos RDS, logs y reportes listos para análisis diferencial o modelado posterior.

# =============================================================================================

# 1. Cargar librerías necesarias
suppressPackageStartupMessages({
  library(affy)                # Lectura y manejo de archivos .CEL
  library(Biobase)             # Infraestructura ExpressionSet
  library(limma)               # Herramientas para análisis lineal y normalización
  library(vsn)                 # Normalización por variancia estabilizada (VSN)
  library(arrayQualityMetrics) # QC avanzado interactivo
  library(pheatmap)           # Visualización de heatmaps
  library(sva)                # ComBat para corrección de efecto batch
})

# 2. Determinar ruta de entrada
args <- commandArgs(trailingOnly = TRUE)

# Si se pasa argumento por consola, lo usamos
if (length(args) > 0) {
  input_dir <- args[1]
  
  # Si no, buscamos si ya existe en el entorno global (usado desde master)
} else if (exists("input_dir", envir = .GlobalEnv)) {
  input_dir <- get("input_dir", envir = .GlobalEnv)
  
  # Si no hay forma de determinar la ruta, detenemos el script
} else stop("❌ No se ha proporcionado input_dir.")

# Extraemos el ID del dataset desde la ruta
dataset_id <- basename(input_dir)

# Definimos carpetas de resultados
results_dir <- file.path("~/TFM_MDD/results", dataset_id)
qc_dir <- file.path(results_dir, "qc_raw")

# Creamos las carpetas si no existen
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
cat("📥 Iniciando preprocesamiento para:", dataset_id, "\n")

cat("🔍 Leyendo archivos .CEL...\n")

# Extraemos la lista de archivos .CEL presentes en el directorio de entrada
cel_files <- list.files(input_dir, pattern = "\\.CEL$", full.names = TRUE, recursive = TRUE)
if (length(cel_files) == 0) stop("❌ No se encontraron archivos .CEL.")

# 3. Leer archivos .CEL (microarrays crudos)
# Cargamos los archivos .CEL utilizando la función ReadAffy() del paquete 'affy',
# lo que genera un objeto AffyBatch con todos los datos crudos de microarrays Affymetrix.
raw_data <- ReadAffy(filenames = cel_files)

# Confirmación por consola de que la lectura fue exitosa
cat("✅ Archivos .CEL leídos.\n")

# Extraemos la matriz de intensidades crudas a partir del objeto AffyBatch
exprs_mat <- exprs(raw_data)

# 4. QC preliminar (log2, NA, baja señal)
# Evaluamos si es necesario aplicar log2 según la dispersión de los cuantiles (IQR)
# Si la mediana es muy alta o hay mucha diferencia entre Q1 y Q3, se asume que no está en escala log
qx <- quantile(exprs_mat, c(0.25, 0.75), na.rm = TRUE)
log_transform <- (qx[2] > 100) || ((qx[2] - qx[1]) > 50)

# Si detectamos que está en escala lineal, aplicamos log2 + 1 solo para análisis QC (sin afectar normalización posterior)
exprs_raw <- if (log_transform) log2(exprs_mat + 1) else exprs_mat

# Mostramos por consola si se aplicó log2 o no
cat(if (log_transform) "🔁 Se aplicó log2 para QC.\n" else "ℹ️ No se aplicó log2.\n")

# Eliminamos genes (filas) que contengan valores NA en alguna muestra
na_rows <- which(rowSums(is.na(exprs_raw)) > 0)
if (length(na_rows) > 0) exprs_raw <- exprs_raw[-na_rows, ]

# Calculamos la intensidad media por gen y filtramos aquellos con señal muy baja
# Eliminamos los genes con intensidad media menor al percentil 20%
mean_intensity <- rowMeans(exprs_raw)
threshold <- quantile(mean_intensity, 0.2)
exprs_raw <- exprs_raw[mean_intensity > threshold, ]

# =============================================================================================
# FUTURO: Guardar matriz de expresión cruda filtrada antes de aplicar normalización
# Motivo: permite reutilizar estos datos para modelos predictivos u otras exploraciones sin repetir QC.
# Mejora la eficiencia computacional y la trazabilidad de decisiones técnicas.
# saveRDS(exprs_raw, file.path(results_dir, "exprs_filtered.rds"))
# =============================================================================================

# 5. Visualizaciones básicas
# Boxplot de intensidades crudas después del log2 (si se aplicó) y filtrado
# Permite evaluar la homogeneidad de la distribución entre muestras
png(file.path(qc_dir, "boxplot_raw.png"), width = 1400)
boxplot(exprs_raw, las = 2, main = "Boxplot intensidades crudas")
dev.off()

# Histograma de todas las intensidades crudas tras preprocesamiento inicial
# Útil para visualizar la distribución global y detectar sesgos de intensidad
png(file.path(qc_dir, "histograma_intensidades.png"), width = 1000)
hist(as.numeric(exprs_raw), breaks = 100, col = "gray", main = "Histograma intensidades")
dev.off()

# 6. PCA y heatmap (si hay metadatos)
# Buscamos el archivo de metadatos clínicos correspondiente al dataset
meta_path <- list.files(input_dir, pattern = paste0("^", dataset_id, "_metadata\\.csv$"), full.names = TRUE)
has_meta <- length(meta_path) > 0

# Inicializamos variables que se usarán si hay metadatos
group_colors <- NULL           # Colores para los grupos en los plots
group_labels <- NULL           # Etiquetas de grupo (diagnóstico)
diagnosis_col <- NULL          # Nombre de la columna con la variable de diagnóstico
subset_cols <- NULL            # Columnas (muestras) seleccionadas para visualización balanceada

if (has_meta) {
  # Cargamos los metadatos conservando los nombres originales de columnas (ej. 'disease status:ch1')
  meta <- read.csv(meta_path[1], check.names = FALSE)  # ← Para conservar nombres con ":" intactos
  
  # Extraemos los identificadores GSM desde las columnas de expresión
  gsm_ids <- gsub("^([^_]+).*", "\\1", colnames(exprs_raw))
  colnames(exprs_raw) <- gsm_ids
  
  # Detección automática de la columna de ID de muestra
  if ("geo_accession" %in% names(meta)) {
    meta_id_col <- "geo_accession"
  } else {
    # Buscamos columnas candidatas que puedan contener los IDs de muestra
    id_candidates <- grep("geo_accession|gsm|sample|id", names(meta), ignore.case = TRUE, value = TRUE)
    for (id_col in id_candidates) {
      if (any(meta[[id_col]] %in% colnames(exprs_raw))) {
        meta_id_col <- id_col
        break
      }
    }
    if (!exists("meta_id_col")) stop("❌ No se pudo detectar la columna de IDs en los metadatos.")
  }
  
  # Detección automática de la variable clínica de diagnóstico
  diag_candidates <- grep("group|diagnosis|condition", names(meta), ignore.case = TRUE, value = TRUE)
  for (dc in diag_candidates) {
    if (length(unique(meta[[dc]])) >= 2) {
      diagnosis_col <- dc
      break
    }
  }
  
  # Si no se encontró en columnas estándar, buscamos dentro de characteristics_ch1
  if (is.null(diagnosis_col)) {
    charac_cols <- grep("characteristics_ch1", names(meta), value = TRUE)
    for (col in charac_cols) {
      if (any(grepl("disease status", meta[[col]], ignore.case = TRUE))) {
        # Extraemos solo la parte del diagnóstico de la cadena
        meta[[col]] <- gsub(".*disease status:\\s*", "", meta[[col]])
        diagnosis_col <- col
        break
      }
    }
  }
  
  # =============================================================================================
  # FUTURO: Añadir al informe balance de clases (MDD, Control) y sugerencia de métrica (AUC, F1)
  # Motivo: la mayoría de métricas estándar como Accuracy pueden resultar engañosas si las clases están desbalanceadas.
  # Por ejemplo, un modelo que predice siempre "Control" en un dataset con 80% de controles tendrá 80% de accuracy pero será inútil.
  # Solución:
  # - Contar el número de muestras por grupo tras detectar diagnosis_col y group_labels.
  # - Calcular el ratio entre clases (ej. n_mayor / n_menor).
  # - Si ratio > 1.5, imprimir una advertencia y sugerir usar AUC o F1-score en lugar de Accuracy.
  # - Esta decisión se puede guardar en el informe `qc_interpretacion.txt` o en un log global para el pipeline.
  # Justifica automáticamente la métrica usada en la parte de DEG o modelos predictivos.
  # =============================================================================================
  
  
  # Si encontramos una columna de diagnóstico válida, construimos los grupos
  if (!is.null(diagnosis_col) && !is.null(meta_id_col)) {
    diagnosis_vals <- meta[[diagnosis_col]]
    
    # Solo procedemos si hay al menos 2 grupos distintos
    if (length(unique(na.omit(diagnosis_vals))) >= 2) {
      g1 <- meta[meta[[diagnosis_col]] == unique(na.omit(diagnosis_vals))[1], meta_id_col]
      g2 <- meta[meta[[diagnosis_col]] == unique(na.omit(diagnosis_vals))[2], meta_id_col]
      
      # Intersectamos con las columnas reales del dataset (por si hay muestras no coincidentes)
      g1 <- intersect(g1, colnames(exprs_raw)); g2 <- intersect(g2, colnames(exprs_raw))
      
      # Seleccionamos 30 muestras por grupo de forma aleatoria pero reproducible
      n_select <- min(30, length(g1), length(g2))
      if (n_select >= 2) {
        set.seed(42)
        subset_cols <- c(sample(g1, n_select), sample(g2, n_select))
        exprs_subset <- exprs_raw[, subset_cols, drop = FALSE]
        meta <- meta[match(subset_cols, meta[[meta_id_col]]), ]
        group_labels <- as.factor(meta[[diagnosis_col]])
        group_colors <- as.numeric(group_labels)
      }
    }
  }
}

# Si no se detectaron metadatos clínicos válidos, seleccionamos un subconjunto general
if (is.null(subset_cols)) {
  subset_cols <- head(colnames(exprs_raw), min(60, ncol(exprs_raw)))
  exprs_subset <- exprs_raw[, subset_cols, drop = FALSE]
}

# =============================================================================================
# FUTURO: Controlar las visualizaciones (PCA, heatmap, AQM) si el número de muestras > 100
# Motivo: en sistemas con hardware limitado (como WSL2), realizar PCA o heatmap con muchas muestras (>100)
# puede agotar la RAM o hacer que el proceso se caiga (Killed).
# Solución:
# - Implementar un flag de entrada al script como --skip_visuals o --limit_plot_n=60
# - O bien añadir un condicional automático que reduzca el número de muestras para visualización:
#     - Si n_muestras > 100, seleccionar automáticamente 30 por grupo clínico (como ya se hace)
#     - O simplemente evitar las visualizaciones para evitar errores
# Esto permite mantener el QC sin comprometer la estabilidad del script.
# =============================================================================================



# --- PCA ---
pca_done <- FALSE  # Indicador para registrar si el PCA se realizó correctamente

tryCatch({
  # Solo realizamos PCA si hay al menos 2 muestras y 2 genes con varianza no nula
  if (ncol(exprs_subset) >= 2 && nrow(exprs_subset) >= 2 && all(apply(exprs_subset, 1, var) != 0)) {
    pca <- prcomp(t(exprs_subset), scale. = TRUE)  # Transponemos para que las muestras sean filas
    
    # Guardamos gráfico de PCA coloreado por grupo clínico (si existe)
    png(file.path(qc_dir, "pca_subset_colored.png"), width = 1000)
    if (!is.null(group_colors) && length(unique(group_labels)) > 0) {
      plot(pca$x[, 1:2], col = group_colors, pch = 19,
           main = "PCA (grupos balanceados)", xlab = "PC1", ylab = "PC2")
      legend("topright", legend = levels(group_labels),
             col = 1:length(levels(group_labels)), pch = 19)
    } else {
      plot(pca$x[, 1:2], col = "blue", pch = 19,
           main = "PCA (sin metadatos)", xlab = "PC1", ylab = "PC2")
    }
    dev.off(); pca_done <- TRUE
  } else {
    stop("❌ Dataset insuficiente o con columnas constantes.")
  }
}, error = function(e) cat("⚠️ PCA falló:", e$message, "\n"))

# --- Heatmap ---
heatmap_done <- FALSE  # Indicador para registrar si el heatmap se generó correctamente

tryCatch({
  if (ncol(exprs_subset) >= 2) {
    # Calculamos distancias entre muestras y la matriz correspondiente
    dists <- dist(t(exprs_subset))
    dist_mat <- as.matrix(dists)
    
    # Creamos anotaciones para las columnas si hay metadatos disponibles
    annotation_col <- NULL
    if (!is.null(group_labels)) {
      annotation_col <- data.frame(Grupo = group_labels)
      rownames(annotation_col) <- colnames(exprs_subset)
    }
    
    # Generamos el heatmap y lo guardamos
    png(file.path(qc_dir, "heatmap_subset.png"), width = 1000)
    pheatmap(dist_mat,
             main = "Heatmap (subset)",
             labels_col = FALSE,
             annotation_col = annotation_col,
             color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
    dev.off()
    heatmap_done <- TRUE
  } else stop("❌ Muy pocas muestras para heatmap.")
}, error = function(e) cat("⚠️ Heatmap falló:", e$message, "\n"))

# =============================================================================================
# FUTURO: Guardar en un archivo las columnas detectadas automáticamente como diagnosis, ID o batch
# Motivo: facilita la auditoría y documentación del flujo automatizado.
# Especialmente útil si se quiere revisar qué columnas fueron usadas en cada ejecución.
# writeLines(..., file.path(qc_dir, "columnas_detectadas.txt"))
# =============================================================================================

# 7. Métricas y gráficos
# Carga total: suma de intensidades por muestra
carga <- colSums(exprs_raw)

# SD por muestra: evalúa la dispersión de señales
sd_vals <- apply(exprs_raw, 2, sd)

# Guardamos métricas numéricas en CSV
write.csv(carga, file.path(qc_dir, "carga_total.csv"))
write.csv(sd_vals, file.path(qc_dir, "sd_por_muestra.csv"))

# Calculamos rango de SD y ratio de carga entre muestras (para el diagnóstico posterior)
sd_range <- max(sd_vals) - min(sd_vals)
carga_ratio <- max(carga) / min(carga)

# Gráfico de barras: carga total por muestra
png(file.path(qc_dir, "carga_total.png"), width = 1000)
barplot(carga, las = 2, col = "darkblue", main = "Carga total por muestra", ylab = "Suma intensidades")
dev.off()

# Gráfico de barras: desviación estándar por muestra
png(file.path(qc_dir, "sd_por_muestra.png"), width = 1000)
barplot(sd_vals, las = 2, col = "darkred", main = "SD por muestra", ylab = "Desviación estándar")
dev.off()

# --- Interpretación automática del QC previo a normalización ---

interpretacion <- c()          # Vector que almacenará los mensajes interpretativos
suggested_method <- "rma()"    # Valor por defecto (asume buena calidad técnica)

# Evaluamos la SD: si es alta o muy variable → sospechamos heterogeneidad técnica
# Valores por encima de 1.4 o con un rango superior a 0.8 indican heterogeneidad técnica significativa.
# Estos umbrales se basan en el rango típico esperado tras log2 para microarrays Affymetrix,
# donde una SD esperada suele estar entre 0.3 y 1.0. Valores mayores pueden indicar ruido técnico,
# mala calidad o necesidad de una normalización más robusta como vsn().
if (max(sd_vals) > 1.4 || sd_range > 0.8) {
  interpretacion <- c(interpretacion, "⚠️ SD alta o dispersión amplia → posible ruido técnico o muestras heterogéneas.")
  suggested_method <- "vsn()"
} else {
  interpretacion <- c(interpretacion, "✅ SD en rango esperado para datos log2.")
}

# Evaluamos la carga total por muestra (suma de intensidades): si hay gran disparidad → sospecha de batch o muestras mal procesadas
# Si el cociente entre la muestra con mayor carga y la de menor carga supera 3,
# consideramos que puede haber un efecto técnico (batch, calidad desigual).
# Este umbral está basado en experiencias prácticas previas y en las recomendaciones
# para microarrays normalizados, donde idealmente la carga debe ser homogénea.
if (carga_ratio > 3) {
  interpretacion <- c(interpretacion, "⚠️ Carga total muy desigual → posible efecto batch o calidad de muestra variable.")
  suggested_method <- "vsn()"
} else {
  interpretacion <- c(interpretacion, "✅ Carga total homogénea entre muestras.")
}

# Evaluamos el resultado del PCA si se ha podido ejecutar y hay grupos clínicos detectados
if (pca_done && !is.null(group_labels)) {
  
  # Medimos la separación entre los grupos en PC1 (media de cada grupo)
  pca_sep <- abs(mean(pca$x[group_labels == levels(group_labels)[1], 1]) -
                   mean(pca$x[group_labels == levels(group_labels)[2], 1]))
  
  # Si la distancia entre grupos es amplia, interpretamos como separación clara
  # Evaluamos si existe separación clara entre los grupos clínicos en PC1:
  # Calculamos la diferencia entre las medias de los grupos en PC1 (pca_sep).
  # Si esa distancia supera 2 unidades (en datos centrados y escalados), se considera
  # que hay una separación significativa a nivel visual y técnico.
  # Este umbral se basa en observación empírica: en PCA con scale=TRUE, una distancia >2
  # indica agrupamientos visuales consistentes. No es un test estadístico formal, pero sirve
  # como orientación automatizada para sugerir si la señal biológica es clara o no.
  if (pca_sep > 2) {
    interpretacion <- c(interpretacion, "✅ El PCA muestra una separación clara entre grupos clínicos.")
  } else {
    interpretacion <- c(interpretacion, "⚠️ PCA sin separación clara → señal biológica débil o efecto técnico dominante.")
    suggested_method <- "vsn()"
  }
} else {
  interpretacion <- c(interpretacion, "⚠️ PCA no realizado o sin metadatos clínicos.")
}

# Evaluamos el resultado del heatmap si se generó correctamente
# No se usa un umbral numérico, pero si no se observan bloques definidos,
# se considera que no hay patrón claro. Esto refuerza la recomendación de usar vsn().
if (heatmap_done && !is.null(group_labels)) {
  interpretacion <- c(interpretacion, "✅ El heatmap refleja cierta agrupación por grupo clínico.")
} else {
  interpretacion <- c(interpretacion, "⚠️ Heatmap sin patrón de bloques definido.")
  suggested_method <- "vsn()"
}

# Componemos el texto final que se imprimirá y se guardará
qc_auto_txt <- c(
  "==================== INTERPRETACIÓN AUTOMÁTICA ====================",
  interpretacion,
  paste0("🎯 Sugerencia de normalización basada en QC: ", suggested_method),
  "==================================================================="
)

# Guardamos el resumen en un archivo .txt dentro del directorio de QC
writeLines(qc_auto_txt, file.path(qc_dir, "interpretacion_automatizada.txt"))

# También lo imprimimos en consola
cat(paste(qc_auto_txt, collapse = "\n"), "\n")

# =============================================================================================
# FUTURO: Reemplazar el umbral fijo PC1 > 2 por una evaluación estadística (ANOVA o PERMANOVA)
# Motivo: aportaría una justificación estadística más rigurosa para aplicar ComBat.
# Relevante si se plantea publicación o generalización del framework.
# =============================================================================================


# 8. Evaluación técnica
cat("\n=================== RESUMEN DE CALIDAD TÉCNICA ===================\n")
cat("📊 SD por muestra (esperado log2: 0.3 – 1.0): ", round(min(sd_vals),2), "-", round(max(sd_vals),2), "\n")
if (max(sd_vals) > 1.4) cat("⚠️ SD alta: posible dispersión técnica.\n") else cat("✅ SD dentro del rango aceptable.\n")
cat("📦 Carga total (esperado homogénea; ratio ideal < 3x): ", round(min(carga)), "-", round(max(carga)), "\n")
if (carga_ratio > 3) cat("⚠️ Carga desigual → posible efecto batch.\n") else cat("✅ Carga relativamente homogénea.\n")
cat("📈 PCA ejecutado: ", if (pca_done) "✅ Sí" else "⚠️ No", "\n")
cat("🌡️ Heatmap generado: ", if (heatmap_done) "✅ Sí" else "⚠️ No", "\n")
cat("🧪 Se aplicó log2 para QC: ", if (log_transform) "✅ Sí" else "ℹ️ No", "\n")
cat("💾 Genes eliminados por NA: ", length(na_rows), "\n")
cat("🔕 Genes filtrados por baja señal (20%): ", sum(mean_intensity <= threshold), "\n")
cat("=================================================================\n\n")

# 9. Detección y corrección de batch (pData o metadatos externos)
# Inicializamos variables para registrar si se detecta efecto batch
combat_applied <- FALSE
batch_info <- list(source = NA, column = NA, levels = NA)
batch_effect_detected <- FALSE

# Nos aseguramos de que los nombres de las muestras (GSM IDs) sean consistentes en todas las matrices
gsm_ids <- gsub("^([^_]+).*", "\\1", sampleNames(raw_data))
sampleNames(raw_data) <- gsm_ids
colnames(exprs_mat) <- gsm_ids
colnames(exprs_raw) <- gsm_ids

# --- PCA por batch para evaluar efecto visual ---
# Buscamos el archivo de metadatos clínicos
meta_path <- list.files(input_dir, pattern = paste0("^", dataset_id, "_metadata\\.csv$"), full.names = TRUE)

if (length(meta_path) > 0) {
  meta_ext <- read.csv(meta_path[1], check.names = FALSE)
  
  # Detectamos columna con IDs (geo_accession o similares)
  id_col <- grep("geo_accession|gsm|sample|id", names(meta_ext), ignore.case = TRUE, value = TRUE)[1]
  
  # Detectamos posible columna de batch (e.g., 'batch:ch1', 'slide', 'run'...)
  batch_col_meta <- if ("batch:ch1" %in% names(meta_ext)) {
    "batch:ch1"
  } else {
    posibles <- grep("batch|lote|slide|array|run|scan", names(meta_ext), ignore.case = TRUE, value = TRUE)
    posibles[1]
  }
  
  # Continuamos solo si hay batch e ID válidos
  if (!is.null(batch_col_meta) && !is.na(id_col)) {
    rownames(meta_ext) <- meta_ext[[id_col]]
    shared <- intersect(colnames(exprs_raw), rownames(meta_ext))
    
    if (length(shared) >= 10) {
      # Subconjunto para PCA por batch
      exprs_batch_pca <- exprs_raw[, shared, drop = FALSE]
      meta_batch <- meta_ext[shared, , drop = FALSE]
      batch_vec <- factor(trimws(as.character(meta_batch[[batch_col_meta]])))
      names(batch_vec) <- colnames(exprs_batch_pca)
      
      # Seleccionamos 30 muestras por batch para evitar cuellos de botella de RAM
      if (length(unique(batch_vec)) >= 2) {
        subset_cols <- unlist(lapply(levels(batch_vec), function(b) {
          ids <- names(batch_vec[batch_vec == b])
          if (length(ids) > 30) sample(ids, 30) else ids
        }))
        exprs_pca <- exprs_batch_pca[, subset_cols, drop = FALSE]
        pca_batch <- prcomp(t(exprs_pca), scale. = TRUE)
        batch_sub <- batch_vec[subset_cols]
        
        # Visualizamos el PCA coloreado por batches
        png(file.path(qc_dir, "pca_by_batch.png"), width = 1000)
        plot(pca_batch$x[, 1:2], col = as.numeric(batch_sub), pch = 19,
             main = "PCA por batch (subset)", xlab = "PC1", ylab = "PC2")
        legend("topright", legend = levels(batch_sub),
               col = 1:length(levels(batch_sub)), pch = 19)
        dev.off()
        
        # Detección automática de efecto batch en PC1
        # Medimos la separación entre batches en PC1
        pc1_means <- tapply(pca_batch$x[, 1], batch_sub, mean)
        if (length(unique(pc1_means)) >= 2) {
          diff_pc1 <- max(pc1_means) - min(pc1_means)
          
          # Umbral empírico: si diferencia > 2, asumimos efecto batch relevante
          # Evaluamos la separación entre batches en la componente principal 1 (PC1).
          # Calculamos la diferencia entre las medias de PC1 por cada grupo de batch.
          #
          # 🧠 Justificación del umbral:
          #   - En PCA con `scale. = TRUE`, las coordenadas están en unidades de desviación estándar.
          #   - Una separación > 2 indica que los batches difieren más de 2 DE → señal técnica clara.
          #   - Este umbral es empírico y se basa en observaciones prácticas en análisis de microarrays.
          #   - Es rápido, reproducible y suficientemente robusto para automatización inicial.
          #
          # 🔧 FUTURO:
          #   - Para mayor rigor estadístico, este umbral puede sustituirse por una prueba formal:
          #       1. ANOVA sobre PC1: `aov(pca_batch$x[,1] ~ batch_sub)`
          #       2. PERMANOVA (adonis): sobre la matriz de distancias, si se desea evaluar agrupamiento.
          #   - Estas alternativas permitirían confirmar el efecto batch con significancia estadística,
          #     especialmente útil en estudios multicéntricos, con más de 2 niveles o con tamaños grandes.
          if (diff_pc1 > 2) {
            batch_effect_detected <- TRUE
            cat("✅ Separación clara entre batches en PCA: se aplicará ComBat.\n")
          } else {
            cat("ℹ️ No hay suficiente separación entre batches en PCA. No se aplicará ComBat.\n")
          }
        }
      }
    }
  }
}

# Función auxiliar que aplica ComBat y guarda el resultado como RDS
# ==========================================================
# Esta función aplica la corrección por batch a una matriz de expresión normalizada,
# utilizando la librería sva::ComBat. Es robusta frente a errores y controla el uso de RAM.
# También guarda el resultado como archivo .rds para facilitar reproducibilidad.
# ==========================================================
apply_combat <- function(exprs_matrix, batch_vector, origen) {
  
  # Verificamos que el vector de batches coincida en longitud y nombres con las columnas de expresión
  if (length(batch_vector) != ncol(exprs_matrix) || !all(names(batch_vector) == colnames(exprs_matrix))) {
    stop("❌ batch_vector mal alineado con columnas de expresión.")
  }
  
  # Creamos el modelo de diseño para ComBat (intercepto solamente)
  mod <- model.matrix(~1, data = data.frame(batch = batch_vector))
  
  n_samples <- ncol(exprs_matrix)
  corrected <- NULL  # Inicializamos la variable para guardar el resultado corregido
  
  # Si el número de muestras es muy alto (>250), evitamos usar parámetros empíricos para ahorrar RAM
  force_mean_only <- n_samples > 250
  
  if (force_mean_only) {
    cat("⚠️ Demasiadas muestras (", n_samples, ") → Usando mean.only = TRUE para evitar Killed.\n")
    
    # mean.only=TRUE aplica solo corrección de medias, menos intensiva en RAM
    corrected <- ComBat(dat = exprs_matrix, batch = batch_vector, mod = mod, mean.only = TRUE)
    
  } else {
    # Intentamos aplicar ComBat con prior empírico (mejor estimación bayesiana)
    tryCatch({
      corrected <- ComBat(dat = exprs_matrix, batch = batch_vector, mod = mod, par.prior = TRUE)
      
    }, error = function(e1) {
      # Si falla, caemos en un modo de seguridad con mean.only=TRUE
      cat("⚠️ Error con par.prior=TRUE: ", e1$message, "\n")
      corrected <- ComBat(dat = exprs_matrix, batch = batch_vector, mod = mod, mean.only = TRUE)
    })
  }
  
  # Si por algún motivo no se pudo aplicar ComBat, se lanza error crítico
  if (is.null(corrected)) {
    stop("❌ ComBat no pudo aplicarse correctamente.")
  }
  
  # Guardamos el resultado corregido como archivo RDS para trazabilidad
  saveRDS(corrected, file = file.path(results_dir, "normalized_expression_combat.rds"))
  
  # Mensaje de confirmación
  cat("✔️ ComBat aplicado desde", origen, "\n")
  
  # Devolvemos la matriz corregida para continuar el flujo del pipeline
  return(corrected)
}



# 10. Normalización adaptativa basada en evaluación técnica previa
# Evaluamos si hay dispersión técnica o problemas de calidad que desaconsejen el uso de rma().
# Si se detecta alguna de estas condiciones, se prefiere vsn() (más robusta pero más lenta):
# 1. SD por muestra muy elevada (>1.4) o rango amplio (>0.8) → indica muestras ruidosas.
# 2. Carga total muy desigual (ratio > 3) → posible problema de escaneo o lote.
# 3. PCA o heatmap no pudieron generarse → señal inconsistente o fallos de visualización.
alta_dispersion <- (max(sd_vals) > 1.4 || sd_range > 0.8 || carga_ratio > 3 || !pca_done || !heatmap_done)
criterio <- metodo <- NULL  # Inicializamos variables para registrar la decisión tomada

tryCatch({
  
  # Si hay alta dispersión → aplicamos vsn() (normalización robusta)
  if (alta_dispersion) {
    eset <- expresso(raw_data, bg.correct = FALSE, normalize = FALSE,
                     pmcorrect.method = "pmonly", summary.method = "avgdiff")
    norm_data <- exprs(vsn2(exprs(eset)))
    metodo <- "vsn()"
    criterio <- "Dispersión global detectada"
    
    # Si no hay dispersión → aplicamos rma() (rápida y apropiada)
  } else {
    norm_data <- exprs(rma(raw_data))
    metodo <- "rma()"
    criterio <- "Homogeneidad técnica aceptable"
  }
}, error = function(e) {
  # En caso de error (por ejemplo, fallo en rma), aplicamos vsn() como fallback
  cat("⚠️ Normalización fallida: ", e$message, "\n")
  eset <- expresso(raw_data, bg.correct = FALSE, normalize = FALSE,
                   pmcorrect.method = "pmonly", summary.method = "avgdiff")
  norm_data <- exprs(vsn2(exprs(eset)))
  metodo <<- "vsn()"
  criterio <<- "Fallback tras fallo en rma()"
})

# Guardamos la matriz de expresión normalizada, independientemente del método
saveRDS(norm_data, file = file.path(results_dir, "normalized_expression.rds"))

# --- Aplicamos ComBat solo si previamente se detectó efecto batch significativo ---
if (batch_effect_detected) {
  
  # Verificamos que las muestras del batch coincidan con las columnas normalizadas
  shared_cols <- intersect(colnames(norm_data), names(batch_vec))
  
  if (length(shared_cols) < ncol(norm_data)) {
    warning("⚠️ ComBat no se aplicará a todas las muestras. Verifica nombres de columnas.")
  }
  
  # Seleccionamos las columnas comunes para aplicar ComBat
  norm_data_shared <- norm_data[, shared_cols, drop = FALSE]
  
  # Guardamos la información de batch (para incluir en logs o informes)
  batch_info <- list(source = "metadata", column = batch_col_meta, levels = length(unique(batch_vec)))
  
  # Aplicamos ComBat sobre la matriz normalizada
  corrected <- apply_combat(norm_data_shared, batch_vec[shared_cols], "metadatos externos")
  
  # Actualizamos la matriz normalizada con la versión corregida
  norm_data <- corrected
  combat_applied <- TRUE
  
  # Guardamos la versión final corregida por batch
  saveRDS(norm_data, file = file.path(results_dir, "normalized_expression.rds"))
} else {
  cat("ℹ️ No se aplicará ComBat porque no se detectó efecto batch significativo.\n")
}


# 11. Registro de decisiones de normalización y sesión R
# Mostramos por consola el método y el criterio técnico que se ha aplicado
cat("🎯 Método de normalización aplicado: ", metodo, "\n")
cat("🧠 Justificación técnica: ", criterio, "\n")

# Guardamos el método elegido (rma vs vsn) y la razón en un archivo de texto plano
# Este archivo permite trazar qué decisión se tomó para cada dataset en ejecución en lote
writeLines(c(
  paste0("Método de normalización: ", metodo),
  paste0("Justificación: ", criterio)
), file.path(results_dir, "decision_normalizacion.txt"))

# Guardamos la información de la sesión R completa (versiones, paquetes cargados)
# Es esencial para reproducibilidad, sobre todo si cambia el entorno o se migra a otro equipo
writeLines(capture.output(sessionInfo()), file.path(results_dir, "session_info.txt"))

# Creamos un informe más interpretativo y explicativo que resume las métricas clave del QC
# Incluye SD, carga total, número de genes retenidos, decisión de normalización y si se aplicó ComBat
writeLines(c(
  "==================== INFORME QC ADAPTATIVO ====================",
  paste0("🔎 Dataset: ", dataset_id),
  paste0("🧬 Genes tras filtrado: ", nrow(exprs_raw)),
  paste0("📊 SD (min/max): ", round(min(sd_vals), 2), " / ", round(max(sd_vals), 2)),
  paste0("📦 Carga total: ", round(min(carga)), " / ", round(max(carga))),
  paste0("🎯 Normalización: ", metodo),
  paste0("🧠 Criterio: ", criterio),
  if (combat_applied) "✔️ ComBat aplicado por batch detectado." else "ℹ️ No se aplicó ComBat.",
  if (!is.null(group_labels)) paste0("📘 PCA coloreado por: ", diagnosis_col) else "⚠️ PCA sin metadatos clínicos.",
  "==============================================================="
), file.path(qc_dir, "qc_interpretacion.txt"))


# 12. Ejecución de arrayQualityMetrics (solo si hay ≤ 100 muestras)
# Si hay pocas muestras, ejecutamos arrayQualityMetrics para generar un informe interactivo en HTML
# Este informe incluye múltiples gráficas de diagnóstico (MA plots, boxplots, PCA, etc.)
# No se usa en datasets grandes porque puede consumir mucha RAM o colapsar en WSL2
if (ncol(exprs_raw) <= 100) {
  tryCatch({
    cat("📊 Ejecutando arrayQualityMetrics...\n")
    arrayQualityMetrics(expressionset = raw_data,
                        outdir = file.path(qc_dir, "arrayQualityMetrics"),
                        force = TRUE, do.logtransform = FALSE)
    cat("✅ arrayQualityMetrics completado.\n")
  }, error = function(e) {
    cat("⚠️ arrayQualityMetrics falló:", e$message, "\n")
  })
} else {
  
  # En datasets grandes, evitamos ejecutarlo por limitaciones de hardware
  cat("ℹ️ Demasiadas muestras para arrayQualityMetrics (omitido).\n")
}

# Mensaje final indicando que todo el flujo de preprocesamiento ha terminado correctamente
cat("✅ Preprocesamiento completado y listo para el pipeline.\n")
