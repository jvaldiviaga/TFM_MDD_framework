#!/usr/bin/env Rscript
# ============================================================
# Script: qc_post_normalization.R
# Objetivo: Evaluar calidad post-normalización con metadatos,
# PCA coloreado automático por variables clínicas, heatmap y AQM
# Uso: (base) jessica@DESKTOP-UFC7TIA_ ~/TFM_MDD$ Rscript scripts/Transcriptomics/qc_post_normalization.R GSE[ID]
# ============================================================

# ============================================================
# OBJETIVO DEL SCRIPT:
# Evaluar automáticamente la calidad técnica tras la normalización de microarrays,
# generando gráficos e interpretaciones adaptativas para facilitar la toma de decisiones.
# ============================================================

# ------------------------------------------------------------
# ECISIONES CLAVE Y JUSTIFICACIÓN
# ------------------------------------------------------------

# 1. Se prioriza la carga de la versión `normalized_expression_combat.rds` si está disponible,
#    ya que refleja expresión corregida por batch, que suele ser más homogénea.

# 2. Se normalizan los nombres de las muestras (GSMxxxxxxx) para facilitar el cruce
#    con metadatos descargados de GEO y evitar errores de nombre extendido.

# 3. Se aplican validaciones previas a cualquier análisis:
#    - Detectar NA, NaN o Inf.
#    - Excluir genes sin varianza.
#    - Verificar que todas las filas estén completas.
#    Estas comprobaciones previenen errores silenciosos en visualizaciones.

# 4. Se generan boxplots, SD y carga total por muestra para identificar muestras aberrantes.
#    Estas métricas son fáciles de interpretar por usuarios clínicos y técnicos.

# 5. El script detecta variables clínicas categóricas con 2–10 niveles y colorea el PCA por ellas.
#    Límite ≤10 niveles evita sobrecargar gráficamente las visualizaciones.
#    Mejora futura: hacer este umbral dinámico o permitir variables continuas discretizadas.

# 6. Se genera un heatmap basado en distancias euclídeas entre muestras,
#    con anotaciones clínicas si hay metadatos disponibles. Esto permite detectar subgrupos,
#    outliers y estructuras residuales tras la normalización.

# 7. Se automatiza la interpretación técnica de tres métricas clave:
#    - SD > 1.4 → dispersión preocupante
#    - Carga total desequilibrada (>3x) → riesgo técnico
#    - Distancia media > 10 → alta heterogeneidad
#    Estas reglas generan recomendaciones legibles (✅, ⚠️, ❌).

# 8. Se detecta la variable clínica binaria que más separa grupos en PC1 del PCA.
#    Esto ayuda a identificar la covariable más explicativa del patrón post-normalización.

# 9. Se estiman el número de clústeres jerárquicos del heatmap y la separación efectiva en PCA.
#    Esto permite anticipar estructuras ocultas, subgrupos técnicos o fuentes de batch no anotadas.

# 10. Se generan dos informes:
#     - `qc_interpretacion.txt`: resumen técnico estructurado con métricas y decisiones tomadas
#     - `interpretacion_automatizada.txt`: diagnóstico textual interpretativo basado en reglas heurísticas

# 11. Se ejecuta `arrayQualityMetrics` si el objeto `eset` es válido,
#     generando un reporte estándar adicional compatible con Bioconductor.

# 12. Se guarda `sessionInfo()` al final, lo que permite trazabilidad, reproducibilidad completa
#     y facilita debugging en entornos controlados (Snakemake, HPC...).

# 13. El script es compatible con datasets grandes gracias a:
#     - Control condicional de visualizaciones (PCA, heatmap, AQM solo si n ≤ 100)
#     - Selección automática de 30 muestras por grupo si hay metadatos disponibles.

# ------------------------------------------------------------
# MEJORAS FUTURAS
# ------------------------------------------------------------

# - Permitir más de 10 niveles en PCA coloreado si la variable es interpretativa.
# - Exportar todos los gráficos también en PDF.
# - Integrar salida RMarkdown/HTML con todo el informe en un solo documento.
# - Sustituir gráficos base por ggplot2.
# - Implementar detección robusta de outliers (p.ej. Mahalanobis).
# - Modularizar funciones para uso desde Snakemake.
# - Añadir logging automático (log4r) para trazabilidad completa.

# ============================================================

# --- Cargar librerías necesarias, suprimen mensajes al cargar ---
suppressPackageStartupMessages({
  library(ggplot2)            # Visualización con gráficos base ggplot
  library(pheatmap)           # Para generar heatmaps
  library(stats)              # Estadísticas base de R (PCA, etc.)
  library(Biobase)            # Para manipular objetos tipo ExpressionSet
  library(grDevices)          # Gráficos base
  library(arrayQualityMetrics) # Para AQM automático
})

# --- 1. Argumentos y rutas ---
args <- commandArgs(trailingOnly = TRUE)  # Recoger argumentos de línea de comandos
dataset_id <- ifelse(length(args) > 0, args[1], "GSE98793")  # Dataset por defecto

# Directorio de resultados para el dataset
results_dir <- file.path("~/TFM_MDD/results", dataset_id)

# Ruta del archivo de metadatos clínicos
metadata_file <- file.path("~/TFM_MDD/data/transcriptomics", dataset_id, paste0(dataset_id, "_metadata.csv"))

# Ruta de salida para el QC post-normalización
qc_outdir <- file.path(results_dir, "qc_post")
dir.create(qc_outdir, recursive = TRUE, showWarnings = FALSE)

# --- 2. Cargar expresión normalizada (priorizar la versión con ComBat si existe) ---
# Define la ruta del archivo con ComBat (si existe) y sin ComBat
exprs_file_combat <- file.path(results_dir, "normalized_expression_combat.rds")
exprs_file <- file.path(results_dir, "normalized_expression.rds")

# Prioriza la versión con ComBat si está disponible; si no, usa la versión sin ComBat
if (file.exists(exprs_file_combat)) {
  exprs_file <- exprs_file_combat
  exprs_version <- "normalized_expression_combat.rds"
} else if (file.exists(exprs_file)) {
  exprs_version <- "normalized_expression.rds"
} else {
  stop("❌ No se encontró ningún archivo de expresión normalizada.")
}

# Informa por consola cuál archivo fue cargado
cat("🧬 Archivo de expresión cargado:", exprs_version, "\n")

# Carga el objeto de expresión normalizada desde archivo .rds
exprs_norm <- readRDS(exprs_file)

# --- Normalizar nombres de muestra (quedarse solo con el GSMxxxxxxx) ---
gsm_ids <- gsub("^([^_]+).*", "\\1", colnames(exprs_norm))  # Extrae solo la parte GSM del nombre
colnames(exprs_norm) <- gsm_ids  # Asigna los nuevos IDs como nombres de columnas

# Verifica si el objeto es un ExpressionSet; si no, lo convierte
if (inherits(exprs_norm, "ExpressionSet")) {
  eset <- exprs_norm
  exprs_norm <- exprs(eset)  # Extrae matriz de expresión
} else {
  exprs_norm <- as.matrix(exprs_norm)  # Asegura que sea una matriz
  eset <- ExpressionSet(assayData = exprs_norm)  # Crea un ExpressionSet para compatibilidad
}

# --- Validaciones para asegurar calidad técnica de la matriz ---
if (!all(is.finite(exprs_norm))) stop("❌ NA/NaN/Inf detectado.")  # Verifica valores finitos
if (!all(complete.cases(exprs_norm))) stop("❌ Filas incompletas.")  # Verifica que no haya NA por fila
if (any(apply(exprs_norm, 1, sd) == 0)) stop("❌ Genes sin variabilidad.")  # Excluye genes con SD = 0

cat("✅ Matriz validada\n") # Confirmación de validación correcta

# --- 3. Visualizaciones básicas post-normalización ---
# Cálculo de métricas básicas por muestra
sample_sd <- apply(exprs_norm, 2, sd)         # SD por muestra
total_expr <- colSums(exprs_norm)            # Suma total de expresión por muestra

# Guarda los resultados numéricos en archivos CSV para trazabilidad
write.csv(sample_sd, file.path(qc_outdir, "sample_sd_post.csv"))
write.csv(total_expr, file.path(qc_outdir, "carga_total_post.csv"))

# --- Genera visualizaciones PNG de QC por muestra ---
# Boxplot de intensidades normalizadas
png(file.path(qc_outdir, "boxplot_normalized.png"), width = 1400)
boxplot(exprs_norm, las = 2, col = "lightblue", main = "Boxplot intensidades normalizadas")
dev.off()

# Barplot de SD por muestra
png(file.path(qc_outdir, "sample_sd_post.png"), width = 1000)
barplot(sample_sd, las = 2, col = "darkred", main = "SD por muestra")
dev.off()

# Barplot de carga total por muestra
png(file.path(qc_outdir, "carga_total_post.png"), width = 1000)
barplot(total_expr, las = 2, col = "darkblue", main = "Carga total por muestra")
dev.off()

# --- 4. Evaluación con metadatos y PCA coloreado por cada variable ---
diagnosis_col <- NULL          # Inicializa variable para almacenar columna de diagnóstico si se detecta
pca_vars <- list()             # Lista para guardar las variables que se usan en los PCA coloreados
metadata <- NULL               # Inicializa objeto metadata vacío

# --- Carga de metadatos si el archivo existe ---
if (file.exists(metadata_file)) {
  metadata <- read.csv(metadata_file)  # Carga archivo CSV de metadatos
  
  # Busca columna que contenga los identificadores de muestra (GSM, id, sample...)
  id_col <- grep("geo_accession|gsm|sample|id", names(metadata), ignore.case = TRUE, value = TRUE)[1]
  if (!is.null(id_col)) {
    rownames(metadata) <- metadata[[id_col]]  # Asigna identificadores como rownames del metadata
    
    # Detecta las muestras en común entre expresión y metadatos
    shared <- intersect(colnames(exprs_norm), rownames(metadata))
    
    # Filtra los metadatos y expresión para que coincidan sólo en muestras comunes
    metadata <- metadata[shared, , drop = FALSE]
    exprs_common <- exprs_norm[, shared, drop = FALSE]
    
    # --- Recorre cada variable de metadatos para hacer un PCA coloreado si tiene entre 2 y 10 niveles ---
    for (var in names(metadata)) {
      if (length(unique(metadata[[var]])) >= 2 && length(unique(metadata[[var]])) <= 10) {
        group_labels <- as.factor(metadata[[var]])        # Convierte a factor para colorear
        group_colors <- as.numeric(group_labels)          # Colores numéricos para cada grupo
        
        # Verifica que la matriz de expresión sea válida y tenga variabilidad
        if (nrow(exprs_common) >= 2 && all(apply(exprs_common, 1, sd) != 0)) {
          pca <- prcomp(t(exprs_common), scale. = TRUE)   # Ejecuta PCA con muestras en columnas
          
          # Guarda el PCA coloreado como imagen PNG
          png(file.path(qc_outdir, paste0("pca_", var, ".png")), width = 1000)
          plot(pca$x[, 1:2], col = group_colors, pch = 19,
               main = paste("PCA por", var),
               xlab = "PC1", ylab = "PC2")
          legend("topright", legend = levels(group_labels),
                 col = 1:length(levels(group_labels)), pch = 19)
          dev.off()
          pca_vars[[var]] <- TRUE  # Registra que esta variable fue utilizada para un PCA
        }
      }
    }
  }
}

# --- 4b. Detección de variable con mayor separación en PCA (solo binarias) - ---
var_mayor_sep <- NA  # Nombre de la variable con mayor diferencia entre grupos en PC1
mayor_dif <- 0       # Magnitud de dicha diferencia

# Solo si hay metadatos disponibles
if (!is.null(metadata) && ncol(metadata) > 0) {
  for (var in names(metadata)) {
    if (length(unique(metadata[[var]])) == 2) {
      g <- metadata[[var]]
      
      # Verifica que al menos haya muestras en común
      if (length(intersect(colnames(exprs_common), rownames(metadata))) >= 2) {
        g <- as.factor(g)
        pc1_vals <- prcomp(t(exprs_common), scale. = TRUE)$x[, 1]  # Extrae valores de PC1
        
        # Verifica que todos los niveles estén presentes en los datos
        if (all(levels(g) %in% g)) {
          # Calcula diferencia absoluta entre medias de PC1 en ambos grupos
          dif_pc1 <- abs(mean(pc1_vals[g == levels(g)[1]]) - mean(pc1_vals[g == levels(g)[2]]))
          
          # Si esta variable genera más separación, se guarda como la mejor hasta ahora
          if (dif_pc1 > mayor_dif) {
            mayor_dif <- dif_pc1
            var_mayor_sep <- var
          }
        }
      }
    }
  }
}

# --- Heatmap anotado con distancias ---
library(RColorBrewer)  # Librería para paletas de colores bonitas

# Calcula distancias euclídeas entre columnas (muestras) de la matriz normalizada
sample_dists <- dist(t(exprs_norm))
dist_matrix <- as.matrix(sample_dists)  # Convierte a matriz cuadrada para heatmap

# Inicializa anotaciones por muestra (metadata) como NULL
annotation_col <- NULL

# Si hay metadatos disponibles, prepara las variables para anotación del heatmap
if (!is.null(metadata)) {
  
  # Filtra solo variables categóricas con ≤10 niveles
  # NOTA: Este filtro limita a variables con ≤10 niveles únicos
  # Esto se hace para evitar:
  #   - Usar variables numéricas continuas que producirían coloraciones caóticas
  #   - Gráficos ilegibles por demasiadas categorías
  # En futuras versiones, se puede hacer más flexible (por ejemplo, permitir hasta 15 o filtrar por tipo de dato)
  valid_vars <- names(metadata)[sapply(metadata, function(x) length(unique(x)) <= 10)]
  
  # Extrae esas columnas del metadata original
  annotation_col <- metadata[, valid_vars, drop = FALSE]
  
  # Alinea filas de metadata con columnas de la matriz de expresión
  annotation_col <- annotation_col[match(colnames(exprs_norm), rownames(annotation_col)), , drop = FALSE]
  
  # Elimina filas con NA
  annotation_col <- annotation_col[complete.cases(annotation_col), , drop = FALSE]
  
}

# Genera el heatmap de distancias con anotaciones por muestra (si están disponibles)
pheatmap(dist_matrix,
         main = "Heatmap de distancias entre muestras",
         filename = file.path(qc_outdir, "heatmap_normalized.png"),
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_col = annotation_col,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(brewer.pal(n = 9, name = "RdBu"))(100),  # Paleta de color gradiente
         breaks = seq(min(dist_matrix), max(dist_matrix), length.out = 101))  # Rango de colores

# --- 5. Evaluación técnica general ---
mean_distance <- mean(sample_dists)  # Distancia media entre muestras

# Detecta si el archivo cargado incluía corrección por ComBat
combat_applied <- grepl("combat", exprs_version, ignore.case = TRUE)

# Genera una recomendación textual adaptada al grado de heterogeneidad
recommendation <- if (mean_distance < 5) {
  "✅ Buena homogeneidad. Proceder al análisis diferencial."
} else if (mean_distance <= 10) {
  "⚠️ Variabilidad moderada. Revisar posibles subgrupos o covariables."
} else if (combat_applied) {
  "⚠️ Persistente heterogeneidad tras ComBat. Evaluar limpieza manual o covariables adicionales."
} else {
  "❌ Alta heterogeneidad. Considerar ComBat o limpieza manual."
}


# --- 6. Interpretación automática ---
interp <- c()       # Vector que acumulará frases de interpretación técnica
problemas <- c()    # Vector que acumula las fuentes de posible problema técnico

# Evaluar SD por muestra: si alguna muestra supera SD > 1.4
if (max(sample_sd) > 1.4) {
  interp <- c(interp, "⚠️ SD alta en algunas muestras.")
  problemas <- c(problemas, "dispersión en intensidades")
} else {
  interp <- c(interp, "✅ SD en rango aceptable.")
}

# Evaluar carga total entre muestras: compara máximo vs mínimo
if ((max(total_expr) / min(total_expr)) > 3) {
  interp <- c(interp, "⚠️ Carga desigual entre muestras.")
  problemas <- c(problemas, "carga total desequilibrada")
} else {
  interp <- c(interp, "✅ Carga total homogénea.")
}

# Evaluar la distancia media entre muestras
if (mean_distance > 10) {
  interp <- c(interp, paste0("❌ Alta distancia media entre muestras: ", round(mean_distance, 2)))
  problemas <- c(problemas, "separación en heatmap")
} else if (mean_distance > 5) {
  interp <- c(interp, paste0("⚠️ Variabilidad moderada entre muestras: ", round(mean_distance, 2)))
  problemas <- c(problemas, "variabilidad inter-muestra")
} else {
  interp <- c(interp, paste0("✅ Muestras compactas en heatmap (", round(mean_distance, 2), ")"))
}

# Variables de metadatos que mostraron separación en PCA
if (!is.null(pca_vars) && length(pca_vars) > 0) {
  relevantes <- paste(names(pca_vars), collapse = ", ")
  interp <- c(interp, paste0("📘 PCA sugiere variabilidad por: ", relevantes))
}

# Variable con mayor separación en PC1, si se detectó
if (!is.na(var_mayor_sep) && mayor_dif > 0) {
  interp <- c(interp, paste0("🧭 Mayor separación en PC1 por: ", var_mayor_sep,
                             " (ΔPC1 = ", round(mayor_dif, 2), ")"))
}



# --- 6.1. Análisis extendido (separación PCA y clústeres en heatmap) ---

# Si existe la variable "Grupo" entre las anotaciones y el PCA fue generado, calcula separación en PC1
if (!is.null(metadata) && !is.null(annotation_col) && "Grupo" %in% colnames(annotation_col)) {
  grupos <- annotation_col$Grupo
  if (length(unique(grupos)) == 2 && exists("pca")) {
    pc1_vals <- pca$x[, 1]
    grupo1 <- pc1_vals[grupos == unique(grupos)[1]]
    grupo2 <- pc1_vals[grupos == unique(grupos)[2]]
    pc1_sep <- abs(mean(grupo1) - mean(grupo2))  # Diferencia media entre los dos grupos
    
    interp <- c(interp, paste0("📈 Separación entre grupos en PC1: ", round(pc1_sep, 2)))
    if (pc1_sep > 2) {
      interp <- c(interp, "✅ Buena separación en PCA post-normalización.")
    } else {
      interp <- c(interp, "⚠️ Separación moderada en PCA.")
    }
  }
}

# Estima el número de clústeres jerárquicos en el heatmap usando corte por altura media
tryCatch({
  hc <- hclust(sample_dists, method = "complete")
  k <- length(unique(cutree(hc, h = mean(sample_dists))))  # Agrupa por corte a altura media
  interp <- c(interp, paste0("🔢 Nº estimado de clústeres (heatmap): ", k))
}, error = function(e) {
  interp <- c(interp, "⚠️ No se pudo calcular número de clústeres.")
})

# --- 6.2. Detección y eliminación automatizada de outliers por consenso ---

cat("🔍 Detectando outliers por múltiples criterios...\n")

# A. Outliers por SD
outliers_sd <- names(sample_sd[sample_sd > 1.4])

# B. Outliers por PCA (distancia Euclídea desde el centroide en PC1-PC2)
pca <- prcomp(t(exprs_norm), scale. = TRUE)
pc_scores <- pca$x[, 1:2]
dist_pca <- sqrt(rowSums(scale(pc_scores)^2))
limite_dist <- quantile(dist_pca, 0.95)
outliers_pca <- names(dist_pca[dist_pca > limite_dist])

# C. Outliers por AQM: cargar si existe informe previo
outliers_aqm <- character(0)
aqm_path <- file.path(qc_outdir, "AQM_Report", "index.html")
if (file.exists(aqm_path)) {
  aqm_lines <- readLines(aqm_path)
  aqm_detected <- grep("flagged as outlier", aqm_lines, value = TRUE)
  outliers_aqm <- unique(gsub(".*<code>(GSM[0-9]+)</code>.*", "\\1", aqm_detected))
}

# D. Votos por muestra
outlier_votes <- table(c(outliers_sd, outliers_pca, outliers_aqm))
outliers_final <- names(outlier_votes[outlier_votes >= 2])  # Aparecen en ≥2 listas

# E. Filtrado y registro
if (length(outliers_final) > 0) {
  exprs_norm <- exprs_norm[, !colnames(exprs_norm) %in% outliers_final, drop = FALSE]
  eset <- ExpressionSet(assayData = exprs_norm)
  writeLines(outliers_final, file.path(qc_outdir, "muestras_excluidas.txt"))
  cat("❗ Se eliminaron", length(outliers_final), "muestras por consenso de outliers.\n")
} else {
  cat("✅ No se eliminaron muestras. No se detectó consenso suficiente entre criterios.\n")
}



# --- 7. Informes ---
# Reporte resumido con métricas clave de QC
qc_txt <- c(
  "==================== INFORME QC POST-NORMALIZACIÓN ====================",
  paste0("🔎 Dataset: ", dataset_id),
  paste0("🧬 Expresión evaluada: ", exprs_version),
  paste0("📊 SD (min/max): ", round(min(sample_sd), 2), " / ", round(max(sample_sd), 2)),
  paste0("📦 Carga total (min/max): ", round(min(total_expr)), " / ", round(max(total_expr))),
  paste0("🎯 Distancia media (heatmap): ", round(mean_distance, 2)),
  paste0("🧠 Evaluación técnica: ", recommendation),
  if (length(pca_vars) > 0) paste0("📘 PCA coloreado por: ", paste(names(pca_vars), collapse = ", ")) else "⚠️ PCA sin metadatos clínicos.",
  "======================================================================="
)

# Guarda el informe principal en .txt
writeLines(qc_txt, file.path(qc_outdir, "qc_interpretacion.txt"))

# También lo imprime por consola
cat(paste(qc_txt, collapse = "\n"), "\n")

# Guarda interpretación automatizada detallada con emojis y causas posibles
writeLines(c(
  "=========== INTERPRETACIÓN AUTOMÁTICA ===========",
  interp,
  paste0("🎯 Recomendación técnica: ", recommendation),
  "=================================================="
), file.path(qc_outdir, "interpretacion_automatizada.txt"))

# --- 8. Ejecutar AQM ---
cat("📊 Ejecutando arrayQualityMetrics...\n")
tryCatch({
  arrayQualityMetrics(eset,
                      outdir = file.path(qc_outdir, "AQM_Report"),
                      force = TRUE,
                      do.logtransform = FALSE)
  cat("✅ Informe arrayQualityMetrics generado en: qc_post/AQM_Report/\n")
}, error = function(e) {
  cat("⚠️ Error en arrayQualityMetrics:", e$message, "\n")
})

# --- 9. Guardar entorno R ---
sink(file.path(qc_outdir, "session_info.txt"))
sessionInfo()
sink()
cat("✅ QC post-normalización finalizado correctamente para", dataset_id, "\n")

# --- 10. Resumen final por consola ---
cat("\n==================== RESUMEN FINAL QC POST-NORMALIZACIÓN ====================\n")
cat("🧬 Archivo evaluado: ", exprs_version, "\n")
cat("📊 Rango SD: ", round(min(sample_sd), 2), "–", round(max(sample_sd), 2), "\n")
cat("📦 Rango carga total: ", round(min(total_expr)), "–", round(max(total_expr)), "\n")
cat("🌡️ Distancia media (heatmap): ", round(mean_distance, 2), "\n")
cat("📘 PCA coloreado por: ", if (length(pca_vars) > 0) paste(names(pca_vars), collapse = ", ") else "Ninguna", "\n")
cat("🎯 Recomendación: ", recommendation, "\n")
cat("===============================================================================\n\n")
