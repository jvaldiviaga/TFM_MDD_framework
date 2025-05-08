#!/usr/bin/env Rscript
# ==========================================================
# Script: differential_expression.R
# Objetivo: Análisis de expresión diferencial entre grupos clínicos
#           con detección automática, explicaciones, visualizaciones
# ==========================================================

# ============================================================
# OBJETIVO GENERAL:
# Este script constituye el núcleo del análisis de expresión diferencial en el framework.
# Ha sido diseñado para ser reutilizable, automatizable y adaptable a distintos datasets,
# integrando desde la limpieza y alineación de datos hasta la exportación de resultados y visualizaciones.

# FORTALEZAS Y JUSTIFICACIÓN DEL FLUJO

# Este script ha sido diseñado como una herramienta modular, automatizable y robusta,
# capaz de adaptarse a múltiples datasets transcriptómicos sin necesidad de intervención manual.

# - Permite ejecución flexible tanto desde consola (CLI) como desde Snakemake, aceptando argumentos dinámicos.
# - Prioriza automáticamente la versión de expresión corregida por batch (ComBat), mejorando la homogeneidad del modelo.
# - Carga y alinea correctamente los metadatos y la matriz de expresión, incluso si los IDs de muestra vienen en formatos no estándar.
# - Permite analizar cualquier variable clínica de interés, con o sin covariables técnicas o clínicas de ajuste (batch, sexo, edad...).
# - Utiliza `limma`, un paquete de referencia para análisis de expresión diferencial, validado y reproducible.
# - Filtra automáticamente genes sin varianza antes de modelar, previniendo errores silenciosos en `eBayes`.
# - Detecta de forma automática el coeficiente correcto del diseño estadístico, incluso cuando hay múltiples efectos.
# - Exporta resultados en formatos estandarizados (.csv, .rds, .txt), incluyendo subconjuntos diferenciados (sobreexpresados, subexpresados).
# - Genera visualizaciones interpretables (volcano plot, heatmap, PCA) esenciales para validación clínica y exploración biológica.
# - Incluye control de errores y validaciones explícitas en todas las fases: argumentos, metadatos, niveles, diseño, resultados.
# - Toda la ejecución queda registrada con trazabilidad completa: parámetros, fecha, archivos generados y log de comandos.

# Su estructura modular y su adaptabilidad lo hacen idóneo para su integración dentro de un framework multi-ómico escalable.

# ------------------------------------------------------------
#  LIMITACIONES ACTUALES Y MEJORAS PENDIENTES
# ------------------------------------------------------------

# 1. Limitación en la detección automática de variable clínica: solo 2–3 niveles únicos.
#    -Esta restricción afecta únicamente a la detección automática en `differential_expression.R`
#    cuando no se indica `group_col`. Se excluyen variables válidas como región cerebral,
#    etnicidad o tratamientos con >3 niveles.
#    -Sin embargo, esta limitación se solventa completamente mediante el script
#    `run_differential_expression_interactive.R`, que permite seleccionar cualquier variable,
#    sin importar el número de niveles.
#    -Mejora futura: ampliar también la lógica automática interna para mayor autonomía del sistema.

# 2. Falta verificación del balance entre grupos → riesgo de sesgo en casos extremos.

# 3. Las covariables no se validan por NA o baja varianza antes de incluirse en el modelo.

# 4. La anotación del heatmap puede fallar si hay desalineación entre columnas y metadatos.

# 5. PCA y heatmap se hacen con los top N genes sin control de redundancia ni filtrado adicional.

# 6. Covariables numéricas como edad no se estandarizan ni normalizan antes de modelar.

# 7. La fórmula del diseño se construye como string sin validación estructural.

# 8. No se generan resultados vacíos explícitos si no hay genes significativos.

# 9. Solo se soportan comparaciones binarias → sin contrasts para modelos multigrupo.

# 10. El log de ejecución no incluye hash de los archivos ni versión de los datos utilizados.

# ------------------------------------------------------------
# MEJORAS FUTURAS
# ------------------------------------------------------------

# - Detectar automáticamente variables clínicas con más de 3 niveles y ofrecer advertencias.
# - Validar n por grupo antes de ejecutar limma; lanzar alerta si hay gran desequilibrio.
# - Verificar y filtrar covariables por NA, varianza o colinealidad.
# - Modularizar en funciones: carga, validación, modelo, exportación.
# - Incluir `makeContrasts()` para modelos complejos multigrupo.
# - Exportar DEG también como Excel (.xlsx).
# - Añadir informes RMarkdown con tablas, plots y resumen HTML.
# - Guardar `sessionInfo()` y hashes de datos para trazabilidad completa.

# ============================================================

# --- Carga de librerías necesarias, suprimiendo mensajes ---
suppressPackageStartupMessages({
  library(readr)        # Para leer CSVs de metadatos
  library(limma)        # Modelo lineal para análisis de expresión diferencial
  library(stringr)      # Funciones útiles para procesar nombres de columnas
  library(ggplot2)      # Visualización: volcano plot, PCA
  library(pheatmap)     # Generación de heatmaps
  library(FactoMineR)   # PCA avanzada (complementaria)
  library(factoextra)   # Visualización de PCA (complementaria)
  library(dplyr)        # Manipulación de data frames
})

# 1. Captura de argumentos
args <- commandArgs(trailingOnly = TRUE)  # Captura los argumentos pasados al script
if (length(args) < 1) {
  stop("❌ Uso: Rscript differential_expression.R <dataset_id> [group_column] [reference_level] [covariables_coma]")
}

# Asigna argumentos según el orden
dataset_id <- args[1]  # Obligatorio
group_col  <- if (length(args) >= 2) args[2] else NULL  # Variable clínica a comparar
ref_group  <- if (length(args) >= 3) args[3] else NULL  # Grupo de referencia (control)
covariables <- if (length(args) == 4) unlist(strsplit(args[4], ",")) else character(0)  # Covariables separadas por coma


# 2. Rutas y archivos
# Define la carpeta base con los datos transcriptómicos y metadata
base_dir    <- file.path("~/TFM_MDD/data/transcriptomics", dataset_id)

# Ruta al archivo de metadata ya limpiado (si existe)
meta_clean  <- file.path(base_dir, paste0(dataset_id, "_metadata_clean.csv"))


# Rutas a la expresión normalizada, con prioridad a la versión corregida por ComBat
expr_file_combat <- file.path("~/TFM_MDD/results", dataset_id, "normalized_expression_combat.rds")
expr_file_raw    <- file.path("~/TFM_MDD/results", dataset_id, "normalized_expression.rds")
expr_file <- if (file.exists(expr_file_combat)) expr_file_combat else expr_file_raw  # Escoge el que esté disponible

# Añade sufijo a la carpeta de salida en función de si hay covariables o no
suffix <- if (length(covariables) == 0) {
  "no_covariates"  # Si no hay covariables, lo indica explícitamente
} else {
  paste0("with_", paste(covariables, collapse = "_"))  # Si hay covariables, las incluye en el nombre de carpeta
}

# Define el directorio de salida donde se guardarán resultados del análisis diferencial
out_dir <- file.path("~/TFM_MDD/results", dataset_id, "differential_expression", suffix)

# Crea el directorio si no existe
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)



# 3. Carga de metadata y expresión normalizada
# Si no existe el archivo de metadatos limpio...
if (!file.exists(meta_clean)) {
  cat("⚠️ Metadata limpio no encontrado. Usando el crudo y guardando versión limpia.\n")
  
  # Construye la ruta al archivo original
  meta_raw <- file.path(base_dir, paste0(dataset_id, "_metadata.csv"))
  
  # Si tampoco está el archivo original, se detiene
  if (!file.exists(meta_raw)) stop("❌ No se encontró metadata.csv para limpiar: ", meta_raw)
  
  # Carga el archivo original
  metadata <- read_csv(meta_raw, show_col_types = FALSE)
  
  # Guarda una copia como "limpio" (aunque aquí aún no se transforma)
  write_csv(metadata, meta_clean)
  cat("✅ Metadata limpio guardado como:", meta_clean, "\n")
} else {
    # Si existe el metadata limpio, lo carga directamente
  metadata <- read_csv(meta_clean, show_col_types = FALSE)
}

# Limpia nombres de columnas para evitar problemas (espacios, símbolos)
colnames(metadata) <- make.names(colnames(metadata))

# Carga la matriz de expresión normalizada (con o sin ComBat)
exprs <- readRDS(expr_file)

# Asegura que sea una matriz (por si es data.frame u otro objeto)
if (!is.matrix(exprs)) exprs <- as.matrix(exprs)

# --- 💡 Alineación de IDs GSM entre metadatos y expresión ---
# Extrae solo el ID base de las muestras (por ejemplo: "GSM12345_rep1" → "GSM12345")
gsm_ids <- gsub("^([^_]+).*", "\\1", colnames(exprs))

# Asigna estos IDs como nombres de columna
colnames(exprs) <- gsm_ids

# Hace lo mismo para los rownames del metadata
rownames(metadata) <- gsub("^([^_]+).*", "\\1", rownames(metadata))

# 4. Detección automática de variable clínica
# Si no se indicó ninguna variable clínica por argumento...
if (is.null(group_col)) {
  cat("\n🔍 Buscando columna clínica candidata para comparación...\n")
  cat("ℹ️ Columnas con 2-3 niveles únicos:\n\n")
  
  # Recorre las columnas del metadata buscando candidatas con 2-3 niveles únicos
  for (col in colnames(metadata)) {
    vals <- unique(na.omit(metadata[[col]]))
    if (length(vals) >= 2 && length(vals) <= 3) {
      cat("•", col, "→", paste(vals, collapse = " / "), "\n")
    }
  }
  
  # Se detiene y sugiere usar un script auxiliar para revisar columnas
  stop("\n❌ No se indicó la columna clínica. Por favor, revísala con: Rscript check_metadata_columns.R ", dataset_id)
}

# Limpia el nombre de la variable de grupo (por si tiene símbolos)
group_col_clean <- make.names(group_col)

# Verifica y muestra los niveles existentes en la variable seleccionada
cat("📊 Verificando niveles de la variable:", group_col_clean, "\n")
print(table(metadata[[group_col_clean]]))

# Limpia también nombres de covariables para su uso en fórmulas
cov_clean <- make.names(covariables)

# Si la columna clínica seleccionada no existe en el metadata, se detiene
if (!(group_col_clean %in% colnames(metadata))) {
  stop("❌ Columna clínica no encontrada:", group_col_clean)
}

# 5. Sugerencia automática de nivel de referencia
# Si no se especificó el grupo de referencia, se toma el más frecuente
if (is.null(ref_group)) {
  ref_group <- names(sort(table(metadata[[group_col_clean]]), decreasing = TRUE))[1]
  cat("ℹ️ Nivel de referencia no especificado. Usando por defecto el más frecuente: '", ref_group, "'\n")
}

# Verifica que el grupo de referencia exista en los niveles de la variable
if (!(ref_group %in% unique(metadata[[group_col_clean]]))) {
  stop("❌ Nivel de referencia no válido. Valores posibles: ", paste(unique(metadata[[group_col_clean]]), collapse = ", "))
}


# 6. Diseño experimental
# Asegura que metadata y expresión estén alineadas correctamente por IDs
# Extrae los IDs GSM base desde los nombres de columna de la matriz de expresión
gsm_ids <- gsub("^([^_]+).*", "\\1", colnames(exprs))
# Reasigna los nombres de columna con solo los IDs
colnames(exprs) <- gsm_ids

# Reasigna los rownames del metadata también a IDs base
rownames(metadata) <- gsm_ids

# Intersección de muestras compartidas entre expresión y metadatos
shared_samples <- intersect(colnames(exprs), rownames(metadata))

# Filtra la matriz de expresión y los metadatos para que tengan las mismas muestras
exprs <- exprs[, shared_samples]
metadata <- metadata[shared_samples, , drop = FALSE]

# Convierte la variable de grupo en factor y define el nivel de referencia
group <- factor(metadata[[group_col_clean]])
group <- relevel(group, ref = ref_group)

# Construye la fórmula de diseño, incluyendo grupo y covariables (si existen)
design_formula <- paste("~", paste(c(group_col_clean, cov_clean), collapse = " + "))

# Genera la matriz de diseño del modelo lineal
design <- model.matrix(as.formula(design_formula), data = metadata)

# 7. Análisis con limma
# Filtra genes sin varianza (evita errores en el modelo lineal)
exprs <- exprs[apply(exprs, 1, sd) > 0, ]

# Ajusta el modelo lineal por gen
fit <- lmFit(exprs, design)

# Aplica eBayes para moderar las estadísticas de prueba
fit <- eBayes(fit)

# Muestra los coeficientes disponibles en el modelo ajustado
cat("🧪 Coeficientes disponibles en el modelo:\n")
print(colnames(fit$coefficients))

# Busca el coeficiente que corresponde a la variable de grupo
coef_candidates <- grep(paste0("^", group_col_clean), colnames(fit$coefficients), value = TRUE)

# Si no se encuentra ningún coeficiente relacionado con la variable clínica, aborta
if (length(coef_candidates) == 0) {
  stop("❌ No se encontró ningún coeficiente para la variable: ", group_col_clean,
       "\n📌 Coeficientes disponibles: ", paste(colnames(fit$coefficients), collapse = ", "))
}

# Si hay varios coeficientes, usa el primero
coef_name <- coef_candidates[1]
cat("✅ Usando coeficiente:", coef_name, "\n")

# Extrae todos los genes con sus estadísticas ajustadas (Benjamini-Hochberg)
deg <- topTable(fit, coef = coef_name, number = Inf, adjust = "BH")
deg <- as.data.frame(deg)

# Verifica que la tabla contiene la columna esperada
cat("🔎 Columnas en topTable:\n")
print(colnames(deg))

# Asegura que se generó la columna de p-valor ajustado
if (!"adj.P.Val" %in% colnames(deg)) {
  stop("❌ La tabla de resultados no contiene 'adj.P.Val'. Algo falló en el modelo.")
}


# 8. Guardar resultados
# Guarda todos los resultados
write.csv(deg, file = file.path(out_dir, "DEG_results.csv"))

# Guarda el top 100 por p-valor
write.csv(head(deg, 100), file = file.path(out_dir, "DEG_top100.csv"))


# 8.1 Guardar genes sobreexpresados y subexpresados
# Genes con logFC positivo y p-valor ajustado < 0.05
deg_up <- deg %>% filter(adj.P.Val < 0.05 & logFC > 0)

# Genes con logFC negativo y p-valor ajustado < 0.05
deg_down <- deg %>% filter(adj.P.Val < 0.05 & logFC < 0)

# Guarda ambos subconjuntos
write_csv(deg_up, file.path(out_dir, "DEG_upregulated.csv"))
write_csv(deg_down, file.path(out_dir, "DEG_downregulated.csv"))

# Informa al usuario del número de genes significativos
cat("🟢 Genes sobreexpresados (logFC > 0):", nrow(deg_up), "\n")
cat("🔴 Genes subexpresados (logFC < 0):", nrow(deg_down), "\n")


# 9. Volcano plot
# Añade nombres de gen y flag de significancia
deg$gene <- rownames(deg)
deg$signif <- deg$adj.P.Val < 0.05

# Genera volcano plot y guarda como imagen PNG
png(file.path(out_dir, "volcano_plot.png"), width = 1000)
ggplot(deg, aes(logFC, -log10(P.Value), color = signif)) +
  geom_point(alpha = 0.6) + theme_minimal() +
  scale_color_manual(values = c("grey", "red")) +
  labs(title = "Volcano plot", x = "Log Fold Change", y = "-log10(P-value)")
dev.off()


# 10. Heatmap

# Selecciona los 50 genes más significativos
top_genes <- rownames(head(deg, 50))

# Extrae la expresión de esos genes
heat_expr <- exprs[top_genes, , drop = FALSE]

# Asegura que las columnas están alineadas con los metadatos
heat_expr <- heat_expr[, intersect(colnames(heat_expr), rownames(metadata)), drop = FALSE]

# Convierte a matriz válida
heat_expr <- as.matrix(heat_expr)

# Elimina genes o muestras con NA o valores infinitos
heat_expr <- heat_expr[complete.cases(heat_expr), ]
heat_expr <- heat_expr[, colSums(is.na(heat_expr) | is.infinite(heat_expr)) == 0, drop = FALSE]

# Si quedan muy pocos genes o muestras, se omite el heatmap
if (nrow(heat_expr) < 2 || ncol(heat_expr) < 2) {
  cat("⚠️ No hay suficientes datos para generar el heatmap\n")
} else {
  # Extrae la variable de grupo para anotación de columnas
  annotation_col <- as.data.frame(metadata[, group_col_clean, drop = FALSE])
  rownames(annotation_col) <- rownames(metadata)
  
  # Alinea las columnas del heatmap con la anotación
  annotation_col <- annotation_col[colnames(heat_expr), , drop = FALSE]
  
  # Genera el heatmap
  png(file.path(out_dir, "heatmap_top50.png"), width = 1000, height = 800)
  pheatmap(heat_expr,
           show_rownames = TRUE,
           show_colnames = FALSE,
           annotation_col = annotation_col,
           main = "Heatmap Top 50 DEGs")
  dev.off()
}


# 11. PCA top 100 genes
# Calcula PCA de los 100 genes más significativos
pca_data <- prcomp(t(exprs[rownames(head(deg, 100)), ]), scale. = TRUE)

# Dibuja gráfico de dispersión PCA coloreado por grupo
png(file.path(out_dir, "pca_top100.png"), width = 1000)
plot(pca_data$x[, 1:2], col = as.factor(metadata[[group_col_clean]]), pch = 19,
     main = "PCA - Top 100 DEGs", xlab = "PC1", ylab = "PC2")

# Añade leyenda por grupo
legend("topright", legend = levels(group), col = 1:length(levels(group)), pch = 19)
dev.off()

# 12. Resumen textual
# Total de genes analizados
n_total <- nrow(deg)

# Genes con p-valor ajustado < 0.05 (significativos)
n_sig <- sum(deg$adj.P.Val < 0.05)

# Genes sobreexpresados (logFC > 0 y significativos)
n_up <- sum(deg$adj.P.Val < 0.05 & deg$logFC > 0)

# Genes subexpresados (logFC < 0 y significativos)
n_down <- sum(deg$adj.P.Val < 0.05 & deg$logFC < 0)

# Construye resumen textual amigable con emojis
resumen <- paste0(
  "🧠 Análisis de expresión diferencial para ", dataset_id, "\n",
  "📊 Comparación: ", group_col, " (referencia: ", ref_group, ")\n",
  "🧬 Genes evaluados: ", n_total, "\n",
  "✅ Genes significativos (adj.P.Val < 0.05): ", n_sig, "\n",
  "🔺 Sobreexpresados: ", n_up, " | 🔻 Subexpresados: ", n_down, "\n"
)

# Guarda el resumen en un archivo de texto
writeLines(resumen, file.path(out_dir, "summary_DEG.txt"))

# También lo imprime en consola
cat(resumen)


# 13. Guardar DEG como .rds y log del comando
# Guarda el data frame DEG completo como archivo binario RDS para análisis downstream
saveRDS(deg, file = file.path(out_dir, "DEG_results.rds"))

# Construye un registro detallado del comando ejecutado (reproducibilidad)
command_log <- paste0(
  "🧬 Análisis de expresión diferencial\n",
  "📁 Dataset: ", dataset_id, "\n",
  "📊 Variable de grupo: ", group_col, "\n",
  "⚙️  Referencia: ", ref_group, "\n",
  "🧩 Covariables: ", ifelse(length(covariables) == 0, "Ninguna", paste(covariables, collapse = ", ")), "\n",
  "📅 Fecha de ejecución: ", Sys.time(), "\n"
)

# Guarda ese log en un archivo de texto
writeLines(command_log, file.path(out_dir, "log_comando.txt"))

# Informa por consola que el log ha sido guardado
cat("\n📝 Log del análisis guardado en log_comando.txt\n")
