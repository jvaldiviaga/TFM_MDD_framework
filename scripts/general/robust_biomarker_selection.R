#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(glmnet)
  library(randomForest)
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Uso: Rscript robust_biomarker_selection.R <dataset_id>")
dataset_id <- args[1]
base_dir <- file.path(Sys.getenv("HOME"), "TFM_MDD")
data_dir <- file.path(base_dir, "data", "transcriptomics", dataset_id)
results_dir <- file.path(base_dir, "results", dataset_id, "biomarker_selection")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# === Buscar archivos
expr_path <- list.files(path = c(data_dir, file.path(base_dir, "results", dataset_id)),
                        pattern = ".*norm.*\\.rds$", full.names = TRUE, recursive = TRUE)[1]
if (is.na(expr_path)) stop("❌ No se encontró expresión normalizada.")
meta_path <- list.files(path = data_dir,
                        pattern = paste0(dataset_id, "_metadata_(clean|filtered)\\.csv$"),
                        full.names = TRUE)[1]
if (is.na(meta_path)) stop("❌ No se encontró metadata limpio ni filtrado.")
deg_path <- list.files(path = file.path(base_dir, "results", dataset_id, "differential_expression"),
                       pattern = "DEG_results_mapped.*\\.csv$", full.names = TRUE, recursive = TRUE)[1]
if (is.na(deg_path)) stop("❌ No se encontró archivo DEG mapeado.")

# === Cargar datos
expr <- readRDS(expr_path)
meta <- read.csv(meta_path, check.names = FALSE)
deg  <- read.csv(deg_path)

# === Detectar columna de SYMBOL
symbol_col <- grep("^symbol$", names(deg), ignore.case = TRUE, value = TRUE)[1]
if (is.na(symbol_col)) stop("❌ No hay columna SYMBOL mapeada en DEG.")

# === Detección de ID de muestra
id_col <- if ("SampleID" %in% names(meta)) "SampleID" else if ("geo_accession" %in% names(meta)) "geo_accession" else {
  gsm_guess <- grep("^GSM", as.character(unlist(meta[1, ])), value = TRUE)
  if (length(gsm_guess) == 0) stop("❌ No se pudo detectar ID de muestra.")
  idx <- which(meta[1, ] == gsm_guess[1])[1]
  colnames(meta)[idx]
}

meta <- meta[!duplicated(meta[[id_col]]), ]
rownames(meta) <- toupper(gsub("\\.", "", meta[[id_col]]))
colnames(expr) <- toupper(gsub("\\.", "", colnames(expr)))
common_samples <- intersect(colnames(expr), rownames(meta))
if (length(common_samples) < 10) stop("❌ Menos de 10 muestras comunes entre expresión y metadatos.")
expr <- expr[, common_samples]
meta <- meta[common_samples, ]

# === DETECCIÓN LITERAL de variable de grupo (versión definitiva) ===
if ("disease status:ch1" %in% names(meta)) {
  group_var <- "disease status:ch1"
} else if ("subject group:ch1" %in% names(meta)) {
  group_var <- "subject group:ch1"
} else if ("X.disease" %in% names(meta)) {
  group_var <- "X.disease"
} else {
  stop("❌ No se encontró ninguna de las columnas de grupo esperadas.")
}

message("✅ Usando variable de grupo: ", group_var)

group <- as.factor(toupper(trimws(as.character(meta[[group_var]]))))
group <- droplevels(group)

if (length(unique(group)) != 2) {
  message("⚠️ Valores detectados en ", group_var, ": ", paste(levels(group), collapse = ", "))
  stop("❌ La variable de grupo no tiene dos niveles válidos.")
}

# === Preparación para modelado
x <- as.matrix(t(expr))  # Transponer: genes como variables, muestras como observaciones
y <- group

# === LASSO ROBUSTO
message("🔁 Ejecutando LASSO con 50 semillas...")
lasso_genes_list <- list()
for (i in 1:50) {
  set.seed(i)
  lasso_cv <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 5)
  coef_lasso <- coef(lasso_cv, s = "lambda.min")
  genes_lasso_i <- rownames(coef_lasso)[which(coef_lasso[, 1] != 0)]
  lasso_genes_list[[i]] <- genes_lasso_i[genes_lasso_i != "(Intercept)"]
}
# Contar frecuencia de aparición
lasso_all_genes <- unlist(lasso_genes_list)
freq_lasso_table <- table(lasso_all_genes)
freq_lasso_df <- data.frame(gene = names(freq_lasso_table), freq_LASSO = as.integer(freq_lasso_table))
genes_lasso <- freq_lasso_df$gene[freq_lasso_df$freq_LASSO >= 10]  # umbral ajustable
write.csv(freq_lasso_df, file.path(results_dir, "top_genes_lasso.csv"), row.names = FALSE)

# === RANDOM FOREST ROBUSTO
message("🌳 Ejecutando Random Forest con 20 iteraciones...")
rf_genes_list <- list()
for (i in 1:20) {
  set.seed(100 + i)
  rf <- randomForest(x, y, ntree = 500, importance = TRUE)
  rf_imp <- importance(rf, type = 2)
  rf_genes_list[[i]] <- sort(rf_imp[, 1], decreasing = TRUE)[1:50]
}
rf_all_genes <- unlist(lapply(rf_genes_list, names))
rf_freq_table <- table(rf_all_genes)
rf_df <- data.frame(gene = names(rf_freq_table), freq_RF = as.integer(rf_freq_table))
genes_rf <- rf_df$gene[rf_df$freq_RF >= 5]  # umbral ajustable
write.csv(rf_df, file.path(results_dir, "top_genes_rf.csv"), row.names = FALSE)

# === DEGs con SYMBOL válido
deg[[symbol_col]] <- as.character(deg[[symbol_col]])
deg_valid <- deg[!is.na(deg[[symbol_col]]) & !grepl("^ILMN_", deg[[symbol_col]]), ]
if ("adj.P.Val" %in% names(deg_valid)) {
  genes_deg <- deg_valid[[symbol_col]][deg_valid$adj.P.Val < 0.05]
  deg_significativos <- length(genes_deg)
} else {
  genes_deg <- character(0)
  deg_significativos <- 0
}

# === Panel por intersección adaptativa (≥2 métodos)
in_deg   <- toupper(deg[[symbol_col]])
in_lasso <- toupper(genes_lasso)
in_rf    <- toupper(genes_rf)

# Unión total
genes_union <- union(union(in_lasso, in_rf), in_deg)

# Contar en cuántos métodos aparece cada gen
gene_counts <- table(unlist(list(in_lasso, in_rf, in_deg)))
genes_common <- names(gene_counts)[gene_counts >= 2]

# Determinar estrategia
if (length(genes_common) > 0) {
  strategy <- "Genes seleccionados por ≥2 métodos (DEG, LASSO, RF)"
} else if (length(in_lasso) > 0) {
  genes_common <- in_lasso
  strategy <- "Fallback: solo LASSO disponible"
} else if (length(in_rf) > 0) {
  genes_common <- in_rf
  strategy <- "Fallback: solo RF disponible"
} else {
  genes_common <- character(0)
  strategy <- "❌ Sin genes seleccionables"
}

# === Crear mapeo robusto SYMBOL <-> ID
deg$gene <- toupper(trimws(as.character(deg$gene)))
deg[[symbol_col]] <- toupper(trimws(as.character(deg[[symbol_col]])))

symbol_map <- deg %>%
  filter(!is.na(gene), !is.na(!!sym(symbol_col))) %>%
  distinct(gene, .data[[symbol_col]]) %>%
  rename(id_plataforma = gene, SYMBOL = !!sym(symbol_col))

# Usar solo IDs de plataforma para las comparaciones
in_deg_ids   <- toupper(deg$gene[deg$adj.P.Val < 0.05])
in_lasso_ids <- toupper(genes_lasso)
in_rf_ids    <- toupper(genes_rf)

# Unión total y recuento
genes_union_ids <- union(union(in_lasso_ids, in_rf_ids), in_deg_ids)
gene_counts <- table(unlist(list(in_lasso_ids, in_rf_ids, in_deg_ids)))
genes_common_ids <- names(gene_counts)[gene_counts >= 2]

# Determinar estrategia
if (length(genes_common_ids) > 0) {
  strategy <- "Genes seleccionados por ≥2 métodos (DEG, LASSO, RF)"
} else if (length(in_lasso_ids) > 0) {
  genes_common_ids <- in_lasso_ids
  strategy <- "Fallback: solo LASSO disponible"
} else if (length(in_rf_ids) > 0) {
  genes_common_ids <- in_rf_ids
  strategy <- "Fallback: solo RF disponible"
} else {
  genes_common_ids <- character(0)
  strategy <- "❌ Sin genes seleccionables"
}

# === Guardar panel_robusto_validado.csv
panel_df <- data.frame(id_plataforma = genes_common_ids)
panel_df <- left_join(panel_df, symbol_map, by = "id_plataforma")
write.csv(panel_df, file.path(results_dir, "panel_robusto_validado.csv"), row.names = FALSE)

# === Guardar panel_robusto_expandido.csv
expandido <- data.frame(id_plataforma = genes_union_ids)
expandido <- expandido %>%
  mutate(
    in_LASSO = id_plataforma %in% in_lasso_ids,
    in_RF    = id_plataforma %in% in_rf_ids,
    in_DEG   = id_plataforma %in% in_deg_ids
  ) %>%
  left_join(symbol_map, by = "id_plataforma")
write.csv(expandido, file.path(results_dir, "panel_robusto_expandido.csv"), row.names = FALSE)

# === Guardar resumen
summary_lines <- c(
  paste0("Selección de biomarcadores para ", dataset_id),
  paste0("Genes DEG significativos: ", sum(deg$adj.P.Val < 0.05, na.rm = TRUE)),
  paste0("Genes LASSO: ", length(genes_lasso)),
  paste0("Genes RF: ", length(genes_rf)),
  paste0("Genes seleccionados (≥2 métodos): ", length(genes_common_ids)),
  paste0("Estrategia aplicada: ", strategy)
)
cat(paste(summary_lines, collapse = "\n"), "\n")
writeLines(summary_lines, file.path(results_dir, "resumen_biomarcadores.txt"))
