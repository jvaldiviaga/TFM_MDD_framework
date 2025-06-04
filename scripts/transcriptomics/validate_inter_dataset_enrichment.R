#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
})

# === Rutas base ===
base_dir <- file.path(Sys.getenv("HOME"), "TFM_MDD")
results_dir <- file.path(base_dir, "results")
output_dir <- file.path(results_dir, "transcriptomics", "VALIDATION_INTER_ENRICHMENT")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# === Cargar todos los resúmenes de validación (por método + dirección logFC)
summary_files <- list.files(results_dir, pattern = "panel_validation_overlap_summary_.*\\.csv$", recursive = TRUE, full.names = TRUE)
all_summaries <- lapply(summary_files, function(f) {
  df <- read_csv(f, show_col_types = FALSE)
  df$dataset <- str_extract(f, "GSE\\d+")
  df$panel_type <- ifelse(str_detect(f, "expandido"), "expandido", "validado")
  return(df)
})
summary_df <- bind_rows(all_summaries)
summary_df <- summary_df %>%
  mutate(SYMBOL = toupper(SYMBOL))

# === Validación por dirección (logFC) desde validation_logfc.csv
validation_files <- list.files(results_dir, pattern = "validation_logfc\\.csv$", recursive = TRUE, full.names = TRUE)
genes_validated <- unique(unlist(lapply(validation_files, function(f) {
  df <- read_csv(f, show_col_types = FALSE)
  df %>% filter(direction_match == TRUE) %>% pull(gene) %>% toupper()
})))

# === Genes del panel principal
main_panel_path <- file.path(results_dir, "GSE98793", "biomarker_selection", "panel_robusto_validado.csv")
main_genes <- if (file.exists(main_panel_path)) {
  read_csv(main_panel_path, show_col_types = FALSE) %>% pull(SYMBOL) %>% toupper()
} else character(0)

# === Enriquecimiento funcional (si existe)
enrichment_dir <- list.dirs(path = file.path(results_dir, "GSE98793", "enrichment"), recursive = TRUE)
enrichment_files <- list.files(enrichment_dir, pattern = ".*_enrichment\\.csv$", full.names = TRUE)
enriched_genes <- unique(unlist(lapply(enrichment_files, function(f) {
  df <- read_csv(f, show_col_types = FALSE)
  if ("geneID" %in% names(df)) unlist(str_split(df$geneID, "/")) else character(0)
}))) %>% toupper()

# === Construir resumen por SYMBOL
summary_global <- summary_df %>%
  group_by(SYMBOL) %>%
  summarise(
    datasets_detected = n_distinct(dataset),
    methods_detected = sum(LASSO + RF + DEG, na.rm = TRUE),
    validated_directions = sum(SYMBOL %in% genes_validated),
    in_main_panel = any(SYMBOL %in% main_genes),
    in_enrichment = any(SYMBOL %in% enriched_genes),
    .groups = "drop"
  ) %>%
  mutate(
    total_score = methods_detected +
                  validated_directions +
                  as.numeric(in_main_panel) +
                  as.numeric(in_enrichment)
  ) %>%
  arrange(desc(total_score))

# === Guardar CSV
write_csv(summary_global, file.path(output_dir, "inter_dataset_robustness_summary.csv"))

# === Generar figura interpretativa
top_genes <- summary_global %>% slice_max(total_score, n = 20)
png(file.path(output_dir, "inter_dataset_robust_genes_plot.png"), width = 1000, height = 600)
ggplot(top_genes, aes(x = reorder(SYMBOL, total_score), y = total_score, fill = in_enrichment)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#2c7fb8", "FALSE" = "#cccccc")) +
  labs(
    title = "Genes robustos entre datasets transcriptómicos",
    subtitle = "Basado en coincidencia de métodos, validación logFC, rutas y panel principal",
    x = "Gen (SYMBOL)", y = "Puntuación total de robustez"
  ) +
  theme_minimal()
dev.off()

# === Log explicativo
log_lines <- c(
  "=== VALIDACIÓN CRUZADA INTER-DATASET DE TRANSCRIPTÓMICA ===",
  paste0("🧬 Genes detectados en ≥1 dataset: ", nrow(summary_global)),
  paste0("🔁 Genes validados por dirección logFC: ", sum(summary_global$validated_directions > 0)),
  paste0("📌 Genes en panel principal (GSE98793): ", sum(summary_global$in_main_panel)),
  paste0("🧠 Genes implicados en rutas GO/KEGG: ", sum(summary_global$in_enrichment)),
  paste0("📄 Resumen global guardado en: inter_dataset_robustness_summary.csv"),
  paste0("🖼️ Imagen de resumen: inter_dataset_robust_genes_plot.png")
)
writeLines(log_lines, file.path(output_dir, "log_validation_inter_dataset.txt"))
cat(paste(log_lines, collapse = "\n"), "\n")
