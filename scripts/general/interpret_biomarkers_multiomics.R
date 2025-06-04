#!/usr/bin/env Rscript

# =============================
# INTERPRETACIÓN MULTI-ÓMICA
# =============================
# Autor: Jessica Valdivia
# Objetivo: generar interpretaciones claras y visuales de biomarcadores en transcriptómica y MEG

suppressMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(gprofiler2)
  library(stringr)
  library(ggplot2)
})

# ============
# RUTAS
# ============
base_dir <- "~/TFM_MDD/results"
output_dir <- file.path(base_dir, "interpretation_multiomics")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============
# ARCHIVOS DE ENTRADA
# ============
genes_file <- file.path(base_dir, "transcriptomics/VALIDATION_INTER_ENRICHMENT/inter_dataset_robustness_summary.csv")
features_file <- file.path(base_dir, "ds005356/stats/lasso_selected_features.csv")
regression_file <- file.path(base_dir, "ds005356/stats/regression_results_linear.csv")
meg_dict_file <- file.path(output_dir, "mechanism_dictionary_meg.csv")

# ============
# TRANSCRIPTÓMICA
# ============
cat("🔬 Interpretando biomarcadores transcriptómicos...\n")

genes_df <- read_csv(genes_file)
gene_list <- unique(genes_df$SYMBOL)

# Enriquecimiento funcional con gprofiler2
gostres <- gost(gene_list, organism = "hsapiens", sources = c("GO:BP", "KEGG", "REAC"))

if (!is.null(gostres$result)) {
  enriched <- gostres$result %>%
    select(term_name, source, term_id, intersection_size, query_size) %>%
    mutate(mecanismo = case_when(
      str_detect(term_name, "synap|plasticity") ~ "neuroplasticity",
      str_detect(term_name, "immune|inflamm") ~ "inflammation",
      str_detect(term_name, "neurotransmit|synapse|transmission") ~ "neurotransmission",
      str_detect(term_name, "stress|oxidative") ~ "oxidative_stress",
      TRUE ~ "other"
    ))

  write_csv(enriched, file.path(output_dir, "transcriptomics_enriched_mechanisms.csv"))

  summary_genes <- enriched %>%
    count(mecanismo, name = "n_rutas") %>%
    arrange(desc(n_rutas))

  write_csv(summary_genes, file.path(output_dir, "summary_mecanismos_transcriptomica.csv"))

  cat("✅ Interpretación transcriptómica completada.\n")
} else {
  cat("⚠️ No se encontraron rutas enriquecidas para los genes.\n")
}

# ============
# MEG – INTERPRETACIÓN FEATURES
# ============
cat("🧠 Interpretando biomarcadores funcionales MEG...\n")

features_df <- read_csv(features_file)
regression_df <- read_csv(regression_file)
meg_dict <- read_csv(meg_dict_file)

# Filtrar solo features significativas seleccionadas
features_sign <- regression_df %>%
  filter(feature %in% features_df$feature & padj < 0.05)

features_annotated <- features_sign %>%
  left_join(meg_dict, by = "feature") %>%
  mutate(explanation = paste0("La feature funcional ", feature, " está asociada con ", mechanism, ". ", justification))

write_csv(features_annotated, file.path(output_dir, "meg_features_interpreted.csv"))

cat("✅ Interpretación MEG completada.\n")

# ============
# INTEGRACIÓN MULTI-ÓMICA
# ============
cat("🔗 Integrando mecanismos compartidos...\n")

# Genes por mecanismo
genes_mechs <- enriched %>%
  filter(mecanismo != "other") %>%
  distinct(mecanismo)

# Features por mecanismo
features_mechs <- features_annotated %>%
  filter(!is.na(mechanism)) %>%
  distinct(mechanism)

# Integración
all_mechs <- union(genes_mechs$mecanismo, features_mechs$mechanism)

summary_overlap <- data.frame(
  mechanism = all_mechs,
  genes = all_mechs %in% genes_mechs$mecanismo,
  features = all_mechs %in% features_mechs$mechanism
)

write_csv(summary_overlap, file.path(output_dir, "mecanismos_comunes_multiomica.csv"))

# Visual
p <- ggplot(summary_overlap, aes(x = features, y = genes, label = mechanism)) +
  geom_text(aes(color = mechanism), size = 4) +
  labs(title = "Comparativa de mecanismos entre ómicas", x = "Presente en MEG", y = "Presente en Transcriptómica") +
  theme_minimal()

ggsave(file.path(output_dir, "mecanismos_overlap_plot.pdf"), p, width = 7, height = 5)

cat("✅ Integración finalizada. Resultados guardados en:\n", output_dir, "\n")
