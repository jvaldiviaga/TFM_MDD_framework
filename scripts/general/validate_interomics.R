#!/usr/bin/env Rscript

# ===============================================
# generate_biomarker_report.R
# Genera informe de biomarcadores multiómicos del TFM
# Incluye interpretación funcional y fisiopatológica
# ===============================================

suppressMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(DOSE)
  library(enrichplot)
  library(ggplot2)
  library(knitr)
  library(rmarkdown)
})

# ====================
# 1. CARGA DE ARCHIVOS
# ====================

cat("📥 Cargando datos reales...\n")

genes <- read_csv("~/TFM_MDD/results/transcriptomics/VALIDATION_INTER_ENRICHMENT/inter_dataset_robustness_summary.csv", show_col_types = FALSE)
meg_lasso <- read_csv("~/TFM_MDD/results/neuroimaging/lasso_selected_features.csv", show_col_types = FALSE)
meg_regression <- read_csv("~/TFM_MDD/results/neuroimaging/regression_results_linear.csv", show_col_types = FALSE)

# ============
# 2. FILTRADO
# ============

top_genes <- genes %>%
  arrange(desc(total_score)) %>%
  slice_head(n = 5) %>%
  pull(SYMBOL)

significant_meg <- meg_regression %>%
  filter(P_value < 0.05, grepl("rwp|p300|gfp", Feature, ignore.case = TRUE)) %>%
  arrange(P_value)

# ============================
# 3. ANÁLISIS FUNCIONAL GENES
# ============================

gene_ids <- bitr(top_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
go <- enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
kegg <- enrichKEGG(gene = gene_ids$ENTREZID, organism = 'hsa')

# ================================
# 4. CLASIFICACIÓN FUNCIONAL (GO)
# ================================

# Diccionario simple: palabra clave -> mecanismo
mecanismos <- list(
  "Neurotransmisión" = c("synapse", "neurotransmission", "GABA", "dopamine", "glutamate"),
  "Inflamación / inmunidad" = c("immune", "inflammatory", "cytokine", "response to virus", "interleukin"),
  "Neuroplasticidad" = c("plasticity", "axon", "neuron projection", "learning", "memory"),
  "Mitocondrial / estrés oxidativo" = c("mitochondrion", "oxidative", "respiratory chain"),
  "Circadiano" = c("circadian", "clock", "sleep"),
  "Señalización intracelular" = c("MAPK", "signaling pathway", "PI3K", "kinase"),
  "Procesamiento cognitivo / MEG" = c("p300", "rwp", "attention", "GFP")
)

mapear_mecanismo <- function(texto) {
  texto <- tolower(texto)
  for (mec in names(mecanismos)) {
    for (patron in mecanismos[[mec]]) {
      if (grepl(patron, texto)) return(mec)
    }
  }
  return("Otro")
}

# Clasificamos genes y rutas
go_df <- as.data.frame(go) %>% select(Description, pvalue) %>%
  mutate(Mecanismo = sapply(Description, mapear_mecanismo))

# ================================
# 5. INTERPRETACIÓN DE MEG
# ================================

interpretacion_meg <- significant_meg %>%
  mutate(Descripcion = case_when(
    grepl("rwp", Feature, ignore.case = TRUE) ~ "Duración de la respuesta prolongada (RWP), relacionada con procesamiento sostenido",
    grepl("p300", Feature, ignore.case = TRUE) ~ "Latencia del componente P300, vinculada a procesamiento atencional",
    grepl("gfp", Feature, ignore.case = TRUE) ~ "Actividad global del campo eléctrico (GFP)",
    TRUE ~ "No interpretado"
  ))

# =========================
# 6. GENERAR TEXTO RESUMEN
# =========================

cat("\n📝 RESUMEN POR PANTALLA:\n")
cat("Top genes implicados:\n")
print(top_genes)

cat("\nFeatures MEG significativas:\n")
print(interpretacion_meg$Feature)

# ========================
# 7. GENERACIÓN DE INFORME
# ========================

output_dir <- "~/TFM_MDD/results/interpretation_multiomics"
rmd_file <- file.path(output_dir, "biomarker_summary_report.Rmd")
pdf_file <- file.path(output_dir, "biomarker_summary_report.pdf")

# Creamos archivo Rmd
rmd_contenido <- c(
  "---",
  "title: \"Resumen de biomarcadores multiómicos\"",
  "output: pdf_document",
  "---",
  "",
  "# Genes principales identificados",
  "",
  paste0("- **", g, "**") %>% purrr::map_chr(~ paste0(.x, ": participa en rutas GO relacionadas con ", 
                                                      paste(go_df$Mecanismo[grepl(.x, tolower(go_df$Description))][1], collapse = ", "), ".")) %>%
    paste(collapse = "\n"),
  "",
  "# Rutas GO enriquecidas",
  "",
  paste0("- ", go_df$Description, " (", go_df$Mecanismo, ")") %>%
    paste(collapse = "\n"),
  "",
  "# Características funcionales MEG significativas",
  "",
  paste0("## ", interpretacion_meg$Feature, "\n\n", interpretacion_meg$Descripcion) %>%
    paste(collapse = "\n\n")
)

writeLines(rmd_contenido, rmd_file)

# Compila PDF
cat("📄 Generando PDF...\n")
tryCatch({
  rmarkdown::render(rmd_file, output_file = pdf_file, quiet = TRUE)
}, error = function(e) {
  cat("❌ No se pudo compilar el PDF:", e$message, "\n")
})

cat(paste0("✅ Informe generado en: ", pdf_file, "\n"))
