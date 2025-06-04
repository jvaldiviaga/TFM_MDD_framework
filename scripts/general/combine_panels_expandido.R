#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))

# === Rutas base ===
base_dir <- file.path(Sys.getenv("HOME"), "TFM_MDD", "results")
datasets <- c("GSE98793", "GSE44593", "GSE39653")
panel_paths <- file.path(base_dir, datasets, "biomarker_selection", "panel_robusto_expandido.csv")

# === Leer y limpiar paneles ===
leer_panel <- function(path, dataset_id) {
  df <- read_csv(path, show_col_types = FALSE)
  df <- df %>%
    filter(!is.na(SYMBOL), SYMBOL != "") %>%
    select(id_plataforma, SYMBOL, in_LASSO, in_RF, in_DEG) %>%
    rename_with(~ paste0(., "_", dataset_id), -SYMBOL)
    return(df)
}

paneles <- map2(panel_paths, datasets, leer_panel)
panel_combinado <- reduce(paneles, full_join, by = "SYMBOL")

# === Sustituir NAs por FALSE
bool_cols <- grep("^in_", names(panel_combinado), value = TRUE)
panel_combinado[bool_cols] <- panel_combinado[bool_cols] %>%
  mutate(across(everything(), ~ replace_na(., FALSE)))

# === Calcular n_methods
panel_combinado$n_methods <- rowSums(panel_combinado[bool_cols])

# === Añadir logFC desde DEG_results_mapped priorizando carpetas con "with_"
for (dataset in datasets) {
  # Buscar todas las rutas con archivo DEG_results_mapped.csv
  deg_files <- list.files(
    path = file.path(base_dir, dataset, "differential_expression"),
    pattern = "DEG_results_mapped\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  # Prioriza rutas con 'with_' en el nombre
  if (length(deg_files) > 1) {
    deg_files <- deg_files[order(!grepl("with_", deg_files), grepl("no_covariates", deg_files))]
  }
  
  deg_file <- deg_files[1]
  
  if (!is.na(deg_file) && file.exists(deg_file)) {
    deg_df <- read_csv(deg_file, show_col_types = FALSE)
    
    # Detectar columnas aunque estén en minúsculas
    colnames_lower <- tolower(names(deg_df))
    symbol_col <- names(deg_df)[colnames_lower == "symbol"]
    logfc_col  <- names(deg_df)[colnames_lower == "logfc"]
    adjp_col   <- names(deg_df)[colnames_lower == "adj.p.val"]
    
    if (all(c(length(symbol_col), length(logfc_col), length(adjp_col)) == 1)) {
      # Agrupar por símbolo y quedarse con la sonda más significativa (menor adj.P.Val)
      deg_df_resumido <- deg_df %>%
        filter(!is.na(.data[[symbol_col]]), !is.na(.data[[logfc_col]]), !is.na(.data[[adjp_col]])) %>%
        group_by(SYMBOL = .data[[symbol_col]]) %>%
        slice_min(.data[[adjp_col]], n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        select(SYMBOL, logFC = all_of(logfc_col)) %>%
        rename(!!paste0("logFC_", dataset) := logFC)
      
      panel_combinado <- left_join(panel_combinado, deg_df_resumido, by = "SYMBOL")
    } else {
      warning(paste0("⚠️ El archivo DEG de ", dataset, " no contiene columnas SYMBOL, logFC y adj.P.Val válidas."))
    }
  } else {
    warning(paste0("⚠️ No se encontró archivo DEG para ", dataset))
  }
}


# === Guardar resultados
output_path <- file.path(base_dir, "panel_robusto_expandido_unificado.csv")
write_csv(panel_combinado, output_path)
cat("✅ Panel combinado guardado en:", output_path, "\n")

# === Log resumen
log_path <- file.path(base_dir, "log_panel_combinado.txt")

totales <- map2_int(paneles, datasets, ~ nrow(.x))
lineas_log <- c(
  "=== Resumen panel combinado ===",
  paste0("Genes totales por dataset:"),
  paste0("• ", datasets, ": ", totales),
  ""
)

n_2metodos <- sum(panel_combinado$n_methods >= 2)
n_3metodos <- sum(panel_combinado$n_methods == 3)
lineas_log <- c(
  lineas_log,
  paste0("Genes presentes en ≥2 métodos: ", n_2metodos),
  paste0("Genes presentes en los 3 métodos: ", n_3metodos),
  ""
)

top_genes <- panel_combinado %>%
  arrange(desc(n_methods)) %>%
  select(SYMBOL, n_methods) %>%
  head(10)

lineas_log <- c(
  lineas_log,
  "Top 10 genes con mayor n_methods:",
  apply(top_genes, 1, function(row) paste0("• ", row["SYMBOL"], ": ", row["n_methods"]))
)

writeLines(lineas_log, log_path)
cat("📄 Log guardado en:", log_path, "\n")
