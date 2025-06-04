#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tools)
  library(ggplot2)
  library(ggrepel)
})

# === Rutas base ===
base_dir <- file.path(Sys.getenv("HOME"), "TFM_MDD")
results_dir <- file.path(base_dir, "results")

# === Buscar archivo de validación
validation_file <- list.files(
  path = results_dir,
  pattern = "validation_logfc\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)[1]

if (is.na(validation_file)) stop("❌ No se encontró ningún archivo validation_logfc.csv")
message("📄 Archivo de validación encontrado: ", validation_file)
validation <- read_csv(validation_file, show_col_types = FALSE)

# Verificación
if (!all(c("gene", "direction_match") %in% colnames(validation))) {
  stop("❌ El archivo debe tener columnas 'gene' y 'direction_match'")
}
validated <- validation %>% filter(direction_match == TRUE) %>% mutate(gene = toupper(gene))

# === Detectar paneles robustos y validados
panel_files <- list.files(path = results_dir, pattern = "panel_robusto_(expandido|validado)\\.csv$",
                          recursive = TRUE, full.names = TRUE)
if (length(panel_files) == 0) stop("❌ No se encontraron paneles robustos.")

log_lines <- c("=== VALIDACIÓN PANEL + DIRECCIÓN DE logFC (direction_match == TRUE) ===")
resumen_global <- list()

for (panel_path in panel_files) {
  dataset <- basename(dirname(dirname(panel_path)))
  tipo_panel <- ifelse(grepl("expandido", panel_path), "expandido", "validado")
  log_lines <- c(log_lines, paste0("\n📂 Dataset: ", dataset, " [", tipo_panel, "]"))
  
  panel <- read_csv(panel_path, show_col_types = FALSE)
  if (!("id_plataforma" %in% names(panel)) || !("SYMBOL" %in% names(panel))) {
    log_lines <- c(log_lines, "⚠️ Faltan columnas necesarias ('id_plataforma', 'SYMBOL').")
    next
  }
  
  panel <- panel %>% mutate(id_plataforma = toupper(id_plataforma))
  matched <- panel %>% filter(id_plataforma %in% validated$gene)
  
  log_lines <- c(log_lines, paste0("🔢 Coincidencias validadas por ID de plataforma: ", nrow(matched)))
  
  if (nrow(matched) == 0) {
    log_lines <- c(log_lines, "⚠️ No hay coincidencias validadas.")
    next
  }
  
  resumen <- matched %>%
    group_by(SYMBOL) %>%
    summarise(
      n_sondas      = n(),
      n_validated   = sum(id_plataforma %in% validated$gene),
      LASSO         = sum(ifelse("in_LASSO" %in% names(.), in_LASSO, 0)),
      RF            = sum(ifelse("in_RF"    %in% names(.), in_RF, 0)),
      DEG           = sum(ifelse("in_DEG"   %in% names(.), in_DEG, 0)),
      .groups = "drop"
    ) %>%
    mutate(total_evidence = n_validated + LASSO + RF + DEG) %>%
    arrange(desc(total_evidence))
  
  # Carpeta de salida por dataset
  out_dir <- file.path(results_dir, dataset, "validate_panel_overlap")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  out_csv <- file.path(out_dir, paste0("panel_validation_overlap_summary_", tipo_panel, ".csv"))
  write_csv(resumen, out_csv)
  log_lines <- c(log_lines, paste0("✅ Resumen guardado: ", out_csv))
  
  resumen$dataset <- dataset
  resumen$panel_type <- tipo_panel
  resumen_global[[length(resumen_global)+1]] <- resumen
  
  # === Figura interpretativa
  fig_path <- file.path(out_dir, paste0("panel_validation_overlap_plot_", tipo_panel, ".png"))
  top_genes <- resumen %>% arrange(desc(total_evidence)) %>% head(min(15, nrow(resumen)))
  
  png(fig_path, width = 1000, height = 700)
  print(
    ggplot(top_genes, aes(x = n_validated, y = LASSO + RF + DEG, label = SYMBOL, color = total_evidence)) +
      geom_point(size = 4) +
      ggrepel::geom_text_repel(max.overlaps = 20) +
      scale_color_gradient(low = "mediumorchid4", high = "violetred1") +
      labs(
        title = paste("Validación cruzada por gen -", dataset, tipo_panel),
        x = "# de sondas validadas (direction_match == TRUE)",
        y = "# de métodos que lo detectaron (LASSO + RF + DEG)"
      ) +
      theme_minimal()
  )
  dev.off()
  log_lines <- c(log_lines, paste0("🖼️ Figura generada: ", fig_path))
}

# === Guardar log
log_path <- file.path(results_dir, "log_panel_validation_overlap_full.txt")
writeLines(log_lines, log_path)
cat(paste(log_lines, collapse = "\n"), "\n")

# === Guardar resumen global
if (length(resumen_global) > 0) {
  resumen_total <- bind_rows(resumen_global)
  out_global <- file.path(results_dir, "transcriptomics", "panel_overlap_validation_summary.csv")
  dir.create(dirname(out_global), recursive = TRUE, showWarnings = FALSE)
  write_csv(resumen_total, out_global)
}

