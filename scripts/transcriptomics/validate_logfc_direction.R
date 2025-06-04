#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tools)
})

# === Directorios base
base_dir <- file.path(Sys.getenv("HOME"), "TFM_MDD", "results")
all_datasets <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)
deg_paths <- list.files(base_dir, pattern = "DEG_results_mapped\\.csv$", recursive = TRUE, full.names = TRUE)

# === Cargar todos los DEG y contar significativos
deg_list <- list()
sig_counts <- c()

for (path in deg_paths) {
  deg <- tryCatch(read_csv(path, show_col_types = FALSE), error = function(e) return(NULL))
  if (!is.null(deg) && "adj.P.Val" %in% colnames(deg)) {
    deg <- deg %>% filter(!is.na(logFC), !is.na(gene), !is.na(adj.P.Val))
    dataset_id <- strsplit(path, "/")[[1]]
    dataset_id <- dataset_id[which(dataset_id %in% all_datasets)][1]
    deg_list[[dataset_id]] <- deg
    sig_counts[dataset_id] <- sum(deg$adj.P.Val < 0.05, na.rm = TRUE)
  }
}

# === Detectar dataset principal: mayor número de DEG significativos
main_dataset <- names(which.max(sig_counts))
cat("📌 Dataset principal detectado:", main_dataset, "\n")

# === Comparar con todos los demás
secondary_datasets <- setdiff(names(deg_list), main_dataset)
main_deg <- deg_list[[main_dataset]] %>%
  filter(adj.P.Val < 0.05) %>%
  select(gene, SYMBOL = symbol, logFC_main = logFC)

log_summary <- c()

for (sec in secondary_datasets) {
  sec_deg <- deg_list[[sec]] %>% select(gene, logFC_val = logFC, symbol)
  
  # Cruce por ID de plataforma
  merged <- inner_join(main_deg, sec_deg, by = "gene")
  if (nrow(merged) == 0) {
    warning(paste("❌ No hay genes en común entre", main_dataset, "y", sec))
    next
  }
  
  merged <- merged %>%
    mutate(
      Same_Direction = ifelse(sign(logFC_main) == sign(logFC_val), "✅", "❌"),
      abs_diff = abs(logFC_main - logFC_val)
    )
  
  # Resumen
  n_total <- nrow(merged)
  n_match <- sum(merged$Same_Direction == "✅")
  percent_match <- round(100 * n_match / n_total, 1)
  
  cat("\n🧪 Comparando", sec, "vs", main_dataset, ":\n")
  cat("   Genes comparables:", n_total, "\n")
  cat("   Coincidencia dirección:", n_match, "(", percent_match, "% )\n")
  
  log_summary <- c(log_summary,
                   paste0("=== ", sec, " vs ", main_dataset, " ==="),
                   paste0("Genes comparables: ", n_total),
                   paste0("Coincidencia dirección: ", n_match, " (", percent_match, "%)"),
                   ""
  )
  
  # Top 10 coincidencias
  top_same <- merged %>%
    filter(Same_Direction == "✅") %>%
    arrange(desc(abs(logFC_main))) %>%
    select(SYMBOL, logFC_main, logFC_val) %>%
    head(10)
  
  cat("🔝 Top 10 genes coincidentes en dirección:\n")
  print(top_same)
  log_summary <- c(log_summary, "Top 10 coincidentes:", capture.output(print(top_same)), "")
  
  # Top 10 discordancias
  top_diff <- merged %>%
    filter(Same_Direction == "❌") %>%
    arrange(desc(abs_diff)) %>%
    select(SYMBOL, logFC_main, logFC_val) %>%
    head(10)
  
  cat("⚠️ Top 10 genes con dirección discordante:\n")
  print(top_diff)
  log_summary <- c(log_summary, "Top 10 discordantes:", capture.output(print(top_diff)), "")
  
  # Guardar CSV individual
  out_dir <- file.path(base_dir, sec, "validation_from", main_dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(merged, file.path(out_dir, "validation_logfc_comparison.csv"))
  
  # Gráfico scatter
  png(file.path(out_dir, "logFC_comparison.png"), width = 1000)
  ggplot(merged, aes(logFC_main, logFC_val, color = Same_Direction)) +
    geom_point(alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_minimal() +
    labs(title = paste("Comparación logFC:", main_dataset, "vs", sec),
         x = paste("logFC -", main_dataset),
         y = paste("logFC -", sec)) +
    scale_color_manual(values = c("✅" = "blue", "❌" = "red")) +
    theme(legend.position = "bottom")
  dev.off()
}

# === Guardar log resumen global
log_path <- file.path(base_dir, "log_validation_summary.txt")
writeLines(log_summary, log_path)
cat("\n📄 Log resumen guardado en:", log_path, "\n")
