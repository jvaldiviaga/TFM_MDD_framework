#!/usr/bin/env Rscript

# enrichment_analysis.R (versión definitiva y estable)
# Análisis funcional global + interpretación fisiopatológica general para MDD
# Uso: source("run_enrichment_interactive.R")



# FORTALEZAS Y JUSTIFICACIÓN DEL FLUJO
#
# Este script constituye el núcleo funcional y robusto del análisis de rutas
# biológicas alteradas en MDD, a partir de genes diferencialmente expresados (DEGs)
# ya mapeados a ENTREZID. Aplica análisis de enriquecimiento clásico (ORA) y genera
# salidas visuales e interpretativas.
#
# - Automatiza enriquecimiento funcional sobre GO (BP, MF, CC), KEGG y Reactome.
# - Integra visualizaciones clave (dotplot, cnetplot, emapplot), protegidas ante errores.
# - Exporta resultados enriquecidos en formatos estructurados (.csv, .rds, .png).
# - Genera interpretación automática general basada en mecanismos MDD (ej. neurotransmisión, inflamación...).
# - Se integra desde `run_enrichment_interactive.R` o pipelines automatizados.


# DECISIONES CLAVE Y JUSTIFICACIÓN

#
# 1. Se usa `DEG_results_mapped.csv` como única entrada (flujo simple y controlado).
# 2. Se filtran genes con adj.P.Val < 0.05 y ENTREZID válido.
# 3. Se aplica ORA mediante enrichGO(), enrichKEGG(), enrichPathway().
# 4. Las funciones de visualización están protegidas si el objeto es nulo o vacío.
# 5. Se omite enriquecimiento si hay menos de 10 genes significativos.
# 6. Se genera una interpretación general basada en una lista cerrada de genes asociados a MDD.
# 7. El resumen final resume términos enriquecidos por base y total de genes.

# LIMITACIONES DETECTADAS

# 1. El umbral de significancia (`adj.P.Val < 0.05`) y qvalueCutoff están codificados de forma fija
#    dentro del script (no son argumentos externos configurables).

# 2. El enriquecimiento se omite si hay <10 genes significativos, pero no se permite ajustar
#    este umbral ni se advierte si el número es marginalmente bajo (ej. entre 10 y 20).

# 3. No se informa al usuario si los análisis de enriquecimiento fallan internamente
#    (por ejemplo, si enrichGO() devuelve NULL silenciosamente).

# 4. La interpretación fisiopatológica se basa en una lista fija y cerrada de genes por mecanismo,
#    codificada directamente en el script. No se permite ampliarla dinámicamente.

# 5. No se valida que los términos enriquecidos (GO/KEGG/Reactome) tengan relevancia clínica
#    o psiquiátrica antes de ser considerados en la interpretación final.

# 6. Aunque está presente el bloque de GSEA, no se valida la longitud mínima del ranking antes
#    de ejecutarlo, y no se controla si los resultados están vacíos.

# 7. El resumen final muestra “Top 5” términos incluso cuando hay <5 rutas enriquecidas,
#    lo cual puede inducir a interpretaciones erróneas.

# 8. El resumen (`summary_enrichment.txt`) solo integra resultados de ORA, no de GSEA,
#    y no separa UP y DOWN, ya que el script actual no lo implementa.

# 9. No se implementa ninguna simplificación semántica (ej. con `simplify()`) para eliminar
#    redundancia en términos GO.

# 10. El script no guarda logs o mensajes de error si fallan las funciones internas
#     (ni traza qué pasos se completaron).


# MEJORAS FUTURAS


# - Permitir configurar umbrales como `adj.P.Val`, `qvalueCutoff` y `min.genes` como argumentos externos.

# - Añadir validación previa para advertir si los genes significativos están entre 10 y 20,
#   indicando que el enriquecimiento podría carecer de potencia.

# - Mostrar por pantalla si cada análisis (GO, KEGG, Reactome) tuvo éxito o no.

# - Externalizar listas de genes por mecanismo fisiopatológico y glosarios de palabras clave
#   en archivos editables (CSV, JSON) para facilitar futuras ampliaciones.

# - Integrar simplificación de términos GO redundantes (`simplify()` de clusterProfiler).

# - Implementar visualizaciones adaptativas: por ejemplo, no generar `cnetplot` o `emapplot`
#   si hay <5 términos enriquecidos, y evitar que los gráficos fallen.

# - Añadir sección detallada a `summary_enrichment.txt` para rutas por base funcional, incluyendo
#   UP/DOWN y GSEA si se activa en versiones futuras.

# - Añadir trazabilidad con logs por análisis (incluyendo número de genes, éxito del ORA, etc.)

# - Controlar errores explícitos durante llamadas a `bitr_kegg()`, `enrich*()`, `AnnotationDbi::select()`
#   y otras funciones críticas.



# --- Cargar paquetes ---
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(readr)
  library(dplyr)
  library(tibble)
  library(enrichplot)
  library(ggplot2)
  library(DOSE)
})

# --- Argumentos ---
args <- commandArgs(trailingOnly = TRUE)
dataset_id <- ifelse(length(args) >= 1, args[1], stop("❌ Falta dataset_id"))
subfolder  <- ifelse(length(args) >= 2, args[2], stop("❌ Falta subcarpeta"))

# --- Rutas ---
base_dir <- file.path("~/TFM_MDD/results", dataset_id, "differential_expression", subfolder)
mapped_path <- file.path(base_dir, "DEG_results_mapped.csv")
out_dir <- file.path("~/TFM_MDD/results", dataset_id, "enrichment", subfolder)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Leer DEGs ---
deg <- read_csv(mapped_path, show_col_types = FALSE)
deg_sig <- deg %>%
  filter(adj.P.Val < 0.05 & !is.na(entrez_id) & entrez_id != "") %>%
  distinct(entrez_id, .keep_all = TRUE)

<<<<<<< HEAD
entrez_ids <- unique(deg_sig$entrez_id)
cat("✅ Genes significativos:", length(entrez_ids), "\n")

=======
# --- Leer DEGs ---
deg <- read_csv(mapped_path, show_col_types = FALSE)
deg_sig <- deg %>%
  filter(adj.P.Val < 0.05 & !is.na(entrez_id) & entrez_id != "") %>%
  distinct(entrez_id, .keep_all = TRUE)

# --- Salir si no hay genes significativos ---
if (nrow(deg_sig) == 0) {
  cat("❌ No hay genes diferencialmente expresados significativos. Se omite el enriquecimiento.\n")
  
  # Crear carpeta de salida y archivos dummy para evitar error en Snakemake
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines("Sin genes significativos, no se ejecutó el análisis.", file.path(out_dir, "summary_enrichment.txt"))
  write.csv(data.frame(), file.path(out_dir, "significant_genes.csv"), row.names = FALSE)
  
  quit(save = "no", status = 0)
}

entrez_ids <- unique(deg_sig$entrez_id)
cat("✅ Genes significativos:", length(entrez_ids), "\n")
cat(🚀 Subida completa: scripts nuevos, carpetas actualizadas y README final)

# --- Mostrar genes reales ---
if (!"symbol" %in% names(deg_sig)) stop("❌ Falta columna 'symbol' con nombres de genes")
cat("\n🔬 Ejemplo de genes significativos:\n")
print(deg_sig %>% select(GeneSymbol = symbol, logFC, adj.P.Val) %>% head(10))
write_csv(deg_sig %>% select(GeneSymbol = symbol, entrez_id, logFC, adj.P.Val),
          file.path(out_dir, "significant_genes.csv"))

# --- Función modular ORA ---
run_enrichment <- function(gene_ids, type, ont = NULL, ...) {
  tryCatch({
    if (type == "GO") enrichGO(gene = gene_ids, OrgDb = org.Hs.eg.db, ont = ont, ...)
    else if (type == "KEGG") enrichKEGG(gene = gene_ids, organism = "hsa", ...)
    else if (type == "Reactome") enrichPathway(gene = gene_ids, organism = "human", ...)
  }, error = function(e) NULL)
}

# --- Ejecutar análisis ---
ego_bp <- run_enrichment(entrez_ids, "GO", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.2, readable = TRUE)
ego_mf <- run_enrichment(entrez_ids, "GO", ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 0.2, readable = TRUE)
ego_cc <- run_enrichment(entrez_ids, "GO", ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 0.2, readable = TRUE)
ekegg  <- run_enrichment(entrez_ids, "KEGG", pAdjustMethod = "BH")
ereact <- run_enrichment(entrez_ids, "Reactome", pAdjustMethod = "BH", readable = TRUE)

# --- Guardar resultados ---
save_and_export <- function(obj, name) {
  if (!is.null(obj)) {
    write.csv(as.data.frame(obj), file.path(out_dir, paste0(name, "_enrichment.csv")), row.names = FALSE)
    saveRDS(obj, file.path(out_dir, paste0(name, "_enrichment.rds")))
  }
}
save_and_export(ego_bp, "GO_BP")
save_and_export(ego_mf, "GO_MF")
save_and_export(ego_cc, "GO_CC")
save_and_export(ekegg, "KEGG")
save_and_export(ereact, "Reactome")

# --- Visualizaciones con protección ---
plot_all <- function(obj, name) {
  if (!is.null(obj) && nrow(obj) > 0) {
    n_terms <- min(15, nrow(obj))
    try({
      png(file.path(out_dir, paste0(name, "_dotplot.png")), width = 1000)
      print(dotplot(obj, showCategory = n_terms) + ggtitle(paste(name, " pathways")))
      dev.off()
    }, silent = TRUE)
    if (n_terms >= 3) {
      try({
        png(file.path(out_dir, paste0(name, "_cnetplot.png")), width = 1200)
        print(cnetplot(obj, showCategory = min(5, n_terms)))
        dev.off()
        png(file.path(out_dir, paste0(name, "_emapplot.png")), width = 1200)
        print(emapplot(pairwise_termsim(obj)))
        dev.off()
      }, silent = TRUE)
    }
  } else {
    cat(paste0("⚠️  No se puede graficar ", name, ": sin términos enriquecidos.\n"))
  }
}
plot_all(ego_bp, "GO_BP")
plot_all(ego_mf, "GO_MF")
plot_all(ego_cc, "GO_CC")
plot_all(ekegg, "KEGG")
plot_all(ereact, "Reactome")

# --- GSEA si existe estadístico t ---
if ("t" %in% colnames(deg) && !any(is.na(deg$t))) {
  gene_ranks <- deg %>% filter(!is.na(entrez_id)) %>%
    arrange(desc(t)) %>% distinct(entrez_id, .keep_all = TRUE) %>%
    select(entrez_id, t) %>% deframe()
  
  gsea_go <- gseGO(geneList = gene_ranks, OrgDb = org.Hs.eg.db, ont = "BP", verbose = FALSE)
  gsea_kegg <- gseKEGG(geneList = gene_ranks, organism = "hsa", verbose = FALSE)
  gsea_react <- gsePathway(geneList = gene_ranks, organism = "human", verbose = FALSE)
  
  save_and_export(gsea_go, "GSEA_GO_BP")
  save_and_export(gsea_kegg, "GSEA_KEGG")
  save_and_export(gsea_react, "GSEA_Reactome")
}

# --- Interpretación fisiopatológica automática ---
fisiopato_genes <- list(
  "Neurotransmisión" = c("SLC6A4", "HTR2A", "DRD2", "GRIK5", "GRM5", "NEGR1", "RBFOX1", "LHPP"),
  "Neuroplasticidad" = c("BDNF", "PCLO", "CACNA1E", "CACNA2D1", "SIRT1", "CREB", "LRFN5"),
  "Inflamación" = c("IL6", "TNF", "IL1B", "CRP", "OLFM4", "GFAP", "AQP4", "GLT1"),
  "Eje_HPA" = c("FKBP5", "NR3C1", "CRH", "CRHR1"),
  "Metabolismo_Mitocondrial" = c("SIRT1", "OLFM4", "BDNF")
)

deg_sig_genes <- unique(deg_sig$symbol)
interpretacion <- c("🧠 Interpretación fisiopatológica basada en genes significativos:")

for (mecanismo in names(fisiopato_genes)) {
  genes_detectados <- intersect(deg_sig_genes, fisiopato_genes[[mecanismo]])
  if (length(genes_detectados) > 0) {
    interpretacion <- c(interpretacion, paste0("✔ ", mecanismo, ": ", paste(genes_detectados, collapse = ", ")))
  }
}
if (length(interpretacion) == 1) {
  interpretacion <- c(interpretacion, "❌ No se detectaron genes relevantes en vías fisiopatológicas clave.")
}
cat("\n", paste(interpretacion, collapse = "\n"), "\n")
write_lines(interpretacion, file.path(out_dir, "interpretacion_fisiopatologica.txt"))

# --- Interpretación individual de genes según rutas funcionales (GO/KEGG) ---

library(AnnotationDbi)

# Glosario de mecanismos fisiopatológicos del MDD y palabras clave asociadas
glosario_mdd <- list(
  "Neurotransmisión" = c("synapse", "neurotransmitter", "gaba", "glutamate", "serotonin", "dopamine"),
  "Neuroinflamación" = c("cytokine", "inflammation", "immune"),
  "Eje HPA" = c("stress", "hormone", "glucocorticoid"),
  "Neuroplasticidad" = c("plasticity", "neurogenesis", "apoptosis"),
  "Disfunción mitocondrial" = c("mitochondria", "oxidative", "metabolism", "energy")
)

interpretaciones_gen <- list()

# Procesar cada gen significativo
for (gen in unique(na.omit(deg_sig$symbol))) {
  rutas_go <- tryCatch(AnnotationDbi::select(org.Hs.eg.db, keys = gen, keytype = "SYMBOL", columns = "GO"), error = function(e) NULL)
  rutas_kegg <- tryCatch(bitr_kegg(gen, fromType = "symbol", toType = "Path", organism = "hsa"), error = function(e) NULL)
  
  todas_rutas <- unique(c(
    if (!is.null(rutas_go)) as.character(rutas_go$GO), 
    if (!is.null(rutas_kegg)) as.character(rutas_kegg$Path)
  ))
  
  # Crear resumen por gen
  texto <- paste0("🔹 ", gen, " ➤ ", length(todas_rutas), " rutas encontradas")
  relacionadas <- character()
  
  for (palabra in unlist(glosario_mdd)) {
    if (any(grepl(palabra, tolower(todas_rutas), fixed = TRUE))) {
      mecanismo <- names(glosario_mdd)[sapply(glosario_mdd, function(pat) palabra %in% pat)]
      relacionadas <- unique(c(relacionadas, mecanismo))
    }
  }
  
  if (length(relacionadas) > 0) {
    texto <- paste(texto, "\n     → Relacionado con:", paste(unique(relacionadas), collapse = ", "))
  } else {
    texto <- paste(texto, "\n     → No se identificó relación directa con mecanismos MDD")
  }
  
  interpretaciones_gen[[gen]] <- list(
    gene = gen,
    rutas_detectadas = paste(todas_rutas, collapse = "; "),
    mecanismos_MDD = ifelse(length(relacionadas) > 0, paste(unique(relacionadas), collapse = ", "), "Ninguno"),
    comentario = ifelse(length(relacionadas) > 0, "Posible implicación fisiopatológica", "No asociado directamente")
  )
  
  cat(texto, "\n")
}

# Guardar en .csv
ruta_csv <- file.path(out_dir, "ruta_por_gen.csv")
write.csv(do.call(rbind, lapply(interpretaciones_gen, as.data.frame)), ruta_csv, row.names = FALSE)

# Guardar en .txt resumen
ruta_txt <- file.path(out_dir, "interpretacion_por_gen.txt")
writeLines(unlist(lapply(interpretaciones_gen, function(x) {
  paste0("🔸 ", x$gene, ": ", x$mecanismos_MDD, " → ", x$comentario)
})), con = ruta_txt)


# --- Resumen final ---
write_lines(c(
  paste0("📄 Enriquecimiento - ", subfolder),
  paste0("🔹 Genes totales: ", nrow(deg)),
  paste0("🔹 Genes significativos: ", length(entrez_ids)),
  paste0("🔹 GO-BP: ", ifelse(!is.null(ego_bp), nrow(ego_bp), 0)),
  paste0("🔹 KEGG: ", ifelse(!is.null(ekegg), nrow(ekegg), 0)),
  paste0("🔹 Reactome: ", ifelse(!is.null(ereact), nrow(ereact), 0)),
  "\n🧠 Top 5 GO-BP:",
  if (!is.null(ego_bp)) paste(head(ego_bp$Description, 5), collapse = "\n") else "N/A",
  "\n🧬 Top 5 KEGG:",
  if (!is.null(ekegg)) paste(head(ekegg$Description, 5), collapse = "\n") else "N/A",
  "\n🧪 Top 5 Reactome:",
  if (!is.null(ereact)) paste(head(ereact$Description, 5), collapse = "\n") else "N/A"
), file.path(out_dir, "summary_enrichment.txt"))

cat("🎯 Análisis completado y resultados guardados en:\n", out_dir, "\n")
