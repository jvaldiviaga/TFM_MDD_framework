#!/usr/bin/env Rscript
# ==========================================================
# enrichment_analysis.R (versión definitiva)
# Análisis funcional + interpretación automática para MDD
# ==========================================================

# ============================================================

# FORTALEZAS Y JUSTIFICACIÓN DEL FLUJO

# Este script constituye el núcleo del análisis funcional e interpretación automática
# del pipeline transcriptómico. Ha sido diseñado para convertir conjuntos de genes
# diferencialmente expresados (DEGs) en hipótesis fisiopatológicas interpretables
# y en evidencia sobre mecanismos biológicos subyacentes en la depresión mayor.

# - Automatiza análisis de enriquecimiento por ORA sobre GO (BP, MF, CC), KEGG y Reactome.
# - Ejecuta enriquecimientos por separado para genes UP y DOWN, diferenciando mecanismos activados e inhibidos.
# - Integra visualizaciones clave (dotplot, cnetplot, emapplot) protegidas ante errores o resultados vacíos.
# - Ejecuta GSEA si existe estadístico t, mejorando sensibilidad frente a análisis tradicionales.
# - Exporta todos los resultados en formatos estructurados (CSV, RDS, PNG).
# - Genera interpretación automática basada en genes significativos vs listas fisiopatológicas clave (MDD).
# - Cruza cada gen con rutas funcionales y detecta relaciones con mecanismos como neurotransmisión, neuroinflamación, eje HPA, etc.
# - Se integra de forma directa desde `run_enrichment_interactive.R` o desde un pipeline automatizado.

# ------------------------------------------------------------
# DECISIONES CLAVE Y JUSTIFICACIÓN
# ------------------------------------------------------------

# 1. Se usan genes con adj.P.Val < 0.05 y ENTREZID válido como entrada para enriquecimiento funcional.
# 2. Se detectan los DEGs mapeados automáticamente, eliminando duplicados.
# 3. La función `run_enrichment()` modulariza GO, KEGG y Reactome con control de errores interno.
# 4. El enriquecimiento UP/DOWN separa funciones activadas/inhibidas para una lectura más fina.
# 5. `plot_all()` protege cada visualización ante errores, generando solo si hay al menos 3 términos.
# 6. El ranking de GSEA se basa en el estadístico t cuando está disponible, generando `gseGO`, `gseKEGG`, `gsePathway`.
# 7. Se implementa interpretación fisiopatológica automática cruzando genes con mecanismos clave del MDD.
# 8. Se analizan las rutas funcionales por gen (GO, KEGG) y se vinculan con glosarios fisiopatológicos.
# 9. Se genera un resumen final con genes totales, DEG significativos y top 5 rutas por cada base.
# 10. Toda la información se guarda en formato estructurado y legible por humanos.

# ------------------------------------------------------------
# LIMITACIONES DETECTADAS
# ------------------------------------------------------------

# 1. El umbral de significancia (`adj.P.Val < 0.05`) y qvalueCutoff están codificados de forma fija.

# 2. No hay control de mínimo de genes necesarios para ejecutar análisis ni advertencias si n < 10.
#    -Si el conjunto de genes de entrada es demasiado pequeño (por ejemplo, <10),
#    el análisis funcional pierde potencia estadística y puede devolver resultados inestables,
#    poco reproducibles o vacíos.
#    Esto afecta tanto a ORA como a GSEA, ya que no se pueden calcular distribuciones esperadas
#    de forma robusta ni representar rutas relevantes con solo 2–3 genes.
#    - Posible mejora: implementar un umbral mínimo configurable, emitir advertencias
#    explícitas al usuario y evitar generar resultados engañosos cuando n es insuficiente.

# 3. No se informa al usuario si los análisis de enriquecimiento fallan o devuelven 0 términos.

# 4. La interpretación fisiopatológica se basa en una lista cerrada de genes y un glosario de palabras clave.

# 5. No se valida que los términos GO/KEGG encontrados sean realmente relevantes antes de asignar mecanismos.

# 6. El GSEA se lanza sin verificar la longitud mínima del ranking ni controla si los resultados están vacíos.

# 7. El resumen final indica “Top 5” incluso si no hay suficientes términos para mostrar.

# 8. No se integran aún los resultados de GSEA ni UP/DOWN en el `summary_enrichment.txt`.

# ------------------------------------------------------------
# MEJORAS FUTURAS
# ------------------------------------------------------------

# - Permitir configurar el umbral de adj.P.Val y qvalueCutoff como argumentos externos.
# - Añadir validación previa al análisis (e.g., mínimo de genes, existencia de resultados).
# - Informar explícitamente qué análisis fueron exitosos o fallaron (ORA, GSEA, UP, DOWN).
# - Externalizar las listas de genes fisiopatológicos y glosarios como archivos editables.
# - Implementar visualizaciones más flexibles según número de rutas enriquecidas.
# - Añadir sección detallada en `summary_enrichment.txt` para UP, DOWN y GSEA.
# - Mejorar trazabilidad con logs por análisis y verificación de integridad del input.
# - Añadir control de warnings o errores durante las llamadas a `bitr_kegg` y `AnnotationDbi`.

# ============================================================

# --- Cargar paquetes necesarios ---
suppressPackageStartupMessages({
  library(clusterProfiler)  # Paquete principal para ORA y GSEA
  library(org.Hs.eg.db)     # Anotaciones genómicas para humano
  library(ReactomePA)       # Enriquecimiento sobre la base de Reactome
  library(readr)            # Lectura de archivos CSV
  library(dplyr)            # Manipulación de datos
  library(tibble)           # Compatibilidad con dataframes modernos
  library(enrichplot)       # Visualizaciones enriquecimiento
  library(ggplot2)          # Gráficos generales
  library(DOSE)             # Requerido por clusterProfiler para análisis de enfermedad
})

# --- Leer argumentos del script ---
args <- commandArgs(trailingOnly = TRUE)
dataset_id <- ifelse(length(args) >= 1, args[1], stop("❌ Falta dataset_id"))     # ID del dataset
subfolder  <- ifelse(length(args) >= 2, args[2], stop("❌ Falta subcarpeta"))     # Subcarpeta (covariables, modelo, etc.)

# --- Definir rutas de trabajo ---
base_dir <- file.path("~/TFM_MDD/results", dataset_id, "differential_expression", subfolder)
mapped_path <- file.path(base_dir, "DEG_results_mapped.csv")  # Input con ENTREZID y SYMBOL
out_dir <- file.path("~/TFM_MDD/results", dataset_id, "enrichment", subfolder)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)    # Crear directorio si no existe

# --- Leer DEGs mapeados ---
deg <- read_csv(mapped_path, show_col_types = FALSE)

# Filtrar genes con p.adj < 0.05, con ENTREZID válido y sin duplicados
deg_sig <- deg %>%
  filter(adj.P.Val < 0.05 & !is.na(entrez_id) & entrez_id != "") %>%
  distinct(entrez_id, .keep_all = TRUE)

# Extraer IDs únicos para el enriquecimiento
entrez_ids <- unique(deg_sig$entrez_id)
cat("✅ Genes significativos:", length(entrez_ids), "\n")

# --- Mostrar ejemplo de genes significativos ---
# Requiere columna 'symbol' para interpretación legible
if (!"symbol" %in% names(deg_sig)) stop("❌ Falta columna 'symbol' con nombres de genes")

cat("\n🔬 Ejemplo de genes significativos:\n")
print(deg_sig %>% select(GeneSymbol = symbol, logFC, adj.P.Val) %>% head(10))

# Exportar los genes significativos para revisión o downstream
write_csv(deg_sig %>% select(GeneSymbol = symbol, entrez_id, logFC, adj.P.Val),
          file.path(out_dir, "significant_genes.csv"))

# --- Función modular para ORA ---
run_enrichment <- function(gene_ids, type, ont = NULL, ...) {
  tryCatch({
    if (type == "GO") enrichGO(gene = gene_ids, OrgDb = org.Hs.eg.db, ont = ont, ...)
    else if (type == "KEGG") enrichKEGG(gene = gene_ids, organism = "hsa", ...)
    else if (type == "Reactome") enrichPathway(gene = gene_ids, organism = "human", ...)
  }, error = function(e) NULL)
}
# Esta función permite reutilizar el mismo llamado para GO, KEGG y Reactome, controlando errores internamente.

# --- Ejecutar análisis ORA para cada base ---
ego_bp <- run_enrichment(entrez_ids, "GO", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.2, readable = TRUE)
ego_mf <- run_enrichment(entrez_ids, "GO", ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 0.2, readable = TRUE)
ego_cc <- run_enrichment(entrez_ids, "GO", ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 0.2, readable = TRUE)
ekegg  <- run_enrichment(entrez_ids, "KEGG", pAdjustMethod = "BH")
ereact <- run_enrichment(entrez_ids, "Reactome", pAdjustMethod = "BH", readable = TRUE)

# --- Guardar resultados si existen ---
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

# --- Leer genes up/down si existen ---
up_path <- file.path(base_dir, "DEG_upregulated.csv")       # Ruta a genes sobreexpresados
down_path <- file.path(base_dir, "DEG_downregulated.csv")   # Ruta a genes subexpresados

# Leer archivos si existen, en otro caso asignar NULL
deg_up <- if (file.exists(up_path)) read_csv(up_path, show_col_types = FALSE) else NULL
deg_down <- if (file.exists(down_path)) read_csv(down_path, show_col_types = FALSE) else NULL

# --- Enriquecimiento separado para genes UP (sobreexpresados) ---
if (!is.null(deg_up)) {
  entrez_up <- unique(na.omit(deg_up$entrez_id))  # IDs válidos
  cat("🔺 Enriquecimiento para genes sobreexpresados:", length(entrez_up), "\n")
  
  # Ejecuta análisis ORA para cada base de conocimiento
  up_go_bp <- run_enrichment(entrez_up, "GO", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.2, readable = TRUE)
  up_kegg  <- run_enrichment(entrez_up, "KEGG", pAdjustMethod = "BH")
  up_react <- run_enrichment(entrez_up, "Reactome", pAdjustMethod = "BH", readable = TRUE)
  
  # Guardar resultados
  save_and_export(up_go_bp, "UP_GO_BP")
  save_and_export(up_kegg,  "UP_KEGG")
  save_and_export(up_react, "UP_Reactome")
  
  # Visualizaciones automáticas (dotplot, cnetplot, emapplot)
  plot_all(up_go_bp, "UP_GO_BP")
  plot_all(up_kegg,  "UP_KEGG")
  plot_all(up_react, "UP_Reactome")
}

# --- Enriquecimiento separado para genes DOWN (subexpresados) ---
if (!is.null(deg_down)) {
  entrez_down <- unique(na.omit(deg_down$entrez_id))  # IDs válidos
  cat("🔻 Enriquecimiento para genes subexpresados:", length(entrez_down), "\n")
  
  # Ejecuta análisis ORA
  down_go_bp <- run_enrichment(entrez_down, "GO", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.2, readable = TRUE)
  down_kegg  <- run_enrichment(entrez_down, "KEGG", pAdjustMethod = "BH")
  down_react <- run_enrichment(entrez_down, "Reactome", pAdjustMethod = "BH", readable = TRUE)
  
  # Guardar resultados
  save_and_export(down_go_bp, "DOWN_GO_BP")
  save_and_export(down_kegg,  "DOWN_KEGG")
  save_and_export(down_react, "DOWN_Reactome")
  
  # Visualizaciones
  plot_all(down_go_bp, "DOWN_GO_BP")
  plot_all(down_kegg,  "DOWN_KEGG")
  plot_all(down_react, "DOWN_Reactome")
}

# --- Visualización robusta con protección de errores ---
plot_all <- function(obj, name) {
  if (!is.null(obj) && nrow(obj) > 0) {
    n_terms <- min(15, nrow(obj))  # Límite visual de categorías enriquecidas
    
    # Dotplot simple
    try({
      png(file.path(out_dir, paste0(name, "_dotplot.png")), width = 1000)
      print(dotplot(obj, showCategory = n_terms) + ggtitle(paste(name, " pathways")))
      dev.off()
    }, silent = TRUE)
    
    # cnetplot y emapplot solo si hay suficientes términos
    if (n_terms >= 3) {
      try({
        png(file.path(out_dir, paste0(name, "_cnetplot.png")), width = 1200)
        print(cnetplot(obj, showCategory = min(5, n_terms)))
        dev.off()
        
        png(file.path(out_dir, paste0(name, "_emapplot.png")), width = 1200)
        print(emapplot(pairwise_termsim(obj)))  # Similaridad semántica entre términos
        dev.off()
      }, silent = TRUE)
    }
  } else {
    # Aviso si no hay datos suficientes para graficar
    cat(paste0("⚠️  No se puede graficar ", name, ": sin términos enriquecidos.\n"))
  }
}

# Llamar a visualizaciones para análisis globales
plot_all(ego_bp, "GO_BP")
plot_all(ego_mf, "GO_MF")
plot_all(ego_cc, "GO_CC")
plot_all(ekegg, "KEGG")
plot_all(ereact, "Reactome")

# --- GSEA si existe columna con estadístico t ---
if ("t" %in% colnames(deg) && !any(is.na(deg$t))) {
  # Prepara el ranking de genes para GSEA (ordenado por estadístico t)
  gene_ranks <- deg %>%
    filter(!is.na(entrez_id)) %>%
    arrange(desc(t)) %>%
    distinct(entrez_id, .keep_all = TRUE) %>%
    select(entrez_id, t) %>%
    deframe()  # Convierte a vector con nombres = ENTREZID
  
  # Ejecuta GSEA sobre GO, KEGG y Reactome
  gsea_go <- gseGO(geneList = gene_ranks, OrgDb = org.Hs.eg.db, ont = "BP", verbose = FALSE)
  gsea_kegg <- gseKEGG(geneList = gene_ranks, organism = "hsa", verbose = FALSE)
  gsea_react <- gsePathway(geneList = gene_ranks, organism = "human", verbose = FALSE)
  
  # Guardar resultados
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
# ↪ Diccionario con genes relevantes para mecanismos clave del MDD

deg_sig_genes <- unique(deg_sig$symbol)
interpretacion <- c("🧠 Interpretación fisiopatológica basada en genes significativos:")

# Compara genes significativos con los mecanismos definidos
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

library(AnnotationDbi)

# Glosario para detectar mecanismos MDD por palabras clave en rutas GO/KEGG
glosario_mdd <- list(
  "Neurotransmisión" = c("synapse", "neurotransmitter", "gaba", "glutamate", "serotonin", "dopamine"),
  "Neuroinflamación" = c("cytokine", "inflammation", "immune"),
  "Eje HPA" = c("stress", "hormone", "glucocorticoid"),
  "Neuroplasticidad" = c("plasticity", "neurogenesis", "apoptosis"),
  "Disfunción mitocondrial" = c("mitochondria", "oxidative", "metabolism", "energy")
)

interpretaciones_gen <- list()

# Para cada gen significativo
for (gen in unique(na.omit(deg_sig$symbol))) {
  rutas_go <- tryCatch(AnnotationDbi::select(org.Hs.eg.db, keys = gen, keytype = "SYMBOL", columns = "GO"), error = function(e) NULL)
  rutas_kegg <- tryCatch(bitr_kegg(gen, fromType = "symbol", toType = "Path", organism = "hsa"), error = function(e) NULL)
  
  todas_rutas <- unique(c(
    if (!is.null(rutas_go)) as.character(rutas_go$GO), 
    if (!is.null(rutas_kegg)) as.character(rutas_kegg$Path)
  ))
  
  # Resumen del gen
  texto <- paste0("🔹 ", gen, " ➤ ", length(todas_rutas), " rutas encontradas")
  relacionadas <- character()
  
  # Buscar coincidencias entre rutas y palabras clave
  for (palabra in unlist(glosario_mdd)) {
    if (any(grepl(palabra, tolower(todas_rutas), fixed = TRUE))) {
      mecanismo <- names(glosario_mdd)[sapply(glosario_mdd, function(pat) palabra %in% pat)]
      relacionadas <- unique(c(relacionadas, mecanismo))
    }
  }
  
  # Añadir interpretación
  if (length(relacionadas) > 0) {
    texto <- paste(texto, "\n     → Relacionado con:", paste(unique(relacionadas), collapse = ", "))
  } else {
    texto <- paste(texto, "\n     → No se identificó relación directa con mecanismos MDD")
  }
  
  # Registrar resultado
  interpretaciones_gen[[gen]] <- list(
    gene = gen,
    rutas_detectadas = paste(todas_rutas, collapse = "; "),
    mecanismos_MDD = ifelse(length(relacionadas) > 0, paste(unique(relacionadas), collapse = ", "), "Ninguno"),
    comentario = ifelse(length(relacionadas) > 0, "Posible implicación fisiopatológica", "No asociado directamente")
  )
  
  cat(texto, "\n")
}

# Guardar resultados por gen en .csv
ruta_csv <- file.path(out_dir, "ruta_por_gen.csv")
write.csv(do.call(rbind, lapply(interpretaciones_gen, as.data.frame)), ruta_csv, row.names = FALSE)

# Guardar resumen por texto
ruta_txt <- file.path(out_dir, "interpretacion_por_gen.txt")
writeLines(unlist(lapply(interpretaciones_gen, function(x) {
  paste0("🔸 ", x$gene, ": ", x$mecanismos_MDD, " → ", x$comentario)
})), con = ruta_txt)


# --- Función para interpretar fisiopatología según lista de genes ---
interpretar_fisiopato <- function(gene_list, etiqueta) {
  deg_sig_genes <- unique(na.omit(gene_list))  # Elimina NA y duplicados
  resumen <- c(paste0("🧠 Interpretación fisiopatológica para genes ", etiqueta, ":"))
  
  for (mecanismo in names(fisiopato_genes)) {
    genes_detectados <- intersect(deg_sig_genes, fisiopato_genes[[mecanismo]])
    if (length(genes_detectados) > 0) {
      resumen <- c(resumen, paste0("✔ ", mecanismo, ": ", paste(genes_detectados, collapse = ", ")))
    }
  }
  
  # Si no se encontró ningún gen asociado
  if (length(resumen) == 1) {
    resumen <- c(resumen, "❌ No se detectaron genes relevantes en vías fisiopatológicas clave.")
  }
  return(resumen)
}

# --- Aplicar función si existen archivos DEG UP/DOWN ---
if (!is.null(deg_up)) {
  up_summary <- interpretar_fisiopato(deg_up$symbol, "sobreexpresados (UP)")
  writeLines(up_summary, file.path(out_dir, "interpretacion_fisiopatologica_UP.txt"))
}
if (!is.null(deg_down)) {
  down_summary <- interpretar_fisiopato(deg_down$symbol, "subexpresados (DOWN)")
  writeLines(down_summary, file.path(out_dir, "interpretacion_fisiopatologica_DOWN.txt"))
}

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
