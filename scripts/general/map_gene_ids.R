#!/usr/bin/env Rscript
# ==========================================================
# Script: map_gene_ids.R (multi-subcarpeta)
# Objetivo: Mapear a ENTREZID/SYMBOL en cada comparación por separado
# Uso: Rscript scripts/map_gene_ids.R GSE98793
# En mi caso: (base) jessica@DESKTOP-UFC7TIA:~/TFM_MDD/scripts/transcriptomics$ Rscript ~/TFM_MDD/scripts/transcriptomics/map_gene_ids.R GSE[ID]
# ==========================================================

# ============================================================

# OBJETIVO GENERAL:
# Este script mapea automáticamente los genes detectados como diferencialmente expresados
# (DEG) a identificadores estándar (ENTREZID y SYMBOL) compatibles con análisis funcionales
# posteriores (GO, KEGG, Reactome), iterando sobre cada configuración del análisis diferencial.

# FORTALEZAS Y JUSTIFICACIÓN DEL FLUJO:
# - Automatiza el proceso de mapeo en todos los análisis diferenciales realizados.
# - Usa múltiples tipos de ID (SYMBOL, ENSEMBL, PROBEID...) para maximizar el éxito del mapeo.
# - Aplica un fallback automático con anotación GPL570 si falla `org.Hs.eg.db`.
# - Exporta los resultados mapeados por configuración, junto con un resumen interpretativo.
# - Se integra sin intervención manual dentro del pipeline general del framework.

# ------------------------------------------------------------
#  LIMITACIONES Y MEJORAS PENDIENTES
# ------------------------------------------------------------

# 1. Solo se mapea `DEG_results.csv`, ignorando `DEG_upregulated.csv` y `DEG_downregulated.csv`,
#    que son esenciales para análisis de enriquecimiento separados.

# 2. Se asume que la primera columna del CSV contiene los IDs de gen.
#    Esto puede no cumplirse y debe corregirse con detección automática del identificador real.

# 3. El mapeo se considera exitoso si al menos 5 genes tienen ENTREZID,
#    lo cual es un umbral arbitrario y débil, sin justificación metodológica.

# 4. El fallback a la anotación GPL570 solo es válido para datasets Affymetrix específicos.
#    No sirve para RNA-seq, otras plataformas microarray ni datos de metilación.

# ------------------------------------------------------------
# MEJORAS FUTURAS
# ------------------------------------------------------------

# - Iterar sobre los subconjuntos `DEG_upregulated.csv` y `DEG_downregulated.csv` también.
# - Detectar automáticamente la columna que contiene los identificadores de gen (o usar rownames).
# - Sustituir el umbral fijo `≥ 5 genes` por uno porcentual o configurable.
# - Detectar automáticamente la plataforma del dataset (p. ej., GPL10558, ENSEMBL, etc.).
# - Permitir especificar manualmente el tipo de identificador a mapear como argumento.
# - Añadir soporte para IDs propios de RNA-seq (ENSG...) sin asumir plataforma fija.
# - Generar un log único por dataset que resuma los mapeos totales y parciales.
# - Añadir soporte para `multiVals = list` para detectar genes con múltiples anotaciones.

# ============================================================

# --- Carga de librerías necesarias ---
suppressPackageStartupMessages({
  library(readr)           # Leer archivos CSV
  library(dplyr)           # Manipulación de dataframes
  library(AnnotationDbi)   # Framework de anotación de Bioconductor
  library(org.Hs.eg.db)    # Base de datos de anotaciones para Homo sapiens
  library(fs)              # Manejo de rutas y sistema de archivos
})

# --- Argumento ---
args <- commandArgs(trailingOnly = TRUE)
dataset_id <- ifelse(length(args) > 0, args[1], stop("❌ Indica el dataset (ej. GSE98793)"))

# --- Buscar subcarpetas dentro de differential_expression ---
base_dir <- file.path("~/TFM_MDD/results", dataset_id, "differential_expression")

# Lista las subcarpetas dentro del directorio de resultados de expresión diferencial
subdirs <- dir_ls(base_dir, type = "directory")

# Si no hay subcarpetas, se detiene el script
if (length(subdirs) == 0) stop("❌ No hay subcarpetas en:", base_dir)

# Itera sobre cada subcarpeta (diferente configuración de covariables, etc.)
for (subdir in subdirs) {
  cat("\n📁 Procesando:", path_file(subdir), "\n")
  input_path <- file.path(subdir, "DEG_results.csv")           # Archivo de entrada
  output_path <- file.path(subdir, "DEG_results_mapped.csv")   # Archivo de salida
  summary_path <- file.path(subdir, "mapping_summary.txt")     # Resumen de mapeo
  
  
  if (!file.exists(input_path)) {
    cat("⚠️  Archivo no encontrado:", input_path, "\n")
    next  # Salta a la siguiente subcarpeta
  }
  
  # Lee el archivo de DEG y renombra la primera columna como 'gene'
  deg <- read_csv(input_path, show_col_types = FALSE)
  colnames(deg)[1] <- "gene"
  gene_ids <- deg$gene
  total_genes <- length(gene_ids)
  success <- FALSE  # Flag de éxito en el mapeo
  
  # --- Mapear usando org.Hs.eg.db ---
  # Intenta varios tipos de identificadores (SYMBOL, ENSEMBL, etc.)
  for (type in c("SYMBOL", "PROBEID", "AFFY", "ENSEMBL", "ENTREZID")) {
    message("🔍 Probando tipo de ID: ", type)
    try({
      symbol_map <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL", keytype = type, multiVals = "first")
      entrez_map <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "ENTREZID", keytype = type, multiVals = "first")
      
      # Considera que el mapeo fue exitoso si al menos 5 genes tienen ENTREZID
      if (sum(!is.na(entrez_map)) >= 5) {
        deg$symbol <- symbol_map[deg$gene]
        deg$entrez_id <- entrez_map[deg$gene]
        success <- TRUE
        break
      }
    }, silent = TRUE)
  }
  
  # --- Alternativa: anotación GPL570 ---
  # Si fallan todos los métodos anteriores, intenta con la tabla de anotación de GEO
  if (!success) {
    message("⚠️ No se pudo mapear vía org.Hs.eg.db. Intentando con anotación GPL570...")
    gpl_url <- "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"
    gpl_file <- tempfile(fileext = ".gz")
    download.file(gpl_url, gpl_file)
    
    lines <- readLines(gzfile(gpl_file))
    table_start <- grep("^ID\\t", lines)[1]
    if (is.na(table_start)) next
    
    tmp_txt <- tempfile(fileext = ".txt")
    writeLines(lines[table_start:length(lines)], tmp_txt)
    annot <- read.delim(tmp_txt, header = TRUE, sep = "\t", quote = "")
    
    # Detecta las columnas relevantes para ID, símbolo y ENTREZID
    cols <- tolower(colnames(annot))
    id_col     <- colnames(annot)[which(cols %in% c("id", "id_ref"))[1]]
    symbol_col <- colnames(annot)[which(cols %in% c("gene.symbol", "gene symbol", "gene title"))[1]]
    entrez_col <- colnames(annot)[which(cols %in% c("gene.id", "entrez gene", "entrezgene", "entrez_gene_id"))[1]]
    
    # Si encuentra todas las columnas necesarias, hace el join con los DEGs
    if (!any(is.na(c(id_col, symbol_col, entrez_col)))) {
      map <- annot[, c(id_col, symbol_col, entrez_col)]
      colnames(map) <- c("gene", "symbol", "entrez_id")
      deg_tmp <- left_join(deg, map, by = "gene")
      if ("entrez_id" %in% colnames(deg_tmp) && sum(!is.na(deg_tmp$entrez_id)) >= 5) {
        deg <- deg_tmp
        success <- TRUE
      }
    }
  }
  
  # --- Alternativa: anotación GPL10558 (Illumina HumanHT-12) ---
  if (!success) {
    message("⚠️ Intentando con anotación Illumina GPL10558...")
    gpl_url <- "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz"
    gpl_file <- tempfile(fileext = ".gz")
    download.file(gpl_url, gpl_file, quiet = TRUE)
    
    lines <- readLines(gzfile(gpl_file))
    table_start <- grep("^ID\\t", lines)[1]
    if (!is.na(table_start)) {
      tmp_txt <- tempfile(fileext = ".txt")
      writeLines(lines[table_start:length(lines)], tmp_txt)
      annot <- read.delim(tmp_txt, header = TRUE, sep = "\t", quote = "")
      
      # Detecta columnas relevantes
      cols <- tolower(colnames(annot))
      id_col     <- colnames(annot)[which(cols %in% c("id", "id_ref"))[1]]
      symbol_col <- colnames(annot)[which(cols %in% c("gene.symbol", "gene symbol", "gene title"))[1]]
      entrez_col <- colnames(annot)[which(cols %in% c("gene.id", "entrez gene", "entrezgene", "entrez_gene_id"))[1]]
      
      if (!any(is.na(c(id_col, symbol_col, entrez_col)))) {
        map <- annot[, c(id_col, symbol_col, entrez_col)]
        colnames(map) <- c("gene", "symbol", "entrez_id")
        deg <- left_join(deg, map, by = "gene")
        if (sum(!is.na(deg$entrez_id)) >= 5) {
          success <- TRUE
        }
      }
    }
  }
  
  
  if (success) {
    # Filtra genes que realmente tienen un ENTREZID válido
    deg <- deg %>% filter(!is.na(entrez_id) & entrez_id != "")
    mapped_genes <- nrow(deg)
    lost_genes <- total_genes - mapped_genes
    
    # Guarda el archivo mapeado
    write_csv(deg, output_path)
    cat("✅ Archivo guardado en:", output_path, "\n")
    cat("🧬 Genes totales:", total_genes, " | Mapeados:", mapped_genes, " | Perdidos:", lost_genes, "\n")
    
    # Genera un resumen interpretativo
    resumen <- c(
      paste0("🧬 Comparación: ", path_file(subdir)),
      paste0("🧪 Genes originales: ", total_genes),
      paste0("🔍 Genes mapeados a ENTREZID: ", mapped_genes),
      paste0("⚠️ Genes perdidos: ", lost_genes),
      if (mapped_genes < 50) "❗ Mapeo insuficiente (<50 genes). Verificar IDs o tecnología." else "✅ Mapeo suficiente para enriquecimiento.",
      paste0("🗓️ Fecha: ", Sys.time())
    )
    writeLines(resumen, summary_path)
    cat(paste(resumen, collapse = "\n"), "\n")
    
  } else {
    cat("❌ No se pudo mapear ningún identificador válido\n")
  }
}