#!/usr/bin/env Rscript

suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(tidyverse)
  library(ggplot2)
  library(rmarkdown)
  library(igraph)
  library(ggraph)
})

cat("📥 Cargando datos reales...\n")

# Archivos de entrada
genes_file <- "~/TFM_MDD/results/transcriptomics/VALIDATION_INTER_ENRICHMENT/inter_dataset_robustness_summary.csv"
lasso_file <- "~/TFM_MDD/results/ds005356/stats/lasso_selected_features.csv"
regression_file <- "~/TFM_MDD/results/ds005356/stats/regression_results_linear.csv"
output_dir <- "~/TFM_MDD/results/interpretation_multiomics/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ================= CARGA GENES ===================
genes <- read_csv(genes_file, show_col_types = FALSE) %>% filter(total_score >= 3)
gene_ids <- bitr(genes$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene_entrez <- gene_ids$ENTREZID

# ================= ENRIQUECIMIENTO FUNCIONAL ===================
ego <- enrichGO(gene = gene_entrez, OrgDb = org.Hs.eg.db, ont="BP", pAdjustMethod="fdr", readable=TRUE)

mecanismos_dict <- c(
  "plasticity|neurogenesis|synaptic" = "Neuroplasticidad",
  "serotonin|dopamin|GABA|glycine" = "Neurotransmisión",
  "inflamm|immune|cytokine" = "Inflamación / inmunidad",
  "oxidative|mitochondr" = "Estrés oxidativo / mitocondrial",
  "signaling|pathway|transduction" = "Señalización intracelular",
  "circadian|rhythm" = "Disfunción circadiana"
)

classify_mechanism <- function(desc) {
  for (key in names(mecanismos_dict)) {
    if (grepl(key, desc, ignore.case = TRUE)) return(mecanismos_dict[key])
  }
  return(NA)
}

ego@result$Mecanismo <- sapply(ego@result$Description, classify_mechanism)

ruta_por_gen <- ego@result %>%
  filter(!is.na(Mecanismo)) %>%
  separate_rows(geneID, sep="/") %>%
  group_by(geneID) %>%
  summarise(
    Rutas = paste(unique(Description), collapse="; "),
    Mecanismos = paste(unique(Mecanismo), collapse="; ")
  ) %>%
  left_join(genes, by = c("geneID" = "SYMBOL"))

# ================== MEG AGRUPADO ====================
lasso <- read_csv(lasso_file, show_col_types = FALSE)
reg <- read_csv(regression_file, show_col_types = FALSE) %>% filter(P_value < 0.05)
sig_meg_raw <- semi_join(reg, lasso, by = "Feature") %>% filter(str_detect(Feature, "rwp|p300|gfp|slope"))

meg_label_dict <- tibble::tribble(
  ~codigo, ~descripcion,
  "rwp", "Reward Positivity (RWP): indicador temprano del procesamiento de recompensa",
  "p300", "P300: potencial tardío vinculado a la atención y control cognitivo",
  "gfp", "Global Field Power (GFP): medida general de actividad cerebral sincronizada",
  "slope", "Pendiente temporal de RWP: velocidad del procesamiento de recompensa"
)

condiciones_var_dict <- tibble::tribble(
  ~cond, ~descripcion,
  "cond_05", "fase de anticipación de recompensa",
  "cond_07", "respuesta a recompensa recibida",
  "cond_03", "fase temprana de expectativa (cue onset)",
  "cond_04", "predicción de feedback positivo"
)

sig_meg <- sig_meg_raw %>%
  mutate(
    tipo = str_extract(Feature, "rwp|p300|gfp|slope"),
    cond = str_extract(Feature, "cond_\\d+")
  ) %>%
  left_join(meg_label_dict, by = c("tipo" = "codigo")) %>%
  left_join(condiciones_var_dict, by = "cond")

sig_meg_grouped <- sig_meg %>%
  group_by(tipo, descripcion.x) %>%
  summarise(
    Condiciones = paste(unique(descripcion.y), collapse=", "),
    Min_p = min(P_value),
    Ejemplos = paste(unique(Feature), collapse=", ")
  )

# ================== TEXTO INTERPRETATIVO AGRUPADO ====================
texto_interpretativo <- c("# 🧬 Informe Clínico Interpretativo de Biomarcadores del MDD\n")

for(i in 1:nrow(ruta_por_gen)) {
  texto_interpretativo <- c(texto_interpretativo,
                            paste0("🧬 El gen ", ruta_por_gen$geneID[i], " participa en las rutas: ", ruta_por_gen$Rutas[i], 
                                   ". Estas están implicadas en ", ruta_por_gen$Mecanismos[i], 
                                   ", mecanismos conocidos por estar alterados en MDD."))
}

texto_interpretativo <- c(texto_interpretativo, "\n# 🧠 Biomarcadores MEG significativos\n")

for(i in 1:nrow(sig_meg_grouped)) {
  texto_interpretativo <- c(texto_interpretativo,
                            paste0("🧠 Variable ", sig_meg_grouped$tipo[i], " (p≈", round(sig_meg_grouped$Min_p[i],3), "): ", 
                                   sig_meg_grouped$descripcion.x[i], ". Observada en condiciones: ", sig_meg_grouped$Condiciones[i], "."))
}

writeLines(texto_interpretativo, file.path(output_dir, "biomarker_summary_report.txt"))

# PDF y RMarkdown
template_rmd <- file.path(output_dir, "biomarker_summary_report.Rmd")
writeLines(c(
  "---",
  "title: 'Informe Biomarcadores MDD'",
  "output:",
  "  pdf_document:",
  "    latex_engine: xelatex",
  "---",
  "",
  texto_interpretativo
), template_rmd)
rmarkdown::render(input = template_rmd, output_file = file.path(output_dir, "biomarker_summary_report.pdf"))

# ================== VISUALIZACIONES ====================
cat("📊 Generando visualizaciones...\n")

# 1. genes_significativos
png(file.path(output_dir, "genes_significativos.png"), width=1000, height=800)
ggplot(genes, aes(x=reorder(SYMBOL, total_score), y=total_score)) +
  geom_col(fill="gray30") + coord_flip() +
  labs(title="Genes significativos en MDD", x="Gen", y="Score") +
  theme_minimal(base_size=14)
dev.off()

# 2. features_meg_significativas
png(file.path(output_dir, "features_meg_significativas.png"), width=1000, height=800)
ggplot(sig_meg_raw, aes(x=reorder(Feature, -log10(P_value)), y=-log10(P_value))) +
  geom_col(fill="blue") + coord_flip() +
  labs(title="Features MEG significativas en MDD", x="Feature", y="-log10(p-value)") +
  theme_minimal(base_size=14)
dev.off()

# 3. dotplot_enrichment_go
png(file.path(output_dir, "dotplot_enrichment_go.png"), width=1200, height=800)
plot(dotplot(ego, showCategory=10, title="GO Enrichment - Genes MDD"))
dev.off()

# 4. biomarkers_by_mechanism_barplot (único válido)
gene_mech <- ego@result %>%
  filter(!is.na(Mecanismo)) %>%
  separate_rows(geneID, sep="/") %>%
  distinct(geneID, Mecanismo) %>%
  mutate(tipo = "Transcriptómica") %>%
  rename(biomarcador = geneID)

meg_mech <- sig_meg_raw %>%
  mutate(Mecanismo = "Procesamiento atencional / cognitivo", tipo = "MEG") %>%
  rename(biomarcador = Feature)

graph_mech <- bind_rows(gene_mech, meg_mech)

png(file.path(output_dir, "biomarkers_by_mechanism_barplot.png"), width=1200, height=800)
ggplot(graph_mech, aes(x=Mecanismo, fill=tipo)) +
  geom_bar() +
  labs(title="Biomarcadores por mecanismo fisiopatológico y capa ómica", x="Mecanismo fisiopatológico", y="Nº de biomarcadores") +
  theme_minimal(base_size=14) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()

# 5. arcplot_genes_meg_mecanisms
genes_nodos <- ruta_por_gen$geneID
mecanismos_nodos <- unique(unlist(str_split(ruta_por_gen$Mecanismos, "; ")))
meg_nodos <- unique(sig_meg_grouped$tipo)
nodo_funcional <- "Procesamiento atencional / cognitivo"

nodes <- tibble(name = unique(c(genes_nodos, mecanismos_nodos, meg_nodos, nodo_funcional)))

edges1 <- ruta_por_gen %>%
  separate_rows(Mecanismos, sep="; ") %>%
  transmute(from = geneID, to = Mecanismos)

edges2 <- sig_meg_grouped %>%
  transmute(from = tipo, to = nodo_funcional)

edges <- bind_rows(edges1, edges2) %>% filter(from != to)

# Validación cruzada: asegura que todos los nodos usados están en 'nodes'
edges <- edges %>% filter(from %in% nodes$name & to %in% nodes$name)

g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

arc_file <- file.path(output_dir, "arcplot_genes_meg_mecanisms.png")
png(arc_file, width=1600, height=1000)
ggraph(g, layout = 'linear', circular = TRUE) +
  geom_edge_arc(alpha=0.4, color="grey50") +
  geom_node_point(size=3) +
  geom_node_text(aes(label=name), size=3, repel=TRUE) +
  theme_void()
dev.off()

cat("✅ Informe generado exitosamente en: ", output_dir, "\n")
