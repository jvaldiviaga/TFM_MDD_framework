#!/usr/bin/env Rscript
# ==========================================================
# Comparación de logFC entre dos datasets (con/sin covariables)
# ==========================================================

# ============================================================

# FORTALEZAS Y JUSTIFICACIÓN DEL FLUJO

# Este script permite evaluar de forma cuantitativa y visual la robustez de los genes
# diferencialmente expresados (DEGs) al comparar la dirección del efecto (logFC) entre dos análisis,
# típicamente con y sin ajuste por covariables. Actúa como un módulo de validación cruzada
# que ayuda a priorizar resultados confiables.

# - Carga automáticamente los resultados DEG mapeados de ambos análisis.
# - Filtra y cruza los genes por `entrez_id`, asegurando trazabilidad biológica.
# - Solo compara genes significativos en el análisis base (adj.P.Val < 0.05), centrando la validación.
# - Calcula si los signos de los logFC coinciden, identificando genes robustos vs inestables.
# - Muestra los 10 genes más concordantes y discordantes, priorizados por magnitud.
# - Exporta la tabla comparativa con anotaciones de coincidencia y diferencia absoluta.
# - Genera un gráfico de dispersión logFC_base vs logFC_val coloreado por concordancia, lo cual
#   facilita la interpretación visual y la identificación de cambios inducidos por covariables.
# - Se integra de forma transparente con el resto del framework (estructura de carpetas, naming, etc.)

# ------------------------------------------------------------
# LIMITACIONES DETECTADAS
# ------------------------------------------------------------

# 1. El análisis es unidireccional: solo compara los genes significativos del análisis base.
#    -No detecta genes significativos en el análisis de validación si no lo eran en el base.

# 2. No hay validación del número mínimo de genes cruzados.
#    -Si n < 5, la concordancia puede ser poco representativa o engañosa.
#    -Debería emitir advertencia o evitar continuar si n es insuficiente.

# 3. Solo evalúa si el signo coincide, sin considerar magnitud o significancia en la validación.
#    -Mejora: calcular también la correlación de logFC (Pearson/Spearman).

# 4. El script no evalúa si los genes validados siguen siendo significativos en el segundo análisis,
#    ya que su objetivo es validar la dirección del efecto, no la significancia estadística.
#    Esta es una decisión metodológica deliberada centrada en robustez biológica más que en umbrales arbitrarios.


# 5. No permite establecer umbrales de cambio para marcar genes inestables.
#    -Ejemplo: flag si abs(logFC_base - logFC_val) > 1.

# 6. El gráfico es estático (PNG) y no permite identificar genes por hover o selección.
#    -Mejora futura: versión interactiva (plotly, Shiny).

# ------------------------------------------------------------
# 🚀 MEJORAS FUTURAS
# ------------------------------------------------------------

# - Añadir opción de usar genes significativos en ambos análisis, no solo en el base,
#   o permitir validar cualquier subconjunto predefinido (e.g., genes enriquecidos).
#
# - Calcular la correlación entre logFC_base y logFC_val como métrica continua
#   de estabilidad (Pearson, Spearman), además de la coincidencia de signo.
#
# - Reportar el número total de genes cruzados, el porcentaje con concordancia
#   de dirección y la distribución por magnitud de cambio.
#
# - Añadir flags de cambio significativo en dirección y magnitud,
#   por ejemplo: marcar como inestable si abs(logFC_base - logFC_val) > 1.
#
# - Exportar un ranking de genes que cambian de dirección con mayor impacto,
#   útil para priorización en análisis funcional o visualización.
#
# - Añadir tabla de contingencia que desglose los cambios UP→UP, UP→DOWN, etc.,
#   para entender mejor los efectos del ajuste por covariables.
#
# - Incorporar anotación funcional (vías GO/KEGG) para destacar si los genes
#   inestables están implicados en mecanismos clave o son biológicamente relevantes.
#
# - Permitir comparar múltiples modelos y visualizar resultados en red o matriz
#   de concordancia entre configuraciones (no solo 1 vs 1).
#
# - Generar gráficos interactivos (e.g., plotly) que permitan identificar genes
#   al hacer hover y explorar patrones de cambio más finos.
#
# - Incluir la evaluación de significancia en el segundo análisis,
#   aunque esta no es prioritaria si el objetivo principal es validar robustez
#   en la dirección del efecto.

# ============================================================


# Cargar solo los paquetes necesarios para manipulación y visualización
suppressPackageStartupMessages({
  library(dplyr)     # Manipulación de data frames
  library(readr)     # Lectura y escritura de CSV
  library(ggplot2)   # Visualización (scatter plot)
})

# Captura de argumentos desde consola: dataset base y validación
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Uso: Rscript validate_logfc_direction.R <dataset_base> <subfolder_base> <dataset_val> <subfolder_val>")
}

# Extraer los valores pasados
base_id <- args[1]
base_sub <- args[2]
val_id  <- args[3]
val_sub <- args[4]

# Función para cargar los resultados DEG mapeados de un análisis
load_deg <- function(ds_id, subfolder) {
  path <- file.path("~/TFM_MDD/results", ds_id, "differential_expression", subfolder, "DEG_results_mapped.csv")
  if (!file.exists(path)) stop("❌ No se encontró:", path)
  
  df <- read_csv(path, show_col_types = FALSE)
  
  # Filtrar genes con logFC y ENTREZID válidos
  df <- df %>% filter(!is.na(logFC), !is.na(entrez_id)) %>%
    distinct(entrez_id, .keep_all = TRUE)
  
  return(df)
}

# Cargar DEG para el análisis base (ej. sin covariables) y validación (ej. con covariables)
deg_base <- load_deg(base_id, base_sub)
deg_val  <- load_deg(val_id, val_sub)

# Filtrar solo genes significativos en el análisis base
deg_base_sig <- deg_base %>% filter(adj.P.Val < 0.05)

# Cruzar genes por ENTREZID con el análisis de validación
merged <- inner_join(
  deg_base_sig %>% select(entrez_id, GeneSymbol = symbol, logFC_base = logFC),
  deg_val %>% select(entrez_id, logFC_val = logFC),
  by = "entrez_id"
)

# Si no hay genes en común, abortar
if (nrow(merged) == 0) stop("❌ No hay genes en común para comparar.")

# Evaluar si la dirección del logFC coincide
merged <- merged %>%
  mutate(
    Same_Direction = ifelse(sign(logFC_base) == sign(logFC_val), "✅", "❌"),
    abs_diff = abs(logFC_base - logFC_val)  # Diferencia absoluta entre logFC
  )


# Resumen e impresión
n_total <- nrow(merged)
n_match <- sum(merged$Same_Direction == "✅", na.rm = TRUE)

# Imprimir resumen general
cat("🧬 Genes significativos validados:", n_total, "\n")
cat("🔁 Concordancia de dirección:", n_match, "(", round(100 * n_match / n_total, 1), "% )\n\n")

# Mostrar 10 genes con dirección consistente y mayor logFC
c("🔝 Top 10 genes coincidentes en dirección:\n")
print(merged %>% filter(Same_Direction == "✅") %>%
        arrange(desc(abs(logFC_base))) %>%
        select(GeneSymbol, logFC_base, logFC_val) %>% head(10))

# Mostrar 10 genes discordantes con mayor diferencia
cat("\n⚠️ Top 10 genes con dirección discordante:\n")
print(merged %>% filter(Same_Direction == "❌") %>%
        arrange(desc(abs_diff)) %>%
        select(GeneSymbol, logFC_base, logFC_val) %>% head(10))

# Guardar resultados
# Crear carpeta de salida (dataset validado / validado desde base)
out_dir <- file.path("~/TFM_MDD/results", val_id, "validation_from", base_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Guardar tabla comparada con entrez_id, logFCs y concordancia
write_csv(merged, file.path(out_dir, paste0("validation_logfc_", base_sub, "_vs_", val_sub, ".csv")))

# Gráfico de dispersión coloreando por coincidencia o discrepancia de dirección
png(file.path(out_dir, paste0("logFC_comparison_", base_sub, "_vs_", val_sub, ".png")), width = 1000)
ggplot(merged, aes(logFC_base, logFC_val, color = Same_Direction)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +  # Ejes de referencia
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = paste("Comparación logFC:", base_id, "vs", val_id),
       x = paste("logFC -", base_sub), y = paste("logFC -", val_sub)) +
  scale_color_manual(values = c("❌" = "red", "✅" = "blue")) +
  theme(legend.position = "bottom")
dev.off()

