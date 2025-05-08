#!/bin/bash
# ==========================================================
# Script: prepare_geo_dataset.sh (versión generalizada)
# Objetivo: Automatizar la preparación de datasets GEO para cualquier tipo de ómica (.CEL, .counts, .fastq, .idat, etc.)
# Uso: (base) jessica@DESKTOP-UFC7TIA:~/TFM_MDD$ bash scripts/general/prepare_geo_dataset.sh GSE[44593[ID]
# ==========================================================

# ==========================================================
# FORTALEZAS Y JUSTIFICACIÓN DEL FLUJO:
# - Automatiza la preparación completa de cualquier dataset descargado desde GEO, sin importar su tipo.
# - Detecta automáticamente el tipo de archivo (.CEL, .counts, .fastq, .idat, .vcf, .nii.gz...) y lo clasifica.
# - Genera el archivo `dataset_type.txt` que guía al preprocesamiento según la ómica detectada.
# - Implementa detección progresiva: descarga primero vía curl y, si falla, usa R (GEOquery).
# - Organiza los datos en carpetas tipo `final/` para evitar mezclar archivos temporales.
# - Compatible con todo tipo de datos ómicos en un framework multiómico unificado.

# ==========================================================

# 1. Captura del ID del dataset
# Obtenemos el primer argumento pasado por línea de comandos (ej. GSE98793)
GSE_ID=$1

# Si no se pasó argumento, mostramos error y salimos
if [ -z "$GSE_ID" ]; then
  echo "❌ Debes proporcionar un ID GEO como GSE98793"
  exit 1
fi

# 2. Definición de rutas base
# Carpeta base donde se almacenarán los archivos
BASE_DIR=~/TFM_MDD/data/transcriptomics/$GSE_ID

# Ruta esperada del archivo .tar que contiene los datos crudos
RAW_TAR="$BASE_DIR/${GSE_ID}_RAW.tar"

# Carpeta final donde moveremos los archivos útiles (descomprimidos y clasificados)
FINAL_DIR="$BASE_DIR/final"

# Creamos la carpeta final si no existe
mkdir -p "$FINAL_DIR"

# LIMITACIÓN ACTUAL:
# Todos los archivos, sin importar el tipo (RNA-seq, idat, .vcf, imagen...), se almacenan dentro de data/transcriptomics/<GSE_ID>
# MEJORA FUTURA:
# Tras detectar el tipo de dato, mover los archivos a su subdirectorio adecuado:
#   - data/transcriptomics/         → para .CEL, .counts
#   - data/epigenomics/             → para .idat
#   - data/genomics/                → para .vcf, .bed
#   - data/neuroimaging/            → para .nii.gz, .fif
# Esto permitirá una estructura multiómica real y facilitará la automatización del pipeline.


# 3. Intento de descarga vía FTP
# URL directa al archivo .tar usando estructura estándar de GEO FTP
FTP_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/${GSE_ID:0:6}nnn/$GSE_ID/suppl/${GSE_ID}_RAW.tar"

# Intentamos descargar con curl, con reintento en caso de corte (-C -)
echo "🌐 Intentando descarga desde FTP con curl..."
curl -C - -L "$FTP_URL" -o "$RAW_TAR"

# 4. Verificación de archivo descargado
# Verificamos si el archivo existe y pesa más de 50 MB (umbral de validez)
if [ -f "$RAW_TAR" ] && [ $(stat -c%s "$RAW_TAR") -gt 50000000 ]; then
  echo "✅ Descarga completada correctamente."
  cd "$BASE_DIR"
else
  # Si falla, recurrimos al script en R que usa GEOquery como alternativa
  echo "⚠️ Descarga incompleta o fallida. Intentando con GEOquery..."
  Rscript ~/TFM_MDD/scripts/general/download_raw_geo.R "$GSE_ID"
  
  # Ajustamos ruta por si el R script creó una subcarpeta adicional
  RAW_TAR="$BASE_DIR/$GSE_ID/${GSE_ID}_RAW.tar"
  cd "$BASE_DIR/$GSE_ID" 2>/dev/null || cd "$BASE_DIR" || exit
fi

# 5. Validación final
# Si el archivo aún no existe, abortamos el script
if [ ! -f "$RAW_TAR" ]; then
  echo "❌ No se pudo obtener un archivo .tar válido."
  exit 1
fi

# 6. Extracción del .tar
echo "📦 Extrayendo contenido del .tar..."
tar -xvf "$RAW_TAR"

# 7. Detección automática de tipo de archivo (Generalizado)
# Este bloque detecta qué tipo de datos contiene el .tar y los clasifica automáticamente
# Además, guarda esta clasificación en un archivo de texto para que el pipeline la reconozca

# Cada bloque detecta un tipo distinto de archivo y actúa en consecuencia

FILE_LIST=$(ls)

if ls *.CEL.gz 1> /dev/null 2>&1; then
  echo "🧬 Detectados archivos .CEL.gz (Microarrays)."
  echo "Tipo de datos: Microarrays (.CEL)" > "$BASE_DIR/dataset_type.txt"
  
  echo "🧬 Descomprimiendo archivos .CEL.gz..."
  gunzip -f *.CEL.gz

  echo "📂 Moviendo archivos .CEL al directorio final..."
  mv *.CEL "$FINAL_DIR"

elif ls *.counts.txt 1> /dev/null 2>&1 || ls *.csv 1> /dev/null 2>&1; then
  echo "🧬 Detectados archivos de conteo (RNA-seq counts)."
  echo "Tipo de datos: RNA-seq (matrices de conteo)" > "$BASE_DIR/dataset_type.txt"

  echo "📂 Moviendo archivos de conteo al directorio final..."
  mv *.counts.txt "$FINAL_DIR" 2>/dev/null
  mv *.csv "$FINAL_DIR" 2>/dev/null

elif ls *.txt.gz 1> /dev/null 2>&1; then
  echo "🧬 Detectados archivos .txt.gz (RNA-seq por muestra, comprimidos)."
  echo "Tipo de datos: RNA-seq (por muestra .txt.gz)" > "$BASE_DIR/dataset_type.txt"

  echo "🧩 Descomprimiendo archivos .txt.gz..."
  gunzip -f *.txt.gz

  echo "📂 Moviendo archivos .txt descomprimidos al directorio final..."
  mv *.txt "$FINAL_DIR"

elif ls *.fastq.gz 1> /dev/null 2>&1; then
  echo "🧬 Detectados archivos .fastq.gz (RNA-seq crudo)."
  echo "Tipo de datos: RNA-seq crudo (.fastq)" > "$BASE_DIR/dataset_type.txt"

  echo "📂 Moviendo archivos .fastq.gz al directorio final..."
  mv *.fastq.gz "$FINAL_DIR"

elif ls *.idat 1> /dev/null 2>&1; then
  echo "🧬 Detectados archivos .idat (Metilación)."
  echo "Tipo de datos: Epigenómica (.idat)" > "$BASE_DIR/dataset_type.txt"

  echo "📂 Moviendo archivos .idat al directorio final..."
  mv *.idat "$FINAL_DIR"

elif ls *.vcf 1> /dev/null 2>&1 || ls *.bed 1> /dev/null 2>&1; then
  echo "🧬 Detectados archivos .vcf o .bed (Genómica/GWAS)."
  echo "Tipo de datos: Genómica (GWAS)" > "$BASE_DIR/dataset_type.txt"

  echo "📂 Moviendo archivos .vcf/.bed al directorio final..."
  mv *.vcf "$FINAL_DIR" 2>/dev/null
  mv *.bed "$FINAL_DIR" 2>/dev/null

elif ls *.fif 1> /dev/null 2>&1 || ls *.nii.gz 1> /dev/null 2>&1; then
  echo "🧬 Detectados archivos .fif o .nii.gz (Neuroimagen/MEG/fMRI)."
  echo "Tipo de datos: Imagen funcional (MEG/fMRI)" > "$BASE_DIR/dataset_type.txt"

  echo "📂 Moviendo archivos de imagen al directorio final..."
  mv *.fif "$FINAL_DIR" 2>/dev/null
  mv *.nii.gz "$FINAL_DIR" 2>/dev/null

else
  echo "⚠️ Tipo de datos no detectado automáticamente."
  echo "Tipo de datos: Desconocido" > "$BASE_DIR/dataset_type.txt"
fi

# 8. Descomprimir automáticamente todos los archivos .gz restantes
echo "🧩 Descomprimiendo automáticamente cualquier archivo .gz restante..."
find "$FINAL_DIR" -name "*.gz" -exec gunzip -f {} \;
find "$BASE_DIR" -name "*.gz" -exec gunzip -f {} \;

# 9. Mensaje de cierre
echo "✅ Dataset $GSE_ID preparado correctamente en: $FINAL_DIR"
echo "🗒️ Información de tipo de datos registrada en: dataset_type.txt"

