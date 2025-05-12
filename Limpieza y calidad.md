# Limpieza y Calidad 🧹

El primer paso será limpiar nuestras secuencias de adaptadores, primers y otras secuencias
que no necesitamos. Para esto utilizaremos **Cutadapt** 

Algunas de sus funciones son:  

1. Eliminar adaptadores de secuencias de lectura. 
2. Recortar secuencias de baja calidad en los extremos. 
3. Filtrar lecturas por longitud mínima/máxima. 
4. Eliminar secuencias contaminantes no deseadas.

Si aún no lo tienes instalado, puedes revisar este [link](https://cutadapt.readthedocs.io/en/stable/installation.html) 

Para esto crearemos dos directorios, uno de entrada con nuestros datos crudos (archivos FASTQ) y otro de salida de los archivos recortados. 

```bash
mkdir –p  data  # Contiene los archivos fastq originales 
mkdir –p  01.cutadapt  # Carpeta para los archivos recortados 
```

Directorio de entrada y salida

RAW_DIR="C:/Users/amplicones/data" # Directorio de entrada con archivos FASTQ
OUT_DIR="C:/Users/amplicones/01.cutadapt" #Directorio de salida con archivos recortados

Ahora definiremos las secuencias de nuestros adaptadores 

```bash
FORWARD_ADAPTER="CCTACGGGNGGCWGCAG" 

REVERSE_ADAPTER="GACTACHVGGGTATCTAATCC"
```

Configuraremos los parámetros de cutadapt

```bash
ERROR_RATE="0.1" 

MINIMUM_LENGTH="1" 

OVERLAP="3" 

QUALITY_CUTOFF="0" 

QUALITY_BASE="33" 

CORES="1" 
```
Ejecutaremos Cutadapt 

```bash
for FILE in "$RAW_DIR"/*_R1_001.fastq.gz; do
base=$(basename "$FILE" _R1_001.fastq.gz)
R2_FILE="$RAW_DIR/${base}_R2_001.fastq.gz"

if [ -f "$R2_FILE" ]; then
echo "Procesando: $base"
cutadapt \
-g $FORWARD_ADAPTER \
-G $REVERSE_ADAPTER \
--error-rate $ERROR_RATE \
--minimum-length $MINIMUM_LENGTH \
--overlap $OVERLAP \
--pair-adapters \
-o "$OUT_DIR/${base}_R1_001.fastq.gz" \
-p "$OUT_DIR/${base}_R2_001.fastq.gz" \
--quality-cutoff $QUALITY_CUTOFF \
--quality-base $QUALITY_BASE \
-j $CORES \
"$FILE" "$R2_FILE"
else
 echo "Advertencia: No se encontró el archivo R2 para ${base}_R1.fastq.gz"
fi
done
```

Para utilizar este código podemos guardarlo como un script en Bash

```bash
cutadapt.sh
```

Y lo corremos desde la terminal de Bash

```bash
bash cutadapt.sh
```
