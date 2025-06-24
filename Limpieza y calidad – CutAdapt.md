# Limpieza y calidad – CutAdapt🧹

Para iniciar con el flujo de trabajo, deberás tener tus archivos FASTQ (si no los has descargado, puedes encontrarlos los archivos de prueba aquí).

El primer paso será limpiar nuestras secuencias de adaptadores, primers y otras secuencias
que no necesitamos. Para esto utilizaremos **Cutadapt** 

Algunas de sus funciones son:  

1. Eliminar adaptadores de secuencias de lectura. 
2. Recortar secuencias de baja calidad en los extremos. 
3. Filtrar lecturas por longitud mínima/máxima. 
4. Eliminar secuencias contaminantes no deseadas.

Si aún no lo tienes instalado, puedes revisar este [link](https://cutadapt.readthedocs.io/en/stable/installation.html) 

Primero posicionate en tu directorio utilizando los comandos

```bash
cd  ruta/a/tu/repositorio # Ingresa la ruta de acceso donde tengas tus archivos
ls # Enlista el contenido dentro de la ruta 
```

Si tu archivo requiere ser descomprimido utiliza el siguiente comando:

```bash
tar -xvzf archivo.tar.gz
```

Después de descomprimir, la estructura del directorio será la siguiente:

```bash
# Este es un ejemplo de la estructura del directorio una vez que se haya descomprimido el archivo.
    archivo/
    ├── data/
    │   ├── 24h_R2.fastq.gz
    │   ├── 48h_R2.fastq.gz
    │   ├── dbs/
    │   │   └── classifier_silva_138_trained.qza
    │   ├── 24h_R1.fastq.gz
    │   ├── 0h_R1.fastq.gz
    │   ├── 9h_R2.fastq.gz
    │   ├── 9h_R1.fastq.gz
    │   ├── 48h_R1.fastq.gz
    │   ├── 0h_R2.fastq.gz
    ├── results/
    ├── src/
    │   └── create_manifest.sh
```

Confirma que han sido descomprimidos

```bash
cd "c/ruta/a/tu/repositorio/archivo/data"
ls # Se enlistan los fastq.gz
```

Posteriormente crearemos y guardaremos el siguiente código como un script en Bash

```bash
nano run_cutadapt.sh
```

Que contendrá los siguientes puntos

```bash
# Crearemos dos directorios, uno de entrada con nuestros datos crudos (archivos FASTQ) y otro de salida de los archivos recortados. 

mkdir –p  data  # Contiene los archivos fastq originales 
mkdir –p  01.cutadapt  # Carpeta para los archivos recortados 

# Directorio de entrada y salida

RAW_DIR="C:/Users/amplicones/data" # Directorio de entrada con archivos FASTQ
OUT_DIR="C:/Users/amplicones/01.cutadapt" #Directorio de salida con archivos recortados

# Ahora definiremos las secuencias de nuestros adaptadores 

FORWARD_ADAPTER="CCTACGGGNGGCWGCAG" 
REVERSE_ADAPTER="GACTACHVGGGTATCTAATCC"

# Configuraremos los parámetros de cutadapt

ERROR_RATE="0.1" 
MINIMUM_LENGTH="1" 
OVERLAP="3" 
QUALITY_CUTOFF="0" 
QUALITY_BASE="33" 
CORES="1" 

# Proceso para cada par R1/R2 - Ejecución de Cutadapt 

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

Guardalo con ctrl + O; luego Enter y Sal con ctrl + X

Da el permiso de ejecutar

```bash
chmod +x run_cutadapt.sh
```

Y correrlo desde la terminal de Bash

```bash
bash run_cutadapt.sh
```

Después de ejecutarlo, revisa que la carpeta 01.cutadapt/ tenga los archivos recortados

```bash
ls 01.cutadapt/
```

Dando lo siguiente:

```
0h_R1.fastq.gz  24h_R1.fastq.gz  48h_R1.fastq.gz  9h_R1.fastq.gz
0h_R2.fastq.gz  24h_R2.fastq.gz  48h_R2.fastq.gz  9h_R2.fastq.gz
```


















 
