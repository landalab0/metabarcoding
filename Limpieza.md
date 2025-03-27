El primer paso será limpiar nuestras secuencias de adaptadores, primers y otras secuencias
que no necesitamos. Para esto utilizaremos **Cutadapt** 

Algunas de sus funciones son:  

Eliminar adaptadores de secuencias de lectura. 
Recortar secuencias de baja calidad en los extremos. 
Filtrar lecturas por longitud mínima/máxima. 
Eliminar secuencias contaminantes no deseadas. 

Para esto crearemos dos directorios, uno de entrada con los datos crudos FASTQ y otro de salida de los archivos recortados. 

```bash
mkdir –p  data  # Contiene los archivos fastq originales 
mkdir –p  Out_dir  # Carpeta para los archivos recortados 
```

Ahora definiremos las secuencias de nuestros adaptadores 

```bash
FORWARD_ADAPTER="CCTACGGGNGGCWGCAG" 

REVERSE_ADAPTER="GACTACHVGGGTATCTAATCC"
```

Configuraremos los parámetros de cutadapt

```bash
#ERROR_RATE="0.1" 

#MINIMUM_LENGTH="1" 

#OVERLAP="3" 

#QUALITY_CUTOFF="0" 

#QUALITY_BASE="33" 

#CORES="1" 
```

