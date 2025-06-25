## Control de calidad con DADA2 🔍

Ahora, nos iremos a R para trabajar con dada2.

Tip: si requieres actualizarlo, ejecuta el siguiente comando

``` r
install.packages("installr")
library(installr)
updateR()
``` 

Instalaremos el paquete dada2

``` r
install.packages("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.16") 
``` 

Cargaremos las librerías necesarias

``` r
library(dada2) # pipeline principal
library(Biostrings) # manejo de secuencias
library(ShortRead) # lecturas de archivos FASTQ
```

Establece el directorio de trabajo si aún no está configurado y verificalo

```r
setwd (ruta/de/tu/directorio)
getwd()
```
Definamos la ruta de nuestros archivos recortados (los resultados del paso anterior)

```r
path.cut <- file.path("01.cutadapt/")
if(!dir.exists(path.cut)) dir.create(path.cut)
```
Listaremos los archivos R1 y R2

```r
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz$", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz$", full.names = TRUE))
```

Para obtener los nombres de nuestras muestras:

```r
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
```

Verificamos

```r
cutFs
cutRs
head(sample.names)
 ``` 

Y finalmente graficaremos la calidad.

```r
plotQualityProfile(cutFs[1:2]) 
plotQualityProfile(cutRs[1:2])
```

## Filtrado y calidad ✂️

Creamos la carpeta de salida para nuestros nuevos archivos filtrados

```r
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
```

La función que usaremos a continuación *filterAndTrim* tiene distintos parámetros que podrás cambiar a tu conveniencia (ten en cuenta la longitud y calidad de tus lecturas), algunos son:

	• truncLen: corta las lecturas a esa longitud.
	• maxEE=c(2,2): controla la calidad permitida en términos de errores esperados. Puedes ser más estricto con c(1,1) si es necesario.
	• truncQ=2: corta la lectura en la primera base con calidad < Q2 (puede dejarse así).
	• maxN=0: descarta cualquier lectura con bases ambiguas (N).


Ahora filtraremos: 

```r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(260,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 

head(out)
```

Resultado

```
                reads.in reads.out
0h_R1.fastq.gz     2e+05    118317
24h_R1.fastq.gz    2e+05    124353
48h_R1.fastq.gz    2e+05    124133
9h_R1.fastq.gz     2e+05    123311
```

Siendo: 
reads.in: número de secuencias que tenía el archivo al inicio (≈ 200,000 en este caso) y reads.out: número de secuencias que pasaron el filtro de calidad y trimming.

Por ejemplo; En 0h_R1.fastq.gz había ~200,000 lecturas originales → 118,317 pasaron los 	filtros, (~59%–62%) común en datos reales.
  
## Modelado y Error  📈
El siguiente paso es modelar los errores de secuenciación.
DADA2 aprende el perfil de errores del secuenciador a partir de los reads.
Esto permite distinguir errores reales de variaciones biológicas reales (ASVs)

```r
errF <- learnErrors(filtFs, multithread = TRUE) # Aprender errores para las lecturas forward (R1)
```
Resultado: 127429640 total bases in 490114 reads from 4 samples will be used for learning the error rates (15 min).

```r
errR <- learnErrors(filtRs, multithread = TRUE) # Aprender errores para las lecturas reverse (R2)
```
Resultado: 98022800 total bases in 490114 reads from 4 samples will be used for learning the error rates (15 min), donde se observan menos bases (~98 millones) porque las lecturas reverse suelen ser más cortas.

Genera un gráfico con las tasas de error

```r
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
```

Dereplicar lecturas

```r
derepFs <- derepFastq(filtFs)
names(derepFs) <- sample.names
```
Con el modelo de errores se comparan las lecturas observadas vs las esperadas por error, y posteriormente identificamos las secuencias reales (ASVs) en cada muestra. Este paso reemplaza el antiguo enfoque de OTUs (clustering).

Inferencia de secuencias y merge con DADA para forward y reverse

```r
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)
```

Resultado: 
```
Sample - lecturas despues del filtrado - secuencias unicas detectadas
Sample 1 - 118317 reads in 25821 unique sequences.
Sample 2 - 124353 reads in 48299 unique sequences.
Sample 3 - 124133 reads in 52299 unique sequences.
Sample 4 - 123311 reads in 44430 unique sequences.
Sample 1 - 118317 reads in 25769 unique sequences.
Sample 2 - 124353 reads in 51070 unique sequences.
Sample 3 - 124133 reads in 48749 unique sequences.
Sample 4 - 123311 reads in 44051 unique sequences.
```

Una vez detectadas las ASVs en forward y reverse, ahora fusionamos los pares, DADA2 busca un solapamiento confiable entre cada par R1 y R2 para reconstruir el amplicón completo.

Fusiona los pares

```r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, minOverlap = 12)
```

Resultado

```
109886 paired-reads (in 2409 unique pairings) successfully merged out of 115858 (in 5288 pairings) input.
103372 paired-reads (in 9139 unique pairings) successfully merged out of 120441 (in 19792 pairings) input.
101310 paired-reads (in 10182 unique pairings) successfully merged out of 119420 (in 21311 pairings) input.
104586 paired-reads (in 7289 unique pairings) successfully merged out of 119821 (in 15948 pairings) input.
```

Construye tabla de secuencias (ASVs por muestra) y exporta las dimensiones de la tabla resultante

```r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# RESULTADO: [1]     4 22352
```

## Remoción secuencias de Quimeras
Las quimeras son artefactos comunes en PCR (combinaciones erróneas de fragmentos reales), DADA2 las identifica y elimina comparando los ASVs con otros más abundantes.

```r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```
Resultado: Identified 21987 bimeras out of 22352 input sequences (15 a 20 min).

Ver distribución de longitudes de ASVs después de quitar quimeras

```r
table(nchar(getSequences(seqtab.nochim)))
```

Resultado

```
260 280 281 382 383 387 400 401 402 403 404 405 406 407 408 409 410 411 414 415
1   2   3   1   2   2   2   1  16   3   5   5   1  24   5   2   1   2   1   2
419 420 421 422 423 424 425 426 427 428 429 442 443 448
2   1   1  17   1   1  13   8 163  54   2   3   2  16
```

Tabla de seguimiento de lecturas únicas que contiene datos sin ruidos

```
getN <- function(x) sum(getUniques(x))

track <- cbind(out_O, sapply(dadaFs, getN), rowSums(seqtab.nochim))
#track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN),
               rowSums(seqtab.nochim))
```

Ajustar nombres de columnas

```
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
#colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
#rownames(track) <- sample.names
head(track)
```
