## Control de calidad con DADA2 🔍

Ahora, nos iremos a R para trabajar con dada2.

Instalaremos el paquete dada2

``` r
install.packages("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.16") 
``` 

Cargaremos las librerías necesarias

``` r
library(dada2)
library(Biostrings)
library(ShortRead)
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
head(sample.names)
```

Y finalmente graficaremos la calidad 

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

Filtramos: 

```r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(260,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 

head(out)
```

## Modelado y Error  📈

```r
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)


#Genera un gráfico con las tasas de error
```r
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
```

#Dereplicar lecturas
```r
derepFs <- derepFastq(filtFs)
names(derepFs) <- sample.names
```

#Inferencia de secuencias 

```r
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)
```
