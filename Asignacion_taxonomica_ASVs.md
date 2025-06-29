## 🧬 Preparación del clasificador taxonómico Greengenes2

En esta sección del flujo de trabajo, se descarga y procesa un clasificador taxonómico basado en **Greengenes2 (2022.10)**, que será utilizado con `assignTaxonomy()` de DADA2 para la asignación de taxonomía.

Descarga de archivos del clasificador

Se obtienen los archivos necesarios desde el repositorio oficial de Greengenes2:

```r
download.file("http://ftp.microbio.me/greengenes_release/current/2022.10.backbone.full-length.fna.qza", 
              "2022.10.backbone.full-length.fna.qza")
download.file("http://ftp.microbio.me/greengenes_release/current/2022.10.backbone.tax.qza", 
              "2022.10.backbone.tax.qza")
```

Descompresión de archivos

```r
unzip("2022.10.backbone.full-length.fna.qza")
unzip("2022.10.backbone.tax.qza")
```

Se definen las rutas de los archivos extraídos

```r
fn <- "a53d9300-5c5c-4774-a2e8-a5e23904f1ae/data/dna-sequences.fasta"
txfn <- "c16a953c-f24d-4d14-927c-40d90ced395e/data/taxonomy.tsv"
```

Se leen las secuencias y la tabla de taxonomía. 

Luego se verifica que las secuencias y sus nombres coincidan con los identificadores de la taxonomía

```r
sq <- getSequences(fn) # Carga las secuencias de referencia
tdf <- read.csv(txfn, sep = "\t", header = TRUE) # Lee la tabla de taxonomía (formato TSV)
```

Asigna los nombres de la taxonomía a las secuencias

```r
tax <- tdf[, 2] # 
names(tax) <- tdf[, 1]
```

Verifica que las secuencias y las taxonomías coinciden en orden

```r
identical(names(sq), names(tax))  # Debe devolver TRUE
```

## 🧬 Asignación Taxonómica con Clasificador Personalizado Greengenes2 (DADA2)

Este bloque de código realiza dos tareas principales:

**1. Limpieza y reestructuración de las anotaciones taxonómicas provenientes de Greengenes2 para generar un clasificador en formato `.fa.gz`.**

**2. Asignación de taxonomía a las ASVs utilizando `assignTaxonomy()` de `dada2`.**

# 1.Limpieza de niveles taxonómicos y construcción del clasificador

El objetivo será corregir y estandarizar los niveles taxonómicos antes de entrenar el clasificador:

Primero elimina el nombre del género duplicado en la especie

```r
taxes <- strsplit(tax, "; ")
for(i in seq(length(taxes))) {
  gen <- taxes[[i]][[6]]
  gen <- substr(gen, 4, nchar(gen)) # Elimina el prefijo "g__"
  taxes[[i]][[7]] <- gsub(gen, "", taxes[[i]][[7]]) # Elimina repeticiones del género
  taxes[[i]][[7]] <- gsub("__ ", "__", taxes[[i]][[7]]) # Limpia espacios
}
```

Después, detecta niveles no asignados

```r
tax_pre <- c("d__", "p__", "c__", "o__", "f__", "g__", "s__") 
is.unassigned <- sapply(taxes, function(tx) {
  tx == tax_pre
}) |> t()
```
table(tax.depth)
```
Determina profundidad de clasificación válida por secuencia

```r
tax.depth <- apply(is.unassigned, 1, function(isu) {
  min(which(isu) - 1L, 7L)
})

table(tax.depth)
```

Genera ID taxonómicos

```r
tax.ids <- sapply(seq_along(taxes), function(i) {
  td <- tax.depth[[i]]
  id.str <- paste(taxes[[i]][1:td], collapse = ";")
  paste0(id.str, ";")
})
names(tax.ids) <- names(taxes)
```

Crea el archivo FASTA con anotaciones corregidas compatible con DADA2

Comienza asociando las secuencias con su respectiva etiqueta y guardando el archivo para el uso con assignTaxonomy

```r
sq.out <- sq
names(sq.out) <- tax.ids
writeFasta(sq.out, "greengenes2_trainset.fa.gz", compress = TRUE)
```

Validación rápida del clasificador (opcional)

```r
dada2:::tax.check("greengenes2_trainset.fa.gz",
                  fn.test = system.file("extdata", "ten_16s.100.fa.gz", package = "dada2"))
```

# 2. Asignación de taxonomía a las ASVs

Cargar objetos generados por DADA2

```r
load("errF, errR, dadaFs, dadaRs, seqtab.nochim")
```

Ejecuta

```
classifier <- "greengenes2_trainset.fa.gz" # Clasificador limpio y entrenado para assignTaxonomy()

taxa <- assignTaxonomy(seqtab.nochim,
                       classifier,
                       multithread = TRUE,
                       tryRC = TRUE)
```

## Exploración de resultados

Elimina nombres de fila para impresión amigable

```r
taxa_print <- taxa
rownames(taxa_print) <- NULL
```

Primeras asignaciones taxonómicas

```r
head(taxa_print)
```

Dimensiones de la tabla, filas ASV's y columnas nivel taxonómico

```r
dim(taxa_print)
```

Exporta los objetos generados durante el preprocesamiento para uso posterior

```r
save(errF, dadaFs, dadaRs, seqtab.nochim, taxa,
     file = "data_ok.RData") # Archivo binario con objetos R para reutilización
```

Por ultimo, exporta los resultados en formato CSV

```r
write.csv(taxa, "taxonomy.csv") # Taxonomía asignada a cada ASV
write.csv(seqtab.nochim, "table.csv") # Matriz ASV por muestra (sin quimeras)
write.csv(track, "stats.csv") # Estadísticas de procesamiento por muestra
```
