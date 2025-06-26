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
sq <- getSequences(fn)
tdf <- read.csv(txfn, sep = "\t", header = TRUE)
tax <- tdf[, 2]
names(tax) <- tdf[, 1]

identical(names(sq), names(tax))  # Debe devolver TRUE
```

Profundidad taxonómica

Se evalúa la cantidad de niveles taxonómicos por entrada (Greengenes2 maneja 7 niveles: dominio a especie):

```r
taxes <- strsplit(tax, "; ")
tax.depth <- sapply(taxes, length)
table(tax.depth)
```

## 🧬 Asignación Taxonómica con Clasificador Personalizado Greengenes2 (DADA2)

Este bloque de código realiza dos tareas principales:

**1. Limpieza y reestructuración de las anotaciones taxonómicas provenientes de Greengenes2 para generar un clasificador en formato `.fa.gz`.**

**2. Asignación de taxonomía a las ASVs utilizando `assignTaxonomy()` de `dada2`.**

# 1.Limpieza de niveles taxonómicos y construcción del clasificador

El objetivo es corregir y estandarizar los niveles taxonómicos antes de entrenar el clasificador:

Eliminar el nombre del género duplicado en la especie

```r
for(i in seq(length(taxes))) {
  gen <- taxes[[i]][[6]]
  gen <- substr(gen, 4, nchar(gen))
  taxes[[i]][[7]] <- gsub(gen, "", taxes[[i]][[7]])
  taxes[[i]][[7]] <- gsub("__ ", "__", taxes[[i]][[7]])
}
```

Detectar niveles no asignados

```r
tax_pre <- c("d__", "p__", "c__", "o__", "f__", "g__", "s__")
is.unassigned <- sapply(taxes, function(tx) {
  tx == tax_pre
}) |> t()
```

Determinar profundidad de clasificación válida

```r
tax.depth <- apply(is.unassigned, 1, function(isu) {
  min(which(isu) - 1L, 7L)
})
```

Generar ID taxonómicos válidos

```r
tax.ids <- sapply(seq_along(taxes), function(i) {
  td <- tax.depth[[i]]
  id.str <- paste(taxes[[i]][1:td], collapse = ";")
  paste0(id.str, ";")
})
names(tax.ids) <- names(taxes)
```

Crear archivo FASTA con anotaciones corregidas

```r
sq.out <- sq
names(sq.out) <- tax.ids
writeFasta(sq.out, "greengenes2_trainset.fa.gz", compress = TRUE)
```

Validación rápida del clasificador

```r
dada2:::tax.check("greengenes2_trainset.fa.gz",
                  fn.test = system.file("extdata", "ten_16s.100.fa.gz", package = "dada2"))
```

# 2. Asignación de taxonomía a las ASVs

```r
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

Dimensiones de la tabla

```r
dim(taxa_print)
```

Exportar objetos generados durante el preprocesamiento para uso posterior

```r
save(errF, dadaFs, dadaRs, seqtab.nochim, taxa,
     file = "data_ok.RData") # Archivo binario con objetos R para reutilización
```

Exportar resultados en formato CSV

```r
write.csv(taxa, "taxonomy.csv") # Taxonomía asignada a cada ASV
write.csv(seqtab.nochim, "table.csv") # Matriz ASV por muestra (sin quimeras)
write.csv(track, "stats.csv") # Estadísticas de procesamiento por muestra
```

