## Este pipeline describe como construir un objeto phyloseq a partir de datos de secuenciación metagenómica ya preprocesados con DADA2.

Comenzamos cargando las librerias para el análisis filogenético y estructural de las comunidades microbianas.

```r
library(tidyverse)
library(vegan)
library(ggplot2)
library(DECIPHER)
library(phangorn)
library(phyloseq)
```

Cargaremos/importaremos los objetos guardados después del pipeline de DADA2 (ASVs, tablas y taxonomía)

```r
load(file = "data_ok.RData")
```

Leemos la metadata asociada a las muestras.

```r
metadata <- read_csv("Metadata.csv") %>%
  mutate_at(vars(depth), as.character) %>%
  column_to_rownames("seqR1")
```

Generamos un árbol filogenético a partir de las secuencias ASV.

```r
ASVs.nochim <- DNAStringSet(colnames(seqtab.nochim))
names(ASVs.nochim) <- paste0("ASV", 1:ncol(seqtab.nochim))

alignment <- AlignSeqs(ASVs.nochim, anchor = NA, processors = 30)
phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)

fit <- pml(treeNJ, data = phang.align)
fitGTR_green <- update(fit, k = 4, inv = 0.2)

save(fitGTR_green, file = "tree.green.ok.RData")
```

Construiremos el objeto phyloseq integrando secuencias, taxonomía, árbol y metadata.

```r
load("tree.green.ok.RData")

rownames(seqtab.nochim) <- gsub(".fastq.gz", "", rownames(seqtab.nochim))
tmp.seqtab <- seqtab.nochim
colnames(tmp.seqtab) <- names(ASVs.nochim)

tmp.taxa <- taxa
rownames(tmp.taxa) <- names(ASVs.nochim)

ps.green.nochim <- phyloseq(
  otu_table(tmp.seqtab, taxa_are_rows = FALSE),
  sample_data(metadata),
  tax_table(tmp.taxa),
  refseq(ASVs.nochim),
  phy_tree(fitGTR_green$tree)
)

save(ps.green.nochim, file = "phyloseq_green_ok.RDATA")
```

Enraizamos el árbol filogenético para análisis de diversidad filogenética.

```r
load("phyloseq_green_ok.RDATA")
ps <- ps.green.nochim

set.seed(1)
phy_tree(ps) <- root(phy_tree(ps), sample(taxa_names(ps), 1), resolve.root = TRUE)

is.rooted(phy_tree(ps))  # Verifica si el árbol está enraizado

save(ps, file = "phyloseq_green_root_ok.RDATA")
```
















