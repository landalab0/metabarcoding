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

# Contiene: errF, errR, dadaFs, dadaRs, seqtab.nochim, taxa
```

Leemos la metadata asociada a las muestras.

```r
metadata <- read_csv("Metadata.csv") %>%
  mutate_at(vars(depth), as.character) %>%
  column_to_rownames("seqR1")
```

🌳 Generamos un árbol filogenético a partir de las secuencias ASV.

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

# Asegurarse que los nombres de muestra coincidan con el metadata
rownames(seqtab.nochim) <- gsub(".fastq.gz", "", rownames(seqtab.nochim))

# Renombrar columnas (ASVs) y filas (taxonomía) con IDs cortos
tmp.seqtab <- seqtab.nochim
colnames(tmp.seqtab) <- names(ASVs.nochim)

tmp.taxa <- taxa
rownames(tmp.taxa) <- names(ASVs.nochim)

# Construir el objeto phyloseq
ps.green.nochim <- phyloseq(
  otu_table(tmp.seqtab, taxa_are_rows = FALSE),
  sample_data(metadata),
  tax_table(tmp.taxa),
  refseq(ASVs.nochim),
  phy_tree(fitGTR_green$tree)
)

# Guardalo en el objeto phyloseq
save(ps.green.nochim, file = "phyloseq_green_ok.RDATA")
```

🌱 Enraizamos el árbol filogenético para análisis de diversidad filogenética.

```r
load("phyloseq_green_ok.RDATA")
ps <- ps.green.nochim

set.seed(1)
phy_tree(ps) <- root(phy_tree(ps), sample(taxa_names(ps), 1), resolve.root = TRUE)

is.rooted(phy_tree(ps))  # Verifica si el árbol está enraizado

save(ps, file = "phyloseq_green_root_ok.RDATA")
```

Visualiza de la composición y curado del objeto phyloseq

Abundancia relativa por Familia

```r
otu_tab <- otu_table(seqtab.nochim, taxa_are_rows = FALSE)
tax_tab <- tax_table(taxa)
ps <- phyloseq(otu_tab, tax_tab)

colnames(tax_table(ps)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

ps_family <- tax_glom(ps, taxrank = "Family")

plot_bar(ps_family, fill = "Family") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  ggtitle("Abundancia relativa por Familia") +
  scale_y_continuous(labels = scales::percent_format(scale = 1))
```

Filtra ASVs no asignados o poco abundantes

```r
ps_filt <- subset_taxa(ps, !is.na(Phylum) & Phylum != "")
ps_filt <- filter_taxa(ps_filt, function(x) sum(x) > 10, TRUE)
```

Verifica reducción

```r
ps
ps_filt
```
Archivo de metadata experimental

```
sample.names <- c("0h", "24h", "48h", "9h")

metadata <- data.frame(
  SampleID = sample.names,
  Tiempo = c(0, 24, 48, 9),
  Tratamiento = c("Control", "Fermentación", "Fermentación", "Fermentación")
)
rownames(metadata) <- metadata$SampleID

write.csv(metadata, "metadata.csv", row.names = TRUE)
```

Integra metadata al objeto phyloseq

```r
metadata <- read.csv("metadata.csv", row.names = 1)
sample_data(ps_filt) <- sample_data(metadata)

ps_filt # Confirmar
```

Grafica abundancia por grupo (Phylum)

plot_bar(ps_filt, fill = "Phylum") +
  facet_wrap(~ Tratamiento) +
  theme_minimal() +
  labs(title = "Composición por Tratamiento (Phylum)")
