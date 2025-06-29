Calculemos la diversidad Alfa (α)

```r
alpha_div <- estimate_richness(ps_filt, measures = c("Observed", "Shannon", "Simpson"))
# Agregar metadata si es necesario
alpha_div$Tratamiento <- sample_data(ps_filt)$Tratamiento
```
Graficar diversidad Shannon por grupo

```r
library(ggplot2)
ggplot(alpha_div, aes(x = Tratamiento, y = Shannon, fill = Tratamiento)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Diversidad Alfa (Shannon)", x = "Tratamiento", y = "Índice de Shannon") +
  theme(legend.position = "none")
```

Calculemos la diversidad Beta (β)

Transforma a proporciones relativas

```r
ps_rel <- transform_sample_counts(ps_filt, function(x) x / sum(x))
```

Calcula distancia Bray-Curtis

```r
ordu <- ordinate(ps_rel, method = "PCoA", distance = "bray")
```

Visualizar PCoA

```r
plot_ordination(ps_rel, ordu, color = "Tratamiento") +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "Diversidad Beta (Bray-Curtis - PCoA)")
```

Prueba estadística PERMANOVA (opcional)

Matriz de distancia
```
library(vegan)
dist_bc <- phyloseq::distance(ps_filt, method = "bray")

meta <- data.frame(sample_data(ps_filt)) # Metadata

adonis2(dist_bc ~ Tratamiento, data = meta)
```
