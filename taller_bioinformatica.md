## TIPS y comando extras que se utilizaron durante el taller

# Script para guardar modelos de error DADA2 (errF, errR)
Ruta donde guardar

```r
ruta_guardado <- "modelos_errores.RData"
```

Guardar los objetos

```r
save(errF, errR, file = ruta_guardado)
cat("✅ Modelos errF y errR guardados en:", ruta_guardado, "\n")
```

# Script para cargarlos errF, errR previamente guardados
Ruta del archivo errF, errR en .RData (ajusta la ruta si es necesario)

```r
ruta_archivo <- "modelos_errores.RData"
```

Cargar el archivo

```r
load(ruta_archivo)
```

Verifica que se cargaron

```r
print(ls())
cat("✅ Modelos cargados: errF y errR disponibles para usar\n")
```

Utilizar para el paso siguiente
#dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
#dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

# Para guardar los gráficos

```r
png("plotErrors_forward.png", width = 800, height = 600)
plotErrors(errF, nominalQ = TRUE)
dev.off()

png("plotErrors_reverse.png", width = 800, height = 600)
plotErrors(errR, nominalQ = TRUE)
dev.off()
```

# Script para cargarlos dadaFs y dadaRs previamente guardados

Ruta del archivo con los objetos

```r
ruta_dada <- "salida_dadaFs_dadaRs.RData"
```

Cargar los objetos

```r
load(ruta_dada)
```

Confirmar contenido

```r
print(ls())
cat("✅ Objetos cargados: dadaFs y dadaRs están listos para usar\n")
```

# Script para cargarlos ruta_nochim previamente guardada (la tabla sin quimeras)

```r
ruta_nochim <- "tabla_no_quimeras.RData"
```

Cargar el objeto

```r
load(ruta_nochim)
```

Verifica que está cargado

```r
print(ls())
cat("✅ Objeto 'seqtab.nochim' cargado y listo para usar\n")
```

# Para revisión rápida o uso en otro software se puede exportar como .csv

```r
write.csv(seqtab.nochim, "tabla_no_quimeras.csv")
```



