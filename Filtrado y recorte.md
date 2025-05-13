
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
