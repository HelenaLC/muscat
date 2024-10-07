<img src="inst/extdata/muscat.png" width="200" align="right"/> 

**`muscat` (**Mu**lti-sample **mu**lti-group **sc**RNA-seq **a**nalysis **t**ools )**

...provides methods for *Differential State* (DS) analyses in scRNA-seq data  
with multiple samples, groups, and (cell)-subpopulations, as elaborated in:

> Crowell HL, Soneson C\*, Germain P-L\*,  
Calini D, Collin L, Raposo C, Malhotra D & Robinson MD:  
"*muscat* detects subpopulation-specific state transitions from  
multi-sample multi-condition single-cell transcriptomics data"  
*Nature Communications* **11**, 6077 (2020)  
[DOI: 10.1038/s41467-020-19894-4](https://doi.org/10.1038/s41467-020-19894-4)

*These authors contributed equally.

### installation

`muscat` is available through Bioconductor, and
can be installed using the following commands:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("muscat")
```

### quick guide

Let `sce` be a [`SingleCellExperiment`](https://www.bioconductor.org/packages/SingleCellExperiment.html) object with cell metadata (`colData`) columns

1. `"sample_id"` specifying unique sample identifiers (e.g., PeterPan1, Nautilus7, ...)
2. `"group_id"` specifying each sample's experimental condition (e.g., reference/stimulated, healthy/diseased, ...)
3. `"cluster_id"` specifying subpopulation (cluster) assignments (e.g., B cells, dendritic cells, ...)

Aggregation-based methods come down to the following simple commands: 

```r
# compute pseudobulks (sum of counts)
pb <- aggregateData(sce, 
    assay = "counts", fun = "sum",
    by = c("cluster_id", "sample_id"))
    
# run pseudobulk (aggregation-based) DS analysis
ds_pb <- pbDS(pb, method = "edgeR")
```

Mixed models can be run directly on cell-level measurements, e.g.:

```r
ds_mm <- mmDS(sce, method = "dream")
```

For details, please see the package vignettes.

### differential detection

`muscat` also supports testing for differential detection as proposed in 

> Gilis J, Perin L, Malfait M, Van den Berge K,  
Assefa AT, Verbist B, Risso D, and Clement L:  
Differential detection workflows for  
multi-sample single-cell RNA-seq data.  
*bioRxiv* (2023). [DOI: 10.1101/2023.12.17.572043](https://doi.org/10.1101/2023.12.17.572043)

Key alterations to the commands above are highlighted below (!!!), 
however, we recommend users consult the corresponding publication 
and package vignette for more details.

```r
# sum binarized counts
pb <- aggregateData(sce, 
    assay = "counts", 
    fun = "num.detected", # !!!
    by = c("cluster_id", "sample_id"))
# test for differential detection
dd <- pbDD(pb) # or..
dd <- pbDS(pb, method = "DD")
```