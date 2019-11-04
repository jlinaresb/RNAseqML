DOI: https://doi.org/10.1007/978-3-030-15628-2_3


# Abstract

Transcriptome analysis, as a tool for the characterization and understanding of phenotypic alterations in molecular biology, plays an integral role in the understanding of complex, multi-factorial and heterogeneous diseases such as cancer. Profiling of transcriptome is used for searching the genes that show differences in their expression level associated with a particular response. RNA-seq data allows researchers to study millions of short reads derived from an RNA sample using next-generation sequencing (NGS) methods. In general terms, such amount of data is difficult to understand and there is no optimal analysis pipeline for every single analysis. Classical statistical approaches are provided in different R packages (i.e. DESeq or edgeR packages). In medicine, a Machine Learning algorithm can be used for the differential expression analysis of a particular response (i.e. sick versus healthy patients) selecting the genes that are more relevant for discriminating both health outcomes, considering biological pathway information, gene relations or using integrative approaches in order to include all the information available from different curated data-sources. The main aim of our proposal is to practically address Machine Learning based approach for gene expression analysis using RNA-seq data for cancer research within the R framework and to compare it with a classical gene expression analysis approach.

# Required packages

```{r}
install.packages(c("ggplot2", "mlr", "dplyr", "vegan"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("tweeDEseq", edgeR))


library(devtools)
install_github("vqv/ggbiplot")
```


```{r }
require(tweeDEseq)
require(edgeR)
require(vegan)
require(ggplot2)
require(ggbiplot)
require(mlr)
require(dplyr)
```


# How to run

```{bash }
Rscript ~/RNAseqML/R/run-classical-analysis.r
Rscript ~/RNAseqML/R/run-MachineLearning-analysis.r
```