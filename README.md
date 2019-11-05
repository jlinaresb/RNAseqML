# Differential Gene Expression Analysis of RNA-seq Data Using Machine Learning for Cancer Research

This repository includes code based on the following chapter book:

*Liñares Blanco J., Gestal M., Dorado J., Fernandez-Lozano C. (2019) Differential Gene Expression Analysis of RNA-seq Data Using Machine Learning for Cancer Research. In: Tsihrintzis G., Virvou M., Sakkopoulos E., Jain L. (eds) Machine Learning Paradigms. Learning and Analytics in Intelligent Systems, vol 1. Springer, Cham*

DOI: https://doi.org/10.1007/978-3-030-15628-2_3

## TCGA data: BRCA and LUAD

RNASeq v2 data from BRCA and LUAD cohort have been download from 'http://gdac.broadinstitute.org/runs/stddata__latest/' URL. TCGA2STAT package version 1.2 have been used for this issue. For convenience, the objects have been stored in the folder /data.  

For the filtering of the genes, those related to the pathways of breast cancer, lung cancer as well as housekeeping genes and cell cycle genes were chosen. This gene signature was obtained through 'importing-genes.r' script.

## Abstract

Transcriptome analysis, as a tool for the characterization and understanding of phenotypic alterations in molecular biology, plays an integral role in the understanding of complex, multi-factorial and heterogeneous diseases such as cancer. Profiling of transcriptome is used for searching the genes that show differences in their expression level associated with a particular response. RNA-seq data allows researchers to study millions of short reads derived from an RNA sample using next-generation sequencing (NGS) methods. In general terms, such amount of data is difficult to understand and there is no optimal analysis pipeline for every single analysis. Classical statistical approaches are provided in different R packages (i.e. DESeq or edgeR packages). In medicine, a Machine Learning algorithm can be used for the differential expression analysis of a particular response (i.e. sick versus healthy patients) selecting the genes that are more relevant for discriminating both health outcomes, considering biological pathway information, gene relations or using integrative approaches in order to include all the information available from different curated data-sources. The main aim of our proposal is to practically address Machine Learning based approach for gene expression analysis using RNA-seq data for cancer research within the R framework and to compare it with a classical gene expression analysis approach.

## Getting Started

There are two RNAseq analysis in this repo, a classical and a machine-learning based analysis. Each has its own R script. But before you run the code, make sure you have the following:

### Prerequisites:

The packages we've used:

```{r}
install.packages(c("ggplot2", "mlr", "dplyr", "vegan", "parallelMap"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("tweeDEseq", "edgeR"))


library(devtools)
install_github("vqv/ggbiplot")
```
so, make sure that the following works!

```{r}
require(tweeDEseq)
require(edgeR)
require(vegan)
require(ggplot2)
require(ggbiplot)
require(mlr)
require(dplyr)
require(parallelMap)
```
#### Some system dependencies

Please, be sure that you have the following packages in Ubuntu:

```sh
sudo apt-get update
sudo apt-get install libxml2-dev
sudo apt-get install r-cran-xml
sudo apt-get install libssl-dev libcurl4-openssl-dev
```

### How to run

Clone this repository:

```{bash}
git clone https://github.com/jlinaresb/RNAseqML.git
```

navigate to the RNAseqML/R folder and execute the following command (choose your OS):

#### In linux environments

```{bash}
Rscript run-classical-analysis.r
Rscript run-MachineLearning-analysis.r
```

#### In windows environments

```{bash}
'C:\Program Files\R\R-3.5.2\bin\Rscript.exe' run-classical-analysis.r
'C:\Program Files\R\R-3.5.2\bin\Rscript.exe' run-MachineLearning-analysis.r
```

## Citation:

If you have found our code and chapter book useful, we kindly ask you to cite our work.

```tex
@Inbook{LiñaresBlanco2019,
author="Li{\~{n}}ares Blanco, Jose
and Gestal, Marcos
and Dorado, Juli{\'a}n
and Fernandez-Lozano, Carlos",
editor="Tsihrintzis, George A.
and Virvou, Maria
and Sakkopoulos, Evangelos
and Jain, Lakhmi C.",
title="Differential Gene Expression Analysis of RNA-seq Data Using Machine Learning for Cancer Research",
bookTitle="Machine Learning Paradigms: Applications of Learning and Analytics in Intelligent Systems",
year="2019",
publisher="Springer International Publishing",
address="Cham",
pages="27--65",
isbn="978-3-030-15628-2",
doi="10.1007/978-3-030-15628-2_3",
url="https://doi.org/10.1007/978-3-030-15628-2_3"
}
```

## Questions?

If you have any questions, please feel free to contact (j.linares@udc.es).

