# GENE EXPRESSION CORRELATES OF CANCER MUTATION SIGNATURES

The aim of this project is to discover the sources and/or consequences of unknown mutational signatures. 

## prerequisites

This project was completed using R version 4.0.2


```bash
library(vcfR)
library(limma)
library(statmod)
library(devtools)
library(MAST)
library(dplyr)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(plyr)
library(AnnotationFilter)
library(fgsea)
library(data.table)
library(ggplot2)
library(BiocParallel)
library(formattable)
library(tidyr)
library(gridExtra)
library(xtable)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(ggbiplot)
library(readxl)
```



## concatVCF.R

`concatVCF.R` is a function that reads in several bgzipped VCF files from a folder and concatenates them so they may be input into SigProfilerJulia

[SigProfilerJulia]: https://bitbucket.org/bbglab/sigprofilerjulia/src/master/


## gsea.R

`gsea.R` performs linear regression to find correlation between genes and the extracted mutational signatures, and gene set enrichment analysis on the output of this linear regression

## plots.R

'plots.R' provides code to create Principal Component Analysis plots

## TCGAbiolinks.R

'TCGAbiolinks.R' provides code that allows the user to use the package TCGAbiolinks, which performs integrative analysis with GDC data such as the TCGA data used in this project

## GSEAplot.R

'GSEAplot.R' provides code allowing the user to perform gene set enrichment, and create a dumbbell plot of the top ten overexpressed and top ten underexpressed pathways 

#### Installing git

git config --global user.name "sheridanemily"
git config --global user.email "sheridanemilyanne@gmail.com"

git init
git config user.name "sheridanemily"
git config user.email "sheridanemilyanne@gmail.com"
git status



#### Stage and commit files to the local repo

(base) a@e:~/MSc_Project$ git log
commit a55ee6ac0095877739fc23610c37b0817e9f9f7c (HEAD -> master)
Author: sheridanemily <sheridanemilyanne@gmail.com>
Date:   Tue Nov 26 13:38:48 2019 +0000

    Add initial version of thesis code.
    



https://github.com/sheridanemily/MSc_Project.git


