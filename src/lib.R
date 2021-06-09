#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("tidyr")
#install.packages("tibble")

setwd('C:/Users/sereb/Documents/Bioinformatics/minor_project/src')
getwd()

library(ggplot2)
library(dplyr)
library(tidyr)   # replace_na
library(tibble)  # column_to_rownames

###

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ChIPseeker")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", force=TRUE)
#BiocManager::install("clusterProfiler")
#BiocManager::install("GenomicFeatures", force=TRUE)

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)
library(GenomicFeatures)
library(org.Hs.eg.db)

###

DATA_DIR <- '../data/'
OUT_DIR <- '../images/'
