#!/usr/bin/env R
#BioCManager.R
if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

BiocManager::install("metagenomeSeq")
