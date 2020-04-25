#!/usr/bin/env Rscript
# de.R
library(tximport)
library(readr)
library(DESeq2)
tx2gene <- read.csv("tx2gene.csv")
head(tx2gene)
samples <- read.csv("/scratch/SampleDataFiles/Samples.csv", header=TRUE)
head(samples)
files <- file.path("quant", samples$Sample, "quant.sf")
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ Menthol + Vibrio)
dds$Vibrio <- relevel(dds$Vibrio, ref = "Control")
dds$Menthol <- relevel(dds$Menthol, ref = "Control")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
padj <- .05
minLog2FoldChange <- .5
dfAll <- data.frame()
# Get all DE results except Intercept, and "flatten" into a single file.
for (result in resultsNames(dds)){
        if(result != 'Intercept'){
                    res <- results(dds, alpha=.05, name=result)
        dfRes <- as.data.frame(res)
                dfRes <- subset(subset(dfRes, select=c(log2FoldChange, padj)))
                dfRes$Factor <- result
                        dfAll <- rbind(dfAll, dfRes)
                    }
}
head(dfAll)
write.csv(dfAll, file="dfAll.csv")

data <- read.csv("dfAll.csv",header=T)
colnames(data) <- c('ko','log2FoldChange','padj','Factor')
data <- subset(data, padj < 0.05)

pathway <- read.table("/scratch/SampleDataFiles/Annotation/path.txt", sep ='\t')
colnames(pathway) <- c("ko","pathway")

dfAll_pathway <- merge(data,pathway)
descr <- read.table("/scratch/SampleDataFiles/Annotation/ko",sep='\t')
colnames(descr) <- c("pathway","description")
pathway_desc <- merge(dfAll_pathway, descr)
write.csv(pathway_desc, file="deAnnotated.csv")

#```{R}
library(knitr)
de_anno <-read.csv("deAnnotated.csv",header=T)
kable(de_anno)
#```
