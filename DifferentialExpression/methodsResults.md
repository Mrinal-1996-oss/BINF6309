## Methods

To perform differential expression analysis, we conducted two steps;
estimate the relative abundance of transcripts, and apply statistical
models to test for differential expression between treatment groups. To
estimate the relative abundance, we determine the numbers NGS reads
match a given gene within a genome. To accomplish, we utilized three
packages Salmon, DESeq2, and tximport.

The first process in differential expression analysis starts with
quantifying transcript abundance from RNA-seq reads (Patro et al). As a
transcriptome-wide quantifier to adjusts fragment GC contents bias,
which helps the accuracy of abundance estimates and following
differential expression analysis (Patro et al). First, we took de-novo
transcriptome and built a salmon index of 25 kmer length. Second, we
align Aip samples to AipIndex using Salmon.

The second process utilizes tximport package to import the Salmon
abundance estimate (Sonsen et al). In order to utilize Salmon input in
tximport, we create a table mapping transcripts to genes. Here, used
annotation files form /scratch/SampleDataFiles/Annotation to build the
table in R.

The third process uses DESeq2 to conduct statistical tests to identify
differentially expressed genes (love et al). This is carried out by
importing salmon alignments into DESeq2 and perform differential
expression analysis.

## Results

``` r
library(knitr)
de_anno <- read.csv("deAnnotated.csv", header= T)
kable(de_anno)
```

|  X | pathway      | ko        | log2FoldChange |      padj | Factor                        | description                               |
| -: | :----------- | :-------- | -------------: | --------: | :---------------------------- | :---------------------------------------- |
|  1 | path:ko03015 | ko:K04354 |     \-2.068948 | 0.0273493 | Menthol\_Menthol\_vs\_Control | mRNA surveillance pathway                 |
|  2 | path:ko04071 | ko:K04354 |     \-2.068948 | 0.0273493 | Menthol\_Menthol\_vs\_Control | Sphingolipid signaling pathway            |
|  3 | path:ko04111 | ko:K04354 |     \-2.068948 | 0.0273493 | Menthol\_Menthol\_vs\_Control | Cell cycle - yeast                        |
|  4 | path:ko04144 | ko:K17920 |       1.118743 | 0.0000011 | Menthol\_Menthol\_vs\_Control | Endocytosis                               |
|  5 | path:ko04151 | ko:K04354 |     \-2.068948 | 0.0273493 | Menthol\_Menthol\_vs\_Control | PI3K-Akt signaling pathway                |
|  6 | path:ko04152 | ko:K04354 |     \-2.068948 | 0.0273493 | Menthol\_Menthol\_vs\_Control | AMPK signaling pathway                    |
|  7 | path:ko04261 | ko:K04354 |     \-2.068948 | 0.0273493 | Menthol\_Menthol\_vs\_Control | Adrenergic signaling in cardiomyocytes    |
|  8 | path:ko04390 | ko:K04354 |     \-2.068948 | 0.0273493 | Menthol\_Menthol\_vs\_Control | Hippo signaling pathway                   |
|  9 | path:ko04391 | ko:K04354 |     \-2.068948 | 0.0273493 | Menthol\_Menthol\_vs\_Control | Hippo signaling pathway - fly             |
| 10 | path:ko04530 | ko:K04354 |     \-2.068948 | 0.0273493 | Menthol\_Menthol\_vs\_Control | Tight junction                            |
| 11 | path:ko04728 | ko:K04354 |     \-2.068948 | 0.0273493 | Menthol\_Menthol\_vs\_Control | Dopaminergic synapse                      |
| 12 | path:ko05142 | ko:K04354 |     \-2.068948 | 0.0273493 | Menthol\_Menthol\_vs\_Control | Chagas disease (American trypanosomiasis) |
| 13 | path:ko05160 | ko:K04354 |     \-2.068948 | 0.0273493 | Menthol\_Menthol\_vs\_Control | Hepatitis C                               |
| 14 | path:ko05165 | ko:K04354 |     \-2.068948 | 0.0273493 | Menthol\_Menthol\_vs\_Control | Human papillomavirus infection            |

## References

<div id="refs" class="references">

<div id="ref-Love">

Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. “Moderated
Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2.”
*Genome Biology* 15 (12): 550–50.
<https://doi.org/10.1186/s13059-014-0550-8>.

</div>

<div id="ref-Patro">

Patro, Rob, Geet Duggal, Michael I. Love, Rafael A. Irizarry, and Carl
Kingsford. 2017. “Salmon Provides Fast and Bias-Aware Quantification of
Transcript Expression.” *Nature Methods* 14 (4): 417–19.
<https://doi.org/10.1038/nmeth.4197>.

</div>

<div id="ref-Soneson">

Soneson, Charlotte, Michael I. Love, and Mark D. Robinson. 2015.
“Differential Analyses for RNA-Seq: Transcript-Level Estimates Improve
Gene-Level Inferences.” *F1000Research* 4 (December): 1521–1.
<https://www.ncbi.nlm.nih.gov/pubmed/26925227>.

</div>

</div>
