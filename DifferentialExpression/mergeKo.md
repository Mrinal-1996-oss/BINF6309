``` r
#!/usr/bin/env Rscript
#mergeKo.R
# Load the kable library to display formatted tables
library(knitr)
# Load BLAST results as a table
blastFile <- "/scratch/SampleDataFiles/Annotation/transcriptBlast.txt"
keggFile <- "/scratch/SampleDataFiles/Annotation/kegg.txt"
koFile <- "/scratch/SampleDataFiles/Annotation/ko.txt"

blast <- read.table(blastFile, sep="\t", header=FALSE)
# Set column names to match fields selected in BLAST
colnames(blast) <- c("trans", "sp", "qlen", "slen", "bitscore", 
                                          "length", "nident", "pident", "evalue", "ppos")
# Calculate the percentage of identical matches relative to subject length
blast$cov <- blast$nident/blast$slen
# Filter for at least 50% coverage of subject(SwissProt) sequence
blast <- subset(blast, cov > .5)
# Check the blast table
kable(head(blast))
```

|      | trans                      | sp     | qlen | slen | bitscore | length | nident | pident | evalue |  ppos |       cov |
| ---- | :------------------------- | :----- | ---: | ---: | -------: | -----: | -----: | -----: | -----: | ----: | --------: |
| 114  | TRINITY\_DN26\_c0\_g1\_i1  | P62489 |  462 |  172 |      261 |    135 |    121 |  89.63 |      0 | 97.04 | 0.7034884 |
| 115  | TRINITY\_DN26\_c0\_g1\_i1  | Q7ZW41 |  462 |  172 |      259 |    135 |    120 |  88.89 |      0 | 97.04 | 0.6976744 |
| 126  | TRINITY\_DN26\_c0\_g2\_i1  | P62489 |  460 |  172 |      261 |    135 |    121 |  89.63 |      0 | 97.04 | 0.7034884 |
| 127  | TRINITY\_DN26\_c0\_g2\_i1  | Q7ZW41 |  460 |  172 |      259 |    135 |    120 |  88.89 |      0 | 97.04 | 0.6976744 |
| 1529 | TRINITY\_DN122\_c1\_g1\_i1 | Q5R988 |  412 |  183 |      223 |    135 |    102 |  75.56 |      0 | 91.11 | 0.5573770 |
| 1530 | TRINITY\_DN122\_c1\_g1\_i1 | Q8BU31 |  412 |  183 |      219 |    135 |    101 |  74.81 |      0 | 88.89 | 0.5519126 |

``` r
# Load SwissProt to KEGG as a table
kegg <- read.table(keggFile, sep="\t", header=FALSE)
# Set the Swissprot to KEGG column names
colnames(kegg) <- c("sp", "kegg")
# Remove the up: prefix from sp column
kegg$sp <- gsub("up:", "", kegg$sp)
# Check the kegg table
kable(head(kegg))
```

| sp     | kegg       |
| :----- | :--------- |
| P62489 | rno:117017 |
| Q7ZW41 | dre:324088 |
| Q8BU31 | mmu:72065  |
| P61227 | rno:170923 |
| P10114 | hsa:5911   |
| P28074 | hsa:5693   |

``` r
# Merge BLAST and SwissProt-to-KEGG
blastKegg <- merge(blast, kegg)
# Check the merged table
kable(head(blastKegg))
```

| sp     | trans                       | qlen | slen | bitscore | length | nident | pident | evalue |  ppos |       cov | kegg                  |
| :----- | :-------------------------- | ---: | ---: | -------: | -----: | -----: | -----: | -----: | ----: | --------: | :-------------------- |
| A0AJB9 | TRINITY\_DN9562\_c0\_g1\_i1 | 1944 |  399 |      401 |    376 |    208 |  55.32 |      0 | 73.94 | 0.5213033 | lwe:lwe1683           |
| A0AKK8 | TRINITY\_DN8143\_c2\_g1\_i1 | 1805 |  295 |      348 |    278 |    165 |  59.35 |      0 | 77.70 | 0.5593220 | lwe:lwe2122           |
| A0ALL5 | TRINITY\_DN9485\_c0\_g1\_i1 | 1847 |  504 |      571 |    502 |    286 |  56.97 |      0 | 73.90 | 0.5674603 | lwe:lwe2479           |
| A0AQ71 | TRINITY\_DN9157\_c0\_g1\_i2 | 1369 |  206 |      213 |    202 |    104 |  51.49 |      0 | 70.79 | 0.5048544 | dsi:Dsimw501\_GD21160 |
| A0AQ71 | TRINITY\_DN9157\_c0\_g1\_i1 | 1393 |  206 |      213 |    202 |    104 |  51.49 |      0 | 70.79 | 0.5048544 | dsi:Dsimw501\_GD21160 |
| A0B562 | TRINITY\_DN7764\_c0\_g1\_i1 |  496 |  142 |      162 |    139 |     79 |  56.83 |      0 | 73.38 | 0.5563380 | mtp:Mthe\_0033        |

``` r
# Load KEGG to KO as a table
ko <- read.table(koFile, sep="\t", header=FALSE)
# Set column names
colnames(ko) <- c("kegg", "ko")
# Check the ko table
kable(head(ko))
```

| kegg           | ko        |
| :------------- | :-------- |
| mu:66105       | ko:K06689 |
| lsv:3772835    | ko:K02704 |
| nma:NMA0736    | ko:K04043 |
| mmi:MMAR\_4089 | ko:K02111 |
| hsa:8337       | ko:K11251 |
| ath:AT1G04480  | ko:K02894 |

``` r
# Merge KOs
blastKo <- merge(blastKegg, ko)
# Check the blast ko table
kable(head(blastKo))
```

| kegg         | sp     | trans                       | qlen | slen | bitscore | length | nident | pident | evalue |  ppos |       cov | ko        |
| :----------- | :----- | :-------------------------- | ---: | ---: | -------: | -----: | -----: | -----: | -----: | ----: | --------: | :-------- |
| aae:aq\_1065 | O67161 | TRINITY\_DN9495\_c0\_g1\_i2 | 1244 |  342 |      316 |    342 |    172 |  50.29 |      0 | 70.76 | 0.5029240 | ko:K00134 |
| aae:aq\_484  | O66778 | TRINITY\_DN9573\_c0\_g1\_i1 | 1574 |  426 |      425 |    434 |    216 |  49.77 |      0 | 69.35 | 0.5070423 | ko:K01689 |
| aae:aq\_679  | O66907 | TRINITY\_DN9485\_c0\_g1\_i1 | 1847 |  503 |      612 |    502 |    296 |  58.96 |      0 | 77.09 | 0.5884692 | ko:K02111 |
| aae:aq\_996  | O67118 | TRINITY\_DN8020\_c0\_g1\_i1 | 2247 |  632 |      633 |    605 |    329 |  54.38 |      0 | 74.55 | 0.5205696 | ko:K04043 |
| aag:5564556  | Q1HQU2 | TRINITY\_DN9453\_c0\_g1\_i1 | 1110 |  297 |      393 |    272 |    197 |  72.43 |      0 | 86.03 | 0.6632997 | ko:K02932 |
| aag:5564556  | Q1HQU2 | TRINITY\_DN8136\_c0\_g1\_i1 |  999 |  297 |      386 |    294 |    209 |  71.09 |      0 | 82.99 | 0.7037037 | ko:K02932 |

``` r
tx2gene <- unique(subset(blastKo, select=c(trans, ko)))
# Check the tx2gene table
kable(head(tx2gene))
```

| trans                       | ko        |
| :-------------------------- | :-------- |
| TRINITY\_DN9495\_c0\_g1\_i2 | ko:K00134 |
| TRINITY\_DN9573\_c0\_g1\_i1 | ko:K01689 |
| TRINITY\_DN9485\_c0\_g1\_i1 | ko:K02111 |
| TRINITY\_DN8020\_c0\_g1\_i1 | ko:K04043 |
| TRINITY\_DN9453\_c0\_g1\_i1 | ko:K02932 |
| TRINITY\_DN8136\_c0\_g1\_i1 | ko:K02932 |

``` r
# Write as a csv file, excluding row.names
write.csv(tx2gene, file="tx2gene.csv", row.names=FALSE)
```
