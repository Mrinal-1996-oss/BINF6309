Lung Metagenomic Analysis
================

``` r
library(metagenomeSeq)
data(lungData)
```

``` r
library(metagenomeSeq)
#Data Preparation
##Reading in a biom file
library(biomformat)
biom_file <- system.file("extdata", "min_sparse_otu_table.biom",
                          package = "biomformat")
b <- read_biom(biom_file)
biom2MRexperiment(b)
```

    ## MRexperiment (storageMode: environment)
    ## assayData: 5 features, 6 samples 
    ##   element names: counts 
    ## protocolData: none
    ## phenoData: none
    ## featureData: none
    ## experimentData: use 'experimentData(object)'
    ## Annotation:

``` r
### MRexperiment (storageMode: environment)
### assayData: 5 features, 6 samples
### element names: counts
### protocolData: none
### phenoData: none
### featureData: none
### experimentData: use 'experimentData(object)'
##Loading Count Data
dataDirectory <- system.file("extdata", package = "metagenomeSeq")
lung = loadMeta(file.path(dataDirectory, "CHK_NAME.otus.count.csv"))
dim(lung$counts)
```

    ## [1] 1000   78

``` r
### [1] 1000   78
##Loading Taxonomy
taxa = read.delim(file.path(dataDirectory, "CHK_otus.taxonomy.csv"), stringsAsFactors = FALSE)
##Loading Metadata
clin = loadPhenoData(file.path(dataDirectory, "CHK_clinical.csv"), tran = TRUE)
ord = match(colnames(lung$counts), rownames(clin))
clin = clin[ord, ]
head(clin[1:2, ])
```

    ##                                          SampleType          SiteSampled
    ## CHK_6467_E3B11_BRONCH2_PREWASH_V1V2 Bronch2.PreWash Bronchoscope.Channel
    ## CHK_6467_E3B11_OW_V1V2                           OW           OralCavity
    ##                                     SmokingStatus
    ## CHK_6467_E3B11_BRONCH2_PREWASH_V1V2        Smoker
    ## CHK_6467_E3B11_OW_V1V2                     Smoker

``` r
###                                          SampleType
### CHK_6467_E3B11_BRONCH2_PREWASH_V1V2 Bronch2.PreWash
### CHK_6467_E3B11_OW_V1V2                           OW
###                                              SiteSampled
### CHK_6467_E3B11_BRONCH2_PREWASH_V1V2 Bronchoscope.Channel
### CHK_6467_E3B11_OW_V1V2                        OralCavity
###                                     SmokingStatus
### CHK_6467_E3B11_BRONCH2_PREWASH_V1V2        Smoker
### CHK_6467_E3B11_OW_V1V2                     Smoker
##Creating an MRexperiment object
phenotypeData = AnnotatedDataFrame(clin)
phenotypeData
```

    ## An object of class 'AnnotatedDataFrame'
    ##   rowNames: CHK_6467_E3B11_BRONCH2_PREWASH_V1V2 CHK_6467_E3B11_OW_V1V2
    ##     ... CHK_6467_E3B09_BAL_A_V1V2 (78 total)
    ##   varLabels: SampleType SiteSampled SmokingStatus
    ##   varMetadata: labelDescription

``` r
### An object of class 'AnnotatedDataFrame'
###   rowNames: CHK_6467_E3B11_BRONCH2_PREWASH_V1V2
###     CHK_6467_E3B11_OW_V1V2 ...
###     CHK_6467_E3B09_BAL_A_V1V2 (78 total)
###   varLabels: SampleType SiteSampled SmokingStatus
###   varMetadata: labelDescription
OTUdata = AnnotatedDataFrame(taxa)
OTUdata
```

    ## An object of class 'AnnotatedDataFrame'
    ##   rowNames: 1 2 ... 1000 (1000 total)
    ##   varLabels: OTU Taxonomy ... strain (10 total)
    ##   varMetadata: labelDescription

``` r
### An object of class 'AnnotatedDataFrame'
###   rowNames: 1 2 ... 1000 (1000 total)
###   varLabels: OTU Taxonomy ... strain (10 total)
###   varMetadata: labelDescription
obj = newMRexperiment(lung$counts,phenoData=phenotypeData,featureData=OTUdata)
### Links to a paper providing further details can be included optionally.
### experimentData(obj) = annotate::pmid2MIAME("21680950")
obj
```

    ## MRexperiment (storageMode: environment)
    ## assayData: 1000 features, 78 samples 
    ##   element names: counts 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: CHK_6467_E3B11_BRONCH2_PREWASH_V1V2
    ##     CHK_6467_E3B11_OW_V1V2 ... CHK_6467_E3B09_BAL_A_V1V2 (78 total)
    ##   varLabels: SampleType SiteSampled SmokingStatus
    ##   varMetadata: labelDescription
    ## featureData
    ##   featureNames: 1 2 ... 1000 (1000 total)
    ##   fvarLabels: OTU Taxonomy ... strain (10 total)
    ##   fvarMetadata: labelDescription
    ## experimentData: use 'experimentData(object)'
    ## Annotation:

``` r
### MRexperiment (storageMode: environment)
### assayData: 1000 features, 78 samples
###   element names: counts
### protocolData: none
### phenoData
###   sampleNames: CHK_6467_E3B11_BRONCH2_PREWASH_V1V2
###     CHK_6467_E3B11_OW_V1V2 ...
###     CHK_6467_E3B09_BAL_A_V1V2 (78 total)
###   varLabels: SampleType SiteSampled SmokingStatus
###   varMetadata: labelDescription
### featureData
###   featureNames: 1 2 ... 1000 (1000 total)
###   fvarLabels: OTU Taxonomy ... strain (10 total)
###   fvarMetadata: labelDescription
### experimentData: use 'experimentData(object)'
### Annotation:
##Example Datasets
data(lungData)
lungData
```

    ## MRexperiment (storageMode: environment)
    ## assayData: 51891 features, 78 samples 
    ##   element names: counts 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: CHK_6467_E3B11_BRONCH2_PREWASH_V1V2
    ##     CHK_6467_E3B11_OW_V1V2 ... CHK_6467_E3B09_BAL_A_V1V2 (78 total)
    ##   varLabels: SampleType SiteSampled SmokingStatus
    ##   varMetadata: labelDescription
    ## featureData
    ##   featureNames: 1 2 ... 51891 (51891 total)
    ##   fvarLabels: taxa
    ##   fvarMetadata: labelDescription
    ## experimentData: use 'experimentData(object)'
    ## Annotation:

``` r
### MRexperiment (storageMode: environment)
### assayData: 51891 features, 78 samples
###   element names: counts
### protocolData: none
### phenoData
###   sampleNames: CHK_6467_E3B11_BRONCH2_PREWASH_V1V2
###     CHK_6467_E3B11_OW_V1V2 ...
###     CHK_6467_E3B09_BAL_A_V1V2 (78 total)
###   varLabels: SampleType SiteSampled SmokingStatus
###   varMetadata: labelDescription
### featureData
###   featureNames: 1 2 ... 51891 (51891 total)
###   fvarLabels: taxa
###   fvarMetadata: labelDescription
### experimentData: use 'experimentData(object)'
### Annotation:
##Useful Commands
phenoData(obj)
```

    ## An object of class 'AnnotatedDataFrame'
    ##   sampleNames: CHK_6467_E3B11_BRONCH2_PREWASH_V1V2
    ##     CHK_6467_E3B11_OW_V1V2 ... CHK_6467_E3B09_BAL_A_V1V2 (78 total)
    ##   varLabels: SampleType SiteSampled SmokingStatus
    ##   varMetadata: labelDescription

``` r
### An object of class 'AnnotatedDataFrame'
###   sampleNames: CHK_6467_E3B11_BRONCH2_PREWASH_V1V2
###     CHK_6467_E3B11_OW_V1V2 ...
###     CHK_6467_E3B09_BAL_A_V1V2 (78 total)
###   varLabels: SampleType SiteSampled SmokingStatus
###   varMetadata: labelDescription
head(pData(obj), 3)
```

    ##                                          SampleType          SiteSampled
    ## CHK_6467_E3B11_BRONCH2_PREWASH_V1V2 Bronch2.PreWash Bronchoscope.Channel
    ## CHK_6467_E3B11_OW_V1V2                           OW           OralCavity
    ## CHK_6467_E3B08_OW_V1V2                           OW           OralCavity
    ##                                     SmokingStatus
    ## CHK_6467_E3B11_BRONCH2_PREWASH_V1V2        Smoker
    ## CHK_6467_E3B11_OW_V1V2                     Smoker
    ## CHK_6467_E3B08_OW_V1V2                  NonSmoker

``` r
###                                          SampleType
### CHK_6467_E3B11_BRONCH2_PREWASH_V1V2 Bronch2.PreWash
### CHK_6467_E3B11_OW_V1V2                           OW
### CHK_6467_E3B08_OW_V1V2                           OW
###                                              SiteSampled
### CHK_6467_E3B11_BRONCH2_PREWASH_V1V2 Bronchoscope.Channel
### CHK_6467_E3B11_OW_V1V2                        OralCavity
### CHK_6467_E3B08_OW_V1V2                        OralCavity
###                                     SmokingStatus
### CHK_6467_E3B11_BRONCH2_PREWASH_V1V2        Smoker
### CHK_6467_E3B11_OW_V1V2                     Smoker
### CHK_6467_E3B08_OW_V1V2                  NonSmoker
featureData(obj)
```

    ## An object of class 'AnnotatedDataFrame'
    ##   featureNames: 1 2 ... 1000 (1000 total)
    ##   varLabels: OTU Taxonomy ... strain (10 total)
    ##   varMetadata: labelDescription

``` r
### An object of class 'AnnotatedDataFrame'
###   featureNames: 1 2 ... 1000 (1000 total)
###   varLabels: OTU Taxonomy ... strain (10 total)
###   varMetadata: labelDescription
head(fData(obj)[, -c(2, 10)], 3)
```

    ##   OTU superkingdom         phylum                  class             order
    ## 1   1     Bacteria Proteobacteria  Epsilonproteobacteria Campylobacterales
    ## 2   2         <NA>           <NA>                   <NA>              <NA>
    ## 3   3     Bacteria Actinobacteria Actinobacteria (class)   Actinomycetales
    ##               family         genus                  species
    ## 1 Campylobacteraceae Campylobacter     Campylobacter rectus
    ## 2               <NA>          <NA>                     <NA>
    ## 3   Actinomycetaceae   Actinomyces Actinomyces radicidentis

``` r
###   OTU superkingdom         phylum                  class
### 1   1     Bacteria Proteobacteria  Epsilonproteobacteria
### 2   2         <NA>           <NA>                   <NA>
### 3   3     Bacteria Actinobacteria Actinobacteria (class)
###               order             family         genus
### 1 Campylobacterales Campylobacteraceae Campylobacter
### 2              <NA>               <NA>          <NA>
### 3   Actinomycetales   Actinomycetaceae   Actinomyces
###                    species
### 1     Campylobacter rectus
### 2                     <NA>
### 3 Actinomyces radicidentis
head(MRcounts(obj[, 1:2]))
```

    ##   CHK_6467_E3B11_BRONCH2_PREWASH_V1V2 CHK_6467_E3B11_OW_V1V2
    ## 1                                   0                      0
    ## 2                                   0                      0
    ## 3                                   0                      0
    ## 4                                   0                      0
    ## 5                                   0                      0
    ## 6                                   0                      0

``` r
###   CHK_6467_E3B11_BRONCH2_PREWASH_V1V2
### 1                                   0
### 2                                   0
### 3                                   0
### 4                                   0
### 5                                   0
### 6                                   0
###   CHK_6467_E3B11_OW_V1V2
### 1                      0
### 2                      0
### 3                      0
### 4                      0
### 5                      0
### 6                      0
featuresToKeep = which(rowSums(obj) >= 100)
samplesToKeep = which(pData(obj)$SmokingStatus == "Smoker")
obj_smokers = obj[featuresToKeep, samplesToKeep]
### MRexperiment (storageMode: environment)
### assayData: 1 features, 33 samples
###   element names: counts
### protocolData: none
### phenoData
###   sampleNames: CHK_6467_E3B11_BRONCH2_PREWASH_V1V2
###     CHK_6467_E3B11_OW_V1V2 ...
###     CHK_6467_E3B09_BAL_A_V1V2 (33 total)
###   varLabels: SampleType SiteSampled SmokingStatus
###   varMetadata: labelDescription
### featureData##   featureNames: 570
###   fvarLabels: OTU Taxonomy ... strain (10 total)
###   fvarMetadata: labelDescription
### experimentData: use 'experimentData(object)'
### Annotation:
head(pData(obj_smokers), 3)
```

    ##                                          SampleType          SiteSampled
    ## CHK_6467_E3B11_BRONCH2_PREWASH_V1V2 Bronch2.PreWash Bronchoscope.Channel
    ## CHK_6467_E3B11_OW_V1V2                           OW           OralCavity
    ## CHK_6467_E3B11_BAL_A_V1V2                     BAL.A                 Lung
    ##                                     SmokingStatus
    ## CHK_6467_E3B11_BRONCH2_PREWASH_V1V2        Smoker
    ## CHK_6467_E3B11_OW_V1V2                     Smoker
    ## CHK_6467_E3B11_BAL_A_V1V2                  Smoker

``` r
###                                          SampleType
### CHK_6467_E3B11_BRONCH2_PREWASH_V1V2 Bronch2.PreWash
### CHK_6467_E3B11_OW_V1V2                           OW
### CHK_6467_E3B11_BAL_A_V1V2                     BAL.A
###                                              SiteSampled
### CHK_6467_E3B11_BRONCH2_PREWASH_V1V2 Bronchoscope.Channel
### CHK_6467_E3B11_OW_V1V2                        OralCavity
### CHK_6467_E3B11_BAL_A_V1V2                           Lung
###                                     SmokingStatus
### CHK_6467_E3B11_BRONCH2_PREWASH_V1V2        Smoker
### CHK_6467_E3B11_OW_V1V2                     Smoker
### CHK_6467_E3B11_BAL_A_V1V2                  Smoker
head(normFactors(obj))
```

    ##                                     [,1]
    ## CHK_6467_E3B11_BRONCH2_PREWASH_V1V2   NA
    ## CHK_6467_E3B11_OW_V1V2                NA
    ## CHK_6467_E3B08_OW_V1V2                NA
    ## CHK_6467_E3B07_BAL_A_V1V2             NA
    ## CHK_6467_E3B11_BAL_A_V1V2             NA
    ## CHK_6467_E3B09_OP_V1V2                NA

``` r
###                                     [,1]
### CHK_6467_E3B11_BRONCH2_PREWASH_V1V2   NA
### CHK_6467_E3B11_OW_V1V2                NA
### CHK_6467_E3B08_OW_V1V2                NA
### CHK_6467_E3B07_BAL_A_V1V2             NA
### CHK_6467_E3B11_BAL_A_V1V2             NA
### CHK_6467_E3B09_OP_V1V2                NA
normFactors(obj) <- rnorm(ncol(obj))
head(normFactors(obj))
```

    ## CHK_6467_E3B11_BRONCH2_PREWASH_V1V2              CHK_6467_E3B11_OW_V1V2 
    ##                         -2.02640688                         -1.19028031 
    ##              CHK_6467_E3B08_OW_V1V2           CHK_6467_E3B07_BAL_A_V1V2 
    ##                          0.07345299                          1.32868132 
    ##           CHK_6467_E3B11_BAL_A_V1V2              CHK_6467_E3B09_OP_V1V2 
    ##                          0.64671372                          1.26433095

``` r
### CHK_6467_E3B11_BRONCH2_PREWASH_V1V2
###                           1.3709584
###              CHK_6467_E3B11_OW_V1V2
###                          -0.5646982
###              CHK_6467_E3B08_OW_V1V2
###                           0.3631284
###           CHK_6467_E3B07_BAL_A_V1V2
###                           0.6328626
###           CHK_6467_E3B11_BAL_A_V1V2
###                           0.4042683
###              CHK_6467_E3B09_OP_V1V2
###                          -0.1061245
head(libSize(obj))
```

    ##                                     [,1]
    ## CHK_6467_E3B11_BRONCH2_PREWASH_V1V2    0
    ## CHK_6467_E3B11_OW_V1V2                16
    ## CHK_6467_E3B08_OW_V1V2                 1
    ## CHK_6467_E3B07_BAL_A_V1V2              2
    ## CHK_6467_E3B11_BAL_A_V1V2            118
    ## CHK_6467_E3B09_OP_V1V2                 5

``` r
###                                     [,1]
### CHK_6467_E3B11_BRONCH2_PREWASH_V1V2    0
### CHK_6467_E3B11_OW_V1V2                16
### CHK_6467_E3B08_OW_V1V2                 1
### CHK_6467_E3B07_BAL_A_V1V2              2
### CHK_6467_E3B11_BAL_A_V1V2            118
### CHK_6467_E3B09_OP_V1V2                 5
libSize(obj) <- rnorm(ncol(obj))
head(libSize(obj))
```

    ## CHK_6467_E3B11_BRONCH2_PREWASH_V1V2              CHK_6467_E3B11_OW_V1V2 
    ##                          1.26303819                          0.01764207 
    ##              CHK_6467_E3B08_OW_V1V2           CHK_6467_E3B07_BAL_A_V1V2 
    ##                         -1.19778706                          0.82826731 
    ##           CHK_6467_E3B11_BAL_A_V1V2              CHK_6467_E3B09_OP_V1V2 
    ##                         -0.89837499                         -1.16201418

``` r
### CHK_6467_E3B11_BRONCH2_PREWASH_V1V2
###                         -0.8857763011
###              CHK_6467_E3B11_OW_V1V2
###                         -1.09978090
###              CHK_6467_E3B08_OW_V1V2
###                          1.51270701
###           CHK_6467_E3B07_BAL_A_V1V2
###                          0.25792144
###           CHK_6467_E3B11_BAL_A_V1V2
###                          0.08844023
###              CHK_6467_E3B09_OP_V1V2
###                         -0.12089654
#Normalization
data(lungData)
p = cumNormStatFast(lungData)
lungData = cumNorm(lungData, p = p)
##Exporting Data
mat = MRcounts(lungData, norm = TRUE, log = TRUE)[1:5, 1:5]
exportMat(mat, file = file.path(dataDirectory, "tmp.tsv"))
###Default value being used
head(read.csv(file = file.path(dataDirectory, "tmp.tsv"), sep = "\t"))
```

    ##   Taxa.and.Samples CHK_6467_E3B11_BRONCH2_PREWASH_V1V2 CHK_6467_E3B11_OW_V1V2
    ## 1                1                                   0                      0
    ## 2                2                                   0                      0
    ## 3                3                                   0                      0
    ## 4                4                                   0                      0
    ## 5                5                                   0                      0
    ##   CHK_6467_E3B08_OW_V1V2 CHK_6467_E3B07_BAL_A_V1V2 CHK_6467_E3B11_BAL_A_V1V2
    ## 1                      0                         0                         0
    ## 2                      0                         0                         0
    ## 3                      0                         0                         0
    ## 4                      0                         0                         0
    ## 5                      0                         0                         0

``` r
###                               Subject Scaling.factor
### 1 CHK_6467_E3B11_BRONCH2_PREWASH_V1V2             67
### 2              CHK_6467_E3B11_OW_V1V2           2475
### 3              CHK_6467_E3B08_OW_V1V2           2198
### 4           CHK_6467_E3B07_BAL_A_V1V2            836
### 5           CHK_6467_E3B11_BAL_A_V1V2           1008
###   Quantile.value Number.of.identified.features Library.size
### 1              2                            60          271
### 2              1                          3299         7863
### 3              1                          2994         8360
### 4              1                          1188         5249
### 5              1                          1098         3383
#Statistical Testing
##  Example using fitFeatureModel for differential abundance testing
data(lungData)
lungData = lungData[, -which(is.na(pData(lungData)$SmokingStatus))] 
lungData = filterData(lungData, present = 30, depth = 1)
lungData <- cumNorm(lungData, p = 0.5)
pd <- pData(lungData)
mod <- model.matrix(~1 + SmokingStatus, data = pd)
lungres1 = fitFeatureModel(lungData, mod)
head(MRcoefs(lungres1))
```

    ##           logFC        se      pvalues   adjPvalues
    ## 3465  -4.824949 0.5697511 0.000000e+00 0.000000e+00
    ## 35827 -4.304266 0.5445548 2.664535e-15 1.079137e-13
    ## 2817   2.320656 0.4324661 8.045793e-08 1.629273e-06
    ## 2735   2.260203 0.4331098 1.803341e-07 2.921412e-06
    ## 5411   1.748296 0.3092461 1.572921e-08 4.246888e-07
    ## 48745 -1.645805 0.3293117 5.801451e-07 7.831959e-06

``` r
###           logFC        se      pvalues   adjPvalues
### 3465  -4.824949 0.5697511 0.000000e+00 0.000000e+00
### 35827 -4.304266 0.5445548 2.664535e-15 1.079137e-13
### 2817   2.320656 0.4324661 8.045793e-08 1.629273e-06
### 2735   2.260203 0.4331098 1.803341e-07 2.921412e-06
### 5411   1.748296 0.3092461 1.572921e-08 4.246888e-07
### 48745 -1.645805 0.3293117 5.801451e-07 7.831959e-06
##Example using fitZig for differential abundance testing
data(lungData)
controls = grep("Extraction.Control", pData(lungData)$SampleType)
lungTrim = lungData[, -controls]
rareFeatures = which(rowSums(MRcounts(lungTrim) > 0) < 10)
lungTrim = lungTrim[-rareFeatures, ]
lungp = cumNormStat(lungTrim, pFlag = TRUE, main = "Trimmed lung data")
```

    ## Default value being used.

![](lungData_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
###Default value being used
lungTrim = cumNorm(lungTrim, p = lungp)

smokingStatus = pData(lungTrim)$SmokingStatus
bodySite = pData(lungTrim)$SampleType
normFactor = normFactors(lungTrim)
normFactor = log2(normFactor/median(normFactor) + 1)
mod = model.matrix(~smokingStatus + bodySite + normFactor)
settings = zigControl(maxit = 10, verbose = TRUE)
fit = fitZig(obj = lungTrim, mod = mod, useCSSoffset = FALSE, 
             control = settings)
```

    ## it= 0, nll=88.42, log10(eps+1)=Inf, stillActive=1029
    ## it= 1, nll=93.56, log10(eps+1)=0.06, stillActive=261
    ## it= 2, nll=93.46, log10(eps+1)=0.05, stillActive=120
    ## it= 3, nll=93.80, log10(eps+1)=0.05, stillActive=22
    ## it= 4, nll=93.94, log10(eps+1)=0.03, stillActive=3
    ## it= 5, nll=93.93, log10(eps+1)=0.00, stillActive=1
    ## it= 6, nll=93.90, log10(eps+1)=0.00, stillActive=1
    ## it= 7, nll=93.87, log10(eps+1)=0.00, stillActive=1
    ## it= 8, nll=93.86, log10(eps+1)=0.00, stillActive=1
    ## it= 9, nll=93.85, log10(eps+1)=0.00, stillActive=1

``` r
### it= 0, nll=88.42, log10(eps+1)=Inf, stillActive=1029
### it= 1, nll=93.56, log10(eps+1)=0.06, stillActive=261
### it= 2, nll=93.46, log10(eps+1)=0.05, stillActive=120
### it= 3, nll=93.80, log10(eps+1)=0.05, stillActive=22
### it= 4, nll=93.94, log10(eps+1)=0.03, stillActive=3
### it= 5, nll=93.93, log10(eps+1)=0.00, stillActive=1
### it= 6, nll=93.90, log10(eps+1)=0.00, stillActive=1
### it= 7, nll=93.87, log10(eps+1)=0.00, stillActive=1
### it= 8, nll=93.86, log10(eps+1)=0.00, stillActive=1
### it= 9, nll=93.85, log10(eps+1)=0.00, stillActive=1
### The default, useCSSoffset = TRUE, automatically includes
### the CSS scaling normalization factor.
##Multiple Groups
###maxit is for demonstration purposes
# maxit=1 is for demonstration purposes
settings = zigControl(maxit = 1, verbose = FALSE) 
mod = model.matrix(~bodySite)
colnames(mod) = levels(bodySite)
# fitting the ZIG model
res = fitZig(obj = lungTrim, mod = mod, control = settings) 
# The output of fitZig contains a list of various useful
# items. hint: names(res). Probably the most useful is the 
# limma 'MLArrayLM' object called fit.
zigFit = slot(res, "fit")
finalMod = slot(res, "fit")$design
contrast.matrix = makeContrasts(BAL.A - BAL.B, OW - PSB, 
                                levels = finalMod) 
fit2 = contrasts.fit(zigFit, contrast.matrix)
fit2 = eBayes(fit2)
topTable(fit2)
```

    ##       BAL.A...BAL.B  OW...PSB   AveExpr         F      P.Value  adj.P.Val
    ## 18531    0.37318792  2.075648 0.7343081 12.715105 5.359780e-05 0.02813711
    ## 6291    -0.10695735  1.658829 0.4671470 12.956898 5.482439e-05 0.02813711
    ## 37977   -0.37995461  2.174071 0.4526060 12.528733 8.203239e-05 0.02813711
    ## 6901     0.17344138  1.466113 0.2435881 12.018652 1.335806e-04 0.03212047
    ## 40291    0.06892926  1.700238 0.2195735 11.803380 1.560761e-04 0.03212047
    ## 36117   -0.28665883  2.233996 0.4084024 10.571931 3.012092e-04 0.05013569
    ## 7343    -0.22859078  1.559465 0.3116465 10.090602 3.931844e-04 0.05013569
    ## 7342     0.59882970  1.902346 0.5334647  9.410984 4.901651e-04 0.05013569
    ## 1727     1.09837459 -2.160466 0.7780167  9.346013 5.027597e-04 0.05013569
    ## 40329   -0.07145998  1.481582 0.2475735  9.700136 5.259032e-04 0.05013569

``` r
###       BAL.A...BAL.B  OW...PSB   AveExpr         F
### 18531    0.37318792  2.075648 0.7343081 12.715105
### 6291    -0.10695735  1.658829 0.4671470 12.956898
### 37977   -0.37995461  2.174071 0.4526060 12.528733
### 6901     0.17344138  1.466113 0.2435881 12.018652
### 40291    0.06892926  1.700238 0.2195735 11.803380
### 36117   -0.28665883  2.233996 0.4084024 10.571931
### 7343    -0.22859078  1.559465 0.3116465 10.090602
### 7342     0.59882970  1.902346 0.5334647  9.410984
### 1727     1.09837459 -2.160466 0.7780167  9.346013
### 40329   -0.07145998  1.481582 0.2475735  9.700136
###            P.Value  adj.P.Val
### 18531 5.359780e-05 0.02813711
### 6291  5.482439e-05 0.02813711
### 37977 8.203239e-05 0.02813711
### 6901  1.335806e-04 0.03212047
### 40291 1.560761e-04 0.03212047
### 36117 3.012092e-04 0.05013569
### 7343  3.931844e-04 0.05013569
### 7342  4.901651e-04 0.05013569
### 1727  5.027597e-04 0.05013569
### 40329 5.259032e-04 0.05013569
##Exporting Fits
taxa = sapply(strsplit(as.character(fData(lungTrim)$taxa), split = ";"),  function(i) {
    i[length(i)]
  })
head(MRcoefs(fit, taxa = taxa, coef = 2))
```

    ##                                   smokingStatusSmoker      pvalues   adjPvalues
    ## Neisseria polysaccharea                     -4.031612 3.927097e-11 2.959194e-08
    ## Neisseria meningitidis                      -3.958899 5.751592e-11 2.959194e-08
    ## Prevotella intermedia                       -2.927686 4.339587e-09 8.930871e-07
    ## Porphyromonas sp. UQD 414                   -2.675306 1.788697e-07 1.357269e-05
    ## Prevotella paludivivens                      2.575672 1.360718e-07 1.272890e-05
    ## Leptotrichia sp. oral clone FP036            2.574172 3.544957e-04 1.414122e-03

``` r
###                                   smokingStatusSmoker
### Neisseria polysaccharea                     -4.031612
### Neisseria meningitidis                      -3.958899
### Prevotella intermedia                       -2.927686
### Porphyromonas sp. UQD 414                   -2.675306
### Prevotella paludivivens                      2.575672
### Leptotrichia sp. oral clone FP036            2.574172
###                                        pvalues   adjPvalues
### Neisseria polysaccharea           3.927097e-11 2.959194e-08
### Neisseria meningitidis            5.751592e-11 2.959194e-08
### Prevotella intermedia             4.339587e-09 8.930871e-07
### Porphyromonas sp. UQD 414         1.788697e-07 1.357269e-05
### Prevotella paludivivens           1.360718e-07 1.272890e-05
### Leptotrichia sp. oral clone FP036 3.544957e-04 1.414122e-03
###Log Normal Permutation Test
coeffOfInterest = 2
res = fitLogNormal(obj = lungTrim, mod = mod, useCSSoffset = FALSE,
                   B = 10, coef = coeffOfInterest)
### extract p.values and adjust for multiple testing res$p are
### the p-values calculated through permutation
adjustedPvalues = p.adjust(res$p, method = "fdr")
### extract the absolute fold-change estimates
foldChange = abs(res$fit$coef[, coeffOfInterest])
### determine features still significant and order by the
sigList = which(adjustedPvalues <= 0.05)
sigList = sigList[order(foldChange[sigList])]
### view the top taxa associated with the coefficient of interest.
head(taxa[sigList])
```

    ## [1] "Burkholderia sp. LSB13479"    "Actinomyces sp. P31/96/2"    
    ## [3] "Veillonella montpellierensis" "Veillonella montpellierensis"
    ## [5] "Campylobacter curvus"         "Anaeroglobus geminatus"

``` r
### [1] "Veillonella montpellierensis"
### [2] "Veillonella sp. oral clone VeillI7"
### [3] "Listeria grayi"
### [4] "Megasphaera micronuciformis"
### [5] "Prevotella intermedia"
### [6] "Campylobacter curvus"
#Feature Specific
head(MRtable(fit, coef = 2, taxa = 1:length(fData(lungTrim)$taxa)))
```

    ##     +samples in group 0 +samples in group 1 counts in group 0 counts in group 1
    ## 63                   24                   6              1538                11
    ## 779                  23                   7              1512                22
    ## 358                  24                   1               390                 1
    ## 499                  21                   2               326                 2
    ## 25                   15                  26               162              1893
    ## 928                   2                  11                 4                91
    ##     smokingStatusSmoker      pvalues   adjPvalues
    ## 63            -4.031612 3.927097e-11 2.959194e-08
    ## 779           -3.958899 5.751592e-11 2.959194e-08
    ## 358           -2.927686 4.339587e-09 8.930871e-07
    ## 499           -2.675306 1.788697e-07 1.357269e-05
    ## 25             2.575672 1.360718e-07 1.272890e-05
    ## 928            2.574172 3.544957e-04 1.414122e-03

``` r
patients = sapply(strsplit(rownames(pData(lungTrim)), split = "_"), 
  function(i) {
    i[3] 
  })
pData(lungTrim)$patients = patients
classIndex = list(smoker = which(pData(lungTrim)$SmokingStatus ==
"Smoker"))
classIndex$nonsmoker = which(pData(lungTrim)$SmokingStatus ==
"NonSmoker") 
otu = 779
# plotOTU
plotOTU(lungTrim, otu = otu, classIndex, main = "Neisseria meningitidis")
```

![](lungData_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
# Now multiple OTUs annotated similarly
x = fData(lungTrim)$taxa[otu]
otulist = grep(x, fData(lungTrim)$taxa)
# plotGenus
plotGenus(lungTrim, otulist, classIndex, labs = FALSE, main = "Neisseria meningit") 
lablist <- c("S", "NS")
axis(1, at = seq(1, 6, by = 1), labels = rep(lablist, times = 3))
```

![](lungData_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->

``` r
#Citing metagenomeSeq
citation("metagenomeSeq")
```

    ## 
    ## Please cite the top for the original statistical method and
    ## normalization method implemented in metagenomeSeq and the bottom for
    ## the software/vignette guide. Time series analysis/function is described
    ## in the third citation.
    ## 
    ##   JN Paulson, OC Stine, HC Bravo, M Pop.  Differential abundance
    ##   analysis for microbial marker-gene surveys. Nat Meth Accepted
    ## 
    ##   JN Paulson, H Talukder, M Pop, HC Bravo. metagenomeSeq: Statistical
    ##   analysis for sparse high-throughput sequencing. Bioconductor package:
    ##   1.26.3. http://cbcb.umd.edu/software/metagenomeSeq
    ## 
    ##   H Talukder*, JN Paulson*, HC Bravo. Longitudinal differential
    ##   abundance analysis of marker-gene surveys. Submitted
    ## 
    ## To see these entries in BibTeX format, use 'print(<citation>,
    ## bibtex=TRUE)', 'toBibtex(.)', or set
    ## 'options(citation.bibtex.max=999)'.

``` r
### To cite the original statistical method and
### normalization method implemented in metagenomeSeq use
### Paulson JN, Stine OC, Bravo HC, Pop M (2013).
### "Differential abundance analysis for microbial
### marker-gene surveys." _Nat Meth_,*advance online
### publication*. doi: 10.1038/nmeth.2658 (URL:
### https://doi.org/10.1038/nmeth.2658), <URL:
### http://www.nature.com/nmeth/journal/vaop/ncurrent/abs/nmeth.2658.html>.
### To cite the metagenomeSeq software/vignette guide use
### Paulson JN, Olson ND, Braccia DJ, Wagner J, Talukder
### H, Pop M, Bravo HC (2013). _metagenomeSeq:
### Statistical analysis for sparse high-throughput
### sequncing._. Bioconductor package, <URL:
### http://www.cbcb.umd.edu/software/metagenomeSeq>.
### To cite time series analysis/function fitTimeSeries
### use
### Paulson*JN, Talukder*H, Bravo HC (2017).
### "Longitudinal differential abundance analysis of
### marker-gene surveys using smoothing splines."
### _biorxiv_. doi: 10.1101/099457 (URL:
### https://doi.org/10.1101/099457), <URL:
### https://www.biorxiv.org/content/10.1101/099457v1>.
### To see these entries in BibTeX format, use
### 'print(<citation>, bibtex=TRUE)', 'toBibtex(.)', or
### set 'options(citation.bibtex.max=999)'.
sessionInfo()
```

    ## R version 3.6.0 (2019-04-26)
    ## Platform: x86_64-redhat-linux-gnu (64-bit)
    ## Running under: CentOS Linux 7 (Core)
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ## [1] biomformat_1.12.0    metagenomeSeq_1.26.3 RColorBrewer_1.1-2  
    ## [4] glmnet_3.0-2         Matrix_1.2-17        limma_3.40.6        
    ## [7] Biobase_2.44.0       BiocGenerics_0.30.0 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.4         plyr_1.8.6         compiler_3.6.0     bitops_1.0-6      
    ##  [5] iterators_1.0.12   tools_3.6.0        digest_0.6.25      rhdf5_2.28.1      
    ##  [9] jsonlite_1.6.1     evaluate_0.14      lattice_0.20-38    rlang_0.4.5       
    ## [13] IHW_1.12.0         foreach_1.5.0      lpsymphony_1.12.0  yaml_2.2.1        
    ## [17] xfun_0.12          stringr_1.4.0      knitr_1.28         gtools_3.8.2      
    ## [21] caTools_1.18.0     locfit_1.5-9.4     grid_3.6.0         Wrench_1.2.0      
    ## [25] fdrtool_1.2.15     rmarkdown_2.1      gdata_2.18.0       Rhdf5lib_1.6.3    
    ## [29] magrittr_1.5       gplots_3.0.3       codetools_0.2-16   htmltools_0.4.0   
    ## [33] matrixStats_0.56.0 shape_1.4.4        KernSmooth_2.23-15 stringi_1.4.6     
    ## [37] slam_0.1-47

``` r
### R version 3.6.2 (2019-12-12)
### Platform: x86_64-pc-linux-gnu (64-bit)
### Running under: Ubuntu 18.04.3 LTS
### Matrix products: default
### BLAS:   /home/biocbuild/bbs-3.10-bioc/R/lib/libRblas.so
### LAPACK: /home/biocbuild/bbs-3.10-bioc/R/lib/libRlapack.so
### locale:
###  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
###  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=C
###  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
###  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
###  [9] LC_ADDRESS=C               LC_TELEPHONE=C
### [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
### attached base packages:
### [1] parallel  stats     graphics  grDevices utils
### [6] datasets  methods   base
### other attached packages:
###  [1] biomformat_1.14.0    gss_2.1-1
###  [3] metagenomeSeq_1.28.2 RColorBrewer_1.1-2
###  [5] glmnet_3.0-2         Matrix_1.2-18
###  [7] limma_3.42.2         Biobase_2.46.0
###  [9] BiocGenerics_0.32.0  knitr_1.27.2
### loaded via a namespace (and not attached):
###  [1] Rcpp_1.0.3         compiler_3.6.2
###  [3] formatR_1.7        highr_0.8
###  [5] plyr_1.8.5         bitops_1.0-6
###  [7] iterators_1.0.12   tools_3.6.2
###  [9] rhdf5_2.30.1       jsonlite_1.6.1
### [11] evaluate_0.14      lattice_0.20-38
### [13] IHW_1.14.0         foreach_1.4.7
### [15] lpsymphony_1.14.0  xfun_0.12
### [17] stringr_1.4.0      gtools_3.8.1
### [19] caTools_1.18.0     locfit_1.5-9.1
### [21] grid_3.6.2         Wrench_1.4.0
### [23] fdrtool_1.2.15     gdata_2.18.0
### [25] Rhdf5lib_1.8.0     magrittr_1.5
### [27] gplots_3.0.1.2     codetools_0.2-16
### [29] matrixStats_0.55.0 shape_1.4.4
### [31] KernSmooth_2.23-16 stringi_1.4.5
### [33] slam_0.1-47
```
