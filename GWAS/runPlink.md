``` r
library(knitr)
```

``` bash
#Get basic summary about the file 
plink --file hapmap1
#Make a binary PED file
plink --file hapmap1 --make-bed --out hapmap1
plink --file hapmap1 --make-bed --mind 0.05 --out highgeno
plink --bfile hapmap1
#Generate a summary of the statistics
plink --bfile hapmap1 --missing --out miss_stat
```

    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to plink.log.
    ## Options in effect:
    ##   --file hapmap1
    ## 
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## Scanning .ped file... 0%1%2%3%4%5%6%7%8%10%11%12%13%14%15%16%17%19%20%21%22%23%24%25%26%28%29%30%31%32%33%34%35%37%38%39%40%41%42%43%44%46%47%48%49%50%51%52%53%55%56%57%58%59%60%61%62%64%65%66%67%68%69%70%71%73%74%75%76%77%78%79%80%82%83%84%85%86%87%88%89%91%92%93%94%95%96%97%98%100%.ped scan complete (for binary autoconversion).
    ## Performing single-pass .bed write (83534 variants, 89 people).
    ## 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%--file: plink.bed + plink.bim + plink.fam written.
    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to hapmap1.log.
    ## Options in effect:
    ##   --file hapmap1
    ##   --make-bed
    ##   --out hapmap1
    ## 
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## Scanning .ped file... 0%1%2%3%4%5%6%7%8%10%11%12%13%14%15%16%17%19%20%21%22%23%24%25%26%28%29%30%31%32%33%34%35%37%38%39%40%41%42%43%44%46%47%48%49%50%51%52%53%55%56%57%58%59%60%61%62%64%65%66%67%68%69%70%71%73%74%75%76%77%78%79%80%82%83%84%85%86%87%88%89%91%92%93%94%95%96%97%98%100%.ped scan complete (for binary autoconversion).
    ## Performing single-pass .bed write (83534 variants, 89 people).
    ## 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%--file: hapmap1-temporary.bed + hapmap1-temporary.bim + hapmap1-temporary.fam
    ## written.
    ## 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.99441.
    ## 83534 variants and 89 people pass filters and QC.
    ## Among remaining phenotypes, 44 are cases and 45 are controls.
    ## --make-bed to hapmap1.bed + hapmap1.bim + hapmap1.fam ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.
    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to highgeno.log.
    ## Options in effect:
    ##   --file hapmap1
    ##   --make-bed
    ##   --mind 0.05
    ##   --out highgeno
    ## 
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## Scanning .ped file... 0%1%2%3%4%5%6%7%8%10%11%12%13%14%15%16%17%19%20%21%22%23%24%25%26%28%29%30%31%32%33%34%35%37%38%39%40%41%42%43%44%46%47%48%49%50%51%52%53%55%56%57%58%59%60%61%62%64%65%66%67%68%69%70%71%73%74%75%76%77%78%79%80%82%83%84%85%86%87%88%89%91%92%93%94%95%96%97%98%100%.ped scan complete (for binary autoconversion).
    ## Performing single-pass .bed write (83534 variants, 89 people).
    ## 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%--file: highgeno-temporary.bed + highgeno-temporary.bim +
    ## highgeno-temporary.fam written.
    ## 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values loaded from .fam.
    ## 0 people removed due to missing genotype data (--mind).
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.99441.
    ## 83534 variants and 89 people pass filters and QC.
    ## Among remaining phenotypes, 44 are cases and 45 are controls.
    ## --make-bed to highgeno.bed + highgeno.bim + highgeno.fam ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.
    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to plink.log.
    ## Options in effect:
    ##   --bfile hapmap1
    ## 
    ## Warning: No output requested.  Exiting.
    ## 
    ##   plink [input flag(s)...] {command flag(s)...} {other flag(s)...}
    ##   plink --help {flag name(s)...}
    ## 
    ## Commands include --make-bed, --recode, --flip-scan, --merge-list,
    ## --write-snplist, --list-duplicate-vars, --freqx, --missing, --test-mishap,
    ## --hardy, --mendel, --ibc, --impute-sex, --indep-pairphase, --r2, --show-tags,
    ## --blocks, --distance, --genome, --homozyg, --make-rel, --make-grm-gz,
    ## --rel-cutoff, --cluster, --pca, --neighbour, --ibs-test, --regress-distance,
    ## --model, --bd, --gxe, --logistic, --dosage, --lasso, --test-missing,
    ## --make-perm-pheno, --tdt, --qfam, --annotate, --clump, --gene-report,
    ## --meta-analysis, --epistasis, --fast-epistasis, and --score.
    ## 
    ## 'plink --help | more' describes all functions (warning: long).
    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to miss_stat.log.
    ## Options in effect:
    ##   --bfile hapmap1
    ##   --missing
    ##   --out miss_stat
    ## 
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.99441.
    ## --missing: Sample missing data report written to miss_stat.imiss, and
    ## variant-based missing data report written to miss_stat.lmiss.

``` r
file1 <- read.table("miss_stat.lmiss", header = TRUE, sep="\t")
kable(head(file1))
```

| CHR………SNP…N\_MISS…N\_GENO…F\_MISS |
| :-------------------------------- |
| 1 rs6681049 0 89 0                |
| 1 rs4074137 0 89 0                |
| 1 rs7540009 0 89 0                |
| 1 rs1891905 0 89 0                |
| 1 rs9729550 0 89 0                |
| 1 rs3813196 0 89 0                |

``` r
file2 <- read.table("miss_stat.imiss", header  = TRUE , sep="\t")
kable(head(file2))
```

| FID..IID.MISS\_PHENO…N\_MISS…N\_GENO…F\_MISS |
| :------------------------------------------- |
| HCB181 1 N 671 83534 0.008033                |
| HCB182 1 N 1156 83534 0.01384                |
| HCB183 1 N 498 83534 0.005962                |
| HCB184 1 N 412 83534 0.004932                |
| HCB185 1 N 329 83534 0.003939                |
| HCB186 1 N 1233 83534 0.01476                |

``` bash
#Perform statistical analysis of chromosomes
plink --bfile hapmap1 --chr 1 --out res1 --missing
plink --bfile hapmap1 --chr 2 --out res2 --missing
#Make a statistics summary of allele frequency
plink --bfile hapmap1 --freq --out freq_stat
plink --bfile hapmap1 --freq --within pop.phe --out freq_stat
```

    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to res1.log.
    ## Options in effect:
    ##   --bfile hapmap1
    ##   --chr 1
    ##   --missing
    ##   --out res1
    ## 
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## 6762 out of 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.995545.
    ## --missing: Sample missing data report written to res1.imiss, and variant-based
    ## missing data report written to res1.lmiss.
    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to res2.log.
    ## Options in effect:
    ##   --bfile hapmap1
    ##   --chr 2
    ##   --missing
    ##   --out res2
    ## 
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## 7667 out of 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.99694.
    ## --missing: Sample missing data report written to res2.imiss, and variant-based
    ## missing data report written to res2.lmiss.
    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to freq_stat.log.
    ## Options in effect:
    ##   --bfile hapmap1
    ##   --freq
    ##   --out freq_stat
    ## 
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.99441.
    ## --freq: Allele frequencies (founders only) written to freq_stat.frq .
    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to freq_stat.log.
    ## Options in effect:
    ##   --bfile hapmap1
    ##   --freq
    ##   --out freq_stat
    ##   --within pop.phe
    ## 
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values loaded from .fam.
    ## --within: 2 clusters loaded, covering a total of 89 people.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.99441.
    ## --freq: Cluster-stratified allele frequencies (founders only) written to
    ## freq_stat.frq.strat .

``` r
file3 <- read.table("freq_stat.frq.strat", header = TRUE, sep="\t")
kable(head(file3))
```

| CHR………SNP…..CLST…A1…A2……MAF….MAC..NCHROBS |
| :---------------------------------------- |
| 1 rs6681049 1 1 2 0.2333 21 90            |
| 1 rs6681049 2 1 2 0.1932 17 88            |
| 1 rs4074137 1 1 2 0.1 9 90                |
| 1 rs4074137 2 1 2 0.05682 5 88            |
| 1 rs7540009 1 0 2 0 0 90                  |
| 1 rs7540009 2 0 2 0 0 88                  |

``` bash
#For single specific SNP
plink --bfile hapmap1 --snp rs1891905 --freq --within pop.phe --out snp1_frq_stat
#Perform basic association analysis for all SNPis
plink --bfile hapmap1 --assoc --out as1
```

    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to snp1_frq_stat.log.
    ## Options in effect:
    ##   --bfile hapmap1
    ##   --freq
    ##   --out snp1_frq_stat
    ##   --snp rs1891905
    ##   --within pop.phe
    ## 
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## 1 out of 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values loaded from .fam.
    ## --within: 2 clusters loaded, covering a total of 89 people.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## --freq: Cluster-stratified allele frequencies (founders only) written to
    ## snp1_frq_stat.frq.strat .
    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to as1.log.
    ## Options in effect:
    ##   --assoc
    ##   --bfile hapmap1
    ##   --out as1
    ## 
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.99441.
    ## 83534 variants and 89 people pass filters and QC.
    ## Among remaining phenotypes, 44 are cases and 45 are controls.
    ## Writing C/C --assoc report to as1.assoc ... 0%1%2%3%4%6%7%8%9%10%11%12%14%15%16%17%18%19%20%22%23%25%26%27%28%29%31%32%33%34%35%36%38%39%40%41%42%43%45%46%47%48%49%50%51%53%54%55%56%58%59%60%61%62%63%65%66%67%68%69%70%72%73%74%76%77%78%79%80%81%83%84%85%87%88%89%90%92%93%94%95%96%97%98%99%done.

``` r
file4 <- read.table("as1.assoc", header = TRUE, sep="\t")
kable(head(file4))
```

| CHR………SNP………BP…A1……F\_A……F\_U…A2……..CHISQ…………P………..OR   |
| :------------------------------------------------------ |
| 1 rs6681049 1 1 0.1591 0.2667 2 3.067 0.07991 0.5203    |
| 1 rs4074137 2 1 0.07955 0.07778 2 0.001919 0.9651 1.025 |
| 1 rs7540009 3 0 0 0 2 NA NA NA                          |
| 1 rs1891905 4 1 0.4091 0.4 2 0.01527 0.9017 1.038       |
| 1 rs9729550 5 1 0.1705 0.08889 2 2.631 0.1048 2.106     |
| 1 rs3813196 6 1 0.03409 0.02222 2 0.2296 0.6318 1.553   |

``` bash
sort --key=7 -nr as1.assoc | head
plink --bfile hapmap1 --assoc --adjust --out as2
```

    ##    9    rs999510      47206    1   0.4091   0.3864    2      0.09488       0.7581          1.1 
    ##    9    rs999484      49016    1  0.02273  0.02222    2    0.0005167       0.9819        1.023 
    ##    9    rs999398      46425    1   0.1591   0.1889    2       0.2747       0.6002       0.8124 
    ##    9    rs998226      47266    1    0.375   0.4778    2        1.921       0.1657       0.6558 
    ##    9    rs997540      49756    1   0.3977   0.4333    2       0.2322       0.6299       0.8636 
    ##    9   rs9969732      48079    0        0        0    2           NA           NA           NA 
    ##    9   rs9969724      46038    0        0        0    2           NA           NA           NA 
    ##    9   rs9969710      47755    1   0.1023   0.1556    2        1.123       0.2893       0.6184 
    ##    9    rs995923      47709    1   0.2045   0.2444    2       0.4066       0.5237       0.7948 
    ##    9    rs995903      47398    1  0.04545      0.1    2        1.955        0.162       0.4286 
    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to as2.log.
    ## Options in effect:
    ##   --adjust
    ##   --assoc
    ##   --bfile hapmap1
    ##   --out as2
    ## 
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.99441.
    ## 83534 variants and 89 people pass filters and QC.
    ## Among remaining phenotypes, 44 are cases and 45 are controls.
    ## Writing C/C --assoc report to as2.assoc ... 0%1%2%3%4%6%7%8%9%10%11%12%14%15%16%17%18%19%20%22%23%25%26%27%28%29%31%32%33%34%35%36%38%39%40%41%42%43%45%46%47%48%49%50%51%53%54%55%56%58%59%60%61%62%63%65%66%67%68%69%70%72%73%74%76%77%78%79%80%81%83%84%85%87%88%89%90%92%93%94%95%96%97%98%99%done.
    ## --adjust: Genomic inflation est. lambda (based on median chisq) = 1.25377.
    ## 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%--adjust values (68727 variants) written to as2.assoc.adjusted .

``` r
file5 <- read.table("as2.assoc.adjusted", header = TRUE, sep="\t")
kable(head(file5))
```

| CHR………SNP……UNADJ………GC…….BONF…….HOLM…SIDAK\_SS…SIDAK\_SD…..FDR\_BH…..FDR\_BY |
| :-------------------------------------------------------------------------- |
| 13 rs9585021 5.586e-06 4.994e-05 0.3839 0.3839 0.3188 0.3188 0.09719 1      |
| 2 rs2222162 5.918e-06 5.232e-05 0.4068 0.4067 0.3342 0.3342 0.09719 1       |
| 9 rs10810856 7.723e-06 6.483e-05 0.5308 0.5308 0.4118 0.4118 0.09719 1      |
| 2 rs4675607 8.05e-06 6.703e-05 0.5533 0.5533 0.4249 0.4249 0.09719 1        |
| 2 rs4673349 8.485e-06 6.994e-05 0.5832 0.5831 0.4419 0.4419 0.09719 1       |
| 2 rs1375352 8.485e-06 6.994e-05 0.5832 0.5831 0.4419 0.4419 0.09719 1       |

``` bash
#Calculate association statistics based of genotypic and association models
plink --bfile hapmap1 --pheno pop.phe --assoc --adjust --out as3
plink --bfile hapmap1 --model --snp rs2222162 --out mod1
plink --bfile hapmap1 --model --cell 0 --snp rs2222162 --out mod2
#Perform stratification analysis
plink --bfile hapmap1 --cluster --mc 2 --ppc 0.05 --out str1
```

    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to as3.log.
    ## Options in effect:
    ##   --adjust
    ##   --assoc
    ##   --bfile hapmap1
    ##   --out as3
    ##   --pheno pop.phe
    ## 
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values present after --pheno.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.99441.
    ## 83534 variants and 89 people pass filters and QC.
    ## Among remaining phenotypes, 44 are cases and 45 are controls.
    ## Writing C/C --assoc report to as3.assoc ... 0%1%2%3%4%6%7%8%9%10%11%12%14%15%16%17%18%19%20%22%23%25%26%27%28%29%31%32%33%34%35%36%38%39%40%41%42%43%45%46%47%48%49%50%51%53%54%55%56%58%59%60%61%62%63%65%66%67%68%69%70%72%73%74%76%77%78%79%80%81%83%84%85%87%88%89%90%92%93%94%95%96%97%98%99%done.
    ## --adjust: Genomic inflation est. lambda (based on median chisq) = 1.78854.
    ## 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%--adjust values (68727 variants) written to as3.assoc.adjusted .
    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to mod1.log.
    ## Options in effect:
    ##   --bfile hapmap1
    ##   --model
    ##   --out mod1
    ##   --snp rs2222162
    ## 
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## 1 out of 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## 1 variant and 89 people pass filters and QC.
    ## Among remaining phenotypes, 44 are cases and 45 are controls.
    ## Writing --model report to mod1.model ... 0%done.
    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to mod2.log.
    ## Options in effect:
    ##   --bfile hapmap1
    ##   --cell 0
    ##   --model
    ##   --out mod2
    ##   --snp rs2222162
    ## 
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## 1 out of 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## 1 variant and 89 people pass filters and QC.
    ## Among remaining phenotypes, 44 are cases and 45 are controls.
    ## Writing --model report to mod2.model ... 0%done.
    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to str1.log.
    ## Options in effect:
    ##   --bfile hapmap1
    ##   --cluster
    ##   --mc 2
    ##   --out str1
    ##   --ppc 0.05
    ## 
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values loaded from .fam.
    ## Using up to 47 threads (change this with --threads).
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.99441.
    ## 83534 variants and 89 people pass filters and QC.
    ## Among remaining phenotypes, 44 are cases and 45 are controls.
    ## 1152 markers complete.2304 markers complete.3456 markers complete.4608 markers complete.5760 markers complete.6912 markers complete.8064 markers complete.9216 markers complete.10368 markers complete.11520 markers complete.12672 markers complete.13824 markers complete.14976 markers complete.16128 markers complete.17280 markers complete.18432 markers complete.19584 markers complete.20736 markers complete.21888 markers complete.23040 markers complete.24192 markers complete.25344 markers complete.26496 markers complete.27648 markers complete.28800 markers complete.29952 markers complete.31104 markers complete.32256 markers complete.33408 markers complete.34560 markers complete.35712 markers complete.36864 markers complete.38016 markers complete.39168 markers complete.40320 markers complete.41472 markers complete.42624 markers complete.43776 markers complete.44928 markers complete.46080 markers complete.47232 markers complete.48384 markers complete.49536 markers complete.50688 markers complete.51840 markers complete.52992 markers complete.54144 markers complete.55296 markers complete.56448 markers complete.57600 markers complete.58752 markers complete.59904 markers complete.61056 markers complete.62208 markers complete.63360 markers complete.64512 markers complete.65664 markers complete.66816 markers complete.67968 markers complete.69120 markers complete.70272 markers complete.71424 markers complete.72576 markers complete.73728 markers complete.74880 markers complete.76032 markers complete.77184 markers complete.78336 markers complete.79488 markers complete.80640 markers complete.81792 markers complete.82944 markers complete.83534 markers complete.IBD calculations complete.  
    ## Clustering... [sorting IBS values]Clustering... done.                        
    ## Writing cluster solution... 0%1%2%3%4%5%6%7%8%10%11%12%13%14%15%16%17%19%20%21%22%23%24%25%26%28%29%30%31%32%33%34%35%37%38%39%40%41%42%43%44%46%47%48%49%50%51%52%53%55%56%57%58%59%60%61%62%64%65%66%67%68%69%70%71%73%74%75%76%77%78%79%80%82%83%84%85%86%87%88%89%91%92%93%94%95%96%97%98%100%Cluster solution written to str1.cluster1 , str1.cluster2 , and str1.cluster3 .

``` r
file6 <- read.table("str1.cluster1", header = TRUE, sep="\t")
kable(head(file6))
```

| SOL.0 | HCB181\_1.JPT260\_1 |
| :---- | :------------------ |
| SOL-1 | HCB182\_1 HCB225\_1 |
| SOL-2 | HCB183\_1 HCB194\_1 |
| SOL-3 | HCB184\_1 HCB202\_1 |
| SOL-4 | HCB185\_1 HCB217\_1 |
| SOL-5 | HCB186\_1 HCB201\_1 |
| SOL-6 | HCB187\_1 HCB189\_1 |

``` bash
#Association analysis, accounting for clusters
plink --bfile hapmap1 --mh --within str1.cluster2 --adjust --out aac1
```

    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to aac1.log.
    ## Options in effect:
    ##   --adjust
    ##   --bfile hapmap1
    ##   --mh
    ##   --out aac1
    ##   --within str1.cluster2
    ## 
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values loaded from .fam.
    ## --within: 45 clusters loaded, covering a total of 89 people.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.99441.
    ## 83534 variants and 89 people pass filters and QC.
    ## Among remaining phenotypes, 44 are cases and 45 are controls.
    ## --mh/--bd: 21 valid clusters, with a total of 21 cases and 21 controls.
    ## Writing report to aac1.cmh ... 0%0%1%1%2%2%3%3%4%4%5%5%6%6%7%7%8%8%9%9%10%10%11%11%12%12%13%13%14%14%15%15%16%16%17%17%18%18%19%19%20%20%21%21%22%22%23%23%24%24%25%25%26%26%27%27%28%28%29%29%30%30%31%31%32%32%33%33%34%34%35%35%36%36%37%37%38%38%39%39%40%40%41%41%42%42%43%43%44%44%45%45%46%46%47%47%48%48%49%50%50%51%51%52%52%53%53%54%54%55%55%56%56%57%57%58%58%59%59%60%60%61%61%62%62%63%63%64%64%65%65%66%66%67%67%68%68%69%69%70%70%71%71%72%72%73%73%74%74%75%75%76%76%77%77%78%78%79%79%80%80%81%81%82%82%83%83%84%84%85%85%86%86%87%87%88%88%89%89%90%90%91%91%92%92%93%93%94%94%95%95%96%96%97%97%98%98%99%done.
    ## --adjust: Genomic inflation est. lambda (based on median chisq) = 1.07656.
    ## 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%--adjust values (66852 variants) written to aac1.cmh.adjusted .

``` r
file7 <- read.table("aac1.cmh.adjusted", header = TRUE, sep ="\t")
kable(head(file7))
```

| CHR………SNP……UNADJ………GC…….BONF…….HOLM…SIDAK\_SS…SIDAK\_SD…..FDR\_BH…..FDR\_BY |
| :-------------------------------------------------------------------------- |
| 13 rs9585021 1.906e-06 4.418e-06 0.1274 0.1274 0.1196 0.1196 0.1274 1       |
| 21 rs3017432 2.209e-05 4.332e-05 1 1 0.7716 0.7716 0.7384 1                 |
| 2 rs2222162 4.468e-05 8.353e-05 1 1 0.9496 0.9495 0.8734 1                  |
| 17 rs3829612 7.177e-05 0.0001299 1 1 0.9918 0.9918 0.8734 1                 |
| 2 rs4673349 9.617e-05 0.0001707 1 1 0.9984 0.9984 0.8734 1                  |
| 2 rs1375352 9.617e-05 0.0001707 1 1 0.9984 0.9984 0.8734 1                  |

``` bash
plink --bfile hapmap1 --cluster --cc --ppc 0.01 --out version2
plink --bfile hapmap1 --mh --within version2.cluster2 --adjust --out aac2
plink --bfile hapmap1 --cluster --K 2 --out version3
```

    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to version2.log.
    ## Options in effect:
    ##   --bfile hapmap1
    ##   --cc
    ##   --cluster
    ##   --out version2
    ##   --ppc 0.01
    ## 
    ## Note: --cc flag deprecated.  Use '--cluster cc'.
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values loaded from .fam.
    ## Using up to 47 threads (change this with --threads).
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.99441.
    ## 83534 variants and 89 people pass filters and QC.
    ## Among remaining phenotypes, 44 are cases and 45 are controls.
    ## 1152 markers complete.2304 markers complete.3456 markers complete.4608 markers complete.5760 markers complete.6912 markers complete.8064 markers complete.9216 markers complete.10368 markers complete.11520 markers complete.12672 markers complete.13824 markers complete.14976 markers complete.16128 markers complete.17280 markers complete.18432 markers complete.19584 markers complete.20736 markers complete.21888 markers complete.23040 markers complete.24192 markers complete.25344 markers complete.26496 markers complete.27648 markers complete.28800 markers complete.29952 markers complete.31104 markers complete.32256 markers complete.33408 markers complete.34560 markers complete.35712 markers complete.36864 markers complete.38016 markers complete.39168 markers complete.40320 markers complete.41472 markers complete.42624 markers complete.43776 markers complete.44928 markers complete.46080 markers complete.47232 markers complete.48384 markers complete.49536 markers complete.50688 markers complete.51840 markers complete.52992 markers complete.54144 markers complete.55296 markers complete.56448 markers complete.57600 markers complete.58752 markers complete.59904 markers complete.61056 markers complete.62208 markers complete.63360 markers complete.64512 markers complete.65664 markers complete.66816 markers complete.67968 markers complete.69120 markers complete.70272 markers complete.71424 markers complete.72576 markers complete.73728 markers complete.74880 markers complete.76032 markers complete.77184 markers complete.78336 markers complete.79488 markers complete.80640 markers complete.81792 markers complete.82944 markers complete.83534 markers complete.IBD calculations complete.  
    ## Clustering... [sorting IBS values]Clustering... done.                        
    ## Writing cluster solution... 0%1%2%3%4%5%6%7%8%10%11%12%13%14%15%16%17%19%20%21%22%23%24%25%26%28%29%30%31%32%33%34%35%37%38%39%40%41%42%43%44%46%47%48%49%50%51%52%53%55%56%57%58%59%60%61%62%64%65%66%67%68%69%70%71%73%74%75%76%77%78%79%80%82%83%84%85%86%87%88%89%91%92%93%94%95%96%97%98%100%Cluster solution written to version2.cluster1 , version2.cluster2 , and
    ## version2.cluster3 .
    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to aac2.log.
    ## Options in effect:
    ##   --adjust
    ##   --bfile hapmap1
    ##   --mh
    ##   --out aac2
    ##   --within version2.cluster2
    ## 
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values loaded from .fam.
    ## --within: 5 clusters loaded, covering a total of 89 people.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.99441.
    ## 83534 variants and 89 people pass filters and QC.
    ## Among remaining phenotypes, 44 are cases and 45 are controls.
    ## --mh/--bd: 5 valid clusters, with a total of 44 cases and 45 controls.
    ## Writing report to aac2.cmh ... 0%0%1%1%2%2%3%3%4%4%5%5%6%6%7%7%8%8%9%9%10%10%11%11%12%12%13%13%14%14%15%15%16%16%17%17%18%18%19%19%20%20%21%21%22%22%23%23%24%24%25%25%26%26%27%27%28%28%29%29%30%30%31%31%32%32%33%33%34%34%35%35%36%36%37%37%38%38%39%39%40%40%41%41%42%42%43%43%44%44%45%45%46%46%47%47%48%48%49%50%50%51%51%52%52%53%53%54%54%55%55%56%56%57%57%58%58%59%59%60%60%61%61%62%62%63%63%64%64%65%65%66%66%67%67%68%68%69%69%70%70%71%71%72%72%73%73%74%74%75%75%76%76%77%77%78%78%79%79%80%80%81%81%82%82%83%83%84%84%85%85%86%86%87%87%88%88%89%89%90%90%91%91%92%92%93%93%94%94%95%95%96%96%97%97%98%98%99%done.
    ## --adjust: Genomic inflation est. lambda (based on median chisq) = 1.027.
    ## 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%--adjust values (68727 variants) written to aac2.cmh.adjusted .
    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to version3.log.
    ## Options in effect:
    ##   --K 2
    ##   --bfile hapmap1
    ##   --cluster
    ##   --out version3
    ## 
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values loaded from .fam.
    ## Using up to 47 threads (change this with --threads).
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.99441.
    ## 83534 variants and 89 people pass filters and QC.
    ## Among remaining phenotypes, 44 are cases and 45 are controls.
    ## 960 markers complete.1920 markers complete.2880 markers complete.3840 markers complete.4800 markers complete.5760 markers complete.6720 markers complete.7680 markers complete.8640 markers complete.9600 markers complete.10560 markers complete.11520 markers complete.12480 markers complete.13440 markers complete.14400 markers complete.15360 markers complete.16320 markers complete.17280 markers complete.18240 markers complete.19200 markers complete.20160 markers complete.21120 markers complete.22080 markers complete.23040 markers complete.24000 markers complete.24960 markers complete.25920 markers complete.26880 markers complete.27840 markers complete.28800 markers complete.29760 markers complete.30720 markers complete.31680 markers complete.32640 markers complete.33600 markers complete.34560 markers complete.35520 markers complete.36480 markers complete.37440 markers complete.38400 markers complete.39360 markers complete.40320 markers complete.41280 markers complete.42240 markers complete.43200 markers complete.44160 markers complete.45120 markers complete.46080 markers complete.47040 markers complete.48000 markers complete.48960 markers complete.49920 markers complete.50880 markers complete.51840 markers complete.52800 markers complete.53760 markers complete.54720 markers complete.55680 markers complete.56640 markers complete.57600 markers complete.58560 markers complete.59520 markers complete.60480 markers complete.61440 markers complete.62400 markers complete.63360 markers complete.64320 markers complete.65280 markers complete.66240 markers complete.67200 markers complete.68160 markers complete.69120 markers complete.70080 markers complete.71040 markers complete.72000 markers complete.72960 markers complete.73920 markers complete.74880 markers complete.75840 markers complete.76800 markers complete.77760 markers complete.78720 markers complete.79680 markers complete.80640 markers complete.81600 markers complete.82560 markers complete.83520 markers complete.83534 markers complete.Distance matrix calculation complete.
    ## Clustering... [sorting IBS values]Clustering... done.                        
    ## Writing cluster solution... 0%1%2%3%4%5%6%7%8%10%11%12%13%14%15%16%17%19%20%21%22%23%24%25%26%28%29%30%31%32%33%34%35%37%38%39%40%41%42%43%44%46%47%48%49%50%51%52%53%55%56%57%58%59%60%61%62%64%65%66%67%68%69%70%71%73%74%75%76%77%78%79%80%82%83%84%85%86%87%88%89%91%92%93%94%95%96%97%98%100%Cluster solution written to version3.cluster1 , version3.cluster2 , and
    ## version3.cluster3 .

``` r
file8<- read.table("aac2.cmh.adjusted", header = TRUE, sep="\t")
kable(head(file8))
```

| CHR………SNP……UNADJ………GC…….BONF…….HOLM…SIDAK\_SS…SIDAK\_SD…..FDR\_BH…..FDR\_BY          |
| :----------------------------------------------------------------------------------- |
| 2 rs2222162 1.474e-08 2.276e-08 0.001013 0.001013 0.001013 0.001013 0.001013 0.01187 |
| 2 rs4675607 1.134e-05 1.479e-05 0.7794 0.7794 0.5413 0.5413 0.2779 1                 |
| 13 rs9585021 1.213e-05 1.58e-05 0.8338 0.8338 0.5656 0.5656 0.2779 1                 |
| 9 rs7046471 3.372e-05 4.278e-05 1 1 0.9015 0.9015 0.4296 1                           |
| 2 rs4673349 3.75e-05 4.746e-05 1 1 0.924 0.924 0.4296 1                              |
| 2 rs1375352 3.75e-05 4.746e-05 1 1 0.924 0.924 0.4296 1                              |

``` bash
plink --bfile hapmap1 --mh --within pop.phe --adjust --out aac3
plink --bfile hapmap1 --cluster --matrix --out ibd_view
```

    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to aac3.log.
    ## Options in effect:
    ##   --adjust
    ##   --bfile hapmap1
    ##   --mh
    ##   --out aac3
    ##   --within pop.phe
    ## 
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values loaded from .fam.
    ## --within: 2 clusters loaded, covering a total of 89 people.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.99441.
    ## 83534 variants and 89 people pass filters and QC.
    ## Among remaining phenotypes, 44 are cases and 45 are controls.
    ## --mh/--bd: 2 valid clusters, with a total of 44 cases and 45 controls.
    ## Writing report to aac3.cmh ... 0%0%1%1%2%2%3%3%4%4%5%5%6%6%7%7%8%8%9%9%10%10%11%11%12%12%13%13%14%14%15%15%16%16%17%17%18%18%19%19%20%20%21%21%22%22%23%23%24%24%25%25%26%26%27%27%28%28%29%29%30%30%31%31%32%32%33%33%34%34%35%35%36%36%37%37%38%38%39%39%40%40%41%41%42%42%43%43%44%44%45%45%46%46%47%47%48%48%49%50%50%51%51%52%52%53%53%54%54%55%55%56%56%57%57%58%58%59%59%60%60%61%61%62%62%63%63%64%64%65%65%66%66%67%67%68%68%69%69%70%70%71%71%72%72%73%73%74%74%75%75%76%76%77%77%78%78%79%79%80%80%81%81%82%82%83%83%84%84%85%85%86%86%87%87%88%88%89%89%90%90%91%91%92%92%93%93%94%94%95%95%96%96%97%97%98%98%99%done.
    ## --adjust: Genomic inflation est. lambda (based on median chisq) = 1.
    ## 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%--adjust values (68727 variants) written to aac3.cmh.adjusted .
    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to ibd_view.log.
    ## Options in effect:
    ##   --bfile hapmap1
    ##   --cluster
    ##   --matrix
    ##   --out ibd_view
    ## 
    ## Note: --matrix flag deprecated.  Migrate to '--distance ibs flat-missing',
    ## '--r2 square', etc.
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values loaded from .fam.
    ## Using up to 47 threads (change this with --threads).
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.99441.
    ## 83534 variants and 89 people pass filters and QC.
    ## Among remaining phenotypes, 44 are cases and 45 are controls.
    ## 960 markers complete.1920 markers complete.2880 markers complete.3840 markers complete.4800 markers complete.5760 markers complete.6720 markers complete.7680 markers complete.8640 markers complete.9600 markers complete.10560 markers complete.11520 markers complete.12480 markers complete.13440 markers complete.14400 markers complete.15360 markers complete.16320 markers complete.17280 markers complete.18240 markers complete.19200 markers complete.20160 markers complete.21120 markers complete.22080 markers complete.23040 markers complete.24000 markers complete.24960 markers complete.25920 markers complete.26880 markers complete.27840 markers complete.28800 markers complete.29760 markers complete.30720 markers complete.31680 markers complete.32640 markers complete.33600 markers complete.34560 markers complete.35520 markers complete.36480 markers complete.37440 markers complete.38400 markers complete.39360 markers complete.40320 markers complete.41280 markers complete.42240 markers complete.43200 markers complete.44160 markers complete.45120 markers complete.46080 markers complete.47040 markers complete.48000 markers complete.48960 markers complete.49920 markers complete.50880 markers complete.51840 markers complete.52800 markers complete.53760 markers complete.54720 markers complete.55680 markers complete.56640 markers complete.57600 markers complete.58560 markers complete.59520 markers complete.60480 markers complete.61440 markers complete.62400 markers complete.63360 markers complete.64320 markers complete.65280 markers complete.66240 markers complete.67200 markers complete.68160 markers complete.69120 markers complete.70080 markers complete.71040 markers complete.72000 markers complete.72960 markers complete.73920 markers complete.74880 markers complete.75840 markers complete.76800 markers complete.77760 markers complete.78720 markers complete.79680 markers complete.80640 markers complete.81600 markers complete.82560 markers complete.83520 markers complete.83534 markers complete.Distance matrix calculation complete.
    ## Writing... 1%Writing... 2%Writing... 3%Writing... 4%Writing... 5%Writing... 6%Writing... 7%Writing... 8%Writing... 10%Writing... 11%Writing... 12%Writing... 13%Writing... 14%Writing... 15%Writing... 16%Writing... 17%Writing... 19%Writing... 20%Writing... 21%Writing... 22%Writing... 23%Writing... 24%Writing... 25%Writing... 26%Writing... 28%Writing... 29%Writing... 30%Writing... 31%Writing... 32%Writing... 33%Writing... 34%Writing... 35%Writing... 37%Writing... 38%Writing... 39%Writing... 40%Writing... 41%Writing... 42%Writing... 43%Writing... 44%Writing... 46%Writing... 47%Writing... 48%Writing... 49%Writing... 50%Writing... 51%Writing... 52%Writing... 53%Writing... 55%Writing... 56%Writing... 57%Writing... 58%Writing... 59%Writing... 60%Writing... 61%Writing... 62%Writing... 64%Writing... 65%Writing... 66%Writing... 67%Writing... 68%Writing... 69%Writing... 70%Writing... 71%Writing... 73%Writing... 74%Writing... 75%Writing... 76%Writing... 77%Writing... 78%Writing... 79%Writing... 80%Writing... 82%Writing... 83%Writing... 84%Writing... 85%Writing... 86%Writing... 87%Writing... 88%Writing... 89%Writing... 91%Writing... 92%Writing... 93%Writing... 94%Writing... 95%Writing... 96%Writing... 97%Writing... 98%IBS matrix written to ibd_view.mibs , and IDs to ibd_view.mibs.id .
    ## Clustering... [sorting IBS values]Clustering... done.                        
    ## Writing cluster solution... 0%1%2%3%4%5%6%7%8%10%11%12%13%14%15%16%17%19%20%21%22%23%24%25%26%28%29%30%31%32%33%34%35%37%38%39%40%41%42%43%44%46%47%48%49%50%51%52%53%55%56%57%58%59%60%61%62%64%65%66%67%68%69%70%71%73%74%75%76%77%78%79%80%82%83%84%85%86%87%88%89%91%92%93%94%95%96%97%98%100%Cluster solution written to ibd_view.cluster1 , ibd_view.cluster2 , and
    ## ibd_view.cluster3 .

![](runPlink_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` bash
#Quantitative trait association analysis
plink --bfile hapmap1 --assoc --pheno qt.phe --out quant1
```

    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to quant1.log.
    ## Options in effect:
    ##   --assoc
    ##   --bfile hapmap1
    ##   --out quant1
    ##   --pheno qt.phe
    ## 
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values present after --pheno.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.99441.
    ## 83534 variants and 89 people pass filters and QC.
    ## Phenotype data is quantitative.
    ## Writing QT --assoc report to quant1.qassoc ... 0%1%2%3%4%6%7%8%9%10%11%12%14%15%16%17%18%19%20%22%23%25%26%27%28%29%31%32%33%34%35%36%38%39%40%41%42%43%45%46%47%48%49%50%51%53%54%55%56%58%59%60%61%62%63%65%66%67%68%69%70%72%73%74%76%77%78%79%80%81%83%84%85%87%88%89%90%92%93%94%95%96%97%98%99%done.

``` bash
plink --bfile hapmap1 --assoc --pheno qt.phe --perm --within str1.cluster2 --out quant2
```

    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to quant2.log.
    ## Options in effect:
    ##   --assoc
    ##   --bfile hapmap1
    ##   --out quant2
    ##   --perm
    ##   --pheno qt.phe
    ##   --within str1.cluster2
    ## 
    ## Note: --perm flag deprecated.  Use e.g. '--model perm'.
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values present after --pheno.
    ## --within: 45 clusters loaded, covering a total of 89 people.
    ## Using up to 47 threads (change this with --threads).
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.99441.
    ## 83534 variants and 89 people pass filters and QC.
    ## Phenotype data is quantitative.
    ## Writing QT --assoc report to quant2.qassoc ... 0%1%2%3%4%6%7%8%9%10%11%12%14%15%16%17%18%19%20%22%23%25%26%27%28%29%31%32%33%34%35%36%38%39%40%41%42%43%45%46%47%48%49%50%51%53%54%55%56%58%59%60%61%62%63%65%66%67%68%69%70%72%73%74%76%77%78%79%80%81%83%84%85%87%88%89%90%92%93%94%95%96%97%98%99%done.
    ## 512 permutations complete.1024 permutations complete.1536 permutations complete.2048 permutations complete.2560 permutations complete.3072 permutations complete.3584 permutations complete.4096 permutations complete.4608 permutations complete.5120 permutations complete.5632 permutations complete.6144 permutations complete.6656 permutations complete.7168 permutations complete.7680 permutations complete.8192 permutations complete.8704 permutations complete.9216 permutations complete.9728 permutations complete.10240 permutations complete.10752 permutations complete.11264 permutations complete.11776 permutations complete.12288 permutations complete.12800 permutations complete.13312 permutations complete.13824 permutations complete.14336 permutations complete.14848 permutations complete.15360 permutations complete.15872 permutations complete.16384 permutations complete.16896 permutations complete.17408 permutations complete.17920 permutations complete.18432 permutations complete.18944 permutations complete.19456 permutations complete.19968 permutations complete.20480 permutations complete.20992 permutations complete.21504 permutations complete.22016 permutations complete.22528 permutations complete.23040 permutations complete.23552 permutations complete.24064 permutations complete.24576 permutations complete.25088 permutations complete.25600 permutations complete.26112 permutations complete.26624 permutations complete.27136 permutations complete.27648 permutations complete.28160 permutations complete.28672 permutations complete.29184 permutations complete.29696 permutations complete.30208 permutations complete.30720 permutations complete.31232 permutations complete.31744 permutations complete.32256 permutations complete.32768 permutations complete.33280 permutations complete.33792 permutations complete.34304 permutations complete.34816 permutations complete.35328 permutations complete.35840 permutations complete.36352 permutations complete.36864 permutations complete.37376 permutations complete.37888 permutations complete.38400 permutations complete.38912 permutations complete.39424 permutations complete.39936 permutations complete.40448 permutations complete.40960 permutations complete.41472 permutations complete.41984 permutations complete.42496 permutations complete.43008 permutations complete.43520 permutations complete.44032 permutations complete.44544 permutations complete.45056 permutations complete.45568 permutations complete.46080 permutations complete.46592 permutations complete.47104 permutations complete.47616 permutations complete.48128 permutations complete.48640 permutations complete.49152 permutations complete.49664 permutations complete.50176 permutations complete.50688 permutations complete.51200 permutations complete.51712 permutations complete.52224 permutations complete.52736 permutations complete.53248 permutations complete.53760 permutations complete.54272 permutations complete.54784 permutations complete.55296 permutations complete.55808 permutations complete.56320 permutations complete.56832 permutations complete.57344 permutations complete.57856 permutations complete.58368 permutations complete.58880 permutations complete.59392 permutations complete.59904 permutations complete.60416 permutations complete.60928 permutations complete.61440 permutations complete.61952 permutations complete.62464 permutations complete.62976 permutations complete.63488 permutations complete.64000 permutations complete.64512 permutations complete.65024 permutations complete.65536 permutations complete.66048 permutations complete.66560 permutations complete.67072 permutations complete.67584 permutations complete.68096 permutations complete.68608 permutations complete.69120 permutations complete.69632 permutations complete.70144 permutations complete.70656 permutations complete.71168 permutations complete.71680 permutations complete.72192 permutations complete.72704 permutations complete.73216 permutations complete.73728 permutations complete.74240 permutations complete.74752 permutations complete.75264 permutations complete.75776 permutations complete.76288 permutations complete.76800 permutations complete.77312 permutations complete.77824 permutations complete.78336 permutations complete.78848 permutations complete.79360 permutations complete.79872 permutations complete.80384 permutations complete.80896 permutations complete.81408 permutations complete.81920 permutations complete.82432 permutations complete.82944 permutations complete.83456 permutations complete.83968 permutations complete.84480 permutations complete.84992 permutations complete.85504 permutations complete.86016 permutations complete.86528 permutations complete.87040 permutations complete.87552 permutations complete.88064 permutations complete.88576 permutations complete.89088 permutations complete.89600 permutations complete.90112 permutations complete.90624 permutations complete.91136 permutations complete.91648 permutations complete.92160 permutations complete.92672 permutations complete.93184 permutations complete.93696 permutations complete.94208 permutations complete.94720 permutations complete.95232 permutations complete.95744 permutations complete.96256 permutations complete.96768 permutations complete.97280 permutations complete.97792 permutations complete.98304 permutations complete.98816 permutations complete.99328 permutations complete.99840 permutations complete.100352 permutations complete.100864 permutations complete.101376 permutations complete.101888 permutations complete.102400 permutations complete.102912 permutations complete.103424 permutations complete.103936 permutations complete.104448 permutations complete.104960 permutations complete.105472 permutations complete.105984 permutations complete.106496 permutations complete.107008 permutations complete.107520 permutations complete.108032 permutations complete.108544 permutations complete.109056 permutations complete.109568 permutations complete.110080 permutations complete.110592 permutations complete.111104 permutations complete.111616 permutations complete.112128 permutations complete.112640 permutations complete.113152 permutations complete.113664 permutations complete.114176 permutations complete.114688 permutations complete.115200 permutations complete.115712 permutations complete.116224 permutations complete.116736 permutations complete.117248 permutations complete.117760 permutations complete.118272 permutations complete.118784 permutations complete.119296 permutations complete.119808 permutations complete.120320 permutations complete.120832 permutations complete.121344 permutations complete.121856 permutations complete.122368 permutations complete.122880 permutations complete.123392 permutations complete.123904 permutations complete.124416 permutations complete.124928 permutations complete.125440 permutations complete.125952 permutations complete.126464 permutations complete.126976 permutations complete.127488 permutations complete.128000 permutations complete.128512 permutations complete.129024 permutations complete.129536 permutations complete.130048 permutations complete.130560 permutations complete.131072 permutations complete.131584 permutations complete.132096 permutations complete.132608 permutations complete.133120 permutations complete.133632 permutations complete.134144 permutations complete.134656 permutations complete.135168 permutations complete.135680 permutations complete.136192 permutations complete.136704 permutations complete.137216 permutations complete.137728 permutations complete.138240 permutations complete.138752 permutations complete.139264 permutations complete.139776 permutations complete.140288 permutations complete.140800 permutations complete.141312 permutations complete.141824 permutations complete.142336 permutations complete.142848 permutations complete.143360 permutations complete.143872 permutations complete.144384 permutations complete.144896 permutations complete.145408 permutations complete.145920 permutations complete.146432 permutations complete.146944 permutations complete.147456 permutations complete.147968 permutations complete.148480 permutations complete.148992 permutations complete.149504 permutations complete.150016 permutations complete.150528 permutations complete.151040 permutations complete.151552 permutations complete.152064 permutations complete.152576 permutations complete.153088 permutations complete.153600 permutations complete.154112 permutations complete.154624 permutations complete.155136 permutations complete.155648 permutations complete.156160 permutations complete.156672 permutations complete.157184 permutations complete.157696 permutations complete.158208 permutations complete.158720 permutations complete.159232 permutations complete.159744 permutations complete.160256 permutations complete.160768 permutations complete.161280 permutations complete.161792 permutations complete.162304 permutations complete.162816 permutations complete.163328 permutations complete.163840 permutations complete.164352 permutations complete.164864 permutations complete.165376 permutations complete.165888 permutations complete.166400 permutations complete.166912 permutations complete.167424 permutations complete.167936 permutations complete.168448 permutations complete.168960 permutations complete.169472 permutations complete.169984 permutations complete.170496 permutations complete.171008 permutations complete.171520 permutations complete.172032 permutations complete.172544 permutations complete.173056 permutations complete.173568 permutations complete.174080 permutations complete.174592 permutations complete.175104 permutations complete.175616 permutations complete.176128 permutations complete.176640 permutations complete.177152 permutations complete.177664 permutations complete.178176 permutations complete.178688 permutations complete.179200 permutations complete.179712 permutations complete.180224 permutations complete.180736 permutations complete.181248 permutations complete.181760 permutations complete.182272 permutations complete.182784 permutations complete.183296 permutations complete.183808 permutations complete.184320 permutations complete.184832 permutations complete.185344 permutations complete.185856 permutations complete.186368 permutations complete.186880 permutations complete.187392 permutations complete.187904 permutations complete.188416 permutations complete.188928 permutations complete.189440 permutations complete.189952 permutations complete.190464 permutations complete.190976 permutations complete.191488 permutations complete.192000 permutations complete.192512 permutations complete.193024 permutations complete.193536 permutations complete.194048 permutations complete.194560 permutations complete.195072 permutations complete.195584 permutations complete.196096 permutations complete.196608 permutations complete.197120 permutations complete.197632 permutations complete.198144 permutations complete.198656 permutations complete.199168 permutations complete.199680 permutations complete.200192 permutations complete.200704 permutations complete.201216 permutations complete.201728 permutations complete.202240 permutations complete.202752 permutations complete.203264 permutations complete.203776 permutations complete.204288 permutations complete.204800 permutations complete.205312 permutations complete.205824 permutations complete.206336 permutations complete.206848 permutations complete.207360 permutations complete.207872 permutations complete.208384 permutations complete.208896 permutations complete.209408 permutations complete.209920 permutations complete.210432 permutations complete.210944 permutations complete.211456 permutations complete.211968 permutations complete.212480 permutations complete.212992 permutations complete.213504 permutations complete.214016 permutations complete.214528 permutations complete.215040 permutations complete.215552 permutations complete.216064 permutations complete.216576 permutations complete.217088 permutations complete.217600 permutations complete.218112 permutations complete.218624 permutations complete.219136 permutations complete.219648 permutations complete.220160 permutations complete.220672 permutations complete.221184 permutations complete.221696 permutations complete.222208 permutations complete.222720 permutations complete.223232 permutations complete.223744 permutations complete.224256 permutations complete.224768 permutations complete.225280 permutations complete.225792 permutations complete.226304 permutations complete.226816 permutations complete.227328 permutations complete.227840 permutations complete.228352 permutations complete.228864 permutations complete.229376 permutations complete.229888 permutations complete.230400 permutations complete.230912 permutations complete.231424 permutations complete.231936 permutations complete.232448 permutations complete.232960 permutations complete.233472 permutations complete.233984 permutations complete.234496 permutations complete.235008 permutations complete.235520 permutations complete.236032 permutations complete.236544 permutations complete.237056 permutations complete.237568 permutations complete.238080 permutations complete.238592 permutations complete.239104 permutations complete.239616 permutations complete.240128 permutations complete.240640 permutations complete.241152 permutations complete.241664 permutations complete.242176 permutations complete.242688 permutations complete.243200 permutations complete.243712 permutations complete.244224 permutations complete.244736 permutations complete.245248 permutations complete.245760 permutations complete.246272 permutations complete.246784 permutations complete.247296 permutations complete.247808 permutations complete.248320 permutations complete.248832 permutations complete.249344 permutations complete.249856 permutations complete.250368 permutations complete.250880 permutations complete.251392 permutations complete.251904 permutations complete.252416 permutations complete.252928 permutations complete.253440 permutations complete.253952 permutations complete.254464 permutations complete.254976 permutations complete.255488 permutations complete.256000 permutations complete.256512 permutations complete.257024 permutations complete.257536 permutations complete.258048 permutations complete.258560 permutations complete.259072 permutations complete.259584 permutations complete.260096 permutations complete.260608 permutations complete.261120 permutations complete.261632 permutations complete.262144 permutations complete.262656 permutations complete.263168 permutations complete.263680 permutations complete.264192 permutations complete.264704 permutations complete.265216 permutations complete.265728 permutations complete.266240 permutations complete.266752 permutations complete.267264 permutations complete.267776 permutations complete.268288 permutations complete.268800 permutations complete.269312 permutations complete.269824 permutations complete.270336 permutations complete.270848 permutations complete.271360 permutations complete.271872 permutations complete.272384 permutations complete.272896 permutations complete.273408 permutations complete.273920 permutations complete.274432 permutations complete.274944 permutations complete.275456 permutations complete.275968 permutations complete.276480 permutations complete.276992 permutations complete.277504 permutations complete.278016 permutations complete.278528 permutations complete.279040 permutations complete.279552 permutations complete.280064 permutations complete.280576 permutations complete.281088 permutations complete.281600 permutations complete.282112 permutations complete.282624 permutations complete.283136 permutations complete.283648 permutations complete.284160 permutations complete.284672 permutations complete.285184 permutations complete.285696 permutations complete.286208 permutations complete.286720 permutations complete.287232 permutations complete.287744 permutations complete.288256 permutations complete.288768 permutations complete.289280 permutations complete.289792 permutations complete.290304 permutations complete.290816 permutations complete.291328 permutations complete.291840 permutations complete.292352 permutations complete.292864 permutations complete.293376 permutations complete.293888 permutations complete.294400 permutations complete.294912 permutations complete.295424 permutations complete.295936 permutations complete.296448 permutations complete.296960 permutations complete.297472 permutations complete.297984 permutations complete.298496 permutations complete.299008 permutations complete.299520 permutations complete.300032 permutations complete.300544 permutations complete.301056 permutations complete.301568 permutations complete.302080 permutations complete.302592 permutations complete.303104 permutations complete.303616 permutations complete.304128 permutations complete.304640 permutations complete.305152 permutations complete.305664 permutations complete.306176 permutations complete.306688 permutations complete.307200 permutations complete.307712 permutations complete.308224 permutations complete.308736 permutations complete.309248 permutations complete.309760 permutations complete.310272 permutations complete.310784 permutations complete.311296 permutations complete.311808 permutations complete.312320 permutations complete.312832 permutations complete.313344 permutations complete.313856 permutations complete.314368 permutations complete.314880 permutations complete.315392 permutations complete.315904 permutations complete.316416 permutations complete.316928 permutations complete.317440 permutations complete.317952 permutations complete.318464 permutations complete.318976 permutations complete.319488 permutations complete.320000 permutations complete.320512 permutations complete.321024 permutations complete.321536 permutations complete.322048 permutations complete.322560 permutations complete.323072 permutations complete.323584 permutations complete.324096 permutations complete.324608 permutations complete.325120 permutations complete.325632 permutations complete.326144 permutations complete.326656 permutations complete.327168 permutations complete.327680 permutations complete.328192 permutations complete.328704 permutations complete.329216 permutations complete.329728 permutations complete.330240 permutations complete.330752 permutations complete.331264 permutations complete.331776 permutations complete.332288 permutations complete.332800 permutations complete.333312 permutations complete.333824 permutations complete.334336 permutations complete.334848 permutations complete.335360 permutations complete.335872 permutations complete.336384 permutations complete.336896 permutations complete.337408 permutations complete.337920 permutations complete.338432 permutations complete.338944 permutations complete.339456 permutations complete.339968 permutations complete.340480 permutations complete.340992 permutations complete.341504 permutations complete.342016 permutations complete.342528 permutations complete.343040 permutations complete.343552 permutations complete.344064 permutations complete.344576 permutations complete.345088 permutations complete.345600 permutations complete.346112 permutations complete.346624 permutations complete.347136 permutations complete.347648 permutations complete.348160 permutations complete.348672 permutations complete.349184 permutations complete.349696 permutations complete.350208 permutations complete.350720 permutations complete.351232 permutations complete.351744 permutations complete.352256 permutations complete.352768 permutations complete.353280 permutations complete.353792 permutations complete.354304 permutations complete.354816 permutations complete.355328 permutations complete.355840 permutations complete.356352 permutations complete.356864 permutations complete.357376 permutations complete.357888 permutations complete.358400 permutations complete.358912 permutations complete.359424 permutations complete.359936 permutations complete.360448 permutations complete.360960 permutations complete.361472 permutations complete.361984 permutations complete.362496 permutations complete.363008 permutations complete.363520 permutations complete.364032 permutations complete.364544 permutations complete.365056 permutations complete.365568 permutations complete.366080 permutations complete.366592 permutations complete.367104 permutations complete.367616 permutations complete.368128 permutations complete.368640 permutations complete.369152 permutations complete.369664 permutations complete.370176 permutations complete.370688 permutations complete.371200 permutations complete.371712 permutations complete.372224 permutations complete.372736 permutations complete.373248 permutations complete.373760 permutations complete.374272 permutations complete.374784 permutations complete.375296 permutations complete.375808 permutations complete.376320 permutations complete.376832 permutations complete.377344 permutations complete.377856 permutations complete.378368 permutations complete.378880 permutations complete.379392 permutations complete.379904 permutations complete.380416 permutations complete.380928 permutations complete.381440 permutations complete.381952 permutations complete.382464 permutations complete.382976 permutations complete.383488 permutations complete.384000 permutations complete.384512 permutations complete.385024 permutations complete.385536 permutations complete.386048 permutations complete.386560 permutations complete.387072 permutations complete.387584 permutations complete.388096 permutations complete.388608 permutations complete.389120 permutations complete.389632 permutations complete.390144 permutations complete.390656 permutations complete.391168 permutations complete.391680 permutations complete.392192 permutations complete.392704 permutations complete.393216 permutations complete.393728 permutations complete.394240 permutations complete.394752 permutations complete.395264 permutations complete.395776 permutations complete.396288 permutations complete.396800 permutations complete.397312 permutations complete.397824 permutations complete.398336 permutations complete.398848 permutations complete.399360 permutations complete.399872 permutations complete.400384 permutations complete.400896 permutations complete.401408 permutations complete.401920 permutations complete.402432 permutations complete.402944 permutations complete.403456 permutations complete.403968 permutations complete.404480 permutations complete.404992 permutations complete.405504 permutations complete.406016 permutations complete.406528 permutations complete.407040 permutations complete.407552 permutations complete.408064 permutations complete.408576 permutations complete.409088 permutations complete.409600 permutations complete.410112 permutations complete.410624 permutations complete.411136 permutations complete.411648 permutations complete.412160 permutations complete.412672 permutations complete.413184 permutations complete.413696 permutations complete.414208 permutations complete.414720 permutations complete.415232 permutations complete.415744 permutations complete.416256 permutations complete.416768 permutations complete.417280 permutations complete.417792 permutations complete.418304 permutations complete.418816 permutations complete.419328 permutations complete.419840 permutations complete.420352 permutations complete.420864 permutations complete.421376 permutations complete.421888 permutations complete.422400 permutations complete.422912 permutations complete.423424 permutations complete.423936 permutations complete.424448 permutations complete.424960 permutations complete.425472 permutations complete.425984 permutations complete.426496 permutations complete.427008 permutations complete.427520 permutations complete.428032 permutations complete.428544 permutations complete.429056 permutations complete.429568 permutations complete.430080 permutations complete.430592 permutations complete.431104 permutations complete.431616 permutations complete.432128 permutations complete.432640 permutations complete.433152 permutations complete.433664 permutations complete.434176 permutations complete.434688 permutations complete.435200 permutations complete.435712 permutations complete.436224 permutations complete.436736 permutations complete.437248 permutations complete.437760 permutations complete.438272 permutations complete.438784 permutations complete.439296 permutations complete.439808 permutations complete.440320 permutations complete.440832 permutations complete.441344 permutations complete.441856 permutations complete.442368 permutations complete.442880 permutations complete.443392 permutations complete.443904 permutations complete.444416 permutations complete.444928 permutations complete.445440 permutations complete.445952 permutations complete.446464 permutations complete.446976 permutations complete.447488 permutations complete.448000 permutations complete.448512 permutations complete.449024 permutations complete.449536 permutations complete.450048 permutations complete.450560 permutations complete.451072 permutations complete.451584 permutations complete.452096 permutations complete.452608 permutations complete.453120 permutations complete.453632 permutations complete.454144 permutations complete.454656 permutations complete.455168 permutations complete.455680 permutations complete.456192 permutations complete.456704 permutations complete.457216 permutations complete.457728 permutations complete.458240 permutations complete.458752 permutations complete.459264 permutations complete.459776 permutations complete.460288 permutations complete.460800 permutations complete.461312 permutations complete.461824 permutations complete.462336 permutations complete.462848 permutations complete.463360 permutations complete.463872 permutations complete.464384 permutations complete.464896 permutations complete.465408 permutations complete.465920 permutations complete.466432 permutations complete.466944 permutations complete.467456 permutations complete.467968 permutations complete.468480 permutations complete.468992 permutations complete.469504 permutations complete.470016 permutations complete.470528 permutations complete.471040 permutations complete.471552 permutations complete.472064 permutations complete.472576 permutations complete.473088 permutations complete.473600 permutations complete.474112 permutations complete.474624 permutations complete.475136 permutations complete.475648 permutations complete.476160 permutations complete.476672 permutations complete.477184 permutations complete.477696 permutations complete.478208 permutations complete.478720 permutations complete.479232 permutations complete.479744 permutations complete.480256 permutations complete.480768 permutations complete.481280 permutations complete.481792 permutations complete.482304 permutations complete.482816 permutations complete.483328 permutations complete.483840 permutations complete.484352 permutations complete.484864 permutations complete.485376 permutations complete.485888 permutations complete.486400 permutations complete.486912 permutations complete.487424 permutations complete.487936 permutations complete.488448 permutations complete.488960 permutations complete.489472 permutations complete.489984 permutations complete.490496 permutations complete.491008 permutations complete.491520 permutations complete.492032 permutations complete.492544 permutations complete.493056 permutations complete.493568 permutations complete.494080 permutations complete.494592 permutations complete.495104 permutations complete.495616 permutations complete.496128 permutations complete.496640 permutations complete.497152 permutations complete.497664 permutations complete.498176 permutations complete.498688 permutations complete.499200 permutations complete.499712 permutations complete.500224 permutations complete.500736 permutations complete.501248 permutations complete.501760 permutations complete.502272 permutations complete.502784 permutations complete.503296 permutations complete.503808 permutations complete.504320 permutations complete.504832 permutations complete.505344 permutations complete.505856 permutations complete.506368 permutations complete.506880 permutations complete.507392 permutations complete.507904 permutations complete.508416 permutations complete.508928 permutations complete.509440 permutations complete.509952 permutations complete.510464 permutations complete.510976 permutations complete.511488 permutations complete.512000 permutations complete.512512 permutations complete.513024 permutations complete.513536 permutations complete.514048 permutations complete.514560 permutations complete.515072 permutations complete.515584 permutations complete.516096 permutations complete.516608 permutations complete.517120 permutations complete.517632 permutations complete.518144 permutations complete.518656 permutations complete.519168 permutations complete.519680 permutations complete.520192 permutations complete.520704 permutations complete.521216 permutations complete.521728 permutations complete.522240 permutations complete.522752 permutations complete.523264 permutations complete.523776 permutations complete.524288 permutations complete.524800 permutations complete.525312 permutations complete.525824 permutations complete.526336 permutations complete.526848 permutations complete.527360 permutations complete.527872 permutations complete.528384 permutations complete.528896 permutations complete.529408 permutations complete.529920 permutations complete.530432 permutations complete.530944 permutations complete.531456 permutations complete.531968 permutations complete.532480 permutations complete.532992 permutations complete.533504 permutations complete.534016 permutations complete.534528 permutations complete.535040 permutations complete.535552 permutations complete.536064 permutations complete.536576 permutations complete.537088 permutations complete.537600 permutations complete.538112 permutations complete.538624 permutations complete.539136 permutations complete.539648 permutations complete.540160 permutations complete.540672 permutations complete.541184 permutations complete.541696 permutations complete.542208 permutations complete.542720 permutations complete.543232 permutations complete.543744 permutations complete.544256 permutations complete.544768 permutations complete.545280 permutations complete.545792 permutations complete.546304 permutations complete.546816 permutations complete.547328 permutations complete.547840 permutations complete.548352 permutations complete.548864 permutations complete.549376 permutations complete.549888 permutations complete.550400 permutations complete.550912 permutations complete.551424 permutations complete.551936 permutations complete.552448 permutations complete.552960 permutations complete.553472 permutations complete.553984 permutations complete.554496 permutations complete.555008 permutations complete.555520 permutations complete.556032 permutations complete.556544 permutations complete.557056 permutations complete.557568 permutations complete.558080 permutations complete.558592 permutations complete.559104 permutations complete.559616 permutations complete.560128 permutations complete.560640 permutations complete.561152 permutations complete.561664 permutations complete.562176 permutations complete.562688 permutations complete.563200 permutations complete.563712 permutations complete.564224 permutations complete.564736 permutations complete.565248 permutations complete.565760 permutations complete.566272 permutations complete.566784 permutations complete.567296 permutations complete.567808 permutations complete.568320 permutations complete.568832 permutations complete.569344 permutations complete.569856 permutations complete.570368 permutations complete.570880 permutations complete.571392 permutations complete.571904 permutations complete.572416 permutations complete.572928 permutations complete.573440 permutations complete.573952 permutations complete.574464 permutations complete.574976 permutations complete.575488 permutations complete.576000 permutations complete.576512 permutations complete.577024 permutations complete.577536 permutations complete.578048 permutations complete.578560 permutations complete.579072 permutations complete.579584 permutations complete.580096 permutations complete.580608 permutations complete.581120 permutations complete.581632 permutations complete.582144 permutations complete.582656 permutations complete.583168 permutations complete.583680 permutations complete.584192 permutations complete.584704 permutations complete.585216 permutations complete.585728 permutations complete.586240 permutations complete.586752 permutations complete.587264 permutations complete.587776 permutations complete.588288 permutations complete.588800 permutations complete.589312 permutations complete.589824 permutations complete.590336 permutations complete.590848 permutations complete.591360 permutations complete.591872 permutations complete.592384 permutations complete.592896 permutations complete.593408 permutations complete.593920 permutations complete.594432 permutations complete.594944 permutations complete.595456 permutations complete.595968 permutations complete.596480 permutations complete.596992 permutations complete.597504 permutations complete.598016 permutations complete.598528 permutations complete.599040 permutations complete.599552 permutations complete.600064 permutations complete.600576 permutations complete.601088 permutations complete.601600 permutations complete.602112 permutations complete.602624 permutations complete.603136 permutations complete.603648 permutations complete.604160 permutations complete.604672 permutations complete.605184 permutations complete.605696 permutations complete.606208 permutations complete.606720 permutations complete.607232 permutations complete.607744 permutations complete.608256 permutations complete.608768 permutations complete.609280 permutations complete.609792 permutations complete.610304 permutations complete.610816 permutations complete.611328 permutations complete.611840 permutations complete.612352 permutations complete.612864 permutations complete.613376 permutations complete.613888 permutations complete.614400 permutations complete.614912 permutations complete.615424 permutations complete.615936 permutations complete.616448 permutations complete.616960 permutations complete.617472 permutations complete.617984 permutations complete.618496 permutations complete.619008 permutations complete.619520 permutations complete.620032 permutations complete.620544 permutations complete.621056 permutations complete.621568 permutations complete.622080 permutations complete.622592 permutations complete.623104 permutations complete.623616 permutations complete.624128 permutations complete.624640 permutations complete.625152 permutations complete.625664 permutations complete.626176 permutations complete.626688 permutations complete.627200 permutations complete.627712 permutations complete.628224 permutations complete.628736 permutations complete.629248 permutations complete.629760 permutations complete.630272 permutations complete.630784 permutations complete.631296 permutations complete.631808 permutations complete.632320 permutations complete.632832 permutations complete.633344 permutations complete.633856 permutations complete.634368 permutations complete.634880 permutations complete.635392 permutations complete.635904 permutations complete.636416 permutations complete.636928 permutations complete.637440 permutations complete.637952 permutations complete.638464 permutations complete.638976 permutations complete.639488 permutations complete.640000 permutations complete.640512 permutations complete.641024 permutations complete.641536 permutations complete.642048 permutations complete.642560 permutations complete.643072 permutations complete.643584 permutations complete.644096 permutations complete.644608 permutations complete.645120 permutations complete.645632 permutations complete.646144 permutations complete.646656 permutations complete.647168 permutations complete.647680 permutations complete.648192 permutations complete.648704 permutations complete.649216 permutations complete.649728 permutations complete.650240 permutations complete.650752 permutations complete.651264 permutations complete.651776 permutations complete.652288 permutations complete.652800 permutations complete.653312 permutations complete.653824 permutations complete.654336 permutations complete.654848 permutations complete.655360 permutations complete.655872 permutations complete.656384 permutations complete.656896 permutations complete.657408 permutations complete.657920 permutations complete.658432 permutations complete.658944 permutations complete.659456 permutations complete.659968 permutations complete.660480 permutations complete.660992 permutations complete.661504 permutations complete.662016 permutations complete.662528 permutations complete.663040 permutations complete.663552 permutations complete.664064 permutations complete.664576 permutations complete.665088 permutations complete.665600 permutations complete.666112 permutations complete.666624 permutations complete.667136 permutations complete.667648 permutations complete.668160 permutations complete.668672 permutations complete.669184 permutations complete.669696 permutations complete.670208 permutations complete.670720 permutations complete.671232 permutations complete.671744 permutations complete.672256 permutations complete.672768 permutations complete.673280 permutations complete.673792 permutations complete.674304 permutations complete.674816 permutations complete.675328 permutations complete.675840 permutations complete.676352 permutations complete.676864 permutations complete.677376 permutations complete.677888 permutations complete.678400 permutations complete.678912 permutations complete.679424 permutations complete.679936 permutations complete.680448 permutations complete.680960 permutations complete.681472 permutations complete.681984 permutations complete.682496 permutations complete.683008 permutations complete.683520 permutations complete.684032 permutations complete.684544 permutations complete.685056 permutations complete.685568 permutations complete.686080 permutations complete.686592 permutations complete.687104 permutations complete.687616 permutations complete.688128 permutations complete.688640 permutations complete.689152 permutations complete.689664 permutations complete.690176 permutations complete.690688 permutations complete.691200 permutations complete.691712 permutations complete.692224 permutations complete.692736 permutations complete.693248 permutations complete.693760 permutations complete.694272 permutations complete.694784 permutations complete.695296 permutations complete.695808 permutations complete.696320 permutations complete.696832 permutations complete.697344 permutations complete.697856 permutations complete.698368 permutations complete.698880 permutations complete.699392 permutations complete.699904 permutations complete.700416 permutations complete.700928 permutations complete.701440 permutations complete.701952 permutations complete.702464 permutations complete.702976 permutations complete.703488 permutations complete.704000 permutations complete.704512 permutations complete.705024 permutations complete.705536 permutations complete.706048 permutations complete.706560 permutations complete.707072 permutations complete.707584 permutations complete.708096 permutations complete.708608 permutations complete.709120 permutations complete.709632 permutations complete.710144 permutations complete.710656 permutations complete.711168 permutations complete.711680 permutations complete.712192 permutations complete.712704 permutations complete.713216 permutations complete.713728 permutations complete.714240 permutations complete.714752 permutations complete.715264 permutations complete.715776 permutations complete.716288 permutations complete.716800 permutations complete.717312 permutations complete.717824 permutations complete.718336 permutations complete.718848 permutations complete.719360 permutations complete.719872 permutations complete.720384 permutations complete.720896 permutations complete.721408 permutations complete.721920 permutations complete.722432 permutations complete.722944 permutations complete.723456 permutations complete.723968 permutations complete.724480 permutations complete.724992 permutations complete.725504 permutations complete.726016 permutations complete.726528 permutations complete.727040 permutations complete.727552 permutations complete.728064 permutations complete.728576 permutations complete.729088 permutations complete.729600 permutations complete.730112 permutations complete.730624 permutations complete.731136 permutations complete.731648 permutations complete.732160 permutations complete.732672 permutations complete.733184 permutations complete.733696 permutations complete.734208 permutations complete.734720 permutations complete.735232 permutations complete.735744 permutations complete.736256 permutations complete.736768 permutations complete.737280 permutations complete.737792 permutations complete.738304 permutations complete.738816 permutations complete.739328 permutations complete.739840 permutations complete.740352 permutations complete.740864 permutations complete.741376 permutations complete.741888 permutations complete.742400 permutations complete.742912 permutations complete.743424 permutations complete.743936 permutations complete.744448 permutations complete.744960 permutations complete.745472 permutations complete.745984 permutations complete.746496 permutations complete.747008 permutations complete.747520 permutations complete.748032 permutations complete.748544 permutations complete.749056 permutations complete.749568 permutations complete.750080 permutations complete.750592 permutations complete.751104 permutations complete.751616 permutations complete.752128 permutations complete.752640 permutations complete.753152 permutations complete.753664 permutations complete.754176 permutations complete.754688 permutations complete.755200 permutations complete.755712 permutations complete.756224 permutations complete.756736 permutations complete.757248 permutations complete.757760 permutations complete.758272 permutations complete.758784 permutations complete.759296 permutations complete.759808 permutations complete.760320 permutations complete.760832 permutations complete.761344 permutations complete.761856 permutations complete.762368 permutations complete.762880 permutations complete.763392 permutations complete.763904 permutations complete.764416 permutations complete.764928 permutations complete.765440 permutations complete.765952 permutations complete.766464 permutations complete.766976 permutations complete.767488 permutations complete.768000 permutations complete.768512 permutations complete.769024 permutations complete.769536 permutations complete.770048 permutations complete.770560 permutations complete.771072 permutations complete.771584 permutations complete.772096 permutations complete.772608 permutations complete.773120 permutations complete.773632 permutations complete.774144 permutations complete.774656 permutations complete.775168 permutations complete.775680 permutations complete.776192 permutations complete.776704 permutations complete.777216 permutations complete.777728 permutations complete.778240 permutations complete.778752 permutations complete.779264 permutations complete.779776 permutations complete.780288 permutations complete.780800 permutations complete.781312 permutations complete.781824 permutations complete.782336 permutations complete.782848 permutations complete.783360 permutations complete.783872 permutations complete.784384 permutations complete.784896 permutations complete.785408 permutations complete.785920 permutations complete.786432 permutations complete.786944 permutations complete.787456 permutations complete.787968 permutations complete.788480 permutations complete.788992 permutations complete.789504 permutations complete.790016 permutations complete.790528 permutations complete.791040 permutations complete.791552 permutations complete.792064 permutations complete.792576 permutations complete.793088 permutations complete.793600 permutations complete.794112 permutations complete.794624 permutations complete.795136 permutations complete.795648 permutations complete.796160 permutations complete.796672 permutations complete.797184 permutations complete.797696 permutations complete.798208 permutations complete.798720 permutations complete.799232 permutations complete.799744 permutations complete.800256 permutations complete.800768 permutations complete.801280 permutations complete.801792 permutations complete.802304 permutations complete.802816 permutations complete.803328 permutations complete.803840 permutations complete.804352 permutations complete.804864 permutations complete.805376 permutations complete.805888 permutations complete.806400 permutations complete.806912 permutations complete.807424 permutations complete.807936 permutations complete.808448 permutations complete.808960 permutations complete.809472 permutations complete.809984 permutations complete.810496 permutations complete.811008 permutations complete.811520 permutations complete.812032 permutations complete.812544 permutations complete.813056 permutations complete.813568 permutations complete.814080 permutations complete.814592 permutations complete.815104 permutations complete.815616 permutations complete.816128 permutations complete.816640 permutations complete.817152 permutations complete.817664 permutations complete.818176 permutations complete.818688 permutations complete.819200 permutations complete.819712 permutations complete.820224 permutations complete.820736 permutations complete.821248 permutations complete.821760 permutations complete.822272 permutations complete.822784 permutations complete.823296 permutations complete.823808 permutations complete.824320 permutations complete.824832 permutations complete.825344 permutations complete.825856 permutations complete.826368 permutations complete.826880 permutations complete.827392 permutations complete.827904 permutations complete.828416 permutations complete.828928 permutations complete.829440 permutations complete.829952 permutations complete.830464 permutations complete.830976 permutations complete.831488 permutations complete.832000 permutations complete.832512 permutations complete.833024 permutations complete.833536 permutations complete.834048 permutations complete.834560 permutations complete.835072 permutations complete.835584 permutations complete.836096 permutations complete.836608 permutations complete.837120 permutations complete.837632 permutations complete.838144 permutations complete.838656 permutations complete.839168 permutations complete.839680 permutations complete.840192 permutations complete.840704 permutations complete.841216 permutations complete.841728 permutations complete.842240 permutations complete.842752 permutations complete.843264 permutations complete.843776 permutations complete.844288 permutations complete.844800 permutations complete.845312 permutations complete.845824 permutations complete.846336 permutations complete.846848 permutations complete.847360 permutations complete.847872 permutations complete.848384 permutations complete.848896 permutations complete.849408 permutations complete.849920 permutations complete.850432 permutations complete.850944 permutations complete.851456 permutations complete.851968 permutations complete.852480 permutations complete.852992 permutations complete.853504 permutations complete.854016 permutations complete.854528 permutations complete.855040 permutations complete.855552 permutations complete.856064 permutations complete.856576 permutations complete.857088 permutations complete.857600 permutations complete.858112 permutations complete.858624 permutations complete.859136 permutations complete.859648 permutations complete.860160 permutations complete.860672 permutations complete.861184 permutations complete.861696 permutations complete.862208 permutations complete.862720 permutations complete.863232 permutations complete.863744 permutations complete.864256 permutations complete.864768 permutations complete.865280 permutations complete.865792 permutations complete.866304 permutations complete.866816 permutations complete.867328 permutations complete.867840 permutations complete.868352 permutations complete.868864 permutations complete.869376 permutations complete.869888 permutations complete.870400 permutations complete.870912 permutations complete.871424 permutations complete.871936 permutations complete.872448 permutations complete.872960 permutations complete.873472 permutations complete.873984 permutations complete.874496 permutations complete.875008 permutations complete.875520 permutations complete.876032 permutations complete.876544 permutations complete.877056 permutations complete.877568 permutations complete.878080 permutations complete.878592 permutations complete.879104 permutations complete.879616 permutations complete.880128 permutations complete.880640 permutations complete.881152 permutations complete.881664 permutations complete.882176 permutations complete.882688 permutations complete.883200 permutations complete.883712 permutations complete.884224 permutations complete.884736 permutations complete.885248 permutations complete.885760 permutations complete.886272 permutations complete.886784 permutations complete.887296 permutations complete.887808 permutations complete.888320 permutations complete.888832 permutations complete.889344 permutations complete.889856 permutations complete.890368 permutations complete.890880 permutations complete.891392 permutations complete.891904 permutations complete.892416 permutations complete.892928 permutations complete.893440 permutations complete.893952 permutations complete.894464 permutations complete.894976 permutations complete.895488 permutations complete.896000 permutations complete.896512 permutations complete.897024 permutations complete.897536 permutations complete.898048 permutations complete.898560 permutations complete.899072 permutations complete.899584 permutations complete.900096 permutations complete.900608 permutations complete.901120 permutations complete.901632 permutations complete.902144 permutations complete.902656 permutations complete.903168 permutations complete.903680 permutations complete.904192 permutations complete.904704 permutations complete.905216 permutations complete.905728 permutations complete.906240 permutations complete.906752 permutations complete.907264 permutations complete.907776 permutations complete.908288 permutations complete.908800 permutations complete.909312 permutations complete.909824 permutations complete.910336 permutations complete.910848 permutations complete.911360 permutations complete.911872 permutations complete.912384 permutations complete.912896 permutations complete.913408 permutations complete.913920 permutations complete.914432 permutations complete.914944 permutations complete.915456 permutations complete.915968 permutations complete.916480 permutations complete.916992 permutations complete.917504 permutations complete.918016 permutations complete.918528 permutations complete.919040 permutations complete.919552 permutations complete.920064 permutations complete.920576 permutations complete.921088 permutations complete.921600 permutations complete.922112 permutations complete.922624 permutations complete.923136 permutations complete.923648 permutations complete.924160 permutations complete.924672 permutations complete.925184 permutations complete.925696 permutations complete.926208 permutations complete.926720 permutations complete.927232 permutations complete.927744 permutations complete.928256 permutations complete.928768 permutations complete.929280 permutations complete.929792 permutations complete.930304 permutations complete.930816 permutations complete.931328 permutations complete.931840 permutations complete.932352 permutations complete.932864 permutations complete.933376 permutations complete.933888 permutations complete.934400 permutations complete.934912 permutations complete.935424 permutations complete.935936 permutations complete.936448 permutations complete.936960 permutations complete.937472 permutations complete.937984 permutations complete.938496 permutations complete.939008 permutations complete.939520 permutations complete.940032 permutations complete.940544 permutations complete.941056 permutations complete.941568 permutations complete.942080 permutations complete.942592 permutations complete.943104 permutations complete.943616 permutations complete.944128 permutations complete.944640 permutations complete.945152 permutations complete.945664 permutations complete.946176 permutations complete.946688 permutations complete.947200 permutations complete.947712 permutations complete.948224 permutations complete.948736 permutations complete.949248 permutations complete.949760 permutations complete.950272 permutations complete.950784 permutations complete.951296 permutations complete.951808 permutations complete.952320 permutations complete.952832 permutations complete.953344 permutations complete.953856 permutations complete.954368 permutations complete.954880 permutations complete.955392 permutations complete.955904 permutations complete.956416 permutations complete.956928 permutations complete.957440 permutations complete.957952 permutations complete.958464 permutations complete.958976 permutations complete.959488 permutations complete.960000 permutations complete.960512 permutations complete.961024 permutations complete.961536 permutations complete.962048 permutations complete.962560 permutations complete.963072 permutations complete.963584 permutations complete.964096 permutations complete.964608 permutations complete.965120 permutations complete.965632 permutations complete.966144 permutations complete.966656 permutations complete.967168 permutations complete.967680 permutations complete.968192 permutations complete.968704 permutations complete.969216 permutations complete.969728 permutations complete.970240 permutations complete.970752 permutations complete.971264 permutations complete.971776 permutations complete.972288 permutations complete.972800 permutations complete.973312 permutations complete.973824 permutations complete.974336 permutations complete.974848 permutations complete.975360 permutations complete.975872 permutations complete.976384 permutations complete.976896 permutations complete.977408 permutations complete.977920 permutations complete.978432 permutations complete.978944 permutations complete.979456 permutations complete.979968 permutations complete.980480 permutations complete.980992 permutations complete.981504 permutations complete.982016 permutations complete.982528 permutations complete.983040 permutations complete.983552 permutations complete.984064 permutations complete.984576 permutations complete.985088 permutations complete.985600 permutations complete.986112 permutations complete.986624 permutations complete.987136 permutations complete.987648 permutations complete.988160 permutations complete.988672 permutations complete.989184 permutations complete.989696 permutations complete.990208 permutations complete.990720 permutations complete.991232 permutations complete.991744 permutations complete.992256 permutations complete.992768 permutations complete.993280 permutations complete.993792 permutations complete.994304 permutations complete.994816 permutations complete.995328 permutations complete.995840 permutations complete.996352 permutations complete.996864 permutations complete.997376 permutations complete.997888 permutations complete.998400 permutations complete.998912 permutations complete.999424 permutations complete.999936 permutations complete.1000000 (adaptive) permutations complete.
    ## Permutation test report written to quant2.qassoc.perm .

``` r
file10 <- read.table("quant2.qassoc.perm", header = TRUE, sep = "\t")
kable(head(file10))
```

| CHR………SNP………EMP1………..NP |
| :---------------------- |
| 1 rs6681049 0.3676 67   |
| 1 rs4074137 0.5926 26   |
| 1 rs7540009 NA NA       |
| 1 rs1891905 0.875 7     |
| 1 rs9729550 0.05092 706 |
| 1 rs3813196 0.7059 16   |

``` bash
plink --bfile hapmap1 --assoc --pheno qt.phe --mperm 1000 --within str1.cluster2 --out quant3
plink --bfile hapmap1 --pheno qt.phe --gxe --covar pop.phe --snp rs2222162 --out quant3
```

    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to quant3.log.
    ## Options in effect:
    ##   --assoc
    ##   --bfile hapmap1
    ##   --mperm 1000
    ##   --out quant3
    ##   --pheno qt.phe
    ##   --within str1.cluster2
    ## 
    ## Note: --mperm flag deprecated.  Use e.g. '--model mperm=[value]'.
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values present after --pheno.
    ## --within: 45 clusters loaded, covering a total of 89 people.
    ## Using up to 47 threads (change this with --threads).
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## Total genotyping rate is 0.99441.
    ## 83534 variants and 89 people pass filters and QC.
    ## Phenotype data is quantitative.
    ## Writing QT --assoc report to quant3.qassoc ... 0%1%2%3%4%5%6%8%9%10%11%12%13%15%16%17%18%19%20%21%23%25%26%27%28%29%31%32%33%34%35%36%38%39%40%41%42%43%44%46%47%48%49%50%51%52%54%55%56%57%59%60%61%62%63%64%65%66%67%68%69%70%72%73%74%76%77%78%79%80%81%82%83%84%85%87%88%89%90%91%92%93%94%95%96%97%98%99%done.
    ## 512 permutations complete.1000 max(T) permutations complete.
    ## Permutation test report written to quant3.qassoc.mperm .
    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Logging to quant3.log.
    ## Options in effect:
    ##   --bfile hapmap1
    ##   --covar pop.phe
    ##   --gxe
    ##   --out quant3
    ##   --pheno qt.phe
    ##   --snp rs2222162
    ## 
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## 1 out of 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values present after --pheno.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## --covar: 1 case/control covariate loaded for --gxe.
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## 1 variant and 89 people pass filters and QC.
    ## Phenotype data is quantitative.
    ## Writing --gxe report to quant3.qassoc.gxe ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.

``` bash
#Extracting a SNP of interest
plink --bfile hapmap1 --snp rs2222162 --recodeAD --out rec_snp1
d <- read.table("rec_snp1.recode.raw" , header=T)
summary(glm(PHENOTYPE-1 ~ rs2222162_A, data=d, family="binomial"))
```

    ## PLINK v1.90b4.9 64-bit (13 Oct 2017)           www.cog-genomics.org/plink/1.9/
    ## (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3
    ## Note: --recodeAD flag deprecated.  Use 'recode AD ...'.
    ## Logging to rec_snp1.log.
    ## Options in effect:
    ##   --bfile hapmap1
    ##   --out rec_snp1
    ##   --recode AD
    ##   --snp rs2222162
    ## 
    ## 193230 MB RAM detected; reserving 96615 MB for main workspace.
    ## 1 out of 83534 variants loaded from .bim file.
    ## 89 people (89 males, 0 females) loaded from .fam.
    ## 89 phenotype values loaded from .fam.
    ## Using 1 thread (no multithreaded calculations invoked).
    ## Before main variant filters, 89 founders and 0 nonfounders present.
    ## Calculating allele frequencies... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99% done.
    ## 1 variant and 89 people pass filters and QC.
    ## Among remaining phenotypes, 44 are cases and 45 are controls.
    ## --recode AD to rec_snp1.raw ... 0%1%2%3%4%5%6%7%8%9%10%11%12%13%14%15%16%17%18%19%20%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40%41%42%43%44%45%46%47%48%49%50%51%52%53%54%55%56%57%58%59%60%61%62%63%64%65%66%67%68%69%70%71%72%73%74%75%76%77%78%79%80%81%82%83%84%85%86%87%88%89%90%91%92%93%94%95%96%97%98%99%done.
    ## bash: -c: line 2: syntax error near unexpected token `('
    ## bash: -c: line 2: `d <- read.table("rec_snp1.recode.raw" , header=T)'
