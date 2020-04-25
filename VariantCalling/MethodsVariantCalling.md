## Variant Calling

The first step towards Variant Calling is alignment.For alignment, the
tools that are being used are BWA(Burrow Wheel
Alligner),Speedseq,Bowtie2,BFAST,Illumina Isaak WGS,etc.In this code, we
are basically using the BWA aligner. BWA is a software package for
mapping low-divergent sequences against a large reference genome, such
as the human genome. It consists of three algorithms: BWA-backtrack,
BWA-SW and BWA-MEM.The first algorithm is designed for Illumina sequence
reads upto 100bp,while the rest two for longer sequences ranged from
70bp to 1Mbp.

For the second step of getting the reads,we use the fastq-dump function
which allows to split files SRR6808334 and gets the specific .log and
.err function.

For the third step , Somatic variant callingis performed.Somatic variant
calling includes steps like SNP , Indels, CNVs and SVs.After getting the
reads from SRR6808334, we are trying to trim thte reads using
Trimmomatic.Trimmomatic helps to trim the read files and distinguish
between the paired and unpaired fastq files. Trimmomatic as a more
flexible and efficient preprocessing tool, which could correctly handle
paired-end data. The value of NGS read preprocessing is demonstrated for
both reference-based and reference-free tasks. Trimmomatic is shown to
produce output that is at least competitive with, and in many cases
superior to, that produced by other tools, in all scenarios tested .

For the fourth step, the genomes are subsequently indexed. Germline
variant callingis performed.A germline variant caller generally has a
ploidy-based genotyping algorithm built in to part of the
algorithm/pipeline.Germline variant caller uses SNPs,Indels,CNVs and SVs
as well.

For the fifth step, reads are aligned using BWA mem and are then sorted
by the same algorithm .The reads are then indexed and a VCF file is
produced using DeepVariantThe variants are filtered and annotated using
the tools like Snpedd,snpsift,annovar,Frequency db,gnomad, dbsnp and
RareVariantVis.

Finally, the data is visualised using tools that include GV ,UCSC and
RareVariantVis.The plots are available in heatmaps,circos plots,scatter
plots, log2 plots,etc

## References

<div id="refs" class="references">

<div id="ref-Bolger2014">

Bolger, Anthony M., Marc Lohse, and Bjoern Usadel. 2014. “Trimmomatic: A
Flexible Trimmer for Illumina Sequence Data.” *Bioinformatics (Oxford,
England)* 30 (15): 2114–20.
<https://doi.org/10.1093/bioinformatics/btu170>.

</div>

<div id="ref-Li2009">

Li, Heng, and Richard Durbin. 2009. “Fast and Accurate Short Read
Alignment with Burrows-Wheeler Transform.” *Bioinformatics (Oxford,
England)* 25 (14): 1754–60.
<https://doi.org/10.1093/bioinformatics/btp324>.

</div>

<div id="ref-McKenna2010">

McKenna, Aaron, Matthew Hanna, Eric Banks, Andrey Sivachenko, Kristian
Cibulskis, Andrew Kernytsky, Kiran Garimella, et al. 2010. “The Genome
Analysis Toolkit: A MapReduce Framework for Analyzing Next-Generation
DNA Sequencing Data.” *Genome Research* 20 (9): 1297–1303.
<https://doi.org/10.1101/gr.107524.110>.

</div>

<div id="ref-Poplin2018">

Poplin, Ryan, Pi-Chuan Chang, David Alexander, Scott Schwartz, Thomas
Colthurst, Alexander Ku, Dan Newburger, et al. 2018. “A Universal SNP
and Small-Indel Variant Caller Using Deep Neural Networks.” *Nature
Biotechnology* 36 (10): 983–87. <https://doi.org/10.1038/nbt.4235>.

</div>

</div>
