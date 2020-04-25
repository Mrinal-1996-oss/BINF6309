#!/usr/bin/env bash
outDir="quant/"
sample="/scratch/SampleDataFiles/Paired/"
R1=".R1.paired.fastq"
R2=".R2.paired.fastq"

function alignAll {
    for fastqs in $sample*Aip*$R1
    do
        pathRemove="${fastqs/$sample/}"
        sampleName="${pathRemove/$R1/}"
        
        salmon quant -l IU \
            -1 $sample$sampleName$R1 \
            -2 $sample$sampleName$R2 \
            -i AipIndex \
            --validateMappings \
            -o $outDir$sampleName
    done    
}

alignAll 1>alignAll.log 2>alignAll.err 
