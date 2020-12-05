#!/bin/bash 

# Define directory and sample names
refdir=reference_GENCODE/hisat2_index/index   # reference directory
samples=(A{1..3} B{1..3})                     # sample names
outdir=hisat2_output                          # output directory
indir=rawdata                                 # input directory

mkdir $outdir 


for read in ${samples[*]} 

do 

    hisat2 -q -p 8 --seed 23 -x $refdir -1 $indir/${read}_1.fastq -2 $indir/${read}_2.fastq -S $outdir/$read.sam 

done 

cd ..
