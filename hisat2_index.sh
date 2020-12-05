#!/bin/bash

# Define directory path or file 
refdir=reference_GENCODE    # reference directory
genome=*genome.fa           # reference genome
outdir=hisat2_index         # ouput directory storing index files

mkdir $refdir/$outdir

hisat2-build -f -o 4 -p 8 --seed 67 $refdir/$genome $refdir/$outdir/index 

cd ..

# Basic command line: hisat2-build [options] <reference_in> <ht2_base>
# -f ~ --seed: options
# "reference_GENCODE/*genome.fa": <reference_in> 
# "index": <ht2_base> 
