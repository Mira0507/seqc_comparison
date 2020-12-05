#!/bin/bash

cd reference_GENCODE

indir=salmon_index

# Assign location of transcripts and genome reference files 
transcripts=*transcripts.fa 
genome=*genome.fa

# Concatenate them and save as.gentrome.fa
cat $transcripts $genome > $indir/gentrome.fa
