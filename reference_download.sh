#!/bin/bash

todir=reference_GENCODE

mkdir $todir

# Assign location of transcriptome file
transcriptome=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.transcripts.fa.gz
genome=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.primary_assembly.genome.fa.gz
gtf=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.primary_assembly.annotation.gtf.gz


# Download human transcriptome and genome files 
wget -c $transcriptome -O $todir/gencode.v36.transcripts.fa.gz
wget -c $genome -O $todir/GRCh38.v36.genome.fa.gz
wget -c $gtf -O $todir/gencode.v36.primary_assembly.annotation.gtf.gz

# Unzip the reference files manually!!!
