#!/bin/bash

# sra toolkit doc: https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=prefetch

# create a directory to store rawdata

cd rawdata



for x in SRR9500{78..85}
do 
    fastq-dump --split-files $x/$x.sra 
done

cd ..
