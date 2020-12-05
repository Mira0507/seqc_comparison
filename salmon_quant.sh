#!/bin/bash

# Define name of input/output directory
in=rawdata
out=salmon_output


# Define index file directory
ind=reference_GENCODE/salmon_index/gencode_index

# Define file names 
samples=(A{1..3} B{1..3})

mkdir $out

cd $in

for read in ${samples[*]}
do
    salmon quant -i ../$ind -l A --gcBias --seqBias -1 ${read}_1.fastq -2 ${read}_2.fastq -p 8 --validateMappings -o ../$out/${read}.salmon_quant
done

cd ..
