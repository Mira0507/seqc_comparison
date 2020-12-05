#!/bin/bash

cd hisat2_output

input=(B{1..3})  

for f in ${input[*]}

do
    samtools view -bS $f.sam > $f.bam 
done 


cd ..
