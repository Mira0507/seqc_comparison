#!/bin/bash

cd hisat2_output

input=(A{1..3} B{1..3})  

for f in ${input[*]}

do
    samtools sort $f.bam -o $f.sorted.bam

done 


cd ..

# Delete SAM files after samtools run
