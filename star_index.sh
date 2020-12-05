
#!/bin/bash 


# Define reference directy 
ref=reference_GENCODE

# Define output directory
index=star_index


# Define reference file names
genome=*genome.fa
gtf=*.gtf

mkdir $ref/$index  

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $ref/$index --genomeFastaFiles $ref/$genome --sjdbGTFfile $ref/$gtf 
