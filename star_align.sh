
#!/bin/bash

# Define directory names 
outdir=star_output   # directory storing output files
indir=rawdata        # input directory (fastq files)
refdir=/home/mira/Documents/programming/Bioinformatics/SEQC/reference_GENCODE     # reference directory (absolute path needed!)
indexdir=star_index    # index directory
genome=*genome.fa      # reference file
gtf=*.gtf              # GTF file


samples=(A{1..3} B{1..3})



for read in ${samples[*]} 

do 
    STAR --runThreadN 8 --runMode alignReads --genomeDir $refdir/$indexdir --sjdbGTFfile $refdir/$gtf -sjdbOverhang 100 --readFilesIn $indir/${read}_1.fastq $indir/${read}_2.fastq --outFileNamePrefix $read --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --outSAMunmapped None --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --twopassMode Basic --chimOutType Junctions

done


# --readFilesIn: For paired-end reads, use comma separated list for read1, followed by space, followed by comma separated list for read2
# --readFilesCommand: required when the input files are .gzip format (e.g. --readFilesCommand zcat, --readFilesCommand gunzip -c, or --readFilesCommand bunzip2)
# Note: the output files are generated in the current directory
