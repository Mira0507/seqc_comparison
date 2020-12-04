## Analysis of SEQC datasets 

### 1. Rawdata 

- GEO number: [GSE49712](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49712)
- SRA number: [SRP028705](https://www.ncbi.nlm.nih.gov/sra?term=SRP028705)
- Reference paper: [Genome Biol. 2013;14(9):R95. doi: 10.1186/gb-2013-14-9-r95.](https://pubmed.ncbi.nlm.nih.gov/24020486)
- Sample preparation: The Sequencing Quality Control Consortium generated two datasets from two reference RNA samples in order to evaluate transcriptome profiling by next-generation sequencing technology. Each sample contains one of the reference RNA source and a set of synthetic RNAs from the External RNA Control Consortium (ERCC) at known concentrations. Group A contains 5 replicates of the Strategene Universal Human Reference RNA (UHRR), which is composed of total RNA from 10 human cell lines, with 2% by volume of ERCC mix 1. Group B includes 5 replicate samples of the Ambion Human Brain Reference RNA (HBRR) with 2% by volume of ERCC mix 2. The ERCC spike-in control is a mixture of 92 synthetic polyadenylated oligonucleotides of 250-2000 nucleotides long that are meant to resemble human transcripts.      
- Raw data was downloaded with [SRA toolkit](https://github.com/Mira0507/using_SRA/blob/master/README.md)

#### 1-1. Downloading SRA files 

- rawdata_download.sh

```bash
#!/bin/bash


# Numbers of SRA data:
# SRR950078: A1
# SRR950079: B2
# SRR950080: A2
# SRR950081: B2
# SRR950082: A3
# SRR950083: B3


cd rawdata


for x in {78..83}
do
    prefetch SRR9500$x

done

cd .. 
```

#### 1-2. Converting to FASTQ files 

- sra_fastq.sh

```bash
#!/bin/bash

# sra toolkit doc: https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=prefetch

# create a directory to store rawdata

cd rawdata



for x in SRR9500{78..83}
do 
    fastq-dump --split-files $x/$x.sra 
done

cd ..
```


### 2. Conda environment

- Conda doc: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html

#### 2-1. Conde env for alignment/mapping

- HISAT2: http://daehwankimlab.github.io/hisat2
- Samtools: http://www.htslib.org/doc/#manual-pages
- salmon: https://salmon.readthedocs.io/en/latest
- bedtools: https://bedtools.readthedocs.io/en/latest
- gawk: https://www.gnu.org/software/gawk/manual/gawk.html
- STAR: https://github.com/alexdobin/STAR

```yml

name: mapping 
channels:
  - conda-forge
  - bioconda 
  - defaults 
dependencies:
  - hisat2=2.2.1
  - samtools=1.11
  - salmon=1.4.0
  - bedtools=2.29.2 
  - gawk=5.1.0 
  - star=2.7.6a

```

#### 2-2. Conda env for counting and differential expression (DE) analysis

- DESeq2: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
- Tximport: http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
- Rsubread: https://bioconductor.org/packages/release/bioc/html/Rsubread.html

```yml

name: r
channels:
  - conda-forge
  - bioconda 
  - defaults 
dependencies:
  - r-base =4.0.3
  - r-tidyverse
  - r-data.table
  - r-ggplot2
  - r-markdown
  - r-upsetr
  - r-ggrepel
  - r-ggplotify
  - r-gridextra
  - r-pheatmap
  - bioconductor-deseq2=1.30.0
  - bioconductor-annotationhub=2.22.0
  - bioconductor-tximport=1.18.0
  - bioconductor-rsubread=2.4.0
  - bioconductor-seqc=1.24.0
```

### 3. Reference files

- Version: [GENCODE](https://www.gencodegenes.org/human) GRCh38 (hg38) release 36 (v36)
- Genome: Genome sequence, primary assembly (GRCh38)
- Transcriptome: Transcript sequences
- GTF: Transcript sequences (PRI)
- reference_download.sh

```bash
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
```

### 4. Salmon mapping 

#### 4-1. Indexing

- **Creating decoys.txt file**: salmon_decoy_txt.sh

```bash
#!/bin/bash

indir=salmon_index

# Move to the indexing directory
cd reference_GENCODE

# Save reference names into decoys.txt
grep "^>" < *.genome.fa | cut -d " " -f 1 > $indir/decoys.txt

# Replace ">" to "" in the decoys.txt file
sed -i.bak -e 's/>//g' $indir/decoys.txt

cd .. 
```


- **Creating gentrome.fa file**: salmon_decoy_gentrome.sh

```bash
#!/bin/bash

cd reference_GENCODE

indir=salmon_index

# Assign location of transcripts and genome reference files 
transcripts=*transcripts.fa 
genome=*genome.fa

# Concatenate them and save as.gentrome.fa
cat $transcripts $genome > $indir/gentrome.fa
```


- **Indexing**: salmon_index.sh

```bash
#!/bin/bash

# Move to the directory where gentrome.fa and decoys.txt files have been created 
cd reference_GENCODE/salmon_index 

salmon index -t gentrome.fa -d decoys.txt -p 8 -i gencode_index --gencode

cd ../..
```

#### 4-2. Mapping

- Flags used in this analysis is decribed in [my previous Salmon workflow](https://github.com/Mira0507/salmon_test/blob/master/workflow.md) except **-r**
- Since this is a paird-end analysis, **-1** and **-2** are used instead of **-r** 


```bash
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
```

### 5. STAR alignment 

- Flags are described in [my previous STAR workflow](https://github.com/Mira0507/star_test/blob/master/workflow_ens.md).

#### 5-1. Indexing

- star_index.sh

```bash 

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
``` 

#### 5-2. Alignment

- star_align.sh

```bash

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
```


### 6. HISAT2 

- Flags are described in [my previous HISAT2 workflow](https://github.com/Mira0507/hisat2_test/blob/master/workflow.md)

#### 6-1. Indexing

- hisat2_index.sh

```bash
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
```

#### 6-2. Alignment

- hisat2_align.sh

```bash
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
```

#### 6-3. Converting SAM to BAM 

- uses samtools
- hisat2_samtobam.sh

```bash
#!/bin/bash

cd hisat2_output

input=(A{1..3} B{1..3})  

for f in ${input[*]}

do
    samtools view -bS $f.sam > $f.bam 
done 


cd ..
```


#### 6-4. Sorting 

- uses samtools
- hisat2_sort.sh

```bash
#!/bin/bash

cd hisat2_output

input=(A{1..3} B{1..3})  

for f in ${input[*]}

do
    samtools sort $f.bam -o $f.sorted.bam

done 


cd ..

# Delete SAM files after samtools run
```
