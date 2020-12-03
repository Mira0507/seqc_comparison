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
- Command line: **conda activate <env name>**

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
```

### 3. Reference files

- Version: [GENCODE](https://www.gencodegenes.org/human) GRCh38(hg38) release 36 (v36)
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

- Creating decoys.txt file
- salmon_decoy_txt.sh

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


- Creating gentrome.fa file
- salmon_decoy_gentrome.sh

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


- Indexing 
- salmon_index.sh

```bash
#!/bin/bash

# Move to the directory where gentrome.fa and decoys.txt files have been created 
cd reference_GENCODE/salmon_index 

salmon index -t gentrome.fa -d decoys.txt -p 8 -i gencode_index --gencode

cd ../..
```

#### 4-2. Mapping
