## Analysis of SEQC datasets 

### 1. Rawdata 

- GEO number: [GSE49712](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49712)
- SRA number: [SRP028705](https://www.ncbi.nlm.nih.gov/sra?term=SRP028705)
- Reference paper: [Genome Biol. 2013;14(9):R95. doi: 10.1186/gb-2013-14-9-r95.](https://pubmed.ncbi.nlm.nih.gov/24020486)
- Sample preparation: The Sequencing Quality Control Consortium generated two datasets from two reference RNA samples in order to evaluate transcriptome profiling by next-generation sequencing technology. Each sample contains one of the reference RNA source and a set of synthetic RNAs from the External RNA Control Consortium (ERCC) at known concentrations. Group A contains 5 replicates of the Strategene Universal Human Reference RNA (UHRR), which is composed of total RNA from 10 human cell lines, with 2% by volume of ERCC mix 1. Group B includes 5 replicate samples of the Ambion Human Brain Reference RNA (HBRR) with 2% by volume of ERCC mix 2. The ERCC spike-in control is a mixture of 92 synthetic polyadenylated oligonucleotides of 250-2000 nucleotides long that are meant to resemble human transcripts.      
- Raw data was downloaded with [SRA toolkit](https://github.com/Mira0507/using_SRA/blob/master/README.md) 

#### 1-1. Downloading SRA files 

- **rawdata_download.sh**

```bash
#!/bin/bash


# Numbers of SRA data:
# SRR950078: A1
# SRR950079: B2
# SRR950080: A2
# SRR950081: B2
# SRR950082: A3
# SRR950083: B3
# SRR950084: A4
# SRR950085: B4



cd rawdata


for x in {78..85}
do
    prefetch SRR9500$x

done

cd .. 
```

#### 1-2. Converting to FASTQ files 

- ** sra_fastq.sh**

```bash
#!/bin/bash

# sra toolkit doc: https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=prefetch

# create a directory to store rawdata

cd rawdata



for x in SRR9500{78..85}
do 
    fastq-dump --split-files $x/$x.sra 
done

cd ..
```
