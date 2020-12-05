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
