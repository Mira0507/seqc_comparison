#!/bin/bash

cd rawdata 

pref=SRR9500
suff=.fastq

mv ${pref}78_1${suff} A1_1$suff
mv ${pref}78_2${suff} A1_2$suff
mv ${pref}79_1${suff} B1_1$suff
mv ${pref}79_2${suff} B1_2$suff
mv ${pref}80_1${suff} A2_1$suff
mv ${pref}80_2${suff} A2_2$suff
mv ${pref}81_1${suff} B2_1$suff
mv ${pref}81_2${suff} B2_2$suff
mv ${pref}82_1${suff} A3_1$suff
mv ${pref}82_2${suff} A3_2$suff
mv ${pref}83_1${suff} B3_1$suff
mv ${pref}83_2${suff} B3_2$suff

cd ..
