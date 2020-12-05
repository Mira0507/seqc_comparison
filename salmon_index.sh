#!/bin/bash

# Move to the directory where gentrome.fa and decoys.txt files have been created 
cd reference_GENCODE/salmon_index 

salmon index -t gentrome.fa -d decoys.txt -p 8 -i gencode_index --gencode

cd ../..
