#!/bin/bash

ls -lah

# cat *fastq.csv >> merged.csv

cat merged.csv > final.csv 
cat *.fastq.csv >> final.csv 

grep -v "transcript" final.csv > a.tmp && mv a.tmp final.csv 

grep -v "delins" final.csv > a.tmp && mv a.tmp final.csv 


# Corrects for ribosomal slippage.
python3 /Users/uwvirongs/Documents/Michelle/lava/bin/ribosomal_slippage.py final.csv proteins.csv

awk NF final.csv > a.tmp && mv a.tmp final.csv

cat *.reads.csv > reads.csv 

cat *.log > complex.log
# TODO error handling @ line 669-683 of lava.py 

 python3 /Users/uwvirongs/Documents/Michelle/lava/bin/genome_protein_plots.py visualization.csv proteins.csv reads.csv . "Plot"

mkdir vcf_files
mv *.vcf vcf_files
