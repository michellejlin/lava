#!/bin/bash

python3 /Users/uwvirongs/Documents/Michelle/lava/bin/Annotate_complex_mutations.py Example1_file1.fastq.csv 0	

mv complex.log Example1_file1.fastq.complex.log
mv reads.csv Example1_file1.fastq.reads.csv
