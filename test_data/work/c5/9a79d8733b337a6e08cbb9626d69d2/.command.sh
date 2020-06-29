#!/bin/bash

python3 /Users/uwvirongs/Documents/Michelle/lava/bin/Annotate_complex_mutations.py Example1_file2.fastq.csv 1	

mv complex.log Example1_file2.fastq.complex.log

mv reads.csv Example1_file2.fastq.reads.csv
