import subprocess 
import argparse
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio import Entrez
import shutil
import sys
import time
from datetime import datetime
import re
import os.path
import pandas as pd
import sys

def read_fasta(fasta_file_loc):
    strain_list = []
    genome_list = []
    dna_string = ''

    # Opens the reference fasta file.
    for line in open(fasta_file_loc):
    	# Grabs strain names from fasta headers, delineated by ">".
        if line[0] == '>':
            strain_list.append(line[1:].split()[0])
            # For every new strain, add existing dna_string to genome list and continue.
            if dna_string != '':
                genome_list.append(dna_string)
                dna_string = ''
        # Grabs sequence with removed whitespace.
        else:
            dna_string += line.strip()
    # Replaces Us with Ts (causes problems downstream).
    genome_list.append(dna_string.replace('U', 'T'))
    return strain_list, genome_list

g = open('consensus.fasta', 'w')
ignore, genome_ref = read_fasta('consensus.fasta')
g.write(">lava\\n")
g.write(genome_ref[0])
g.close()



gene_loc_list = []
gene_product_list = []
allow_one = False

# Pulls CDS annotations from downloaded genbank file 
for line in open('Lava_ref.gbk'):
    if 'CDS' in line and '..' in line:
        gene_loc_list.append(re.findall(r'\\d+', line))
        allow_one = True
    if '/product="' in line and allow_one:
        allow_one = False
        px = line.split('=')[1][1:-2]
        px = px.replace(' ', '_')
        gene_product_list.append(px)

# Writes a new file for MAFFT input 
z = open('aligner.fasta', 'w')
fe = open('consensus.fasta')
for line in fe:
    z.write(line)
fe.close()
ge = open('lava_ref.fasta')
z.write('\\n')
for line in ge:
    z.write(line)
ge.close()
z.close()