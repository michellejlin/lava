# Checks for complex mutations and prints a warning message.

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

print(sys.argv[1])
print(sys.argv[2])

complex_list = []
line_list = []
last = ''
last_line = ''
add_a_complex = False
seen_one_complex = False
complex_log_name = 'complex.log'
complex_log = open(complex_log_name, 'w')
for line in open(sys.argv[1]):
    # Make sure we only read lines that contain data 
    if ',' in line:
        # Detects complex mutations by seeing if the amino acid residue is the same as the last one we read.
        # This assumes that merged.csv is in sorted order, which we guarantee that it is.
        # Prints a warning message if complex mutation is detected.
        if line.split(',')[4][:-1] == last:
            if not seen_one_complex:
                print('WARNING: Complex mutation detected (multiple nucleotide changes within the same codon)! Sample=' + line.split(',')[0] + 
                    ' Change1=' + line.split(',')[4] + ' and Change2=' + last_line.split(',')[4] + ' this may happen again but future warnings will be suppressed. You can find the complete list in the output folder under complex.log.')
                seen_one_complex = True
            
            complex_log.write('Sample=' + line.split(',')[0] + ' Change1=' + line.split(',')[4] + ' and Change2=' + last_line.split(',')[4])
            complex_list.append(last_line)
            add_a_complex = True

        last = line.split(',')[4][:-1]
        
        # If we don't have time data in the metadata, add it here.
        if line.strip().split(',')[10] == '':
            line = line.strip() + str(sys.argv[2]) + '\n'
        if add_a_complex:
            complex_list.append(line)
            add_a_complex = False

        line_list.append(line)
        last_line = line

g = open(sys.argv[1], 'w')

# Rewrites file with 'complex' as change type.
for line in line_list:
    if line in complex_list:
        g.write(','.join(line.split(',')[:8]) + ',complex,' + ','.join(line.split(',')[9:]))
    else:
        g.write(line)
g.close()