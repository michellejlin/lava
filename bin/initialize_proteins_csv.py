# Creates list of protein names and locations (proteins.csv) based on GFF annotations.

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
    
# Pulls protein information from GFF into proteins.csv. 
# $12 = protein name, $4 = beginning nucleotide, $5 = ending nucleotide.
# Sorts by increasing order of numbers.
subprocess.call('grep "ID=transcript:" ' + 'lava_ref.gff' + ' | awk -F\'[\t;:]\' \'{print $12 "," $4 "," $5}\' | sort -t \',\' -k2 -n > ' + 
    'proteins.csv', shell=True)
    
# Checks if protein list has duplicates, and prints error message.
df = pd.read_csv("proteins.csv", names = ["Protein","First", "Last"])
col = list(df.Protein)
for i in col:
    if col.count(i) > 1:
        sys.exit("ERROR: Duplicate protein names found. Use custom GFF file without duplicates.")