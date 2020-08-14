# This script pulls reference fasta and reference GenBank file given accession number specified by --GENBANK.

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

Entrez.email = 'vpeddu@uw.edu'
# Pulls reference GenBank file from Entrez
record = Entrez.read(Entrez.esearch(db='nucleotide', term= sys.argv[1]))
h2 = Entrez.efetch(db='nucleotide', id=record["IdList"][0], rettype='gb', retmode='text')
e = open('lava_ref.gbk', 'w')
e.write(h2.read())
e.close()

# Pulls reference fasta from Entrez
h = Entrez.efetch(db='nucleotide', id=record["IdList"][0], rettype='fasta', retmode='text')
d = open('lava_ref.fasta', 'w')
d.write(h.read())
d.close()

#lava.py -q NC_039477.1 Example2_file1.fastq Example2_metadata.csv -o Example2_output
