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

tmp_mat_peptides = open("mat_peptides_additions.txt", "w")

# Original mat_peptides.txt in format: name, start, end
for mat_peptide in open("mat_peptides.txt"):
    mat_peptide_start = int(mat_peptide.split(',')[1])
    mat_peptide_name = mat_peptide.split(',')[0]

    for line in open("proteins.csv"):
        protein_start = int(line.split(',')[1])
        protein_end = int(line.split(',')[2])

        # Adds fourth column in mat_peptides_additions.txt that is the number of nucleotides
        # the mature peptide is from the original protein's start, which we will fix in ribosomal_slippage.py.
        if mat_peptide_start >= protein_start and mat_peptide_start <= protein_end:
            correction_number = mat_peptide_start - protein_start
            if "_rib_" in mat_peptide_name:
                before_slip = int(mat_peptide_name.split("_rib_")[1])
                correction_number = correction_number - before_slip
            tmp_mat_peptides.write(mat_peptide.rstrip() + "," + str(correction_number) + "\n")