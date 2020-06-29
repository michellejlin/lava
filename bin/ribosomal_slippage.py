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

correction_number = 0
residue_correction_number = 0
# -1 for coronavirus, HIV
slippage_number = -1

# grabs the start of the CDS for the original protein before ribosomal slip
slippage_cd_start = open("ribosomal_start.txt").read()

for line in open("proteins.csv"):
    if '_ribosomal_slippage' in line:
        slip_site = line.split(',')[1]
        correction_number = int(slip_site) - int(slippage_cd_start) - slippage_number
        residue_correction_number = correction_number / 3

        # For some reason, nucleotide counting is off but residue number is correct.
        # For now, adding back protein start is vaguely correct.
        correction_number = correction_number + int(slippage_cd_start)


# Writes corrected lines in new file.
temp = open('final_corrected_slippage.csv', 'w')
temp2 = open('visualization.csv', 'w')

correction_number = int(correction_number)
residue_correction_number = int(residue_correction_number)

with open("final.csv") as f:
    header = f.readline()
    temp.write(header)
    temp2.write(header.rstrip() + ",NucCorrect,AminoCorrect" + '\n')

    next(f)

    for row in f:
        line = row.rstrip()
        # Adds the number of where ribosomal slippage happens to where Annovar begins.
        nuc = line.split(',')[6]
        if 'del' in nuc:
            nuc_num = int(nuc[0:-4])
            print("del nuc_num" + str(nuc_num))
        else:
            nuc_num = int(nuc[1:-1])
        amino = line.split(',')[4]
        amino_num = int(amino[1:-1])

        # position = line.split(',')[2]
            # if position < int(start_num_list[index]):
            #     index = index + 1


        if '_ribosomal_slippage' in line:
            if 'del' in nuc:
                del_nuc_replacement = str(nuc_num + correction_number) + nuc[-4:]
                new_line = line.replace(nuc, del_nuc_replacement)
                temp2_line = new_line + "," + del_nuc_replacement
            else:
                nuc_replacement = nuc[0] + str(nuc_num + correction_number) + nuc[-1]
                new_line = line.replace(nuc, nuc_replacement)
                temp2_line = new_line + "," + nuc_replacement

            amino_replacement = amino[0] + str(amino_num + residue_correction_number) + amino[-1]
            new_line = new_line.replace(amino, amino_replacement)

            temp2_line = temp2_line.replace(amino, amino_replacement)
            temp2_line = temp2_line + "," + amino_replacement

            temp.write(new_line + '\n')
            temp2.write(temp2_line + '\n')
        else:
            temp.write(line + '\n')
            temp2.write(line + "," + nuc + "," + amino + '\n')
