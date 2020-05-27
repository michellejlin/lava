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

start_num = ""
# Grabs where ribosomal slippage happens from proteins.csv.
start_num_list = []
for line in open("proteins.csv"):
    if '_ribosomal_slippage' in line:
        start_num = line.split(',')[1]
        if not start_num in start_num_list:
            start_num_list.append(start_num)

# Writes corrected lines in new file.
temp = open('temp.csv', 'w')
index = 0
if start_num != "":
    for line in open("merged.csv", 'r'):
        # Adds the number of where ribosomal slippage happens to where Annovar begins.
        if '_ribosomal_slippage' in line:
            nuc = line.split(',')[6]
            if 'del' in nuc:
                nuc_num = int(nuc[0:-4])
                print("del nuc_num" + str(nuc_num))
            else:
                nuc_num = int(nuc[1:-1])
            amino = line.split(',')[4]
            amino_num = int(amino[1:-1])
            position = line.split(',')[2]
            if position < int(start_num_list[index]):
                index = index + 1
            if 'del' in nuc:
                new_line = line.replace(nuc, str(nuc_num + int(start_num_list[index])) + nuc[-4:])
            else:
                new_line = line.replace(nuc, nuc[0] + str(nuc_num + int(start_num_list[index])) + nuc[-1])
            amino_num2 = int(start_num_list[index])/3
            new_line = new_line.replace(amino, amino[0] + str(amino_num + amino_num2) + amino[-1])
            temp.write(new_line)
        else:
            temp.write(line)