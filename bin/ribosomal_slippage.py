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
import math

correction_number = 0
residue_correction_number = 0

# # -1 for coronavirus, HIV
# slippage_number = -1
#
# # Grabs the start of the CDS for the original protein before ribosomal slip
# slippage_cd_start = open("ribosomal_start.txt").read()
#
# for line in open("proteins.csv"):
#     if '_ribosomal_slippage' in line:
#         slip_site = line.split(',')[1]
#         residue_correction_number = int(slip_site) - int(slippage_cd_start) - slippage_number
#
#         # For some reason, nucleotide counting is off but residue number is correct.
#         # For now, adding back protein start is vaguely correct.
#         correction_number = int(slip_site) - 1


# Ribosomal_corrected has corrected ribosomal slippage annotations.
# Visualization is what will be fed into Bokeh, which will include old and new annotations.
ribosomal_corrected = open('final_corrected_slippage.csv', 'w')
visualization = open('visualization.csv', 'w')

correction_number = int(correction_number)
residue_correction_number = int(residue_correction_number)

# Looping through all mutations we found
with open("final.csv") as f:
    # Adds headers as appropriate to each file
    header = f.readline()
    ribosomal_corrected.write(header.rstrip() + ",MatPeptide" + '\n')
    visualization.write(header.rstrip() + ",NucCorrect,AminoCorrect,MatPeptide" + '\n')

    for row in f:
        line = row.rstrip()
        # Finds the nucleotide and amino acid numbers that need to be changed.
        # Formatting is different for deletions because of extra 'del.'
        nuc = line.split(',')[6]
        if 'del' in nuc:
            nuc_num = int(nuc[0:-4])
        else:
            nuc_num = int(nuc[1:-1])
        amino = line.split(',')[4]
        amino_num = int(amino[1:-1])
        position = int(line.split(',')[2])

        # position = line.split(',')[2]
            # if position < int(start_num_list[index]):
            #     index = index + 1

        mat_peptide = "-"
        # Going through and checking which mature peptide it falls under
        for mature_peptide in open("mat_peptides_additions.txt"):
            mat_name = mature_peptide.split(',')[0]
            mat_start = int(mature_peptide.split(',')[1])
            mat_end = int(mature_peptide.split(',')[2])
            mat_correction = int(mature_peptide.split(',')[3])

            # Check to see if mutation falls within this mature peptide
            if position >= mat_start and position <= mat_end:
                # Subtracts the difference between mature peptide and start of protein
                mat_nuc_num = position - mat_start + 1
                if (mat_name == "RNA-dependent_RNA_polymerase_rib_26"):
                    mat_nuc_num = mat_nuc_num + 27
                mat_aa_num = math.ceil(mat_nuc_num/3)

                # Grabs correct nucleotide annotation.
                if 'del' in nuc:
                    mat_nuc = str(mat_nuc_num) + nuc[-4:]
                else:
                    mat_nuc = nuc[0] + str(mat_nuc_num) + nuc[-1]

                # Grabs correct amino acid mutation.
                mat_aa = amino[0] + str(mat_aa_num) + amino[-1]

                # Writes full mature peptide annotation.
                mat_peptide = mat_name + ": " + mat_aa + "; " + mat_nuc

        # Corrects for ribosomal slippage by adding correction_number to
        # original nucleotide/residue number.
        if '_ribosomal_slippage' in line:
            if 'del' in nuc:
                del_nuc_replacement = str(nuc_num + correction_number) + nuc[-4:]
                new_line = line.replace(nuc, del_nuc_replacement)
                visualization_line = new_line + "," + del_nuc_replacement
            else:
                nuc_replacement = nuc[0] + str(nuc_num + correction_number) + nuc[-1]
                new_line = line.replace(nuc, nuc_replacement)
                visualization_line = new_line + "," + nuc_replacement

            # Round up if decimal, which should only be if ribosomal slippage happens?
            amino_replacement = int(nuc_num) + residue_correction_number
            amino_replacement = math.ceil(amino_replacement / 3)
            amino_replacement = amino[0] + str(amino_replacement) + amino[-1]
            new_line = new_line.replace(amino, amino_replacement)

            visualization_line = visualization_line.replace(amino, amino_replacement)
            visualization_line = visualization_line + "," + amino_replacement

            ribosomal_corrected.write(new_line + "," + mat_peptide + '\n')
            visualization.write(visualization_line + "," + mat_peptide + '\n')
        else:
            ribosomal_corrected.write(line + "," + mat_peptide + '\n')
            visualization.write(line + "," + nuc + "," + amino + "," + mat_peptide + '\n')
