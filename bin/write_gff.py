# Creates a GFF (lava_ref.gff) for our consensus fasta with adjusted numbers according to reference,
# using CDS annotations from our reference GenBank file.
# Also grabs ribosomal slippage information and mature peptide information if applicable.

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

# Reads in a fasta and returns fasta name, fasta sequence with Ts instead of Us.
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

# Takes our consensus fasta and a reference fasta. 
# Builds arrays for both so we can adjust numberings downstream.
def build_num_arrays(our_seq, ref_seq):
	ref_count = 0
	our_count = 0
	ref_num_array = []
	our_num_array = []

	# Loops through the GenBank reference fasta sequence.
	for x in range(0, len(ref_seq)):
		if ref_seq[x] != '-':
			ref_count += 1
			ref_num_array.append(ref_count)
		else:
			ref_num_array.append(-1)

		if our_seq[x] != '-':
			our_count += 1
			our_num_array.append(our_count)
		else:
			our_num_array.append(-1)

	return our_num_array, ref_num_array

# Takes number arrays and adjusts genome position in our consensus fasta
# to match corresponding position in reference fasta.
def adjust(given_num, our_num_array, ref_num_array, genome):
	
	# Handles gene lengths that go off the end of the genome
	if given_num == len(genome):
		return len(genome)

	# Goes through our number array and searches for the number of interest
	if our_num_array[given_num] == '-1':
		in_dex = given_num
		while our_num_array[in_dex != '-1']:
			in_dex += 1
			break
		return str(our_num_array[in_dex])

	else:
		found = False
		for x in range(0, len(our_num_array)):
			if ref_num_array[x] == given_num:
				index = x
				found = True
				break

	# Now index is the absolute location of what we want
	if found:
		return str(our_num_array[index])
	else:
		return str(len(genome))    

# Reads in consensus fasta.
ignore, genome_ref = read_fasta('consensus.fasta')
g = open('consensus.fasta', 'w')
# Replaces consensus fasta header (pulled from GenBank) with >lava and writes the fasta sequence (with Ts instead of Us).
g.write(">lava\n")
g.write(genome_ref[0])
g.close()

gene_loc_list = []
gene_product_list = []
allow_one = False

mat_peptide_loc_list = []
mat_peptide_product_list = []
mat_peptide_current_loc = []
allow_one_mat = False

# Pulls CDS annotations from downloaded GenBank file.
for line in open('lava_ref.gbk'):
	# Grabs locations of proteins from CDS annotations.
	if 'CDS' in line and '..' in line:
		gene_loc_list.append(re.findall(r'\d+', line))
		allow_one = True
	# Grabs protein names from CDS annotations.
	if '/product="' in line and allow_one:
		allow_one = False
		# Example product annotation: /product="ORF1ab polyprotein"
		# Grabs name of gene in between quotes and replaces spaces with underscores to prevent issues downstream.
		px = line.split('=')[1][1:-2]
		px = px.replace(' ', '_')
		gene_product_list.append(px)

	# Pulls mature peptide annotations (for coronavirus: long polyproteins).
	if 'mat_peptide' in line and '..' in line:
		mat_peptide_current_loc = re.findall(r'\d+', line)
		allow_one_mat = True
	if '/product="' in line and allow_one_mat:
		allow_one_mat = False
		px = line.split('=')[1][1:-2]
		px = px.replace(' ', '_')
		# Mature peptides can repeat in GenBank files, don't want that, so limit to one.
		if px not in mat_peptide_product_list:
			mat_peptide_product_list.append(px)
			mat_peptide_loc_list.append(mat_peptide_current_loc)
	
# Writes a new file for MAFFT input 
z = open('aligner.fasta', 'w')
fe = open('consensus.fasta')
for line in fe:
	z.write(line)
fe.close()
ge = open('lava_ref.fasta')
z.write('\n')
for line in ge:
	z.write(line)
ge.close()
z.close()

# Aligns our conensus fasta to the reference fasta with MAFFT.
nullo, genome = read_fasta('consensus.fasta')
subprocess.call('/usr/local/miniconda/bin/mafft --quiet ' + 'aligner.fasta > ' + 'lava.ali', shell=True)
ali_list, ali_genomes = read_fasta('lava.ali')
ref_seq = ali_genomes[1]
our_seq = ali_genomes[0]
our_seq_num_array, ref_seq_num_array = build_num_arrays(our_seq, ref_seq)

# Goes through every annotation and adjusts the coordinates from the reference to the consensus sequence
for entry in range(0, len(gene_loc_list)):
	for y in range(0, len(gene_loc_list[entry])):
		gene_loc_list[entry][y] = adjust(int(gene_loc_list[entry][y]), our_seq_num_array, ref_seq_num_array, genome)   

# Adjusts for mature peptides too.
for entry in range(0, len(mat_peptide_loc_list)):
	for y in range(0, len(mat_peptide_loc_list[entry])):
		mat_peptide_loc_list[entry][y] = adjust(int(mat_peptide_loc_list[entry][y]), our_seq_num_array, ref_seq_num_array, genome)   

# Writes mature peptide names and locations to new file.
mat_peptide_list = open('mat_peptides.txt', 'w')
for x in range(0, len(mat_peptide_product_list)):
	# Some peptides have ribosomal slippage too, 
	# grab the last two locations (after slippage join)
	if len(mat_peptide_loc_list[x]) == 4:
		peptide_start = str(mat_peptide_loc_list[x][2])
		peptide_end = str(mat_peptide_loc_list[x][3])
		before_slip = int(mat_peptide_loc_list[x][1])-int(mat_peptide_loc_list[x][0])
		before_slip = int(before_slip)
		# Append the length of the protein before ribosomal slippage to name of mat_peptide,
		# to be parsed later by ribosomal_slippage.py
		mat_peptide_list.write(mat_peptide_product_list[x] + "_rib_" + str(before_slip) + 
		"," + peptide_start + "," + peptide_end + "\n")

	# Otherwise, just grab the locations of the mat peptide
	else:
		peptide_start = str(mat_peptide_loc_list[x][0])
		peptide_end = str(mat_peptide_loc_list[x][1])
		mat_peptide_list.write(mat_peptide_product_list[x] + "," + peptide_start + "," + peptide_end + "\n")
mat_peptide_list.close()

# Writes a new .gff file for annovar to use later 
print("Writing GFF file...")
g = open('lava_ref.gff', 'w')
g.write('##gff-version 3\n##source-version geneious 9.1.7\n')
name= 'lava'

slip_start = ""

for x in range(0, len(gene_product_list)):
# Hard coded for the bad HPIV3 reference
	gene_name = gene_product_list[x].replace("potein", "protein")
	gene_start = str(gene_loc_list[x][0])
	gene_end = gene_loc_list[x][1]
	if(gene_name == "hemagglutinin-neuraminidase"):
		gene_end = int(gene_end) - 6
		gene_name = "HN"
	elif(gene_name == "phosphoprotein"):
		gene_name = "P"
	elif(gene_name == "nucleocapsid_protein"):
		gene_name = "NP"
	elif(gene_name == "matrix_protein"):
		gene_name = "M"
	elif(gene_name == "large_protein"):
		gene_name = "L"
	elif(gene_name == "fusion_protein"):
		gene_name = "F"
	gene_end = str(gene_end)

	## For ribosomal slippage, creates fake new protein to get the second set of values
	if len(gene_loc_list[x]) == 4:
		g.write(name + '\tLAVA\tgene\t' + str(gene_loc_list[x][2]) + '\t' + str(gene_loc_list[x][3]) + '\t.\t+\t.\tID=gene:' + gene_product_list[x] + '_ribosomal_slippage;biotype=protein_coding\n')
		g.write(name + '\tLAVA\tCDS\t' +  str(gene_loc_list[x][2]) + '\t' + str(gene_loc_list[x][3]) + '\t.\t+\t0\tID=CDS:' + gene_product_list[x] + '_ribosomal_slippage;Parent=transcript:' + gene_product_list[x] + '_ribosomal_slippage;biotype=protein_coding\n')
		g.write(name + '\tLAVA\ttranscript\t' + str(gene_loc_list[x][2]) + '\t' + str(gene_loc_list[x][3]) + '\t.\t+\t.\tID=transcript:' + gene_product_list[x] + '_ribosomal_slippage;Parent=gene:' + gene_product_list[x] + '_ribosomal_slippage;biotype=protein_coding\n')

		slip_start = str(gene_loc_list[x][0])
	
	elif(gene_name != "C_protein" and gene_name != "D_protein"):
		g.write(name + '\tLAVA\tgene\t' + gene_start + '\t' + gene_end + '\t.\t+\t.\tID=gene:' + gene_name + ';biotype=protein_coding\n')
		g.write(name + '\tLAVA\tCDS\t' +  gene_start + '\t' + gene_end + '\t.\t+\t0\tID=CDS:' + gene_name + ';Parent=transcript:' + gene_name + ';biotype=protein_coding\n')
		g.write(name + '\tLAVA\ttranscript\t' + gene_start + '\t' + gene_end + '\t.\t+\t.\tID=transcript:' + gene_name + ';Parent=gene:' + gene_name + ';biotype=protein_coding\n')

# Writes where the the protein containing the ribosomal slippage starts in separate file, to be parsed later.
slip_start_file = open('ribosomal_start.txt', 'w')
slip_start_file.write(slip_start)