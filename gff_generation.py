import subprocess
import re
import argparse
import timeit
import os
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
import platform
import sys
from Bio import Entrez
import time
import datetime
import shutil
Entrez.email = 'uwvirongs@gmail.com'

# Reads in a fasta file that should have strain names for the names of the sequences -  can handle any number of
# sequences. Also strips leading and trailing Ns or ?s from the provided sequence. Returns two lists with the names of
# the strains in the first one and the genomes as strings in the second list, also changes U's to T's
def read_fasta(fasta_file_loc):
    strain_list = []
    genome_list = []
    dna_string = ''

    for line in open(fasta_file_loc):
        if line[0] == '>':
            strain_list.append(line[1:].split()[0])
            if dna_string != '':
                # strip leading and trailing Ns or ?s because there's no reason to submit them
                xip = 0
                while dna_string[xip] == 'N' or dna_string[xip] == '?':
                    xip += 1

                y = len(dna_string)
                while dna_string[y-1] == 'N' or dna_string[y-1] == '?':
                    y -= 1

                dna_string = dna_string[xip:y]

                genome_list.append(dna_string)
                dna_string = ''
        else:
            dna_string += line.strip()
    # Just to make sure all our sequences are on the same page
    genome_list.append(dna_string.replace('U', 'T'))
    return strain_list, genome_list


# Spell checking functionality provided by Entrez
def spell_check(query_string):
    print('in spellcheck')
    handle = Entrez.espell(term=query_string)
    record = Entrez.read(handle)
    corrected_query = record["CorrectedQuery"]
    # Entrez returns blank strings for numerals or things that are spelled correctly
    # Since this is based on NCBI's spell checking protein names are included and correct
    # However this won't correct SUPER messed up words or made up words
    if corrected_query != '':
        print('Checking spelling on ' + query_string)
        print(query_string + ' was corrected to: ' + corrected_query)
        return corrected_query
    else:
        return query_string


# This function takes the strain name and a location of the individual fasta file we saved earlier and runs a blast
# Search saving the top 35 hits - then the top hit is found that is a complete genome and the fasta and .gbk of that
# are saved - then we run alignment on the two and return two strings one of them our sequence with '-' and the
# other of our new reference sequence
def align(our_fasta_loc, ref_fasta, ref_gb):
    #makes aligner.fasta for MAAFT which has the two fastas together
    s = 'cat ' + ref_fasta + ' > aligner.fasta'
    subprocess.call(s, shell=True)
    s = 'cat ' + our_fasta_loc + ' >> aligner.fasta'
    subprocess.call(s, shell=True)
	
    # Windows
    if SLASH == '\\':
        s = 'mafft-win\\mafft.bat --quiet aligner.fasta > aligner.ali'
        subprocess.call(s, shell=True)
    else:
        try:
            subprocess.call('mafft --quiet aligner.fasta > aligner.ali',
                    shell=True)
        except:
            print('Running on a non windows system, which means you need to install mafft and put it on the sys path '
                  'yourself.\nI suggest using brew or apt')
    ali_list, ali_genomes = read_fasta('aligner.ali')

    # SWAPPED THESE DURING DEBUGGING
    our_seq = ali_genomes[1]
    ref_seq = ali_genomes[0]

    return our_seq, ref_seq


# Takes in two sequences with gaps inserted inside of them and returns arrays that have a -1 in the gap locations and
# count up from 1 in the nucleotide areas - This data structure allows for extremely rapid conversion between relative
# locations in the two sequences although does assume that these genes are of uniform length
# NOTE: This means that when we have reads that like don't have the start codons of the first gene or something we'll
# get a -1 for the start location on our annotation
def build_num_arrays(our_seq, ref_seq):
    ref_count = 0
    our_count = 0
    ref_num_array = []
    our_num_array = []

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
	

# Takes a gene start index relative to an unaligned reference sequence and then returns the location of the same start
# area on the unaligned sequence that we're annotating using the number arrays to finish
def adjust(given_num, our_num_array, ref_num_array):

    # Go through our number array and search for the number of interest
    if our_num_array[given_num] == '-1':
        in_dex = given_num
        while our_num_array[in_dex != '-1']:
            in_dex += 1
            break
        print(in_dex)
        print(str(our_num_array[in_dex]))
        return str(our_num_array[in_dex])

    else:
        found = False
        for x in range(0, len(our_num_array)):
            if ref_num_array[x] == given_num:
                index = x
                print("x:",x)
                found = True
                break
    # now index is the absolute location of what we want
    return str(our_num_array[index])


# this opens up the reference .gbk file and pulls all of the annotations, it then adjusts the annotations to the
# relative locations that they should appear on our sequence
def pull_correct_annotations(our_seq, ref_seq, ref_gb):
    # Read the reference gbk file and extract lists of all of the protein locations and annotations!
    gene_loc_list = []
    gene_product_list = []
    allow_one = False

    for line in open(ref_gb):
        if ' CDS ' in line:
            # this is now going to be a list of numbers, start-stop start-stop
            # this line simply makes sure we read in reversed start-stops in the true reversed direction
            if 'complement' in line:
                whack = re.findall(r'\d+', line)
                whack.reverse()
                gene_loc_list.append(whack)
            else:
                gene_loc_list.append(re.findall(r'\d+', line))
            allow_one = True

        if '/gene="' in line and allow_one:
            allow_one = False
            # Inconsistent naming of protein products
            px = line.split('=')[1][1:-2]
            if args.no_spell_check:
                px_word_list = px.split()
                for word in px_word_list:
                    if '1' or '2' or '3' or '4' or '5' or '6' or '7' or '8' or '9' or '0' not in word: 
                        word = spell_check(word)

                px = ' '.join(px_word_list)
            if px == 'phospho protein':
                px = 'phoshoprotein'
            gene_product_list.append(px)
            
    our_seq_num_array, ref_seq_num_array = build_num_arrays(our_seq, ref_seq)
    
    
    entries_to_delete=[]
    to_delete = ["'", "C", "Y"]
    for entry in range(0, len(gene_product_list)):
    	protein = gene_product_list[entry]
    	if("P/V/C" in protein):
    		continue
    	elif any(c in protein for c in to_delete):
    		entries_to_delete.append(entry)
    for i in sorted(entries_to_delete, reverse=True):
    	del(gene_product_list[i])
    	del(gene_loc_list[i])
    	
    # Adjust every locus so that we actually put in correct annotations	
    for entry in range(0, len(gene_loc_list)):
        for y in range(0, len(gene_loc_list[entry])):
            gene_loc_list[entry][y] = adjust(int(gene_loc_list[entry][y]), our_seq_num_array, ref_seq_num_array)
    return gene_loc_list, gene_product_list

# Creates the gff for the new transferred annotated proteins
def write_gff (gene_loc_list, gene_product_list, our_fasta_loc):
	name = our_fasta_loc.split(".fasta")[0]
	w = open(name + '.gff', "w")
	w.write('##gff-version 3\n##source-version lava\n')
	for entry in range(0, len(gene_loc_list)):
		start_num = gene_loc_list[entry][0]
		end_num = gene_loc_list[entry][1]
		w.write(name + '\tLava\ttranscript\t' + start_num + '\t' + end_num + '\t.\t+\t.\tID=transcript:' + 
			gene_product_list[entry] + ';Parent=gene:' + gene_product_list[entry] + ';biotype=protein_coding\n')
		w.write(name + '\tLava\tCDS\t' + start_num + '\t' + end_num + '\t.\t+\t0\tID=CDS:' + 
			gene_product_list[entry] + ';Parent=transcript:' + gene_product_list[entry] + ';biotype=protein_coding\n')
		w.write(name + '\tLava\tgene\t' + start_num + '\t' + end_num + '\t.\t+\t.\tID=gene:' + 
			gene_product_list[entry] + ';biotype=protein_coding\n')

# quick check to make sure slashes go the right way on both Windows and Mac/Linux
def check_os():
    if platform.system() == 'Linux' or platform.system() == 'Darwin':
        return '/'
    else:
        return '\\'


if __name__ == '__main__':

    SLASH = check_os()
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('fasta_file', help='Input file in .fasta format containing complete or near complete '
                                           'genomes for all the viruses that you want to have annotated')
    parser.add_argument('ref_fasta', help='Reference fasta pulled off of Genbank')
    parser.add_argument('ref_gb', help='Reference Genbank annotations pulled off of Genbank')
    parser.add_argument('--no_spell_check', action='store_false', help='Turn off the automatic spellchecking for protein annoations ')
    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    our_fasta_loc = args.fasta_file
    ref_fasta = args.ref_fasta
    ref_gb = args.ref_gb
    our_seq, ref_seq = align(our_fasta_loc, ref_fasta, ref_gb)
    gene_loc_list, gene_product_list = pull_correct_annotations(our_seq, ref_seq, ref_gb)
    write_gff(gene_loc_list, gene_product_list, our_fasta_loc)