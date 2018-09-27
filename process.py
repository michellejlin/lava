import argparse
import sys
from Bio import Entrez
import re
import subprocess
Entrez.email = 'uwvirongs@gmail.com'


def read_fasta(fasta_file_loc):
    strain_list = []
    genome_list = []
    dna_string = ''

    for line in open(fasta_file_loc):
        if line[0] == '>':
            strain_list.append(line[1:].split()[0])
            if dna_string != '':
                genome_list.append(dna_string)
                dna_string = ''
        else:
            dna_string += line.strip()
    genome_list.append(dna_string.replace('U', 'T'))
    return strain_list, genome_list

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

def adjust(given_num, our_num_array, ref_num_array, genome):

    # Handles gene lengths that go off the end of the genome
    if given_num == len(genome):
        return len(genome)

    # Go through our number array and search for the number of interest
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

    # now index is the absolute location of what we want
    if found:
        return str(our_num_array[index])
    else:
        return str(len(genome))


if __name__ =='__main__':
    parser = argparse.ArgumentParser(description='Internal tool to pull reference gff and fasta')
    parser.add_argument('ref', help='provide a Genbank acession number to use as a reference')
    parser.add_argument('fasta', help='Reference Fasta')
    
    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)
    fasta = args.fasta
    ref_seq_gb = args.ref
    record = Entrez.read(Entrez.esearch(db='nucleotide', term=ref_seq_gb))
    h2 = Entrez.efetch(db='nucleotide', id=record["IdList"][0], rettype='gb', retmode='text')
    e = open('lava_ref.gbk', 'w')
    e.write(h2.read())
    e.close()  
    
    h = Entrez.efetch(db='nucleotide', id=record["IdList"][0], rettype='fasta', retmode='text')
    d = open('lava_ref.fasta', 'w')
    d.write(h.read())
    d.close()
    
    gene_loc_list = []
    gene_product_list = []
    allow_one = False
    
    for line in open('lava_ref.gbk'):
        if ' CDS ' in line:
            gene_loc_list.append(re.findall(r'\d+', line))
            allow_one = True

        if '/product="' in line and allow_one:
            allow_one = False
            px = line.split('=')[1][1:-2]
            px = px.split()[0]
            gene_product_list.append(px)
	# once we reach this point we have a list of our reference genes and locations which we need to align to the provided reference sequence
	z = open('aligner.fasta', 'w')
    fe = open(fasta)
    for line in fe:
        z.write(line)
    fe.close()
    ge = open('lava_ref.fasta')
    z.write('\n')
    for line in ge:
        z.write(line)
    ge.close()
    z.close()
    nullo, genome = read_fasta(fasta)
    subprocess.call('mafft --quiet aligner.fasta > lava.ali', shell=True)
    ali_list, ali_genomes = read_fasta('lava.ali')
    ref_seq = ali_genomes[1]
    our_seq = ali_genomes[0]
    our_seq_num_array, ref_seq_num_array = build_num_arrays(our_seq, ref_seq)                
    
    for entry in range(0, len(gene_loc_list)):
        for y in range(0, len(gene_loc_list[entry])):
            gene_loc_list[entry][y] = adjust(int(gene_loc_list[entry][y]), our_seq_num_array, ref_seq_num_array, genome)   
             
    g = open('lava_ref.gff', 'w')
    g.write('##gff-version 3\n##source-version geneious 9.1.7\n')
    name= 'lava'
    for x in range(0, len(gene_product_list)):
        g.write(name + '\tGeneious\tgene\t' + gene_loc_list[x][0] + '\t' + gene_loc_list[x][1] + '\t.\t+\t.\tID=gene:' + gene_product_list[x] + ';biotype=protein_coding\n')
        g.write(name + '\tGeneious\tCDS\t' +  gene_loc_list[x][0] + '\t' + gene_loc_list[x][1] + '\t.\t+\t0\tID=CDS:' + gene_product_list[x] + ';Parent=transcript:' + gene_product_list[x] + ';biotype=protein_coding\n')
        g.write(name + '\tGeneious\ttranscript\t' + gene_loc_list[x][0] + '\t' + gene_loc_list[x][1] + '\t.\t+\t.\tID=transcript:' + gene_product_list[x] + ';Parent=gene:' + gene_product_list[x] + ';biotype=protein_coding\n')
