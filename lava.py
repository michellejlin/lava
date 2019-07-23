## lava Version 1.02
## Longitudinal Analysis of Viral Alleles 

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
from install_check import check_picard, check_gatk, check_varscan

Entrez.email = 'uwvirongs@gmail.com'
VERSION = 'v1.02'

# Takes a file path pointing to a fasta file and returns two lists.
# The first list is a list of all the fasta headers and the second is a list of all the sequences
# for each record. These will always be in the same order. Also replaces Us with Ts.
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


# Takes two MAFFT-aligned sequences and returns two arrays of numbers corresponding
# to relative genome position for each. Any position that doesn't align will have a -1.
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


# Takes two number arrays built by build_num_arrays and adjusts a given nucleotide position from a reference 
# to a given genome, then returns this adjusted index as a string. 
# If the adjustment fails (for example off the edges of the alignments) a -1 will be returned.
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

# Takes a location to a two column metadata.csv sheet with names of each sample we want to analyze and time data. 
# Returns two lists, one of sample/file names and the other of their relative temporal relationship.
def read_metadata(filepath):
	sample_list = []
	sample_time_list = []
	first = True
	for line in open(filepath):
		if first:
			first = False
		# Protects against empty lines
		elif ',' in line:
			#print(line)
			sample_list.append(line.split(',')[0])
			sample_time_list.append(line.split(',')[1])

	if not os.path.isfile(sample_list[0]):
		print('METADATA file is incorrectly formatted!')
		print(sample_list[0] + ' was the first sample but could not be found in the working directory.')
		print('Please move your samples to the directory you are running lava from, this is:')
		subprocess.call('pwd', shell=True)
		print('Also make sure your samples are pointing to the correct directory and include the file extension .fastq.')
		sys.exit(1)
	return sample_list, sample_time_list


# Takes a fastq file, GenBank Accession number and output directory. Aligns fastq file to genbank accession number, 
# and extracts majority consensus as well as transfers annotations from Genbank to newly extracted consensus. 
# Returns a filepath to a fasta of the majority consensus and a filepath to a corresponding 
# .gff file containing protein annotations. Also writes these files into the output directory as well as all intermediate files. 
def process(ref_seq_gb, fastq, new_dir):

	# Pulls whole genbank reference from Entrez
	record = Entrez.read(Entrez.esearch(db='nucleotide', term=ref_seq_gb))
	h2 = Entrez.efetch(db='nucleotide', id=record["IdList"][0], rettype='gb', retmode='text')
	e = open(new_dir + '/lava_ref.gbk', 'w')
	e.write(h2.read())
	e.close()

	# Pulls reference fasta from Entrez
	h = Entrez.efetch(db='nucleotide', id=record["IdList"][0], rettype='fasta', retmode='text')
	d = open(new_dir + '/lava_ref.fasta', 'w')
	d.write(h.read())
	d.close()

	# Indexes reference fasta
	subprocess.call('bwa index ' + new_dir + '/lava_ref.fasta 2> ' + new_dir + '/lava.log', shell=True)

	## Deleted the -v flag to try to fix the dang consensus generation 
	# Aligns first fartq to downloaded reference fasta.
	subprocess.call('bwa mem -M ' + new_dir + '/lava_ref.fasta ' + fastq + ' > ' + new_dir + '/aln.sam 2>> ' + new_dir + '/lava.log', shell=True)
	# Converts sam to bam.
	subprocess.call('samtools view -S -b ' + new_dir + '/aln.sam > ' + new_dir + '/aln.bam 2>> '+ new_dir + '/lava.log', shell=True)
	# Sorts the bam file.
	subprocess.call('samtools sort ' + new_dir + '/aln.bam -o ' + new_dir + '/aln.sorted.bam 2>> ' + new_dir + '/lava.log', shell=True)
	
	# Creates vcf and indexes.
	# Added max-depth to consider the first 50,000 reads at each position and reduced the prior by 7 orders of magnitude 
	subprocess.call('bcftools mpileup --max-depth 500000 -P 1.1e-100 -Ou -f ' + new_dir + '/lava_ref.fasta ' + new_dir + '/aln.sorted.bam | bcftools call -m -Oz -o ' + 
		new_dir + '/calls.vcf.gz 2>> ' + new_dir + '/lava.log', shell=True)
	subprocess.call('tabix ' + new_dir + '/calls.vcf.gz 2>> ' + new_dir + '/lava.log', shell=True)

	# Filters first so consensus actually takes most the common allele.
	subprocess.call('gunzip ' + new_dir + '/calls.vcf.gz', shell=True)
	subprocess.call('bcftools filter -i \'(DP4[0]+DP4[1]) < (DP4[2]+DP4[3]) && ((DP4[2]+DP4[3]) > 0)\' ' + new_dir + '/calls.vcf -o ' 
		+ new_dir + '/calls2.vcf 2>>' + new_dir + '/lava.log', shell=True)
	subprocess.call('bgzip ' + new_dir + '/calls2.vcf 2>> ' + new_dir + '/lava.log', shell=True)
	subprocess.call('tabix ' + new_dir + '/calls2.vcf.gz 2>> ' + new_dir + '/lava.log', shell=True)

	# Creates the consensus "lava_ref.fasta" and writes it to the output folder 
	subprocess.call('cat ' + new_dir + '/lava_ref.fasta | bcftools consensus ' + new_dir + '/calls2.vcf.gz > ' + new_dir + '/consensus.fasta', shell=True)

	# Rewrites fasta with lava as the sequence header so that annovar can match this with the reference .gff
	fasta = new_dir + '/consensus.fasta'
	ignore, genome_ref = read_fasta(fasta)
	g = open(fasta, 'w')
	g.write('>lava\n')
	g.write(genome_ref[0])
	g.close()

	gene_loc_list = []
	gene_product_list = []
	allow_one = False

	# Pulls CDS annotations from downloaded genbank file 
	for line in open(new_dir + '/lava_ref.gbk'):
		if 'CDS' in line and '..' in line:
			gene_loc_list.append(re.findall(r'\d+', line))
			allow_one = True
		if '/product="' in line and allow_one:
			allow_one = False
			px = line.split('=')[1][1:-2]
			px = px.replace(' ', '_')
			gene_product_list.append(px)

	# Writes a new file for MAFFT input 
	z = open(new_dir + '/aligner.fasta', 'w')
	fe = open(fasta)
	for line in fe:
		z.write(line)
	fe.close()
	ge = open(new_dir + '/lava_ref.fasta')
	z.write('\n')
	for line in ge:
		z.write(line)
	ge.close()
	z.close()

	# Uses MAFFT to align consensus fasta with downloaded fasta 
	nullo, genome = read_fasta(fasta)
	subprocess.call('mafft --quiet ' + new_dir + '/aligner.fasta > ' + new_dir + '/lava.ali', shell=True)
	ali_list, ali_genomes = read_fasta(new_dir + '/lava.ali')
	ref_seq = ali_genomes[1]
	our_seq = ali_genomes[0]
	our_seq_num_array, ref_seq_num_array = build_num_arrays(our_seq, ref_seq)

	# Goes through every annotation and adjusts the coordinates from the reference to the consensus sequence
	for entry in range(0, len(gene_loc_list)):
		for y in range(0, len(gene_loc_list[entry])):
			gene_loc_list[entry][y] = adjust(int(gene_loc_list[entry][y]), our_seq_num_array, ref_seq_num_array, genome)   

	# Writes a new .gff file for annovar to use later 
	print("Writing GFF file...")
	g = open(new_dir + '/lava_ref.gff', 'w')
	g.write('##gff-version 3\n##source-version geneious 9.1.7\n')
	name= 'lava'
	for x in range(0, len(gene_product_list)):

		g.write(name + '\tLAVA\tgene\t' + str(gene_loc_list[x][0]) + '\t' + str(gene_loc_list[x][1]) + '\t.\t+\t.\tID=gene:' + gene_product_list[x] + ';biotype=protein_coding\n')
		g.write(name + '\tLAVA\tCDS\t' +  str(gene_loc_list[x][0]) + '\t' + str(gene_loc_list[x][1]) + '\t.\t+\t0\tID=CDS:' + gene_product_list[x] + ';Parent=transcript:' + gene_product_list[x] + ';biotype=protein_coding\n')
		g.write(name + '\tLAVA\ttranscript\t' + str(gene_loc_list[x][0]) + '\t' + str(gene_loc_list[x][1]) + '\t.\t+\t.\tID=transcript:' + gene_product_list[x] + ';Parent=gene:' + gene_product_list[x] + ';biotype=protein_coding\n')
		
		## For ribosomal slippage, creates fake new protein to get the second set of values
		# if len(gene_loc_list[x]) == 4:
		# 	g.write(name + '\tLAVA\tgene\t' + str(gene_loc_list[x][2]) + '\t' + str(gene_loc_list[x][3]) + '\t.\t+\t.\tID=gene:' + gene_product_list[x] + '_ribosomal_slippage;biotype=protein_coding\n')
		# 	g.write(name + '\tLAVA\tCDS\t' +  str(gene_loc_list[x][2]) + '\t' + str(gene_loc_list[x][3]) + '\t.\t+\t0\tID=CDS:' + gene_product_list[x] + '_ribosomal_slippage;Parent=transcript:' + gene_product_list[x] + '_ribosomal_slippage;biotype=protein_coding\n')
		# 	g.write(name + '\tLAVA\ttranscript\t' + str(gene_loc_list[x][2]) + '\t' + str(gene_loc_list[x][3]) + '\t.\t+\t.\tID=transcript:' + gene_product_list[x] + '_ribosomal_slippage;Parent=gene:' + gene_product_list[x] + '_ribosomal_slippage;biotype=protein_coding\n')


	# Returns the filepaths for the new fasta and .gff file as strings 
	return fasta, new_dir + '/lava_ref.gff'

# Goes through the merged.csv file and makes sure that the temporal metadata has been added correctly to the output.
# Scans for complex mutations and prints a helpful warning message to the terminal listing what samples and their changes.
# Re-writes the file with complex mutations being annotated as 'complex'.
def add_passage(sample, passage, the_dir):
	complex_list = []
	line_list = []
	last = ''
	last_line = ''
	add_a_complex = False
	seen_one_complex = False
	complex_log = open(the_dir + '/complex.log', 'w')
	for line in open(sample):
		# Make sure we only read lines that contain data 
		if ',' in line:
			# Detects complex mutations by seeing if the amino acid residue is the same as the last one we read.
			# This assumes that merged.csv is in sorted order, which we guaruntee that it is.
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
				line = line.strip() + str(passage) + '\n'
			if add_a_complex:
				complex_list.append(line)
				add_a_complex = False

			line_list.append(line)
			last_line = line

	g = open(sample, 'w')

	# Rewrites file with 'complex' as change type.
	for line in line_list:
		if line in complex_list:
			g.write(','.join(line.split(',')[:8]) + ',complex,' + ','.join(line.split(',')[9:]))
		else:
			g.write(line)
	g.close()

# Handles ribosomal slippage by adding a second protein after ribosomal slippage.
## Currently not in use, because of wonky nucleotide numbers by Annovar.
def add_ribosomal_slippage(new_dir):
	start_num = ""
	# Grabs where ribosomal slippage happens from proteins.csv.
	start_num_list = []
	for line in open(new_dir + "/proteins.csv"):
		if '_ribosomal_slippage' in line:
			start_num = line.split(',')[1]
			if not start_num in start_num_list:
				start_num_list.append(start_num)

	# Writes corrected lines in new file.
	temp = open(new_dir + '/temp.csv', 'w')
	index = 0
	if start_num != "":
		for line in open(new_dir + "/merged.csv", 'r'):
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
	## subprocess.call('mv ' + new_dir + '/temp.csv ' + new_dir + '/merged.csv', shell=True)


# Removes the bigger intermediate files that LAVA creates.
def clean_up(new_dir):
	subprocess.call('rm ' + new_dir + '/*.bam', shell=True)
	subprocess.call('rm ' + new_dir + '/*.sam', shell=True)
	subprocess.call('rm ' + new_dir + '/*.avinput', shell=True)
	subprocess.call('rm ' + new_dir + '/*.vcf*', shell=True)
	subprocess.call('rm ' + new_dir + '/*.pileup', shell=True)
	subprocess.call('rm ' + new_dir + '/*.bai', shell=True)
	subprocess.call('rm ' + new_dir + '/*variant_function', shell=True)
	subprocess.call('rm ' + new_dir + '/*.txt', shell=True)


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Version ' + VERSION + '\nLongitudinal analysis of minor variants across whole viral'
												 ' genomes. Use either -g and -f together or use -q to annotate control fastq.')
	parser.add_argument('-f', help='Specify a reference fasta with the majority consensus of the control fastq. This option must be used '
								   'with the -g flag to specify the protein annotations relative to the start of this fasta.')
	parser.add_argument('-g', help='Specify a reference gff file with the protein annotations for the reference fasta supplied with the -f flag.'
								   ' This option must be paired with the -f flag.')
	parser.add_argument('-q', help='Provide a Genbank accession number. This record will be used to generate a majority consensus from the '
								   'control fastq, and this consensus will be annotated from the downloaded genbank record as well.')
	parser.add_argument('-nuc', action='store_true', help='Results are listed as nucleotide changes not amino acid changes. Do not use with -png.')
	parser.add_argument('-af', help='Specify an allele frequency percentage to cut off - with a minimum of 1 percent - in whole numbers.')
	parser.add_argument('-png', action='store_true', help='Output results as a png. Do not use with -nuc.')
	parser.add_argument('-dedup', action='store_true', help='Optional flag, will perform automatic removal of PCR duplicates via DeDup.')
	parser.add_argument('control_fastq', help='Required argument: The fastq reads for the first sample in your longitudinal analysis')
	parser.add_argument('metadata', help='Required argument: A two column csv - the first column is the name of all the fastqs you wish '
										 'to include in your analysis. All fastqs that you want to include need to be specified in this '
										 'file AND be located in the folder from which you are running lava. The second column is the '
										 'temporal seperation between the samples.  This is unitless so you can input passage number, '
										 'days, or whatever condition your experiment happens to have.')
	parser.add_argument('-o', help='Optional flag to name the output folder that lava will stuff output into. If a name isn\'t provided '
								   'folder will be named lava-date.')
	parser.add_argument('-save', action='store_true', help='Optional argument to save intermediate alignment files (sams, bams, vcfs, ect) '
		'LAVA\'s default behavior is to remove these after use to save disk footprint.')

	# Checks for argument sanity.
	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(1)

	# If an output directory is not provided by -o, pulls the date and time in isoformat with no colons.
	if args.o != None:
		new_dir = args.o
	else:
		new_dir = str(datetime.now().isoformat())
		new_dir = '_'.join(new_dir.split(':'))

	# Creates flags for passing to genome_protein_plots.py 
	if args.nuc:
		nuc_flag = '-nuc'
	else:
		nuc_flag = ''

	if args.png:
		png_flag = '-png'
	else:
		png_flag = ''

	user_af = -1
	if args.af:
		user_af = '-af ' + args.af
	else:
		user_af = ''

	# Writes new output folder.
	subprocess.call('mkdir -p ' + new_dir, shell=True)

	# Setting variables
	dir_path = os.path.dirname(os.path.realpath(__file__))
	plot_title = new_dir
	metadata_location = args.metadata
	control_fastq = args.control_fastq
	sample_list, sample_time_list = read_metadata(metadata_location)

	# Checks for picard, gatk, and varscan, and exits if we can't find them.
	PICARD = check_picard(dir_path)
	GATK = check_gatk(dir_path)
	VARSCAN = check_varscan(dir_path)

	# Copies fastqs to new directory - slow but protects user data from the insanity that is about to occur.
	subprocess.call('cp -p ' + control_fastq + ' ' + new_dir + '/', shell=True)
	control_fastq = new_dir + '/' + control_fastq

	# Makes sure that we've got a way of pulling annotations, if the user gives -f -g and -q then we only use -f and -g.
	# If user-provided fasta and gff, copies into directory.
	if args.f != None and args.g != None:
		print('Using -f and -g flags to annotate control fastq, -q flag will be ignored.')
		print('This method of reference generation assumes that your fasta and .gff file are formatted correctly.')
		print('If you are using this method and LAVA is crashing or producing whack output verify these files. A helpful guide '
			'is available in the README.')
		reference_fasta = args.f
		reference_gff = args.g
		subprocess.call('cp -p ' + reference_fasta + ' ' + new_dir + '/', shell=True)
		subprocess.call('cp -p ' + reference_gff + ' ' + new_dir + '/', shell=True)
		reference_fasta = new_dir + '/' + reference_fasta
		reference_gff = new_dir + '/' + reference_gff

	# Aligns specified fastq file with Genbank accession number and extracts majority consensus with correct annotations. 
	elif args.q != None:
		print('Using -q flag to automatically create reference from control fastq and genbank record...')
		reference_fasta, reference_gff = process(args.q, control_fastq, new_dir)

	# Exits out if users fail to provide either of the two methods: -f + -g, or -q.
	else:
		print('Improper arguments for proper reference generation. Either use -f and -g to specify a fasta and gff file '
			  'respectively, or use -q to automatically pull a genbank record to use as a reference for your control fastq')
		sys.exit(1)

	# Moves users provided files into new output directory, slow but protects original input.
	sample_path_list = []
	for sample in sample_list:
		sample_path_list.append(new_dir + '/' + sample)
		subprocess.call('cp ' + sample + ' ' + new_dir + '/',shell=True)


	# Setup complete! Now onto the pipeline.
	# Indexes the reference fasta.
	print('Indexing reference...')
	subprocess.call('bwa index ' + reference_fasta + ' 2> ' + new_dir + '/lava.log', shell=True)
	print('Done indexing.')
	subprocess.call('samtools faidx ' + reference_fasta + ' 2>> ' + new_dir + '/lava.log', shell=True)
	subprocess.call(GATK + ' CreateSequenceDictionary -R ' + reference_fasta + ' --VERBOSITY ERROR 2>> ' + new_dir + '/lava.log', shell=True)

	# For every sequence, aligns to consensus fasta and creates .bam and subsequent pileup files.
	# Removes PCR duplicates if -dedup specified, and extracts genome coverage data.
	for sample in sample_path_list:

		# Aligns sequence to consensus fasta, generating .sam files.
		print('Aligning reads for sample ' + sample + '...')
		subprocess.call('bwa mem -M -R \'@RG\\tID:group1\\tSM:' + sample + '\\tPL:illumina\\tLB:lib1\\tPU:unit1\' -p -t 6 -L [17,17] ' + 
			reference_fasta + ' ' + sample + ' > ' + sample + '.sam' + ' 2>> ' + new_dir + '/lava.log', shell=True)

		# Sorts and converts alignments from .sam to .bam.
		subprocess.call('java -jar ' + PICARD + ' SortSam INPUT=' + sample + '.sam OUTPUT=' + sample + 
			'.bam SORT_ORDER=coordinate VERBOSITY=ERROR 2>> ' + new_dir + '/lava.log', shell=True)

		# If -dedup specified, removes PCR duplicates with PICARD MarkDuplicates.
		if args.dedup:
			print('Removing PCR duplicates from sample ' + sample + '...')
			# Tags and removes duplicates using Picard.
			subprocess.call('java -jar ' + PICARD + ' MarkDuplicates INPUT=' + sample + '.bam OUTPUT=' + sample + 
				'_dedup.bam METRICS_FILE=metrics.txt VERBOSITY=ERROR REMOVE_DUPLICATES=true 2>> ' + new_dir + '/lava.log',shell=True)
			# Copies over the _dedup bam into the actual bam.
			subprocess.call('cat ' + sample + '_dedup.bam > ' + sample + '.bam' + ' 2>> ' + new_dir + '/lava.log', shell=True)
			print('Done removing PCR duplicates from sample ' + sample + '.')
		subprocess.call('java -jar ' + PICARD + ' BuildBamIndex INPUT=' + sample + '.bam VERBOSITY=ERROR 2>> ' + new_dir + '/lava.log', shell=True)
		
		# Extracts genome coverage information into separate .genomecov files for each sample.
		subprocess.call('echo sample\tposition\tcov > ' + sample + '.genomecov', shell=True)
		subprocess.call('bedtools genomecov -d -ibam ' + sample + '.bam >> ' + sample +'.genomecov' + ' 2>> ' + new_dir + '/lava.log', shell=True)
		
		# Creates pileups from .bam files.
		subprocess.call('samtools mpileup -f ' + reference_fasta + ' ' + sample + '.bam > ' + sample + '.pileup' + ' 2>> ' + 
			new_dir + '/lava.log', shell=True)

	# Checks for these files in both the current directory and one directory up, if not present print helpful error message. 
	if os.path.isfile(dir_path + '/gff3ToGenePred'):
		subprocess.call(dir_path + '/gff3ToGenePred ' + reference_gff + ' ' + new_dir + '/AT_refGene.txt -warnAndContinue -useName -allowMinimalGenes 2>> ' + new_dir 
			+ '/lava.log', shell=True)
	else:
		print('gff3ToGenePred not found - ')

	if os.path.isfile(dir_path + '/retrieve_seq_from_fasta.pl'):
		subprocess.call(dir_path + '/retrieve_seq_from_fasta.pl --format refGene --seqfile ' + reference_fasta + ' ' + new_dir + '/AT_refGene.txt --out AT_refGeneMrna.fa 2>> ' 
			+ new_dir + '/lava.log', shell=True)
	else:
		print('retrieve_seq_from_fasta.pl not found, to fix this move these files from ANNOVAR into the main lava directory. For more information and help check out the readme.')
		sys.exit(1)

	# Creates annovar /db from reference gff.
	subprocess.call('mkdir ' + new_dir + '/db/', shell=True)
	subprocess.call('mv AT_refGeneMrna.fa ' + new_dir + '/db/', shell=True)
	subprocess.call('mv ' + new_dir + '/AT_refGene.txt ' + new_dir + '/db/', shell=True)

	# Initializes merged.csv.
	subprocess.call('echo "Sample,Amino Acid Change,Position,AF,Change,Protein,NucleotideChange,LetterChange,Syn,Depth,Passage" > ' + 
		new_dir + '/merged.csv', shell=True)

	# Pulls protein information from GFF into proteins.csv. 
	# $12 = protein name, $4 = beginning nucleotide, $5 = ending nucleotide.
	# Sorts by increasing order of numbers.
	subprocess.call('grep "ID=transcript:" ' + reference_gff + ' | awk -F\'[\t;:]\' \'{print $12 "," $4 "," $5}\' | sort -t \',\' -k2 -n > ' + 
		new_dir + '/proteins.csv', shell=True)
		
	# Checks if protein list has duplicates, and prints error message.
	df = pd.read_csv(new_dir + "/proteins.csv", names = ["Protein","First", "Last"])
	col = list(df.Protein)
	for i in col:
		if col.count(i) > 1:
			sys.exit("ERROR: Duplicate protein names found. Use custom GFF file without duplicates.")


	# For all samples, extracts minor allele variants and detects amino acid changes.
	# Extracts relevant information from each sample to merged.csv: amino acid change, protein, residue number,
	# nucleotide letter change, type of mutation, depth, and metadata information.
	ref_done = False
	for x in range(0, len(sample_path_list)):
		sample = sample_path_list[x]
		passage = sample_time_list[x]

		## RYAN IS CHANGING THIS TO ALWAYS EXECUTE
		## if sample != control_fastq:
		if 1 == 1:
			print('Analyzing variants in sample ' + sample + '...')

			# Makes vcf using VarScan.
			# Current variables:
			# --validation 1 = outputs all compared positions even if non-variant (to get whole genome)
			# --min-coverage 2 = minimum coverage in both samples to call variant
			subprocess.call('java -jar ' +  VARSCAN + ' somatic ' + control_fastq + '.pileup ' + sample + '.pileup ' + sample + 
				'.vcf --validation 1 --output-vcf 1 --min-coverage 2 2>> ' + new_dir + '/lava.log', shell=True)

			# Cleans up vcf files to fix ploidy issues
			subprocess.call('mv ' + sample + '.vcf.validation ' + sample + '.vcf', shell=True)
			subprocess.call('awk -F $\'\t\' \'BEGIN {FS=OFS="\t"}{gsub("0/0","0/1",$10)gsub("0/0","1/0",$11)}1\' ' + 
							sample + '.vcf > ' + sample + '_p.vcf', shell=True)

			# Removes ambiguous calls (Y, R, M, K, S, W)
			subprocess.call('awk \'$1 ~ /^#/ {print $0;next} {if ($4 ~ /A|C|T|G/ && $5 ~ /.|A|C|T|G/) print $0}\' ' + sample + '_p.vcf > ' + sample + '_filtered.vcf', shell=True)

			# Converts .vcfs to files annovar can accept (.avinput).
			subprocess.call(dir_path + '/convert2annovar.pl -withfreq -format vcf4 -includeinfo ' + sample + '_filtered.vcf > ' + sample + '.avinput 2>> ' + new_dir + '/lava.log', shell=True)
			subprocess.call(dir_path + '/annotate_variation.pl -outfile ' + sample + ' -v -buildver AT ' + sample + '.avinput ' + new_dir + '/db/ 2>> ' + new_dir + '/lava.log', shell=True)

			# For reference sequence, grabs information (location of this information differs for reference).
			if not ref_done:

				# Grabs mutations that are SNVs with AF > 1% from the reference vcf into a separate file.
				# $18 = reference allele frequency
				if user_af=='':
					subprocess.call('awk -F":" \'($18+0)>=1{print}\' ' + sample + '.exonic_variant_function > ' + new_dir + '/ref.txt', shell=True)
				# Filters data by user specified allele frequency.
				else:
					subprocess.call('awk -v af=' + args.af + ' -F":" \'($18+0)>=' + args.af + '{print}\' ' + sample + '.exonic_variant_function > ' + new_dir + '/ref.txt', shell=True)
				subprocess.call('grep "SNV" ' + new_dir + '/ref.txt > a.tmp && mv a.tmp ' + new_dir + '/ref.txt', shell = True)

				# Grabs data from columns into csv
				# sample name (printed), protein+aa change, position, af, aa change, protein, nuc residue change,
				# nuc letter change, type of mutation, depth, 0 (because it's the reference)
				subprocess.call('awk -v ref=' + control_fastq + ' -F \'[\t:,]\' \'{print ref,","$6" "substr($9,3)","$12","$39+0","substr($9,3)","$6","substr($8,3)","substr($8,3,1)" to "substr($8,length($8))","$2","$36",0"}\' ' 
					+ new_dir + '/ref.txt > ' + new_dir + '/ref.csv', shell=True)
				
				## removed this line to attempt to fix duplication of control tabs - RCS
				## subprocess.call('cat ' + new_dir + '/ref.csv >> ' + new_dir + '/merged.csv', shell=True)
				subprocess.call('printf ' + control_fastq + '"," > ' + new_dir + '/reads.csv', shell=True)

				# Gets reads data for ref sample.
				subprocess.call('samtools flagstat ' + control_fastq + '.bam | awk \'NR==1{printf $1","} NR==5{printf $1","} NR==5{print substr($5,2)}\' '
					'>> ' + new_dir + '/reads.csv  2>> ' + new_dir + '/lava.log', shell=True)

				# Gets genome coverage data for ref sample.
				subprocess.call('echo sample	position	cov > ' + sample + '.genomecov', shell=True)
				subprocess.call('bedtools genomecov -d -ibam ' + sample + '.bam >> ' + sample + '.genomecov 2>> ' + new_dir + '/lava.log', shell=True)

				ref_done = True

			# Gets genome coverage data for each sample.
			subprocess.call('echo \'sample	position	cov\' > ' + sample + '.genomecov', shell=True)
			subprocess.call('bedtools genomecov -d -ibam ' + sample + '.bam >> ' + sample + '.genomecov 2>> ' + new_dir + '/lava.log', shell=True)

			# Gets reads data for sample.
			subprocess.call('printf ' + sample + '"," >> ' + new_dir + '/reads.csv', shell=True)
			subprocess.call('samtools flagstat ' + sample + '.bam | awk \'NR==1{printf $1","} NR==5{printf $1","} NR==5{print substr($5,2)}\' >> ' 
				+ new_dir + '/reads.csv 2>> ' + new_dir + '/lava.log', shell=True)

			# Grabs variants with AF > 1%.
			if user_af=='':
				subprocess.call('awk -F":" \'($24+0)>=1{print}\' ' + sample + '.exonic_variant_function > ' + sample + '.txt', shell=True)
			# Filters data by user specified allele frequency.
			else:
				subprocess.call('awk -v af=' + args.af + ' -F":" \'($24+0)>=' + args.af+ '{print}\' ' + sample + '.exonic_variant_function > ' + sample + '.txt', shell=True)
			
			# Grabs SNVs and stopgains/stoplosses.
			subprocess.call('grep "SNV" ' + sample + '.txt > a.tmp ', shell=True)
			subprocess.call('grep "stop" ' + sample + '.txt >> a.tmp', shell=True)
			subprocess.call('mv a.tmp ' + sample + '.txt', shell=True)

			# Gets sample metadata from metadata.csv.
			subprocess.call('SAMPLE="$(awk -F"," -v name=' + sample + ' \'$1==name {print $2}\' ' + metadata_location + ')" ', shell=True)

			# Grabs data from columns into csv: 
			# sample name (printed), protein+aa change, position, af, aa change, protein, nuc residue change,
			# nuc letter change, type of mutation, depth, sample metadata
			subprocess.call('awk -v name=' + sample + ' -v sample=$SAMPLE -F\'[\t:,]\' \'{print name","$6" "substr($9,3)","$12","$49+0","substr($9,3)","$6","substr($8,3)","substr($8,3,1)" to "substr($8,length($8))","$2","$46","sample}\' ' + 
				sample + '.txt > ' + sample + '.csv', shell = True)

			# Goes through and make sure that we have passage data in the preliminary csv files.
			add_passage(sample + '.csv',passage, new_dir)

			# Adds preliminary csv file into giant merged.csv file that combines all sample data.
			subprocess.call('cat ' + sample + '.csv >> ' + new_dir + '/merged.csv', shell=True)

			# Gets rid of "transcript"s in merged.csv.
			subprocess.call('grep -v "transcript" ' + new_dir + '/merged.csv > a.tmp && mv a.tmp ' + new_dir + '/merged.csv', shell=True)

			# Gets rid of "delins" in merged.csv
			subprocess.call('grep -v "delins" ' + new_dir + '/merged.csv > a.tmp && mv a.tmp ' + new_dir + '/merged.csv', shell=True)

	# Corrects for ribosomal slippage.
	## add_ribosomal_slippage(new_dir)

	# Gets rid of underscores in protein names.
	subprocess.call('awk \'{gsub("_"," "); print}\' ' + new_dir + '/proteins.csv > a.tmp && mv a.tmp ' + new_dir + '/proteins.csv', shell=True)

	# Removes intermediate files (default behavior), unless otherwise specified by -save.
	print('Cleaning up...')
	if not args.save:
		clean_up(new_dir)

	# Detects and prints warning if more than 5000 lines in final output folder.
	num_lines = sum(1 for line in open(new_dir + '/merged.csv'))
	if num_lines > 5000:
		print('WARNING: greater than 5,000 lines detected in the final output folder. This may cause browsers to crash when using depth or allele frequncy sliders. If this happens you may need to perform more QC on your reads or on the merged.csv file.')
	
	for line in open(new_dir + '/lava.log'):
		if 'WARNING: Unable to retrieve regions at' in line and 'due to lack of sequence information' in line:
			print('LAVA has crashed due to improper formatting of gff/fasta pair. Make sure that your fasta header has the same name as the first column of the gff file. For more information look at the gff formatting guide on the README')
			sys.exit()
		if 'A total of' in line and 'sequences will be ignored due to lack of correct ORF annotation' in line:
			print('Some (or all) of your feature annotations have failed to improper formatting of gff/fasta pair. Make sure that your gff file is formatted as specified in the README. Also check to make sure that your gff protein locations are correct and contain valid open reading frames. If using an online reference, make sure it is correct, references that have different protein lengths than your sample will produce this error.')

		if 'Warning: skipping: invalid strand' in line:
			print('One (or more) of the features on your gff file is improperly formatted, please check the gff guide on the readme for more help. Other coding features should still be annotated correctly.')
	# Pipeline is done! Now on to the visualization.
	print('Generating visualization...')

	# Calls visualizer genome_protein_plots.py depending on what directory the user is in.
	## merged, proteins, reads 
	if os.path.isfile(dir_path + '/genome_protein_plots.py'):
		subprocess.call('python3 ' + dir_path + '/genome_protein_plots.py ' + nuc_flag + ' ' + png_flag + ' ' + user_af + ' ' +
			new_dir + '/merged.csv ' + new_dir + '/proteins.csv ' + new_dir + '/reads.csv ' + new_dir + ' ' + plot_title, shell=True)
	elif os.path.isfile('genome_protein_plots.py'):
		subprocess.call('python3 genome_protein_plots.py ' + nuc_flag + ' ' + png_flag + ' ' + user_af + ' ' +
			+ new_dir + '/merged.csv ' + new_dir + '/proteins.csv ' + new_dir + '/reads.csv ' + new_dir+ ' ' + plot_title, shell=True)
	else:
		print('Genome_protein_plots could not be found. Output will not be visualized - go to the README for help.')
		sys.exit(1)
