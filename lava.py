# lava Version 0.9rs 
# Longitudinal Analysis of Viral Alleles 

import subprocess 
import argparse
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio import Entrez
import shutil
import sys
import timeit
import time
from datetime import datetime
import re
import os.path

Entrez.email = 'uwvirongs@gmail.com'

VERSION = 'v0.9rs'

def check_picard():
	if os.path.isfile('./picard.jar'):
		return './picard.jar'
	elif os.path.isfile('../picard.jar'):
		return '../picard.jar'
	else:
		print('Picard not found - lava is being executed from : ')
		subprocess.call('pwd', shell=True)
		print('LAVA checked for picard in the above folder and the main lava folder')
		print('to fix this error download picard and unzip it into the main lava directory - for more indepth help check out the readme')
		sys.exit(1)

def check_gatk():
	if os.path.isfile('./gatk-4.0.11.0/gatk'):
		return './gatk-4.0.11.0/gatk'
	elif os.path.isfile('../gatk-4.0.11.0/gatk'):
		return '../gatk-4.0.11.0/gatk'
	else:
		print('GATK not found - lava is being executed from : ')
		subprocess.call('pwd', shell=True)
		print('LAVA checked for GATK in the above folder and the main lava folder')
		print('to fix this error download GATK and unzip it into the main lava directory - for more indepth help check out the readme')
		sys.exit(1)

def check_varscan():
	if os.path.isfile('./VarScan'):
		return './VarScan'
	elif os.path.isfile('../VarScan'):
		return '../VarScan'
	else:
		print('VarScan not found - lava is being executed from : ')
		subprocess.call('pwd', shell=True)
		print('LAVA checked for VarScan in the above folder and the main lava folder')
		print('to fix this error download VarScan and unzip it into the main lava directory NOTE: the jar file needs to be named VarScan - for more indepth help check out the readme')
		sys.exit(1)


# Takes a file path pointing to a fasta file and returns two lists, the first is a list of all the fasta headers and the second is 
# A list of all the sequence for each record, these will always be in the same order. Also replaces Us with Ts
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


# Takes two MAFFT aligned sequences and returns two arrays of numbers coresponding to relative genome position for each
# any position that doesn't align will have a -1 
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


# takes two number arrays built by build_num_arrays and adjusts a given nucleotide position from a reference to a given genome 
# returns adjusted index as a string. If the adjustment fails (for example off the edges of the alignments) a -1 will be returned
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

# Takes a location to a two column csv metadata sheet with names of each sample we want to analyze and time data. Returns two lists
# One of sample/file names and the other of their relative temporal relationship
def read_metadata(filepath):
	sample_list = []
	sample_time_list = []
	first = True
	for line in open(filepath):
		if first:
			first = False
		# protects against empty lines
		elif ',' in line:
			print(line)
			sample_list.append(line.split(',')[0])
			sample_time_list.append(line.split(',')[1])

	if not os.path.isfile(sample_list[0]):
		print('METADATA file is incorecctly formatted!')
		print(sample_list[0] + ' was the first sample but could not be found in the working direcotry')
		print('Please move your samples to the directory you are running lava from, this is:')
		subprocess.call('pwd', shell=True)
		sys.exit(1)
	return sample_list, sample_time_list


# Takes a fastq file, GenBank Accession number and output directory. Aligns fastq file to genbank acsession, and extracts majority consensus as well as 
# transfers annotations from Genbank to newly extracted consensus. Returns a filepath to a fasta of the majority consensus and a filepath to a coresponding 
# .gff file containing protein annotations. Also writes these files into the output directory as well as all intermediate files 
def process(ref_seq_gb, fastq, new_dir):

	# Pull whole genbank reference from Entrez
	record = Entrez.read(Entrez.esearch(db='nucleotide', term=ref_seq_gb))
	h2 = Entrez.efetch(db='nucleotide', id=record["IdList"][0], rettype='gb', retmode='text')
	e = open(new_dir + '/lava_ref.gbk', 'w')
	e.write(h2.read())
	e.close()

	# pull reference fasta from Entrez
	h = Entrez.efetch(db='nucleotide', id=record["IdList"][0], rettype='fasta', retmode='text')
	d = open(new_dir + '/lava_ref.fasta', 'w')
	d.write(h.read())
	d.close()

	# index reference fasta
	subprocess.call('bwa index ' + new_dir + '/lava_ref.fasta 2> ' + new_dir + '/lava.log', shell=True)

	# align passed fastq reads to the downloaded reference fasta create vcf file and index it, sterr is sent to lava.log 
	subprocess.call('bwa mem -M ' + new_dir + '/lava_ref.fasta ' + fastq + ' > ' + new_dir + '/aln.sam 2>> ' + new_dir + '/lava.log', shell=True)
	subprocess.call('samtools view -S -b ' + new_dir + '/aln.sam > ' + new_dir + '/aln.bam 2>> '+ new_dir + '/lava.log', shell=True)
	subprocess.call('samtools sort ' + new_dir + '/aln.bam -o ' + new_dir + '/aln.sorted.bam 2>> ' + new_dir + '/lava.log', shell=True)
	subprocess.call('bcftools mpileup -Ou -f ' + new_dir + '/lava_ref.fasta ' + new_dir + '/aln.sorted.bam | bcftools call -mv -Oz -o ' + 
		new_dir + '/calls.vcf.gz 2>> ' + new_dir + '/lava.log', shell=True)
	subprocess.call('tabix ' + new_dir + '/calls.vcf.gz 2>> ' + new_dir + '/lava.log', shell=True)
	# create the consensus and write it to the output folder 
	subprocess.call('cat ' + new_dir + '/lava_ref.fasta | bcftools consensus ' + new_dir + '/calls.vcf.gz > ' + new_dir + '/consensus.fasta', shell=True)


	# rewrite fasta with lava as the sequence header so that annovar can match this with the reference .gff
	fasta = new_dir + '/consensus.fasta'
	ignore, genome_ref = read_fasta(fasta)
	g = open(fasta, 'w')
	g.write('>lava\n')
	g.write(genome_ref[0])
	g.close()

	gene_loc_list = []
	gene_product_list = []
	allow_one = False

	# pull CDS annotations from downloaded genbank file 
	for line in open(new_dir + '/lava_ref.gbk'):
		if ' CDS ' in line:
			gene_loc_list.append(re.findall(r'\d+', line))
			allow_one = True

		if '/product="' in line and allow_one:
			allow_one = False
			px = line.split('=')[1][1:-2]
			px = px.split()[0]
			gene_product_list.append(px)

	# write a new file for mafft input 
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

	# use maaft to align consensus fasta with downloaded fasta 
	nullo, genome = read_fasta(fasta)
	subprocess.call('mafft --quiet ' + new_dir + '/aligner.fasta > ' + new_dir + '/lava.ali', shell=True)
	ali_list, ali_genomes = read_fasta(new_dir + '/lava.ali')
	ref_seq = ali_genomes[1]
	our_seq = ali_genomes[0]
	our_seq_num_array, ref_seq_num_array = build_num_arrays(our_seq, ref_seq)

	# go through every annotation and adjust the coordinates from the reference to the consensus sequence
	for entry in range(0, len(gene_loc_list)):
		for y in range(0, len(gene_loc_list[entry])):
			gene_loc_list[entry][y] = adjust(int(gene_loc_list[entry][y]), our_seq_num_array, ref_seq_num_array, genome)   

	# write a new .gff file for annovar to use later 
	g = open(new_dir + '/lava_ref.gff', 'w')
	g.write('##gff-version 3\n##source-version geneious 9.1.7\n')
	name= 'lava'
	for x in range(0, len(gene_product_list)):
		g.write(name + '\tGeneious\tgene\t' + gene_loc_list[x][0] + '\t' + gene_loc_list[x][1] + '\t.\t+\t.\tID=gene:' + 
			gene_product_list[x] + ';biotype=protein_coding\n')
		g.write(name + '\tGeneious\tCDS\t' +  gene_loc_list[x][0] + '\t' + gene_loc_list[x][1] + '\t.\t+\t0\tID=CDS:' + 
			gene_product_list[x] + ';Parent=transcript:' + gene_product_list[x] + ';biotype=protein_coding\n')
		g.write(name + '\tGeneious\ttranscript\t' + gene_loc_list[x][0] + '\t' + gene_loc_list[x][1] + '\t.\t+\t.\tID=transcript:' + 
			gene_product_list[x] + ';Parent=gene:' + gene_product_list[x] + ';biotype=protein_coding\n')
	# return the filepaths for the new fasta and .gff file as strings 
	return fasta, new_dir + '/lava_ref.gff'

# goes through the merged.csv file and makes sure that the temporal metadata has been added correctly to the output
# also scans for complex mutations and prints a helpfull warning message to the terminal listing what samples and the changes
# also re-writes the file with complex mutations being annotated as 'complex'
def add_passage(sample, passage):
	complex_list = []
	line_list = []
	last = ''
	last_line = ''
	add_a_complex = False
	for line in open(sample):
		# make sure we only read lines that contain data 
		if ',' in line:
			# if the amino acid reside is the same as the last one we read than we have a complex mutaiton, this assumes that merged.csv is in sorted order
			# which we garuntee that it is
			if line.split(',')[4][:-1] == last:
				print(line)
				print('WARNING: Complex mutation detected (multiple nucleotide changes within the same codon)! Sample=' + line.split(',')[0] + 
				  	' Change1=' + line.split(',')[4] + ' and Change2=' + last_line.split(',')[4])
				complex_list.append(last_line)
				add_a_complex = True

			last = line.split(',')[4][:-1]
			
			# if we don't have time data add it here 
			if line.strip().split(',')[10] == '':
				line = line.strip() + str(passage) + '\n'
			if add_a_complex:
				complex_list.append(line)
				add_a_complex = False

			line_list.append(line)
			last_line = line

	g = open(sample, 'w')

	# re-write file with 'complex' as change type 
	for line in line_list:
		if line in complex_list:
			# re-write line but with complex instead of whatever change was previously listed 
			g.write(','.join(line.split(',')[:8]) + ',complex,' + ','.join(line.split(',')[9:]))
		else:
			g.write(line)
	g.close()


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Version ' + VERSION + '\nLongitudinal analysis of minor variants across whole viral'
												 ' genomes. Use either -g and -f together or use -q to annotate control fastq')
	parser.add_argument('-f', help='Specify a reference fasta with the majority consensus of the control fastq. This option must be used '
								   'with the -g flag to specify the protein annotations relative to the start of this fasta.')
	parser.add_argument('-g', help='Specify a refernce gff file with the protein annotations for the reference fasta supplied with the -f flag.'
								   ' This option must be paired with the -f flag')
	parser.add_argument('-q', help='Provide a Genbank accession number. This record will be used to generate a majority consensus from the '
								   'control fastq, this consensus will be annotated from the downloaded genbank record as well')
	parser.add_argument('-nuc', action='store_true', help='Results are listed as nucleotide changes not amino acid changes. Do not use with -png')
	parser.add_argument('-png', action='store_true', help='Output results as a png. Do not use with -nuc')
	parser.add_argument('-dedup', action='store_true', help='Optional flag, will perform automatic removal of PCR duplicates via DeDup')
	parser.add_argument('control_fastq', help='Required argument: The fastq reads for the first sample in your longitudinal analysis')
	parser.add_argument('metadata', help='Required argument: A two column csv the first column is the name of all the fastqs you wish '
										 'to include in your analysis. All fastqs that you want to include need to be specified in this '
										 'file AND be located in the folder from which you are running lava. The second column is the '
										 'temporal seperation between the samples.  This is unitless so you can input passage number, '
										 'days, or whatever condition your experiment happens to have.')
	parser.add_argument('-o', help='Optional flag to name the output folder that lava will stuff output into. If a name isn\'t provided '
								   'folder will be named lava-date')

	# check for argument sanity
	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)

	# if an output directory is not provided pull the date and time in isoformat with no colons
	if args.o != None:
		new_dir = args.o
	else:
		new_dir = str(datetime.now().isoformat())
		new_dir = '_'.join(new_dir.split(':'))

	# create flags for passing to genome_protein_plots.py 
	if args.nuc:
		nuc_flag = '-nuc'
	else:
		nuc_flag = ''

	if args.png:
		png_flag = '-png'
	else:
		png_flag = ''

	# check for picard, gatk, and varscan exit if we can't find them 
	PICARD = check_picard()
	GATK = check_gatk()
	VARSCAN = check_varscan()
	# write new output folder 
	subprocess.call('mkdir -p ' + new_dir, shell=True)
	# debuggin console print
	subprocess.call('pwd', shell=True)
	subprocess.call('pwd', shell=True)
	metadata_location = args.metadata
	control_fastq = args.control_fastq

	sample_list, sample_time_list = read_metadata(metadata_location)
	# debugging print
	print(sample_list)
	print(sample_time_list)
	# copy fastqs to new directory - slow but protects user data from the insanity that is about to occur 
	subprocess.call('cp ' + control_fastq + ' ' + new_dir + '/', shell=True)
	control_fastq = new_dir + '/' + control_fastq

	# make sure that we've got a way of pulling annotations, if the user gives -f -g and -q then we only use -f and -g 
	if args.f != None and args.g != None:
		print('Using -f and -g flags to annotate control fastq, -q flag will be ignored.')
		print('This method of reference generation assumes that your fasta and .gff file are formated correctly')
		print('If you are using this method and lava is crashing or producing whack output verify these files. A helpful guide is avalible in the README')
		reference_fasta = args.f
		reference_gff = args.g
		subprocess.call('cp ' + reference_fasta + ' ' + new_dir + '/', shell=True)
		subprocess.call('cp ' + reference_gff + ' ' + new_dir + '/', shell=True)
		reference_fasta = new_dir + '/' + reference_fasta
		reference_gff = new_dir + '/' + reference_gff

	elif args.q != None:
		print('Using -q flag to automatically create reference from control fastq and genbank record.')
		reference_fasta, reference_gff = process(args.q, control_fastq, new_dir)

	else:
		print('Improper arguments for proper reference generation. Either use -f and -g to specify a fasta and gff file, '
			  'respectively or use -q to automatically pull a genbank record to use as a reference for your control fastq')
		sys.exit(1)

	
	# move users provided files into new output directory, slow but protects original input 
	sample_path_list = []
	for sample in sample_list:
		sample_path_list.append(new_dir + '/' + sample)
		subprocess.call('cp ' + sample + ' ' + new_dir + '/',shell=True)
	

	print('Indexing reference...')
	subprocess.call('bwa index ' + reference_fasta + ' 2> ' + new_dir + '/lava.log', shell=True)
	print('Done indexing.')

	subprocess.call('samtools faidx ' + reference_fasta + ' 2>> ' + new_dir + '/lava.log', shell=True)
	subprocess.call(GATK + ' CreateSequenceDictionary -R ' + reference_fasta + ' --VERBOSITY ERROR 2>> ' + new_dir + '/lava.log', shell=True)
	# debuggin print command 
	print(sample_path_list)
	for sample in sample_path_list:

		print('Aligning reads for sample ' + sample)
		subprocess.call('bwa mem -M -R \'@RG\\tID:group1\\tSM:' + sample + '\\tPL:illumina\\tLB:lib1\\tPU:unit1\' -p -t 6 -L [17,17] ' + 
			reference_fasta + ' ' + sample + ' > ' + sample + '.sam' + ' 2>> ' + new_dir + '/lava.log', shell=True)

		subprocess.call('java -jar ' + PICARD + ' SortSam INPUT=' + sample + '.sam OUTPUT=' + sample + 
			'.bam SORT_ORDER=coordinate VERBOSITY=ERROR 2>> ' + new_dir + '/lava.log', shell=True)

		if args.dedup:
			print('Removing PCR duplicates from sample ' + sample)
			subprocess.call('java -jar ' + PICARD + ' MarkDuplicates INPUT=' + sample + '.bam OUTPUT=' + sample + 
				'_dedup.bam METRICS_FILE=metrics.txt VERBOSITY=ERROR 2>> ' + new_dir + '/lava.log',shell=True)
			subprocess.call('cat ' + name + '_dedup.bam > ' + name + '.bam' + ' 2>> ' + new_dir + '/lava.log', shell=True)
			print('done removing PCR duplicates from sample ' + sample)

		subprocess.call('java -jar ' + PICARD + ' BuildBamIndex INPUT=' + sample + '.bam VERBOSITY=ERROR 2>> ' + new_dir + '/lava.log', shell=True)
		subprocess.call('echo sample\tposition\tcov > ' + sample + '.genomecov', shell=True)
		subprocess.call('bedtools genomecov -d -ibam ' + sample + '.bam >> ' + sample +'.genomecov' + ' 2>> ' + new_dir + '/lava.log', shell=True)
		subprocess.call('samtools mpileup -f ' + reference_fasta + ' ' + sample + '.bam > ' + sample + '.pileup' + ' 2>> ' + 
			new_dir + '/lava.log', shell=True)
	# create annovar /db from reference gff 
	subprocess.call('gff3ToGenePred ' + reference_gff + ' ' + new_dir + '/AT_refGene.txt -warnAndContinue -useName -allowMinimalGenes 2>> ' + new_dir 
		+ '/lava.log', shell=True)
	# check for these files in both the current directory and one directory up, if not present print helpful error message 
	if os.path.isfile('../retrieve_seq_from_fasta.pl'):
		subprocess.call('../retrieve_seq_from_fasta.pl --format refGene --seqfile ' + reference_fasta + ' ' + new_dir + '/AT_refGene.txt --out AT_refGeneMrna.fa 2>> ' 
			+ new_dir + '/lava.log', shell=True)
	elif os.path.isfile('./retrieve_seq_from_fasta.pl'):
		subprocess.call('./retrieve_seq_from_fasta.pl --format refGene --seqfile ' + reference_fasta + ' ' + new_dir + '/AT_refGene.txt --out AT_refGeneMrna.fa 2>> ' 
			+ new_dir + '/lava.log', shell=True)
	else:
		print('retrieve_seq_from_fasta.pl not found, to fix this move these files from ANNOVAR into the main lava directory. For more information and help check out the readme')
		sys.exit(1)

	subprocess.call('mkdir ' + new_dir + '/db/', shell=True)
	subprocess.call('mv AT_refGeneMrna.fa ' + new_dir + '/db/', shell=True)
	subprocess.call('mv ' + new_dir + '/AT_refGene.txt ' + new_dir + '/db/', shell=True)
	# TODO: Re-write these in python 
	# intialzie merged.csv 
	subprocess.call('echo "Sample,Amino Acid Change,Position,AF,Change,Protein,NucleotideChange,LetterChange,Syn,Depth,Passage" > ' + 
		new_dir + '/merged.csv', shell=True)

	# pull annotations for plotting 
	subprocess.call('grep "ID=transcript:" ' + reference_gff + ' | awk -F\'[\t;:]\' \'{print $12 "," $4 "," $5}\' | sort -t \',\' -k2 -n > ' + 
		new_dir + '/proteins.csv', shell=True)

	ref_done = False
	for x in range(0, len(sample_path_list)):
		sample = sample_path_list[x]
		passage = sample_time_list[x]
		# RYAN IS CHANGING THIS TO ALWAYS EXECUTE
		#if sample != control_fastq:
		if 1 == 1:
			print('Analyzing variants in sample ' + sample)
			subprocess.call('java -jar ' +  VARSCAN + ' somatic ' + control_fastq + '.pileup ' + sample + '.pileup ' + sample + 
				'.vcf --validation 1 --output-vcf 1 --min-coverage 2 2>> ' + new_dir + '/lava.log', shell=True)


			subprocess.call('mv ' + sample + '.vcf.validation ' + sample + '.vcf', shell=True)
			subprocess.call('awk -F $\'\t\' \'BEGIN {FS=OFS="\t"}{gsub("0/0","0/1",$10)gsub("0/0","1/0",$11)}1\' ' + 
							sample + '.vcf > ' + sample + '_p.vcf', shell=True)

			subprocess.call('../convert2annovar.pl -withfreq -format vcf4 -includeinfo ' + sample + '_p.vcf > ' + sample + '.avinput 2>> ' + new_dir + '/lava.log', shell=True)
			#annotates all mutations with codon changes
			subprocess.call('../annotate_variation.pl -outfile ' + sample + ' -v -buildver AT ' + sample + '.avinput ' + new_dir + '/db/ 2>> ' + new_dir + '/lava.log', shell=True)

			if not ref_done:
				subprocess.call('awk -F":" \'($18+0)>=5{print}\' ' + sample + '.exonic_variant_function > ' + new_dir + '/ref.txt', shell=True)
				subprocess.call('grep "SNV" ' + new_dir + '/ref.txt > a.tmp && mv a.tmp ' + new_dir + '/ref.txt', shell = True)
				subprocess.call('awk -v ref=' + control_fastq + ' -F \'[\t:,]\' \'{print ref,","$6" "substr($9,3)","$12","$39+0","substr($9,3)","$6","substr($8,3)","substr($8,3,1)" to "substr($8,length($8))","$2","$36",0"}\' ' 
					+ new_dir + '/ref.txt > ' + new_dir + '/ref.csv', shell=True)
				# removed this line to attempt to fix duplication of control tabs - RCS
				#subprocess.call('cat ' + new_dir + '/ref.csv >> ' + new_dir + '/merged.csv', shell=True)
				subprocess.call('printf ' + control_fastq + '"," > ' + new_dir + '/reads.csv', shell=True)

				subprocess.call('samtools flagstat ' + control_fastq + '.bam | awk \'NR==1{printf $1","} NR==5{printf $1","} NR==5{print substr($5,2)}\' '
					'>> ' + new_dir + '/reads.csv  2>> ' + new_dir + '/lava.log', shell=True)
				subprocess.call('echo sample	position	cov > ' + sample + '.genomecov', shell=True)
				subprocess.call('bedtools genomecov -d -ibam ' + sample + '.bam >> ' + sample + '.genomecov 2>> ' + new_dir + '/lava.log', shell=True)

				ref_done = True


			subprocess.call('echo \'sample	position	cov\' > ' + sample + '.genomecov', shell=True)
			subprocess.call('bedtools genomecov -d -ibam ' + sample + '.bam >> ' + sample + '.genomecov 2>> ' + new_dir + '/lava.log', shell=True)
			subprocess.call('printf ' + sample + '"," >> ' + new_dir + '/reads.csv', shell=True)
			subprocess.call('samtools flagstat ' + sample + '.bam | awk \'NR==1{printf $1","} NR==5{printf $1","} NR==5{print substr($5,2)}\' >> ' 
				+ new_dir + '/reads.csv 2>> ' + new_dir + '/lava.log', shell=True)
			subprocess.call('awk -F":" \'($24+0)>=5{print}\' ' + sample + '.exonic_variant_function > ' + sample + '.txt', shell=True)

			subprocess.call('grep "SNV" ' + sample + '.txt > a.tmp ', shell=True)
			subprocess.call('grep "stop" ' + sample + '.txt >> a.tmp', shell=True)
			subprocess.call('mv a.tmp ' + sample + '.txt', shell=True)
			subprocess.call('SAMPLE="$(awk -F"," -v name=' + sample + ' \'$1==name {print $2}\' ' + metadata_location + ')" ', shell=True)
			subprocess.call('awk -v name=' + sample + ' -v sample=$SAMPLE -F\'[\t:,]\' \'{print name","$6" "substr($9,3)","$12","$49+0","substr($9,3)","$6","substr($8,3)","substr($8,3,1)" to "substr($8,length($8))","$2","$46","sample}\' ' + 
				sample + '.txt > ' + sample + '.csv', shell = True)
			# go through and make sure that we have passage data in the preliminary csv files 
			add_passage(sample + '.csv',passage )
			subprocess.call('cat ' + sample + '.csv >> ' + new_dir + '/merged.csv', shell=True)

	print('generating Vizualisation')
	if os.path.isfile('ngls_test.html'):
		shutil.copy('ngls_test.html', new_dir)
	elif os.path.isfile('../ngls_test.html'):
		subprocess.call('cp ../ngls_test.html .', shell=True)
	else:
		print('ngls_test.html file not found - 3d crystal structure will not be integrated')

	# merged, proteins, reads 
	if os.path.isfile('../genome_protein_plots.py'):
		subprocess.call('pythonw ../genome_protein_plots.py ' + new_dir + '/merged.csv ' + new_dir + '/proteins.csv ' + new_dir + '/reads.csv ' + 
			new_dir + ' ' + nuc_flag + ' ' + png_flag, shell=True)
	elif os.path.isfile('genome_protein_plots.py'):
		subprocess.call('pythonw genome_protein_plots.py ' + new_dir + '/merged.csv ' + new_dir + '/proteins.csv ' + new_dir + '/reads.csv ' + 
			new_dir+ ' ' + nuc_flag + ' ' + png_flag, shell=True)
	else:
		print('Genome_protein_plots could not be found. Output will not be visualized - go to XXXX for help')
		sys.exit(1)



