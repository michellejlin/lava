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


PICARD='/Users/uwvirongs/downloads/picard-2.18.7/picard/build/libs/picard.jar'
GATK='/Users/uwvirongs/downloads/gatk-4.0.5.1/gatk'
VARSCAN='/Users/uwvirongs/Downloads/VarScan.v2.3.9.jar'
ANNOVAR='/Users/uwvirongs/Downloads/annovar'

Entrez.email = 'uwvirongs@gmail.com'

VERSION = 'v0.9rs'


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

# Takes a location to a two column csv metadata sheet with names of each sample we want to analyze and time data. Returns two lists
# One of sample/file names and the other of their relative temporal relationship
def read_metadata(filepath):
	sample_list = []
	sample_time_list = []
	first = True
	for line in open(filepath):
		if first:
			first = False
		else:
			sample_list.append(line.split(',')[0])
			sample_time_list.append(line.split(',')[1])

	if not os.path.isfile(sample_list[0]):
		print('METADATA file is incorecctly formatted!')
		print(sample_list[0] + ' was the first sample but could not be found in the working direcotry')
		print('Please move your samples to the directory you are running lava from, this is:')
		subprocess.call('pwd', shell=True)
		sys.exit(1)
	return sample_list, sample_time_list


# Essentially import the entire process.py file into this function 
def process(ref_seq_gb, fastq, new_dir):

	record = Entrez.read(Entrez.esearch(db='nucleotide', term=ref_seq_gb))
	h2 = Entrez.efetch(db='nucleotide', id=record["IdList"][0], rettype='gb', retmode='text')
	e = open(new_dir + '/lava_ref.gbk', 'w')
	e.write(h2.read())
	e.close()

	h = Entrez.efetch(db='nucleotide', id=record["IdList"][0], rettype='fasta', retmode='text')
	d = open(new_dir + '/lava_ref.fasta', 'w')
	d.write(h.read())
	d.close()


	subprocess.call('bwa index ' + new_dir + '/lava_ref.fasta 2> ' + new_dir + '/lava.log', shell=True)
	subprocess.call('bwa mem -M ' + new_dir + '/lava_ref.fasta ' + fastq + ' > ' + new_dir + '/aln.sam 2>> ' + new_dir + '/lava.log', shell=True)
	subprocess.call('samtools view -S -b ' + new_dir + '/aln.sam > ' + new_dir + '/aln.bam 2>> '+ new_dir + '/lava.log', shell=True)
	subprocess.call('samtools sort ' + new_dir + '/aln.bam -o ' + new_dir + '/aln.sorted.bam 2>> ' + new_dir + '/lava.log', shell=True)
	subprocess.call('bcftools mpileup -Ou -f ' + new_dir + '/lava_ref.fasta ' + new_dir + '/aln.sorted.bam | bcftools call -mv -Oz -o ' + 
		new_dir + '/calls.vcf.gz 2>> ' + new_dir + '/lava.log', shell=True)
	subprocess.call('tabix ' + new_dir + '/calls.vcf.gz 2>> ' + new_dir + '/lava.log', shell=True)
	subprocess.call('cat ' + new_dir + '/lava_ref.fasta | bcftools consensus ' + new_dir + '/calls.vcf.gz > ' + new_dir + '/consensus.fasta', shell=True)

	fasta = new_dir + '/consensus.fasta'
	ignore, genome_ref = read_fasta(fasta)
	g = open(fasta, 'w')
	g.write('>lava\n')
	g.write(genome_ref[0])
	g.close()

	gene_loc_list = []
	gene_product_list = []
	allow_one = False

	for line in open(new_dir + '/lava_ref.gbk'):
		if ' CDS ' in line:
			gene_loc_list.append(re.findall(r'\d+', line))
			allow_one = True

		if '/product="' in line and allow_one:
			allow_one = False
			px = line.split('=')[1][1:-2]
			px = px.split()[0]
			gene_product_list.append(px)

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

	nullo, genome = read_fasta(fasta)
	subprocess.call('mafft --quiet ' + new_dir + '/aligner.fasta > ' + new_dir + '/lava.ali', shell=True)
	ali_list, ali_genomes = read_fasta(new_dir + '/lava.ali')
	ref_seq = ali_genomes[1]
	our_seq = ali_genomes[0]
	our_seq_num_array, ref_seq_num_array = build_num_arrays(our_seq, ref_seq)

	for entry in range(0, len(gene_loc_list)):
		for y in range(0, len(gene_loc_list[entry])):
			gene_loc_list[entry][y] = adjust(int(gene_loc_list[entry][y]), our_seq_num_array, ref_seq_num_array, genome)   

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
	return fasta, new_dir + '/lava_ref.gff'

def add_passage(sample, passage):

	line_list = []
	last = ''
	last_line = ''
	for line in open(sample):
		if line.split(',')[4][:-1] == last:
			print('WARNING: Complex mutation detected (multiple nucleotide changes within the same codon)! Sample=' + line.split(',')[0] + ' Change1=' + line.split(',')[4] + ' and Change2=' + last_line.split(',')[4])
		last = line.split(',')[4][:-1]
		last_line = line
		if line.strip().split(',')[10] == '':
			line = line.strip() + str(passage)
		line_list.append(line)

	g = open(sample, 'w')
	for line in line_list:
		g.write(line)
	g.close()


if __name__ == '__main__':

	start_time = timeit.default_timer()

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

	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(0)

	if args.o != None:
		new_dir = args.o
	else:
		new_dir = str(datetime.now().isoformat())
		new_dir = '_'.join(new_dir.split(':'))

	if args.nuc:
		nuc_flag = '-nuc'
	else:
		nuc_flag = ''

	if args.png:
		png_flag = '-png'
	else:
		png_flag = ''
		
	subprocess.call('mkdir -p ' + new_dir, shell=True)

	metadata_location = args.metadata
	control_fastq = args.control_fastq

	sample_list, sample_time_list = read_metadata(metadata_location)

	subprocess.call('cp ' + control_fastq + ' ' + new_dir + '/', shell=True)
	control_fastq = new_dir + '/' + control_fastq

	if args.f != None and args.g != None:
		print('Using -f and -g flags to annotate control fastq, -q flag will be ignored.')
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

	

	sample_path_list = []
	for sample in sample_list:
		sample_path_list.append(new_dir + '/' + sample)

		subprocess.call('cp ' + sample + ' ' + new_dir + '/', shell=True)
	

	print('Indexing reference...')
	subprocess.call('bwa index ' + reference_fasta + ' 2> ' + new_dir + '/lava.log', shell=True)
	print('Done indexing.')

	subprocess.call('samtools faidx ' + reference_fasta + ' 2>> ' + new_dir + '/lava.log', shell=True)
	subprocess.call(GATK + ' CreateSequenceDictionary -R ' + reference_fasta + ' --VERBOSITY ERROR 2>> ' + new_dir + '/lava.log', shell=True)

	for sample in sample_path_list:

		print('Aligning reads for sample ' + sample)
		subprocess.call('bwa mem -M -R \'@RG\\tID:group1\\tSM:' + sample + '\\tPL:illumina\\tLB:lib1\\tPU:unit1\' -p -t 6 -L [17,17] ' + 
			reference_fasta + ' ' + sample + ' > ' + sample + '.sam' + ' 2>> ' + new_dir + '/lava.log', shell=True)

		subprocess.call('java -jar ' + PICARD + ' SortSam INPUT=' + sample + '.sam OUTPUT=' + sample + 
			'.bam SORT_ORDER=coordinate VERBOSITY=ERROR 2>> ' + new_dir + '/lava.log', shell=True)

		if args.dedup:
			subprocess.call('java -jar ' + PICARD + ' MarkDuplicates INPUT=' + sample + '.bam OUTPUT=' + sample + 
				'_dedup.bam METRICS_FILE=metrics.txt VERBOSITY=ERROR 2>> ' + new_dir + '/lava.log',shell=True)
			subprocess.call('cat ' + name + '_dedup.bam > ' + name + '.bam' + ' 2>> ' + new_dir + '/lava.log', shell=True)

		subprocess.call('java -jar ' + PICARD + ' BuildBamIndex INPUT=' + sample + '.bam VERBOSITY=ERROR 2>> ' + new_dir + '/lava.log', shell=True)
		subprocess.call('echo sample\tposition\tcov > ' + sample + '.genomecov', shell=True)
		subprocess.call('bedtools genomecov -d -ibam ' + sample + '.bam >> ' + sample +'.genomecov' + ' 2>> ' + new_dir + '/lava.log', shell=True)
		subprocess.call('samtools mpileup -f ' + reference_fasta + ' ' + sample + '.bam > ' + sample + '.pileup' + ' 2>> ' + 
			new_dir + '/lava.log', shell=True)

	subprocess.call('gff3ToGenePred ' + reference_gff + ' ' + new_dir + '/AT_refGene.txt -warnAndContinue -useName -allowMinimalGenes 2>> ' + new_dir 
		+ '/lava.log', shell=True)
	if os.path.isfile('../retrieve_seq_from_fasta.pl'):
		subprocess.call('../retrieve_seq_from_fasta.pl --format refGene --seqfile ' + reference_fasta + ' ' + new_dir + '/AT_refGene.txt --out AT_refGeneMrna.fa 2>> ' 
			+ new_dir + '/lava.log', shell=True)
	elif os.path.isfile('./retrieve_seq_from_fasta.pl'):
		subprocess.call('./retrieve_seq_from_fasta.pl --format refGene --seqfile ' + reference_fasta + ' ' + new_dir + '/AT_refGene.txt --out AT_refGeneMrna.fa 2>> ' 
			+ new_dir + '/lava.log', shell=True)
	else:
		print('retrieve_seq_from_fasta.pl not found, to fix this move these files from ANNOVAR into the main lava directory')
		sys.exit(1)

	subprocess.call('mkdir ' + new_dir + '/db/', shell=True)
	subprocess.call('mv AT_refGeneMrna.fa ' + new_dir + '/db/', shell=True)
	subprocess.call('mv ' + new_dir + '/AT_refGene.txt ' + new_dir + '/db/', shell=True)
	# TODO: Re-write these in python 
	subprocess.call('echo "Sample,Amino Acid Change,Position,AF,Change,Protein,NucleotideChange,LetterChange,Syn,Depth,Passage" > ' + 
		new_dir + '/merged.csv', shell=True)
	subprocess.call('grep "ID=transcript:" ' + reference_gff + ' | awk -F\'[\t;:]\' \'{print $12 "," $4 "," $5}\' | sort -t \',\' -k2 -n > ' + 
		new_dir + '/proteins.csv', shell=True)

	ref_done = False
	for x in range(0, len(sample_path_list)):
		sample = sample_path_list[x]
		passage = sample_time_list[x]
		if sample != control_fastq:
			subprocess.call('java -jar ' +  VARSCAN + ' somatic ' + control_fastq + '.pileup ' + sample + '.pileup ' + sample + 
				'.vcf --validation 1 --output-vcf 1 --min-coverage 2 2>> ' + new_dir + '/lava.log', shell=True)


			subprocess.call('mv ' + sample + '.vcf.validation ' + sample + '.vcf', shell=True)
			subprocess.call('awk -F $\'\t\' \'BEGIN {FS=OFS="\t"}{gsub("0/0","0/1",$10)gsub("0/0","1/0",$11)}1\' ' + 
							sample + '.vcf > ' + sample + '_p.vcf', shell=True)

			subprocess.call('export PATH=$PATH:/Users/uwvirongs/Downloads/annovar;export PATH="/Users/uwvirongs/downloads/annovar/:$PATH";convert2annovar.pl '
				'-withfreq -format vcf4 -includeinfo ' + sample + '_p.vcf > ' + sample + '.avinput 2>> ' + new_dir + '/lava.log', shell=True)
			#annotates all mutations with codon changes
			subprocess.call('export PATH=$PATH:/Users/uwvirongs/Downloads/annovar;export PATH="/Users/uwvirongs/downloads/annovar/:$PATH";annotate_variation.pl '
				'-outfile ' + sample + ' -v -buildver AT ' + sample + '.avinput ' + new_dir + '/db/ 2>> ' + new_dir + '/lava.log', shell=True)

			if not ref_done:
				subprocess.call('awk -F":" \'($18+0)>=5{print}\' ' + sample + '.exonic_variant_function > ' + new_dir + '/ref.txt', shell=True)
				subprocess.call('grep "SNV" ' + new_dir + '/ref.txt > a.tmp && mv a.tmp ' + new_dir + '/ref.txt', shell = True)
				subprocess.call('awk -v ref=' + control_fastq + ' -F \'[\t:,]\' \'{print ref,","$6" "substr($9,3)","$12","$39+0","substr($9,3)","$6","substr($8,3)","substr($8,3,1)" to "substr($8,length($8))","$2","$36",0"}\' ' 
					+ new_dir + '/ref.txt > ' + new_dir + '/ref.csv', shell=True)
				subprocess.call('cat ' + new_dir + '/ref.csv >> ' + new_dir + '/merged.csv', shell=True)
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



