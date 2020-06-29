#!/bin/bash

set -e 

   #for logging

echo Example1_ref.fasta

   ls -latr 

   #Entrez fetch function

if [[ Example1_ref.fasta == "NO_FILE" ]]
	then
		python3 /Users/uwvirongs/Documents/Michelle/lava/bin/pull_entrez.py False
	else 
		mv Example1_ref.fasta lava_ref.fasta
		mv Example1_ref.gff lava_ref.gff

		#Creates empty txt file
		touch ribosomal_start.txt
fi


echo Example1_file1.fastq.gz

   /usr/local/miniconda/bin/bwa index lava_ref.fasta

   /usr/local/miniconda/bin/bwa mem -t !{task.cpus} -M lava_ref.fasta Example1_file1.fastq.gz | /usr/local/miniconda/bin/samtools view -Sb - > aln.bam

/usr/local/miniconda/bin/samtools sort aln.bam -o aln.sorted.bam 

   /usr/local/miniconda/bin/bcftools mpileup --max-depth 500000 -P 1.1e-100 -Ou -f lava_ref.fasta aln.sorted.bam | /usr/local/miniconda/bin/bcftools call -m -Oz -o calls.vcf.gz 

   /usr/local/miniconda/bin/tabix calls.vcf.gz

   gunzip calls.vcf.gz 

   /usr/local/miniconda/bin/bcftools filter -i '(DP4[0]+DP4[1]) < (DP4[2]+DP4[3]) && ((DP4[2]+DP4[3]) > 0)' calls.vcf -o calls2.vcf

   /usr/local/miniconda/bin/bgzip calls2.vcf

   /usr/local/miniconda/bin/tabix calls2.vcf.gz 

   cat lava_ref.fasta | /usr/local/miniconda/bin/bcftools consensus calls2.vcf.gz > consensus.fasta

if [[ Example1_ref.fasta == "NO_FILE" ]]
	then
		python3 /Users/uwvirongs/Documents/Michelle/lava/bin/write_gff.py
fi


 # Avoiding filename colision during run_pipeline process 
 mv Example1_file1.fastq.gz CONTROL.fastq
