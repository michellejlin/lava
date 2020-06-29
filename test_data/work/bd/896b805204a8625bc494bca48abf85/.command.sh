#!/bin/bash

echo aligning "Example1_file1.fastq"


/usr/local/miniconda/bin/bwa mem -t 10 -M -R '@RG\tID:group1\tSM:Example1_file1.fastq\tPL:illumina\tLB:lib1\tPU:unit1' -p -t 10 -L [17,17] consensus.fasta Example1_file1.fastq > Example1_file1.fastq.sam


java -jar /usr/bin/picard.jar SortSam INPUT=Example1_file1.fastq.sam OUTPUT=Example1_file1.fastq.bam SORT_ORDER=coordinate VERBOSITY=ERROR 

if false 
	then
		echo "Deduplicating Example1_file1.fastq"
		java -jar /usr/bin/picard.jar MarkDuplicates INPUT=${R1}.bam OUTPUT=${R1}_dedup.bam METRICS_FILE=metrics.txt VERBOSITY=ERROR REMOVE_DUPLICATES=true
		cat ${R1}_dedup.bam > ${R1}.bam
fi

java -jar /usr/bin/picard.jar BuildBamIndex INPUT=Example1_file1.fastq.bam VERBOSITY=ERROR

echo sample	position	cov > Example1_file1.fastq.genomecov

/usr/local/miniconda/bin/bedtools genomecov -d -ibam  Example1_file1.fastq.bam >> Example1_file1.fastq.genomecov

/usr/local/miniconda/bin/bcftools mpileup -f consensus.fasta Example1_file1.fastq.bam > Example1_file1.fastq.pileup


	if [[ "`basename /Users/uwvirongs/Documents/Michelle/lava/test_data/Example1_file1.fastq`" == "`basename Example1_file1.fastq`" ]]
	then 
		echo `basename /Users/uwvirongs/Documents/Michelle/lava/test_data/Example1_file1.fastq` found
		mv Example1_file1.fastq.bam FIRSTFILE.bam

	else 
		echo "not first file"
fi
