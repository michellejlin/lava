#!/bin/bash

echo "Sample,Amino Acid Change,Position,AF,Change,Protein,NucleotideChange,LetterChange,Syn,Depth,Passage" > merged.csv

python3 /Users/uwvirongs/Documents/Michelle/lava/bin/initialize_merged_csv.py


# Creating pileup for control fastq here 

/usr/local/miniconda/bin/bwa mem -t !{task.cpus}  -M -R '@RG\tID:group1\tSM:CONTROL.fastq\tPL:illumina\tLB:lib1\tPU:unit1' -p -t 6 -L [17,17] consensus.fasta CONTROL.fastq > CONTROL.fastq.sam

java -jar /usr/bin/picard.jar SortSam INPUT=CONTROL.fastq.sam OUTPUT=CONTROL.fastq.bam SORT_ORDER=coordinate VERBOSITY=ERROR 

if false 
	then
		echo Deduplicating CONTROL.fastq
		java -jar /usr/bin/picard.jar MarkDuplicates INPUT=CONTROL.fastq.bam OUTPUT=CONTROL.fastq_dedup.bam METRICS_FILE=metrics.txt VERBOSITY=ERROR REMOVE_DUPLICATES=true
		cat CONTROL.fastq_dedup.bam > CONTROL.fastq.bam
fi

java -jar /usr/bin/picard.jar BuildBamIndex INPUT=CONTROL.fastq.bam VERBOSITY=ERROR

echo sample	position	cov > CONTROL.fastq.genomecov

/usr/local/miniconda/bin/bedtools genomecov -d -ibam  CONTROL.fastq.bam >> CONTROL.fastq.genomecov

/usr/local/miniconda/bin/bcftools mpileup -f consensus.fasta CONTROL.fastq.bam > CONTROL.fastq.pileup
