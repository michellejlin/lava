#!/bin/bash
echo Example1_file1.fastq

	if [[ "`basename /Users/uwvirongs/Documents/Michelle/lava/test_data/Example1_file1.fastq`" == "`basename Example1_file1.fastq`" ]]
	then 
		echo `basename /Users/uwvirongs/Documents/Michelle/lava/test_data/Example1_file1.fastq` found
		echo first file found ending process execution for Example1_file1.fastq
		exit 0

	else 
		echo "not first file"
fi

echo "continuing execution for ${R1}"

echo 'sample	position	cov' > Example1_file1.fastq.genomecov 

/usr/local/miniconda/bin/bedtools genomecov -d -ibam FIRSTFILE.bam.bam >> Example1_file1.fastq.genomecov

# reads.csv from all processes will be merged together at end 
 printf Example1_file1.fastq"," > reads.csv

/usr/local/miniconda/bin/samtools flagstat FIRSTFILE.bam.bam | 	awk 'NR==1{printf $1","} NR==5{printf $1","} NR==5{print substr($5,2)}' >> reads.csv

awk -F":" '($24+0)>=1{print}' blank.exonic_variant_function.samp> Example1_file1.fastq.txt

 grep "SNV" Example1_file1.fastq.txt > a.tmp
 grep "stop" Example1_file1.fastq.txt >> a.tmp
 mv a.tmp Example1_file1.fastq.txt

SAMPLE="$(awk -F"," -v name=Example1_file1.fastq '$1==name {print $2}' Example1_metadata.csv)"

echo $SAMPLE

 awk -v name=Example1_file1.fastq -v sample=0 -F'[	:,]' '{print name","$6" "substr($9,3)","$12","$46+0","substr($9,3)","$6","substr($8,3)","substr($8,3,1)" to "substr($8,length($8))","$2","$43","sample}' Example1_file1.fastq.txt > Example1_file1.fastq.csv
