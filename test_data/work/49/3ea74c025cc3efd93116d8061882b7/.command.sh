#!/bin/bash

echo Example1_file1.fastq
echo NO_VAL

if [[ "NO_VAL" == "NO_VAL" ]]
	then
		awk -F":" '($18+0)>=1{print}' Example1_file1.fastq.exonic_variant_function > ref.txt

	else	
		awk -v af=NO_VAL -F":" '($18+0)>=NO_VAL{print}' Example1_file1.fastq.exonic_variant_function > ref.txt
fi

grep "SNV" ref.txt > a.tmp && mv a.tmp ref.txt 

awk -v ref=CONTROL.fastq -F '[	:,]' 	'{print ref,","$6" "substr($9,3)","$12","$39+0","substr($9,3)","$6","substr($8,3)","substr($8,3,1)" 	to "substr($8,length($8))","$2","$36",0"}' ref.txt > ref.csv

printf Example1_file1.fastq"," >> reads.csv

/usr/local/miniconda/bin/samtools flagstat FIRSTFILE.bam | 	awk 'NR==1{printf $1","} NR==5{printf $1","} NR==5{print substr($5,2)}' >> reads.csv

 echo sample	position	cov > Example1_file1.fastq.genomecov
 /usr/local/miniconda/bin/bedtools genomecov -d -ibam FIRSTFILE.bam >> Example1_file1.fastq.genomecov



if [[ "NO_VAL" == "NO_VAL" ]]
	then
		awk -F":" '($24+0)>=1{print}' Example1_file1.fastq.exonic_variant_function > Example1_file1.fastq.txt 
	else
		awk -v af=NO_VAL -F":" '($24+0)>=NO_VAL{print}' Example1_file1.fastq.exonic_variant_function > Example1_file1.fastq.txt 
fi

 grep "SNV" Example1_file1.fastq.txt > a.tmp
 grep "stop" Example1_file1.fastq.txt >> a.tmp
 mv a.tmp Example1_file1.fastq.txt

SAMPLE="$(awk -F"," -v name=Example1_file1.fastq '$1==name {print $2}' Example1_metadata.csv)"

 awk -v name=Example1_file1.fastq -v sample=0 -F'[	:,]' '{print name","$6" "substr($9,3)","$12","$46+0","substr($9,3)","$6","substr($8,3)","substr($8,3,1)" to "substr($8,length($8))","$2","$43","sample}' Example1_file1.fastq.txt > Example1_file1.fastq.csv
