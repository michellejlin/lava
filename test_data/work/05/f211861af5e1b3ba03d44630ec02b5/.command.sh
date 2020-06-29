#!/bin/bash

ls -latr

echo Analyzing variants in sample Example1_file2.fastq

# here for file passthrough (input -> output)
mv Example1_file2.fastq.bam Example1_file2.fastq.bam.bam 

java -jar /usr/local/bin/VarScan somatic CONTROL.fastq.pileup Example1_file2.fastq.pileup Example1_file2.fastq.vcf --validation 1 --output-vcf 1 --min-coverage 2

mv Example1_file2.fastq.vcf.validation Example1_file2.fastq.vcf

awk -F $'	' 'BEGIN {FS=OFS="	"}{gsub("0/0","0/1",$10)gsub("0/0","1/0",$11)gsub("1/1","0/1",$10)gsub("1/1","1/0",$11)}1' Example1_file2.fastq.vcf > Example1_file2.fastq_p.vcf

file="Example1_file2.fastq""_p.vcf"


sed 's/NC_039477.1/lava/g' Example1_file2.fastq_p.vcf > test.txt 

mv test.txt Example1_file2.fastq_p.vcf

#convert2annovar.pl -withfreq -format vcf4 -includeinfo Example1_file2.fastq_p.vcf > Example1_file2.fastq.avinput 
convert2annovar.pl -withfreq -format vcf4old -includeinfo Example1_file2.fastq_p.vcf > Example1_file2.fastq.avinput 

annotate_variation.pl -outfile Example1_file2.fastq -v -buildver AT Example1_file2.fastq.avinput .


ls -lah 

if [[ "`basename /Users/uwvirongs/Documents/Michelle/lava/test_data/Example1_file1.fastq`" == "`basename Example1_file2.fastq`" ]]
	then 
		echo `basename /Users/uwvirongs/Documents/Michelle/lava/test_data/Example1_file1.fastq` found
		touch blank.exonic_variant_function.samp

	else 
		echo "not first file"
		echo `basename Example1_file2.fastq` `basename /Users/uwvirongs/Documents/Michelle/lava/test_data/Example1_file1.fastq`
		mv Example1_file2.fastq.exonic_variant_function Example1_file2.fastq.exonic_variant_function.samp
fi
