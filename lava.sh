#!/bin/bash

ref=$1
control=$2
gff=$3

#sets all the paths, horribly
export PICARD=/Users/uwvirongs/downloads/picard-2.18.7/picard/build/libs/picard.jar
export GATK=/Users/uwvirongs/downloads/gatk-4.0.5.1/gatk
export VARSCAN=/Users/uwvirongs/Downloads/VarScan.v2.3.9.jar
export PATH=$PATH:/Users/uwvirongs/Downloads/annovar 
export SNPEFFPATH=/Users/uwvirongs/snpEff/
export RATT_HOME=/Users/uwvirongs/Downloads/ratt-code/
export PATH="/Users/uwvirongs/downloads/annovar/:$PATH"
export PATH="/Users/uwvirongs/downloads/MAC-master/:$PATH"
export PATH="/Users/uwvirongs/Downloads/MUMmer3.232/scripts/:$PATH"

# esearch -db nucleotide -query "parainfluenza 1 complete genome refseq" -sort relevance| efetch -format fasta | sed -n '1,/>/p' | sed \$d > test.fasta
# reference=$(awk 'NR==1{print substr($1,2)}' test.fasta)
# esearch -db nucleotide -query $reference | efetch -format gb > test.gb
# python gff_generation.py 5S5.fasta test.fasta test.gb --no_spell_check

# gzip test.gb
# genbank2gff3.pl test.gb.gz --split
# esearch -db nucleotide -query "parainfluenza 1 complete genome refseq" -sort relevance | efetch -format gb | sed -n '1,/LOCUS/p' | sed \$d | sed \$d > test.gb
# genbank2embl.pl test.gb test.embl
# $RATT_HOME/start.ratt.sh ./ 5S5.fasta ref.fasta Assembly
# 
# replace Us with Ts
# sed -i -e 's/U/T/g' $ref
# 
# #indexes $reference
# bwa index $ref
# 
# #indexes $reference for GATK use
# samtools faidx $ref
# $GATK CreateSequenceDictionary -R $ref --VERBOSITY WARNING
# 
# for sample in *.fastq
# do
# 	name=$(basename "$sample" .fastq)
# 	
# 	aligns
# 	bwa mem -M -R '@RG\tID:group1\tSM:'"$name"'\tPL:illumina\tLB:lib1\tPU:unit1' -p -t 6 $ref $sample > $name.sam
# 
# 	converts sam to bam and sorts
# 	java -jar $PICARD SortSam \
# 		INPUT=$name.sam \
# 		OUTPUT=$name.bam \
# 		SORT_ORDER=coordinate \
# 		VERBOSITY=WARNING
# 
# 	gets rid of pcr duplicates
# 	java -jar $PICARD MarkDuplicates \
# 		INPUT=$name.bam \
# 		OUTPUT=$name'_dedup'.bam \
# 		METRICS_FILE=metrics.txt \
# 		VERBOSITY=WARNING
# 
# 	indexes bam file
# 	java -jar $PICARD BuildBamIndex INPUT=$name'_dedup'.bam VERBOSITY=WARNING
# 	
# 	creates pileup
# 	samtools mpileup -f $ref $name'_dedup'.bam > $name'_dedup'.pileup
# done

# samtools mpileup -uf test.fasta $name'_dedup'.bam | bcftools call -c | vcfutils.pl vcf2fq > consensus.fq
# seqtk seq -A consensus.fq > consensus.fasta

# Change the .gff file here to make it acceptable
# make sure none of the other things have more annotations than the transcript

# awk 'BEGIN{FS=OFS="\t"} $3=="cds" {$3="CDS"}1;' $gff > a.tmp && mv a.tmp $gff
# awk 'BEGIN{FS=OFS="\t"} $3=="CDS" {$8=0}1;' $gff > a.tmp && mv a.tmp $gff
# awk 'BEGIN{IGNORECASE=1} {gsub(/mrna/,"transcript")}1' $gff > a.tmp && mv a.tmp $gff	
# awk -F" |=" 'BEGIN{OFS=":"}{t=$2; $2=$3; $3=t; print;}' $gff > a.tmp && mv a.tmp $gff
# awk '{gsub(/Name:/,"ID=")}1' $gff > a.tmp && mv a.tmp $gff
# awk 'BEGIN{FS=OFS="\t"} $3=="transcript" {$9=$9";Parent=gene:REPLACEME;biotype=protein_coding"}1' $gff > a.tmp && mv a.tmp $gff
# awk 'BEGIN{FS=OFS="\t"} $3=="CDS" {$9=$9";Parent=transcript:REPLACEME;biotype=protein_coding"}1' $gff > a.tmp && mv a.tmp $gff
# grep -v "regulatory" $gff > a.tmp && mv a.tmp $gff
# awk -F":|;" '{gsub(/REPLACEME/,$2)}1' $gff > a.tmp && mv a.tmp $gff
# awk 'BEGIN{FS=OFS="\t"} $3=="gene" {$9=$9";biotype=protein_coding"}1' $gff > a.tmp && mv a.tmp $gff
# sed "1s/.*/##gff-version 3/" $gff > a.tmp && mv a.tmp $gff
# sed "2s/.*/##source-version geneious 9.1.7/" $gff > a.tmp && mv a.tmp $gff

#converts the gff into something annovar can use
gff3ToGenePred $gff AT_refGene.txt -warnAndContinue -useName -allowMinimalGenes
retrieve_seq_from_fasta.pl --format refGene --seqfile $ref AT_refGene.txt --out AT_refGeneMrna.fa
mkdir db/
mv AT_refGeneMrna.fa db/
mv AT_refGene.txt db/

#created merged header
echo "Sample,Amino Acid Change,Position,AF,Change,Protein,NucleotideChange,Syn,Depth,Passage" > merged.csv
#extract out protein names and sequences
grep "ID=transcript:" $gff | awk -F'[\t;:]' '{print $12 "," $4 "," $5}' | sort -t ',' -k2 -n > proteins.csv

REF_DONE=false

for sample in *.fastq
do
	name=$(basename "$sample" .fastq)
	con=$(basename "$control" .fastq)
	if [ $name != $con ]
	then
		#creates vcf of all bases with reference and alternate alleles
		java -jar $VARSCAN somatic $con'_dedup'.pileup $name'_dedup'.pileup $name.vcf --validation 1 --output-vcf 1 --min-coverage 2
		mv $name.vcf.validation $name.vcf
		awk -F $'\t' 'BEGIN {FS=OFS="\t"}{gsub("0/0","0/1",$10)gsub("0/0","1/0",$11)}1' $name.vcf > $name_p.vcf
		#annotates all mutations with codon changes
		convert2annovar.pl -withfreq -format vcf4 -includeinfo $name_p.vcf>$name.avinput
		annotate_variation.pl -outfile $name -v -buildver AT $name.avinput db/
		if [ $REF_DONE = false ]
		then
			#grep -v '0:0%' $name.exonic_variant_function | awk -F":" '($18+0)>=5{print}' > ref.txt
			awk -F":" '($18+0)>=5{print}' $name.exonic_variant_function > ref.txt
			grep "SNV" ref.txt > a.tmp && mv a.tmp ref.txt
			awk -v ref=$con -F '[\t:,]' '{print ref,","$6" "substr($9,3)","$12","$39+0","substr($9,3)","$6","substr($8,3,1)" to "substr($8,length($8))","$2","$36",0"}' ref.txt > ref.csv
			cat ref.csv >> merged.csv
			printf $con"," > reads.csv
			samtools flagstat $con'_dedup'.bam | awk 'NR==1{printf $1","} NR==5{printf $1","} NR==5{print substr($5,2)}' >> reads.csv
			REF_DONE=true
		fi
		printf $name"," >> reads.csv
		samtools flagstat $name'_dedup'.bam | awk 'NR==1{printf $1","} NR==5{printf $1","} NR==5{print substr($5,2)}' >> reads.csv
		
		awk -F":" '($24+0)>=5{print}' $name.exonic_variant_function >$name.txt 
		grep "SNV" $name.txt > a.tmp 
		grep "stop" $name.txt >> a.tmp
		mv a.tmp $name.txt 
		#gets Passage number from metadata.csv 
		SAMPLE="$(awk -F"," -v name=$name '$1==name {print $2}' metadata.csv)" 
		awk -v name=$name -v sample=$SAMPLE -F'[\t:,]' '{print name","$6" "substr($9,3)","$12","$49+0","substr($9,3)","$6","substr($8,3,1)" to "substr($8,length($8))","$2","$46","sample}' $name.txt > $name.csv 
		cat $name.csv >> merged.csv
	fi
done	

genome_protein_plots

#Graveyard :(





# mkdir $SNPEFFPATH/data/$gff2
# cp $gff2.gff $SNPEFFPATH/data/$gff2
# cp $gff2.fa $SNPEFFPATH/data/$gff2
# cd $SNPEFFPATH/data/$gff2
# mv $gff2.gff genes.gff
# mv $gff2.fa sequences.fa
# cd $SNPEFF_PATH
# echo "$gff2.genome : Para1" >> snpEff.config
# echo "$gff2.genome : Para1" >> snpEff.config
# java -Xmx4g -jar $SNPEFFPATH/snpEff.jar build -gff3 -v $gff2
# 
# bcftools csq -f $gff.fasta -g genes.gff -l -p m test2.vcf > test -o out.bcf
# bgzip test.vcf
# tabix test.vcf.gz
# bcftools norm -m-both -o test2.vcf test.vcf.gz
# bcftools norm -f $gff.fasta -o test3.vcf test2.vcf

#region-based
#annotate_variation.pl -regionanno -outfile test -dbtype gff3 -gff3dbfile $gff -v -buildver AF test.avinput ./

#Doesn't seem to be filtering for AF -> gets blank output
#$GATK SelectVariants -R testref.fasta -V added.vcf -select "AF > 0.05" --select-type-to-exclude MNP -O addedtest.vcf

#no difference between changing ploidy and not
#$GATK Mutect2 -R testref.fasta -I control_dedup.bam -I cancer_dedup.bam -tumor sample1 -normal sample1 -O test2.vcf -ploidy 100	

#zipping and indexing for bcf tools
#bgzip -c added.vcf >added.vcf.gz
#tabix -p vcf added.vcf.gz

#calls variants, germline, messes up allele frequencies
#$GATK HaplotypeCaller -R para1.fasta \
#	-I test_dedup.bam \
#	-O test_variants.vcf
#	--sample-ploidy 1

#allele frequencies are messed up
#samtools index cancer.bam
#freebayes -f testref.fasta cancer_dedup.bam > out.vcf -p 1 --pooled-continuous --min-alternate-fraction 0.05 --haplotype-length 0 --min-alternate-count 1

#invalid sequence genome??
#library(deepSNV())
#regions <- data.frame(chr="", start, end)
#deepSNV("test.bam", "cancer.bam", regions=regions, q=quality)

##WORKING CODE FOR CALLING ONLY VARIANTS

# 	calls variants
# 	$GATK Mutect2 \
# 		-R $ref \
# 		-I $name'_dedup'.bam \
# 		-tumor $name \
# 		-O $name.vcf \
# 		--verbosity WARNING
# 		
# 	$GATK Mutect2 -R testref.fasta -I control_	dedup.bam -I cancer_dedup.bam -tumor sample1 -normal sample1 -O test.vcf
# 	maybe --output-mode EMIT_ALL_SITES
# 
# 	only takes allele frequencies > 5% (takes multiallelic sites separately because VariantFiltration can't handle them)
# 	$GATK SelectVariants \
# 		-R $ref \
# 		-V $name.vcf \
# 		-O $name'_split'.vcf \
# 		--restrict-alleles-to BIALLELIC \
# 		--verbosity WARNING
# 	$GATK VariantFiltration \
# 		-R $ref \
# 		-V $name.vcf \
# 		-O $name'_filtered'.vcf \
# 		--genotype-filter-expression "AF < 0.05" \
# 		--genotype-filter-name "Low" \
# 		--verbosity WARNING
# 		
# 	grep -v "Low" $name'_filtered'.vcf >> $name'_final'.vcf 
# 
# 	$GATK SelectVariants \
# 		-R $ref \
# 		-V $name.vcf \
# 		-O $name'_split'.vcf \
# 		--restrict-alleles-to MULTIALLELIC \
# 		--verbosity WARNING
# 		
# 	grep "^[^#;]" $name'_split'.vcf >> $name'_final'.vcf
# 
# 	prints to nicerish table
# 	$GATK VariantsToTable \
# 		-R $ref \
# 		-V $name'_final'.vcf \
# 		-F CHROM -F POS -F REF -F ALT -F QUAL -GF AF \
# 		-O $name'_final'.txt \
# 		--verbosity WARNING
# 	
# 	replaces useless CHROM tag with sample name
# 	awk -v var="$name" -F $'\t' 'BEGIN {FS=OFS="\t"}{$1=var;print}' $name'_final'.txt >> combined_final.txt