#sets all the paths
export PICARD=/Users/uwvirongs/downloads/picard-2.18.7/picard/build/libs/picard.jar
export GATK=/Users/uwvirongs/downloads/gatk-4.0.5.1/gatk
export VARSCAN=/Users/uwvirongs/Downloads/VarScan.v2.3.9.jar
export PATH=$PATH:/Users/uwvirongs/Downloads/annovar 
export SNPEFFPATH=/Users/uwvirongs/snpEff/
export PATH="/Users/uwvirongs/downloads/annovar/:$PATH"
dedup=false
script_path=$(dirname "$0")

while getopts ":f::q::g::h::d" opt; do
  case $opt in
    f)
    ref=$OPTARG
      if ! [[ $OPTARG =~ \.fasta$ ]]
      then
        echo 'ERROR: Not a fasta file. Please ensure file is in valid format and try again.'
        exit 1
      elif ! [ -f $OPTARG ]
      then
        echo 'ERROR: File does not exist.'
        exit 1
      fi
      ;;
    q)
      query=$OPTARG
      echo $query "was the query." >&2
      ;;
    g)
    gff=$OPTARG
    if ! [[ $OPTARG =~ \.gff$ ]]
      then
        echo 'ERROR: Not a gff file. Please ensure file is in valid format and try again.'
        exit 1
      elif ! [ -f $OPTARG ]
      then
        echo 'ERROR: File does not exist.'
        exit 1
      fi
      ;;
    h)
      echo "Usage: lava.sh [options] control.fastq"
      echo "	-h	Display this help message."
      echo "	-q	User query for searching Genbank for suitable reference."
      echo "			NOT to be mixed with -g. Recommended way of getting gff file."
      echo "	-g	User-provided .gff file pulled off Geneious."
      echo "			Only for specialized cases. NOT to be mixed with -q."
      echo "	-f	User-provided .fasta file. Only for specialized cases."
      echo "			Intended for use with user-provided .gff file."
      exit 0
      ;;
    d)
    dedup=true
    echo "lava is running in dedup mode, PCR duplicates will be removed"
      ;;
    \?)
      echo "ERROR: Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "ERROR: Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

control=${@:OPTIND:1}

if test -z "$control"
then
	echo "ERROR: Control fastq not given. Please provide fastq at the end of the command and try again."
	exit 1
fi 

if [[ ! -z "$query" && ! -z "$gff" ]]
then
	echo "ERROR: Please only specify one method of generating a gff."
	exit 1
fi

if [[ -z "$query" && -z "$gff" ]]
then
	echo "ERROR: Please specify at least one method of generating a gff. The flags -q and -g cannot be mixed."
	echo "-q is highly recommended for automatic generation of a gff. Please try again."
	exit 1
fi

if [[ -z "$query" && ! (( ! -z "$gff" && ! -z "$ref" )) ]]
then
	echo "ERROR: Please specify both a gff file and a fasta file if providing those files manually."
	exit 1
fi 

if [[ ! -z "$query" ]]
then
	python ../process.py $query $ref
	gff="lava_ref.gff"
	echo "Created fasta and gff files."
else
	gff_check=$(grep "Parent=gene:" $gff)
	if [[ -z "$gff_check" ]]
	then
		echo "Processing geneious gff..."
		awk 'BEGIN{FS=OFS="\t"} $3=="cds" {$3="CDS"}1;' $gff > a.tmp && mv a.tmp $gff
		awk 'BEGIN{FS=OFS="\t"} $3=="CDS" {$8=0}1;' $gff > a.tmp && mv a.tmp $gff
		awk 'BEGIN{IGNORECASE=1} {gsub(/mrna/,"transcript")}1' $gff > a.tmp && mv a.tmp $gff	
		awk -F" |=" 'BEGIN{OFS=":"}{t=$2; $2=$3; $3=t; print;}' $gff > a.tmp && mv a.tmp $gff
		awk '{gsub(/Name=/,"ID=")}1' $gff > a.tmp && mv a.tmp $gff
		awk 'BEGIN{FS=OFS="\t"} $3=="transcript" {$9=$9";Parent=gene:REPLACEME;biotype=protein_coding"}1' $gff > a.tmp && mv a.tmp $gff
		awk 'BEGIN{FS=OFS="\t"} $3=="CDS" {$9=$9";Parent=transcript:REPLACEME;biotype=protein_coding"}1' $gff > a.tmp && mv a.tmp $gff
		grep -v "regulatory" $gff > a.tmp && mv a.tmp $gff
		awk -F":|;" '{gsub(/REPLACEME/,$2)}1' $gff > a.tmp && mv a.tmp $gff
		awk 'BEGIN{FS=OFS="\t"} $3=="gene" {$9=$9";biotype=protein_coding"}1' $gff > a.tmp && mv a.tmp $gff
		sed "1s/.*/##gff-version 3/" $gff > a.tmp && mv a.tmp $gff
		sed "2s/.*/##source-version geneious 9.1.7/" $gff > a.tmp && mv a.tmp $gff
	fi
fi

echo "ref is "$ref
echo "gff is "$gff

# replace Us with Ts
sed -i -e 's/U/T/g' $ref

#indexes $reference
bwa index $ref

#indexes $reference for GATK use
samtools faidx $ref
$GATK CreateSequenceDictionary -R $ref --VERBOSITY WARNING

for sample in *.fastq
do
	name=$(basename "$sample" .fastq)
	
	#aligns with updaed end error masking penalty - RCS
	bwa mem -M -R '@RG\tID:group1\tSM:'"$name"'\tPL:illumina\tLB:lib1\tPU:unit1' -p -t 6 -L [17,17] $ref $sample > $name.sam

	#converts sam to bam and sorts
	java -jar $PICARD SortSam \
		INPUT=$name.sam \
		OUTPUT=$name.bam \
		SORT_ORDER=coordinate \
		VERBOSITY=WARNING
		
	# optional flag
	if $dedup; then
		echo 'Removing PCR duplicates'
		java -jar $PICARD MarkDuplicates \
			INPUT=$name.bam \
			OUTPUT=$name'_dedup'.bam \
			METRICS_FILE=metrics.txt \
			VERBOSITY=WARNING
			
		cat $name'_dedup'.bam > $name.bam
	fi

	#indexes bam file
	java -jar $PICARD BuildBamIndex INPUT=$name.bam VERBOSITY=WARNING
	# Generate by nt coverage
	echo 'sample	position	cov' > $name.genomecov
	bedtools genomecov -d -ibam $name.bam >> $name.genomecov
	
	#creates pileup
	samtools mpileup -f $ref $name.bam > $name.pileup
done

#converts the gff into something annovar can use
gff3ToGenePred $gff AT_refGene.txt -warnAndContinue -useName -allowMinimalGenes
retrieve_seq_from_fasta.pl --format refGene --seqfile $ref AT_refGene.txt --out AT_refGeneMrna.fa
mkdir db/
mv AT_refGeneMrna.fa db/
mv AT_refGene.txt db/

#created merged header
echo "Sample,Amino Acid Change,Position,AF,Change,Protein,NucleotideChange,LetterChange,Syn,Depth,Passage" > merged.csv
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
		java -jar $VARSCAN somatic $con.pileup $name.pileup $name.vcf --validation 1 --output-vcf 1 --min-coverage 2
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
			awk -v ref=$con -F '[\t:,]' '{print ref,","$6" "substr($9,3)","$12","$39+0","substr($9,3)","$6","substr($8,3)","substr($8,3,1)" to "substr($8,length($8))","$2","$36",0"}' ref.txt > ref.csv
			cat ref.csv >> merged.csv
			printf $con"," > reads.csv
			samtools flagstat $con.bam | awk 'NR==1{printf $1","} NR==5{printf $1","} NR==5{print substr($5,2)}' >> reads.csv
			# Coverage Generation
			echo 'sample	position	cov' > $name.genomecov
			bedtools genomecov -d -ibam $name.bam >> $name.genomecov
			
			REF_DONE=true
		fi
		printf $name"," >> reads.csv
		
		# More coverage code, one of these is most likely reduntant 
		echo 'sample	position	cov' > $name.genomecov
		bedtools genomecov -d -ibam $name.bam >> $name.genomecov
		
		samtools flagstat $name.bam | awk 'NR==1{printf $1","} NR==5{printf $1","} NR==5{print substr($5,2)}' >> reads.csv
		
		awk -F":" '($24+0)>=5{print}' $name.exonic_variant_function >$name.txt 
		grep "SNV" $name.txt > a.tmp 
		grep "stop" $name.txt >> a.tmp
		mv a.tmp $name.txt 
		#gets Passage number from metadata.csv 
		SAMPLE="$(awk -F"," -v name=$name '$1==name {print $2}' metadata.csv)" 
		awk -v name=$name -v sample=$SAMPLE -F'[\t:,]' '{print name","$6" "substr($9,3)","$12","$49+0","substr($9,3)","$6","substr($8,3)","substr($8,3,1)" to "substr($8,length($8))","$2","$46","sample}' $name.txt > $name.csv 
		cat $name.csv >> merged.csv
	fi
done	

$script_path/genome_protein_plots.py

# name=$(basename "$control" .fastq)
# mv genome_protein_plots.html $name'_genome_protein_plots'.html
