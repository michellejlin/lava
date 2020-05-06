/*
 * Define the parameters used in the processes below
 */

/*
 * Define the processes used in this workflow
 */




process CreateGFF { 
    
    container "quay.io/vpeddu/lava_image:latest"

    //errorStrategy 'retry'
    //maxRetries 3
    // Define the input files
    input:
      file PULL_ENTREZ
      val(GENBANK)
      file CONTROL_FASTQ
      file MAFFT_PREP
      file GFF_WRITE
    // Define the output files
    output: 
      file "lava_ref.gbk"
      file "lava_ref.fasta"
      file "consensus.fasta"
      file "lava_ref.gff"
	  file "CONTROL.fastq"

	  //file "*.pileup"
    // Code to be executed inside the task
    script:
    """
    #!/bin/bash
    
    #for logging
    
    ls -latr 

    #Entrez fetch function
    python3 ${PULL_ENTREZ} ${GENBANK}

	echo ${CONTROL_FASTQ}

    /usr/local/miniconda/bin/bwa index lava_ref.fasta

    /usr/local/miniconda/bin/bwa mem -M lava_ref.fasta ${CONTROL_FASTQ} | /usr/local/miniconda/bin/samtools view -Sb - > aln.bam

	/usr/local/miniconda/bin/samtools sort aln.bam -o aln.sorted.bam 

    /usr/local/miniconda/bin/bcftools mpileup --max-depth 500000 -P 1.1e-100 -Ou -f lava_ref.fasta aln.sorted.bam | /usr/local/miniconda/bin/bcftools call -m -Oz -o calls.vcf.gz 

    /usr/local/miniconda/bin/tabix calls.vcf.gz

    gunzip calls.vcf.gz 

    /usr/local/miniconda/bin/bcftools filter -i '(DP4[0]+DP4[1]) < (DP4[2]+DP4[3]) && ((DP4[2]+DP4[3]) > 0)' calls.vcf -o calls2.vcf

    /usr/local/miniconda/bin/bgzip calls2.vcf

    /usr/local/miniconda/bin/tabix calls2.vcf.gz 

     cat lava_ref.fasta | /usr/local/miniconda/bin/bcftools consensus calls2.vcf.gz > consensus.fasta

     #python3 ${MAFFT_PREP}

     #mafft --quiet aligner.fasta > lava.ali

     python3 ${GFF_WRITE}

	 # Avoiding filename colision during run_pipeline process 
	 mv ${CONTROL_FASTQ} CONTROL.fastq

    """
}


process Alignment_prep { 
    
    container "quay.io/vpeddu/lava_image:latest"

    //errorStrategy 'retry'
    //maxRetries 3
    // Define the input files
    input:
      file "lava_ref.gbk"
      file "lava_ref.fasta"
      file "consensus.fasta"
      file "lava_ref.gff"
	  //file "*.pileup"
    // Define the output files
    output: 
	tuple file('consensus.fasta.amb'), file('consensus.fasta.bwt'), file('consensus.fasta.sa'), file('consensus.fasta'), file('consensus.fasta.ann'), file('consensus.fasta.pac')
	file "AT_refGene.txt"
	file "AT_refGeneMrna.fa"

    // Code to be executed inside the task
    script:
    """
    #!/bin/bash

    /usr/local/miniconda/bin/bwa index consensus.fasta

	/usr/local/miniconda/bin/samtools faidx consensus.fasta 

	gatk CreateSequenceDictionary -R consensus.fasta  --VERBOSITY ERROR --QUIET true


	# Annovar db build step
	gff3ToGenePred lava_ref.gff AT_refGene.txt -warnAndContinue -useName -allowMinimalGenes

	retrieve_seq_from_fasta.pl --format refGene --seqfile consensus.fasta AT_refGene.txt --out AT_refGeneMrna.fa 

    """
}

process Align_samples { 

   container "quay.io/vpeddu/lava_image:latest"

    //errorStrategy 'retry'
    //maxRetries 3
    // Define the input files

    input:
	tuple file(R1), val(PASSAGE)
	tuple file('consensus.fasta.amb'), file('consensus.fasta.bwt'), file('consensus.fasta.sa'), file('consensus.fasta'), file('consensus.fasta.ann'), file('consensus.fasta.pac')
	val(INPUT)
	tuple val(FIRST_FILE), val(PASSAGE)

	output: 
	tuple file(R1), file("*.pileup")
	file "FIRSTFILE.bam" optional true



	shell:
	'''
	#!/bin/bash

	echo aligning "!{R1}"


	/usr/local/miniconda/bin/bwa mem -M -R \'@RG\\tID:group1\\tSM:!{R1}\\tPL:illumina\\tLB:lib1\\tPU:unit1\' -p -t 6 -L [17,17] consensus.fasta !{R1} > !{R1}.sam


	java -jar /usr/bin/picard.jar SortSam INPUT=!{R1}.sam OUTPUT=!{R1}.bam SORT_ORDER=coordinate VERBOSITY=ERROR 

	# TODO NEED TO FIX THIS PARAMETER. ITS PERMANENTLY FALSE FOR NOW 
	if ${params.DEDUPLICATE} 
		then
			echo Deduplicating ${R1}
			java -jar /usr/bin/picard.jar MarkDuplicates INPUT=${R1}.bam OUTPUT=${R1}_dedup.bam METRICS_FILE=metrics.txt VERBOSITY=ERROR REMOVE_DUPLICATES=true
			cat ${R1}_dedup.bam > ${R1}.bam
	fi

	java -jar /usr/bin/picard.jar BuildBamIndex INPUT=!{R1}.bam VERBOSITY=ERROR

	echo sample\tposition\tcov > !{R1}.genomecov
	
	/usr/local/miniconda/bin/bedtools genomecov -d -ibam  !{R1}.bam >> !{R1}.genomecov

	/usr/local/miniconda/bin/samtools mpileup -f consensus.fasta !{R1}.bam > !{R1}.pileup


		if [[ "`basename !{FIRST_FILE}`" == "`basename !{R1}`" ]]
		then 
			echo `basename !{FIRST_FILE}` found
			mv !{R1}.bam FIRSTFILE.bam

		else 
			echo "not first file"
	fi


	'''

}

process Pipeline_prep { 

    //errorStrategy 'retry'
    //maxRetries 3
    // Define the input files

	container "quay.io/vpeddu/lava_image:latest"

	input: 
		file blank_ignore
		file "lava_ref.gff"
		file INITIALIZE_MERGED_CSV
		file CONTROL_FASTQ	
		tuple file('consensus.fasta.amb'), file('consensus.fasta.bwt'), file('consensus.fasta.sa'), file('consensus.fasta'), file('consensus.fasta.ann'), file('consensus.fasta.pac')


	output: 
		file 'merged.csv'
		file 'proteins.csv'
		file "${CONTROL_FASTQ}.pileup"
		file "${CONTROL_FASTQ}.bam"

	script:
	"""
	#!/bin/bash

	echo "Sample,Amino Acid Change,Position,AF,Change,Protein,NucleotideChange,LetterChange,Syn,Depth,Passage" > merged.csv

	python3 ${INITIALIZE_MERGED_CSV} 


	# Creating pileup for control fastq here 

	/usr/local/miniconda/bin/bwa mem -M -R \'@RG\\tID:group1\\tSM:${CONTROL_FASTQ}\\tPL:illumina\\tLB:lib1\\tPU:unit1\' -p -t 6 -L [17,17] consensus.fasta ${CONTROL_FASTQ} > ${CONTROL_FASTQ}.sam

	java -jar /usr/bin/picard.jar SortSam INPUT=${CONTROL_FASTQ}.sam OUTPUT=${CONTROL_FASTQ}.bam SORT_ORDER=coordinate VERBOSITY=ERROR 

	# TODO NEED TO FIX THIS PARAMETER. ITS PERMANENTLY FALSE FOR NOW 
	if ${params.DEDUPLICATE} 
		then
			echo Deduplicating ${CONTROL_FASTQ}
			java -jar /usr/bin/picard.jar MarkDuplicates INPUT=${CONTROL_FASTQ}.bam OUTPUT=${CONTROL_FASTQ}_dedup.bam METRICS_FILE=metrics.txt VERBOSITY=ERROR REMOVE_DUPLICATES=true
			cat ${CONTROL_FASTQ}_dedup.bam > ${CONTROL_FASTQ}.bam
	fi

	java -jar /usr/bin/picard.jar BuildBamIndex INPUT=${CONTROL_FASTQ}.bam VERBOSITY=ERROR

	echo sample\tposition\tcov > ${CONTROL_FASTQ}.genomecov
	
	/usr/local/miniconda/bin/bedtools genomecov -d -ibam  ${CONTROL_FASTQ}.bam >> ${CONTROL_FASTQ}.genomecov

	/usr/local/miniconda/bin/samtools mpileup -f consensus.fasta ${CONTROL_FASTQ}.bam > ${CONTROL_FASTQ}.pileup

	

	"""
}

process Create_VCF { 

    //errorStrategy 'retry'
    //maxRetries 3
    // Define the input files

	container "quay.io/vpeddu/lava_image:latest"

	input: 
		file CONTROL_FASTQ
		file CONTROL_PILEUP
		tuple file(R1), file(R1_PILEUP)
		file ATREF
		tuple val(FIRST_FILE), val(PASSAGE)
		file ATREF_MRNA


	output: 
		file "*exonic_variant_function" optional true
	shell:
	'''
	#!/bin/bash
## 
	echo Analyzing variants in sample !{R1}

	java -jar /usr/local/bin/VarScan somatic !{CONTROL_PILEUP} !{R1_PILEUP} !{R1}.vcf --validation 1 --output-vcf 1 --min-coverage 2

	mv !{R1}.vcf.validation !{R1}.vcf

	awk -F $\'\t\' \'BEGIN {FS=OFS="\t"}{gsub("0/0","0/1",$10)gsub("0/0","1/0",$11)gsub("1/1","0/1",$10)gsub("1/1","1/0",$11)}1\' !{R1}.vcf > !{R1}_p.vcf

	file="!{R1}""_p.vcf"


	sed 's/NC_039477.1/lava/g' !{R1}_p.vcf > test.txt 

	mv test.txt !{R1}_p.vcf

	#convert2annovar.pl -withfreq -format vcf4 -includeinfo !{R1}_p.vcf > !{R1}.avinput 
	convert2annovar.pl -withfreq -format vcf4old -includeinfo !{R1}_p.vcf > !{R1}.avinput 

	annotate_variation.pl -outfile !{R1} -v -buildver AT !{R1}.avinput .


	ls -lah 
	
	if [[ "`basename !{FIRST_FILE}`" == "`basename !{R1}`" ]]
		then 
			echo `basename !{FIRST_FILE}` found

		else 
			echo "not first file"
			echo `basename !{R1}` `basename !{FIRST_FILE}`
			rm *exonic_variant_function
	fi
			
	'''
}

process Ref_done { 

    //errorStrategy 'retry'
    //maxRetries 3
    // Define the input files

	container "quay.io/vpeddu/lava_image:latest"

	input: 
		tuple file(FIRST_FILE), val(PASSAGE)
		val(AF)
		file exonic_variant_function
		file CONTROL_FASTQ
		file CONTROL_BAM
		file FIRSTBAM


	output: 
	shell:
	'''
	#!/bin/bash

	echo !{FIRST_FILE}

	# Michelle what the hell is with all of this awk 

	#TODO NO OPTION FOR ALLELE FILTERING RIGHT NOW
	awk -F":" \'($18+0)>=1{print}\' !{FIRST_FILE}.exonic_variant_function > ref.txt

	grep "SNV" ref.txt > a.tmp && mv a.tmp ref.txt 

	awk -v ref=!{CONTROL_FASTQ} -F '[\t:,]' \
	'{print ref,","$6" "substr($9,3)","$12","$39+0","substr($9,3)","$6","substr($8,3)","substr($8,3,1)" \
	to "substr($8,length($8))","$2","$36",0"}' ref.txt > ref.csv

	printf !{FIRST_FILE}"," >> reads.csv

	/usr/local/miniconda/bin/samtools flagstat !{FIRSTBAM} | \
	awk 'NR==1{printf $1","} NR==5{printf $1","} NR==5{print substr($5,2)}' >> reads.csv

	 echo sample	position	cov > !{FIRST_FILE}.genomecov
	 /usr/local/miniconda/bin/bedtools genomecov -d -ibam !{FIRSTBAM} >> !{FIRST_FILE}.genomecov


	 # TODO AF FLAG HERE 

	 awk -F":" '($24+0)>=1{print}' !{FIRST_FILE}.exonic_variant_function > !{FIRST_FILE}.txt 

	 grep "SNV" !{FIRST_FILE}.txt > a.tmp
	 grep "stop" !{FIRST_FILE}.txt >> a.tmp
	 mv a.tmp !{FIRST_FILE}.txt

	 awk -v name=!{FIRST_FILE} -v \
	 sample=$SAMPLE -F'[\t:,]' '{print name","$6" "substr($9,3)","$12","$49+0","substr($9,3)","$6","substr($8,3)","substr($8,3,1)" to "substr($8,length($8))","$2","$46","sample}'\
	!{FIRST_FILE}.txt > !{FIRST_FILE}.csv

	'''
}


