// Uses accession number specified by --GENBANK to create our own GFF (lava_ref.gff) for a consensus fasta
// generated from the alignment of "Passage 0" sample to reference fasta.
process CreateGFF { 
    container "quay.io/vpeddu/lava_image:latest"

	// Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 1
	
    input:
      val(GENBANK)
      file CONTROL_FASTQ
	  file PULL_ENTREZ
	  file WRITE_GFF
	  //file FASTA
	  //file GFF

    output: 
      file "lava_ref.fasta"
      file "consensus.fasta"
      file "lava_ref.gff"
	  file "CONTROL.fastq"
	  file "ribosomal_start.txt"
	  file "mat_peptides.txt"

    script:
    """
    #!/bin/bash
    
	set -e 
	echo ${CONTROL_FASTQ}
    
	# Pulls reference fasta and GenBank file using accession number specified by --GENBANK.
	python3 ${PULL_ENTREZ} ${GENBANK}

	# Indexes and aligns "Sample 0" fastq to reference fasta
    /usr/local/miniconda/bin/bwa index lava_ref.fasta
    /usr/local/miniconda/bin/bwa mem -t ${task.cpus} -M lava_ref.fasta ${CONTROL_FASTQ} | /usr/local/miniconda/bin/samtools view -Sb - > aln.bam

	# Generates new consensus fasta from aligned "Sample 0" and reference fasta.
	/usr/local/miniconda/bin/samtools sort -@ ${task.cpus} aln.bam -o aln.sorted.bam 
    /usr/local/miniconda/bin/bcftools mpileup --max-depth 500000 -P 1.1e-100 -Ou -f lava_ref.fasta aln.sorted.bam | /usr/local/miniconda/bin/bcftools call -m -Oz -o calls.vcf.gz 
    /usr/local/miniconda/bin/tabix calls.vcf.gz
    gunzip calls.vcf.gz 
    /usr/local/miniconda/bin/bcftools filter -i '(DP4[0]+DP4[1]) < (DP4[2]+DP4[3]) && ((DP4[2]+DP4[3]) > 0)' calls.vcf -o calls2.vcf
    /usr/local/miniconda/bin/bgzip calls2.vcf
    /usr/local/miniconda/bin/tabix calls2.vcf.gz 
    cat lava_ref.fasta | /usr/local/miniconda/bin/bcftools consensus calls2.vcf.gz > consensus.fasta

	# Creates a GFF (lava_ref.gff) for our consensus fasta per CDS annotations from our reference GenBank file.
	python3 ${WRITE_GFF}

	 # Avoiding filename collision during run_pipeline process 
	 mv ${CONTROL_FASTQ} CONTROL.fastq
    """
}

// Prep for downstream steps.
// Indexes and prepares consensus fasta, and generates prerequisite files for Annovar.
process Alignment_prep { 
    container "quay.io/vpeddu/lava_image:latest"

    errorStrategy 'retry'
    maxRetries 3

    input:
      file "lava_ref.fasta"
      file "consensus.fasta"
      file "lava_ref.gff"

    output: 
	tuple file('consensus.fasta.amb'), file('consensus.fasta.bwt'), file('consensus.fasta.sa'), file('consensus.fasta'), file('consensus.fasta.ann'), file('consensus.fasta.pac')
	file "AT_refGene.txt"
	file "AT_refGeneMrna.fa"

    script:
    """
    #!/bin/bash

	# Indexes and prepares consensus fasta for downstream steps
    /usr/local/miniconda/bin/bwa index consensus.fasta
	/usr/local/miniconda/bin/samtools faidx consensus.fasta 
	gatk CreateSequenceDictionary -R consensus.fasta  --VERBOSITY ERROR --QUIET true

	# Preparatory steps for Annovar downstream
	# Creates Annovar database from our consensus fasta and generated GFF
	gff3ToGenePred lava_ref.gff AT_refGene.txt -warnAndContinue -useName -allowMinimalGenes
	retrieve_seq_from_fasta.pl --format refGene --seqfile consensus.fasta AT_refGene.txt --out AT_refGeneMrna.fa 

    """
}

// Aligns all samples to consensus fasta and removes duplicates if --DEDUPLICATE specified.
// Also generates genomecov files and pileups.
process Align_samples { 
   container "quay.io/vpeddu/lava_image:latest"

    errorStrategy 'retry'
    maxRetries 3

    input:
	tuple file(R1), val(PASSAGE)
	tuple file('consensus.fasta.amb'), file('consensus.fasta.bwt'), file('consensus.fasta.sa'), file('consensus.fasta'), file('consensus.fasta.ann'), file('consensus.fasta.pac')
	tuple val(FIRST_FILE), val(NULL)
	val DEDUPLICATE

	output: 
	tuple file(R1), file("*.pileup"), file("*.bam"), val(PASSAGE)
	file "FIRSTFILE.bam" optional true
	file "${R1}.genomecov"

	shell:
	'''
	#!/bin/bash

	echo Aligning" !{R1}"

	# Align each sample to consensus fasta.
	/usr/local/miniconda/bin/bwa mem -t !{task.cpus} -M -R \'@RG\\tID:group1\\tSM:!{R1}\\tPL:illumina\\tLB:lib1\\tPU:unit1\' -p -L [17,17] consensus.fasta !{R1} > !{R1}.sam

	# Sorts SAM.
	java -jar /usr/bin/picard.jar SortSam INPUT=!{R1}.sam OUTPUT=!{R1}.bam SORT_ORDER=coordinate VERBOSITY=ERROR 

	# Removes duplicates (e.g. from library construction using PCR) if --DEDUPLICATE flag specified.
	if !{DEDUPLICATE} 
		then
			echo "Deduplicating !{R1}"
			java -jar /usr/bin/picard.jar MarkDuplicates INPUT=${R1}.bam OUTPUT=${R1}_dedup.bam METRICS_FILE=metrics.txt VERBOSITY=ERROR REMOVE_DUPLICATES=true
			cat ${R1}_dedup.bam > ${R1}.bam
	fi

	# Creates genomecov file from BAM so we can generate coverage graphs later.
	echo sample\tposition\tcov > !{R1}.genomecov	
	/usr/local/miniconda/bin/bedtools genomecov -d -ibam  !{R1}.bam >> !{R1}.genomecov

	java -jar /usr/bin/picard.jar BuildBamIndex INPUT=!{R1}.bam VERBOSITY=ERROR
	# Generates pileup that VCF can be called off of later.
	/usr/local/miniconda/bin/samtools mpileup --max-depth 500000 -f consensus.fasta !{R1}.bam > !{R1}.pileup

	# Renames first file to avoid file name collision.
		if [[ "`basename !{FIRST_FILE}`" == "`basename !{R1}`" ]]
		then 
			echo `basename !{FIRST_FILE}` found
			mv !{R1}.bam FIRSTFILE.bam

		else 
			echo "not first file"
	fi

	'''
}

// Initializes proteins.csv - list of protein names and locations - from our generated GFF.
// Also does same thing as Align_samples for our "Sample 0" file.
process Pipeline_prep { 

    errorStrategy 'retry'
    maxRetries 3

	container "quay.io/vpeddu/lava_image:latest"

	input: 
		file blank_ignore
		file "lava_ref.gff"
		file CONTROL_FASTQ	
		tuple file('consensus.fasta.amb'), file('consensus.fasta.bwt'), file('consensus.fasta.sa'), file('consensus.fasta'), file('consensus.fasta.ann'), file('consensus.fasta.pac')
		file INITIALIZE_PROTEINS_CSV

	output: 
		file 'merged.csv'
		file 'proteins.csv'
		file "${CONTROL_FASTQ}.pileup"
		file "${CONTROL_FASTQ}.bam"

	script:
	"""
	#!/bin/bash

	# Creates header for final csv.
	echo "Sample,Amino Acid Change,Position,AF,Change,Protein,NucleotideChange,LetterChange,Syn,Depth,Passage" > merged.csv

	# Creates list of protein names and locations (proteins.csv) based on GFF annotations.
	python3 ${INITIALIZE_PROTEINS_CSV}

	# Aligns and generates genomecov and pileup for "Sample 0" fastq.
	/usr/local/miniconda/bin/bwa mem -t ${task.cpus}  -M -R \'@RG\\tID:group1\\tSM:${CONTROL_FASTQ}\\tPL:illumina\\tLB:lib1\\tPU:unit1\' -p -L [17,17] consensus.fasta ${CONTROL_FASTQ} > ${CONTROL_FASTQ}.sam
	java -jar /usr/bin/picard.jar SortSam INPUT=${CONTROL_FASTQ}.sam OUTPUT=${CONTROL_FASTQ}.bam SORT_ORDER=coordinate VERBOSITY=ERROR 

	# Removes duplicates if specified.
	if ${params.DEDUPLICATE} 
		then
			echo Deduplicating ${CONTROL_FASTQ}
			java -jar /usr/bin/picard.jar MarkDuplicates INPUT=${CONTROL_FASTQ}.bam OUTPUT=${CONTROL_FASTQ}_dedup.bam METRICS_FILE=metrics.txt VERBOSITY=ERROR REMOVE_DUPLICATES=true
			cat ${CONTROL_FASTQ}_dedup.bam > ${CONTROL_FASTQ}.bam
	fi

	java -jar /usr/bin/picard.jar BuildBamIndex INPUT=${CONTROL_FASTQ}.bam VERBOSITY=ERROR

	# Creates genome coverage file for "Sample 0" fastq.
	echo sample\tposition\tcov > ${CONTROL_FASTQ}.genomecov
	/usr/local/miniconda/bin/bedtools genomecov -d -ibam  ${CONTROL_FASTQ}.bam >> ${CONTROL_FASTQ}.genomecov

	# Creates pileup for "Sample 0" fastq.
	/usr/local/miniconda/bin/samtools mpileup --max-depth 500000 -f consensus.fasta ${CONTROL_FASTQ}.bam > ${CONTROL_FASTQ}.pileup
	"""
}

// Generates VCF for all the samples and converts to .avinput for Annovar.
process Create_VCF { 
    errorStrategy 'retry'
    maxRetries 3

	container "quay.io/vpeddu/lava_image:latest"

	input: 
		file CONTROL_FASTQ
		file CONTROL_PILEUP
		tuple file(R1), file(R1_PILEUP), file(BAM), val(PASSAGE)
		file ATREF
		tuple val(FIRST_FILE), val(NULL)
		file ATREF_MRNA


	output: 
		file "*exonic_variant_function" optional true
		tuple file(R1), file("*.bam"), file( "*.exonic_variant_function.samp"), val(PASSAGE)
		file "${R1}.vcf"

	shell:
	'''
	#!/bin/bash

	ls -latr
	echo Analyzing variants in sample !{R1}

	# here for file passthrough (input -> output)
	mv !{BAM} !{BAM}.bam 

	# Generates VCF outputting all bases with a min coverage of 2.
	java -jar /usr/local/bin/VarScan somatic !{CONTROL_PILEUP} !{R1_PILEUP} !{R1}.vcf --validation 1 --output-vcf 1 --min-coverage 2
	mv !{R1}.vcf.validation !{R1}.vcf

	# Fixes ploidy issues.
	awk -F $\'\t\' \'BEGIN {FS=OFS="\t"}{gsub("0/0","0/1",$10)gsub("0/0","1/0",$11)gsub("1/1","0/1",$10)gsub("1/1","1/0",$11)}1\' !{R1}.vcf > !{R1}_p.vcf

	# Converts VCF to .avinput for Annovar.
	file="!{R1}""_p.vcf"
	#convert2annovar.pl -withfreq -format vcf4 -includeinfo !{R1}_p.vcf > !{R1}.avinput 
	convert2annovar.pl -withfreq -format vcf4old -includeinfo !{R1}_p.vcf > !{R1}.avinput 
	annotate_variation.pl -outfile !{R1} -v -buildver AT !{R1}.avinput .

	ls -lah 
	
	if [[ "`basename !{FIRST_FILE}`" == "`basename !{R1}`" ]]
		then 
			echo `basename !{FIRST_FILE}` found
			touch blank.exonic_variant_function.samp

		else 
			echo "not first file"
			echo `basename !{R1}` `basename !{FIRST_FILE}`
			mv !{R1}.exonic_variant_function !{R1}.exonic_variant_function.samp
	fi
			
	'''
}

// Extract variants for "Passage 0" sample.
process Ref_done { 
    errorStrategy 'retry'
    maxRetries 3

	container "quay.io/vpeddu/lava_image:latest"

	input: 
		tuple file(FIRST_FILE), val(PASSAGE)
		val(ALLELE_FREQ)
		file exonic_variant_function
		file CONTROL_FASTQ
		file CONTROL_BAM
		file FIRSTBAM
		file METADATA


	output:
		tuple file("reads.csv"), val(PASSAGE), file(FIRST_FILE), file("${FIRST_FILE}.csv")
		

	shell:
	'''
	#!/bin/bash

	echo !{FIRST_FILE}
	echo !{ALLELE_FREQ}
	
	# Filters by specified allele frequency; otherwise, outputs all variants with greater than 1% AF.
	if [[ "!{ALLELE_FREQ}" == "NO_VAL" ]]
		then
			awk -F":" '($18+0)>=1{print}' !{FIRST_FILE}.exonic_variant_function > ref.txt

		else	
			awk -v af=!{ALLELE_FREQ} -F":" '($18+0)>=!{ALLELE_FREQ}{print}' !{FIRST_FILE}.exonic_variant_function > ref.txt
	fi

	# Filters by only single nucleotide mutations.
	grep "SNV" ref.txt > a.tmp && mv a.tmp ref.txt 

	# Grabs columns matching our final header.
	awk -v ref=!{CONTROL_FASTQ} -F '[\t:,]' \
	'{print ref,","$6" "substr($9,3)","$12","$39+0","substr($9,3)","$6","substr($8,3)","substr($8,3,1)" \
	to "substr($8,length($8))","$2","$36",0"}' ref.txt > ref.csv

	printf !{FIRST_FILE}"," >> reads.csv

	/usr/local/miniconda/bin/samtools flagstat !{FIRSTBAM} | \
	awk 'NR==1{printf $1","} NR==5{printf $1","} NR==5{print substr($5,2)}' >> reads.csv

	 echo sample	position	cov > !{FIRST_FILE}.genomecov
	 /usr/local/miniconda/bin/bedtools genomecov -d -ibam !{FIRSTBAM} >> !{FIRST_FILE}.genomecov

	if [[ "!{ALLELE_FREQ}" == "NO_VAL" ]]
		then
			awk -F":" '($24+0)>=1{print}' !{FIRST_FILE}.exonic_variant_function > !{FIRST_FILE}.txt 
		else
			awk -v af=!{ALLELE_FREQ} -F":" '($24+0)>=!{ALLELE_FREQ}{print}' !{FIRST_FILE}.exonic_variant_function > !{FIRST_FILE}.txt 
	fi

	 grep "SNV" !{FIRST_FILE}.txt > a.tmp
	 grep "stop" !{FIRST_FILE}.txt >> a.tmp
	 mv a.tmp !{FIRST_FILE}.txt

	SAMPLE="$(awk -F"," -v name=!{FIRST_FILE} '$1==name {print $2}' !{METADATA})"
	
	 awk -v name=!{FIRST_FILE} -v sample=!{PASSAGE} -F'[\t:,]' '{print name","$6" "substr($9,3)","$12","$46+0","substr($9,3)","$6","substr($8,3)","substr($8,3,1)" to "substr($8,length($8))","$2","$43","sample}' !{FIRST_FILE}.txt > !{FIRST_FILE}.csv

	'''
}

// Extract variants for all other samples.
process Extract_variants { 
    errorStrategy 'retry'
    maxRetries 3

	container "quay.io/vpeddu/lava_image:latest"

	input: 
		tuple val(FIRST_FILE), val(NULL)
		tuple file(R1), file(BAM), file(EXONICVARIANTS), val(PASSAGE)
		file METADATA
	output:
		tuple file("${R1}.csv"), val(PASSAGE), file("reads.csv"), file(R1) optional true
		tuple file(R1), val(PASSAGE) optional true
	shell:

	'''
	#!/bin/bash
	echo !{R1}

		if [[ "`basename !{FIRST_FILE}`" == "`basename !{R1}`" ]]
		then 
			echo `basename !{FIRST_FILE}` found
			echo first file found ending process execution for !{R1}
			exit 0

		else 
			echo "not first file"
	fi

	echo "continuing execution for !{R1}"

	# Creates genomecov files for genome coverage graphs later.
	echo 'sample	position	cov' > !{R1}.genomecov 
	/usr/local/miniconda/bin/bedtools genomecov -d -ibam !{BAM} >> !{R1}.genomecov

	# reads.csv from all processes will be merged together at end 
	 printf !{R1}"," > reads.csv

	/usr/local/miniconda/bin/samtools flagstat !{BAM} | \
	awk 'NR==1{printf $1","} NR==5{printf $1","} NR==5{print substr($5,2)}' >> reads.csv

	awk -F":" '($24+0)>=1{print}' !{EXONICVARIANTS}> !{R1}.txt

	 grep "SNV" !{R1}.txt > a.tmp
	 grep "stop" !{R1}.txt >> a.tmp
	 mv a.tmp !{R1}.txt

	SAMPLE="$(awk -F"," -v name=!{R1} '$1==name {print $2}' !{METADATA})"

	echo $SAMPLE
	
	 awk -v name=!{R1} -v sample=!{PASSAGE} -F'[\t:,]' '{print name","$6" "substr($9,3)","$12","$46+0","substr($9,3)","$6","substr($8,3)","substr($8,3,1)" to "substr($8,length($8))","$2","$43","sample}' !{R1}.txt > !{R1}.csv
	'''
}

// Checks for multi-nucleotide mutations and prints out warning message.
// Currently LAVA does not handle complex mutations and instead annotates it as such for manual review.
process Annotate_complex { 
    errorStrategy 'retry'
    maxRetries 3

	container "quay.io/vpeddu/lava_image:latest"

	input: 
		tuple file(SAMPLE_CSV), val(PASSAGE), file("reads.csv"), file(R1)
		file ANNOTATE_COMPLEX_MUTATIONS
	output:
		file R1
		file "${R1}.complex.log"
		file "${R1}.reads.csv"
		file SAMPLE_CSV

	script:
	"""
	#!/bin/bash

	# Checks for complex mutations and prints a warning message.
	python3 ${ANNOTATE_COMPLEX_MUTATIONS} ${SAMPLE_CSV} ${PASSAGE}	

	# Renaming files to avoid file collision
	mv complex.log ${R1}.complex.log
	mv reads.csv ${R1}.reads.csv
	"""
}

// Checks for multi-nucleotide mutations in first file and prints out warning message.
process Annotate_complex_first_passage { 
    errorStrategy 'retry'
    maxRetries 3

	container "quay.io/vpeddu/lava_image:latest"

	input: 
		tuple file("reads.csv"), val(PASSAGE), file(FIRST_FILE), file(FIRST_FILE_CSV)
		file ANNOTATE_COMPLEX_MUTATIONS

	output:
		tuple file(FIRST_FILE), val(PASSAGE), file("*.complex.log"), file("*.reads.csv"), file(FIRST_FILE_CSV)
	script:

	"""
	#!/bin/bash

	# Checks for complex mutations and prints a warning message.
	python3 ${ANNOTATE_COMPLEX_MUTATIONS} ${FIRST_FILE_CSV} ${PASSAGE}	
	mv complex.log ${FIRST_FILE}.complex.log
	mv reads.csv ${FIRST_FILE}.reads.csv
	"""
}

process Generate_output { 
    errorStrategy 'retry'
    maxRetries 3

	container "quay.io/vpeddu/lava_image:latest"

	input: 
		tuple file(FIRST_R1), val(FIRST_PASSAGE), file(FIRST_COMPLEX_LOG), file(FIRST_READS_CSV), file(FIRST_SAMPLE_CSV)
		file R1 
		file COMPLEX_LOG
		file READS_CSV
		file SAMPLE_CSV
		file MERGED_CSV
		file PROTEINS_CSV
		file GENOMECOV
		file VCF
		file RIBOSOMAL_LOCATION
		file MAT_PEPTIDE_LOCATIONS
		file MAT_PEPTIDE_ADDITION
		file RIBOSOMAL_SLIPPAGE
		file GENOME_PROTEIN_PLOTS
		file PALETTE

	output:
		file "*.html"
		file "*.log"
		file "final.csv"
		file "*.csv"
		file "vcf_files"
		file "genomecov"
		file "all_files"
	script:

	// Assumes that "Passage" info given in metadata file is a numerical value.
	if (params.CATEGORICAL == 'false') {
		"""
		#!/bin/bash

		ls -lah

		# cat *fastq.csv >> merged.csv
		cat merged.csv > final.csv 

		# Takes fastq.gz and fastq
		# if [[ gzip -t \$${R1} ]]
		if ls *.gz &>/dev/null
		then
			cat *.fastq.gz.csv >> final.csv
		else
			cat *.fastq.csv >> final.csv
		fi

		# Gets rid of non-SNPs
		grep -v "transcript" final.csv > a.tmp && mv a.tmp final.csv 
		grep -v "delins" final.csv > a.tmp && mv a.tmp final.csv 

		# Sorts by beginning of mat peptide
		sort -k2 -t, -n mat_peptides.txt > a.tmp && mv a.tmp mat_peptides.txt
		# Adds mature peptide differences from protein start.
		python3 ${MAT_PEPTIDE_ADDITION}
		rm mat_peptides.txt

		# Corrects for ribosomal slippage.
		python3 ${RIBOSOMAL_SLIPPAGE} final.csv proteins.csv

		awk NF final.csv > a.tmp && mv a.tmp final.csv

		cat *.reads.csv > reads.csv 

		cat *.log > complex.log
		# TODO error handling @ line 669-683 of lava.py 
		python3 ${GENOME_PROTEIN_PLOTS} visualization.csv proteins.csv reads.csv . "Plot"

		mkdir vcf_files
		mv *.vcf vcf_files

		mkdir genomecov
		mv *.genomecov genomecov

		mkdir all_files

		cp -r *.txt all_files
		"""
	} 
	// Interprets "Passage" info given in metadata file as categorical, and will rename samples in visualization 
	// based on that categorical value.
	else {
		"""
		#!/bin/bash

		ls -lah

		# cat *fastq.csv >> merged.csv

		head ${PALETTE}

		cat merged.csv > final.csv 

		#Takes fastq.gz and fastq
		# if [[ gzip -t \$${R1} ]]
		if ls *.gz &>/dev/null
		then
			cat *.fastq.gz.csv >> final.csv
		else
			cat *.fastq.csv >> final.csv
		fi

		grep -v "transcript" final.csv > a.tmp && mv a.tmp final.csv 

		grep -v "delins" final.csv > a.tmp && mv a.tmp final.csv 

		# Sorts by beginning of mat peptide
		sort -k2 -t, -n mat_peptides.txt > a.tmp && mv a.tmp mat_peptides.txt
		# Adds mature peptide differences from protein start.
		python3 ${MAT_PEPTIDE_ADDITION}
		rm mat_peptides.txt

		# Corrects for ribosomal slippage.
		python3 ${RIBOSOMAL_SLIPPAGE} final.csv proteins.csv


		awk NF final.csv > a.tmp && mv a.tmp final.csv

		cat *.reads.csv > reads.csv 

		cat *.log > complex.log
		# TODO error handling @ line 669-683 of lava.py 
		python3 ${GENOME_PROTEIN_PLOTS} visualization.csv proteins.csv reads.csv . "Plot" -categorical

		mkdir vcf_files
		mv *.vcf vcf_files

		mkdir genomecov
		mv *.genomecov genomecov

		mkdir all_files

		cp -r *.txt all_files
		"""
	}
} 
