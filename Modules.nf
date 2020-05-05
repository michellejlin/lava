/*
 * Define the parameters used in the processes below
 */

/*
 * Define the processes used in this workflow
 */

process setup{ 

    //Retry at most 3 times
    errorStrategy 'retry'
    maxRetries 3
    
    // Define the Docker container used for this step
    // TODO FILL THIS IN 
    container "docker pull quay.io/vpeddu/lava_image:latest"

    // Define the input files
    input:
      file R1
    // Define the output files

    script: 
    """
    #!/usr/bin/env bash

    ls 


    """
}


process CreateConsensus { 
    //Retry at most 3 times
    //errorStrategy 'retry'
    //maxRetries 3
    
    // Define the Docker container used for this step
    // TODO FILL THIS IN 
    container "quay.io/vpeddu/lava_image:latest"

    // Define the input files
    input:
      file METADATA_FILE
      file R1
      val(GENBANK_ACC)
      file GFF_FILE
      file FASTA_FILE
      val(DEDUPLICATE)
    // Define the output files
    output:
      file 
    // Code to be executed inside the task
    script:
    """
    #!/usr/bin/env bash

    echo "Indexing reference"

    /usr/local/miniconda/bin/bwa index ${FASTA_FILE}

    echo "Finished indexing reference"

    /usr/local/miniconda/bin/samtools faidx ${FASTA_FILE}

    GATK CreateSequenceDictionary -R ${FASTA_FILE} --VERBOSITY ERROR --QUIET true
    """
}

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
	tuple file('lava_ref.fasta.amb'), file('lava_ref.fasta.bwt'), file('lava_ref.fasta.sa'), file('lava_ref.fasta'), file('lava_ref.fasta.ann'), file('lava_ref.fasta.pac')


    // Code to be executed inside the task
    script:
    """
    #!/bin/bash

    /usr/local/miniconda/bin/bwa index lava_ref.fasta

	/usr/local/miniconda/bin/samtools faidx lava_ref.fasta 

	gatk CreateSequenceDictionary -R lava_ref.fasta  --VERBOSITY ERROR --QUIET true


	# Annovar db build step
	gff3ToGenePred lava_ref.gff AT_refGene.txt -warnAndContinue -useName -allowMinimalGenes

	retrieve_seq_from_fasta.pl --format refGene --seqfile lava_ref.fasta AT_refGene.txt --out AT_refGeneMrna.fa 

    """
}

process Align_samples { 

   container "quay.io/vpeddu/lava_image:latest"

    //errorStrategy 'retry'
    //maxRetries 3
    // Define the input files

    input:
	tuple file(R1), val(PASSAGE)
	tuple file('lava_ref.fasta.amb'), file('lava_ref.fasta.bwt'), file('lava_ref.fasta.sa'), file('lava_ref.fasta'), file('lava_ref.fasta.ann'), file('lava_ref.fasta.pac')
	val(INPUT)

	output: 

	script:
	"""
	#!/bin/bash

	echo aligning "${R1}"

	file_loc="${R1}"

	/usr/local/miniconda/bin/bwa mem -M -R \'@RG\\tID:group1\\tSM:${R1}\\tPL:illumina\\tLB:lib1\\tPU:unit1\' -p -t 6 -L [17,17] lava_ref.fasta ${R1} > ${R1}.sam

	java -jar /usr/bin/picard.jar SortSam INPUT=${R1}.sam OUTPUT=${R1}.bam SORT_ORDER=coordinate VERBOSITY=ERROR 

	# TODO NEED TO FIX THIS PARAMETER. ITS FALSE FOR NOW PERMANENTLY
	if ${params.DEDUPLICATE} 
		then
			echo Deduplicating ${R1}
			java -jar /usr/bin/picard.jar MarkDuplicates INPUT=${R1}.bam OUTPUT=${R1}_dedup.bam METRICS_FILE=metrics.txt VERBOSITY=ERROR REMOVE_DUPLICATES=true
			cat ${R1}_dedup.bam > ${R1}.bam
	fi

	java -jar /usr/bin/picard.jar BuildBamIndex INPUT=${R1}.bam VERBOSITY=ERROR

	echo sample\tposition\tcov > ${R1}.genomecov
	
	/usr/local/miniconda/bin/bedtools genomecov -d -ibam  ${R1}.bam >> ${R1}.genomecov

	/usr/local/miniconda/bin/samtools mpileup -f lava_ref.fasta ${R1}.bam > ${R1}.pileup


	"""

}

process Create_Merged_CSV { 

    //errorStrategy 'retry'
    //maxRetries 3
    // Define the input files

	container "quay.io/vpeddu/lava_image:latest"

	input: 
		file blank_ignore
		file "lava_ref.gff"
		file INITIALIZE_MERGED_CSV

	output: 
		file 'merged.csv'
		file 'proteins.csv'
	script:
	"""
	#!/bin/bash

	echo "Sample,Amino Acid Change,Position,AF,Change,Protein,NucleotideChange,LetterChange,Syn,Depth,Passage" > merged.csv

	python3 ${INITIALIZE_MERGED_CSV} 
	"""
}

process Run_pipeline { 

    //errorStrategy 'retry'
    //maxRetries 3
    // Define the input files

	container "quay.io/vpeddu/lava_image:latest"

	input: 
		file R1

	output: 
	shell:
	'''
	#!/bin/bash

	echo Analyzing variants in sample !{R1}

	# java -jar /usr/local/bin/VarScan somatic 

	'''
}