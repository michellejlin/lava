#!/bin/bash

   /usr/local/miniconda/bin/bwa index consensus.fasta

/usr/local/miniconda/bin/samtools faidx consensus.fasta 

gatk CreateSequenceDictionary -R consensus.fasta  --VERBOSITY ERROR --QUIET true


# Annovar db build step
gff3ToGenePred lava_ref.gff AT_refGene.txt -warnAndContinue -useName -allowMinimalGenes

retrieve_seq_from_fasta.pl --format refGene --seqfile consensus.fasta AT_refGene.txt --out AT_refGeneMrna.fa
