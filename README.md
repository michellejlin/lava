# lava: longitudinal analysis and variant annotation
lava is a script that visualizes the minor allele variants in longitudinal sequences. Given a fastq file to serve as a control (usually in longitudinal studies this will be the sequence at Passage 0), lava will annotate the amino acid changes in the control sequence and all subsequent sequences. These changes will then be plotted for ease of viewability of the changes in minor allele variants.

# Input Files

**Mandatory Files**
1. A fastq must be specified as the control. This will be the sample to which lava aligns the rest of the sequences to. Also in the folder include all the sequences to be aligned to the control, but no other.
2. A metadata.csv file must be included in the folder with the fastq. This will contain two columns: Sample and Passage. 

**Optional Files**
1. A fasta of the control may be provided. Otherwise, lava will automatically align the control .fastq to a reference sequence of the virus pulled from Genbank.
2. A gff may also be provided of the annotations of the coding regions of the control. Otherwise, lava will automatically transfer the annotations off a reference sequence pulled from Genbank.

**Note:**
The optional files should not have any special characters in them, particularly underscores. Simply renaming the fasta or gff will not work - make sure the names within the files themselves do not have special characters as well.

# Usage

The general usage of lava is 

`lava.sh [options] control.fastq`

For more detailed instructions: Make sure lava is in your path. This can be done by the command:

`echo 'export PATH=${PATH}:/path/to/lava' >> ~/.bash_profile`

where path/to/lava is the absolute path to the directory in which lava is downloaded.

Create input files as desired (for most, this may only be the metadata.csv file). 

Move to the folder your sequences are in. Make sure this folder contains the necessary input files: the control fastq, the rest of your sequences in fastq format, the metadata.csv file containing information about your samples, and optionally the fasta and gff file.

Then, run the command: 

`lava.sh -q "[user query for virus]" fastq_file`

with your control fastq specified, as well as what virus the sequence is, or the accession number desired for mapping, for the default lava usage. Running the command:

`lava.sh -f [fasta_file] -g [gff_file] fastq_file`

allows the user to specify their own fasta file or gff file for lava. Please note the -q and -g arguments CANNOT be both specified, as only one option for gff generation should be used. The -g argument is intended only for highly specific cases; otherwise, the -q argument in which lava automatically generates the gff is highly recommended.

# Output Files

lava will output an html file of the visualization that can be easily shareable and downloaded. In addition, lava will also output a summary table of the amino acid changes taking place, called merged.csv. 

# Common Errors

1. `WARNING: A total of 1 sequences will be ignored due to lack of correct ORF annotation`
	
	This error will occur when the open reading frame is judged to be wrong by Annovar. It usually happens when the frame does not end in a stop codon.