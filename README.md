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

Create input files as desired (for most, this may only be the metadata.csv file). 

Move to the folder lava is in. For example, if lava were cloned from Github into downloads:

`cd C:\User\Downloads\lava\`

or the desired file path. 

Then, run the command: 

`lava.sh [fastq_file] -g "[user query for virus]"`

with your control fastq specified, as well as what virus the sequence is, or the accession number desired for mapping, for the default lava usage. Running the command:

`lava.sh [fastq_file] -f [fasta_file] -r [gff_file]`

allows the user to specify their own fasta file or gff file for lava. Please note the -g and -r arguments CANNOT be both specified, as only one option for gff generation should be used. The -r argument is intended only for highly specific cases; otherwise, the -g argument in which lava automatically generates the gff is highly recommended.

# Output Files

lava will output an html file of the visualization that can be easily shareable and downloaded. In addition, lava will also output a summary table of the amino acid changes taking place, called merged.csv. 