# lava: Longitudinal Analysis and Variant Annotation
lava analyzes and visualizes minor allele variants in longitudinal sequence data. lava takes a reference fasta (normally representing the first sample in your longitudinal analysis), fastq files (for every sample in your analysis), and a metadata sheet (providing information on what day or passage each sample was collected). Output will be displayed as an interactive graph in your web browser. lava will only work on Mac and Linux machines. 

# Installation
**Dependencies**
The following dependencies are required to run lava: 
1. Python
2. BioPython, numpy, pandas, seaborn, bokeh modules for python which can be installed with `python -m pip install biopython` and `python -m pip install whateverpackage` 
3. Picard (https://broadinstitute.github.io/picard/)
4. GATK (https://software.broadinstitute.org/gatk/download/)
5. VARSCAN (http://varscan.sourceforge.net/)
6. ANNOVAR (http://annovar.openbioinformatics.org/en/latest/user-guide/download/)
7. snpEff (http://snpeff.sourceforge.net/download.html)
8. maaft (https://mafft.cbrc.jp/alignment/software/) mafft needs to be on your path. Install with something like `brew install maaft` or `apt-get install maaft` and mafft will automatically be placed on your path. 

**Once dependencies are installed:**

1. Download this repository either through git or download and unzip. 
2. Open up lava.sh in a text editor and set the export statements in lines 2-7 to point to where you installed the executable for each of your dependencies. 

You may find it helpful to put lava in your path. This can be done by the command:

`echo 'export PATH=${PATH}:/path/to/lava' >> ~/.bash_profile`

Now you're ready to do some longitudinal analysis of minor alleles! 

# Input

**Mandatory Files**

To run lava you need, at a minimum: (Example files are included in the example folder)
1. fastq files for all of your samples, lava does not perform any adapter or quality trimming so this should be done beforehand. (trimmomatic ect. ). You need at least two samples to perform a meaningful longitudinal analysis. `example1.fastq example2.fastq`
2. A fasta file representing the majority consensus of your first sample. `example_reference.fasta`
3. Either a .gff file with protein annotation for the above reference fasta OR a Genbank accession number pointing to a sample that contains annotations that you would like transferred to your reference fasta `example_reference.gff`
4. A metadata.csv file that must contain two columns: Sample and Passage. Then each the names of every fastq file you want to analyse in the sample column and the passage number or day that the sample on that row was collected. `metadata.csv`


Once you've got all the required files above collected make a new folder and place all the files into this folder. Then execute lava from inside this folder. 

`cd /User/uwvirongs/Downloads/lava/example/` 

`bash ../lava.sh -f example_reference.fasta -g example_reference.gff example1.fastq`

# Usage

To run lava you need to make sure you have placed all the fastq files you want to analyze as well as your metadata.csv file inside a folder. Then you have two choices for running lava

With a reference fasta and a reference gff

`lava.sh -f example_reference.fasta -g example_reference.fasta example1.fastq`

And to pull the reference from Genbank

`lava.sh -f example_reference.fasta -q GENBANK_ACCESSION_NUMBER example1.fastq`

For additional help you can also run 

`lava.sh -h`

# Output Files

Output files will be placed into the same folder you placed all your input in. An interactive graph will be automatically opened on your default browser. This graph is saved as genome_protein_plots.html and sharing is as easy as sending this html file over email (no other files are required once genome_protein_plots.html has been generated).

Additionally you can examine the data more in depth via the merged.csv file which will be created - you can also examine the alignments and read mapping of each of your fastq files be picking the appropriate .bam file. (I.e. if you wanted to see how example1.fastq mapped you can pull example1.bam and examine it yourself.)

# Common Errors
1. `WARNING: A total of 1 sequences will be ignored due to lack of correct ORF annotation`

	This error will occur when annovar finds an error in an open reading frame. This is most often due to an incorrect gff file. (If the start and stop positions for a protein is wrong annovar will detect this and produce this error). Therefore the solution is to make sure the protein start and stops in your gff file are in the correct location nucleotide position should be relative to the matching fasta file. 
	
2. `Can't find annotation record "transcript:3D" referenced by "3D" Parent attribute`

	Once again this error is due to an incorrect gff file. Make sure the gff file has every line with the correct ID= and Parent= , with the correct type (transcript, gene, CDS, etc.) For help formatting your gff files correctly you can look at the included example_reference.gff
	
3. `Exception in thread "main" htsjdk.samtools.SAMFormatException: Error parsing text SAM file. Empty field at position 9 (zero-based); File SC168.sam; Line 1081`

	Make sure there are no special characters (dashes, underscores, etc.) in the reference genome name.````
