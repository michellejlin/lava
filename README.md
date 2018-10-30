# LAVA: Longitudinal Analysis of Viral Alleles
LAVA analyzes and visualizes minor allele variants in longitudinal sequence data. LAVA takes a reference fasta (normally representing the first sample in your longitudinal analysis), fastq files (for every sample in your analysis), and a metadata sheet (providing information on what day or passage each sample was collected). Output will be displayed as an interactive graph in your web browser. LAVA will only work on Mac and Linux machines. 

# Installation
**Before you get started**

Although we have a script that wil do the bulk of the installation for you - this script requries a few things. 

1. Mac or Linux operating system. 
2. An internet connection. 
3. Python - This should come installed by default on most Mac and Linux operating systems. However if for some reason it isn't you can download it from (https://www.python.org/)
4. If you are on a Mac operating system you need to have brew installed. Instructions for downloading and installing brew can be found at (https://brew.sh/)

**Installation Instructions**

To make installation easier (there are a lot of dependencies!), we've provided a short script to install most of these for you. 

1. Ensure you have the requirements listed in the 'Before you get started section' 
2. Download the LAVA repository. If you have git installed on the command line you can open up a terminal window and type `git clone https://github.com/michellejlin/lava.git`. Or if you're not comfortable doing that you can click the green button that says Clone or download and click on download zip. Then unzip the folder to wherever you want on your computer. 

3. Open up a terminal window and change directories to the main lava folder. If you downloaded the zip file to your downloads folder this would look something like: `cd /Users/username/Downloads/lava-master/' 

4. Run the install script by typing into the terminal window `python install.py`. The install script will work for a while and install anything you don't already have on your computer. When installation is complete you'll see this message `Installation Complete! The only remaining step is to download ANNOVAR and put all of the files into this folder! For instructions on how to do this see the README.`

5. Download ANNOVAR (http://www.openbioinformatics.org/annovar/annovar_download_form.php). This program requires registration with a .edu email, or request for access. Once you've recieved acsess to the annovar download, download and unzip it. Then copy and paste all the files inside the ANNOVAR folder into your main lava folder. 

That's it! Now you're ready to do some longitudinal analysis of minor alleles! 

# Input

**Mandatory Files**

To run LAVA you need, at a minimum: (Example files are included in the example folder)
1. fastq files for all of your samples, LAVA does not perform any adapter or quality trimming so this should be done beforehand. (trimmomatic ect. ). You need at least two samples to perform a meaningful longitudinal analysis. `example1.fastq example2.fastq`
2. A fasta file representing the majority consensus of your first sample. `example_reference.fasta`
3. Either a .gff file with protein annotation for the above reference fasta OR a Genbank accession number pointing to a sample that contains annotations that you would like transferred to your reference fasta `example_reference.gff`
4. A metadata.csv file that must contain two columns: Sample and Passage. Then each of the names of every fastq file you want to analyse in the sample column and the passage number or day that the sample on that row was collected. `metadata.csv`


Once you've got all the required files above collected make a new folder and place all the files into this folder. Then execute LAVA from inside this folder. 

`cd /User/uwvirongs/Downloads/LAVA/example/` 

`bash ../lava.sh -f example_reference.fasta -g example_reference.gff example1.fastq`

# Usage

To run LAVA you need to make sure you have placed all the fastq files you want to analyze as well as your metadata.csv file inside a folder. Then you have two choices for running LAVA:

1. With a reference fasta and a reference gff

`lava.sh -f example_reference.fasta -g example_reference.fasta example1.fastq`

2. And to pull the reference from Genbank

`lava.sh -f example_reference.fasta -q GENBANK_ACCESSION_NUMBER example1.fastq`

For additional help you can also run 

`lava.sh -help`

# Output Files

Output files will be placed into the same folder you placed all your input in. An interactive graph will be automatically opened on your default browser. This graph is saved as genome_protein_plots.html and sharing is as easy as sending this html file over email (no other files are required once genome_protein_plots.html has been generated).

Additionally you can examine the data more in depth via the merged.csv file which will be created, which includes information such as position, nucleotide changes, allele frequency, depth, and so forth. You can also examine the alignments and read mapping of each of your fastq files be picking the appropriate .bam file. (I.e. if you wanted to see how example1.fastq mapped you can pull example1.bam and examine it yourself.)

# Common Errors
1. `WARNING: A total of 1 sequences will be ignored due to lack of correct ORF annotation`

	This error will occur when annovar finds an error in an open reading frame. This is most often due to an incorrect gff file. (If the start and stop positions for a protein is wrong annovar will detect this and produce this error). Therefore the solution is to make sure the protein start and stops in your gff file are in the correct location nucleotide position should be relative to the matching fasta file. 
	
2. `Can't find annotation record "transcript:3D" referenced by "3D" Parent attribute`

	Once again this error is due to an incorrect gff file. Make sure the gff file has every line with the correct ID= and Parent= , with the correct type (transcript, gene, CDS, etc.) For help formatting your gff files correctly you can look at the included example_reference.gff.
	
3. `Exception in thread "main" htsjdk.samtools.SAMFormatException: Error parsing text SAM file. Empty field at position 9 (zero-based); File SC168.sam; Line 1081`

	Make sure there are no special characters (dashes, underscores, etc.) in the reference genome name.
