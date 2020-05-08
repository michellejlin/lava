# LAVA: Longitudinal Analysis of Viral Alleles
LAVA analyzes and visualizes minor allele variants in longitudinal sequence data. LAVA takes a reference fasta (normally representing the first sample in your longitudinal analysis), fastq files (for every sample in your analysis), and a metadata sheet (providing information on what day or passage each sample was collected). Output will be displayed as an interactive graph in your web browser. LAVA will only work on Mac and Linux machines.

# Table of Contents
1. [Installation](#Installation)
2. [Input](#Input)
3. [Usage](#Usage)
4. [Output Files](#Output-Files)
5. [Common Errors](#Common-Errors)
6. [GFF Creation Guide](#GFF-Creation-Guide)

## Installation
**Before you get started**

Although we have a script that will do the bulk of the installation for you - this script requires a few things. 

1. Mac or Linux operating system. 
2. An Internet connection. 
3. Python 3 - This should come installed by default on most Mac and Linux operating systems. However if for some reason it isn't you can [download Python](https://www.python.org/). If you have both Python 2.7 and Python 3, LAVA will automatically choose to use Python 3.
4. If you are on a Mac operating system you need to have brew installed. Instructions for downloading and installing brew can be found at https://brew.sh/.
5. You also need a java runtime environment if not already installed on your computer. This can be installed with brew (for Mac) with `brew cask install java` and with apt-get (for Linux) via  `sudo apt-get install openjdk-8-jdk`. Both of these commands can be run simply by opening a terminal window and typing them in. For the Linux installation, make sure the package index is updated first by `sudo apt-get update`.

**Installation Instructions**

To make installation easier (there are a lot of dependencies!), we've provided a short script to install most of these for you. 

1. Ensure you have the requirements listed in the 'Before you get started' section.
2. Download the LAVA repository. If you have git installed on the command line you can open up a terminal window and type `git clone https://github.com/michellejlin/lava.git`. Or if you're not comfortable doing that you can click the green button that says Clone or download and click on download zip. Then unzip the folder to wherever you want on your computer. 

3. Open up a terminal window and change directories to the main LAVA folder. If you downloaded the zip file to your downloads folder this would look something like: `cd /Users/username/Downloads/lava-master/` 

4. [Download ANNOVAR](http://www.openbioinformatics.org/annovar/annovar_download_form.php). This program requires registration with a .edu email, or request for access. Once you've recieved access to the annovar download, download and unzip it. Then copy and paste all the files with extension .pl inside the ANNOVAR folder into your main LAVA folder. (There is no need to copy example or humandb.)

5. Run the install script by typing into the terminal window `python3 install.py -i -c`. The install script will work for a while and install anything you don't already have on your computer. When installation is complete you'll see this message: `All dependencies working properly! Time to do some longitudinal analysis of viral alleles! :DDDDDD.` This means LAVA should work properly!
The dependencies that LAVA currently use are:
	* Python modules: biopython, numpy, pandas, bokeh
	* Picard, GATK, VarScan, Annovar
	* bedtools, samtools, bwa, mafft, bcftools

Now we have to make sure the script can be run from anywhere. Navigate to your main LAVA folder (same as step #3).

6. Type in `pwd`. This will give you the current path. Copy this path.
7. Type in `nano ~/.bashrc`. For MAC OSX users, replace 'bashrc' with 'bash_profile', keeping the rest of the punctuation.
8. This will open up an editor in the terminal window. Scroll down to get to the bottom of this file.
9. Type in ``alias lava.py="python INSERTPATHHERE/lava.py"`` into the terminal, where you replace INSERTPATHHERE with the path you just copied. An example might look like this: `alias lava.py="python /Users/uwvirongs/Downloads/lava/lava.py"`.
10. Hit Ctrl+X to quit out of the editor. Save when prompted (this may mean typing 'Y') and press Enter for any new prompts that show up.
11. Type in `source ~/.bashrc` (or replace 'bashrc' with 'bash_profile' for Mac OSX), to refresh the file.

That's it! Now you're ready to do some longitudinal analysis of minor alleles! 

## Input

**Mandatory Files**

Example files are included in the 'example' folder. It is HIGHLY recommended you examine each of the example files before performing your own analysis.

To run LAVA you need, at a minimum:

1. fastq files for all of your samples. LAVA does not perform any adapter or quality trimming so this should be done beforehand (trimmomatic, etc.). You need at least two samples to perform a meaningful longitudinal analysis. `Example1_file1.fastq Example1_file2.fastq`. You only need to specify the first reference fastq for LAVA to point at.

2. A fasta file representing the majority consensus of your first sample. There are two options: 1) A reference fasta and a .gff file with protein annotation for the above reference fasta. `Example1_ref.fasta Example1_ref.gff` OR 2) a Genbank accession number pointing to a sample that contains annotations that you would like transferred to your reference fasta `NC_039477.1` (This is the Genbank reference for Example 2, included in the example folder.) 

NOTE: The examples provided are mainly to illustrate how to use either method - fasta and gff file, or Genbank accession number - can be used effectively. Example 1 uses a provided fasta and gff file and Example 2 uses `-q MF795094.1` to pull the reference from GenBank. The examples provided is real data that was used in [this paper](https://mbio.asm.org/content/mbio/9/4/e00898-18.full.pdf). Example1_file1 is sample SC1201, and Example1_file2 is CUL1201. The data from Example 2 is from [this study](https://www.biorxiv.org/content/10.1101/473405v1), with samples ST107, ST283, and ST709 being Example2_file1, Example2_file2, and Example2_file3, respectively. I have drastically reduced the number of reads in the fastq files to make downloading and running these examples extremely fast, so the data does differ from what's presented a bit. However, if you wish to run the full analysis, all files used are publically available on SRA. 

3. A metadata.csv file that must contain two columns: Sample and Passage. Then each of the names of every fastq file you want to analyse in the sample column and the passage number or day that the sample on that row was collected. Take a look at the example metadata file to see the formatting. `Example1_metadata.csv` 


Once you've got all the required files above collected make a new folder and place all the files into this folder. (This has obviously already been done for the example files). 

NOTE: The metadata file MUST be in the folder LAVA is executed from, because LAVA takes the sample names from the metadata file to search for the .fastqs. This means the metadata file can be in the folder and contain sample names like "Sample1.fastq, Sample2.fastq" or the sample names can be paths to the different fastq files.

We strongly recommend running the examples to see how it works. To do this, navigate to the example folder inside your downloaded LAVA folder. For example, this may look like: `cd /User/uwvirongs/Downloads/LAVA/example/`.

To run Example 1:

	lava.py -f Example1_ref.fasta -g Example1_ref.gff Example1_file1.fastq Example1_metadata.csv -o Example1_output

To run Example 2:

	lava.py -q NC_039477.1 Example2_file1.fastq Example2_metadata.csv -o Example2_output


## Usage

To run LAVA you need to make sure you have placed all the fastq files you want to analyze as well as your metadata.csv file inside a folder. You can then use the terminal to execute LAVA from this folder. You have two choices for running LAVA:

1. With a reference fasta and a reference gff, with the optional -o argument placing output into a folder named output (as seen in Example 1):
	
	`lava.py -f example_reference.fasta -g example_reference.gff example-P0.fastq metadata.csv -o output`

2. And to pull the reference from Genbank, this will place all output into a folder named the current data and time (as seen in Example 2):

	`lava.py -q GENBANK_ACCESSION_NUMBER example-P0.fastq metadata.csv`

Other optional arguments include:

Examining nucleotide changes by type (A -> C, etc.) with the -nuc argument:

	lava.py -f example_reference.fasta -g example_reference.gff example-P0.fastq metadata.csv -nuc -o output

Removing PCR dupicates from reads with the -dedup argument:

	lava.py -f example_reference.fasta -g example_reference.gff example-P0.fastq metadata.csv -dedup -o output

For additional help you can also run `lava.py -help`.

## Output Files

Output files will be placed into the same folder you placed all your input in. An interactive graph will be automatically opened on your default browser. This graph is saved as genome_protein_plots.html and sharing is as easy as sending this html file over email (no other files are required once genome_protein_plots.html has been generated).

Additionally you can examine the data more in depth via the merged.csv file which will be created, which includes information such as position, nucleotide changes, allele frequency, depth, and so forth. 

You can also examine the alignments and read mapping of each of your fastq files be picking the appropriate .bam file. (i.e. if you wanted to see how example1.fastq mapped you can pull example1.bam and examine it yourself.) These files are automatically deleted by defualt, but to save them run lava with the `-save` option.

## Common Errors
1. `WARNING: A total of 1 sequences will be ignored due to lack of correct ORF annotation`

	This error will occur when annovar finds an error in an open reading frame. This is most often due to an incorrect gff file. (If the start and stop positions for a protein is wrong annovar will detect this and produce this error). Therefore the solution is to make sure the protein start and stops in your gff file are in the correct location nucleotide position should be relative to the matching fasta file. This will also occur when there is a early stop codon.
	
2. `Can't find annotation record "transcript:3D" referenced by "3D" Parent attribute`

	Once again this error is due to an incorrect gff file. Make sure the gff file has every line with the correct ID= and Parent= , with the correct type (transcript, gene, CDS, etc.) For help formatting your gff files correctly you can look at the included example_reference.gff.
	
3. `Exception in thread "main" htsjdk.samtools.SAMFormatException: Error parsing text SAM file. Empty field at position 9 (zero-based); File SC168.sam; Line 1081`

	Make sure there are no special characters (dashes, underscores, etc.) in the reference genome name.

4. `Picard not found - lava is being executed from : ` (or with GATK or VARSCAN replacing picard in this error)

	This error is because LAVA couldn't find PICARD in the LAVA folder. Check to see if you ran the install.py script corectly. You can do this by navigating to the LAVA folder and typing in `python install.py -c` to check if all the dependencies are there.

5. ``import: not authorized `subprocess' @ error/constitute.c/WriteImage/1028.``
	
	This error may be accompanied by your mouse cursor turning into a thick plus, and happens because lava.py is not being run in python. Type in `nano ~/.bashrc` or `nano ~/.bash_profile` for Mac OSX, and make sure the alias statement includes python in front of the path, as in `alias lava.py="python /Users/uwvirongs/Downloads/lava/lava.py"`.

6. ``IOError: [Errno 2] No such file or directory: 'Example1_metadata.csv'``

	Make sure your terminal is based in the directory with your metadata file and fastq files! You can type `pwd` to check where your terminal is located. 
	
7. ``WARNING: Unable to retrieve regions at Example1_ref due to lack of sequence information. WARNING: Cannot identify sequence for transcript:L (starting from Example1_ref:8590). WARNING: Cannot identify sequence for gene:L (starting from Example1_ref:8590)``

	Examine your reference fasta and reference gff files. Make sure that internally, both have the same name as the file name. For example, the first column of Example1_ref.gff must say Example1_ref.
	
8. *My genome_protein_plots and merged.csv file is blank. Help!*

	This can be caused by a number of different factors. Most likely it is a GFF issue. Make sure you are following the GFF creation guide closely. Some of the most common problems that should be fixed are: 
	* Make sure you have gene, CDS, and transcript for all your proteins.
	* CDS should be in all caps.
	* For all CDS, make sure there is a `Parent=transcript:`, and for all transcripts, make sure there is a `Parent=gene:` in the last column. Look at the GFF guide for how to format this.
	* For all CDS, make sure the phase is set to "0". This is the second to last column. Everyything else (genes, transcripts) should have the phase set to "."


## GFF Creation Guide
Perhaps the most difficult aspect of running this program is properly formatting your reference fasta and .gff files. In order to have a longitudinal analysis that makes sense, you need to specify a fasta file containing the majority consensus for the first sample. This allows you to examine minor variants in your first sample properly. If you use a fasta that is not representative of your first sample LAVA will Genbank many mutations at 100% allele frequency in your first sample. 

In order to avoid this issue, we recommend using the `-q` flag to specify a GenBank record that is a reference for your samples. Assuming the selected reference has accurate annotations, LAVA will automatically assemble a working consensus sequence for your first set of reads and use this as the reference. 

However, for situations that are not covered by GenBank references, you would need to manually generate your own .fasta and .gff files.

**Example: Using the Template GFF**

An example of something that would not be covered by GenBank references, and thus would not be recommended to use the `-q` flag, is if you wanted to analyze all Influenza A segments at once.

In this case you need to use your favorite method of generating a consensus fasta for your first set of reads (we mainly use Geneious). Once this is done you need to make your .gff file. However, ANNOVAR requires a VERY strict formatting of these gff files. 

The easiest way of generating a new gff file is to edit gene/CDS/transcript names and locations in the provided `Example1_ref.gff`. We highly recommend this method to avoid lots of formatting issues.

To do this you first need to have your start and end nucleotide postitions for your protein locations - these must be relative to the start of the provided reference fasta. (So if you're using a reference fasta that you added 200 Ns to the start of the sequence, all protein starts and stops would need to be increased by 200). 

Below is a visual guide of what you need to change the template `Example1_ref.gff` for your purposes.

![Visual Guide](https://github.com/michellejlin/lava/blob/master/Gff-editing-guide.png)

The color coded regions are the only portions that need to be changed. (The first two lines of the .gff files are comments and can be safely ignored.)

Proteins are coded by 3 tab separated lines (gene, CDS, transcript). 

* The first column must be the name of your .fasta reference sequence. (So if the first line of your reference .fasta is >example, the first column of each row should read `example`. Here, the name is `WSN_reference2`.
* The fourth column should contain your start and end nucleotide positions of protein locations. Change the numbers to match your fasta file. Make sure to keep the blocks of 3, so that there are correct protein locations for each of gene, CDS, transcript.
* The last column has protein names that need to be replaced. Make sure you are replacing after both the `ID=` and the `Parent=`.

**Specific Requirements**

If you don't want to use the template GFF, or want to troubleshoot any problems with the GFF that may be popping up, here are the specific formatting requirements for each column.

| Fasta Name | Source | Feature | Start | End | Score | Strand | Phase | Attributes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Must match the name of your .fasta reference sequence: both the first line, and the file name. | Anything in this column. | One of 3 things: gene, CDS, transcript. CDS must be in all caps. Each protein MUST have all 3 features. | Beginning position of the protein. | End position of the protein. | Only contains "." | Only contains "+". | Contains a "0" for all CDS lines, and "." for all others. | Contains ID=`feature type`, where `feature type` is one of gene, CDS, or transcript, followed by the protein name. For CDS lines, it must also contain a `Parent=transcript:` identifier, followed by the protein name. For transcript lines, it must also contain a `Parent=gene:` identifier, followed by the protein name. All lines must end with `biotype=protein_coding`. Each of these tags should be separated by semicolons.|

If you experience any difficulties doing this, or have any other questions about LAVA, feel free to email us at uwvirongs@gmail.com and we'll be happy to help you out!

### Example 1 nextflow command 
`nextflow run vpeddu/lava --OUTDIR test_data/example_1_output/ --FASTA test_data/Example1_ref.fasta --GFF test_data/Example1_ref.gff --CONTROL_FASTQ test_data/Example1_file1.fastq --METADATA test_data/Example1_metadata.csv -with-docker ubuntu:18.04`

### Example 2 nextflow command 
`nextflow run vpeddu/lava --OUTDIR test_data/example_2_output/ --GENBANK NC_039477.1 --CONTROL_FASTQ test_data/Example2_file1.fastq --METADATA test_data/Example2_metadata.csv -with-docker ubuntu:18.04`

