# LAVA: Longitudinal Analysis of Viral Alleles
LAVA analyzes and visualizes minor allele variants in longitudinal sequence data. LAVA takes a reference fasta (normally representing the first sample in your longitudinal analysis), fastq files (for every sample in your analysis), and a metadata sheet (providing information on what day or passage each sample was collected). Output will be displayed as an interactive graph in your web browser. LAVA will only work on Mac and Linux machines. 

# Installation
**Before you get started**

Although we have a script that will do the bulk of the installation for you - this script requires a few things. 

1. Mac or Linux operating system. 
2. An Internet connection. 
3. Python - This should come installed by default on most Mac and Linux operating systems. However if for some reason it isn't you can [download Python](https://www.python.org/).
4. If you are on a Mac operating system you need to have brew installed. Instructions for downloading and installing brew can be found at https://brew.sh/.
5. You also need a java runtime environment if not already installed on your computer. This can be installed with brew (for Mac) with `brew cask install java` and with apt-get (for Linux) via  `sudo apt-get install openjdk-8-jdk`. Both of these commands can be run simply by opening a terminal window and typing them in. For the Linux installation, make sure the package index is updated first by `sudo apt-get update`.

**Installation Instructions**

To make installation easier (there are a lot of dependencies!), we've provided a short script to install most of these for you. 

1. Ensure you have the requirements listed in the 'Before you get started' section.
2. Download the LAVA repository. If you have git installed on the command line you can open up a terminal window and type `git clone https://github.com/michellejlin/lava.git`. Or if you're not comfortable doing that you can click the green button that says Clone or download and click on download zip. Then unzip the folder to wherever you want on your computer. 

3. Open up a terminal window and change directories to the main LAVA folder. If you downloaded the zip file to your downloads folder this would look something like: `cd /Users/username/Downloads/lava-master/` 

4. [Download ANNOVAR](http://www.openbioinformatics.org/annovar/annovar_download_form.php). This program requires registration with a .edu email, or request for access. Once you've recieved access to the annovar download, download and unzip it. Then copy and paste all the files with extension .pl inside the ANNOVAR folder into your main LAVA folder. (There is no need to copy example or humandb.)

5. Run the install script by typing into the terminal window `python install.py -i -c`. The install script will work for a while and install anything you don't already have on your computer. When installation is complete you'll see this message: `All dependencies working properly! Time to do some longitudinal analysis of viral alleles! :DDDDDD.` This means LAVA should work properly!

Now we have to make sure the script can be run from anywhere. Navigate to your main LAVA folder (same as step #3).

6. Type in `pwd`. This will give you the current path. Copy this path.
7. Type in `nano ~/.bashrc`. For MAC OSX users, replace 'bashrc' with 'bash_profile', keeping the rest of the punctuation.
8. This will open up an editor in the terminal window. Scroll down to get to the bottom of this file.
9. Copy paste ``alias lava.py="python PATH/lava.py"`` into the terminal, where you replace PATH with the path you just copied. An example might look like this: `alias lava.py="python /Users/uwvirongs/Downloads/lava/lava.py"`.
10. Hit Ctrl+X to quit out of the editor. Save when prompted (this may mean typing 'Y') and press Enter for any new prompts that show up.
11. Type in `source ~/.bashrc` (or replace 'bashrc' with 'bash_profile' for Mac OSX), to refresh the file.

That's it! Now you're ready to do some longitudinal analysis of minor alleles! 

# Input

**Mandatory Files**

Example files are included in the 'example' folder. It is HIGHLY recommended you examine each of the example files before performing your own analysis.

To run LAVA you need, at a minimum:

1. fastq files for all of your samples. LAVA does not perform any adapter or quality trimming so this should be done beforehand (trimmomatic, etc.). You need at least two samples to perform a meaningful longitudinal analysis. `Example1_file1.fastq Example1_file2.fastq`

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


# Usage

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

# Output Files

Output files will be placed into the same folder you placed all your input in. An interactive graph will be automatically opened on your default browser. This graph is saved as genome_protein_plots.html and sharing is as easy as sending this html file over email (no other files are required once genome_protein_plots.html has been generated).

Another file is created which contains an integrated viewer provided by NGL (provided graciously under the MIT licence - which is the same as this software). The source code is available at https://github.com/arose/ngl.

Additionally you can examine the data more in depth via the merged.csv file which will be created, which includes information such as position, nucleotide changes, allele frequency, depth, and so forth. 

You can also examine the alignments and read mapping of each of your fastq files be picking the appropriate .bam file. (i.e. if you wanted to see how example1.fastq mapped you can pull example1.bam and examine it yourself.) These files are automatically deleted by defualt, but to save them run lava with the `-save` option.

# Common Errors
1. `WARNING: A total of 1 sequences will be ignored due to lack of correct ORF annotation`

	This error will occur when annovar finds an error in an open reading frame. This is most often due to an incorrect gff file. (If the start and stop positions for a protein is wrong annovar will detect this and produce this error). Therefore the solution is to make sure the protein start and stops in your gff file are in the correct location nucleotide position should be relative to the matching fasta file. 
	
2. `Can't find annotation record "transcript:3D" referenced by "3D" Parent attribute`

	Once again this error is due to an incorrect gff file. Make sure the gff file has every line with the correct ID= and Parent= , with the correct type (transcript, gene, CDS, etc.) For help formatting your gff files correctly you can look at the included example_reference.gff.
	
3. `Exception in thread "main" htsjdk.samtools.SAMFormatException: Error parsing text SAM file. Empty field at position 9 (zero-based); File SC168.sam; Line 1081`

	Make sure there are no special characters (dashes, underscores, etc.) in the reference genome name.

4. `Picard not found - lava is being executed from : ` (or with GATK or VARSCAN replacing picard in this error)

	This error is because LAVA couldn't find PICARD in the LAVA folder. Check to see if you ran the install.py script corectly. You can do this by navigating to the LAVA folder and typing in `python install.py -c` to check if all the dependencies are there.

5. ``import: not authorized `subprocess' @ error/constitute.c/WriteImage/1028.``
	
	This error may be accompanied by your mouse cursor turning into a thick plus, and happens because lava.py is not being run in python. Type in `nano ~/.bashrc` or `nano ~/.bash_profile` for Mac OSX, and make sure the alias statement includes python in front of the path, as in `alias lava.py="python /Users/uwvirongs/Downloads/lava/lava.py"`.
	
# GFF Creation Guide
Perhaps the most difficult aspect of running this program is properly formatting your reference fasta and .gff files. In order to have a longitudinal analysis that makes sense, you need to specify a fasta file containing the majority consensus for the first sample. This allows you to examine minor variants in your first sample properly. If you use a fasta that is not representative of your first sample LAVA will Genbank many mutations at 100% allele frequency in your first sample. One potential fix for this is to use the `-q` flag and specify a genbank record that is a reference for your samples. When using the -q flag LAVA will automatically assemble a consensus sequence for your first set of reads and use this as the reference. However, for situations that are not covered by genbank references (For example if you wanted to analyze all Influenza A segments at once) you would need to manually generate your .fasta and .gff files.

In this case you need to use your favorite method of generating a consensus fasta for your first set of reads (we mainly use Geneious). Once this is done you need to make your .gff file. However, ANNOVAR (for reasons I don't fully understand) requires a VERY strict formatting of these gff files. Therefore, I find that the easiest way of generating a new gff file is to edit gene/CDS/transcript names and locations in the provided `example.gff`. 

To do this you first need to have your start and end nucleotide postitions for your protein locations - these must be relative to the start of the provided reference fasta. (So if you're using a reference fasta that you added 200 Ns to the start of the sequence, all protein starts and stops would need to be increased by 200). 

The first two lines of the .gff file are comments, you can safely ignore these. Then all proteins are coded by 3 tab seperated lines. The first column must be the name of your fasta reference sequence. (So if the first line of your reference fasta is `>example` the first column of each row should read `example`. Second column doesn't matter. Then change the start/stop and protein name for your proteins in blocks of 3. 

![Visual Guide](https://github.com/michellejlin/lava/blob/master/GFF_editing_guide.png)

See the picture for a clearer explanation. If you experience any difficulties doing this, or have any other questions about LAVA, feel free to email us at uwvirongs@gmail.com and we'll be happy to help you out!
