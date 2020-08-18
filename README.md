## To make sure you have the latest version of LAVA, please see the [Greninger lab repo](https://github.com/greninger-lab/lava) for both the longitudinal (master branch) and unlongitudinal versions of LAVA.

# LAVA: Longitudinal Analysis of Viral Alleles

![LAVA](https://github.com/vpeddu/lava/workflows/LAVA/badge.svg)

LAVA analyzes and visualizes minor allele variants in longitudinal sequence data. LAVA takes FASTQ files (for every sample in your analysis), a metadata sheet (providing info on what day or passage each sample was collected), and a reference genome either by accession number or your own FASTA and GFF. Output will be displayed as an interactive graph in your web browser, and a table of all mutations of all samples.

# Table of Contents
1. [Installation](#Installation)
2. [Input](#Input)
3. [Usage](#Usage)
4. [Output Files](#Output-Files)
5. [Common Errors](#Common-Errors)
6. [GFF Creation Guide](#GFF-Creation-Guide)

## Installation
**Before you get started**

LAVA is written in [`Nextflow`](https://www.nextflow.io/) and uses `Docker` containers in the interest of ease of use and reproducibility. The only dependencies are `Nextflow` and `Docker`: 

1. If you are on OSX, install `Nextflow` via `Homebrew` by running `brew install nextflow`. If you are on Linux/Windows, find instructions on install `Nextflow` [here](https://www.nextflow.io/docs/latest/getstarted.html). 
2. `Docker` can be installed from [here](https://docs.docker.com/get-docker/).

That's it! Now you're ready to do some longitudinal analysis of minor alleles! 

## Input

**Mandatory Files**

Example files are included in the `test_data` folder. It is HIGHLY recommended you examine each of the example files before performing your own analysis.

To run LAVA you need, at a minimum:

1. FASTQ files for all of your samples. LAVA does not perform any adapter or quality trimming so this should be done beforehand (e.g. with trimmomatic). You need at least two samples to perform a meaningful longitudinal analysis, e.g. `Example1_file1.fastq Example1_file2.fastq`. You only need to specify the first reference FASTQ for LAVA to point at.

2. A fasta file representing the majority consensus of your first sample. There are two options: 
* A reference fasta and a .gff file with protein annotation for the above reference fasta, e.g. `Example1_ref.fasta Example1_ref.gff`.
* A Genbank accession number pointing to a genome that contains annotations that you would like transferred to your reference fasta, e.g. `NC_039477.1`. (This is the Genbank reference for Example 2, included in the example folder.) In this case, the majority consensus fasta will be automatically generated for you.

NOTE: The examples provided are mainly to illustrate how to use either method - fasta and .gff file, or Genbank accession number - can be used effectively. Example 1 uses a provided fasta and .gff file and Example 2 uses `--GENBANK MF795094.1` to pull the reference from GenBank. The examples provided are real data that was used in [this paper](https://mbio.asm.org/content/mbio/9/4/e00898-18.full.pdf). `Example1_file1.fastq` is sample SC1201, and `Example1_file2.fastq` is CUL1201. The data from Example 2 is from [this study](https://www.biorxiv.org/content/10.1101/473405v1), with samples ST107, ST283, and ST709 being `Example2_file1.fastq`, `Example2_file2.fastq`, and `Example2_file3.fastq`, respectively. I have drastically reduced the number of reads in the FASTQ files to make downloading and running these examples extremely fast, so the data does differ from what's presented a bit. However, if you wish to run the full analysis, all files used are publically available on SRA. 

3. A metadata.csv file that must contain two column headers: Sample and Passage. Sample values must be the full path to the FASTQ for the sample. Passage must contain passage number, day the sample was collected, or any other numerical categorical variable. Take a look at the example metadata file `Example1_metadata.csv` to see the formatting.

## Usage

To run LAVA you need to make sure you have placed all the FASTQ files you want to analyze as well as your metadata.csv file inside a folder. You can then use the terminal to execute LAVA from this folder. You have two choices for running LAVA:

* If your computer doesn't have at least 4 cores and 6GB of ram, run your LAVA commands with `-profile testing`

1. With a reference fasta and a reference gff (as seen in Example 1):
	
`nextflow run greninger-lab/lava --OUTDIR test_data/example_1_output/ --FASTA test_data/Example1_ref.fasta --GFF test_data/Example1_ref.gff --CONTROL_FASTQ test_data/Example1_file1.fastq --METADATA test_data/Example1_metadata.csv -with-docker ubuntu:18.04`

2. And to pull the reference from Genbank, this will place all output into a folder named the current data and time (as seen in Example 2):

`nextflow run greninger-lab/lava --OUTDIR test_data/example_2_output/ --GENBANK NC_039477.1 --CONTROL_FASTQ test_data/Example2_file1.fastq --METADATA test_data/Example2_metadata.csv -with-docker ubuntu:18.04`

For additional help you can also run `lava.py -help`: 

      * --CONTROL_FASTQ The fastq reads for the first sample in
                        your longitudinal analysis [REQUIRED]

      *  --METADATA     A two column csv - the first column is the path to all the 
	  					fastqs you wish to include in your analysis.
                    	The second column is the temporal seperation
                        between the samples. This is unitless so you can input
                        passage number, days, or whatever condition your experiment
                        happens to have. [REQUIRED]
        
       * --OUTDIR        Output directory [REQUIRED]
        
	   * --FASTA        Specify a reference fasta with the majority consensus of the
                        control fastq. This option must be used with the -g flag to
                        specify the protein annotations relative to the start of this
                        fasta. [REQUIRED IF NOT --GENBANK]

        * --GFF         Specify a reference gff file with the protein annotations for
                        the reference fasta supplied with the -f flag. This option
                        must be paired with the -f flag. [REQUIRED IF NOT GENBANK]

        * --GENBANK     Provide a Genbank accession number. This record will be used
                        to generate a majority consensus from the control fastq, and
                        this consensus will be annotated from the downloaded genbank
                        record as well. [REQUIRED IF NOT --FASTA + --GFF]

        * --NUC         Results are listed as nucleotide changes not amino acid
                        changes. Do not use with -png. Doesn't work currently

        * --ALLELE_FREQ Specify an allele frequency percentage to cut off - with a
                        minimum of 1 percent - in whole numbers.

        * --PNG         Output results as a png. Do not use with -nuc. Doesn't work currently
        
        * --DEDUPLICATE Optional flag, will perform automatic removal of PCR
                        duplicates via DeDup.
	* --CATEGORICAL	Optional flag, will specify that column under "Passage" in metadata file will be read as categorical variables.

## Output Files

Output files will be placed the folder specified by `--OUTDIR`. An interactive graph will be automatically opened on your default browser. This graph is saved as `LAVA_plots.html`. Sharing the output is as easy as sending this html file over email.

Additionally you can examine the data more in depth via the `final.csv` file which will be created, which includes information such as position, nucleotide changes, allele frequency, depth, and so forth. 

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
	
	This error may be accompanied by your mouse cursor turning into a thick plus, and happens because lava.py is not being run in python. Type in `nano ~/.bashrc` or `nano ~/.bash_profile` for Mac OSX, and make sure the alias statement includes python in front of the path, as in `alias lava.py="python3 /Users/uwvirongs/Downloads/lava/lava.py"`.

6. ``IOError: [Errno 2] No such file or directory: 'Example1_metadata.csv'``

	Make sure your terminal is based in the directory with your metadata file and FASTQ files! You can type `pwd` to check where your terminal is located. 
	
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

