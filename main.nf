#!/usr/bin/env nextflow
/*
========================================================================================
                         LAVA
========================================================================================
 Longitudinal Analysis of Viral Alleles
 #### Homepage / Documentation
 https://github.com/vpeddu/LAVA/tree/nextflow
----------------------------------------------------------------------------------------
*/

// Using the Nextflow DSL-2 to account for the logic flow of this workflow
nextflow.preview.dsl=2


def helpMessage() {
    log.info"""
    LAVA: Longitudinal Analysis of Viral Alleles

    Usage:

    An example command for running the pipeline is as follows:

    nextflow run FredHutch/CLOMP \\
        --INPUT_FOLDER reads/ \\
        --PAIRED_END \\
        --HOST_FILTER_FIRST \\
        --OUTDIR output/
        
  
        --CONTROL_FASTQ The fastq reads for the first sample in
                        your longitudinal analysis [REQUIRED]

        --METADATA       Required argument: A two column csv - the first column is the
                        name of all the fastqs you wish to include in your analysis.
                        All fastqs that you want to include need to be specified in
                        this file AND be located in the folder from which you are
                        running lava. The second column is the temporal seperation
                        between the samples. This is unitless so you can input
                        passage number, days, or whatever condition your experiment
                        happens to have. [REQUIRED]
        
        --OUTDIR        Output directory
        
        --FASTA         Specify a reference fasta with the majority consensus of the
                        control fastq. This option must be used with the -g flag to
                        specify the protein annotations relative to the start of this
                        fasta.

        --AF            pecify an allele frequency percentage to cut off - with a minimum of 1 percent - in whole numbers. default = ' '

        --GFF           Specify a reference gff file with the protein annotations for
                        the reference fasta supplied with the -f flag. This option
                        must be paired with the -f flag.

        --GENBANK       Provide a Genbank accession number. This record will be used
                        to generate a majority consensus from the control fastq, and
                        this consensus will be annotated from the downloaded genbank
                        record as well.

        --NUC           Results are listed as nucleotide changes not amino acid
                        changes. Do not use with -png.

        --ALLELE_FREQ   Specify an allele frequency percentage to cut off - with a
                        minimum of 1 percent - in whole numbers.

        --PNG           Output results as a png. Do not use with -nuc.
        
        --DEDUPLICATE   Optional flag, will perform automatic removal of PCR
                        duplicates via DeDup.

        --OUT           Optional flag to name the output folder that lava will stuff
                        output into. If a name isn't provided folder will be named
                        lava-date.

        --SAVE          Optional argument to save intermediate alignment files (sams,
                        bams, vcfs, ect) LAVA's default behavior is to remove these
                        after use to save disk footprint.

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Set defaults
// These values are overridden by what the user specifies (e.g. with --R1)
/*
params.INPUT_FOLDER = false
params.INPUT_SUFFIX = ".fastq.gz"
params.PAIRED_END = false
params.OUTDIR = false
params.SNAP_INDEXES = "s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.00/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.01/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.02/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.03/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.04/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.05/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.06/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.07/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.08/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.09/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.10/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.11/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.12/,s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/full_clomp_db/nt.13/" 
params.SNAP_OPTIONS = "-mrl 65 -d 9 -h 30000 -om 1 -omax 20"
params.HOST_FILTER_FIRST = false
params.SECOND_PASS = false
params.TRIMMOMATIC_OPTIONS = ':2:30:10 HEADCROP:10 SLIDINGWINDOW:4:20 CROP:65 MINLEN:65'
params.BBDUK_TRIM_OPTIONS = 'ktrim=r k=27 hdist=1 edist=0 mink=4 qtrim=rl trimq=6 minlength=65 ordered=t qin=33'
params.TRIMMOMATIC_JAR_PATH = "s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/trimmomatic-0.38.jar"
params.TRIMMOMATIC_ADAPTER_PATH = "s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/adapters.fa"
params.SEQUENCER = 'ILLUMINACLIP:'
params.BWT_DB_PREFIX = "hg38"
params.BWT_DB_LOCATION = "s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/hg38/"
params.BWT_SECOND_PASS_OPTIONS = '-D 35 -R 5 -N 1 -L 19 -i S,1,0.50'
params.BLAST_EVAL = 0.001
params.BLAST_CHECK = false
params.DEDUPE = true
params.BLAST_CHECK_DB = false
params.FILTER_LIST = "[12908,28384,48479]"
params.KRAKEN_DB_PATH = "s3://fh-ctr-public-reference-data/tool_specific_data/CLOMP/kraken_db/"
params.H_STRICT = false
params.H_TAXID = 9606
params.LOGIC = "strict"
params.INCLUSION_TAXID = 2759
params.EXCLUSION_TAXID = 9604
params.ENTREZ_EMAIL = "uwvirongs@gmail.com"
params.BASE_DELIMITER = "_"
params.MIN_READ_CUTOFF = 10
params.SAM_NO_BUILD_LIST = "[2759,77133]"
params.ADD_HOST_FILTERED_TO_REPORT = true
params.HOST_FILTER_TAXID = 9606
params.WRITE_UNIQUES = true
params.EDIT_DISTANCE_OFFSET = 6
params.BUILD_SAMS = false
params.SNAP_BATCHSIZE = 20
params.TIEBREAKING_CHUNKS = 2
*/
params.OUTDIR= false

// Check to make sure that the required parameters have been set
//if (!params.INPUT_FOLDER){ exit 1, "Must provide folder containing input files with --INPUT_FOLDER" }
if (!params.OUTDIR){ exit 1, "Must provide output directory with --OUTDIR" }
if (!params.CONTROL_FASTQ){ exit 1, "Must provide control FASTQ with --ControlFastq" }
if (!params.INPUT_FOLDER){ exit 1, "Must provide input folder with --INPUT_FOLDER" }

params.GENBANK = 'False'
params.GFF = 'False'
params.FASTA = 'NO_FILE'
params.DEDUPLICATE = 'false' 
params.AF = ' '
// Make sure the output directory has a trailing slash
if (!params.OUTDIR.endsWith("/")){
    params.OUTDIR = "${params.OUTDIR}/"
}
//if (params.GENBAKN ==)
// Identify some resource files
METADATA_FILE = file(params.METADATA)

/*
 * Import the processes used in this workflow
 */

//include run_lava from './Modules.nf'
//inlcude fml from './Modules.nf'
include CreateGFF from './Modules.nf'
include Alignment_prep from './Modules.nf'
include Align_samples from './Modules.nf' 
include Pipeline_prep from './Modules.nf'
include Create_VCF from './Modules.nf'
include Ref_done from './Modules.nf'
include Extract_variants from './Modules.nf'
include Annotate_complex from './Modules.nf'
include Annotate_complex_first_passage from './Modules.nf'
include Generate_output from './Modules.nf'



CONTROL_FASTQ = file(params.CONTROL_FASTQ)
FASTA = file(params.FASTA)
 input_read_ch = Channel
 .fromPath("${params.INPUT_FOLDER}*fastq")



input_read_ch = Channel
    .fromPath(METADATA_FILE)
    .splitCsv(header:false)

input_read_ch = Channel
    .fromPath(METADATA_FILE)
    .splitCsv(header:true)
    .map{ row-> tuple(file(row.Sample), (row.Passage)) }

// TODO: logic to check --FASTA and --GFF are in together if no --GENBANK
// Run the workflow
workflow {
        //fml() 

        CreateGFF ( 
            params.GENBANK, 
            CONTROL_FASTQ,
            file(params.FASTA),
            file(params.GFF)
        )
        
        Alignment_prep ( 
            CreateGFF.out[0],
            CreateGFF.out[1],
            CreateGFF.out[2]
        )

        Align_samples ( 
            input_read_ch,
            Alignment_prep.out[0],
            params.INPUT_FOLDER,
            input_read_ch.first()
            
        )

        Pipeline_prep ( 
            Align_samples.out[0].collect(),
            CreateGFF.out[2],
            //INITIALIZE_MERGED_CSV,
            CreateGFF.out[3],
            Alignment_prep.out[0]
        )

        Create_VCF ( 
            CreateGFF.out[3],
            Pipeline_prep.out[2],
            Align_samples.out[0],
            Alignment_prep.out[1],
            input_read_ch.first(),
            Alignment_prep.out[2]
        )

        Ref_done ( 
            input_read_ch.first(),
            params.AF,
            Create_VCF.out[0],
            CreateGFF.out[3],
            Pipeline_prep.out[3],
            Align_samples.out[1],
            METADATA_FILE
        )

        Extract_variants ( 
            input_read_ch.first(),
            Create_VCF.out[1],
            METADATA_FILE

        )

        Annotate_complex( 
            Extract_variants.out[0]//,
            //ANNOTATE_COMPLEX

        )

        Annotate_complex_first_passage( 
            Ref_done.out[0],
            //ANNOTATE_COMPLEX
        )

        Generate_output( 
            Annotate_complex_first_passage.out,
            Annotate_complex.out[0].collect(),
            Annotate_complex.out[1].collect(),
            Annotate_complex.out[2].collect(),
            Annotate_complex.out[3].collect(),
            Pipeline_prep.out[0],
            Pipeline_prep.out[1],
            //GENOME_PROTEIN_PLOTS,
            Align_samples.out[2].collect()
        )
        
    publish:
        Generate_output.out to: "${params.OUTDIR}"
        //filter_human_single.out[1] to: "${params.OUTDIR}/logs/"
}
