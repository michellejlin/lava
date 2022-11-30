#!/usr/bin/env nextflow
/*
========================================================================================
                         LAVA
========================================================================================
 Longitudinal Analysis of Viral Alleles
 #### Homepage / Documentation
https://github.com/greninger-lab/lava
----------------------------------------------------------------------------------------
*/

// DSL2 enable syntax depending on different nextflow version
// nextflow_dsl2_v = '20.07.1'
// if ( nextflow.version.matches(">= $nextflow_dsl2_v") ) {
//     nextflow.enable.dsl=2
// } else {
//     nextflow.preview.dsl=2
// }
nextflow.enable.dsl=2

def helpMessage() {
    log.info"""
    RAVA: Reference-based Analysis of Viral Alleles
    Usage:
    An example command for running the pipeline is as follows:
    nextflow run greninger-lab/lava \\
        --METADATA      Required argument: A two column csv - the first column is the
                        path to all the fastqs you wish to include in your analysis.
                        All fastqs that you want to include need to be specified in
                        this file AND be located in the folder from which you are
                        running lava. The second column is the temporal seperation
                        between the samples. This is unitless so you can input
                        passage number, days, or whatever condition your experiment
                        happens to have. [REQUIRED]

        --OUTDIR        Output directory [REQUIRED]

        --FASTA         Specify a reference fasta with the majority consensus of the
                        control fastq. This option must be used with the -g flag to
                        specify the protein annotations relative to the start of this
                        fasta. [REQUIRED IF NOT --GENBANK]
        --GFF           Specify a reference gff file with the protein annotations for
                        the reference fasta supplied with the -f flag. This option
                        must be paired with the -f flag. [REQUIRED IF NOT GENBANK]
        --GENBANK       Provide a Genbank accession number. This record will be used
                        to generate a majority consensus from the control fastq, and
                        this consensus will be annotated from the downloaded genbank
                        record as well. [REQUIRED IF NOT --FASTA + --GFF]
        --AF            pecify an allele frequency percentage to cut off
                        - with a minimum of 1 percent - in whole numbers. default = ' '
        --NUC           Results are listed as nucleotide changes not amino acid
                        changes. Do not use with -png.
        --ALLELE_FREQ   Specify an allele frequency percentage to cut off - with a
                        minimum of 1 percent - in whole numbers.
        --PNG           Output results as a png. Do not use with -nuc.

        --DEDUPLICATE   Optional flag, will perform automatic removal of PCR
                        duplicates via DeDup.
        --NAME          Changes graph name
    """.stripIndent()
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*          SET UP CONFIGURATION VARIABLES            */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

params.OUTDIR= false
params.GENBANK = 'False'
params.GFF = 'False'
params.FASTA = 'NO_FILE'
params.DEDUPLICATE = 'false'
params.ALLELE_FREQ = 'NO_VAL'

params.NAME = 'false'

METADATA_FILE = file(params.METADATA)

/*
 * Import the processes used in this workflow
 */

include { CreateGFF_Genbank }from './Modules.nf'
include { CreateGFF } from './Modules.nf'
include { Alignment_prep } from './Modules.nf'
include { Align_samples } from './Modules.nf'
include { Pipeline_prep } from './Modules.nf'
include { Create_VCF } from './Modules.nf'
include { Extract_variants } from './Modules.nf'
include { Annotate_complex } from './Modules.nf'
include { Generate_output } from './Modules.nf'


// Staging python scripts
PULL_ENTREZ = file("$workflow.projectDir/bin/pull_entrez.py")
WRITE_GFF = file("$workflow.projectDir/bin/write_gff.py")
WRITE_GFF2 = file("$workflow.projectDir/bin/write_gff2.py")

INITIALIZE_PROTEINS_CSV = file("$workflow.projectDir/bin/initialize_proteins_csv.py")
ANNOTATE_COMPLEX_MUTATIONS = file("$workflow.projectDir/bin/Annotate_complex_mutations.py")
MAT_PEPTIDE_ADDITION = file("$workflow.projectDir/bin/mat_peptide_addition.py")
RIBOSOMAL_SLIPPAGE = file("$workflow.projectDir/bin/ribosomal_slippage.py")
GENOME_PROTEIN_PLOTS = file("$workflow.projectDir/bin/genome_protein_plots.py")
PALETTE = file("$workflow.projectDir/bin/palette.py")
ANNOVAR2 = file("$workflow.projectDir/bin/ANNOVAR_Replacement.py")

//input_read_ch = Channel


// Error handling for input flags
//if OUTDIR not set
if (params.OUTDIR == false) {
    println( "Must provide an output directory with --OUTDIR")
    exit(1)
}
// If --GENBANK and --FASTA or --GFF are specified at the same time
// if(((params.GENBANK != "False") && (params.FASTA != "NO_FILE"))){
//     println("--GENBANK cannot be used with --FASTA or --GFF")
//     exit(1)
// }
// if(((params.GENBANK != "False") && (params.GFF != "False"))){
//     println("--GENBANK cannot be used with --FASTA or --GFF")
//     exit(1)
// }
// // If --FASTA without --GENBANK or vice versa
// if( (params.FASTA != "NO_FILE") && params.GFF == 'False'){
//     println('--GFF needs to be specified with --FASTA')
//     exit(1)
// }
// if( (params.GFF != "False") && params.FASTA == 'NO_FILE'){
//     println('--FASTA needs to be specified with --GFF')
//     exit(1)
// }
// If no flags specified
// if(params.GENBANK == "False"){
//     println('Must provide --GENBANK flag with GenBank accession number.')
//     exit(1)
// }

// Make sure OUTDIR ends with trailing slash
if (!params.OUTDIR.endsWith("/")){
   params.OUTDIR = "${params.OUTDIR}/"
}

input_read_ch = Channel
    .fromPath(METADATA_FILE)
    .splitCsv(header:true)
    .map{ row-> tuple(file(row.Sample), (row.Passage)) }

// Throws exception if paths in METADATA are not valid
Channel
    .fromPath(METADATA_FILE)
    .splitCsv(header:true)
    .map{row-> (file(row.Sample, checkIfExists:true))}.ifEmpty{error "Check metadata file"}
    //.map{row-> (file(row.Sample).isEmpty())}
    //.filter{ it == false}.subscribe{println it}


////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                 RUN THE WORKFLOW                   */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

workflow {
    log.info nfcoreHeader()
    if(params.FASTA == "NO_FILE") {
        CreateGFF_Genbank (
            params.GENBANK,
            PULL_ENTREZ,
            WRITE_GFF
        )

        Alignment_prep (
            CreateGFF_Genbank.out[0],
            CreateGFF_Genbank.out[1],
            CreateGFF_Genbank.out[2],
            CreateGFF_Genbank.out[3],
            CreateGFF_Genbank.out[4]
        )
    }
    else {
        CreateGFF (
            params.GENBANK,
            PULL_ENTREZ,
            WRITE_GFF2,
            file(params.FASTA),
            file(params.GFF)
        )

        Alignment_prep (
            CreateGFF.out[0],
            CreateGFF.out[1],
            CreateGFF.out[2],
            CreateGFF.out[3],
            CreateGFF.out[4]
        )
    }

        Align_samples (
            input_read_ch,
            Alignment_prep.out[0],
            params.DEDUPLICATE

        )

        Pipeline_prep (
            Align_samples.out[0].collect(),
            Alignment_prep.out[3],
            Alignment_prep.out[0],
            INITIALIZE_PROTEINS_CSV
        )

        Create_VCF (
            Align_samples.out[0],
            Alignment_prep.out[1],
            Alignment_prep.out[2],
            ANNOVAR2,
            file(params.FASTA),
            file(params.GFF),
            CreateGFF.out[2]
        )

        Extract_variants (
            Create_VCF.out[1],
            METADATA_FILE

        )

        Annotate_complex(
            Extract_variants.out[0],
            ANNOTATE_COMPLEX_MUTATIONS
        )

        Generate_output(
            Annotate_complex.out[0].collect(),
            Annotate_complex.out[1].collect(),
            Annotate_complex.out[2].collect(),
            Annotate_complex.out[3].collect(),
            Pipeline_prep.out[0],
            Pipeline_prep.out[1],
            Align_samples.out[1].collect(),
            Create_VCF.out[2].collect(),
            Alignment_prep.out[4],
            Alignment_prep.out[5],
            MAT_PEPTIDE_ADDITION,
            RIBOSOMAL_SLIPPAGE,
            GENOME_PROTEIN_PLOTS,
            PALETTE
        )
}

def nfcoreHeader() {

    return """
                       ooO
                     ooOOOo
                   oOOOOOOoooo
                 ooOOOooo  oooo
                /vvv\\
               /V V V\\
              /V  V  V\\
             /         \\            oh wow  look at these alleles
            /           \\          /         /
          /               \\   	  o          o
__       /                 \\     /-   o     /-
/\\     /                     \\  /\\  -/-    /\\
                                    /\\


        `7MM"'"Mq.        db `7MMF'   `7MF' db
          MM   `MM.      ;MM:  `MA     ,V  ;MM:
          MM   ,M9      ,V^MM.  VM:   ,V  ,V^MM.
          MMmmdM9      ,M  `MM   MM.  M' ,M  `MM
          MM  YM.      AbmmmqMA  `MM A'  AbmmmqMA
          MM   `Mb.   A'     VML  :MM;  A'     VML
        .JMML. .JMM..AMA.   .AMMA. VF .AMA.   .AMMA.

                      Version 3
    """.stripIndent()
}
