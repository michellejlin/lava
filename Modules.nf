/*
 * Define the parameters used in the processes below
 */

/*
 * Define the processes used in this workflow
 */


process generate_report {


    //Retry at most 3 times
    errorStrategy 'retry'
    maxRetries 3
    
    // Define the Docker container used for this step
    container "quay.io/fhcrc-microbiome/clomp:v0.1.3"

    // Define the input files
    input:
      tuple val(base), file(kraken_tsv_list), file(unassigned_txt_list), file(assigned_txt_list)
      file BLAST_CHECK_DB
      file "kraken_db/"
    // Define the output files
    output:
      file "${base}.final_report.tsv"
      file "${base}_unassigned.txt"
      file "${base}_assigned.txt"
    // Code to be executed inside the task
    script:
      """


"""
}