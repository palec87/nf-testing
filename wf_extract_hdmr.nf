#!/usr/bin/env nextflow

/* 
The pipeline does the following steps
1. finds the association between the archive name and ref_code in the batch1and2 csv file (probably python)
2. create a correctly named floder and copy results
3. write a log if some of the files are missing.
*/

include { readYAML } from './modules/read_yaml.nf'
include { unzipArchive } from './modules/unzip_archive.nf'
include { mergeI5 } from './modules/merge_i5_chunks.nf'


// list of files
params.files = "inp_files.csv"
params.folder_extracted_tables = "${projectDir}/results-hcmr"

// setup for CCMAR data
// params.archives_root = "/usr/local/scratch/metaGOflow-COMPLETED-results/Batch1and2/CCMAR-data"  // archive folder redi

// setup for HCMR data
params.archives_root = "/usr/local/scratch/metaGOflow-COMPLETED-results/Batch1and2/HCMR-data/data"  // archive folder redi


// I think this is the solution for the renaming problems too, or????
process extractTables {
    conda '/usr/local/scratch/nf-metaGOflow/wf-test/nf-testing/conda.yaml'
    publishDir "results-hcmr", mode: 'copy'
    debug true

    input:
    path root_folder
    path archive_name
    val out_name
    
    output:
    path "${out_name}-tables"
    
    script:
    """
    mkdir -p ${archive_name}-tables
    cp ${archive_name}/results/functional-annotation/*.merged.summary.go ${out_name}-tables
    cp ${archive_name}/results/functional-annotation/*.merged.summary.go_slim ${out_name}-tables
    cp ${archive_name}/results/functional-annotation/*.merged.summary.ips ${out_name}-tables
    cp ${archive_name}/results/functional-annotation/*.merged.summary.ko ${out_name}-tables
    cp ${archive_name}/results/functional-annotation/*.merged.summary.pfam ${out_name}-tables
    # cp ${archive_name}/results/functional-annotation/*.merged_CDS.I5.tsv.gz ${out_name}-tables
    cp ${archive_name}/results/taxonomy-summary/LSU/*.merged_LSU.fasta.mseq.tsv ${out_name}-tables
    cp ${archive_name}/results/taxonomy-summary/SSU/*.merged_SSU.fasta.mseq.tsv ${out_name}-tables
    """
}

process findDatasets {
    conda '/usr/local/scratch/nf-metaGOflow/wf-test/nf-testing/conda.yaml'
    // publishDir "results-tables", mode: 'copy'
    debug true

    input:
    path archive_name
    path target_directory
    
    output:
    path "${target_directory}-tables"
    val true, emit: trigger
    
    script:
    """
    mkdir -p ${target_directory}-tables
    cp ${archive_name}/results/functional-annotation/*.merged.summary.go ${target_directory}-tables
    cp ${archive_name}/results/functional-annotation/*.merged.summary.go_slim ${target_directory}-tables
    cp ${archive_name}/results/functional-annotation/*.merged.summary.ips ${target_directory}-tables
    cp ${archive_name}/results/functional-annotation/*.merged.summary.ko ${target_directory}-tables
    cp ${archive_name}/results/functional-annotation/*.merged.summary.pfam ${target_directory}-tables
    # cp ${archive_name}/results/functional-annotation/*.merged_CDS.I5.tsv.gz ${target_directory}-tables
    cp ${archive_name}/results/taxonomy-summary/LSU/*.merged_LSU.fasta.mseq.tsv ${target_directory}-tables
    cp ${archive_name}/results/taxonomy-summary/SSU/*.merged_SSU.fasta.mseq.tsv ${target_directory}-tables
    """
}

workflow {
    def python_dir = file("/usr/local/scratch/nf-metaGOflow/wf-test/nf-testing/python_src")
    // def python_unzip_script = python_dir.resolve("prepare_data.py")
    // def python_merge_script = python_dir.resolve("mgf_concatenate_I5_chunks.py")
    // def python_combine_script = python_dir.resolve("combine_tables.py")
    // def yaml_file = python_dir.resolve("ro-crate.yaml")

    ch_subfolders = Channel.of("220223_JARVIS_HWLTKDRXY", "220622_JARVIS_HMGW5DSX3", "220712_JARVIS_HMNJKDSX3")
    ch_archives_root = Channel.value(params.archives_root)


    mergeI5(python_merge_script, unzipArchive.out.archive_name)

    extractTables(unzipArchive.out.archive_name, ch_newArchive)

    deleteWorkFiles(extractTables.out.trigger, unzipArchive.out.archive_name)

    // println(params.folder_extracted_tables)
    // combineTables(extractTables.out.trigger, params.folder_extracted_tables, python_combine_script)
    
}