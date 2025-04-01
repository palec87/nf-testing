#!/usr/bin/env nextflow

// PROBLEMATIC FILES
// this one gives problesm, HCFCYDSX5.UDI137.tar.bz2: I still do not know why
// DBB_AAAJOSDA_4_HMGW5DSX3.UDI217.zip, HCMR data missing functional annotation
// HMGW5DSX3.UDI229.zip, HCMR data missing functional annotation
// DBB_AAAMOSDA_4_HMGW5DSX3.UDI241.zip, HCMR, missing chuncks data to sucessfully merge I5 files
// HMGW5DSX3.UDI253 no functional annotation
// HMGW5DSX3.UDI227 no functional annotation

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
    publishDir "results-tables", mode: 'copy'
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

// delete work files
process deleteWorkFiles {
    debug true

    input:
    val ready
    path archive_name

    script:
    """
    rm -rf ${archive_name}
    """
}


// this does not work again because of the paths of inputs and outputs.
process combineTables {
    conda '/usr/local/scratch/nf-metaGOflow/wf-test/nf-testing/conda.yaml'
    debug true

    input:
    val ready
    path folder_path
    path python_path


    
    script:
    """
    python ${python_path} ${folder_path}
    """
}

workflow {
    def python_dir = file("/usr/local/scratch/nf-metaGOflow/wf-test/nf-testing/python_src")
    def python_unzip_script = python_dir.resolve("prepare_data.py")
    def python_merge_script = python_dir.resolve("mgf_concatenate_I5_chunks.py")
    def python_combine_script = python_dir.resolve("combine_tables.py")
    def yaml_file = python_dir.resolve("ro-crate.yaml")

    ch_archives_root = Channel.value(params.archives_root)
    ch_file_path = Channel.fromPath(params.files)
                            .splitCsv()
                            .map { row -> file(row[0]) }
                            .view()

    unzipArchive(python_unzip_script, ch_archives_root, ch_file_path) // Unzip archives
    readYAML(python_dir, unzipArchive.out.archive_name, yaml_file) // Read YAML file

    ch_newArchive = readYAML.out.newArchiveName
        .splitCsv()
        .map { csv -> file(csv[0]) }
        .view()

    mergeI5(python_merge_script, unzipArchive.out.archive_name)

    extractTables(unzipArchive.out.archive_name, ch_newArchive)

    deleteWorkFiles(extractTables.out.trigger, unzipArchive.out.archive_name)

    // println(params.folder_extracted_tables)
    // combineTables(extractTables.out.trigger, params.folder_extracted_tables, python_combine_script)
    
}