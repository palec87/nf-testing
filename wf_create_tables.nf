#!/usr/bin/env nextflow

include { readYAML } from './modules/read_yaml.nf'
include { unzipArchive } from './modules/unzip_archive.nf'


// list of files
params.files = "inp_files.csv"
// Default parameter input
params.archives_root = "/usr/local/scratch/metaGOflow-COMPLETED-results/Batch1and2/CCMAR-data/FILTERS"



// I thinnk this is the solution for the renaming problems too, or????
process extractTables {
    conda '/usr/local/scratch/nf-metaGOflow/wf-test/nf-testing/conda.yaml'
    publishDir "results-tables", mode: 'copy'

    input:
    path archive_name
    path target_directory
    
    output:
    path "${target_directory}-tables"
    
    script:
    """
    mkdir -p ${target_directory}-tables
    cp ${archive_name}/results/functional-annotation/DBB.merged.summary.go ${target_directory}-tables
    cp ${archive_name}/results/functional-annotation/DBB.merged.summary.go_slim ${target_directory}-tables
    cp ${archive_name}/results/functional-annotation/DBB.merged.summary.ips ${target_directory}-tables
    cp ${archive_name}/results/functional-annotation/DBB.merged.summary.ko ${target_directory}-tables
    cp ${archive_name}/results/functional-annotation/DBB.merged.summary.pfam ${target_directory}-tables
    cp ${archive_name}/results/taxonomy-summary/LSU/DBB.merged_LSU.fasta.mseq.tsv ${target_directory}-tables
    cp ${archive_name}/results/taxonomy-summary/SSU/DBB.merged_SSU.fasta.mseq.tsv ${target_directory}-tables
    """
}

workflow {
    def python_dir = file("/usr/local/scratch/nf-metaGOflow/wf-test/nf-testing/python_src")
    def python_unzip_script = python_dir.resolve("prepare_data.py")
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

    extractTables(unzipArchive.out.archive_name, ch_newArchive)
    
}