#!/usr/bin/env nextflow

include { readYAML } from './modules/read_yaml.nf'
include { unzipArchive } from './modules/unzip_archive.nf'


// list of files
params.files = "inp_files.csv"
// Default parameter input
params.archives_root = "/usr/local/scratch/metaGOflow-COMPLETED-results/Batch1and2/CCMAR-data/FILTERS"


// Remove chunks from the I5 files
// Concatenate the I5 files:
// <prefix>.merged_CDS.I5_001.tsv.gz
// <prefix>.merged_CDS.I5_002.tsv.gz
// into:
// <prefix>.merged_CDS.I5.tsv.gz

// original structure, HMNJKDSX3:
// |   |-- functional-annotation
// |   |   |-- DBB.merged_CDS.I5_001.tsv.gz
// |   |   |-- DBB.merged_CDS.I5_002.tsv.gz
// |   |   |-- DBB.merged_CDS.I5.tsv.chunks

process mergeI5 {
    conda '/usr/local/scratch/nf-metaGOflow/wf-test/nf-testing/conda.yaml'

    input:
    path python_path
    path archive_name

    // output:
    // val identifier

    script:
    """
    python ${python_path} -d ${archive_name} 
    """
}


// I think this is the solution for the renaming problems too, or????
process extractTables {
    conda '/usr/local/scratch/nf-metaGOflow/wf-test/nf-testing/conda.yaml'
    publishDir "results-tables", mode: 'copy'

    input:
    path archive_name
    path target_directory
    
    output:
    path "${target_directory}-tables"
    val true, emit: trigger
    
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


// this does not work again because of the paths of inputs and outputs.
process combineTables {
    conda '/usr/local/scratch/nf-metaGOflow/wf-test/nf-testing/conda.yaml'
    publishDir "results-tables", mode: 'move'

    input:
    val ready
    path python_path

    output:
    path 'combined_tables/*'
    
    script:
    """
    python ${python_path}
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

    combineTables(extractTables.out.trigger, python_combine_script)
    
}