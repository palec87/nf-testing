#!/usr/bin/env nextflow

// Default parameter input
params.archives_root = "/usr/local/scratch/metaGOflow-COMPLETED-results/Batch1and2/CCMAR-data/FILTERS"
// params.archive_list = ("HMNJKDSX3.UDI244.tar.bz2")
params.file_path = "HMNJKDSX3.UDI244.tar.bz2"
// unzipArchive process

process unzipArchive {
    publishDir "results/unzip"
    
    input:
    path python_path
    path archives_root
    // path archive_list
    val file_path

    output:
    path 'prepared_archives/*'
    
    script:
    // """
    // for archive in ${archive_list}
    // do
    //     python ${python_path} ${archives_root} \${archive} -d
    // done
    // """
    """
    python ${python_path} ${archives_root} ${file_path} -d
    """
}

// Workflow block
workflow {
    def python_dir = file("/usr/local/scratch/nf-metaGOflow/wf-test/nf-testing/python_src")
    def python_script = python_dir.resolve("prepare_data.py")
    ch_archives_root = Channel.of(params.archives_root)
    // ch_archives = Channel.of(params.archive_list).toList().view() // Create a channel using parameter input
    ch_file_path = Channel.of(params.file_path)
    unzipArchive(python_script, ch_archives_root, ch_file_path) // Unzip archives
}


