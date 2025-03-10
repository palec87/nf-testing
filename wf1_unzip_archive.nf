#!/usr/bin/env nextflow

// Default parameter input
params.archives_root = "/usr/local/scratch/metaGOflow-COMPLETED-results/Batch1and2/CCMAR-data/FILTERS"
params.archive_list = ["HMNJKDSX3.UDI244.tar.bz2"]
// unzipArchive process

process unzipArchive {
    publishDir "results/unzip"
    
    input:
    path python_path
    path archives_root
    val archive

    output:
    path 'unzip_*'
    
    script:
    """
    python ${python_path}/prepare_data.py ${archives_root} ${archive} -d
    """
}

// Workflow block
workflow {
    def python_dir = file("/usr/local/scratch/nf-metaGOflow/wf-test/nf-testing/python_src")
    def python_script = python_dir.resolve("prepare_data.py")
    def archives_root = file(params.archives_root)
    ch_archives = Channel.of(params.archive_list) // Create a channel using parameter input
    unzipArchive(python_script, archives_root, ch_archives) // Unzip archives
}


