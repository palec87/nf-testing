#!/usr/bin/env nextflow

// list of files
params.files = "inp_files.csv"
// Default parameter input
params.archives_root = "/usr/local/scratch/metaGOflow-COMPLETED-results/Batch1and2/CCMAR-data/FILTERS"


// unzipArchive process
process unzipArchive {
    publishDir "results/unzip", mode: 'move'
    
    input:
    path python_path
    path archives_root
    path files

    output:
    path 'prepared_archives/*'
    
    script:
    """
    python ${python_path} ${archives_root} ${files} -d
    """
}

// Workflow block
workflow {
    def python_dir = file("/usr/local/scratch/nf-metaGOflow/wf-test/nf-testing/python_src")
    def python_script = python_dir.resolve("prepare_data.py")
    ch_archives_root = Channel.of(params.archives_root)
    ch_file_path = Channel.fromPath(params.files)
                            .splitCsv()
                            .map { csv -> file(csv[0]) }
                            .view { csv -> "After map: $csv" }

    unzipArchive(python_script, ch_archives_root, ch_file_path) // Unzip archives
}


