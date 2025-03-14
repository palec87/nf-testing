// list of files
params.files = "inp_files.csv"
// Default parameter input
params.archives_root = "/usr/local/scratch/metaGOflow-COMPLETED-results/Batch1and2/CCMAR-data/FILTERS"

process printInput {
    debug true

    input:
    path folder
    path file_path

    script:
    """
    echo ${file_path}
    """
}

workflow {
ch_archives_root = Channel.value(params.archives_root)
ch_file_path = Channel.fromPath(params.files)
                        .splitCsv()
                        .map { row -> file(row[0]) }
                        .view()
printInput(ch_archives_root, ch_file_path)
}
