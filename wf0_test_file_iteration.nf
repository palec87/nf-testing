params.files = "inp.csv"
params.root_folder = "/media/davidp/Data/coding/nf-testing/nf-testing"

process iterateFiles {
    debug true
    publishDir "results/iterate", mode: 'move'
    
    input:
    // path folder
    path files
    
    output:
    path 'file_*'
    
    script:
    """
    cp ${files} file_${files}
    """
}
workflow {
    ch_files = Channel.fromPath(params.files)
                      .view { csv -> "Before splitCsv: $csv" }
                      .splitCsv()
                      .view { csv -> "After splitCsv: $csv" }
                      .map { csv -> file(csv[0]) }
                      .view { csv -> "After map: $csv" }

    iterateFiles(ch_files)
}