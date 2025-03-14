#!/usr/bin/env nextflow

// list of files
params.files = "inp_files.csv"
// Default parameter input
params.archives_root = "/usr/local/scratch/metaGOflow-COMPLETED-results/Batch1and2/CCMAR-data/FILTERS"


// unzipArchive process
process unzipArchive {
    conda '/usr/local/scratch/nf-metaGOflow/wf-test/nf-testing/conda.yaml'
    
    input:
    path python_path
    path archives_root
    path files

    output:
    path 'prepared_archives/*', emit: archive_name
    
    script:
    """
    python ${python_path} ${archives_root} ${files} -d
    """
}

process readYAML {
    conda '/usr/local/scratch/nf-metaGOflow/wf-test/nf-testing/conda.yaml'
    publishDir "results", mode: 'copy'

    input:
    path python_path
    path target_directory
    path yaml_file
    

    output:
    path 'ro-crate-name.csv', emit: newArchiveName
    
    script:
    """
    #!/usr/bin/env python
    import sys
    import os
    from pathlib import Path


    sys.path.append("${python_path}")
    import read_config as rc

    conf = rc.read_yaml("${yaml_file}")
    print(conf)
    run_id = Path("${target_directory}").name

    conf["run_id"] = run_id
    conf = rc.get_ref_code_and_prefix(conf)

    with open("ro-crate-name.csv", "w") as f:
        f.write(conf["source_mat_id"])
    """
}

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
    cp ${archive_name}/results/functional_annotation/DBB.merged.summary.go ${target_directory}-tables
    cp ${archive_name}/results/functional_annotation/DBB.merged.summary.go_slim ${target_directory}-tables
    cp ${archive_name}/results/functional_annotation/DBB.merged.summary.ips ${target_directory}-tables
    cp ${archive_name}/results/functional_annotation/DBB.merged.summary.ko ${target_directory}-tables
    cp ${archive_name}/results/functional_annotation/DBB.merged.summary.pfam ${target_directory}-tables
    cp ${archive_name}/results/taxonomy_summary/LSU/DBB.merged_LSU.fasta.mseq.tsv ${target_directory}-tables
    cp ${archive_name}/results/taxonomy_summary/SSU/DBB.merged_SSU.fasta.mseq.tsv ${target_directory}-tables
    """
}

workflow {
    def python_dir = file("/usr/local/scratch/nf-metaGOflow/wf-test/nf-testing/python_src")
    def python_unzip_script = python_dir.resolve("prepare_data.py")
    // def python_ro_crate_script = python_dir.resolve("create-ro-crate_part1.py")
    // def python_ro_crate_script2 = python_dir.resolve("create-ro-crate_part2.py")
    // def python_ro_crate_script3 = python_dir.resolve("create-ro-crate_part3.py")
    def yaml_file = python_dir.resolve("ro-crate.yaml")

    ch_archives_root = Channel.of(params.archives_root)
    ch_file_path = Channel.fromPath(params.files)
                            .splitCsv()
                            .map { csv -> file(csv[0]) }
                            .view { csv -> "After map: $csv" }

    unzipArchive(python_unzip_script, ch_archives_root, ch_file_path) // Unzip archives
    readYAML(python_dir, unzipArchive.out.archive_name, yaml_file) // Read YAML file

    ch_newArchive = readYAML.out.newArchiveName
        .splitCsv()
        .map { csv -> file(csv[0]) }

    extractTables(unzipArchive.out.archive_name, ch_newArchive)
    
}