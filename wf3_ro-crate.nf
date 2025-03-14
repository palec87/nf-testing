#!/usr/bin/env nextflow

// list of files
params.files = "inp_files.csv"
// Default parameter input
params.archives_root = "/usr/local/scratch/metaGOflow-COMPLETED-results/Batch1and2/CCMAR-data/FILTERS"
params.output_dir = './results'


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

process createRoCrate {
    conda '/usr/local/scratch/nf-metaGOflow/wf-test/nf-testing/conda.yaml'
    publishDir "results", mode: 'copy'
    // publishDir "results/${outFolder}", mode: 'copy'
    
    input:
    path archive_folder
    path python_path
    path yaml_file
    path outFolder

    output:
    path "${archive_folder}/*", emit: folder_path
    path 'metadata_part1.json', emit: metadata1
    
    script:
    """
    python ${python_path} ${archive_folder} ${yaml_file} -d
    """
}


process createRoCrate2 {
    conda '/usr/local/scratch/nf-metaGOflow/wf-test/nf-testing/conda.yaml'
    publishDir "results", mode: 'copy'
    // publishDir "results/${outFolder}", mode: 'copy'
    
    input:
    path archive_folder
    path python_path
    path yaml_file
    path outFolder

    output:
    path "${archive_folder}/*", emit: folder_path
    
    script:
    """
    python ${python_path} -d -u ${archive_folder} ${yaml_file}
    """
}

// comments on nomencalture
// create-ro-crate_partXX.py
// target_directory = folder where the files are, unzipArchive.out.archive_name
// Workflow block
workflow {
    def python_dir = file("/usr/local/scratch/nf-metaGOflow/wf-test/nf-testing/python_src")
    def python_unzip_script = python_dir.resolve("prepare_data.py")
    def python_ro_crate_script = python_dir.resolve("create-ro-crate_part1.py")
    // def python_ro_crate_script2 = python_dir.resolve("create-ro-crate_part2.py")
    def python_ro_crate_script3 = python_dir.resolve("create-ro-crate_part3.py")
    def yaml_file = python_dir.resolve("ro-crate.yaml")

    ch_archives_root = Channel.of(params.archives_root)
    ch_file_path = Channel.fromPath(params.files)
                            .splitCsv()
                            .map { csv -> file(csv[0]) }
                            .view { csv -> "After map: $csv" }
    // sp =  ch_file_path.split(".") // Split the file name and extension
    // ch_archiveName = "${sp[0]}.${sp[1]}" // Get the file name

    unzipArchive(python_unzip_script, ch_archives_root, ch_file_path) // Unzip archives
    readYAML(python_dir, unzipArchive.out.archive_name, yaml_file) // Read YAML file

    readYAML.out.newArchiveName.view { it -> "New Ro-Crate name: ${it}" }
    ch_newArchive = readYAML.out.newArchiveName
        .splitCsv()
        .map { csv -> file(csv[0]) }

    // ro-crate from the unzipped archive
    createRoCrate(
        unzipArchive.out.archive_name,
        python_ro_crate_script,
        yaml_file,
        ch_newArchive,
    )

    createRoCrate.out.folder_path.view { it -> "folder_path variable: ${it}" }
    createRoCrate.out.metadata1.view { it -> "metadata1 variable: ${it}" }

    println("projectDir: $projectDir")
    println("###### PART 2 of the rocrate creation ######")

    createRoCrate2(
        unzipArchive.out.archive_name,
        python_ro_crate_script3,
        yaml_file,
        ch_newArchive,
    )

}