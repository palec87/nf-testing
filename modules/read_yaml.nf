#!/usr/bin/env nextflow

process readYAML {
    debug true
    conda '/usr/local/scratch/nf-metaGOflow/wf-test/nf-testing/conda.yaml'
    publishDir "results", mode: 'copy'
    tag "$target_directory"

    input:
    path python_path
    path target_directory
    path yaml_file
    

    output:
    path "ro-crate-${target_directory}.csv", emit: newArchiveName
    
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

    with open(f"ro-crate-${target_directory}.csv", "w") as f:
        f.write(conf["source_mat_id"])
    """
}