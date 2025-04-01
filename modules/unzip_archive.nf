#!/usr/bin/env nextflow

// unzipArchive process
process unzipArchive {
    conda '/usr/local/scratch/nf-metaGOflow/wf-test/nf-testing/conda.yaml'
    
    input:
    path python_path
    path archives_root
    path file
    val hcmr

    output:
    path 'prepared_archives/*', emit: archive_name
    
    script:
    """
    python ${python_path} ${archives_root} ${file} -d ${hcmr ? '--hcmr' : ''}
    """
}

// python ${python_path} ${archives_root} ${file} -d --hcmr