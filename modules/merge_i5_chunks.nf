#!/usr/bin/env nextflow

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
    debug true

    input:
    path python_path
    path archive_name

    script:
    """
    python ${python_path} -d ${archive_name} 
    """
}