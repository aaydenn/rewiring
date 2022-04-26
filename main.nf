#!/usr/bin/env nextflow

params.sdrf = 'data/intermittent_new.tsv'
params.annotation = 'data/miRNA-4_1-st-v1.annotations.20160922.csv.zip'
params.features = 'data/features.tsv'
params.outdir = 'results'

sdrf_file = file(params.sdrf)
anno_file = file(params.annotation)
feat_file = file(params.features)

process readAndNormalize {
    publishDir params.outdir, mode:'copy'

    output:
    file "*.rda" into data_ch
    file "*.png" into qc_ch

    script:
    """
    Rscript normalize.r -s $sdrf_file -a $anno_file
    """
}

process contrastModel {
    input:
    file norm_data from data_ch

    output:
    file "*.tsv" into de_table_ch
    file "*.rda" into de_ch

    script:
    """
    Rscript model.r -n $norm_data -f $feat_file
    """
}

process jointGraphicalModel {
    input:
    file norm_data from data_ch
    file de_list from de_ch

    output:
    file "*.rda" into fgl_table_ch
    file "*.png" into fgl_ch

    script:
    """
    Rscript jgl.r -n $de_list
    """
}
