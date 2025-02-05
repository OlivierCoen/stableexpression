process MERGE_DATA {

    label 'process_low'

    errorStrategy = {
        if (task.exitStatus == 100) {
            log.error(
                "No count could be found before merging datasets! "
                + "Please check the provided accessions and datasets and run again"
                )
            return 'terminate'
        }
        if (task.exitStatus == 101) {
            log.error(
                "When filtering for genes with at least one count in all datasets, no gene was found! "
                + "Please check the provided accessions and datasets and run again."
                )
            return 'terminate'
        }
    }

    publishDir "${params.outdir}/merged_counts"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f8a5d02e7b31980c887253a9f118da0ef91ead1c7b158caf855199e5c5d5473/data':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    path count_files, stageAs: "?/*"
    path design_files, stageAs: "?/*"
    val nb_candidate_genes

    output:
    path 'all_counts.parquet',                                                                                        emit: all_counts
    path 'candidate_gene_counts.parquet',                                                                             emit: candidate_gene_counts
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions


    script:
    """
    merge_data.py \
        --counts "$count_files" \
        --designs "$design_files" \
        --nb-candidate-genes "$nb_candidate_genes"
    """

}
