process MERGE_COUNTS {
    // debug true
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
        'oras://community.wave.seqera.io/library/polars_python:1636dad2fac97d0b':
        'community.wave.seqera.io/library/polars_python:cab787b788e5eba7' }"

    input:
    path count_files, stageAs: "?/*"
    val(filter_out_zero_counts)
    val(get_means)

    output:
    path 'all_counts.parquet',                                                                                        emit: counts
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('polars'),   eval('python3 -c "import polars; print(polars.__version__)"'),     topic: versions


    script:
    def filter_out_zero_counts_arg = filter_out_zero_counts ? "--filter-out-zero-counts" : ""
    def get_means_arg = get_means ? "--get-means" : ""
    """
    merge_counts.py --counts "$count_files" $filter_out_zero_counts_arg $get_means_arg
    """

}
