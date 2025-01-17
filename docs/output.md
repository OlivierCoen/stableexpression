# nf-core/stableexpression: Output

## Introduction

This document describes the output produced by the pipeline.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Expression Atlas](#expression-atlas): get Expression Atlas accessions and download data
- [Normalisation](#normalisation): normalise raw data (with DESeq2 or EdgeR)
- [gProfiler](#gprofiler-idmapping): map gene IDS to Ensembl IDS
- [Gene Statistics](#gene-statistics): merge all counts, compute gene variation statistics and get the most stable genes
- [MultiQC](#multiqc): generate reports

## Output files

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - MultiQC report file: `multiqc_report.html`.
  - MultiQC data dir: `multiqc_data`.
  - Plots created by MultiQC: `multiqc_plots`.

</details>

### Gene Variation

<details markdown="1">
<summary>Output files</summary>

- `gene_variation/`
  - A list of the most stable genes in `stats_most_stable_genes.csv`.
  - Descriptive statistics for all genes in `stats_all_genes.csv`
  - All normalised counts (for each gene and each sample) in `count_summary.csv`.

</details>

### Expression Atlas

<details markdown="1">
<summary>Output files</summary>

- `expressionatlas/`
  - List of accessions found when querying Expression Atlas: `accessions.txt`.
  - List of count datasets (normalized: `*.normalised.csv` / raw: `*.raw.csv`) and experimental designs (`*.design.csv`) downloaded from Expression Atlas.

</details>

### Normalisation

<details markdown="1">
<summary>Output files</summary>

List of newly normalised datasets in `normalisation/`

- `normalisation/deseq2/` for DESeq2
- `normalisation/edger/` for EdgeR

</details>

### GProfiler IDMapping

<details markdown="1">
<summary>Output files</summary>

- `idmapping/`
  - Count datasets whose gene IDs have been mapped to Ensembl IDs: `*_renamed.csv`.
  - Table associating original gene IDs and Ensembl IDs: `*_mapping.csv`.
  - Ensembl gene metadata (name and description): `*_metadata.csv`.

</details>

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
