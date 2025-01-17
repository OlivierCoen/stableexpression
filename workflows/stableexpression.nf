nextflow.preview.topic = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXPRESSIONATLAS_FETCHDATA              } from '../subworkflows/local/expressionatlas_fetchdata/main'
include { PAIRWISE_GENE_VARIATION                } from '../subworkflows/local/pairwise_gene_variation/main.nf'

include { DESEQ2_NORMALISE                       } from '../modules/local/deseq2/normalise/main'
include { EDGER_NORMALISE                        } from '../modules/local/edger/normalise/main'
include { GPROFILER_IDMAPPING                    } from '../modules/local/gprofiler/idmapping/main'
include { MERGE_COUNTS                           } from '../modules/local/merge_counts/main'
include { GENE_STATISTICS                        } from '../modules/local/gene_statistics/main'
include { MULTIQC                                } from '../modules/nf-core/multiqc/main'

include { parseInputDatasets                     } from '../subworkflows/local/utils_nfcore_stableexpression_pipeline'
include { customSoftwareVersionsToYAML           } from '../subworkflows/local/utils_nfcore_stableexpression_pipeline'
include { validateInputParameters                } from '../subworkflows/local/utils_nfcore_stableexpression_pipeline'
include { methodsDescriptionText                 } from '../subworkflows/local/utils_nfcore_stableexpression_pipeline'
include { paramsSummaryMultiqc                   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { paramsSummaryMap                       } from 'plugin/nf-schema'




/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow STABLEEXPRESSION {

    take:
    ch_raw_datasets
    ch_normalised_datasets


    main:

    ch_multiqc_files = Channel.empty()

    ch_species = Channel.value( params.species.split(' ').join('_') )

    //
    // SUBWORKFLOW: fetching Expression Atlas datasets if needed
    //

    EXPRESSIONATLAS_FETCHDATA(
        ch_species,
        params.eatlas_accessions,
        params.eatlas_keywords,
        params.fetch_eatlas_accessions
    )

    // putting all normalized and raw datasets together (local datasets + Expression Atlas datasets)
    ch_normalised_datasets = ch_normalised_datasets.concat( EXPRESSIONATLAS_FETCHDATA.out.normalised )
    ch_raw_datasets = ch_raw_datasets.concat( EXPRESSIONATLAS_FETCHDATA.out.raw )

    //
    // MODULE: normalisation of raw count datasets (including RNA-seq datasets)
    //

    if ( params.normalisation_method == 'deseq2' ) {
        DESEQ2_NORMALISE( ch_raw_datasets )
        ch_raw_datasets_normalised = DESEQ2_NORMALISE.out.cpm

    } else { // 'edger'
        EDGER_NORMALISE( ch_raw_datasets )
        ch_raw_datasets_normalised = EDGER_NORMALISE.out.cpm
    }

    // putting all normalised count datasets together
    ch_normalised_datasets.concat( ch_raw_datasets_normalised ).set{ ch_all_normalised }


    //
    // MODULE: ID Mapping
    //

    // tries to map gene IDs to Ensembl IDs whenever possible
    GPROFILER_IDMAPPING( ch_all_normalised.combine(ch_species) )
    ch_counts_gene_renamed = GPROFILER_IDMAPPING.out.renamed
    ch_gene_metadata = GPROFILER_IDMAPPING.out.metadata
    ch_renaming_mapping = GPROFILER_IDMAPPING.out.mapping

    //
    // MODULE: Merge count files
    //

    MERGE_COUNTS( ch_counts_gene_renamed.collect() )
    ch_merged_counts = MERGE_COUNTS.out.counts

    //
    // STEP: Calculate gene variation
    //

    if (params.skip_gene_variation_calc) {

        ch_m_measures = Channel.of( 'none' )

    } else {

        if ( params.gene_variation_method == 'pairwise_gene_variation' ) {
            PAIRWISE_GENE_VARIATION ( ch_merged_counts )
            ch_m_measures = PAIRWISE_GENE_VARIATION.out.m_measures
        }

    }


    //
    // MODULE: Gene statistics
    //
    GENE_STATISTICS(
        ch_merged_counts,
        ch_gene_metadata.collect(),
        ch_renaming_mapping.collect(),
        ch_m_measures
    )

    ch_stats_most_stable_genes = GENE_STATISTICS.out.stats_most_stable_genes
    ch_stats_all_genes = GENE_STATISTICS.out.stats_all_genes
    ch_log_count_summary = GENE_STATISTICS.out.log_count_summary

    ch_multiqc_files = ch_multiqc_files
                        .mix( ch_stats_most_stable_genes.collect() )
                        .mix( ch_stats_all_genes.collect() )
                        .mix( ch_log_count_summary.collect() )


    //
    // Collate and save software versions
    // TODO: use the nf-core functions when they are adapted to channel topics
    //

    ch_collated_versions = customSoftwareVersionsToYAML( Channel.topic('versions') )
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_stableexpression_software_mqc_versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //

    // Load MultiQC configuration file
    ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.fromPath("${workflow.projectDir}/docs/images/nf-core-stableexpression_logo_light.png", checkIfExists: true)

    // Prepare the workflow summary
    ch_workflow_summary = Channel.value(
        paramsSummaryMultiqc(
            paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        )
    ).collectFile(name: 'workflow_summary_mqc.yaml')

    // Prepare the methods section
    ch_methods_description = Channel.value(
        methodsDescriptionText(
            params.multiqc_methods_description
                ? file(params.multiqc_methods_description)
                : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        )
    ).collectFile(name: 'methods_description_mqc.yaml')

    // Add summary, versions, and methods to the MultiQC input file list
    ch_multiqc_files = ch_multiqc_files
        .mix(ch_workflow_summary)
        .mix(ch_collated_versions)
        .mix(ch_methods_description)

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )
    ch_multiqc_report = MULTIQC.out.report


    emit:
        stats_most_stable_genes = ch_stats_most_stable_genes
        stats_all_genes = ch_stats_all_genes
        count_summary = ch_log_count_summary
        multiqc_report = ch_multiqc_report


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
