nextflow.preview.topic = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXPRESSIONATLAS_FETCHDATA              } from '../subworkflows/local/expressionatlas_fetchdata/main'
include { PAIRWISE_GENE_VARIATION                } from '../subworkflows/local/pairwise_gene_variation/main.nf'
include { EXPRESSION_NORMALISATION               } from '../subworkflows/local/expression_normalisation/main.nf'

include { GPROFILER_IDMAPPING                    } from '../modules/local/gprofiler/idmapping/main'
include { MERGE_DATA                             } from '../modules/local/merge_data/main'
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
    ch_input_datasets


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

    // putting all datasets together (local datasets + Expression Atlas datasets)
    ch_datasets = ch_input_datasets.concat( EXPRESSIONATLAS_FETCHDATA.out.downloaded_datasets )

    //
    // MODULE: ID Mapping
    //

    ch_gene_metadata = Channel.empty()
    if ( params.gene_metadata ) {
        ch_gene_metadata = Channel.fromPath( params.gene_metadata, checkIfExists: true )
    }

    if ( params.skip_gprofiler ) {

        ch_gene_id_mapping = Channel.empty()
        if ( params.gene_id_mapping ) {
            // the gene id mappings will only be those provided by the user
            ch_gene_id_mapping = Channel.fromPath( params.gene_id_mapping, checkIfExists: true )
        }

    } else {
        // tries to map gene IDs to Ensembl IDs whenever possible
        GPROFILER_IDMAPPING(
            ch_datasets.combine( ch_species ),
            params.gene_id_mapping ? Channel.fromPath( params.gene_id_mapping, checkIfExists: true ) : 'none'
        )
        ch_datasets = GPROFILER_IDMAPPING.out.renamed
        ch_gene_metadata = ch_gene_metadata.mix( GPROFILER_IDMAPPING.out.metadata )
        // the gene id mappings are the sum
        // of those provided by the user and those fetched from g:Profiler
        ch_gene_id_mapping = GPROFILER_IDMAPPING.out.mapping
    }

    //
    // SURBWORKFLOW: normalisation of raw count datasets (including RNA-seq datasets)
    //

    EXPRESSION_NORMALISATION(
        ch_datasets,
        params.normalisation_method
    )

    //
    // MODULE: Merge count files and design files and filter out zero counts
    //

    EXPRESSION_NORMALISATION.out.normalised_counts
                                    .map { meta, file -> [file] }
                                    .collect()
                                    .set { ch_count_files }

    EXPRESSION_NORMALISATION.out.normalised_counts
                                    .map { meta, file -> [meta.design] }
                                    .collect()
                                    .set { ch_design_files }

    MERGE_DATA(
        ch_count_files,
        ch_design_files,
        params.nb_candidates_gene_variation
    )

    ch_merged_counts = MERGE_DATA.out.all_counts

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
        ch_gene_id_mapping.collect(),
        ch_m_measures
    )

    ch_top_stable_genes_summary = GENE_STATISTICS.out.top_stable_genes_summary
    ch_all_statistics = GENE_STATISTICS.out.all_statistics
    ch_log_counts = GENE_STATISTICS.out.log_counts
    ch_top_stable_genes_log_counts = GENE_STATISTICS.out.top_stable_genes_transposed_log_counts

    ch_multiqc_files = ch_multiqc_files
                        .mix( ch_top_stable_genes_summary.collect() )
                        .mix( ch_all_statistics.collect() )
                        .mix( ch_log_counts.collect() )
                        .mix( ch_top_stable_genes_log_counts.collect() )


    //
    // Collate and save software versions
    // TODO: use the nf-core functions when they are adapted to channel topics
    //

    ch_collated_versions = customSoftwareVersionsToYAML( Channel.topic('versions') )
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'stableexpression_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
        top_stable_genes_summary = ch_top_stable_genes_summary
        log_counts = ch_log_counts
        top_stable_genes_log_counts = ch_top_stable_genes_log_counts
        multiqc_report = MULTIQC.out.report.toList()

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
