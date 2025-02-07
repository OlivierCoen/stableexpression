//
// Subworkflow with functionality specific to the nf-core/stableexpression pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DESEQ2_NORMALISE                     } from '../../../modules/local/deseq2/normalise/main'
include { EDGER_NORMALISE                      } from '../../../modules/local/edger/normalise/main'
include { QUANTILE_NORMALISE                   } from '../../../modules/local/quantile_normalisation/main'
include { DATASET_STATISTICS                   } from '../../../modules/local/dataset_statistics/main'

/*
========================================================================================
    SUBWORKFLOW TO NORMALISE AND HARMONISE EXPRESSION DATASETS
========================================================================================
*/

workflow EXPRESSION_NORMALISATION {

    take:
    ch_datasets
    normalisation_method


    main:

    //
    // MODULE: normalisation of raw count datasets (including downloaded RNA-seq datasets)
    // at the same time, removing genes that show only zero counts
    //

    ch_datasets = ch_datasets.branch {
        meta, file ->
            raw: meta.normalised == false
            normalised: meta.normalised == true
        }

    ch_raw_rnaseq_datasets = ch_datasets.raw.filter { meta, file -> meta.platform == 'rnaseq' }

    if ( normalisation_method == 'deseq2' ) {
        DESEQ2_NORMALISE( ch_raw_rnaseq_datasets )
        ch_raw_rnaseq_datasets_normalised = DESEQ2_NORMALISE.out.cpm

    } else { // 'edger'
        EDGER_NORMALISE( ch_raw_rnaseq_datasets )
        ch_raw_rnaseq_datasets_normalised = EDGER_NORMALISE.out.cpm
    }

    //
    // MODULE: Quantile normalisation
    //

    // putting all normalised count datasets together and performing quantile normalisation
    ch_datasets.normalised.concat( ch_raw_rnaseq_datasets_normalised ) | QUANTILE_NORMALISE
    ch_quantile_normalised_datasets = QUANTILE_NORMALISE.out.counts

    //
    // MODULE: Dataset statistics
    //

    DATASET_STATISTICS( ch_quantile_normalised_datasets )

    emit:
    normalised_counts = ch_quantile_normalised_datasets
    dataset_statistics = DATASET_STATISTICS.out.stats

}




