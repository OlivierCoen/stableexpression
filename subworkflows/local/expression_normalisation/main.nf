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
    //

    ch_datasets.branch {
        meta, file ->
            raw: meta.normalised == false
            normalised: meta.normalised == true
        }
        .set { ch_datasets }

    if ( normalisation_method == 'deseq2' ) {
        DESEQ2_NORMALISE( ch_datasets.raw )
        ch_raw_datasets_normalised = DESEQ2_NORMALISE.out.cpm

    } else { // 'edger'
        EDGER_NORMALISE( ch_datasets.raw )
        ch_raw_datasets_normalised = EDGER_NORMALISE.out.cpm
    }

    // updating "normalised" state in metadata (only used for introspection for now)
    ch_raw_datasets_normalised
        .map { meta, file ->
            meta.normalised = true
            [meta, file]
        }
        .set { ch_raw_datasets_normalised }

    //
    // MODULE: Quantile normalisation
    //

    // putting all normalised count datasets together and performing quantile normalisation
    ch_datasets.normalised.concat( ch_raw_datasets_normalised ) | QUANTILE_NORMALISE

    emit:
    normalised_counts = QUANTILE_NORMALISE.out.counts

}




