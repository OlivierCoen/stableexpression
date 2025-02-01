//
// Subworkflow with functionality specific to the nf-core/stableexpression pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DESEQ2_NORMALISE                  } from '../../../modules/local/deseq2/normalise/main'
include { EDGER_NORMALISE                   } from '../../../modules/local/edger/normalise/main'
/*
========================================================================================
    SUBWORKFLOW TO NORMALISE AND HARMONISE EXPRESSION DATASETS
========================================================================================
*/

workflow EXPRESSION_NORMALISATION {

    take:
    ch_datasets


    main:

    //
    // MODULE: normalisation of raw count datasets (including downloaded RNA-seq datasets)
    //

    ch_raw_datasets = ch_datasets.filter { it.normalised == false }

    if ( params.normalisation_method == 'deseq2' ) {
        DESEQ2_NORMALISE( ch_raw_datasets )
        ch_raw_datasets_normalised = DESEQ2_NORMALISE.out.cpm

    } else { // 'edger'
        EDGER_NORMALISE( ch_raw_datasets )
        ch_raw_datasets_normalised = EDGER_NORMALISE.out.cpm
    }

    // putting all normalised count datasets together
    ch_normalised_datasets = ch_normalised_datasets.concat( ch_raw_datasets_normalised )

    // getting count means for all genes
    MERGE_COUNTS(
        ch_normalised_counts,
        false, // DO NOT filter out genes with zero counts
        true // return counts means
        )

    //
    // MODULE: Merge count files
    //
    ch_normalised_datasets
                    .map { meta, file -> [file] }
                    .collect()
                    .set { ch_normalised_counts}

    // merging counts and filtering out zero counts
    MERGE_COUNTS(
        ch_normalised_counts,
        true, // filter out genes with zero counts
        false // return counts
        )

    emit:
    merged_counts = MERGE_COUNTS.out.counts

}




