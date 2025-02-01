//
// Subworkflow with functionality specific to the nf-core/stableexpression pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MAKE_CHUNKS                  } from '../../../modules/local/pairwise_gene_variation/make_chunks/main'
include { CROSS_JOIN                   } from '../../../modules/local/pairwise_gene_variation/cross_join/main'
include { EXPRESSION_RATIO             } from '../../../modules/local/pairwise_gene_variation/expression_ratio/main'
include { RATIO_STANDARD_VARIATION     } from '../../../modules/local/pairwise_gene_variation/ratio_standard_variation/main'
include { COMPUTE_M_MEASURE            } from '../../../modules/local/pairwise_gene_variation/compute_m_measure/main'

/*
========================================================================================
    SUBWORKFLOW TO COMPUTE PAIRWISE GENE VARIATION
========================================================================================
*/

workflow EXPRESSION_HARMONISATION {

    take:
    ch_counts


    main:



    emit:


}




