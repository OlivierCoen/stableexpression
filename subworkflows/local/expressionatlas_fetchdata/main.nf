//
// Subworkflow with functionality specific to the nf-core/stableexpression pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXPRESSIONATLAS_GETACCESSIONS          } from '../../../modules/local/expressionatlas/getaccessions/main'
include { EXPRESSIONATLAS_GETDATA                } from '../../../modules/local/expressionatlas/getdata/main'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
========================================================================================
*/

workflow EXPRESSIONATLAS_FETCHDATA {

    take:
    ch_species
    eatlas_accessions
    eatlas_keywords
    fetch_eatlas_accessions


    main:

    ch_accessions = Channel.fromList( eatlas_accessions.tokenize(',') )

    // fetching Expression Atlas accessions if applicable
    if ( fetch_eatlas_accessions || eatlas_keywords ) {

        //
        // MODULE: Expression Atlas - Get accessions
        //
        ch_eatlas_keywords = Channel.value( eatlas_keywords )

        // getting Expression Atlas accessions given a species name and keywords
        // keywords can be an empty string
        EXPRESSIONATLAS_GETACCESSIONS( ch_species, ch_eatlas_keywords )

        // appending to accessions provided by the user
        // ensures that no accessions is present twice (provided by the user and fetched from E. Atlas)
        ch_accessions = ch_accessions
                            .concat( EXPRESSIONATLAS_GETACCESSIONS.out.txt.splitText() )
                            .unique()
                            .map { it -> it.trim() }
                            .filter { it.startsWith('E-') && !it.startsWith('E-PROT-') }
    }

    //
    // MODULE: Expression Atlas - Get data
    //

    // Downloading Expression Atlas data for each accession in ch_accessions
    EXPRESSIONATLAS_GETDATA( ch_accessions )

    // adding dataset id (accession + data_type) in the file meta
    ch_etlas_output_design = augmentWithDatasetId( EXPRESSIONATLAS_GETDATA.out.design.flatten() )
    ch_eatlas_microarray_normalised_output = augmentWithDatasetId( EXPRESSIONATLAS_GETDATA.out.microarray_normalised.flatten() )
    ch_eatlas_rnaseq_raw_output = augmentWithDatasetId( EXPRESSIONATLAS_GETDATA.out.rnaseq_raw.flatten() )

    // adding design files to the meta their respective count files
    ch_microarray_normalised = groupFilesByDatasetId(
        ch_etlas_output_design,
        ch_eatlas_microarray_normalised_output
    )

    ch_rnaseq_raw = groupFilesByDatasetId(
        ch_etlas_output_design,
        ch_eatlas_rnaseq_raw_output
    )

    emit:
    microarray_normalised = ch_microarray_normalised
    rnaseq_raw = ch_rnaseq_raw
}



/*
========================================================================================
    FUNCTIONS
========================================================================================
*/


//
// Get Expression Atlas Batch ID (accession + data_type) from file stem
//
def augmentWithDatasetId( ch_files ) {
    return ch_files
            .map {
                file ->
                    def meta = [dataset: file.getSimpleName()]
                    [meta, file]
            }
}

//
// Groups design and data files by accession and data_type
//
def groupFilesByDatasetId(ch_design, ch_counts) {
    return ch_design
        .concat( ch_counts ) // puts counts at the end of the resulting channel
        .groupTuple() // groups by dataset ID; design files are necessarily BEFORE count files
        .filter {
            it.get(1).size() == 2 // only groups with two files
        }
        .filter { // only groups with first file as design file and second one as count file
            meta, files ->
                files.get(0).name.endsWith('.design.csv') && !files.get(1).name.endsWith('.design.csv')
        }
        .map { // putting design file in meta
            meta, files ->
                def new_meta = meta + [design: files[0]]
                [new_meta, files[1]]
        }
}
