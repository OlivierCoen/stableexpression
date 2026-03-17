#!/usr/bin/env nextflow

include { GENORM } from '../../../../subworkflows/local/genorm'


workflow {

    main:
    ch_counts = channel.fromPath(params.counts, checkIfExists: true)
                    .map{ file -> [ [dataset: file.name, section: file.name], file] }
    ch_counts.view()
    GENORM(ch_counts)

    publish:
    measures = GENORM.out.m_measures
}

output {
    measures {
    }
}
