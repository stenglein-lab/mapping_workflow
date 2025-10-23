#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { MAPPING_WORKFLOW } from './subworkflows/stenglein-lab/mapping_workflow'

workflow {
    MAPPING_WORKFLOW ()
}

