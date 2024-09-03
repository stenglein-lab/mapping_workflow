#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { COMPETITIVE_MAPPING_WORKFLOW } from './subworkflows/stenglein-lab/competitive_mapping_workflow'

workflow {
    COMPETITIVE_MAPPING_WORKFLOW ()
}

