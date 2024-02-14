/*
  Filter bam files with no mapping reads, which may have file size > 0 because of header lines
 */
process FILTER_EMPTY_BAM {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
    tuple val(meta), path (input)

    output:
    tuple val(meta), path("*.filt.bam"),  emit: bam

    when:
    task.ext.when == null || task.ext.when

    shell:
    '''
	 # use samtools stats to count the # of mapped reads
    mapping_read_count=$(samtools stats !{input} | grep "^SN" | grep "reads mapped:" | cut -f 3)

    # only create a link to bam if there are >0 mapping reads in bam file
    if [[ $mapping_read_count -gt 0 ]]
    then
      ln !{input} !{meta.id}.filt.bam
    else
      touch !{meta.id}.filt.bam
    fi

    '''
}
