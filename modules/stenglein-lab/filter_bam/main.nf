process FILTER_BAM {
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
      path  "versions.yml",            emit: versions

    when:
      task.ext.when == null || task.ext.when

    script:
	   def args  = task.ext.args ?: ''
	   def args2 = task.ext.args2 ?: ''
		// index only in case need to do random retrieval (specifying regions in samtools view)
		def index_command = args2 ? "samtools index $input" : ''

          // --exclude-flags 0x4 \\

      """
      $index_command

      samtools \\
          view \\
          --threads ${task.cpus-1} \\
          $args \\
          --output-fmt bam \\
          -o ${meta.id}.filt.bam \\
          $input \\
          $args2


  
      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
      END_VERSIONS
    """
}
