process PROCESS_MAPPING_OUTPUT {
  label 'lowmem_non_threaded'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.3.2"
  }     

  input:
  path (script_to_run)
  path (input_tsv)
  path (contig_species_map)
  path (metadata)
  path (dataset_sizes_file)
  path (R_script_dir)

  output:
  path "*.txt"         , emit: txt, includeInputs: true

  when:
  task.ext.when == null || task.ext.when

  script:

  def args             = task.ext.args ?: ''

  """
   Rscript ${script_to_run} ${input_tsv} ${contig_species_map} ${metadata} ${dataset_sizes_file} ${R_script_dir}
  """

}
