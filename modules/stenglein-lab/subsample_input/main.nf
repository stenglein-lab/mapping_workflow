process subsample_input {
  label 'lowmem'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3"
  } else {                                                                      
      container "quay.io/biocontainers/seqtk:1.3--h5bf99c6_3"
  }   

  input:
  tuple val (meta), path(input_fastq)
  val(subsample_size)                  



  output:
  tuple val (meta), path("*.ss.fastq.gz) , emit: fastq
  path  "versions.yml"                   , emit: versions


      when:
    task.ext.when == null || task.ext.when


  script:
  def r1 = input_fastq[0]
  def r1_ss = r1.name.replaceAll("fastq", "ss.fastq")
  def r2 = input_fastq[1] 
  def r2_ss = input_fastq[1] ? r2.name.replaceAll("fastq", "ss.fastq") : ""
  def r1_command = "seqtk sample $r1 ${params.subsample_fraction} | gzip > $r1_ss" 
  def r2_command = input_fastq[1] ? "seqtk sample $r2 ${params.subsample_fraction} | gzip > $r2_ss" : ""
  """
  $r1_command
  $r2_command
  """
}

