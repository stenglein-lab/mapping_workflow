// assumes input reads have been preprocessed (adapter/quality trimmed)
include { MARSHALL_FASTQ                             } from '../../subworkflows/stenglein-lab/marshall_fastq'
include { BUILD_GENOME_INDEX                         } from '../../subworkflows/stenglein-lab/build_genome_index'
include { MAP_TO_GENOME                              } from '../../subworkflows/stenglein-lab/map_to_genome'
include { MAPPING_STATS                              } from '../../subworkflows/stenglein-lab/mapping_stats'
// include { TABULATE_MAPPING_STATS                     } from '../../modules/stenglein-lab/tabulate_mapping_stats'
// include { PROCESS_MAPPING_STATS                      } from '../../modules/stenglein-lab/process_species_output'

workflow MAPPING_WORKFLOW {                                                    

  def count_fastq = false
  MARSHALL_FASTQ(params.fastq_dir, params.fastq_pattern, count_fastq, params.subsample_size)

  BUILD_GENOME_INDEX(params.genome_fasta)

  MAP_TO_GENOME(MARSHALL_FASTQ.out.reads, BUILD_GENOME_INDEX.out.index)

  def per_base_coverage = false
  MAPPING_STATS(MAP_TO_GENOME.out.bam, BUILD_GENOME_INDEX.out.fasta, per_base_coverage)

}
