// assumes input reads have been preprocessed (adapter/quality trimmed)
include { MARSHAL_FASTQ                              } from '../../subworkflows/stenglein-lab/marshal_fastq'
include { SAVE_OUTPUT_FILE as SAVE_COUNTS_FILE       } from '../../modules/stenglein-lab/save_output_file'
include { BUILD_GENOME_INDEX                         } from '../../subworkflows/stenglein-lab/build_genome_index'
include { MAP_TO_GENOME                              } from '../../subworkflows/stenglein-lab/map_to_genome'
include { MAPPING_STATS                              } from '../../subworkflows/stenglein-lab/mapping_stats'
include { PREPEND_TSV_WITH_ID                        } from '../../modules/stenglein-lab/prepend_tsv_with_id'

include { SAVE_OUTPUT_FILE                             } from '../../modules/stenglein-lab/save_output_file'
include { SAVE_OUTPUT_FILE as SAVE_COLLECTED_COVERAGE  } from '../../modules/stenglein-lab/save_output_file'
include { SAVE_OUTPUT_FILE as SAVE_COLLECTED_STATS     } from '../../modules/stenglein-lab/save_output_file'
include { SAVE_OUTPUT_FILE as SAVE_COLLECTED_DEPTH     } from '../../modules/stenglein-lab/save_output_file'


workflow MAPPING_WORKFLOW {                                                    

  def count_fastq = true

  // make sure FASTQ in order and count # of reads 
  MARSHAL_FASTQ(params.fastq_dir, params.fastq_pattern, count_fastq, params.subsample_size)
  SAVE_COUNTS_FILE(MARSHAL_FASTQ.out.fastq_counts.collectFile(name: "initial_fastq_counts.txt"))

  // make the genome index
  BUILD_GENOME_INDEX(params.genome_fasta)

  // do mapping
  MAP_TO_GENOME(MARSHAL_FASTQ.out.reads, BUILD_GENOME_INDEX.out.index)

  // tabulate mapping stats
  def per_base_coverage = !params.skip_per_base_coverage
  MAPPING_STATS(MAP_TO_GENOME.out.bam, BUILD_GENOME_INDEX.out.fasta, per_base_coverage)

  // save consolidated output files from mapping stats
  SAVE_COLLECTED_COVERAGE(MAPPING_STATS.out.prepended_coverage.collectFile(name: "collected_per_refseq_coverage.tsv"){it[1]})
  SAVE_COLLECTED_STATS   (MAPPING_STATS.out.prepended_stats.collectFile(name: "collected_stats.tsv"){it[1]})
  SAVE_COLLECTED_DEPTH   (MAPPING_STATS.out.prepended_depth.collectFile(name: "collected_per_base_depth.tsv"){it[1]})

  // prepend samtools output with sample IDs
  PREPEND_TSV_WITH_ID(MAPPING_STATS.out.coverage)
  SAVE_OUTPUT_FILE(PREPEND_TSV_WITH_ID.out.tsv.collectFile(name: "collected_coverage.txt"){it[1]})

}
