// assumes input reads have been preprocessed (adapter/quality trimmed)
include { MARSHALL_FASTQ                             } from '../../subworkflows/stenglein-lab/marshall_fastq'
include { SAVE_OUTPUT_FILE as SAVE_COUNTS_FILE       } from '../../modules/stenglein-lab/save_output_file'
include { BUILD_GENOME_INDEX                         } from '../../subworkflows/stenglein-lab/build_genome_index'
include { MAP_TO_GENOME                              } from '../../subworkflows/stenglein-lab/map_to_genome'
include { MAPPING_STATS                              } from '../../subworkflows/stenglein-lab/mapping_stats'
include { PREPEND_TSV_WITH_ID                        } from '../../modules/stenglein-lab/prepend_tsv_with_id'
include { SAVE_OUTPUT_FILE                           } from '../../modules/stenglein-lab/save_output_file'
include { PROCESS_MAPPING_OUTPUT                     } from '../../modules/stenglein-lab/process_mapping_output'

// metadata
include { COLLECT_METADATA                               } from '../../modules/stenglein-lab/collect_metadata'
include { SAVE_OUTPUT_FILE as SAVE_METADATA_FILE         } from '../../modules/stenglein-lab/save_output_file'


workflow COMPETITIVE_MAPPING_WORKFLOW {                                                    

  def count_fastq = true

  // make sure FASTQ in order and count # of reads 
  MARSHALL_FASTQ(params.fastq_dir, params.fastq_pattern, count_fastq, params.subsample_size)
  SAVE_COUNTS_FILE(MARSHALL_FASTQ.out.fastq_counts.collectFile(name: "fastq_counts.txt"))

  // a directory with additional R scripts
  R_script_dir_ch = Channel.fromPath(params.R_shared_script_dir)

  // directory with metadata files
  metadata_dir_ch = Channel.fromPath(params.metadata_dir)
  COLLECT_METADATA(R_script_dir_ch, metadata_dir_ch)

  // save to main results dir
  SAVE_METADATA_FILE(COLLECT_METADATA.out.collected_metadata)

  // make the genome index
  BUILD_GENOME_INDEX(params.genome_fasta)

  // do mapping
  MAP_TO_GENOME(MARSHALL_FASTQ.out.reads, BUILD_GENOME_INDEX.out.index)

  // tabulate mapping stats
  def per_base_coverage = !params.skip_per_base_coverage
  MAPPING_STATS(MAP_TO_GENOME.out.bam, BUILD_GENOME_INDEX.out.fasta, per_base_coverage)

  // prepend samtools output with sample IDs
  PREPEND_TSV_WITH_ID(MAPPING_STATS.out.coverage)
  SAVE_OUTPUT_FILE(PREPEND_TSV_WITH_ID.out.tsv.collectFile(name: "collected_coverage.txt"){it[1]})

  // call R script to analyze data

  contig_species_map_file_ch = Channel.fromPath(params.contig_species_map_file)
  output_processing_script_ch = Channel.fromPath(params.output_processing_R_script)

  PROCESS_MAPPING_OUTPUT(output_processing_script_ch,
                         SAVE_OUTPUT_FILE.out.file, 
                         contig_species_map_file_ch,
                         COLLECT_METADATA.out.rds,
                         SAVE_COUNTS_FILE.out.file,
                         R_script_dir_ch)

}
