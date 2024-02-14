// assumes input reads have been preprocessed (adapter/quality trimmed)
include { MARSHALL_FASTQ                             } from '../../subworkflows/stenglein-lab/marshall_fastq'
include { BUILD_GENOME_INDEX                         } from '../../subworkflows/stenglein-lab/build_genome_index'
include { MAP_TO_GENOME as MAP_TO_REFERENCE          } from '../../subworkflows/stenglein-lab/map_to_genome'
include { BAM_TO_SAM as BAM_TO_SAM_SPECIES_ID        } from '../../modules/stenglein-lab/bam_to_sam'
include { EXTRACT_SAM_INFO                           } from '../../modules/stenglein-lab/extract_sam_info'
include { TABULATE_MAPPING_STATS                     } from '../../modules/stenglein-lab/tabulate_mapping_stats'
include { PROCESS_SPECIES_OUTPUT                     } from '../../modules/stenglein-lab/process_species_output'

workflow SPECIES_ID {                                                    

  MARSHALL_FASTQ(params.fastq_dir, params.fastq_pattern)

  BUILD_GENOME_INDEX(params.species_id_fasta)

  MAP_TO_REFERENCE(MARSHALL_FASTQ.out.reads, BUILD_GENOME_INDEX.out.index)

  BAM_TO_SAM_SPECIES_ID(MAP_TO_REFERENCE.out.bam.filter{it[1].size() > 0})

  EXTRACT_SAM_INFO(BAM_TO_SAM_SPECIES_ID.out.sam)

  TABULATE_MAPPING_STATS(EXTRACT_SAM_INFO.out.txt.collectFile(name: "all_sam_stats.txt"){it[1]})

  PROCESS_SPECIES_OUTPUT(TABULATE_MAPPING_STATS.out.txt)
}

