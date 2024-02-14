include { BOWTIE2_ALIGN       } from '../../modules/nf-core/bowtie2/align/main'
include { FILTER_BAM          } from '../../modules/stenglein-lab/filter_bam/main'

workflow MAP_TO_GENOME {

 take:
  reads
  index

 main:

  ch_versions     = Channel.empty()                                               

  // map reads 
  def sort_bam = true 
  def save_unaligned = false 
  BOWTIE2_ALIGN (reads, index, save_unaligned, sort_bam)
  ch_versions = ch_versions.mix ( BOWTIE2_ALIGN.out.versions )      

  FILTER_BAM(BOWTIE2_ALIGN.out.bam)
  ch_versions = ch_versions.mix ( FILTER_BAM.out.versions )      

 emit: 
  versions      = ch_versions
  bam           = FILTER_BAM.out.bam

}
