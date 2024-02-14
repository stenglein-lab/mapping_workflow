include { SAMTOOLS_STATS                         } from '../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_COVERAGE                      } from '../../modules/nf-core/samtools/coverage/main'
include { SAMTOOLS_DEPTH                         } from '../../modules/nf-core/samtools/depth/main'

/*
  Calculate basic mapping summary statistics from bam
 */
workflow MAPPING_STATS {

 take:
  bam             // [meta, bam]
  genome_fasta    // [meta, fasta]
  per_base_depth  // boolean: calculate per-base coverage depth values using samtools depth?

 main:

  // define some empty channels for keeping track of stuff
  ch_versions     = Channel.empty()                                               

  // ------------------
  // samtools stats
  // ------------------

  // re: using a fake name for a missing path value:
  // https://github.com/nextflow-io/nextflow/issues/1532
  SAMTOOLS_STATS(bam, genome_fasta)
  ch_versions = ch_versions.mix ( SAMTOOLS_STATS.out.versions )      

  // ------------------
  // samtools coverage
  // ------------------

  // add an extra null argument to end of bam input because samtools sort expects a 3rd index
  SAMTOOLS_COVERAGE(bam)
  ch_versions = ch_versions.mix ( SAMTOOLS_COVERAGE.out.versions )      

  // ------------------
  // samtools depth
  // ------------------
  // this outputs per-based coverage depth values
  // only run if requested to do so

  ch_depth = Channel.empty()
  if (per_base_depth) {
    // the second null input is placeholder for a possible interval bedfile
    SAMTOOLS_DEPTH(bam)
	 ch_depth    = ch_depth.mix    ( SAMTOOLS_DEPTH.out.depth )
    ch_versions = ch_versions.mix ( SAMTOOLS_DEPTH.out.versions )      
  }

 emit: 
  versions      = ch_versions
  stats         = SAMTOOLS_STATS.out.stats
  coverage      = SAMTOOLS_COVERAGE.out.coverage
  depth         = ch_depth

}
