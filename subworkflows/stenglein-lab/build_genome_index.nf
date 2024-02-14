include { BOWTIE2_BUILD               } from '../../modules/nf-core/bowtie2/build/main'

workflow BUILD_GENOME_INDEX {

 take:

  genome_fasta

 main:

  // define some empty channels for keeping track of stuff
  ch_versions     = Channel.empty()                                               

  // prepare a genome index 
  // the nf-core bowtie2 index module has a tuple input val(meta), path(fasta), so make 
  // some a metadata string to make it happy
  ch_index = Channel.value("bowtie2-index-meta").map{it -> [it, file(genome_fasta)]}

  BOWTIE2_BUILD (ch_index)
  ch_versions = ch_versions.mix ( BOWTIE2_BUILD.out.versions )      
  
  
 emit: 
  versions      = ch_versions
  index         = BOWTIE2_BUILD.out.index
  fasta         = ch_index

}
