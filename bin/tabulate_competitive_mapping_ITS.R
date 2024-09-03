#!/usr/bin/env Rscript

#
# tabulate the fractions of reads mapping to different reference 
# sequences that are assigned to a host species
# from a competitive mapping analysis
#
# MDS 6/3/2024
#

library(tidyverse)

# This code block sets up input arguments to either come from the command line
# (if running from the pipeline, so not interactively) or to use expected default values
# (if running in interactive mode, for example in RStudio for troubleshooting
# or development).
#
if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  tsv_input            = args[1]
  contig_host_map_file = args[2]
  metadata_input       = args[3]
  R_lib_dir            = args[4]
  R_script_dir         = args[5]
  
  output_dir           = "./"
} else {
  # if running via RStudio
  tsv_input            = "../results_its/collected_coverage.txt"
  contig_host_map_file = "../refseq/host_contig_map_its.txt"
  metadata_input       = "../../metadata/collected_metadata.rds"
  co1_species_id       = "../../species_id/results_species_id/process/co1_based_species_assignments.txt"
  R_lib_dir            = "../lib/R/"
  R_script_dir         = "../../scripts/"
  output_dir           = "../results_its/"
}


# these libraries are not part of the standard tidyverse, so may have to load it 
# from a specified path
# either from pipeline's R lib dir or from R environment
# if (R_lib_dir != "NA") {
  # library(rstatix, lib.loc=R_lib_dir)
  # library(ggpubr, lib.loc=R_lib_dir)
  # library(patchwork, lib.loc=R_lib_dir)
# } else {
  # in this case assuming these will be installed
  # library(rstatix)
  # library(ggpubr)
  # library(patchwork)
# }

source(paste0(R_script_dir, "/plot_colors.R"))
fancy_color_scale <- c(fresh_frozen_color, experimental_dried_color, old_collection_color)

# read in metadata
metadata <- readRDS(metadata_input)

# read in the ITS-mapping data
df <- read.delim(tsv_input, sep="\t", header=F)

# name columns: this is samtools coverage output format
# see: http://www.htslib.org/doc/samtools-coverage.html
colnames(df) <- c("sample_id", "rname", "startpos", "endpos", "numreads", 
                  "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq")


# read in tab-delimited file that maps refseq contigs to a host species
contig_host_map <- read.delim(contig_host_map_file, sep="\t", header=F)
colnames(contig_host_map) <- c("rname", "contig_species")

# change species names for nicer plotting 
contig_host_map$contig_species <- str_replace(contig_host_map$contig_species, "_", "\n")

# join in contig->species mapping
df <- left_join(df, contig_host_map)

# read in COX1-based species assignments and join in
co1 <- read.delim(co1_species_id, sep="\t", header=T)
co1 <- co1 %>% select(dataset, co1_assigned_species) %>% rename(sample_id = dataset)
df <- left_join(df, co1)

# tabulate # of mapped reads per dataset per species
per_sample_per_species <- df %>% 
  group_by(sample_id, contig_species) %>% 
  summarize(per_species_reads = sum(numreads), .groups="drop")

# total # of reads per dataset
per_sample <- df %>% group_by(sample_id) %>% summarize(total_dataset_reads = sum(numreads))

per_sample_per_species <- left_join(per_sample_per_species, per_sample)

per_sample_per_species <- per_sample_per_species %>% mutate(fraction_reads = per_species_reads / total_dataset_reads)

# calculate Shannon entropy of each dataset (ecological evenness)
per_sample_per_species <- per_sample_per_species %>% mutate(pi_lnpi = if_else(fraction_reads == 0, 0, fraction_reads * log(fraction_reads)))
shannons <- per_sample_per_species %>% group_by(sample_id) %>% summarize(shannon = -1 * sum(pi_lnpi))

# join in Shannon 
per_sample_per_species <- left_join(per_sample_per_species, shannons)

# rejoin in metadata
per_sample_per_species <- left_join(per_sample_per_species, metadata, by="sample_id")

# rejoin in COX1-assigned species
per_sample_per_species <- left_join(per_sample_per_species, co1)

# reorder samples by collection date
# per_sample_per_species$sample_id <- fct_reorder(per_sample_per_species$sample_id, 
                                                # per_sample_per_species$date_collected, .na_rm = F, .desc=F)


# reorder species so they display in order of most to fewest samples
per_sample_per_species$contig_species <- fct_relevel(per_sample_per_species$contig_species,
                                                     "Drosophila\nmelanogaster",
                                                     "Drosophila\nsimulans",
                                                     "Drosophila\nputridfa",
                                                     "Scaptoymza\npallida")

make_species_plot <- function(df_to_plot) {
  ggplot(df_to_plot) +
    geom_tile(aes(y=sample_id_in_paper, x=contig_species, fill=fraction_reads)) +
    geom_text(data=filter(df_to_plot, per_species_reads >0),
              aes(y=sample_id_in_paper, x=contig_species, label=per_species_reads), size=2, angle=0, color="white") +
    scale_fill_gradient(low="white", high="darkslateblue", na.value = "grey95") +
    theme_this_paper(base_size=10) +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.title.x = element_blank()) +
    ylab("") 
}

add_x_legend <- function(p) {
  p + 
    theme(axis.text.x = element_text(),
          axis.title.x = element_text(size=9)) +
    xlab("Number of reads mapping to indicated species rRNA ITS") 
}

oc_mel <- filter(per_sample_per_species, sample_type == "Museum\nsamples" & co1_assigned_species == "Drosophila_melanogaster")
oc_sim <- filter(per_sample_per_species, sample_type == "Museum\nsamples" & co1_assigned_species == "Drosophila_simulans")
oc_put <- filter(per_sample_per_species, sample_type == "Museum\nsamples" & co1_assigned_species == "Drosophila_putrida")
oc_sca <- filter(per_sample_per_species, sample_type == "Museum\nsamples" & co1_assigned_species == "Scaptomyza_pallida")
oc_mel_p <- make_species_plot(oc_mel)
oc_sim_p <- make_species_plot(oc_sim)
oc_put_p <- make_species_plot(oc_put)
oc_sca_p <- make_species_plot(oc_sca)
oc_sca_p <- add_x_legend(oc_sca_p)

# do this so each sub-part of plot is proportional to the # of samples 
n_each <- c(nrow(oc_mel %>% group_by(sample_id) %>% summarize()),
            nrow(oc_sim %>% group_by(sample_id) %>% summarize()),
            nrow(oc_put %>% group_by(sample_id) %>% summarize()),
            nrow(oc_sca %>% group_by(sample_id) %>% summarize()) + 0.25) # extra because x-axis text

oc_mel_p + oc_sim_p + oc_put_p + oc_sca_p + plot_layout(ncol=1, heights=n_each)

ggsave(paste0(output_dir, "Fig_12_ITS_mapping_read_counts.pdf"), units="in", width=7, height=6)

# oc_p <- make_species_plot(filter(per_sample_per_species, sample_type == "Museum\nsamples"))
            
ff   <- filter(per_sample_per_species, sample_type == "Fresh\nfrozen")
ed   <- filter(per_sample_per_species, sample_type == "Experimental\ndried")
ctrl <- filter(per_sample_per_species, str_detect(sample_type,"control"))

n_each <- c(nrow(ff %>% group_by(sample_id) %>% summarize()),
            nrow(ed %>% group_by(sample_id) %>% summarize()),
            nrow(ctrl %>% group_by(sample_id) %>% summarize()))
            
ff_p <-   make_species_plot(ff)
ed_p <-   make_species_plot(ed)
ctrl_p <- make_species_plot(ctrl)

ff_p + ed_p + add_x_legend(ctrl_p) + plot_layout(ncol = 1, heights=n_each)

ggsave(paste0(output_dir, "Fig_S_11_ITS_mapping_read_counts_in_controls.pdf"), units="in", width=7, height=10)





