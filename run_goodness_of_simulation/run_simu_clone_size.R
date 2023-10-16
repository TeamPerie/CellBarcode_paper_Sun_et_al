library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(diptest)

library(CellBarcode)

library(devtools)
load_all("../../lib/CellBarcodeSim/")


## small clone size = 0.5, big clone size = 2.5

clone_size_sd = c(0, 0.5, 1, 1.5, 2, 2.5, 3)
clone_n = lapply(clone_size_sd, function(x) {
    3000 / mean(ceiling(rlnorm(1000, 1, x)))
}) %>% unlist %>% round()

true_v = fread("./data/Mathild_gRNA_wenjie_batch3_2014013/ref_barcode.tsv")[[1]] %>% substring(1, 19)
lib = data.frame(seq = true_v, freq = 1)
write_tsv(lib, "./tmp/barcode_library_gRNA.tsv")

## Simulate a Barcode sequencing experiment
for (i in seq_along(clone_size_sd)) {

    output_prefix = paste0("./tmp/barcode_sim_random_miseq_sd_", clone_size_sd[i])

    simulate_main(
        barcode_library_file = "./tmp/barcode_library_gRNA.tsv",  ## Define the barcode library
        clone_n              = clone_n[i],  ## Define the number of clones
        clone_size_dist      = "lognormal", ## Define the clone size distribution
        clone_size_dist_par  = list(size_mean = 1, size_variant = clone_size_sd[i]),  ## Define the parameters of the clone size distribution
        cycle                = 30,  ## Define the number of cell cycles
        efficiency           = 0.705, ## Define the efficiency of the cell cycle
        error                = 1e-5,  ## Define the error rate of the PCR, mutation per base per cycle
        pcr_read_per_cell    = 10,  ## Define the number of PCR reads per cell (clone_n)
        output_prefix        = output_prefix ,  ## Define the output prefix
        ngs_profile          = "MSv1",  ## Define the NGS profile (refer to ART sequencing simulator manual)
        reads_length         = 100,  ## Define the length of the reads
        is_replicate         = F,  ## Define whether to divide the reads into two replicates
        top_seq              = "AAGGACGAAACACCG",  ## Define the fixed sequence at the 5' end of the reads to be added
        bottom_seq           = "ATCGATCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        sequence_trunk       = 100, ## Define the length of the fixed sequence to be added
        art_bin              = NULL ## Use the default ART binary file
    )
}

