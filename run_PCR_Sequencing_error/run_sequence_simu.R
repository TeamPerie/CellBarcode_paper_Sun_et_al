library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(devtools)
library(CellBarcode)
library(Biostrings)

load_all("../../lib/CellBarcodeSim/")

#' # Simulation

barcode_library_file = system.file("data", "random_barcodes.tsv", package = "CellBarcodeSim")
dir.create("tmp", showWarnings = FALSE)

## Simulate a Barcode sequencing experiment
simulate_main(
    barcode_library_file = barcode_library_file,  ## Define the barcode library
    clone_n              = 20,  ## Define the number of clones
    clone_size_dist      = "uniform", ## Define the clone size distribution
    clone_size_dist_par  = list(size_max = 1000, size_min = 1),  ## Define the parameters of the clone size distribution
    cycle                = 30,  ## Define the number of cell cycles
    efficiency           = 0.705, ## Define the efficiency of the cell cycle
    error                = 1e-6,  ## Define the error rate of the PCR, mutation per base per cycle
    pcr_read_per_cell    = 50,  ## Define the number of PCR reads per cell (clone_n)
    output_prefix        = "./tmp/barcode_sim" ,  ## Define the output prefix
    ngs_profile          = "MSv1",  ## Define the NGS profile (refer to ART sequencing simulator manual)
    reads_length         = 35,  ## Define the length of the reads
    is_replicate         = F,  ## Define whether to divide the reads into two replicates
    top_seq              = "AAAAAAAAAAGGGGG",  ## Define the fixed sequence at the 5' end of the reads to be added
    bottom_seq           = "TTTTTTTTTT",  ## Define the fixed sequence at the 3' end of the reads to be added
    sequence_trunk       = 10, ## Define the length of the fixed sequence to be added
    art_bin              = NULL ## Use the default ART binary file
)


#' # Compare PCR and Sequencing error
pattern = "AAAAAAAAAAGGGGG(.{10})TTTTTTTTTT"

f_fa = "./tmp/barcode_sim_library.fasta"
f_fq = "./tmp/barcode_sim.fq"
ds = Biostrings::readDNAStringSet(f_fa)

d_fa = bc_extract(ds, pattern)
d_fq = bc_extract(f_fq, pattern)@messyBc[[1]]

d = merge(d_fa, d_fq, by = "barcode_seq", all = T)
d[is.na(d)] = 0
d


## true sequence
true_seq = fread("./tmp/barcode_sim_ref.tsv")

hist(log10(d_fa$count), main = "Sequence dist after PCR", xlab = "log10 sequence count", ylim = c(0, 300))
d_fa$count = d_fa$count + runif(nrow(d_fa), 0, 1)
rug(log10(d_fa$count))
rug(log10(d_fa[barcode_seq %in% true_seq$seq, count]), col = "red")


hist(log10(d_fq$count), main = "Sequence dist after NGS", xlab = "log10 sequence count", ylim = c(0, 300))
d_fq$count = d_fq$count + runif(nrow(d_fq), 0, 1)
rug(log10(d_fq$count))
rug(log10(d_fa[barcode_seq %in% true_seq$seq, count]), col = "red")
