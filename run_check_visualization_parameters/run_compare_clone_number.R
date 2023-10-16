library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(CellBarcode)

library(devtools)

load_all("../../lib/CellBarcodeSim/")

#+ eval=F
time_start <- Sys.time()
simulate_main(
    barcode_library_file = "../lib/CellBarcodeSim/example_barcode_library/random_barcodes.tsv",
    clone_size_dist      = "lognormal",
    clone_n              = 300,
    clone_size_dist_par  = list(size_mean = 1.2, size_variant = 1),
    cycle                = 15,
    efficiency           = 0.705,
    error                = 1e-6,
    pcr_read_per_cell    = 50,
    output_prefix        = paste0("./tmp/clone_n_300"),
    ngs_profile          = "HS20",
    reads_length         = 100,
    is_replicate         = F,
    top_seq              = "AAAAAAAAAAGGGGG",
    bottom_seq           = "ATCGATCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
    sequence_trunk       = 110,
    art_bin              = "../lib/art_bin_MountRainier/art_illumina"
)
Sys.time() - time_start


time_start <- Sys.time()
simulate_main(
    barcode_library_file = "../lib/CellBarcodeSim/example_barcode_library/random_barcodes.tsv",
    clone_size_dist      = "lognormal",
    clone_n              = 3000,
    clone_size_dist_par  = list(size_mean = 1.2, size_variant = 1),
    cycle                = 15,
    efficiency           = 0.705,
    error                = 1e-6,
    pcr_read_per_cell    = 50,
    output_prefix        = paste0("./tmp/clone_n_3k"),
    ngs_profile          = "HS20",
    reads_length         = 100,
    is_replicate         = F,
    top_seq              = "AAAAAAAAAAGGGGG",
    bottom_seq           = "ATCGATCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
    sequence_trunk       = 110,
    art_bin              = "../lib/art_bin_MountRainier/art_illumina"
)
Sys.time() - time_start

time_start <- Sys.time()
simulate_main(
    barcode_library_file = "../lib/CellBarcodeSim/example_barcode_library/random_barcodes.tsv",
    clone_size_dist      = "lognormal",
    clone_n              = 30000,
    clone_size_dist_par  = list(size_mean = 1.2, size_variant = 1),
    cycle                = 15,
    efficiency           = 0.705,
    error                = 1e-6,
    pcr_read_per_cell    = 50,
    output_prefix        = paste0("./tmp/clone_n_30k"),
    ngs_profile          = "HS20",
    reads_length         = 100,
    is_replicate         = F,
    top_seq              = "AAAAAAAAAAGGGGG",
    bottom_seq           = "ATCGATCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
    sequence_trunk       = 110,
    art_bin              = "../lib/art_bin_MountRainier/art_illumina"
)
Sys.time() - time_start

#+ eval=T
f = dir("./tmp/", pattern = "fq", full.names = T)

pattern = c("AAAAAAAAAAGGGGG(.*)ATCGATCGTTTTTTTT")
bo = bc_extract(f, pattern)
d = bo@messyBc

d[[1]]$count %>% log10 %>% hist
d[[2]]$count %>% log10 %>% hist
d[[3]]$count %>% log10 %>% hist

d = rbindlist(d, idcol= T)

# ggplot(d) + aes(x = count, color = .id) + stat_ecdf()
# ggplot(d) + aes(x = log(count), color = .id) + stat_ecdf()

ggplot(d) + aes(x = count, color = .id) + geom_histogram() + geom_rug()
ggplot(d) + aes(x = log(count), color = .id) + stat_ecdf()


