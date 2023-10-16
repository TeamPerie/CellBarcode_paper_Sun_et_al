library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

library(CellBarcode)

library(devtools)

load_all("../../lib/CellBarcodeSim/")

vdj_lib = system.file("data", "vdj_barcodes.tsv", package = "CellBarcodeSim")
d = fread(vdj_lib)
d$seq %<>% paste0("CCTCGAGGTCATCGAAGTATCAAG", .)
write_tsv(d, "tmp/vdj_barcodes.tsv")

## number of cells: 50000, 12500, 3125, 781, 195)
simulate_main_umi(
    barcode_library_file = "./tmp/vdj_barcodes.tsv",
    clone_size_dist      = "lognormal",
    clone_n              = 100,
    clone_size_dist_par  = list(size_mean = 1.2, size_variant = 1),
    cycle                = 30,
    efficiency           = 0.705,
    error                = 1e-6,
    pcr_read_per_umi     = 100,
    output_prefix        = "./tmp/vdj_simu",
    ngs_profile          = "MSv1",
    reads_length         = 111,
    top_seq              = "",
    bottom_seq           = "CCGTAGCAAGCTCGAGAGTAGACCTACTGGAATCAGACCGCCACCATGGTGAGCACACGTCTGAACTCCAGTCACTCAGTCAATCTCGTATGCCGTCTTCTGCTTG",
    sequence_trunk       = 110,
    preamp_n             = 10,
    umi_length           = 16,
    umi_tagging_efficiency = 0.02,
)



