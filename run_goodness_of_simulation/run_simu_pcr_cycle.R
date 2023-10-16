library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

library(CellBarcode)

library(devtools)

load_all("../../lib/CellBarcodeSim/")

## barcode sequence and cell clone size

ref = data.table(
    barcode_seq = c(
        "CTCGAGGTCATCGAAGTATCAAGTCCAGTTCTACTATCGTAGCTACTA", 
        "CTCGAGGTCATCGAAGTATCAAGTCCAGTATCGTTACGCTACTA", 
        "CTCGAGGTCATCGAAGTATCAAGTCCAGTCTACTATCGTTACGACAGCTACTA",    
        "CTCGAGGTCATCGAAGTATCAAGTCCAGTTCTACTATCGTTACGAGCTACTA", 
        "CTCGAGGTCATCGAAGTATCAAGTCCATCGTAGCTACTA",  
        "CTCGAGGTCATCGAAGTATCAAGTCCAGTACTGTAGCTACTA",   
        "CTCGAGGTCATCGAAGTATCAAGTCCAGTACTATCGTACTA" 
    ), clone_size = 2^c(1, 7, 3, 4, 5, 6, 2)
)


lib = ref[, .(seq = barcode_seq, freq = clone_size)]
write_tsv(lib, "./data/mef_lib.tsv")

dir("../../data/5290/", "mef", full = T)

## number of cells: 50000, 12500, 3125, 781, 195)
n_cells = round(c(50000, 12500, 3125, 781) / 2)
pre_amp_pcr_cycles = c(0, 2, 4, 6)

n_cells * 2^pre_amp_pcr_cycles

for (i in seq_along(n_cells)) {

    output_prefix = paste0("./tmp/simu_mef_", n_cells[i])

    simulate_main_umi(
        barcode_library_file = "./data/mef_lib.tsv",
        clone_size_dist      = "uniform",
        clone_n              = n_cells[i],
        clone_size_dist_par  = list(size_max = 1, size_min = 1),
        cycle                = pre_amp_pcr_cycles[i] + 30,
        efficiency           = 0.705,
        error                = 1e-5,
        pcr_read_per_umi     = 1000,
        output_prefix        = output_prefix,
        ngs_profile          = "HS20",
        reads_length         = 100,
        top_seq              = "",
        bottom_seq           = "CCGTAGCAAGCTCGAGAGTAGACCTACTGGAATCAGACCGCCACCATGGTGAGCA",
        sequence_trunk       = 60,
        preamp_n             = pre_amp_pcr_cycles[i],
        umi_length           = 12,
        umi_tagging_efficiency = 0.02,
    )
}



