library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

library(CellBarcode)

library(devtools)

load_all("../../lib/CellBarcodeSim/")

ref = readLines("../../data/Libraries/new_lib/Ref_Lib_20bp/LG22_filtered.fa") %>% grep(">", ., invert=T, value=T)
ref[1] %>% nchar
x = table(ref)
ref_tab = data.table(seq = names(x), freq = as.integer(x))
write_tsv(ref_tab, "./tmp/ref_lib.tsv")

simulate_main(
        barcode_library_file = "./tmp/ref_lib.tsv",
        clone_size_dist      = "lognormal",
        clone_n              = 15,
        clone_size_dist_par  = list(size_mean = 1.2, size_variant = 3),
        cycle                = 30,
        efficiency           = 0.705,
        error                = 1e-6,
        pcr_read_per_cell    = 50,
        output_prefix        = "./tmp/lenti_simu",
        ngs_profile          = "HS20",
        reads_length         = 65,
        top_seq              = "ACGGAATGCTAGAACACTCGAGATCAG",
        bottom_seq           = "ATGTGGTATGATGTATCA",
        sequence_trunk       = 60
    )

system('grep -B 1 -A 2 -E "^ACGGAATG" tmp/lenti_simu.fq | grep -v "^--$" > tmp/lenti_sumu_filtered.fq')
