library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

library(CellBarcode)

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

n_cells = round(c(50000, 12500, 3125, 781) / 2)

#' # Simulated sample

pattern = "([ACGT]{12})CTCGAGGTCATCGAAGTATC([ACGT]+)CCGTAGCAAGCTCGAGAGTAGACCTACT"
f_l = dir("./tmp/", pattern = "simu_mef_", full = T) %>% grep(".fq", ., value = T)
bo_s = bc_extract(f_l, pattern, pattern_type = c("UMI" = 1, "barcode" = 2))

data.table(bc_meta(bo_s))[, barcode_read_count / raw_read_count]

bo_s = bc_cure_depth(bo_s)
d_l = bo_s@messyBc
d_l = d_l %>% lapply(function(x) {
    y = sample(x$barcode_seq, 1e5, prob = x$count, replace = T)
    x = data.table(barcode_seq = names(table(y)), count = as.numeric(table(y)))
    x
})

n_cells = round(c(50000, 12500, 3125, 781) / 2)


x = bc_cure_umi(bo_s, 100)[, 1] %>% bc_2dt
plot(y = x$count %>% sort, x = 2^c(1:7), xlab = "Expected clone size", ylab = "Simulated clone size")


