library(CellBarcode)
library(stringr)
library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(ggrepel)
library(ggpubr)

x_known_barcode = data.table(
  clone_id = c("3", "6", "14", "17", "30", "46", "53"),
  barcode_seq = c(
    "AAGTCCAGTATCGTTACGCTACTA",
    "AAGTCCAGTACTATCGTACTA",
    "AAGTCCAGTCTACTATCGTTACGACAGCTACTA",
    "AAGTCCATCGTAGCTACTA",
    "AAGTCCAGTACTGTAGCTACTA",
    "AAGTCCAGTTCTACTATCGTTACGAGCTACTA",
    "AAGTCCAGTTCTACTATCGTAGCTACTA"
    ),
  clone_size_factor = c(7, 2, 3, 5, 6, 4, 1), 
  is_known = T
)

# write_tsv(x_known_barcod, "./tmp/MEF_known_barcode.tsv")
cell_n = c(195, 781, 3125, 12500, 50000)

#+ eval=F
fastq_files = dir("../../data/5290/", "mef", full = T)
pattern = "([ATCG]{12})CTCGAGGTCATCGAAGTATC([ACTG]+)CCGTAGCAAGCTCGAGAGTAGAC"
sample_names = fastq_files %>% basename %>% str_match("BCM_\\d+_mef_mix.") %>% extract(, 1)
bc_obj = bc_extract(fastq_files, pattern = pattern, sample_name = sample_names, pattern_type = c(UMI=1, barcode = 2))

save(bc_obj, file = "tmp/bc_obj.RData")

#+ eval=T


#' # Explore threshold
d_plot = lapply(seq(1, 200, 5), function(i) {
    print(i)
    bc_obj_UMI = bc_cure_umi(bc_obj[, 4], depth = i, isUniqueUMI = F) 
    bc_obj_UMI = bc_cure_depth(bc_obj_UMI, depth = 10)
    d = bc_2dt(bc_obj_UMI) 
    d = d[, .(barcode_seq = unique(barcode_seq)), by = sample_name][, .N, by = sample_name]
    d$filter = i
    d
})

d_plot <- rbindlist(d_plot)
d_plot_plus = d_plot[, .N, by = filter] ## Scale the left and right y-axis

ggplot(d_plot) + aes(x = filter, y = N) +
    geom_point() + theme_bw() +
    scale_y_log10() + labs(x = "UMI depth threshold", y = "Number of barcodes")

#' # Apply threshold
d = bc_cure_umi(bc_obj[, 4], depth = 100) %>% bc_cure_depth(10) %>% bc_2dt()
d
ann = x_known_barcode[, clone_size := 125000 * (2^clone_size_factor) / sum(2^clone_size_factor)]
d[, clone_size_estimated := 125000 * count / sum(count)]
d_plot = merge(d, ann, by = "barcode_seq", all = T)
ggplot(d_plot) + aes(x = clone_size, y = clone_size_estimated) +
    geom_point(position = position_jitter(w = 0.3, h = 0)) +
    geom_smooth(method = "lm", se = F) +
    stat_regline_equation() +
    scale_x_log10() + scale_y_log10() + coord_fixed() + theme_bw() +
    labs(x = "Clone size (Expected)", y = "Clone size (Estimated)")

