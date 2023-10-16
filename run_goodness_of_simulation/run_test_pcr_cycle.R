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

#' Clone sd and mean

lapply(n_cells, function(x) {
    mean = mean(log(ref$clone_size * x / sum(ref$clone_size)))
    sd = sd(log(ref$clone_size * x / sum(ref$clone_size)))
    c(mean = mean, sd = sd)
})

ref$barcode_seq %<>% sub("CTCGAGGTCATCGAAGTATC", "", .)

#' # Real sample

fq_files = dir("../../data/5290/", "mef", full = T)[c(1, 4, 6, 8)]

metadata <- stringr::str_split_fixed(basename(fq_files), "_", 10)[, c(4, 6)]
metadata <- as.data.frame(metadata)
sample_name <- apply(metadata, 1, paste, collapse = "_")
colnames(metadata) = c("cell_number", "replication")
rownames(metadata) = sample_name

qc = bc_seq_qc(fq_files)
bc_plot_seqQc(qc)

bo_e <- bc_extract(
  fq_files,  # fastq file
  pattern = "([ACGT]{12})CTCGAGGTCATCGAAGTATC([ACGT]+)CCGTAGCAAGCTCGAGAGTAGACCTACT", 
  pattern_type = c("UMI" = 1, "barcode" = 2),
  sample_name = sample_name,
  metadata = metadata
)

bo_e = bc_cure_depth(bo_e)

data.table(bc_meta(bo_e))[, barcode_read_count / raw_read_count]

d_l = bo_e@messyBc
d_l = d_l %>% lapply(function(x) {
    y = sample(x$barcode_seq, 1e5, prob = x$count, replace = T)
    x = data.table(barcode_seq = names(table(y)), count = as.numeric(table(y)))
    x
})

n_cells = round(c(50000, 12500, 3125, 781) / 2)
pre_amp_pcr_cycles = c(10, 12, 14, 16)
total_pcr_cycle = pre_amp_pcr_cycles + 20

g_l = lapply(d_l, function(x) {
    ggplot(x, aes(x = log10(count + runif(length(count))), color = barcode_seq %in% ref$barcode_seq)) +
        geom_histogram(bins = 20, fill = "skyblue", color = "black", alpha = 0.7, boundary = 0) +
        geom_rug(linewidth = 0.5, alpha = 0.5, aes(color = barcode_seq %in% ref$barcode_seq)) +  
        labs(x = "log10 sequence count") +
        theme_bw() + scale_color_manual(values = c("black", "red")) + guides(color = 'none')

})


egg::ggarrange(plots = g_l, nrow = 1)



#' # Simulated sample

pattern = "([ACGT]{12})CTCGAGGTCATCGAAGTATC([ACGT]+)CCGTAGCAAGCTCGAGAGTAGACCTACT"
f_l = dir("./tmp/", pattern = "simu_mef_", full = T) %>% grep(".fq", ., value = T)
bo_s = bc_extract(f_l, pattern, pattern_type = c("UMI" = 1, "barcode" = 2))

qc = bc_seq_qc(f_l)
bc_plot_seqQc(qc)

data.table(bc_meta(bo_s))[, barcode_read_count / raw_read_count]

bo_s = bc_cure_depth(bo_s)
d_l = bo_s@messyBc
d_l = d_l %>% lapply(function(x) {
    y = sample(x$barcode_seq, 1e5, prob = x$count, replace = T)
    x = data.table(barcode_seq = names(table(y)), count = as.numeric(table(y)))
    x
})

n_cells = round(c(50000, 12500, 3125, 781) / 2)
pre_amp_pcr_cycles = c(10, 12, 14, 16)
total_pcr_cycle = pre_amp_pcr_cycles + 20

g_l = lapply(d_l, function(x) {
    ggplot(x, aes(x = log10(count + runif(length(count))), color = barcode_seq %in% ref$barcode_seq)) +
        geom_histogram(bins = 20, fill = "skyblue", color = "black", alpha = 0.7, boundary = 0) +
        geom_rug(linewidth = 0.5, alpha = 0.5, aes(color = barcode_seq %in% ref$barcode_seq)) +  
        labs(x = "log10 sequence count") +
        theme_bw() + scale_color_manual(values = c("black", "red")) + guides(color = 'none')
})


egg::ggarrange(plots = g_l, nrow = 1)

#' # Barcode base percentage heatmap
f1 = "./tmp/simu_mef_6250.fq"
f2 = "../../data/5290/5290_4_BCM_12500_mef_mixb_TAACTGC_S8_R1_001.fastq.gz"

x_qc = bc_seq_qc(c(f1, f2))
bc_plot_seqQc(x_qc)

bo = bc_extract(c(f1, f2), pattern, pattern_type = c("UMI" = 1, "barcode" = 2))

x = bc_meta(bo)
x[, 2] / x[, 1]
