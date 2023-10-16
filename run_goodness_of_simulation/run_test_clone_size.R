library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(diptest)

library(CellBarcode)

#' ## gRNA dataset

#' Sample R31: 46h_Brice_bank_LRT R34: 7d_Brice_bank_LRT There are around 1800 barcodes.

f = dir("./data/Mathild_gRNA_wenjie_batch3_2014013/", recursive = T, pattern = "fastq.gz", full=T) %>% grep("trimmed|unde", ., value=T, invert=T)
f = f %>% grep("R31|R34", ., value=T)
f

pattern = "AAGGACGAAACACCG(.{19})"
bo_gRNA = bc_extract(f, pattern, sample_name = basename(f) %>% str_replace(".R1.fastq.gz", ""))


d_l = bo_gRNA@messyBc
true_v = fread("./data/Mathild_gRNA_wenjie_batch3_2014013/ref_barcode.tsv")[[1]] %>% substring(1, 19)

#' Clone SD

d_l[[1]][barcode_seq %in% true_v]$count %>% log %>% sd
d_l[[2]][barcode_seq %in% true_v]$count %>% log %>% sd


#' Visualizaion
d_l = d_l %>% lapply(function(x) {
    y = sample(x$barcode_seq, 1e5, prob = x$count, replace = T)
    x = data.table(barcode_seq = names(table(y)), count = as.numeric(table(y)))
    x
})




g_l = lapply(d_l, function(d) {
    ggplot(d, aes(x = log10(count + runif(length(count))))) +
        geom_histogram(bins = 20, fill = "skyblue", color = "black", alpha = 0.7) +
        geom_rug(aes(color = barcode_seq %in% true_v), linewidth = 0.5, alpha = 0.7) +
        scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +  # Set custom colors
        theme_bw() + guides(color = 'none') +
        labs(x = "log10 sequence count")
})

egg::ggarrange(g_l[[1]], g_l[[2]], ncol = 2, nrow = 1)

#' ## gRNA simulation

#' ### Analyze the simulated dataset

f = dir("./tmp/", recursive = T, pattern = "fq", full=T) %>% grep("miseq_sd", ., value=T)
f = f[c(2, 1, 4, 3, 6, 5, 7)]
f

pattern = "AAGGACGAAACACCG(.{19})"

bo_gRNA_simu = bc_extract(f, pattern, sample_name = basename(f) %>% str_replace(".fq", ""))

bo = bo_gRNA_simu

d_l = bo_gRNA_simu@messyBc
d_l = d_l %>% lapply(function(x) {
    y = sample(x$barcode_seq, 1e5, prob = x$count, replace = T)
    x = data.table(barcode_seq = names(table(y)), count = as.numeric(table(y)))
    x
})

g_l = lapply(d_l, function(d) {
    ggplot(d, aes(x = log10(count + runif(length(count))))) +
        geom_histogram(bins = 20, fill = "skyblue", color = "black", alpha = 0.7) +
        geom_rug(aes(color = barcode_seq %in% true_v), linewidth = 0.5, alpha = 0.7) +
        scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +  # Set custom colors
        theme_bw() + guides(color = 'none') +
        labs(x = "log10 sequence count")
})

egg::ggarrange(plots = g_l,  nrow = 1)



#' # Apply statistics

#' ## Barcode extracting efficiency

d_meta_gRNA = bc_meta(bo_gRNA) %>% data.table(keep.rownames = T)
d_meta_gRNA_simu = bc_meta(bo_gRNA_simu) %>% data.table(keep.rownames = T)

d_meta_gRNA[, .(barcode_ratio = barcode_read_count / raw_read_count, rn)]
d_meta_gRNA_simu[, .(barcode_ratio = barcode_read_count / raw_read_count, rn)]
