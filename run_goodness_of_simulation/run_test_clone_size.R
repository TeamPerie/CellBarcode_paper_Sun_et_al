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

#' Returned barcodes in list

bo_gRNA= bc_cure_depth(bo_gRNA, depth = 100, isUpdate = F)
d_l = bo_gRNA@cleanBc
d_l = d_l %>% lapply(function(x) {
    y = sample(x$barcode_seq, 1e5, prob = x$count, replace = T)
    x = data.table(barcode_seq = names(table(y)), count = as.numeric(table(y)))
    x
})

lapply(d_l, function(x) {
    nrow(x[barcode_seq %in% true_v]) / nrow(x)
})



#' Read median

d_l[[1]][!(barcode_seq %in% true_v)]$count %>% quantile(0.50) / d_l[[1]][barcode_seq %in% true_v]$count %>% quantile(0.50)

d_l[[2]][!(barcode_seq %in% true_v)]$count %>% quantile(0.50) / d_l[[2]][barcode_seq %in% true_v]$count %>% quantile(0.50)


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

#' Read median

d_l[[2]][!(barcode_seq %in% true_v)]$count %>% quantile(0.5) / d_l[[2]][barcode_seq %in% true_v]$count %>% quantile(0.5) 

d_l[[4]][!(barcode_seq %in% true_v)]$count %>% quantile(0.5) / d_l[[4]][barcode_seq %in% true_v]$count %>% quantile(0.5) 

d_l[[5]][!(barcode_seq %in% true_v)]$count %>% quantile(0.5) / d_l[[5]][barcode_seq %in% true_v]$count %>% quantile(0.5) 


#' Returned barcodes in list
bo_gRNA_simu = bc_cure_depth(bo_gRNA_simu, depth = 10, isUpdate = F)
d_l = bo_gRNA_simu@cleanBc
d_l = d_l %>% lapply(function(x) {
    y = sample(x$barcode_seq, 1e5, prob = x$count, replace = T)
    x = data.table(barcode_seq = names(table(y)), count = as.numeric(table(y)))
    x
})

lapply(d_l, function(x) {
    nrow(x[barcode_seq %in% true_v]) / nrow(x)
})

d_l



#' # Apply statistics

#' ## Barcode extracting efficiency

d_meta_gRNA = bc_meta(bo_gRNA) %>% data.table(keep.rownames = T)
d_meta_gRNA_simu = bc_meta(bo_gRNA_simu) %>% data.table(keep.rownames = T)

d_meta_gRNA[, .(barcode_ratio = barcode_read_count / raw_read_count, rn)]
d_meta_gRNA_simu[, .(barcode_ratio = barcode_read_count / raw_read_count, rn)]


#' # Compare reads between simulation and real data

exp_f = "./data/Mathild_gRNA_wenjie_batch3_2014013//L412R31/L412R31.R1.fastq.gz"
sim_f = "./tmp/barcode_sim_random_miseq_sd_1.fq"

exp_b = bc_extract(exp_f, "TTGTGGA(.{50})", sample_name = "exp")
d_input = exp_b@messyBc[[1]]
names(d_input) = c("seq", "freq")

## barcode rate
barcode_rate = bc_extract(list(d_input), pattern) %>% bc_meta() %>% data.table(keep.rownames = T)
barcode_rate[, barcode_read_count / raw_read_count]

bc_extract(sim_f, pattern) %>% bc_meta() %>% data.table(keep.rownames = T) %>% 
    .[, barcode_read_count / raw_read_count]


## Basepair frequency per cycle
# calculate basepair frequency per cycle for Simulation data
sim_b = bc_seq_qc(sim_f) 
d = data.table(sim_b@base_freq_per_cycle)[Cycle <= 50]
d$fileName = "Simulation"
d[, base_num := sum(Count), by = .(Cycle, fileName)]
d[, base_percent := Count / base_num, by = .(Base, Cycle, fileName)]

# Calculate basepair frequency per cycle for Experiment data
expanded_data <- do.call(rbind, apply(d_input, 1, function(row) {
  data.frame(Base = unlist(strsplit(as.character(row["seq"]), "")), 
             Cycle = seq_along(strsplit(as.character(row["seq"]), "")[[1]]), 
             Count = as.numeric(row["freq"]))
}))
d_real = aggregate(Count ~ Base + Cycle, data = expanded_data, FUN = sum)
d_real = data.table(d_real)

d_real$fileName = "Experiment"
d_real[, base_num := sum(Count), by = .(Cycle)]
d_real[, base_percent := Count / base_num, by = .(Base, Cycle, fileName)]

# Merge data
d = rbind(d, d_real)

## Plot & statistcs
d_stat = d[Cycle <= 34, .(Base = Base[which.max(base_percent)], base_percent = max(base_percent)), by = .(Cycle, fileName)]
d = d[Cycle <= 34]

m = dcast(d_stat, Cycle ~ fileName, value.var = "base_percent") 
cor.test(m[[2]], m[[3]])

# Heatmap
ggplot(d, aes(Cycle, fileName, fill = Base, alpha = base_percent)) + 
    geom_tile(color = "black") + labs(y = "Sample Name") + theme_bw()

# Line plot
ggplot(d, aes(x = Cycle, y = base_percent, color = Base)) + 
    geom_line() + facet_wrap(~fileName, ncol = 1) + theme_classic()

# Stack bar plot
ggplot(d, aes(x = Cycle, y = base_percent, fill = Base)) + 
    geom_bar(stat = "identity") + facet_wrap(~fileName, ncol = 1) + theme_classic()
    # scale_fill_manual(values = c("#532026FF", "#BA141EFF", "#E2E3E7FF", "#61829CFF"))
write_tsv(d, "./tmp/Figure_2B.tsv")
ggsave("./tmp/Figure_2B.pdf", width = 5, height = 4)

