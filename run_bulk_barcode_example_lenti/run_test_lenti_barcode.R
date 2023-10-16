library(CellBarcode)
library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(ggpubr)

theme0 <- theme_bw() + theme(
    text = element_text(size = 15),
    line = element_line(linewidth = 1),
    axis.line = element_line(linewidth = 1),
    axis.ticks.length = unit(3, units = "mm"),
    axis.text.x = element_text(
        margin = margin(t = 2, unit = "mm")
        , angle = 60, vjust = 1, size = 15, hjust = 1),
    axis.text.y = element_text(margin = margin(r = 3, l = 5, unit = "mm")),
    legend.position = "right",
) 
theme1 = theme0 + theme(
    axis.text.x = element_text( 
        margin = margin(t = 2, unit = "mm")
        , angle = 0, vjust = 1, size = 12, hjust = 0.5)
)


#' # Reference sequence
ref = readLines("../../data/Libraries/new_lib/Ref_Lib_20bp/LG22_filtered.fa") %>% grep(">", ., invert=T, value=T)

## check the base pair length
ref[1] %>% nchar

## number of unique barcodes
length(ref)
unique(ref) %>% length

## The filtered sample based on plate index
# grep -B 1 -A 2 -E "^ACGGAATG" A1007R99.R1.fastq | grep -v "^--$" > sample1_filtered.fastq
# grep -B 1 -A 2 -E "^ACGGAATG" A1007R100.R1.fastq | grep -v "^--$" > sample2_filtered.fastq
f = c("../../data/A1007/sample1_filtered.fastq", "../../data/A1007/sample2_filtered.fastq")


## Check sequencing quality
x = bc_seq_qc(f)
bc_plot_seqQc(x)
pdf("./tmp/fig_bc_seqQc_rep2.pdf", height = 5, width = 7)
bc_plot_seqQc(x[2])
dev.off()
pdf("./tmp/fig_bc_seqQc_rep1.pdf", height = 5, width = 7)
bc_plot_seqQc(x[1])
dev.off()


#' # index 8: ACGGAATG
pattern = "ACGGAATGCTAGAACACTCGAGATCAG(.{20})ATGTGGTATGATGTATC"

time_start = Sys.time()
bc_obj = bc_extract(f, pattern)
Sys.time() - time_start

## apply filtering
bc_obj1 = bc_cure_depth(bc_obj, depth = 0, isUpdate =F)
bc_obj2 = bc_cure_cluster(bc_obj1, count_threshold = 1e7)

d = bc_2dt(bc_obj2)
d[, count_sum := sum(count), by = sample_name]

## Normalize cell number
d[, cell_count := 13564 / 2 / count_sum * count]

## Prepare for plotting
d_plot_cell_count = dcast(d, barcode_seq ~ sample_name, sep = "", value.var="cell_count")
d_plot_read_count = dcast(d, barcode_seq ~ sample_name, sep = "", value.var="count")
d_plot_cell_count[is.na(d_plot_cell_count)] = 0
d_plot_read_count[is.na(d_plot_read_count)] = 0

d_plot_cell_count$in_ref = "No"
d_plot_cell_count[barcode_seq %in% ref, in_ref:= "Yes"]

d_plot_read_count$in_ref = "No" 
d_plot_read_count[barcode_seq %in% ref, in_ref:= "Yes"]

#' # Log SD
d_plot_cell_count[in_ref == "Yes"][sample1_filtered.fastq > 1][, sample1_filtered.fastq] %>% log %>% sd

ggplot(d_plot_read_count) + aes(x = sample1_filtered.fastq + 1, y = sample2_filtered.fastq + 1, color = in_ref) +
    geom_point() + scale_y_log10() + scale_x_log10() +
    theme_bw() + scale_color_manual(values = c("No" = "red", "Yes" = "black"))

ggplot(d_plot_cell_count) + aes(x = sample1_filtered.fastq, y = sample2_filtered.fastq, color = in_ref) +
    geom_point() + scale_y_log10() + scale_x_log10() +
    geom_abline(slope = 1) + geom_vline(xintercept = 1) + geom_hline(yintercept = 1) +
    #     geom_vline(xintercept = 13564 * 0.0001, color = 'blue') + geom_hline(yintercept = 12564 * 0.0001, color = 'blue') +
    theme_bw() + scale_color_manual(values = c("No" = "red", "Yes" = "black"))
# ggsave("./tmp/fig_clustering_norm.pdf", width = 4, height = 3)


#' # Without clustering
d = bc_2dt(bc_obj1)[count > 1]
d[, count_sum := sum(count), by = sample_name]
d[, cell_count := 13564 / 2 / count_sum * count]
d_plot_cell_count = dcast(d, barcode_seq ~ sample_name, sep = "", value.var="cell_count")
d_plot_read_count = dcast(d, barcode_seq ~ sample_name, sep = "", value.var="count")
d_plot_cell_count[is.na(d_plot_cell_count)] = 0
d_plot_read_count[is.na(d_plot_read_count)] = 0

d_plot_cell_count$in_ref = "No"
d_plot_cell_count[barcode_seq %in% ref, in_ref:= "Yes"]

d_plot_read_count$in_ref = "No" 
d_plot_read_count[barcode_seq %in% ref, in_ref:= "Yes"]

ggplot(d_plot_read_count) + aes(x = sample1_filtered.fastq + 1, y = sample2_filtered.fastq + 1, color = in_ref) +
    geom_point() + scale_y_log10() + scale_x_log10() +
    theme_bw() + scale_color_manual(values = c("No" = "red", "Yes" = "black"))

ggplot(d_plot_cell_count) + aes(x = sample1_filtered.fastq, y = sample2_filtered.fastq, color = in_ref) +
    geom_point() + scale_y_log10() + scale_x_log10() +
    geom_abline(slope = 1) + geom_vline(xintercept = 1) + geom_hline(yintercept = 1) +
    #     geom_vline(xintercept = 13564 * 0.0001, color = 'blue') + geom_hline(yintercept = 12564 * 0.0001, color = 'blue') +
    theme_bw() + scale_color_manual(values = c("No" = "red", "Yes" = "black"))


#' # compare to the reads count in original paper
d = bc_2dt(bc_obj1)[count > 1]
d[, count_sum := sum(count), by = sample_name]
d[, cell_count := 13564 / 2 / count_sum * count]
d_plot_cell_count = dcast(d, barcode_seq ~ sample_name, sep = "", value.var="cell_count")
d_plot_read_count = dcast(d, barcode_seq ~ sample_name, sep = "", value.var="count")
d_plot_cell_count[is.na(d_plot_cell_count)] = 0
d_plot_read_count[is.na(d_plot_read_count)] = 0

d_plot_cell_count$in_ref = "No"
d_plot_cell_count[barcode_seq %in% ref, in_ref:= "Yes"]

d_plot_read_count$in_ref = "No" 
d_plot_read_count[barcode_seq %in% ref, in_ref:= "Yes"]
d_plot_read_count[, sample1_filtered.fastq := sample1_filtered.fastq / sum(sample1_filtered.fastq) * 1e5 + 1]
d_plot_read_count[, sample1_filtered.fastq := sample1_filtered.fastq/ sum(sample1_filtered.fastq) * 1e5 + 1]

d_origin = fread("./data/EPO_1000ngml_exp1_all_filtered_Ctrl_1_2_3_4_5_EPO_6_7.txt")
d_origin = d_origin[, .(tag, var.Totalreads.AE3_HIGH_1_M_a, var.Totalreads.AE3_HIGH_1_M_b)]
d_origin = d_origin[var.Totalreads.AE3_HIGH_1_M_a > 0 | var.Totalreads.AE3_HIGH_1_M_b > 0]

d_plot_compare = merge(
    d_plot_read_count[barcode_seq %in%
        d_plot_cell_count[sample1_filtered.fastq> 1 & in_ref == "Yes", barcode_seq]], 
    d_origin, 
    by.x = "barcode_seq", by.y = "tag", all = T)

d_plot_compare[sample1_filtered.fastq > 1]
d_plot_compare[is.na(d_plot_compare)] = 0

d_plot_compare[sample1_filtered.fastq < 1 & var.Totalreads.AE3_HIGH_1_M_a > 100]

g = ggplot(d_plot_compare[sample1_filtered.fastq >= 1 | var.Totalreads.AE3_HIGH_1_M_a > 0]) + 
    aes(x = var.Totalreads.AE3_HIGH_1_M_a + 1, y = sample1_filtered.fastq + 1) +
    geom_point() + scale_y_log10() + scale_x_log10() +
    geom_abline(slope = 1, color = 'red', alpha = 0.5) +
    theme_bw() + labs(x = "Normalized read count (Eisele et al.)", y = "Normmalized read count (CellBarcode)") + ggpubr::stat_cor(method = "spearman")

g

d_plot_compare[sample1_filtered.fastq > 1 & var.Totalreads.AE3_HIGH_1_M_a > 0]
d_plot_compare[sample1_filtered.fastq <= 1 & var.Totalreads.AE3_HIGH_1_M_a > 0]
d_plot_compare[sample1_filtered.fastq > 1 & var.Totalreads.AE3_HIGH_1_M_a <= 0]

#' # compare the reads count to genBaRcode

library(genBaRcode)

bb = "ACGGAATGCTAGAACACTCGAGATCAGNNNNNNNNNNNNNNNNNNNNATGTGGTATGATGTATC"

dir = "../../data/A1007"

time_start = Sys.time()
BC_data <- processingRawData(file_name = "sample1_filtered.fastq",
    source_dir = dir,
    mismatch = 0,
    label = "test",
    bc_backbone = bb,
    bc_backbone_label = "BC_1",
    min_score = 0,
    min_reads = 0,
    save_it = FALSE,
    seqLogo = FALSE,
    cpus = 4,
    strategy = "sequential",
    full_output = FALSE,
    wobble_extraction = TRUE,
    dist_measure = "hamming")

BC_data_EC <- errorCorrection(BC_dat = BC_data,
    maxDist = 1,
    save_it = FALSE,
    cpus = 1,
    strategy = "sequential",
    m = "hamming",
    type = "standard",
    only_EC_BCs = TRUE,
    EC_analysis = FALSE,
    start_small = T)
Sys.time() - time_start

#' Compare to CellBarcode time

time_start = Sys.time()
f = "../../data/A1007/sample1_filtered.fastq"
pattern = "ACGGAATGCTAGAACACTCGAGATCAG(.{20})ATGTGGTATGATGTATC"
bc_obj = bc_extract(f, pattern)
bc_obj1 = bc_cure_depth(bc_obj, depth = 0, isUpdate =F)
bc_obj2 = bc_cure_cluster(bc_obj1, count_threshold = 1e7)
Sys.time() - time_start

#' Compare the result between genBaRcode and CellBarcode

d_g = BC_data_EC@reads %>% data.table
d_c = bc_2dt(bc_obj2) %>% data.table

d_merge = merge(d_g[read_count > 1], d_c[count > 1], by.x = "barcode", by.y = "barcode_seq") %>%
    data.table

d_merge[, count_g := read_count / sum(read_count) * 13564 / 2]
d_merge[, count_c := count / sum(count) * 13564 / 2]
# d_merge = d_merge[barcode %in% ref]

ggplot(d_merge) + aes(y = count_g, x = count_c) + geom_point() + 
    scale_y_log10() + scale_x_log10() + 
    labs(y = "Normalyzed read count (genBaRcode)", x = "Normalized read count (CellBarcode)") +
    geom_abline(slope = 1, color = 'red', alpha = 0.5) +
    geom_hline(yintercept = 1, color = 'red', alpha = 0.5) +
    geom_vline(xintercept = 1, color = 'red', alpha = 0.5) +
    theme_bw() + ggpubr::stat_cor(method = "spearman")

#' Compare the result between Bartender and CellBarcode

d_b = fread("./1M_barcode_cluster.csv")
d_c = bc_2dt(bc_obj2) %>% data.table

d_merge = merge(d_c, d_b, by.x = "barcode_seq", by.y = "Center", all = T) %>%
    data.table

d_merge[is.na(d_merge)] = 0
d_merge[, count_c := count / sum(count) * 13564 / 2]
d_merge[, count_b := time_point_1 / sum(time_point_1) * 13564 / 2]
# d_merge = d_merge[barcode_seq %in% ref]

ggplot(d_merge) + aes(y = count_b, x = count_c) + 
    geom_point()  + theme1 +
    scale_x_log10() + scale_y_log10() + 
    labs(y = "Normalized read count (Bartender)", x = "Normalized read count (CellBarcode)") +
    geom_abline(slope = 1, color = 'red', alpha = 0.5) +
    geom_hline(yintercept = 1, color = 'red', alpha = 0.5) +
    geom_vline(xintercept = 1, color = 'red', alpha = 0.5) +
    theme_bw() + ggpubr::stat_cor(method = "spearman")


# d_merge[count_p.x != 0 & count_p.y != 0] %>% nrow
# d_merge[count_p.x == 0 & count_p.y != 0] %>% nrow
# d_merge[count_p.x != 0 & count_p.y == 0] %>% nrow


#' # Analysis Simulation results


f = c("./tmp/lenti_sumu_filtered.fq", "../../data/A1007/sample1_filtered.fastq")

x = bc_seq_qc(f, sample_name = c("sample2", "sample1"))
bc_plot_seqQc(x)

bc_obj = bc_extract(f, pattern)

x = bc_meta(bc_obj)
x[, 2] / x[, 1]

d_l = bc_obj@messyBc
d_l = d_l %>% lapply(function(x) {
    y = sample(x$barcode_seq, 1e5, prob = x$count, replace = T)
    x = data.table(barcode_seq = names(table(y)), count = as.numeric(table(y)))
    x
})

g_l = lapply(d_l, function(x) {
    ggplot(x, aes(x = log10(count + runif(length(count))))) +
        geom_histogram(bins = 20, fill = "skyblue", color = "black", alpha = 0.7, boundary = 0) +
        geom_rug(linewidth = 0.5, alpha = 0.5) +  
        labs(x = "log10 sequence count") +
        theme_bw() + guides(color = 'none')
})

egg::ggarrange(plots = g_l, nrow = 1)





sessionInfo()

