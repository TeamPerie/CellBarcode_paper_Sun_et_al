library(CellBarcode)
library(ggpubr)
library(data.table)
library(ggplot2)
library(magrittr)
library(readr)
library(stringr)

ref = data.frame(
    barcode_seq = c(
        "AAGTCCAGTTCTACTATCGTAGCTACTA",	
        "AAGTCCAGTATCGTTACGCTACTA",	
        "AAGTCCAGTCTACTATCGTTACGACAGCTACTA",	
        "AAGTCCAGTTCTACTATCGTTACGAGCTACTA",	
        "AAGTCCATCGTAGCTACTA",	
        "AAGTCCAGTACTGTAGCTACTA",	
        "AAGTCCAGTACTATCGTACTA"	
    ), clone_size = 2^c(1, 7, 3, 4, 5, 6, 2)
)

#' # The log-sd will not change with the scale of the data

ref$clone_size %>% log2() %>% sd
(ref$clone_size * 1000) %>% log2() %>% sd

dir("../../data/5290/", "mef", full = T)
fq_files = dir("../../data/5290/", "mef", full = T)[c(1, 3, 2, 10)]

metadata <- stringr::str_split_fixed(basename(fq_files), "_", 10)[, c(4, 6)]
metadata <- as.data.frame(metadata)
sample_name <- apply(metadata, 1, paste, collapse = "_")
colnames(metadata) = c("cell_number", "replication")
rownames(metadata) = sample_name

#+ eval=T, echo = T
bc_obj <- bc_extract(
  fq_files,  # fastq file
  pattern = "([ACGT]{12})CTCGAGGTCATCGAAGTATC([ACGT]+)CCGTAGCAAGCTCGAGAGTAGACCTACT", 
  pattern_type = c("UMI" = 1, "barcode" = 2),
  sample_name = sample_name,
  metadata = metadata
)
bc_obj

#' ## Check if UMI is unique or not
x = bc_obj@messyBc[[1]]
x[umi_seq == "GCGTGTTACCTC" & count > 100, ]

#' # Unique UMI

bc_obj = bc_cure_umi(bc_obj, depth = 10, isUniqueUMI=T)
bc_plot_mutual(bc_obj, highlight = ref$barcode_seq)
bc_obj = bc_cure_depth(bc_obj, depth = 0)
# bc_obj = bc_cure_cluster(bc_obj, count_threshold = 1e7)

bc_obj@metadata

d = bc_2dt(bc_obj)[sample_name %in% c("195_mixa", "195_mixb"), ]
d = d[, .(count = sum(count)), by = barcode_seq]
d[, count_sum := sum(count)]
d = d[count > count_sum * 0.0001]
d[, count_sum := sum(count)]
d[, cell_count := 195 / count_sum * count]
d[is.na(d)] = 0
names(d) %<>% make.names

d$in_ref = F
d[barcode_seq %in% ref$barcode_seq, in_ref := T]
d = merge(d, ref, by = "barcode_seq", all = T)
d[is.na(d)] = 0

ggplot(d) + aes(x = clone_size, y = cell_count, color = in_ref) +
    geom_point() + scale_y_log10() + scale_x_log10() +
    geom_smooth(method = "lm", col = 'black') +
    stat_regline_equation(col = "black") +
    theme_bw() + scale_color_manual(values = c("black", "red")) 

#' # Non-unique UMI
bc_obj = bc_cure_umi(bc_obj, depth = 10, isUniqueUMI=F)
bc_plot_mutual(bc_obj, highlight = ref$barcode_seq)
bc_obj = bc_cure_depth(bc_obj, depth = 0)
# bc_obj = bc_cure_cluster(bc_obj, count_threshold = 1e7)

bc_obj@metadata

d = bc_2dt(bc_obj)[sample_name %in% c("195_mixa", "195_mixb"), ]
d = d[, .(count = sum(count)), by = barcode_seq]
d[, count_sum := sum(count)]
d = d[count > count_sum * 0.0001]
d[, count_sum := sum(count)]
d[, cell_count := 195 / count_sum * count]
# d = dcast(d, barcode_seq ~ sample_name, sep = "", value.var="cell_count")
d[is.na(d)] = 0
names(d) %<>% make.names

d$in_ref = F
d[barcode_seq %in% ref$barcode_seq, in_ref := T]
d = merge(d, ref, by = "barcode_seq", all = T)
d[is.na(d)] = 0

ggplot(d) + aes(x = clone_size, y = cell_count, color = in_ref) +
    geom_point() + scale_y_log10() + scale_x_log10() +
    geom_smooth(method = "lm", col = 'black') +
    stat_regline_equation(col = "black") +
    theme_bw() + scale_color_manual(values = c("black", "red")) 

#' # Do not use UMI

bc_obj = bc_cure_depth(bc_obj, depth = 0, isUpdate = F)
bc_obj = bc_cure_depth(bc_obj, depth = -1)

# d = bc_2dt(bc_obj)[sample_name %in% c("195_mixa", "195_mixb"), ]
d = d[, .(count = sum(count)), by = barcode_seq]
d[, count_sum := sum(count)]
d = d[count > count_sum * 0.0001]
d[, count_sum := sum(count)]
d[, cell_count := 195 / count_sum * count]
# d = dcast(d, barcode_seq ~ sample_name, sep = "", value.var="cell_count")
d[is.na(d)] = 0
names(d) %<>% make.names

d$in_ref = F
d[barcode_seq %in% ref$barcode_seq, in_ref := T]
d = merge(d, ref, by = "barcode_seq", all = T)
d[is.na(d)] = 0

ggplot(d) + aes(x = clone_size * 2^8, y = cell_count, color = in_ref) +
    geom_point() + scale_y_log10() + scale_x_log10() +
    geom_smooth(method = "lm", col = 'black') +
    stat_regline_equation(col = "black") +
    theme_bw() + scale_color_manual(values = c("black", "red")) 


