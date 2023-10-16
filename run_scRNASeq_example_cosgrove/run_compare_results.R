## Jason Final result
x_jason = fread("./data/top_bc_per_cell_VDJ_barcodes_m534_52_10_both.csv")
x_jason = x_jason[, .(cell_barcode = Cell.bc %>% sub("-1", "", .), barcode_seq = Con.seq %>% sub("CCG$", "", .), JS_UMI_count = Numi.filt)]

## Analysis Anne-Marie's result
x = fread("./data/agrep_10xbc_and_vbc_m534_both.txt", header = F)
y = x[, .N, by = .(V1, V2, V3)]
names(y) = c("cell_barcode", "umi", "barcode_seq", "count")
z = y[count >= 3][, .(barcode_seq = barcode_seq[which.max(count)]), by = .(umi, cell_barcode)]
o = z[, .N, by = .(cell_barcode, barcode_seq)][N >= 1]
o$cell_barcode %<>% sub("-1", "", .)
o$barcode_seq %<>% sub("CCG$", "", .)

o[, .N, by = cell_barcode][N > 1]

## compare with my result
x_my = fread("./data/KDI_2015730_2023-01-31_14-58-15/L438T02.R1.fastq.gz.tsv")
x_my = rename(x_my, c("sample_name" = "cell_barcode"))

x_merge = merge(x_my, o, by = c("cell_barcode", "barcode_seq"), all.x = T, all.y = T)
names(x_merge) = c("cell_barcode", "barcode_seq", "WJ_UMI_count", "AM_UMI_count")
x_merge = merge(x_merge, x_jason, all = T)

x_merge[is.na(x_merge)] = 0
x_merge[JS_UMI_count > 0]
x_merge[WJ_UMI_count > 0]
x_b1 = x_merge[WJ_UMI_count == 0, barcode_seq] %>% unique
x_b2 = x_merge[WJ_UMI_count > 0, barcode_seq] %>% unique
sum(x_b1 %in% x_b2)
sum(x_b2 %in% x_b1)

library(ggplot2)
library(ggpubr)
ggplot(x_merge, aes(x = AM_UMI_count, y = WJ_UMI_count)) + geom_jitter(size = 0.5) + geom_abline(intercept = 0, slope = 1) + xlab("Anne-Marie's result") + ylab("Wenjie's result")
ggplot(x_merge, aes(x = AM_UMI_count, y = JS_UMI_count)) + geom_jitter(size = 0.5) + geom_abline(intercept = 0, slope = 1) + xlab("Anne-Marie's result") + ylab("Jason's result")

ggplot(x_merge, aes(x = WJ_UMI_count, y = JS_UMI_count)) + geom_jitter(size = 0.5) + geom_abline(intercept = 0, slope = 1) + xlab("UMI count (CellBarcode)") + ylab("UMI count (Cosgrove et al.)") +
    stat_cor(method = "spearman", label.round = 2, size = 3)

x_merge[WJ_UMI_count > 0, length(unique(barcode_seq))]
x_merge[AM_UMI_count > 0, length(unique(barcode_seq))]
x_merge[JS_UMI_count > 0, length(unique(barcode_seq))]

x_merge %>% write_tsv("./tmp/compare_result.tsv")
x_merge[, .N, by = cell_barcode]
x_merge
