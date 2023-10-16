library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(CellBarcode)
library(ggrepel)
library(ROCR)
library(ggsci)
library(ggpubr)

theme0 <- theme_bw() + theme(
    text = element_text(size = 15),
    line = element_line(size = 1),
    axis.line = element_line(size = 1),
    axis.ticks.length = unit(3, units = "mm"),
    axis.text.x = element_text(
        margin = margin(t = 2, unit = "mm")
        , angle = 60, vjust = 1, size = 15, hjust = 1),
    axis.text.y = element_text(margin = margin(r = 3, l = 5, unit = "mm")),
    legend.position = "right",
) 

theme0_lt = theme0 + theme(legend.position = "top")
theme1 = theme0 + theme(
    axis.text.x = element_text( 
        margin = margin(t = 2, unit = "mm")
        , angle = 0, vjust = 1, size = 12, hjust = 0.5)
)

#' # Sample info

sample_info = fread("../run_simulation_umi/umi_simu_design_matrix.tsv")
sample_info$simu_id %<>% as.character
setkey(sample_info, "simu_id")

i_level_simu_name = sample_info[order(as.integer(simu_id)), simu_name][c(1, 6)]

#' # Clone size disbribution
true_barcode_l = read_rds("../run_preprocess_simulation/tmp/umi_simulatino_ref.rds")
d_true_barcode = true_barcode_l %>% rbindlist(idcol = "simu_file") 
x = d_true_barcode$simu_file %>% str_match("barcode_simulation_\\d+_simu_seq_(.*)") %>% extract(, c(2))
table(x, useNA = "always")
d_true_barcode = cbind(d_true_barcode, simu_id = x)
# d_true_barcode[simu_name == "hamming" & grepl("simulation_19_", simu_file), .N, by = barcode_seq][N > 1]
# d_true_barcode[simu_name == "vdj_barcode" & grepl("simulation_1_", simu_file), .N, by = barcode_seq][N > 1]
# d_true_barcode[simu_name == "vdj_barcode"]
d_true_barcode$simu_name = sample_info[d_true_barcode$simu_id, simu_name]
table(d_true_barcode$simu_name, useNA = "always")
d_true_barcode = rename(d_true_barcode, c("seq" = "barcode_seq"))

## Barcode frequency
ggplot(d_true_barcode[, .(freq = sum(freq)), by = .(simu_file, barcode_seq, simu_name)]) + aes(x = simu_name, y = freq) +
    geom_boxplot() +
    # geom_jitter(alpha = 0.2) + 
    theme0 + scale_y_log10() +
    scale_x_discrete(limits = i_level_simu_name)


#' # AUC
# NOTE: input
d_auc = fread("./tmp/table_auc_umi.tsv")
d_count_true_barcode = fread("tmp/table_count_true_barcode_umi.tsv")

#' # count - true barcode
ggplot(d_count_true_barcode) + aes(x = simu_name, y = count, color = factor(is_true)) +
    geom_boxplot(alpha = 0.2) +
    theme0_lt + scale_y_log10() + scale_color_npg() +
    scale_x_discrete(limits = i_level_simu_name)


#' # AUC versus simu conditions
ggplot(d_auc) + aes(x = simu_name, y = auc) + geom_boxplot() + theme0 +
    scale_x_discrete(limits = i_level_simu_name) +
    labs(y = "AUC") + lims(y = c(0, 1))

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01,
          0.05, Inf), symbols = c("****", "***", "**",
          "*", "ns"))
ggplot(d_auc) + aes(x = simu_name, y = aucpr) + geom_boxplot() + theme0 +
    scale_x_discrete(limits = i_level_simu_name) +
    labs(y = "P-R AUC") + lims(y = c(0, 1)) + 
    stat_compare_means(method = "wilcox.test", na.rm = T, hide.ns = T, symnum.args = symnum.args)
ggsave("./tmp/fig_aucpr_umi.pdf")


d = fread('./tmp/table_count_true_barcode_umi.tsv')
d[simu_id == 6][is_true == 0]

