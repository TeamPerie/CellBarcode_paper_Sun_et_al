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

sample_info = fread("../run_simulation_no_umi/non_umi_simu_design_matrix.tsv")
sample_info$simu_id %<>% as.character
setkey(sample_info, "simu_id")

i_level_simu_name = sample_info[order(as.integer(simu_id)), simu_name][c(7, 1, 8, 23, 22, 24)]

#' # Clone size disbribution
# true_barcode_l = read_rds("../tmp/non_umi_simulatino_ref_mainfigure.rds")
true_barcode_l = read_rds("../run_preprocess_simulation/tmp/non_umi_simulation_ref.rds")
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
d_auc = fread("./tmp/table_auc_no_umi_depth.tsv")
d_count_true_barcode = fread("./tmp/table_count_true_barcode_depth.tsv")

#' # count - true barcode
ggplot(d_count_true_barcode) + aes(x = simu_name, y = count, color = factor(is_true)) +
    geom_boxplot(alpha = 0.2) +
    theme0_lt + scale_y_log10() + scale_color_npg() +
    scale_x_discrete(limits = i_level_simu_name)


#' # AUC versus simu conditions
ggplot(d_auc) + aes(x = simu_name, y = auc) + geom_boxplot() + theme0 +
    scale_x_discrete(limits = i_level_simu_name) +
    labs(y = "AUC") + lims(y = c(0, 1))

comparisons_l = list(
    c("size_variation_1", "base"),
    c("size_variation_1", "size_variation_3"),
    c("vdj_barcode_size_variation_1", "vdj_barcode"),
    c("vdj_barcode_size_variation_1", "vdj_barcode_size_variation_3")
    )

d_auc_sub = d_auc[simu_name %in% i_level_simu_name]
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01,
          0.05, Inf), symbols = c("****", "***", "**",
          "*", "ns"))
ggplot(d_auc_sub) + aes(x = simu_name, y = aucpr) + geom_boxplot() + theme0 +
    scale_x_discrete(limits = i_level_simu_name) +
    labs(y = "P-R AUC") + lims(y = c(0, 1.5)) + 
    stat_compare_means(method = "wilcox.test", comparisons = comparisons_l, na.rm = T, hide.ns = T, symnum.args = symnum.args)
ggsave("./tmp/fig_aucpr_no_umi_depth.pdf")

#' # Resolusion index
d_count_true_barcode = fread("tmp/table_count_true_barcode_depth.tsv")

# d_count_true_barcode[, count_cat := cut(clone_size, breaks = c(0, 1, 2, 3, 4, 5, 10, 20, Inf))]
# lapply(1:100, function(i) {
# d_count_true_barcode[, sum(count < i & is_true), ,by = .(clone_size, simu_file)]
# })


#' ## Apply the autothreshold strategy
# NOTE: input
d = fread("./tmp/table_no_umi_auto_threshold_predict_rate_depth.tsv")

# ggplot(d) + aes(x = tpr, y = fpr, label = simu_name) + geom_point() +
# geom_label_repel(alpha = 0.5)
d_plot_pr = melt(d, id.vars = "simu_name", measure.vars=c("pre", "rec"), variable.name = "pr", value.name = "value") 
d_plot_tf = melt(d, id.vars = "simu_name", measure.vars=c("tpr", "fpr"), variable.name = "pr", value.name = "value") 

d_plot_pr = d_plot_pr[simu_name %in% i_level_simu_name]
ggplot(d_plot_pr) + aes(x = simu_name, y = value, fill = factor(pr)) + 
    geom_boxplot() + theme0 +
    scale_x_discrete(limits = i_level_simu_name)

#+ fig.height=6
ggplot(d_plot_pr) + aes(x = simu_name, y = value) + 
    geom_boxplot() +
    theme0 + facet_grid(pr ~ .) +
    scale_x_discrete(limits = i_level_simu_name) +
    stat_compare_means(method = "wilcox.test", comparisons = comparisons_l, na.rm = T, hide.ns = T, symnum.args = symnum.args)
ggsave("./tmp/fig_pr_no_umi_depth.pdf")

ggplot(d_plot_tf) + aes(x = simu_name, y = value) + 
    geom_boxplot() +
    theme0 + facet_grid(pr ~ .) +
    scale_x_discrete(limits = i_level_simu_name)

#' ## Apply the threshold of 0.00001 of left reads
# NOTE: input
d = fread("./tmp/table_no_umi_p001_threshold_predict_rate_depth.tsv")

# ggplot(d) + aes(x = tpr, y = fpr, label = simu_name) + geom_point() +
# geom_label_repel(alpha = 0.5)
d_plot_pr = melt(d, id.vars = "simu_name", measure.vars=c("pre", "rec"), variable.name = "pr", value.name = "value") 
d_plot_tf = melt(d, id.vars = "simu_name", measure.vars=c("tpr", "fpr"), variable.name = "pr", value.name = "value") 

ggplot(d_plot_pr) + aes(x = simu_name, y = value, fill = factor(pr)) + 
    geom_boxplot() + theme0 +
    scale_x_discrete(limits = i_level_simu_name)

#+ fig.height=6
d_plot_pr = d_plot_pr[simu_name %in% i_level_simu_name]
ggplot(d_plot_pr) + aes(x = simu_name, y = value) + 
    geom_boxplot() +
    theme0 + facet_grid(pr ~ .) +
    scale_x_discrete(limits = i_level_simu_name) +
    stat_compare_means(method = "wilcox.test", comparisons = comparisons_l, na.rm = T, hide.ns = T, symnum.args = symnum.args)
ggsave("./tmp/fig_pr_no_umi_depth_more.pdf")



ggplot(d_plot_tf) + aes(x = simu_name, y = value) + 
    geom_boxplot() +
    theme0 + facet_grid(pr ~ .) +
    scale_x_discrete(limits = i_level_simu_name)

#' ## Merge threshold
# d1 = fread("./tmp/table_no_umi_auto_threshold_predict_rate_depth.tsv")[, resolution := "T auto"]
# d2 = fread("./tmp/table_no_umi_p01_threshold_predict_rate_depth.tsv")[, resolution := "T 0.01"]
d3 = fread("./tmp/table_no_umi_p001_threshold_predict_rate_depth.tsv")[, resolution := "T 0.001"]
d4 = fread("./tmp/table_no_umi_p0001_threshold_predict_rate_depth.tsv")[, resolution := "T 0.0001"]

d = rbindlist(list(d3, d4))
d = melt(d, id.vars = c("simu_name", "resolution"), measure.vars=c("pre", "rec"), variable.name = "pr", value.name = "value") 
d$resolution %<>% factor(levels = c("T auto", "T 0.01", "T 0.001", "T 0.0001"))

d = d[simu_name %in% i_level_simu_name]
ggplot(d) + aes(x = simu_name, y = value, color = resolution) +
    geom_boxplot() + 
    theme0 + facet_grid(pr ~ .) + scale_color_npg() +
    scale_x_discrete(limits = i_level_simu_name) + 
    theme(strip.background = element_blank()) +
    stat_compare_means(method = "wilcox.test", comparisons = comparisons_l, na.rm = T, hide.ns = T, symnum.args = symnum.args)
ggsave("./tmp/fig_pr_no_umi_depth_more.pdf")




