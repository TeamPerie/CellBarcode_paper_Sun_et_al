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
    line = element_line(size = 1),
    axis.line = element_line(size = 1),
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

#' # Read in data

u_f1 = "../../data/VDJ_mammry_gland_meghan_UMI/UMI_m9250_P0_Lum_Rep1/L345R45.R1.fastq.gz"
u_f2 = "../../data/VDJ_mammry_gland_meghan_UMI/UMI_m9250_P0_Lum_Rep2/L345R46.R1.fastq.gz"

u_f3 = "../../data/VDJ_mammry_gland_meghan_UMI/noUMI_m9250_P0_Lum_Rep1/L412R10.R1.fastq.gz"
u_f4 = "../../data/VDJ_mammry_gland_meghan_UMI/noUMI_m9250_P0_Lum_Rep2/L412R11.R1.fastq.gz"


f_v1 = c(u_f1, u_f2)
f_v2 = c(u_f3, u_f4)

bc_qc = bc_seq_qc(c(f_v1[2], f_v2[2]))
bc_plot_seqQc(bc_qc)
bc_plot_seqQc(bc_qc[1])
bc_plot_seqQc(bc_qc[2])

bc_qc = bc_seq_qc(f_v2)
bc_plot_seqQc(bc_qc)


pattern1 = "(.{16})CCTCGAGGTCATCGAAGTATCAAG(.*)CCGTAGCAAGCTCGAGAGTAGACCTACT"
pattern2 = "CCTCGAGGTCATCGAAGTATCAAG(.*)CCGTAGCAAGCTCGAGAGTAGACCTACT"

time_start = Sys.time()
bc_obj_umi = bc_extract(f_v1, pattern1, sample_name = c("rep1", "rep2"), pattern_type = c(barcode = 2, UMI = 1))
Sys.time() - time_start

time_start = Sys.time()
bc_obj_noumi = bc_extract(f_v2, pattern2, sample_name = c("rep1", "rep2"))
Sys.time() - time_start


#' # Unique UMI

## UMI reads sensitivity test
d_plot = ldply(1:10, function(x){
    y = bc_cure_umi(bc_obj_umi, depth = x, isUniqueUMI = T, doFish= F)
    y = bc_cure_depth(y, depth = 0)
    n = bc_2dt(y)[sample_name == "rep2"][count >= 1] %>% nrow
    c(cutoff = x, n = n)
})
knitr::kable(d_plot)

#' ## Test UMI reads cutoff

d_plot = ldply(1:10, function(x){
    y = bc_cure_umi(bc_obj_umi, depth = x, isUniqueUMI = T, doFish= F)
    y = bc_cure_depth(y, depth = 0)
    n = bc_2dt(y)[sample_name == "rep1"][count >= 1] %>% nrow
    c(cutoff = x, n = n)
})

knitr::kable(d_plot)

ggplot(d_plot) + aes(x = cutoff, y = n) +
    geom_bar(stat = "identity") + theme1

#' ## UMI filtering

time_start = Sys.time()
bc_obj_umi = bc_cure_umi(bc_obj_umi, depth = 10, isUniqueUMI = T, doFish= F)
Sys.time() - time_start
bc_plot_mutual(bc_obj_umi, count_marks = 2)
bc_obj_umi = bc_cure_depth(bc_obj_umi, depth = 0)
bc_2df(bc_obj_umi)

#' ## No UMI filtering

#' ### No filtering

time_start = Sys.time()
bc_obj_noumi = bc_cure_depth(bc_obj_noumi, depth = 0, isUpdate = F)
Sys.time() - time_start
cutoff_x = bc_auto_cutoff(bc_obj_noumi)
bc_plot_mutual(bc_obj_noumi, count_marks = cutoff_x)

#' ### Auto filtering

bc_obj_noumi = bc_cure_depth(bc_obj_noumi, depth = -1, isUpdate = F)
bc_plot_mutual(bc_obj_noumi)

bc_obj_noumi_d = bc_2dt(bc_obj_noumi)
bc_obj_umi_d = bc_2dt(bc_obj_umi)

#' ## Compare UMI and no UMI results

d_merge = merge(bc_obj_noumi_d, bc_obj_umi_d, all = T, by = c("sample_name", "barcode_seq"))
d_plot = d_merge[sample_name == "rep2"]
d_plot[is.na(d_plot)] = 0
d_plot[count.x == 0 | count.y == 0]

## Y is the umi count
t_y = sum(d_plot$count.y) * 0.0001

d_plot[count.y < t_y, count.y := 0]
d_plot = d_plot[count.x > 0 | count.y > 0]

d_plot$count.y = d_plot$count.y / sum(d_plot$count.y)
d_plot$count.x = d_plot$count.x / sum(d_plot$count.x)
ggplot(d_plot) + aes(x = count.x, y = count.y) + 
    geom_point()  + theme1 +
    labs(x = "Reads Count", y = "UMI Count") +
    geom_smooth(method = "lm") + stat_regline_equation()

d_plot[count.x != 0 & count.y != 0]
d_plot[count.x == 0 ]
d_plot[count.y == 0 ]

#' # clone size variance
(d_plot$count.x * 10000 + 1) %>% log %>% sd

#' # Compare Bartender

#' ## Bartender UMI result v.s. CellBarcode UMI result

d_b = fread("./1M_barcode_cluster.csv")
d_b = d_b[, count_p := time_point_1 / sum(time_point_1)]
d_c = bc_obj_umi_d[sample_name == "rep2"]
d_c = d_c[, count_p := count / sum(count)]

d_merge = merge(d_c, d_b, by.x = "barcode_seq", by.y = "Center", all = T)
d_merge[is.na(d_merge)] = 0

ggplot(d_merge) + aes(x = count_p.x, y = count_p.y) + 
    geom_point()  + theme1 +
    scale_x_log10() + scale_y_log10() + 
    geom_smooth(method = "lm") + stat_regline_equation() +
    labs(x = "CellBarcode UMI Count%", y = "Bartender Count%")

d_merge[count_p.x != 0 & count_p.y != 0] %>% nrow
d_merge[count_p.x == 0 & count_p.y != 0] %>% nrow
d_merge[count_p.x != 0 & count_p.y == 0] %>% nrow

#' ## Bartender non-UMI result v.s. CellBarcode no UMI result

d_b = fread("./1M_barcode_non_umi_cluster.csv")
d_b = d_b[, count_p := time_point_1 / sum(time_point_1)]
d_c = bc_obj_noumi_d[sample_name == "rep2"]
d_c = d_c[, count_p := count / sum(count)]

d_merge = merge(d_c, d_b, by.x = "barcode_seq", by.y = "Center", all = T)
d_merge[is.na(d_merge)] = 0

ggplot(d_merge) + aes(x = count_p.x, y = count_p.y) + 
    geom_point()  + theme1 +
    labs(x = "UMI Count", y = "Bartender Count") + 
    scale_x_log10() + scale_y_log10() + 
    geom_smooth(method = "lm") + stat_regline_equation() +
    labs(x = "CellBarcode Reads Count%", y = "Bartender Count%")

d_merge[count_p.x != 0 & count_p.y != 0] %>% nrow
d_merge[count_p.x == 0 & count_p.y != 0] %>% nrow
d_merge[count_p.x != 0 & count_p.y == 0] %>% nrow


#' ## Bartender UMI result v.s. Bartender non-UMI result

d_b_u = fread("./1M_barcode_cluster.csv")
d_b_u = d_b_u[, count_p := time_point_1 / sum(time_point_1)]
d_b_nu = fread("./1M_barcode_non_umi_cluster.csv")
d_b_nu = d_b_nu[, count_p := time_point_1 / sum(time_point_1)]

d_merge = merge(d_b_u, d_b_nu, by.x = "Center", by.y = "Center", all = T)

ggplot(d_merge) + aes(x = count_p.x, y = count_p.y) + 
    geom_point()  + theme1 +
    labs(x = "UMI Count", y = "Bartender Count") + 
    scale_x_log10() + scale_y_log10() + 
    geom_smooth(method = "lm") + stat_regline_equation() +
    labs(x = "UMI Count%", y = "reads Count%")

sessionInfo()
