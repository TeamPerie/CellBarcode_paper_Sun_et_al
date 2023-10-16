library(CellBarcode)
library(ggplot2)
library(ggsci)
library(scales)
library(seewave)
library(data.table)
library(magrittr)
library(stringr)


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

#' # Simulate the data

#+ eval=F

source("../../lib/CellBarcodeSim/lib/lib_simulation.R")
dir.create("./tmp", showWarnings = F)

simulate_main(
    barcode_library_file = "../lib/CellBarcodeSim/example_barcode_library/vdj_barcodes.tsv",
    clone_size_dist      = "lognormal",
    clone_n              = 300,
    clone_size_dist_par  = list(size_mean = 1.2, size_variant = 0.8),
    cycle                = 30,
    efficiency           = 0.705,
    error                = 1e-6,
    pcr_read_per_cell    = 50,
    output_prefix        = paste0("./tmp/small_clone_size_sd_simu"),
    ngs_profile          = "HS20",
    reads_length         = 100,
    is_replicate         = F,
    top_seq              = "AAAAAAAAAAGGGGG",
    bottom_seq           = "ATCGATCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
    sequence_trunk       = 110,
    art_bin              = "../lib/art_bin_MountRainier/art_illumina"
)

fq = c("./tmp/small_clone_size_sd_simu.fq")
sample_name = basename(fq) %>% sub("\\.fq", "", .)
pattern = "AAAAAAAAAAGGGGG(.*)ATCGATCGTTTTTTTTTTTTTTTTTTT"
s_d = bc_extract(fq, pattern = pattern, sample_name = sample_name)@messyBc



#' ## Explore the effect of log transformation

## Without do the log transformation
s_x = s_d[[1]][count > 0, ]
hist((s_x$count / sum(s_x$count) * 1e6 + 1))
rug((s_x$count / sum(s_x$count) * 1e6 + 1) + rnorm(nrow(s_x), 0, 0.1))

s_x = s_d[[1]][count > 0, ]
hist(log10(s_x$count / sum(s_x$count) * 1e6 + 1))
rug(log10(s_x$count / sum(s_x$count) * 1e6 + 1) + rnorm(nrow(s_x), 0, 0.1))


## Simulated data
x1 = rlnorm(10000, 0, 1)
x2 = rlnorm(100, 10, 1)

hist(c(x1, x2))
rug(c(x1, x2))
hist(log10(c(x1, x2)))
rug(log10(c(x1, x2)))


#' ## Explore the effect of removing 1 read sequences


## Without remove the 1 reads case
# e_x = e_d[count > 0, ]
# hist(log10(e_x$count / sum(e_x$count) * 1e6 + 1)) 
# rug(log10(e_x$count / sum(e_x$count) * 1e6 + 1)) #+ rnorm(nrow(e_x), 0, 0.1))
s_x = s_d[[1]][count > 0, ]
hist(log10(s_x$count / sum(s_x$count) * 1e6 + 1))
rug(log10(s_x$count / sum(s_x$count) * 1e6 + 1) + rnorm(nrow(s_x), 0, 0.1))

## Remove the 1 reads case
# e_x = e_d[count > 1, ]
# hist(log10(e_x$count / sum(e_x$count) * 1e6 + 1)) 
# rug(log10(e_x$count / sum(e_x$count) * 1e6 + 1)) #+ rnorm(nrow(e_x), 0, 0.1))
s_x = s_d[[1]][count > 2, ]
hist(log10(s_x$count / sum(s_x$count) * 1e6 + 1))
rug(log10(s_x$count / sum(s_x$count) * 1e6 + 1) + rnorm(nrow(s_x), 0, 0.1))


