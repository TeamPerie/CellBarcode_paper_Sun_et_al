#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

## NOTE: "./tmp/non_umi_simulation/barcode_simulation_id
output_prefix = args[1]

## NOTE: design matrix
design_matrix_file = args[2]

library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(CellBarcode)
library(ggrepel)

source("lib/lib_simulation.R")

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


##################################
#  Parameters of the simulation  #
##################################

#' ```
#' Run the simulation
#'   parameters:
#'     barcode related:
#'       barcode_type: hamming, random, vdj, custom
#'       barcode_library_file: default NULL
#'     clone size related:
#'       clone_size_dist: uniform, lognormal
#'       clone_n: number of clone
#'       clone_size_dist_par:
#'         for lognormal: list(size_mean=1.2, size_variant=2)
#'         for uniform: list(size_max=1000, size_min=1)
#'     pcr related:
#'       cycle: default 30
#'       efficiency: 0.705
#'       error: 1e-6
#'       pcr_read_per_cell: 50
#'     sequence related:
#'       output_prefix: "./seq"
#'       ngs_profile: MSv1, HS20
#'       reads_length: default 110
#' ```

#' ## Simulate Random barcode

#' ```
#' Base condition:
#'   Barcode type: Random barcodes
#'   Clone size:
#'     distribution: lognormal
#'     n: 300
#'     size_mean:  1.2
#'     size_variant:  2
#'   PCR:
#'     cycle: 30
#'     efficiency: 0.705
#'     error: 1e-6
#'     pcr_read_per_cell: 50
#'   sequence:
#'     ngs_profile: HS20
#' Variant condition:
#'   Barcode type: Random barcodes
#'   Clone size:
#'     distribution: lognormal
#'     n: 150, 600, 1200
#'     size_mean:  0.6, 3.0, 6
#'     size_variant:  1, 4, 8
#'   PCR:
#'     cycle: 20, 35
#'     efficiency: 0.5, 0.9
#'     error: 1e-5, 1e7
#'     pcr_read_per_cell: 25, 100
#'   sequence:
#'     ngs_profile: MSv1
#' ```

parameter_design_matrix = fread(design_matrix_file)

#+ eval=T
for (i in 1:nrow(parameter_design_matrix)) {
    output_prefix_each = paste0(output_prefix, "_simu_seq_", parameter_design_matrix[i, simu_id])
    simulate_main_umi(
        barcode_type = parameter_design_matrix[i, barcode_type],
        barcode_library_file = NULL,
        clone_size_dist = parameter_design_matrix[i, clone_size_dist],
        clone_n = parameter_design_matrix[i, clone_n],
        clone_size_dist_par = list(
            size_mean=parameter_design_matrix[i, size_mean],
            size_variant=parameter_design_matrix[i, size_variant]
            ),
        cycle = parameter_design_matrix[i, cycle],
        efficiency = parameter_design_matrix[i, efficiency],
        error = parameter_design_matrix[i, error],
        pcr_read_per_umi = parameter_design_matrix[i, pcr_read_per_umi],
        output_prefix = output_prefix_each,
        ngs_profile = parameter_design_matrix[i, ngs_profile],
        reads_length = 100,
        umi_length = parameter_design_matrix[i, umi_length],
        preamp_n = parameter_design_matrix[i, preamp_cycle],
        umi_tagging_efficiency = parameter_design_matrix[i, umi_tagging_efficiency],
        top_seq = "ATCGATCG",
        bottom_seq = "ATCGATCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        sequence_trunk = parameter_design_matrix[i, barcode_length]
    )
}

