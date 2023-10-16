#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

## NOTE: "./tmp/non_umi_simulation/barcode_simulation_id
output_prefix = args[1]

## NOTE: design matrix
design_matrix_file = args[2]

## Load library
library(data.table)
devtools::load_all("../../lib/CellBarcodeSim")

parameter_design_matrix = fread(design_matrix_file)

random_barcode_library = system.file("data", "random_barcodes.tsv", package = "CellBarcodeSim")
vdj_barcode_library = system.file("data", "vdj_barcodes.tsv", package = "CellBarcodeSim")
hamming_barcode_library = system.file("data", "hamming_barcodes.tsv", package = "CellBarcodeSim")

#+ eval=T
for (i in 1:nrow(parameter_design_matrix)) {
#for (i in 1:28) {

    ## define output prefix
    output_prefix_each = paste0(output_prefix, "_simu_seq_", parameter_design_matrix[i, simu_id])

    ## define barcode library file
    if (parameter_design_matrix[i, barcode_type] == "random") {
        barcode_library_file = random_barcode_library
    } else if (parameter_design_matrix[i, barcode_type] == "vdj") {
        barcode_library_file = vdj_barcode_library
    } else if (parameter_design_matrix[i, barcode_type] == "hamming") {
        barcode_library_file = hamming_barcode_library
    } else {
        barcode_library_file = NULL
    }

    ## run simulation
    simulate_main(
        barcode_library_file = barcode_library_file,
        clone_size_dist = parameter_design_matrix[i, clone_size_dist],
        clone_n = parameter_design_matrix[i, clone_n],
        clone_size_dist_par = list(
            size_mean=parameter_design_matrix[i, size_mean],
            size_variant=parameter_design_matrix[i, size_variant]
            ),
        cycle = parameter_design_matrix[i, cycle],
        efficiency = parameter_design_matrix[i, efficiency],
        error = parameter_design_matrix[i, error],
        pcr_read_per_cell = parameter_design_matrix[i, pcr_read_per_cell],
        output_prefix = output_prefix_each,        
        ngs_profile = parameter_design_matrix[i, ngs_profile],
        reads_length         = 100,
        is_replicate         = F,
        top_seq              = "AAAAAAAAAAGGGGG",
        bottom_seq           = "ATCGATCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        sequence_trunk = parameter_design_matrix[i, barcode_length]
    )
}

