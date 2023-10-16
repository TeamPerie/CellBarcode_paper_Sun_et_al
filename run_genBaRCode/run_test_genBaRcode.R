#' # Test genBaRcode

library(genBaRcode)

#' ## Input files
b_file = "genBaRcode_test.fq"
bb = "AAAAAAAAAAGGGGGNNNNNNNNNNNNNNATCGATCGTTTTTTT"

#' ## Extract barcodes by genBaRcode
t_init = Sys.time()
BC_data <- processingRawData(file_name = b_file,
    source_dir = "./",
    mismatch = 0,
    label = "test",
    bc_backbone = bb,
    bc_backbone_label = "BC_1",
    min_score = 0,
    min_reads = 0,
    save_it = FALSE,
    seqLogo = FALSE,
    cpus = 1,
    strategy = "sequential",
    full_output = FALSE,
    wobble_extraction = TRUE,
    dist_measure = "hamming")
Sys.time() - t_init
gb_d_raw = BC_data@reads %>% data.table

#' ## Error correction by genBaRcode
t_init = Sys.time()
BC_data_EC <- errorCorrection(BC_dat = BC_data,
    maxDist = 1,
    save_it = FALSE,
    cpus = 1,
    strategy = "sequential",
    m = "hamming",
    type = "standard",
    only_EC_BCs = TRUE,
    EC_analysis = FALSE,
    start_small = TRUE)
Sys.time() - t_init

gb_d_correct = BC_data_EC@reads

#' # Test CellBarcode

#' ## Input files
library(CellBarcode)
b_file = "genBaRcode_test.fq"
bb = "AAAAAAAAAAGGGGG([ATCG]{14})ATCGATCGTTTTTTT"

#' ## Extract barcodes by CellBarcode

t_init = Sys.time()
x = bc_extract(as.list(b_file), bb)
Sys.time() - t_init

cb_d_raw = x@messyBc[[1]]

#' ## Test error correction by CellBarcode
t_init = Sys.time()
cb_d_correct = bc_cure_cluster(bc_cure_depth(x))@cleanBc[[1]]
Sys.time() - t_init


#' # Compare genBaRcode and CellBarcode results

d = merge(gb_d_raw, cb_d_raw, by.x = "barcode", by.y = "barcode_seq", all = T)
d[is.na(d)] = 0
plot(d$read_count + 1, d$count + 1, log= "xy", xlab = "genBaRcode", ylab = "CellBarcode")
sum(d$read_count != d$count)
d[read_count != count]

ref_d = fread("genBaRcode_ref.tsv")
d = merge(gb_d_correct, cb_d_correct, by.x = "barcode", by.y = "barcode_seq", all = T)
d = merge(d, ref_d, by.x="barcode", by.y = "seq", all=T)
d[is.na(d)] = 0
plot(d$read_count + 1, d$count + 1, log= "xy", xlab = "genBaRcode", ylab = "CellBarcode")
with(d, cor(count, read_count))
with(d, cor(count, freq))
with(d, cor(read_count, freq))


