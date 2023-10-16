library(data.table)
library(stringdist)
library(devtools)

load_all("../../lib/CellBarcodeSim/")

random_barcode_ref = system.file("data", "random_barcodes.tsv", package = "CellBarcodeSim")
hamming_barcode_ref = system.file("data", "hamming_barcodes.tsv", package = "CellBarcodeSim")
vdj_barcode_ref = system.file("data", "vdj_barcodes.tsv", package = "CellBarcodeSim")

vdj_ref = fread(vdj_barcode_ref)
ran_ref = fread(random_barcode_ref)

x = vdj_ref[order(freq, decreasing = T)][1:500]$seq
y = stringdistmatrix(x, x, method = "hamming")
hist(y, main = "Histgram of hamming distance between top 500 VDJ barcode")

x = ran_ref[order(freq, decreasing = T)][1:500]$seq
y = stringdistmatrix(x, x, method = "hamming")
hist(y, main = "Histgram of hamming distance between top 500 Random barcode")
