x = fread("./tmp/GSE99915_reprogramming.timecourse.dge_barcodes.tsv", head=F)[[1]]

y = str_split_fixed(x, "_", 3) %>% data.table

z = y[V2 == "HF-1" & grepl("-6", V3)][, V3]
z %<>% sub("-6", "-1", .)
write_lines(z, "./tmp/hf1.d15.barcode.txt")
