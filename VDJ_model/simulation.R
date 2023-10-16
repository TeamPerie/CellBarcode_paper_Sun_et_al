library(magrittr)
library(data.table)
library(stringr)
library(DNAscar)
library(parallel)
library(readr)

mo = load_model("final_parms.txt")
ma = load_marginal("final_marginals.txt")


y = simu_sequence(ma, mo, 1e7, freq=T, seed=0)
d = data.table(
    seq = names(y) %>% sub("CGAAGTATCAAG", "", .) %>% sub("CCGTAGCAAGCTCGAGAGTAGACCTACT", "", .),
    freq = as.integer(y))
write_tsv(d, "./vdj_barcodes.tsv")
