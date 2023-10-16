library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(CellBarcode)

## create the output dir
dir.create("./tmp", showWarnings = FALSE)

#' Read in all the files 
## for each simulation, find the files
## read in the data
## add annotation
## output the data

simu_dir = "../run_simulation_umi/tmp/"

fq_file_v = dir(simu_dir, "fq", full=T)
ref_file_v = dir(simu_dir, "ref", full=T)

fq_name = basename(fq_file_v) %>% sub(".fq", "", .)
ref_name = basename(ref_file_v) %>% sub("_ref.tsv", "", .)

ref_l = lapply(ref_file_v, fread)
names(ref_l) = ref_name

sample_info = fread("../run_simulation_umi/umi_simu_design_matrix.tsv")
setkey(sample_info, "simu_id")
simu_id = fq_file_v %>% basename %>% str_match("(\\d+).fq") %>% extract(, 2) %>% as.integer
i_umi8 = sample_info[umi_length == 8, simu_id]
i_umi4 = sample_info[umi_length == 4, simu_id]
i_umi12 = sample_info[umi_length == 12, simu_id]

## umi length of 8
pattern = "^ATCGATCG(.{8})(.*)ATCGATCGTTTT"
fq_name[i_umi8]
bc_obj8 = bc_extract(fq_file_v[i_umi8], pattern = pattern, pattern_type = c(UMI=1, barcode=2), sample_name = fq_name[i_umi8])


## umi length of 12
pattern = "^ATCGATCG(.{12})(.*)ATCGATCGTTTT"
fq_name[i_umi12]
bc_obj_umi12 = bc_extract(as.list(fq_file_v[i_umi12]), pattern = pattern, pattern_type = c(UMI=1, barcode=2), sample_name = fq_name[i_umi12])


## umi length of 4
pattern = "^ATCGATCG(.{4})(.*)ATCGATCGTTTT"
fq_name[i_umi4]
bc_obj_umi4 = bc_extract(as.list(fq_file_v[i_umi4]), pattern = pattern, pattern_type = c(UMI=1, barcode=2), sample_name = fq_name[i_umi4])


write_rds(bc_obj8, "./tmp/umi_simulation_bcobj_umi8.rds")
write_rds(bc_obj_umi12, "./tmp/umi_simulation_bcobj_umi12.rds")
write_rds(bc_obj_umi4, "./tmp/umi_simulation_bcobj_umi4.rds")
write_rds(ref_l, "./tmp/umi_simulatino_ref.rds")

# bc_obj = bc_cure_umi(bc_obj, depth = 5, isUniqueUMI = F, doFish=F)
# bc_obj = bc_cure_depth(bc_obj, depth = 0)
# bc_obj
# 
# x = bc_2dt(bc_obj[, 1])
# x[order(count)]
# y = ref_l[[9]][, .(freq = sum(freq)), by = seq]
# 
# res = merge(x, y, by.x = "barcode_seq", by.y = "seq", all = T)
# res[is.na(freq)]
# res[is.na(count)]
# 

# z = read_lines("./tmp/umi_simulation/barcode_simulation_1_simu_seq_1_library.fasta")
# zz = res[is.na(freq)]$barcode_seq
# res_zz = mclapply(zz, function(i) {
#     print(i)
#     grep(i, z) %>% length
# })

# res[is.na(res)] = 0
# with(res, plot(x = count + 1, y = freq + 1, xlab = "UMI count + 1", ylab = "Clone size + 1", log = "xy"))
# qith(res, cor.test(x = count + 1, y = freq + 1, xlab = "UMI count + 1", ylab = "Clone size + 1", log = "xy"))
#  
# res[is.na(count)]
# res[is.na(freq)]
# 
# bc_2dt(bc_obj[, 1])[barcode_seq %in% ref_l[[1]]$seq]
# 
# (x %in% ref_l[[1]]$seq)  %>% sum
# ref_l[[1]][!(seq %in% x)][, .(freq = sum(freq)), by = seq]
# ref_l[[1]][(seq %in% x)][, .(freq = sum(freq)), by = seq][order(freq, decreasing=T)][1:50]
# ref_l[[1]]$seq %>% unique %>% length
# ref_l[[5]]
# 
# ref_l = lapply(ref_file_v, fread)
# names(ref_l) = ref_name
# ref_l[c(1, 6)]
# 
# write_rds(bc_obj, "./tmp/non_umi_simulation_bcobj.rds")
