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

simu_dir = "../run_simulation_no_umi/tmp/"

fq_file_v = dir(simu_dir, "fq", full=T)
basename(fq_file_v) %>% str_match(., "_simulation_(.*)_simu") %>% extract(, 2) %>% table
ref_file_v = dir(simu_dir, "ref", full=T)

fq_name = basename(fq_file_v) %>% sub(".fq", "", .)
ref_name = basename(ref_file_v) %>% sub("_ref.tsv", "", .)

bc_obj = bc_extract(fq_file_v, pattern = "^AAAAAAAAAAGGGGG(.*)ATCGATCGTTTTTTT", sample_name = fq_name)

ref_l = lapply(ref_file_v, fread)
names(ref_l) = ref_name

write_rds(bc_obj, "./tmp/non_umi_simulation_bcobj.rds")
write_rds(ref_l, "./tmp/non_umi_simulation_ref.rds")

#' Read in part file for main figures

simu_dir = "../run_simulation_no_umi/tmp/"

fq_file_v = dir(simu_dir, "fq", full=T) %>% grep("seq_(1|7|8|22|23|24).fq", ., value=T)
ref_file_v = dir(simu_dir, "ref", full=T) %>% grep("seq_(1|7|8|22|23|24)_ref.tsv", ., value=T)

fq_name = basename(fq_file_v) %>% sub(".fq", "", .)
ref_name = basename(ref_file_v) %>% sub("_ref.tsv", "", .)

bc_obj = bc_extract(fq_file_v, pattern = "^AAAAAAAAAAGGGGG(.*)ATCGATCGTTTTTTT", sample_name = fq_name)

ref_l = lapply(ref_file_v, fread)
names(ref_l) = ref_name

write_rds(bc_obj, "./tmp/non_umi_simulation_bcobj_mainfigure.rds")
write_rds(ref_l, "./tmp/non_umi_simulation_ref_mainfigure.rds")

