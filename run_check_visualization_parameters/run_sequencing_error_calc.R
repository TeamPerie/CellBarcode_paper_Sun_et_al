library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

#' # Sequencing error

d = fread("../../lib/art_bin_MountRainier/Illumina_profiles/HiSeq2kL100R1.txt", skip = 27)

x_score = unlist(d[V2 == 15][1])[-c(1,2)] %>% as.numeric
count = unlist(d[V2 == 15][2])[-c(1,2)] %>% as.numeric

sum(count * 10^(x_score / (-10))) / sum(count)

#' # PCR error

#' Given the PCR cycle of 30 and PCR error of 10e-6

1 - (1 - 10e-6) ^ 30
