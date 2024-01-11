library(CellBarcode)
library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(ggpubr)

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

pattern = "ACGGAATGCTAGAACACTCGAGATCAG(.{20})ATGTGGTATGATGTATC"

f = c("../../data/A1007/sample1_filtered.fastq", "../../data/A1007/sample2_filtered.fastq")
bc_obj = bc_extract(f, pattern)

d = bc_obj@messyBc[[1]]

## Original distribution
ggplot(d, aes(x = log10(count))) + geom_histogram(bins = 15) 
## Remove 1 reads barcode
ggplot(d[count != 1], aes(x = log10(count))) + geom_histogram(bins = 15) 


## Sampling 1e5 reads
d2 = sample(d$barcode_seq, 1e5, replace = T, prob = d$count)
d2 = data.table(barcode_seq = names(table(d2)), count = as.numeric(table(d2)))
ggplot(d2, aes(x = log10(count))) + geom_histogram(bins = 15) 

## Sampling 1e4 reads
d2 = sample(d$barcode_seq, 1e4, replace = T, prob = d$count)
d2 = data.table(barcode_seq = names(table(d2)), count = as.numeric(table(d2)))
ggplot(d2, aes(x = log10(count))) + geom_histogram(bins = 15) 

