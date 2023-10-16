library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)

cycle = 1:60

error_p = c(1e-5, 1e-6, 1e-7)

d = lapply(error_p, function(x){
    base_error = 1 - (1-x)^cycle
    data.table(
        error_rate = x,
        cycle = cycle,
        base_error = base_error
        )
}) %>% rbindlist

ggplot(d) + aes(x = cycle, y = base_error, color = factor(error_rate)) + geom_line() + 
    scale_x_continuous(breaks = seq(0, 60, 5)) + theme_bw() + theme(legend.position = "bottom") +
    labs(color = "PCR Error rate / base / cycle")

ggplot(d) + aes(x = cycle, y = -10 * log10(base_error), color = factor(error_rate)) + geom_line() + 
    scale_x_continuous(breaks = seq(0, 60, 5)) + theme_bw() + theme(legend.position = "bottom") +
    labs(color = "PCR Error rate / base / cycle", y = "Q score")


