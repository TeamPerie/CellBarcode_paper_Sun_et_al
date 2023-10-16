library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(dplyr)

#' # Simulation the range of clone size in log normal distribution

#' Simulation parameter range
lmean = seq(0.5, 3, 0.1)
lsd = c(1, 2, 3)

df = expand.grid(lmean = lmean, lsd = lsd)
d = data.table(df)

#' Calculate the range of clone size
d[, range := qlnorm(0.99, lmean, lsd) / qlnorm(0.01, lmean, lsd) , by = .(lmean, lsd)]


#' Plot the range of clone size
ggplot(d) + aes(x = lmean, y = range, color = factor(lsd)) + geom_smooth(se = TRUE) +
    scale_y_log10() + labs(y = "log normal Quantile 99% / Quantile 1%", x = "log mean", color = "log sd")


x = rlnorm(1000, 1, 1) %>% sort
sum(x[1:50]) / sum(x[51:1000])
