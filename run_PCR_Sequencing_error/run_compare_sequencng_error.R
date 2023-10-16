library(data.table)
library(ggplot2)
library(plyr)
library(magrittr)
library(readr)
library(stringr)
library(gmp)

knitr::opts_chunk$set(
    fig.height = 3, fig.width = 5, fig.align = "center"
)

#' function to read a sequence profile file into a data.table

calculate_median_from_frequency <- function(class, frequency) {
  # Create a data frame
  df <- data.frame(class = class, frequency = frequency)
  
  # Calculate cumulative frequency
  df$cumulative_frequency <- cumsum(df$frequency)
  
  # Find n/2
  n_half <- sum(df$frequency) / 2
  
  # Identify the median class
  median_class <- subset(df, cumulative_frequency >= n_half)[1, ]
  
  # Determine lower boundary of median class
  lower_boundary <- median_class$class - 0.5
  
  # Interpolate to find the median
  median <- lower_boundary + ((n_half - (median_class$cumulative_frequency - median_class$frequency)) / median_class$frequency) * diff(range(median_class$class))
  
  return(median)
}

parse_quality_file = function(file) {
    l = readLines(file)
    q = l[seq(1, length(l), 2)]
    n = l[seq(2, length(l), 2)]

    q_s = str_split(q, "\t")
    n = str_split(n, "\t")

    res = lapply(seq_along(q_s), function(i) {
        q_v = q_s[[i]]
        n_v = n[[i]]
        base = q_v[1]
        cycle = q_v[2]
        data.table(
            base = base,
            cycle = as.integer(cycle),
            q = as.numeric(q_v[-c(1, 2)]),
            n = as.numeric(n_v[-c(1, 2)])
        )
    }) %>% rbindlist

    res[q != ""]

    # d_plot = res[, .(
    #     q_mean = sum(q * n) / sum(n),
    #     q_sd = sqrt(sum(q^2 * n) / sum(n) - (sum(q * n) / sum(n))^2)
    #     ), by = .(base, cycle)]
    #
    # d_plot
}

plot_art_profile = function(f) {
    res = parse_quality_file(f)
    d_plot = res[, .(
        q_mean = sum(q * n) / sum(n),
        q_median = calculate_median_from_frequency(q, n),
        q_sd = sqrt(sum(q^2 * n) / sum(n) - (sum(q * n) / sum(n))^2)
        ), by = .(base, cycle)]
    g = ggplot(d_plot[cycle <= 100 & base == "."]) + aes(x = cycle, y = q_median) +
        geom_line() +
        ylim(0, 45) + theme_bw() + labs(y = "Median Q score")
    g
}

#' # MiSeq V1 250R1

f = "./EmpMiSeq250R1.txt"
plot_art_profile(f) + labs(title = "MiSeq")

#' # HiSeq2000 100R1

f = "./HiSeq2kL100R1.txt"
plot_art_profile(f) + labs(title = "HiSeq")

#' # NovaSeq Adil

f = "./Nova.txt"
plot_art_profile(f) + labs(title = "NovaSeq Adil")

#' # MiSeq Adil

f = "./MiSeq.txt"
plot_art_profile(f) + labs(title = "MiSeq Adil")


