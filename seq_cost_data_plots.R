
# Sequence costs and GenBank deposits through time

suppressPackageStartupMessages({
    library(tidyverse)
})

source('theme_black.R')


seq_cost <- read_csv('seqcost2015_4.csv') %>% 
    filter(!is.na(cost_mb))

moores_law <- function(d) {
    t <- as.numeric(d - seq_cost$date[1], units = 'days')
    seq_cost$cost_mb[1] * 0.5^(t / (365 * 2))
}


ggplot(seq_cost, aes(date, cost_mb)) +
    scale_x_date(NULL, limits = c(min(seq_cost$date), max(genbank$date))) +
    scale_y_continuous('Cost per Mb of DNA', breaks = 10^(-1:4), 
                       labels = c('$0.1', '$1', '$10', '$100', '$1k', '$10k')) +
    stat_function(fun = moores_law, linetype = 3, color = 'white') +
    geom_line(color = 'dodgerblue', size = 1) +
    coord_trans(y = 'log10', limy = c(0.01, 10e3)) +
    annotate(geom = 'text', label = "Moore's Law", hjust = 0.5, vjust = 0,
             x = median(seq_cost$date) + 365, y = 750, color = 'white') +
    theme_black()



genbank <- read_csv('genbank_wgs_stats.csv', col_types = 'iDdd')

ggplot(genbank, aes(date, bases / 1e9)) +
    scale_x_date(NULL, limits = c(min(seq_cost$date), max(genbank$date))) +
    scale_y_log10('Gb on GenBank', breaks = 10^(0:3), position = 'right') +
    geom_line(color = 'red', size = 1) +
    theme_black() +
    theme(panel.background = element_blank(), plot.background = element_blank())

