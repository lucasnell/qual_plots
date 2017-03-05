# This one was never used. See `wolbachia.R` file for the one I used in the meeting

# From the following paper (Table 2):
# Hilgenboecker, K., P. Hammerstein, P. Schlattmann, A. Telschow, and J. H. Werren. 2008.
#   How many species are infected with Wolbachia? – a statistical analysis of current
#   data. FEMS Microbiology Letters 281:215–220.

library(tidyverse)

df <- data_frame(
    n = c("1", "2", "10", "≥10", ">100"),
    s = c(547, 110, 6, 115, 13),
    inf = c(19, 21, 33, 54, 92)) %>%
    mutate(n = factor(n, levels = c("1", "2", "10", "≥10", ">100")))

ggplot(df, aes(n, inf)) +
    theme_classic() +
    geom_bar(stat = 'identity', color = NA, fill = 'dodgerblue') +
    ylab('% species infected') +
    xlab('Sample size per species') +
    geom_text(aes(label = s, y = inf + 0.5), vjust = 0, size = 3)# +
    # labs(caption = "(Hilgenboecker et al. 2008)")
