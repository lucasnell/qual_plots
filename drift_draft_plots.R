# For plots of possible scenarios of genetic drift and draft


suppressPackageStartupMessages({
    library(tidyverse)
})

ar1 <- 0.4
ar2 <- 0.4
df <- data_frame(drift = numeric(100), draft = numeric(100), both = numeric(100),
                 gen = 1:100)
set.seed(1)
for (t in 3:100){
    df$drift[t] <- ar2 * df$drift[(t-2)] + ar1 * df$drift[(t-1)] + rnorm(1)
    df$draft[t] <- ar2 * df$draft[(t-2)] + ar1 * df$draft[(t-1)] + rnorm(1)
    df$both[t] <- ar2 * df$both[(t-2)] + ar1 * df$both[(t-1)] + rnorm(1)
}

df <- df %>% 
    mutate(drift = drift + abs(min(df$drift)) + 1,
           both = both + abs(min(df$both)) + 1, 
           draft = draft + abs(min(df$draft)) + 10) %>% 
    gather(type, N, drift:both)


cols <- c('#66c2a5','#fc8d62','#8da0cb')

df %>% 
    ggplot(aes(gen, N, color = type, linetype = type)) +
    geom_line(size = 1) +
    geom_text(data = df %>% filter(gen == 55), aes(label = type),
              fontface = 'bold', vjust = 0, hjust = 1,
              nudge_y = c(2, 3, -1), nudge_x = c(-0.75, 2, 12)) +
    ylab('Population size') +
    xlab('Time') +
    theme_classic() +
    theme(axis.ticks = element_blank(), axis.text = element_blank(), 
          legend.position = 'none') +
    scale_color_manual(values = c('purple', 'red', 'dodgerblue'))



shift <- 2

df2 <- data_frame(x = seq(-3, 3 + shift, length.out = 100),
                  y = sapply(x, function(z) dnorm(z)),
                  y2 = sapply(x, function(z) dnorm(z - shift, sd = 0.25))) %>% 
    mutate(y2 = (y2 / diff(range(y2))) * (diff(range(y))))



df2 %>% 
    ggplot(aes(x, y)) +
    geom_line(size = 1) +
    theme_classic() +
    ylab('Frequency') +
    scale_x_continuous('Host attribute', limits = c(-3, 3 + shift)) +
    theme(axis.ticks = element_blank(), axis.text = element_blank())

df2 %>% 
    ggplot(aes(x, y)) +
    geom_line(size = 1, linetype = 2) +
    geom_line(aes(y = y2), size = 1) +
    theme_classic() +
    ylab('Frequency') +
    scale_x_continuous('Host attribute', limits = c(-3, 3 + shift)) +
    theme(axis.ticks = element_blank(), axis.text = element_blank())

curve(dnorm(x), -3, 3)















