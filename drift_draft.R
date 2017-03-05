#' ---
#' title: "Drift and draft"
#' author: "Lucas Nell"
#' date: "`r format(Sys.Date())`"
#' output: github_document
#' ---
#' 

suppressPackageStartupMessages({
    library(ggplot2)
    library(tidyr)
    library(dplyr)
    # library(Rcpp)
})


source('theme_black.R')

# Using "Genetic hitchhiking" by Barton (doi 10.1098/rstb.2000.0716), p 1154
#     for how sweeps change u through time
# And "Population genetics: a concise guide" by Gillespie, p 54
#     for change of p due to selection






sim_sweep <- function(N, r, s, N_sim, h = 0.5, w_bar = 0.7, u0 = 0.25, N_gen = 100,
                      diploid = TRUE){
    
    if (diploid) {
        N_chrom <- 2*N
    } else {
        N_chrom <- N
    }
    
    if (!length(N) == 1) {
        stop('N should be of length 1.')
    }
    
    sim_inds <- 1:N_sim + 2
    
    u_inc <- (1 - u0) * N_chrom^(-r / s)
    
    p_mat <- matrix(0, nrow = N_gen + 1, ncol = N_sim + 2)
    p_mat[1,] <- c(1/N_chrom, rep(u0, N_sim+1))
    
    u <- u0
    p <- 1/N_chrom
    q <- 1 - p
    for(t in 2:(N_gen+1)){
        delta_p <- {p * q * s * (p * h + q * {1 - h})} / w_bar
        p <- p + delta_p
        q <- 1 - p
        delta_u <- u_inc * p
        # {(r * (1 - u0)) / (s * p)} * exp(-r * t)
        u <- u0 + delta_u
        indivs_U <- rbinom(1, N_chrom, prob = u)
        # u <- min(c(1, u + (UP_start - UQ_start) * exp(-r * t)))
        p_mat[t,1] <- p  # indivs_P / N_chrom
        p_mat[t,2] <- indivs_U / N_chrom
        p_mat[t, sim_inds] <- sapply(
            sim_inds, function(si) rbinom(1, N_chrom, prob = p_mat[t-1,si]) / N_chrom)
    }
    
    p_df <- p_mat[2:(N_gen+1),] %>%
        as_data_frame %>%
        rename_(.dots = setNames(paste0('V', 1:ncol(p_mat)), paste0(1:ncol(p_mat)))) %>% 
        mutate(gen = 1:N_gen, N_chrom = N_chrom) %>%
        gather_('sim', 'p', paste0(1:ncol(p_mat))) %>% 
        mutate(type = factor(c(rep('sweep', N_gen), rep('draft', N_gen), 
                        rep('drift', N_gen * N_sim)), 
                        levels = c('drift', 'sweep', 'draft')),
               sim = as.integer(sim))
    
    return(p_df)
}



# N = 200; r = 0.01; s = 0.15; N_sim = 5; h = 0.5; w_bar = 0.7; u0 = 0.25; N_gen = 100
# diploid = TRUE
# rm(N, r, s, N_sim, h, w_bar, u0, N_gen, diploid)


set.seed(2)
sw_df <- sim_sweep(N = 200, r = 0.02, s = 0.2, N_sim = 5)


cols <- c('#66c2a5','#fc8d62','#8da0cb')


sw_df %>%
    filter(type == 'drift') %>% 
    ggplot(aes(x = gen, y = p, color = type)) +
    geom_line(aes(group = sim)) +
    geom_text(data = sw_df %>% filter(gen == 75, type == 'drift') %>% 
                  group_by(gen, type) %>% summarize(p = max(p)), 
              aes(label = type), 
              nudge_y = 0.1, nudge_x = 0, vjust = 0, hjust = 1) +
    xlab("Generation") +
    ylab("Allele frequency") +
    ylim(0,1) +
    theme_black() +
    theme(legend.position = 'none') +
    scale_color_manual(values = cols[1])



sw_df %>%
    filter(type != 'draft') %>% 
    ggplot(aes(x = gen, y = p, color = type)) +
    geom_line(aes(group = sim)) +
    geom_text(data = sw_df %>% filter(gen == 75, type != 'draft') %>% 
                  group_by(gen, type) %>% summarize(p = max(p)), 
              aes(label = type), 
              nudge_y = c(0.1,0), nudge_x = c(0,-12.5), vjust = 0, hjust = 1) +
    xlab("Generation") +
    ylab("Allele frequency") +
    ylim(0,1) +
    theme_black() +
    theme(legend.position = 'none') +
    scale_color_manual(values = cols[1:2])



sw_df %>%
    ggplot(aes(x = gen, y = p, color = type)) +
    geom_line(aes(group = sim)) +
    geom_text(data = sw_df %>% filter(gen == 75) %>% group_by(gen, type) %>% 
                  summarize(p = max(p)), aes(label = type), 
              nudge_y = c(0.1, 0, 0.05), nudge_x = c(0,-12.5,0), vjust = 0, hjust = 1) +
    xlab("Generation") +
    ylab("Allele frequency") +
    ylim(0,1) +
    theme_black() +
    theme(legend.position = 'none') +
    scale_color_manual(values = cols)
