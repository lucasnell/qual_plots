# https://statisticalrecipes.blogspot.com/2012/02/simulating-genetic-drift.html
## Simulate Genetic Drift (using Wright-Fisher model)
library(ggplot2)
library(tidyr)
library(dplyr)
library(Rcpp)

source('theme_black.R')


# # Set up parameters
# N <- 20 # number of diploid individuals
# p <- 0.2
# q <- 1 - p
# N_gen <- 100 # number of generations
# N_sim <- 5 # number of simulations


R_sim <- function(N, p, N_gen, N_sim, diploid = TRUE){
    
    if (diploid) {
        N_chrom <- 2*N
    } else {
        N_chrom <- N
    }
    
    if (length(N) == 1) {
        N_chrom <- rep(N_chrom, N_gen)
    } else if (!length(N) == N_gen) {
        stop('N should be of length 1 or N_gen.')
    }
    # Adding an NA at the beginning so it's the same length as p_mat below
    N_chrom <- c(NA, N_chrom)
    p_mat <- matrix(0, nrow = N_gen + 1, ncol = N_sim)
    # initialize p in first generation
    p_mat[1,] <- rep(p, N_sim)
    for(s in 1:N_sim){
        for(t in 2:(N_gen+1)){
            indivs <- rbinom(1, N_chrom[t], prob = p_mat[t-1,s])
            p_mat[t,s] <- indivs / N_chrom[t]
        }
    }
    p_df <- p_mat[2:(N_gen+1),] %>% 
        as_data_frame %>%
        rename_(.dots = setNames(paste0('V', 1:N_sim), paste0('sim', 1:N_sim))) %>% 
        mutate(gen = 1:N_gen, N_chrom = N_chrom[2:(N_gen+1)]) %>%
        gather_('sim', 'p', paste0('sim', 1:N_sim))
    return(p_df)
}





# cppFunction("NumericVector cpprbinom(int n, double size, NumericVector prob) {
#     NumericVector v(n);
#     for (int i=0; i<n; i++) {v[i] = as<double>(rbinom(1, size, prob[i]));}
#     return(v); }")



set.seed(2); R_sim(200, 0.2, 100, 5) %>% 
    ggplot(aes(x = gen, y = p, color = sim)) +
    geom_line() +
    # ggtitle("Simulations of Genetic Drift") +
    xlab("Generation") +
    ylab("Allele frequency") +
    ylim(0,1) +
    theme_black() +
    theme(legend.position = 'none')




# =========================
# =========================
# With linearly decreasing N
# =========================
# =========================

N_vec <- as.integer(sapply(0:99, function(x) 1000 - 5*x))
N_int <- round(1 / (sum(1 / N_vec) * (1/length(N_vec))), 0)

set.seed(9)
sim_data <- list(
    R_sim(N_vec, 0.2, 100, 100) %>%
        mutate(type = 'decreasing'),
    R_sim(N_int, 0.2, 100, 100) %>%
        mutate(type = 'steady')) %>%
    bind_rows

ggplot(sim_data, aes(x = gen, y = p, group = interaction(sim, type), color = type)) +
    geom_line(alpha = 0.4) +
    ggtitle("Simulations of Genetic Drift") +
    xlab("Generation") +
    ylab("Allele Frequency") +
    ylim(0,1) +
    theme_classic()


sim_data %>%
    filter(gen == max(gen)) %>%
    ggplot(aes(type, p)) +
    theme_classic() +
    geom_jitter(alpha = 0.4) +
    stat_summary(fun.data = "mean_cl_boot", size = 1)





# =========================
# =========================
# Now with random autocorrelative time series
# =========================
# =========================


set.seed(1)
a1 <- 0.5
a2 <- 0.2
N_vec <- numeric(100)
for (i in 3:length(N_vec)){
    N_vec[i] <- a1 * N_vec[i-1] + a2 * N_vec[i-2] + rnorm(1, sd = 100)
}
N_vec <- round(N_vec,0) + 1000
# plot(seq_along(N_vec), N_vec, type = 'l')
N_int <- round(1 / (sum(1 / N_vec) * (1/length(N_vec))), 0)

set.seed(9)
sim_data <- list(
    R_sim(N_vec, 0.2, 100, 100) %>%
        mutate(type = 'AR2'),
    R_sim(N_int, 0.2, 100, 100) %>%
        mutate(type = 'steady')) %>%
    bind_rows

ggplot(sim_data, aes(x = gen, y = p, group = interaction(sim, type), color = type)) +
    geom_line(alpha = 0.4) +
    xlab("Generation") +
    ylab("Allele Frequency") +
    ylim(0,1) +
    theme_classic()


sim_data %>%
    filter(gen == max(gen)) %>%
    ggplot(aes(type, p)) +
    theme_classic() +
    geom_jitter(alpha = 0.4) +
    stat_summary(fun.data = "mean_cl_boot", size = 1)







# =========================
# =========================
# Now with cyclical time series
# =========================
# =========================


N_vec <- sapply(0:99, function(x) 900 * cos({x/pi} * 2) + 1000)
N_vec <- round(N_vec,0)
# plot(N_vec, type = 'l')

N_int <- round(1 / (sum(1 / N_vec) * (1/length(N_vec))), 0)




set.seed(9)
sim_data <- list(
    R_sim(N_vec, 0.2, 100, 100) %>%
        mutate(type = 'cyclical'),
    R_sim(N_int, 0.2, 100, 100) %>%
        mutate(type = 'steady')) %>%
    bind_rows

ggplot(sim_data, aes(x = gen, y = p, color = type)) +
    geom_line(aes(group = sim), alpha = 0.4) +
    xlab("Generation") +
    ylab("Allele Frequency") +
    ylim(0,1) +
    theme_classic() +
    facet_grid(type ~ .) +
    theme(legend.position = 'none')

sim_data %>% 
    filter(type == 'cyclical') %>% 
    ggplot(aes(x = gen, y = p)) +
    # geom_line(data = data.frame(gen = 1:100, p = N_vec / max(N_vec)), 
    #           color = 'black') +
    geom_line(aes(group = sim), color = 'dodgerblue', alpha = 0.4) +
    xlab("Generation") +
    ylab("Allele Frequency") +
    ylim(0,1) +
    theme_classic()

    


sim_data %>%
    filter(gen == max(gen)) %>%
    ggplot(aes(type, p)) +
    theme_classic() +
    geom_jitter(alpha = 0.4) +
    stat_summary(fun.data = "mean_cl_boot", size = 1)






