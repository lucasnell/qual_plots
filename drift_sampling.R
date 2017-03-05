# This simulates drift for a diploid mating organism, where mating pairs are explicitly
# simulated


library(ggplot2)
library(tidyr)
library(dplyr)
library(Rcpp)



# ===========================================================================
# ===========================================================================

# Functions

# ===========================================================================
# ===========================================================================


# -----
# A single simulation
# ------

one_sim <- function(N_gen, N_vec, geno_ps, genos) {
    
    p_vec <- numeric(N_gen + 1)
    # Initialize individual genotypes in first generation
    geno_nums <- round(geno_ps * N_vec[1], 0)
    # If rounding makes it not add up to N_vec, then add 1 randomly to one 
    # of the end points
    if (sum(geno_nums) != N_vec[1]){
        rand_i <- sample(1:3, 1)
        geno_nums[rand_i] <- geno_nums[rand_i] + (N_vec[1] - sum(geno_nums))
    }
    # This is the matrix for the previous generation's genotypes
    genos_t0 <- lapply(1:3, function(j) matrix(rep(genos[j,], geno_nums[j]), ncol = 2, 
                                               byrow = TRUE)) %>% 
        do.call(what = rbind, args = .)
    # Shuffle genotypes by row
    genos_t0 <- genos_t0[sample.int(nrow(genos_t0)),]
    # Assign to a matrix by sex
    F_n <- round(N_vec[1]/2, 0)
    F_t0 <- genos_t0[1:F_n,]
    M_t0 <- genos_t0[(F_n+1):N_vec[1],]
    
    for (t in 1:N_gen){
        
        p_vec[t] <- sum(genos_t0) / prod(dim(genos_t0))
        
        repro_F <- sample.int(nrow(F_t0), N_vec[t], TRUE)
        repro_M <- sample.int(nrow(M_t0), N_vec[t], TRUE)
        alleles_F <- sample.int(2, length(repro_F), TRUE)
        alleles_M <- sample.int(2, length(repro_M), TRUE)
        genos_t0 <- matrix(0, nrow = N_vec[t], ncol = 2)
        genos_t0[,1] <- F_t0[cbind(repro_F, alleles_F)]
        genos_t0[,2] <- M_t0[cbind(repro_M, alleles_M)]
        
        # Assign to a matrix by sex
        F_n <- round(N_vec[t]/2, 0)
        F_t0 <- genos_t0[1:F_n,]
        M_t0 <- genos_t0[(F_n+1):N_vec[t],]
        
    }
    p_vec[(N_gen+1)] <- sum(genos_t0) / prod(dim(genos_t0))
    return(p_vec)
}





# -----
# Combining multiple simulations
# ------

R_sim_diploid <- function(N, p, N_gen, N_sim){
    
    if (length(N) == 1) {
        N_vec <- rep(N, N_gen)
    } else if (length(N) != N_gen) {
        stop('N should be of length 1 or N_gen.')
    }
    
    # Assuming HWE
    geno_ps <- c(p^2, 2 * p * (1-p), (1-p)^2)
    genos <- rbind(c(1,1), c(1,0), c(0,0))
    
    p_mat <- sapply(1:N_sim, function(i) one_sim(N_gen, N_vec, geno_ps, genos))
    
    p_df <- p_mat[2:(N_gen+1),] %>% 
        as_data_frame %>%
        rename_(.dots = setNames(paste0('V', 1:N_sim), paste0('sim', 1:N_sim))) %>% 
        mutate(gen = 1:N_gen, N = N_vec) %>%
        gather_('sim', 'p', paste0('sim', 1:N_sim))
    
    return(p_df)
}







# ===========================================================================
# ===========================================================================

# With linearly decreasing N

# ===========================================================================
# ===========================================================================


N_vec <- as.integer(sapply(0:99, function(x) 1000 - 5*x))
N_int <- round(1 / (sum(1 / N_vec) * (1/length(N_vec))), 0)

set.seed(9)
sim_data <- list(
    R_sim_diploid(N_vec, 0.2, 100, 100) %>%
        mutate(type = 'decreasing'),
    R_sim_diploid(N_int, 0.2, 100, 100) %>%
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
    N_vec[i] <- a1 * N_vec[i-1] + a2 * N_vec[i-2] + rnorm(1, sd = 200)
}
N_vec <- round(N_vec,0) + 1000
# plot(N_vec, type = 'l')
N_int <- round(1 / (sum(1 / N_vec) * (1/length(N_vec))), 0)

set.seed(9)
sim_data <- list(
    R_sim_diploid(N_vec, 0.2, 100, 100) %>%
        mutate(type = 'AR2'),
    R_sim_diploid(N_int, 0.2, 100, 100) %>%
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


N_vec <- sapply(0:99, function(x) 900 * cos({x/pi} * 0.5) + 1000)
N_vec <- round(N_vec,0)
plot(N_vec, type = 'l')

N_int <- round(1 / (sum(1 / N_vec) * (1/length(N_vec))), 0)


set.seed(9)
sim_data <- list(
    R_sim_diploid(N_vec, 0.5, 100, 50) %>%
        mutate(type = 'cyclical'),
    R_sim_diploid(N_int, 0.5, 100, 50) %>%
        mutate(type = 'steady')) %>%
    bind_rows

ggplot(sim_data, aes(x = gen, y = p, color = type)) +
    geom_line(aes(group = interaction(sim, type)), alpha = 0.4) +
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
