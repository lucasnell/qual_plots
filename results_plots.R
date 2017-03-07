# Plots of potential results



suppressPackageStartupMessages({
    library(tidyverse)
})



#' <!-- ======================================================================== -->
#' 
#' <!-- Attributes plot -->
#' 
#' <!-- ======================================================================== -->
#' 
#' 
#' Attributes:
#' 
#' 1. Dispersal
#' 2. Population growth
#' 3. Survival at high population densities
#' 

set.seed(1)
spp <- sample(LETTERS[1:10], 10)

set.seed(2)
sp_df <- as_data_frame(rbind(
    # 3 specialists
    c(1,runif(2, 0.05, 0.3)),
    c(runif(1, 0.05, 0.3),1,runif(1, 0.05, 0.3)),
    c(runif(2, 0.05, 0.3),1),
    # 2 generalists
    runif(3, 0.6, 0.8),
    runif(3, 0.6, 0.8),
    # Others
    runif(3, 0.05, 0.8),
    runif(3, 0.05, 0.8),
    runif(3, 0.05, 0.8),
    runif(3, 0.05, 0.8),
    runif(3, 0.05, 0.8)
    )) %>% 
    rename_(.dots = setNames(paste0('V', 1:3), c('disp', 'grow', 'surv'))) %>% 
    mutate(species = factor(spp)) %>% 
    arrange(species) %>% 
    gather('attr', 'score', disp:surv)


# Looking at means (this is not the final plot)
sp_df %>% ggplot(aes(attr, score, color = attr)) +
    theme_bw() +
    theme(strip.background = element_rect(color = NA, fill = NA), 
          strip.text = element_text(face = 'bold'),
          legend.position = 'none', panel.grid = element_blank()) +
    geom_point(size = 2) +
    facet_wrap(~ species, nrow = 2) +
    xlab('Attribute') +
    ylab('Score')


sp_df %>% filter(score == 1)
# Specialists are C, D, and E
# Generalists are B and G







indiv_fun <- function(sp, sim_sd = 0.05, sim_N = 10) {
    attrs <- sp_df %>% filter(species == sp) %>% arrange(attr)
    
    ran_mat <- matrix(rnorm(sim_N * 3, 0, sim_sd), ncol = 3) + 
        matrix(rep(attrs$score, sim_N), ncol = 3, byrow = TRUE)
    
    for (i in 1:3) {
        while (!isTRUE(all.equal(mean(ran_mat[,i]), attrs$score[i]))){
            ran_mat[,i] <- ran_mat[,i] - mean(ran_mat[,i]) + attrs$score[i]
        }
    }
    
    # If it has something <= 0, change it to a really small value > 0
    ran_mat <- ifelse(ran_mat > 0, ran_mat, 0.001)

    ran_df <- as_data_frame(ran_mat)
    colnames(ran_df) <- c('disp', 'grow', 'surv')
    
    return(mutate(ran_df, species = sp))
}


indiv_df <- lapply(unique(sp_df$species), indiv_fun, sim_N = 100, sim_sd = 0.1) %>% 
    bind_rows %>% 
    arrange(species) %>% 
    gather('attr', 'score', disp:surv)


attr_p <- ggplot(sp_df, aes(attr, score, color = attr)) +
    theme_bw() +
    theme(strip.background = element_rect(color = NA, fill = NA), 
          strip.text = element_text(face = 'bold'),
          legend.position = 'none', panel.grid = element_blank()) +
    geom_jitter(data = indiv_df, alpha = 0.3, size = 1, shape = 16) +
    geom_point(size = 2, color = 'black') +
    facet_wrap(~ species, nrow = 2) +
    xlab('Attribute') +
    ylab('Score')
attr_p


# ggsave('~/Desktop/attr_p.pdf', attr_p, width = 6, height = 3, units = 'in')





#' <!-- ======================================================================== -->
#' 
#' <!-- Abundances plot -->
#' 
#' <!-- ======================================================================== -->
#' 
#' 
#' Scenarios:
#' 
#' 1. Specialists - 3 spp
#' 2. Generalist - 1 or more spp
#' 3. Mixture - 2 or more spp
#' 

# Specialists are C, D, and E
# Generalists are B and G

scen_df <- read_csv('species,type,s0,s1,s2,s3
A,other,0.1411,0.0,0.0,0.00
B,generalist,0.0418,0.0,0.3,0.00
C,specialist,0.0586,0.2,0.0,0.15
D,specialist,0.0902,0.5,0.0,0.40
E,specialist,0.143,0.3,0.0,0.25
F,other,0.0318,0.0,0.0,0.00
G,generalist,0.1415,0.0,0.7,0.20
H,other,0.1488,0.0,0.0,0.00
I,other,0.1041,0.0,0.0,0.00
J,other,0.0991,0.0,0.0,0.00', col_types = 'ccdddd') %>% 
    arrange(species)


scen_df <- scen_df %>% gather(scen, abund, s0:s3)

# # To replicate the random numbers under s0:
# set.seed(0)
# x <- runif(10) %>% round(., digits = 4)
# x <- round(x / sum(x), 4)
# i <- sample.int(length(x), 1)
# x[i] <- x[i] - (sum(x) - 1)
# sum(x); min(x)  # should be 1 and >0, respectively



blues <- c('#bdc9e1','#74a9cf','#0570b0')
oranges <- c('#feb24c','#f03b20')
greens <- c('#a1d99b','#74c476','#41ab5d','#238b45','#005a32')

scen_df %>% 
    filter(type != 'other', scen != 's0') %>% 
    ggplot(aes(scen, abund, fill = species)) +
    theme_classic() +
    geom_bar(stat = 'identity', position = 'stack') +
    scale_fill_manual(values = c(oranges[1], blues, oranges[2])) +
    theme(legend.position = 'none') +
    xlab('Scenario') +
    ylab('Abundance')


scen_df %>% filter(scen == 's2', abund > 0)


scen_df %>% filter(scen == 's3', abund > 0)




scen_fun <- function(scenario, sim_sd = 0.05, sim_N = 8) {

    abunds <- scen_df %>% filter(scen == scenario) %>% arrange(species)
    
    mat <- matrix(rep(abunds$abund, sim_N), nrow = sim_N, byrow = TRUE)
    
    ran_mat <- t(apply(mat, 1, function(x){
        y <- x[x > 0]
        y <- y + rnorm(length(y), 0, sim_sd)
        y <- ifelse(y > 0, y, 0)
        # Dealing with zeros separately
        z <- x[x == 0]
        z <- z + (rbinom(length(z), 1, prob = sim_sd) * 0.01)
        # Now make sure it adds to one
        out <- x
        out[x > 0] <- y
        out[x == 0] <- z
        out <- out / sum(out)
        return(out)
    }))
    
    ran_df <- as_data_frame(ran_mat)
    colnames(ran_df) <- unique(abunds$species)
    ran_df <- gather(ran_df, species, abundance)
    
    return(mutate(ran_df, scen = scenario, 
                  cage = rep(1:sim_N, length(unique(abunds$species)))))
}


set.seed(3)
iscen_df <- lapply(unique(scen_df$scen), scen_fun, sim_sd = 0.2) %>% 
    bind_rows %>% 
    arrange(species) %>% 
    mutate(species = factor(species), 
           scen = factor(scen, levels = paste0('s',0:3), 
                         labels = c('Random', 'Specialists', 'Generalists', 'Mix')))

colors <- character(length(unique(scen_df$species)))
colors[(scen_df$type == 'generalist')[1:10]] <- oranges
colors[(scen_df$type == 'specialist')[1:10]] <- blues
colors[(scen_df$type == 'other')[1:10]] <- greens

abund_p <- ggplot(iscen_df, aes(cage, abundance, fill = species)) +
    theme_bw() +
    theme(strip.background = element_rect(color = NA, fill = NA), 
          strip.text = element_text(face = 'bold'),
          legend.position = 'none', panel.grid = element_blank()) +
    geom_bar(stat = 'identity', position = 'stack') +
    facet_wrap(~ scen, nrow = 1) +
    scale_fill_manual(values = colors) +
    xlab('Cage number') +
    ylab('Relative abundance')
abund_p


# ggsave('~/Desktop/abund_p.pdf', abund_p, width = 6, height = 3, units = 'in')



