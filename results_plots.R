#' 
#' # Plots of potential results
#' 



suppressPackageStartupMessages({
    library(tidyverse)
})


#' 
#' <!-- ======================================================================== -->
#' 
#' <!-- Alternative states plot -->
#' 
#' <!-- ======================================================================== -->
#' 
#' 


k1 = k2 = 100

a12 = a21 = 2
m1 = -1 / a12
m2 = -a21

# Non-stable equilibrium coordinates
ns_equil <- c(x = {k1 - a12*k2} / {a12*(m2-m1)}, 
              y = m2 * ({k1 - a12*k2} / {a12*(m2-m1)}) + k2)
# Stable equilibrium coordinates
equil <- list(sp1 = c(x = k1, y = 0), sp2 = c(x = 0, y = k2))

# curve(m1 * x + (k1 / a12), 0, 100, xlim = c(0,k1), ylim = c(0,k1), 
#       xlab = expression(N[1]), ylab = expression(N[2]))
# curve(m2 * x + k2, 0, 100, add = TRUE, col = 'red')
# points(x = ns_equil['x'], y = ns_equil['y'])

# ----------
# Creating data frame for red and blue lines (the "zero isoclines")
# ----------
line_df <- data_frame(x = c(seq(0, 100, length.out = 11))) %>% 
    mutate(sp1 = m1 * x + (k1 / a12), sp2 = m2 * x + k2) %>% 
    gather(spp, y, sp1:sp2) %>% 
    mutate(spp = factor(spp)) %>% 
    filter(y >= 0)




# ----------
# Creating ggplot object with no arrows
# ----------
no_arrow_p <- ggplot(line_df, aes(x, y)) +
    theme_classic() +
    # Red and blue lines
    geom_line(aes(color = spp), size = 1) +
    # Hacky axes
    geom_path(data = data_frame(x = c(0, k1 + 10), y = c(0,0))) +
    geom_path(data = data_frame(x = c(0,0), y = c(0, k2 + 10))) +
    # Equilibrium points, first unstable, then stable
    geom_point(data = data_frame(x = ns_equil['x'], y = ns_equil['y']),
               shape = 21, color= 'black', fill = 'white', size = 3) +
    geom_point(data = data_frame(x = c(k1, 0), y = c(0, k2), 
                                 species = factor(c('sp1', 'sp2'))),
               aes(color = species), size = 3, shape = 16) +
    # Aesthetics
    scale_x_continuous(expression(bold(N[1])), limits = c(0, k1 + 10),
                       breaks = c(k2/a21, k1), 
                       labels = c(expression(K[2] / alpha[21]), expression(K[1]))) +
    scale_y_continuous(expression(bold(N[2])), limits = c(0, k2 + 10),
                       breaks = c(k1 / a12, k2), 
                       labels = c(expression(frac(K[1], alpha[12])), expression(K[2]))) +
    scale_color_manual(values = rgb(c(86, 230), c(180,159), c(233, 0), 
                                    maxColorValue = 255)) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), 
          axis.title = element_text(size = 14),
          axis.text.x = element_text(margin = margin(-8, 0, 0, 0), size = 10, 
                                     color = 'black'), 
          axis.text.y = element_text(margin = margin(0, -10, 0, 0), size = 10,
                                     color = 'black'), 
          legend.position = 'none')


no_arrow_p <- no_arrow_p + 
    geom_text(data = line_df %>% 
                  filter((spp == 'sp1' & x == 70) | (spp == 'sp2' & x == 20)), 
              aes(label = paste('Zero isocline for species', gsub('sp', '', spp)),
                  angle = ifelse(spp == 'sp1', {atan(m1) / (2 * pi)} * 360, 
                                 {atan(m2) / (2 * pi)} * 360), 
                  color = spp, x = x + 2.5, y = y + 2.5), 
              size = 3, fontface = 'bold')
no_arrow_p
#




# ggsave('~/Desktop/alt_states.pdf', no_arrow_p, width = 4, height = 4, units = 'in')



# Not doing below bc it's very time-consuming
# ----------
# Creating data frame for arrows
# ----------
m12 <- {(k2 + k1/a12)/2} / {(k1 + k2/a21)/2}  # mean(c(m1, m2))
b12 <- mean(c(k2, k1/a12))

arrow_fun <- function(x, direction) {
    if (direction == 'down') {
        -m12 * (x - ns_equil['x']) + ns_equil['y']
    }
}

# Arrows going downward from non-stable equilibrium
arrow_df <- data_frame(x = round(ns_equil['x'] + 5), x2 = round(ns_equil['x'] + 11)) %>%
    mutate(y = arrow_fun(x, 'down'), y2 = arrow_fun(x2, 'down'))

# # Arrows going upward from non-stable equilibrium
# arrow_df <- data_frame(x = round(ns_equil['x'] + 5), x2 = round(ns_equil['x'] + 11)) %>%
#     mutate(y = m12*x + b12, y2 = m12*x2 + b12)


p + geom_segment(data = arrow_df, aes(xend = x2, yend = y2),
                 arrow = arrow(length = unit(0.02, "npc"), type = 'closed'))










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
#' 0. Random - No apparent pattern
#' 1. Specialists - 3 spp
#' 2. Generalist - 1 or more spp
#' 3. Alternative states - Switch btwn individual generalists
#' 

# Specialists are C, D, and E
# Generalists are B and G

scen_df <- rbind(c('A','other',0.1411,0.0,0.0),
    c('B','generalist',0.0418,0.0,0.3),
    c('C','specialist',0.0586,0.2,0.0),
    c('D','specialist',0.0902,0.5,0.0),
    c('E','specialist',0.143,0.3,0.0),
    c('F','other',0.0318,0.0,0.0),
    c('G','generalist',0.1415,0.0,0.7),
    c('H','other',0.1488,0.0,0.0),
    c('I','other',0.1041,0.0,0.0),
    c('J','other',0.0991,0.0,0.0)) %>% 
    as_data_frame %>% 
    rename_(.dots = setNames(paste0('V', 1:5), 
                             c('species','type',paste0('s',0:2)))) %>% 
    mutate_at(3:5, as.numeric) %>% 
    arrange(species) %>%
    gather(scen, abund, s0:s2)

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


# Function to add noise to one row (i.e., one simulation/cage) of an abundance matrix
randomize_matrix <- function(x, sim_sd){
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
}


scen_fun <- function(scenario, sim_sd = 0.2, sim_N = 8) {

    abunds <- scen_df %>% filter(scen == scenario) %>% arrange(species)
    
    mat <- matrix(rep(abunds$abund, sim_N), nrow = sim_N, byrow = TRUE)
    
    ran_mat <- t(apply(mat, 1, randomize_matrix, sim_sd = sim_sd))
    
    ran_df <- as_data_frame(ran_mat)
    colnames(ran_df) <- unique(abunds$species)
    ran_df <- gather(ran_df, species, abundance)
    
    return(mutate(ran_df, scen = scenario, 
                  cage = rep(1:sim_N, length(unique(abunds$species)))))
}

sim_alt_states <- function(states, sim_sd = 0.2, sim_N = 8){
    abunds <- scen_df %>% filter(scen %in% states) %>% arrange(scen, species)
    
    ran_states <- sample(states, sim_N, replace = TRUE)
    
    ran_mat <- sapply(
        ran_states,
        function(s) {
            mat <- matrix(scen_df$abund[scen_df$scen == s], nrow = 1)
            randomize_matrix(mat, sim_sd)
        }, USE.NAMES = FALSE)
    ran_mat <- t(ran_mat)
    
    ran_df <- mutate(as_data_frame(ran_mat), cage = 1:sim_N, scen = 's3')
    colnames(ran_df)[1:ncol(ran_mat)] <- unique(abunds$species)
    ran_df <- gather(ran_df, species, abundance, -cage, -scen)
    
    return(ran_df)
}









set.seed(3)
iscen_df <- lapply(unique(scen_df$scen), scen_fun) %>% 
    append(., list(sim_alt_states(c('s1', 's2')))) %>% 
    bind_rows %>% 
    arrange(species) %>% 
    mutate(species = factor(species), 
           scen = factor(scen, levels = paste0('s',0:3), 
                         labels = c('Random','Specialists','Generalists','Alt. States')))



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



