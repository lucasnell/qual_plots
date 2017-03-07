library(readr)
library(tidyr)
library(ggplot2)
library(taxize)
library(ape)
library(phytools)
library(ggtree)
library(ggplot2)

#' 
#' This contains functions to search WOS and loads the following additional packages:
#' 
#' * `dplyr`
#' * `purrr`
#' * `stringr`
#' * `stringi`
#' * `xml2`
#' 

# Made a symlink to access the folder this file is in
source('./wos/wos_search.R')


load('classifications.RData')


# Read CSV of data on sequenced genomes and standardize column names
all_genomes <- read_csv("~/Google Drive/Wisconsin/gwas/genomes_2017-02-12.csv.gz", 
                        col_types = str_c(rep('c', 19), collapse = '')) %>% 
    select(-Replicons) %>% 
    rename_(
        .dots = setNames(
            sprintf('`%s`', colnames(.)), 
            colnames(.) %>% str_replace_all(c('[#()%]' = '', '[/ ]' = '_')) %>% 
                str_to_lower))

# Change from organism names to species names and strain/subspecies names
sp_names <- all_genomes$organism_name %>% 
    str_split(' ') %>% 
    map_chr(~ str_c(.x[1:2], collapse = ' '))

sub_str_names <- all_genomes$organism_name %>% 
    str_split(' ') %>% 
    map_chr(~ ifelse(length(.x) > 2, str_c(.x[c(-1,-2)], collapse = ' '), NA))

all_genomes <- all_genomes %>% 
    select(-organism_name) %>% 
    mutate(sp_name = sp_names, sub_str_name = sub_str_names) %>% 
    select(sp_name, sub_str_name, everything()) %>% 
    filter(!grepl(pattern = 'sp.', x = sp_name), !grepl(pattern = ' x$', x = sp_name))


insect_names <- (all_genomes %>% filter(subgroup == 'Insects'))$sp_name %>% unique

plant_names <- (all_genomes %>% filter(subgroup == 'Land Plants'))$sp_name %>% unique


# lp_tr <- read.tree('land_plants.nwk')
# plot(lp_tr, cex = 0.5, no.margin = TRUE)  #, type = 'fan')




# ======================================================================================
# ======================================================================================

#       Insect families investigated for presence of Wolbachia

# ======================================================================================
# ======================================================================================


wolb_df <- read_csv('./wos/Duron.csv', col_types = paste(
    c(rep('c', 5), 'iii', rep('d', 21)), 
    collapse = '')) %>% 
    na.omit %>% 
    mutate(species = gsub('  ', ' ', species)) %>% 
    mutate(N = N_f + N_m + N_u,
           genus = str_split(species, ' ', simplify = TRUE)[,1],
           species = str_split(species, ' ', simplify = TRUE)[, -1] %>% 
               apply(., 1, function(x) {
                   z <- paste(x, collapse = ' ')
                   return(gsub(' $', '', z))
               }),
           wolb = wolb_o) %>% 
    select(class, order, family, genus, species, N, wolb)

# Summarize by species
wolb_df <- wolb_df %>% 
    group_by(class, order, family, genus, species) %>% 
    summarize(N = sum(N), wolb = mean(wolb)) %>% 
    ungroup

ins_wolb <- wolb_df %>% filter(class == 'Insecta')


# Manually adding two Drosophila species from the study by Mateos et al. (2006)
# [doi: 10.1534/genetics.106.058818]
ins_wolb <- add_row(ins_wolb, class = 'Insecta', order = 'Diptera', 
                    family = 'Drosophilidae', genus = 'Drosophila', 
                    species = c('melanogaster', 'macrospina'), N = c(45,7), 
                    wolb = c(0.75,0))

# ins_taxa <- classification(insect_names, db = 'ncbi')

ins_families <- lapply(ins_taxa, function(z) z$name[z$rank == 'family'][1]) %>% 
    unlist(., use.names = FALSE) %>% 
    ifelse(. == 'Termitoidae', 'Termopsidae', .)

ins_all_fams <- read.tree('insects_family.nwk')
ins_fam_tr <- drop.tip(
    ins_all_fams, 
    ins_all_fams$tip.label[!ins_all_fams$tip.label %in% ins_families]
    )



# Species names that have been sequenced whose families are represented by Wolbachia 
# study
surveyed <- insect_names[ins_families %in% ins_wolb$family]

# Check if a species name is part of a family surveyed for Wolbachia
# and if so, whether it was determined to have been infected
check_wolb <- function(sp_name){
    fam_name <- ins_families[insect_names == sp_name]
    if (!fam_name %in% ins_wolb$family){
        ret <- 1
    } else {
        wolb <- filter(ins_wolb, family == fam_name)$wolb
        if (all(wolb == 0)){
            ret <- 2
        } else {
            ret <- 3
        }
    }
    return(ret)
}

# Same, but by family
check_wolb_fam <- function(fam_name){
    if (!fam_name %in% ins_wolb$family){
        ret <- 1
    } else {
        wolb <- filter(ins_wolb, family == fam_name)$wolb
        if (all(wolb == 0)){
            ret <- 2
        } else {
            ret <- 3
        }
    }
    return(ret)
}


ins_tr <- read.tree('./wos/insects.nwk')
ins_tr$tip.label <- gsub('_', ' ', ins_tr$tip.label)



ins_df <- data_frame(sp = ins_tr$tip.label, 
                     wolb = factor(
                         sapply(ins_tr$tip.label, check_wolb, USE.NAMES = FALSE),
                         levels = 1:3,
                         labels = c('unknown', 'none', 'infected')
                         ))


ins_fam_df <- data_frame(fam = ins_fam_tr$tip.label, 
                     wolb = factor(
                         sapply(ins_fam_tr$tip.label, check_wolb_fam, USE.NAMES = FALSE),
                         levels = 1:3,
                         labels = c('unknown', 'none', 'infected')
                     ))



i <- ggtree(ins_fam_tr)
i$data$x <- i$data$x - max(i$data$x)
ip <- i %<+% ins_fam_df +
    theme_tree2(axis.title.x = element_text(size = 14)) +
    # geom_treescale(x = 20, y = 30, width = 50, linesize = 1, offset = -2) +
    geom_tiplab(aes(color = wolb), size = 2, fontface = 'bold.italic') +
    scale_x_continuous('Time (mya)', limits = c(-500, 50),
                       breaks = seq(-500, 0, 100), labels = seq(500, 0, -100)) +
    scale_color_manual(values = c('gray80','#1b9e77','#d95f02'))# +
ip

# ggsave(ip, filename = '~/Desktop/ip.pdf', width = 6, height = 8, units = 'in')





# i <- ggtree(ins_tr)
# i %<+% ins_df +
#     theme_tree2() +
#     # geom_treescale(x = 20, y = 30, width = 50, linesize = 1, offset = -2) +
#     geom_tiplab(aes(color = wolb), size = 1, fontface = 'bold.italic') +
#     scale_x_continuous(limits = c(0, 550)) +
#     scale_color_manual(values = c('gray80','#1b9e77','#d95f02','#7570b3'))# +
#     # theme(legend.position = c(0.25, 0.75))

# cat(paste(unique(ins_families), collapse = '\n'))






# ======================================================================================
# ======================================================================================

#       Plant families investigated for infection by tobacco mosaic virus (TMV)

# ======================================================================================
# ======================================================================================

plant_tr <- read.tree('./wos/land_plants.nwk')
plant_tr$tip.label <- gsub('_', ' ', plant_tr$tip.label)


tmv_df <- read_csv('./wos/Bald.csv', col_types = 'ccciiiiii') %>% 
    separate(orders_families, c('of1', 'of2'), sep = " ", fill = 'right') %>% 
    gather('w', 'order_family', of1:of2, na.rm = TRUE) %>% 
    select(group, subdivision, order_family, starts_with('tmv_'))

# plant_taxa <- classification(plant_names, db = 'ncbi')


# tmv_taxa <- classification(unique(tmv_df$order_family), db = 'ncbi')
tmv_fams <- tmv_taxa %>% 
    discard(~ is.null(nrow(.x))) %>% 
    names


plant_families <- lapply(plant_taxa, function(z) z$name[z$rank == 'family'][1]) %>% 
    unlist(., use.names = FALSE)
plant_orders <- lapply(plant_taxa, function(z) z$name[z$rank == 'order'][1]) %>% 
    unlist(., use.names = FALSE)

# tmv_fams[tmv_fams %in% c(plant_orders, plant_families)]

# plant_families <- unique(c(order_fams, plant_families))



# Families in the order(s) that've been examined for tmv
order_fams <- plant_taxa %>% 
    keep(~ any(.x$name %in% tmv_fams[tmv_fams %in% plant_orders])) %>% 
    map_chr(~ .x$name[.x$rank == 'family']) %>% 
    unique



plant_all_fams <- read.tree('plants_family.nwk')
plant_fam_tr <- drop.tip(
    plant_all_fams, 
    plant_all_fams$tip.label[!plant_all_fams$tip.label %in% unique(plant_families)]
)




# Check if a species name is part of a family surveyed for TMV
# and if so, whether it was determined to have been infected
check_tmv <- function(sp_name){
    
    fam_name <- plant_families[plant_names == sp_name]
    ord_name <- plant_orders[plant_names == sp_name]
    
    if (!fam_name %in% tmv_fams & !ord_name %in% tmv_fams){
        tmv <- as_data_frame(matrix(rep(NA,3), nrow = 1))
        colnames(tmv) <- c('tmv_systemic', 'tmv_leaf', 'tmv_none')
    } else {
        tmv <- filter(tmv_df, order_family %in% c(fam_name, ord_name)) %>% 
            select(starts_with('tmv_'))
    }
    return(tmv)
}




# Check if a family name has been surveyed for TMV
# and if so, whether it was determined to have been infected
check_tmv_fam <- function(fam_name){
    
    ord_name <- plant_orders[plant_families == fam_name][1]
    
    if (!fam_name %in% tmv_fams & !ord_name %in% tmv_fams){
        tmv <- as_data_frame(matrix(rep(NA,3), nrow = 1))
        colnames(tmv) <- c('tmv_systemic', 'tmv_leaf', 'tmv_none')
    } else {
        tmv <- filter(tmv_df, order_family %in% c(fam_name, ord_name)) %>% 
            select(starts_with('tmv_'))
    }
    return(tmv)
}




plant_df <- lapply(plant_tr$tip.label, check_tmv) %>% 
    bind_rows %>% 
    mutate(species = plant_tr$tip.label,
           inf = 1 - (tmv_none / (tmv_systemic + tmv_leaf + tmv_none)),
           status = factor(ifelse(is.na(inf), 0, ifelse(inf == 0, 1, 2)), 
                           levels = 0:2, 
                           labels = c('unknown', 'uninfected', 'infected'))) %>%
    select(species, everything())

plant_fam_df <- lapply(plant_fam_tr$tip.label, check_tmv_fam) %>% 
    bind_rows %>% 
    mutate(fam = plant_fam_tr$tip.label,
           inf = 1 - (tmv_none / (tmv_systemic + tmv_leaf + tmv_none)),
           status = factor(ifelse(is.na(inf), 0, ifelse(inf == 0, 1, 2)), 
                           levels = 0:2, 
                           labels = c('unknown', 'uninfected', 'infected'))) %>%
    select(fam, everything())







p <- ggtree(plant_fam_tr)
p$data$x <- p$data$x - max(p$data$x)
pp <- p %<+% plant_fam_df +
    theme_tree2(axis.title.x = element_text(size = 14)) +
    geom_tiplab(aes(color = status), size = 2, fontface = 'bold.italic') +
    scale_x_continuous('Time (mya)', limits = c(-550, 50),
                       breaks = seq(-500, 0, 100), labels = seq(500, 0, -100)) +
    scale_color_manual(values = c('gray80','#1b9e77','#d95f02'))
pp

# ggsave(pp, filename = '~/Desktop/pp.pdf', width = 6, height = 8, units = 'in')





# p <- ggtree(plant_tr)
# p %<+% plant_df +
#     theme_tree2() +
#     geom_tiplab(aes(color = status), size = 1, fontface = 'bold.italic') +
#     # geom_tippoint(aes(color = status), size = 3) +
#     scale_x_continuous(limits = c(0, 550)) +
#     scale_color_manual(values = c('gray80','#1b9e77','#7570b3')) +
#     theme(legend.position = c(0.25, 0.75))



# save(ins_taxa, plant_taxa, tmv_taxa, file = 'classifications.RData')




# # all_genomes %>% 
# #     filter(subgroup == 'Insects', is.na(sub_str_name)) %>% 
# #     group_by(sp_name) %>% 
# #     summarize(N = n()) %>% 
# #     filter(N > 1) %>% as.data.frame
# # 
# # 
# # all_genomes %>% 
# #     filter(subgroup == 'Insects', !is.na(sub_str_name)) %>% 
# #     group_by(sp_name, sub_str_name) %>% 
# #     summarize(N = n()) %>% 
# #     as.data.frame
# 
# 
# 
# pub_test <- wos_search(
#     list('topic' = sprintf('"Wolbachia" AND "%s"', insect_names[1])))# %>% 
#     # lapply(., clean_entry)# %>% 
#     # bind_rows %>% 
#     # mutate(url = get_url(doi))
# 
# bind_rows(pub_test)
# 
# 
# 
# 
# 
# ------------------------
# This was run to see how many species came up with the above search. Answer = very few.
# ------------------------

# check_pres <- function(sp_name) {
#     wos_list <- wos_search(list('topic' = sprintf('"Wolbachia" AND "%s"', sp_name)))
#     return(ifelse(length(wos_list) > 0, 1, 0))
# }
# 
# 
# # Checks for whether each insect species has >1 Wolbachia paper
# insect_wols <- sapply(insect_names, check_pres, USE.NAMES = FALSE)
# plant_tmv <- sapply(plant_names, check_pres, USE.NAMES = FALSE)
# 
# insect_names[insect_wols == 1]
# plant_names[plant_tmv == 1]









# ====================================================================================
# ====================================================================================
# ====================================================================================
# ====================================================================================
# ====================================================================================
# ====================================================================================
# ====================================================================================
# ====================================================================================
# ====================================================================================





# paste0('Wolbachia AND (', 
#        paste0('"', insect_names, '"', collapse = ' OR '), ')') %>% cat
# 
# paste0('"tobacco mosaic virus" AND (', 
#        paste0('"', plant_names, '"', collapse = ' OR '), ')') %>% cat
# 
# 
# wol_refs <- read_tsv('year	records	percent
# 1989	2	0.123
# 1991	2	0.123
# 1992	4	0.247
# 1993	9	0.555
# 1994	9	0.555
# 1995	16	0.986
# 1996	12	0.740
# 1997	18	1.110
# 1998	26	1.603
# 1999	25	1.541
# 2000	30	1.850
# 2001	31	1.911
# 2002	47	2.898
# 2003	45	2.774
# 2004	60	3.699
# 2005	56	3.453
# 2006	71	4.377
# 2007	55	3.391
# 2008	60	3.699
# 2009	92	5.672
# 2010	77	4.747
# 2011	116	7.152
# 2012	123	7.583
# 2013	147	9.063
# 2014	157	9.679
# 2015	152	9.371
# 2016	175	10.789
# 2017	5	0.308')
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# tmv_refs <- read_tsv('year	records	percent
# 1966	3	0.222
# 1967	6	0.443
# 1969	2	0.148
# 1970	2	0.148
# 1971	4	0.295
# 1972	5	0.369
# 1973	1	0.074
# 1974	2	0.148
# 1975	4	0.295
# 1976	2	0.148
# 1977	2	0.148
# 1978	3	0.222
# 1979	3	0.222
# 1980	2	0.148
# 1981	8	0.591
# 1982	1	0.074
# 1983	3	0.222
# 1984	4	0.295
# 1985	3	0.222
# 1986	5	0.369
# 1988	2	0.148
# 1989	4	0.295
# 1990	3	0.222
# 1991	32	2.363
# 1992	38	2.806
# 1993	36	2.659
# 1994	34	2.511
# 1995	38	2.806
# 1996	30	2.216
# 1997	37	2.733
# 1998	36	2.659
# 1999	51	3.767
# 2000	42	3.102
# 2001	45	3.323
# 2002	42	3.102
# 2003	54	3.988
# 2004	70	5.170
# 2005	49	3.619
# 2006	49	3.619
# 2007	60	4.431
# 2008	52	3.840
# 2009	56	4.136
# 2010	52	3.840
# 2011	59	4.357
# 2012	63	4.653
# 2013	57	4.210
# 2014	70	5.170
# 2015	60	4.431
# 2016	67	4.948
# 2017	1	0.074')
# 
# 
# wol_N <- (wol_refs %>% filter(year != 2017))$records %>% sum
# 
# wol_refs %>% 
#     filter(year != 2017) %>% 
#     ggplot(aes(year, records)) +
#     labs(y = 'Reference counts', x = NULL, 
#          # title = 'Papers through time',
#          # subtitle = paste('Search terms: "Wolbachia" and any insect species with',
#          #                  'a sequenced genome'),
#          caption = 'Source: Web of Science') +
#     theme_classic() +
#     geom_line(color = 'dodgerblue') +
#     annotate(geom = 'text', label = paste('Total =', 
#                                           format(wol_N, big.mark = ',')),
#              x = 2000, y = 75, vjust = 0)
# 
# tmv_N <- (tmv_refs %>% filter(year != 2017))$records %>% sum
# tmv_refs %>% 
#     filter(year != 2017) %>% 
#     ggplot(aes(year, records)) +
#     labs(y = 'Reference counts', x = 'Year', 
#          # title = 'Papers through time',
#          # subtitle = paste('Search terms: "Tobacco mosaic virus" and any plant species',
#          #                  'with a sequenced genome'),
#          caption = 'Source: Web of Science') +
#     theme_classic() +
#     geom_line(color = 'dodgerblue') +
#     annotate(geom = 'text', label = paste('Total =', 
#                                           format(tmv_N, big.mark = ',')),
#              x = 1989, y = 40, vjust = 0, hjust = 1)








