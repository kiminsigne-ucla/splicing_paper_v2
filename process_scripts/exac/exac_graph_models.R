load_pkgs <- function(pkgs){
    new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, 'Package'])]
    if(length(new_pkgs)) install.packages(new_pkgs)
    for(pkg in pkgs){
        suppressWarnings(suppressMessages(library(pkg, character.only = T)))
    }
}

pkgs <- c('dplyr', 'tidyr', 'ggplot2', 'cowplot')
load_pkgs(pkgs)

options(stringsAsFactors = F, warn = -1, warnings = -1)

plot_format <- '.png'

pr_curve_info <- read.table('../../processed_data/exac/exac_models_pr_curves.txt', 
                           sep = '\t', header = T)

ggplot(pr_curve_info, aes(recall, precision)) + geom_line(aes(color = method)) +
    scale_y_log10() + annotation_logticks(sides = 'l')

ggplot(pr_curve_info, aes(recall, precision)) + geom_line() +
    facet_wrap(~ method) +
    scale_y_log10() + annotation_logticks(sides = 'l')

data_annot <- read.table('../../processed_data/exac/exac_func_annot.txt',
                        sep = '\t', header = T)

hal <- read.table('../../processed_data/exac/exac_HAL_scores.txt',
                   sep = '\t', header = T)

spanr <- read.table('../../processed_data/exac/exac_SPANR_scores_capped.txt',
                    sep = '\t', header = T)

data_annot <- data_annot %>% 
    left_join(select(hal, id, hal_dpsi = DPSI_pred), by = 'id') %>% 
    left_join(select(spanr, id, spanr_dpsi = dpsi_spanr_capped), by = 'id')

data_annot$consequence <- factor(data_annot$consequence, 
                                 levels = c('splice_acceptor_variant', 
                                            'splice_donor_variant',
                                            'stop_gained', 'missense_variant', 
                                            'synonymous_variant', 'intron_variant'))

annot_names <- c(
    'splice_acceptor_variant' = 'splice acceptor variant',
    'splice_donor_variant' = 'splice donor variant',
    'stop_gained' = 'stop gained variant',
    'missense_variant' = 'missense variant',
    'synonymous_variant' = 'synonymous variant',
    'intron_variant' = 'intron variant'
)

gg1 <- data_annot %>% 
    filter(!is.na(consequence)) %>% 
    ggplot(aes(v2_dpsi, noncoding_score)) + geom_point(alpha = 0.50) +
    facet_wrap(~ consequence, labeller = as_labeller(annot_names)) +
    labs(x = expression(paste(Delta, ' inclusion index')),
         y = 'FATHMM non-coding score')

gg2 <- data_annot %>% 
    filter(!is.na(consequence)) %>% 
    ggplot(aes(v2_dpsi, coding_score)) + geom_point(alpha = 0.50) +
    facet_wrap(~ consequence, labeller = as_labeller(annot_names)) +
    labs(x = expression(paste(Delta, ' inclusion index')),
         y = 'FATHMM coding score')

plot_grid(gg1, gg2, nrow = 2)


gg1 <- data_annot %>% 
    filter(!is.na(consequence)) %>% 
    ggplot(aes(v2_dpsi, dann_score)) + geom_point(alpha = 0.50) +
    facet_wrap(~ consequence, labeller = as_labeller(annot_names)) +
    labs(x = expression(paste(Delta, ' inclusion index')),
         y = 'DANN score')

gg2 <- data_annot %>% 
    filter(!is.na(consequence)) %>% 
    ggplot(aes(v2_dpsi, cadd_score)) + geom_point(alpha = 0.50) +
    facet_wrap(~ consequence, labeller = as_labeller(annot_names)) +
    labs(x = expression(paste(Delta, ' inclusion index')),
         y = 'CADD score')

plot_grid(gg1, gg2, nrow = 2)


gg1 <- data_annot %>% 
    filter(!is.na(consequence)) %>% 
    ggplot(aes(v2_dpsi, fitCons_score)) + geom_point(alpha = 0.50) +
    facet_wrap(~ consequence, labeller = as_labeller(annot_names)) +
    labs(x = expression(paste(Delta, ' inclusion index')),
         y = 'fitCons score')


gg2 <- data_annot %>% 
    filter(!is.na(consequence)) %>% 
    ggplot(aes(v2_dpsi, linsight_score)) + geom_point(alpha = 0.50) +
    facet_wrap(~ consequence, labeller = as_labeller(annot_names)) +
    labs(x = expression(paste(Delta, ' inclusion index')),
         y = 'LINSIGHT score')

plot_grid(gg1, gg2, nrow = 2)

