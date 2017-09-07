load_pkgs <- function(pkgs){
    new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, 'Package'])]
    if(length(new_pkgs)) install.packages(new_pkgs)
    for(pkg in pkgs){
        suppressWarnings(suppressMessages(library(pkg, character.only = T)))
    }
}

pkgs <- c('dplyr', 'tidyr', 'ggplot2', 'cowplot', 'grid')
load_pkgs(pkgs)

options(stringsAsFactors = F, warn = -1, warnings = -1)

plot_format <- '.png'

###############################################################################
# precision-recall curves
###############################################################################
pr_curve_all <- read.table('../../processed_data/exac/exac_models_pr_curves_all.txt', 
                           sep = '\t', header = T) %>% 
    filter(method != 'hal') %>% # only scores exonic variants
    mutate(type = 'all')

pr_curve_exon <- read.table('../../processed_data/exac/exac_models_pr_curves_exon.txt', 
                            sep = '\t', header = T) %>% 
    mutate(type = 'exon')

pr_curve_intron <- read.table('../../processed_data/exac/exac_models_pr_curves_intron.txt', 
                              sep = '\t', header = T) %>% 
    mutate(type = 'intron')

pr_curve_info <- bind_rows(pr_curve_all, pr_curve_exon) %>% 
    bind_rows(pr_curve_intron)

# all variants
pr_curve_all %>% 
    filter(method != 'fathmm_noncoding') %>% # pick one FATHMM score method
    mutate(method = factor(method, labels = c('CADD', 'DANN', 'FATHMM-MKL', 'fitCons', 
                                              'LINSIGHT', 'SPANR'))) %>% 
    ggplot(aes(recall, precision)) + 
    geom_line(aes(color = method), size = 1) + 
    scale_y_log10() + annotation_logticks(sides = 'l') +
    labs(color = '') +
    theme(legend.position = c(0.80, 0.85)) +
    geom_hline(yintercept = 3.7, linetype = 'dashed', color = 'grey')

ggsave(paste0('../../figs/exac/exac_pr_curves', plot_format), 
       height = 4, width = 4, units = 'in')

# split by intron/exon
ggplot(pr_curve_info, aes(recall, precision)) + geom_line(aes(color = method)) +
    scale_y_log10() + annotation_logticks(sides = 'l') +
    facet_grid(~ type)

###############################################################################
# ROC curves
###############################################################################
roc_curve_all <- read.table('../../processed_data/exac/exac_models_roc_curves_all.txt', 
                            sep = '\t', header = T) %>% 
    filter(method != 'hal') %>% 
    mutate(type = 'all')

roc_curve_exon <- read.table('../../processed_data/exac/exac_models_roc_curves_exon.txt', 
                             sep = '\t', header = T) %>% 
    mutate(type = 'exon')

roc_curve_intron <- read.table('../../processed_data/exac/exac_models_roc_curves_intron.txt', 
                               sep = '\t', header = T) %>% 
    mutate(type = 'intron')

roc_curve_info <- bind_rows(roc_curve_all, roc_curve_exon) %>% 
    bind_rows(roc_curve_intron)


roc_curve_all %>% 
    filter(method != 'fathmm_noncoding') %>% # pick one FATHMM score method
    mutate(method = factor(method, labels = c('CADD', 'DANN', 'FATHMM', 'fitCons', 
                                              'LINSIGHT', 'SPANR'))) %>% 
    ggplot(aes(false_positive_rate, true_positive_rate)) + 
    geom_line(aes(color = method)) + 
    labs(color = '', x = 'false positive rate', y = 'true positive rate')

# split by intron/exon
ggplot(roc_curve_info, aes(false_positive_rate, true_positive_rate)) + 
    geom_line(aes(color = method)) +
    geom_abline(intercept = 0, linetype = 'dashed') +
    facet_grid(~type) +
    labs(x = 'false positive rate', y = 'true positive rate (recall)')

# ###############################################################################
# scores by annotation category
# ###############################################################################
data_annot <- read.table('../../processed_data/exac/exac_func_annot.txt',
                        sep = '\t', header = T)
# 
# hal <- read.table('../../processed_data/exac/exac_HAL_scores.txt',
#                    sep = '\t', header = T)
# 
# spanr <- read.table('../../processed_data/exac/exac_SPANR_scores_capped.txt',
#                     sep = '\t', header = T)
# 
# data_annot <- data_annot %>% 
#     left_join(select(hal, id, hal_dpsi = DPSI_pred), by = 'id') %>% 
#     left_join(select(spanr, id, spanr_dpsi = dpsi_spanr_capped), by = 'id')
# 
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

# gg1 <- data_annot %>% 
#     filter(!is.na(consequence)) %>% 
#     ggplot(aes(v2_dpsi, noncoding_score)) + geom_point(alpha = 0.50) +
#     facet_wrap(~ consequence, labeller = as_labeller(annot_names)) +
#     labs(x = expression(paste(Delta, ' inclusion index')),
#          y = 'FATHMM non-coding score')
# 
# gg2 <- data_annot %>% 
#     filter(!is.na(consequence)) %>% 
#     ggplot(aes(v2_dpsi, coding_score)) + geom_point(alpha = 0.50) +
#     facet_wrap(~ consequence, labeller = as_labeller(annot_names)) +
#     labs(x = expression(paste(Delta, ' inclusion index')),
#          y = 'FATHMM coding score')
# 
# plot_grid(gg1, gg2, nrow = 2)
# 
# 
# gg1 <- data_annot %>% 
#     filter(!is.na(consequence)) %>% 
#     ggplot(aes(v2_dpsi, dann_score)) + geom_point(alpha = 0.50) +
#     facet_wrap(~ consequence, labeller = as_labeller(annot_names)) +
#     labs(x = expression(paste(Delta, ' inclusion index')),
#          y = 'DANN score')
# 
# gg2 <- data_annot %>% 
#     filter(!is.na(consequence)) %>% 
#     ggplot(aes(v2_dpsi, cadd_score)) + geom_point(alpha = 0.50) +
#     facet_wrap(~ consequence, labeller = as_labeller(annot_names)) +
#     labs(x = expression(paste(Delta, ' inclusion index')),
#          y = 'CADD score')
# 
# plot_grid(gg1, gg2, nrow = 2)
# 
# 
# gg1 <- data_annot %>%
#     filter(!is.na(consequence)) %>%
#     ggplot(aes(v2_dpsi, fitCons_score)) + geom_point(alpha = 0.50) +
#     facet_wrap(~ consequence, labeller = as_labeller(annot_names)) +
#     labs(x = expression(paste(Delta, ' inclusion index')),
#          y = 'fitCons score')
# 
# 
# gg2 <- data_annot %>%
#     filter(!is.na(consequence)) %>%
#     ggplot(aes(v2_dpsi, linsight_score)) + geom_point(alpha = 0.50) +
#     facet_wrap(~ consequence, labeller = as_labeller(annot_names)) +
#     labs(x = expression(paste(Delta, ' inclusion index')),
#          y = 'LINSIGHT score')
#  
# plot_grid(gg1, gg2, nrow = 2)

