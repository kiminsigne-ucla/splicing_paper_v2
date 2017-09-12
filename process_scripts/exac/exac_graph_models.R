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
    geom_line(aes(color = method), size = 1.25) + 
    # scale_y_log10() +
    scale_y_log10(breaks = c(3.6, 10, 100), limits = c(2,100)) +
    annotation_logticks(sides = 'l') +
    labs(x = 'Recall (%)', y = 'Precision (%)', color = '') +
    theme(legend.position = c(0.80, 0.80),
          legend.key = element_rect(size = 5),
          legend.key.size = unit(1.25, 'lines'),
          legend.text = element_text(size = 16), 
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 14, color = 'grey20')) +
    geom_hline(yintercept = 3.6, linetype = 'dashed', color = 'grey40') 

ggsave(paste0('../../figs/exac/exac_fig4E_exac_pr_curves', '.tiff'), 
       height = 4, width = 5, units = 'in', dpi = 600)

# split by intron/exon
ggplot(pr_curve_info, aes(recall, precision)) + geom_line(aes(color = method)) +
    scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100), limits = c(0.001, 10)) + 
    annotation_logticks(sides = 'l') +
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
