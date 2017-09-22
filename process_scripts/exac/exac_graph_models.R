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

dir <- '../../figs/exac/'
plot_format_main <- '.png'
plot_format <- '.png'
hi_res <- 600
lo_res <- 300

color_cadd = '#FF9912'
color_dann = '#941494'
color_fathmm = '#108C44'
color_fitcons = '#6AA5CD'
color_linsight = '#ABB9B9'
color_spanr =  '#ED1E24'
color_hal = '#000080'

###############################################################################
# Precision-recall curves
###############################################################################
pr_curve_all <- 
  read.table('../../processed_data/exac/exac_models_pr_curves_all.txt', 
                           sep = '\t', header = T) %>% 
    filter(method != 'hal') %>% # only scores exonic variants
    mutate(type = 'all')

pr_curve_exon <- 
  read.table('../../processed_data/exac/exac_models_pr_curves_exon.txt', 
                            sep = '\t', header = T) %>% 
    mutate(type = 'exon')

pr_curve_intron <- 
  read.table('../../processed_data/exac/exac_models_pr_curves_intron.txt', 
                              sep = '\t', header = T) %>% 
    mutate(type = 'intron')

pr_curve_info <- bind_rows(pr_curve_all, pr_curve_exon) %>% 
    bind_rows(pr_curve_intron)

# all variants
pr_curve_all %>% 
    filter(method != 'fathmm_noncoding') %>% # pick one FATHMM score method
    mutate(method = factor(method, 
                           labels = c('CADD', 'DANN', 'FATHMM-MKL', 'fitCons', 
                                              'LINSIGHT', 'SPANR'))) %>% 
    ggplot(aes(recall, precision)) + 
    geom_line(aes(color = method), size = 1.25) + 
    # scale_y_log10() +
    scale_y_log10(breaks = c(3.6, 10, 100), limits = c(2,100)) +
    annotation_logticks(sides = 'l') +
    labs(x = 'Recall (%)', y = 'Precision (%)', color = '') +
    theme(
          axis.title.x = element_text(size = 19, vjust = -1.5),
          axis.title.y = element_text(size = 19, vjust = -1.5),
          axis.line.x = element_line(color = 'grey30'),
          axis.line.y = element_line(color = 'grey30'),
          axis.ticks = element_line(color = 'grey30'),
          axis.text = element_text(size = 14, color = 'grey20'),
          legend.justification = 'center',
          legend.direction = 'vertical',
          legend.title = element_blank(),
          legend.text = element_text(size = 16),
          legend.key.size = unit(0.5, 'lines'),
          legend.key.height = unit(0.25, "inch"),
          legend.key.width = unit(0.6, "inch")) +
    geom_hline(yintercept = 3.6, linetype = 'dashed', color = 'grey40') +
    scale_color_manual(values = c(color_cadd, color_dann, color_fathmm, 
                                color_fitcons, color_linsight, color_spanr))

ggsave(paste0(dir, 'exac_fig4E_exac_pr_curves_with_legend', 
              plot_format_main), 
       height = 4, width = 7, units = 'in', dpi = hi_res)
       

pr_curve_all %>% 
  filter(method != 'fathmm_noncoding') %>% # pick one FATHMM score method
  mutate(method = factor(method, 
                         labels = c('CADD', 'DANN', 'FATHMM-MKL', 'fitCons', 
                                            'LINSIGHT', 'SPANR'))) %>% 
  ggplot(aes(recall, precision)) + 
  geom_line(aes(color = method), size = 1.25) + 
  # scale_y_log10() +
  scale_y_log10(breaks = c(3.6, 10, 100), limits = c(2,100)) +
  annotation_logticks(sides = 'l') +
  labs(x = 'Recall (%)', y = 'Precision (%)', color = '') +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14, color = 'grey20'),
        axis.line.x = element_line(color = 'grey30'),
        axis.line.y = element_line(color = 'grey30'),
        axis.ticks = element_line(color = 'grey30'),
        legend.position = 'none') +
  geom_hline(yintercept = 3.6, linetype = 'dashed', color = 'grey40') +
  scale_color_manual(values = c(color_cadd, color_dann, color_fathmm, 
                                color_fitcons, color_linsight, color_spanr))

ggsave(paste0(dir, 'exac_fig4E_exac_pr_curves_no_legend', 
       plot_format_main), height = 4, width = 4, units = 'in', dpi = hi_res)

# split by intron/exon
pr_curve_info %>%
  filter(method != 'fathmm_noncoding') %>%# pick one FATHMM score method
  mutate(method = factor(method, labels = c('CADD', 'DANN', 'FATHMM-MKL', 
                                            'fitCons', 'HAL', 'LINSIGHT', 
                                            'SPANR')),
         type   = factor(type, levels = c('exon', 'intron', 'all'), 
                         labels = c('Exonic SNV', 'Intronic SNV', 'All SNV'))) %>%
    ggplot(aes(recall, precision)) + geom_line(aes(color = method)) +
    # scale_y_log10() +
    scale_y_log10(breaks = c(0.01, 0.1, 1, 3.6, 10, 100)) +
    annotation_logticks(sides = 'l') +
    facet_grid(~ type) +
    scale_color_manual(values = c(color_cadd, color_dann, color_fathmm, 
                                color_fitcons, color_hal, color_linsight, color_spanr)) +
  labs(x = 'Recall (%)', y = 'Precision (%)', color = '') +
  geom_hline(yintercept = 3.6, linetype = 'dashed', color = 'grey40') +
  theme(strip.text = element_text(size = 18.5),
        strip.background = element_rect(fill = "#E8E8E8", color = "white"),
        panel.grid = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.border = element_rect(fill = NA, color = 'grey50'),
        axis.title.y = element_text(size = 19, vjust = 2.75),
        axis.title.x = element_text(size = 19, vjust = -2.25),
        axis.text.y = element_text(size = 14, color = 'grey30'),
        axis.text.x = element_text(size = 14, color = 'grey30'),
        axis.line.y = element_line(color = 'grey30'),
        axis.line.x = element_line(color = 'grey30'),
        axis.ticks = element_line(color = 'grey30'),
        legend.justification = 'center',
        legend.direction = 'vertical',
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, 'lines'),
        legend.key.height = unit(0.25, "inch"),
        legend.key.width = unit(0.6, "inch")) +
  guides(colour = guide_legend(override.aes = list(size = 3)))

ggsave(paste0(dir, 'exac_fig4E_exac_pr_curves_type', plot_format), 
       height = 4.5, width = 12, units = 'in', dpi = lo_res)

###############################################################################
# ROC curves
###############################################################################
roc_curve_all <- 
  read.table('../../processed_data/exac/exac_models_roc_curves_all.txt', 
                            sep = '\t', header = T) %>% 
    filter(method != 'hal') %>% 
    mutate(type = 'all')

roc_curve_exon <- 
  read.table('../../processed_data/exac/exac_models_roc_curves_exon.txt', 
                             sep = '\t', header = T) %>% 
    mutate(type = 'exon')

roc_curve_intron <- 
  read.table('../../processed_data/exac/exac_models_roc_curves_intron.txt', 
                               sep = '\t', header = T) %>% 
    mutate(type = 'intron')

roc_curve_info_main <- bind_rows(roc_curve_all, roc_curve_exon) %>% 
    bind_rows(roc_curve_intron) %>%
    filter(method != 'fathmm_noncoding') %>%
    mutate(method = factor(method, labels = c('CADD', 'DANN', 'FATHMM-MKL', 
                                            'fitCons', 'HAL', 'LINSIGHT', 
                                            'SPANR')),
           type   = factor(type, levels = c('exon', 'intron', 'all'), 
                       labels = c('Exonic SNV', 'Intronic SNV', 'All SNV'))) 
                                                                                     
roc_curve_all %>% 
    filter(method != 'fathmm_noncoding') %>% # pick one FATHMM score method
    mutate(method = factor(method, labels = c('CADD', 'DANN', 'FATHMM-MKL', 
                                              'fitCons', 'LINSIGHT', 'SPANR'))) %>% 
    ggplot(aes(false_positive_rate, true_positive_rate)) + 
    geom_line(aes(color = method)) +
    labs(color = '', x = 'False positive rate', y = 'True positive rate')

# split by intron/exon
ggplot(roc_curve_info_main, aes(false_positive_rate, true_positive_rate)) + 
    geom_line(aes(color = method), size = 1.25) +
    geom_abline(intercept = 0, linetype = 'dashed', color = 'grey50') +
    facet_grid( ~ type) +
    theme(strip.text = element_text(size = 18.5),
          strip.background = element_rect(fill = "#E8E8E8", color = "white"),
          panel.grid = element_blank(),
          panel.spacing = unit(1, "lines"),
          panel.border = element_rect(fill = NA, color = 'grey50'),
          axis.title.y = element_text(size = 19, vjust = 2.75),
          axis.title.x = element_text(size = 19, vjust = -2.25),
          axis.text.y = element_text(size = 14, color = 'grey30'),
          axis.text.x = element_text(size = 14, color = 'grey30'),
          axis.line.y = element_line(color = 'grey30'),
          axis.line.x = element_line(color = 'grey30'),
          axis.ticks = element_line(color = 'grey30'),
          legend.position = 'none',
          legend.justification = 'center',
          legend.direction = 'vertical',
          legend.title = element_blank(),
          legend.text = element_text(size = 16),
          legend.key.size = unit(0.5, 'lines'),
          legend.key.height = unit(0.2, "inch"),
          legend.key.width = unit(0.4, "inch")
    ) +
    labs(x = 'False positive rate', y = 'True positive rate') +
    scale_color_manual(values = c(color_cadd, color_dann, color_fathmm, 
                                  color_fitcons, color_hal, 
                                  color_linsight, color_spanr)) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
    scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
    guides(colour = guide_legend(override.aes = list(size = 3)))
    
ggsave(paste0(dir, 'exac_fig4E_exac_roc_curves_type_no_legend', plot_format_main), 
       height = 4, width = 10, units = 'in', dpi = hi_res)

# Legend
ggplot(roc_curve_info_main, aes(false_positive_rate, true_positive_rate)) + 
  geom_line(aes(color = method), size = 1.25) +
  geom_abline(intercept = 0, linetype = 'dashed', color = 'grey50') +
  facet_grid( ~ type) +
  theme(strip.text = element_text(size = 18.5),
        strip.background = element_rect(fill = "#E8E8E8", color = "white"),
        panel.grid = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.border = element_rect(fill = NA, color = 'grey50'),
        axis.title.y = element_text(size = 19, vjust = 2.75),
        axis.title.x = element_text(size = 19, vjust = -2.75),
        axis.text.y = element_text(size = 14, color = 'grey20'),
        axis.text.x = element_text(size = 14, color = 'grey20'),
        axis.line.y = element_line(color = 'grey30'),
        axis.line.x = element_line(color = 'grey30'),
        axis.ticks = element_line(color = 'grey30'),
        legend.position = 'right',
        legend.justification = 'center',
        legend.direction = 'vertical',
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, 'lines'),
        legend.key.height = unit(0.25, "inch"),
        legend.key.width = unit(0.4, "inch")
  ) +
  labs(x = 'False positive rate', y = 'True positive rate') +
  scale_color_manual(values = c(color_cadd, color_dann, color_fathmm, 
                                color_fitcons, color_hal, 
                                color_linsight, color_spanr)) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  guides(colour = guide_legend(override.aes = list(size = 3)))


ggsave(paste0(dir, 'exac_fig4E_exac_roc_curves_type_legend', plot_format_main),
       height = 4, width = 12, units = 'in', dpi = hi_res)
