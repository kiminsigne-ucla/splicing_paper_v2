###############################################################################
# set-up
###############################################################################
load_pkgs <- function(pkgs){
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, 'Package'])]
  if(length(new_pkgs)) install.packages(new_pkgs)
  for(pkg in pkgs){
    suppressWarnings(suppressMessages(library(pkg, character.only = T)))
  }
}

pkgs <- c('dplyr', 'tidyr', 'ggplot2', 'cowplot', 'gridExtra', 'grid')
load_pkgs(pkgs)

options(stringsAsFactors = F, warn = -1, warnings = -1)

plot_format_main <- '.tiff'
plot_format <- '.png'
hi_res <- 600
lo_res <- 300
jitter_alpha <- 0.50

# custom color palette
source("../color_palette.R")
# specify new color palette
steps <- c("blue2", "cyan", "white", "yellow", "red2")
pal <- color.palette(steps, c(160, 1, 1, 160), space = "rgb")

###############################################################################
# Read in data
###############################################################################
data <- read.table('../../processed_data/splicemod/splicemod_data_clean.txt', 
                   sep = '\t', header = T, 
                   colClasses = c('sub_id' = 'character')) %>% 
  filter(rep_quality == 'high')

# read in re-scored reference file
updated_ref <- read.csv('../../ref/splicemod/splicemod_ref_rescored.txt',
                        header = T, sep = '\t',
                        colClasses = c('sub_id' = 'character'))

data <- data %>% 
  left_join(select(updated_ref, id, exon_seq:correct_don_score), by = 'id')

# calculate change in donor/acceptor score between mutant and wild-type 
data <- data %>% 
  group_by(ensembl_id) %>% 
  filter(any(sub_id == '000')) %>% 
  mutate(correct_acc_score_nat = correct_acc_score[sub_id == '000'],
         correct_don_score_nat = correct_don_score[sub_id == '000'],
         delta_acc_score = correct_acc_score - correct_acc_score_nat,
         delta_don_score = correct_don_score - correct_don_score_nat,
         # calculate fold change in score relative to wild-type
         don_score_fold_change = (correct_don_score - correct_don_score_nat) / 
           abs(correct_don_score_nat),
         acc_score_fold_change = (correct_acc_score - correct_acc_score_nat) / 
           abs(correct_acc_score_nat)) %>%
  ungroup()

# compare ESE changes to random
data %>% 
    mutate(category = ifelse(grepl('ESE', category), 'perturb_ESE', category),
           category = ifelse(grepl('rnd_exon', category), 'rnd_exon', category)) %>% 
    filter(category == 'perturb_ESE' | category == 'rnd_exon') %>% 
    t.test(dpsi_smn1 ~ category, data = ., alternative = 'less')

# compare ESS changes to random
data %>% 
    mutate(category = ifelse(grepl('ESS', category), 'perturb_ESS', category),
           category = ifelse(grepl('rnd_exon', category), 'rnd_exon', category)) %>% 
    filter(category == 'perturb_ESS' | category == 'rnd_exon') %>% 
    t.test(dpsi_smn1 ~ category, data = ., alternative = 'greater')

# strongest ESE
data %>% 
    mutate(category = ifelse(grepl('rnd_exon', category), 'rnd_exon', category)) %>% 
    filter(category == 'clst_Ke2011_ESE' | category == 'rnd_exon') %>% 
    t.test(dpsi_smn1 ~ category, data = ., alternative = 'less')

# strongest ESS
data %>% 
    mutate(category = ifelse(grepl('rnd_exon', category), 'rnd_exon', category)) %>% 
    filter(category == 'clst_Ke2011_ESS' | category == 'rnd_exon') %>% 
    t.test(dpsi_smn1 ~ category, data = ., alternative = 'greater')

# intron changes
data %>% 
    mutate(category = ifelse(grepl('ICS', category), 'intron', category),
           category = ifelse(grepl('rnd_intron', category), 'rnd_intron', category)) %>% 
    filter(category == 'intron' | category == 'rnd_intron') %>% 
    t.test(dpsi_smn1 ~ category, data = .)

# label variants as exonic or intronic
tmp <- data %>% 
    separate(loc, into = c('loc_start', 'loc_end', sep = ':', convert = T)) %>% 
    mutate(exon_start = ,
           exon_end = ,
           exon_overlap = ifelse(intersect(seq(loc_start, loc_end), seq(exon_start, exon_end))))

###############################################################################
# MaxEnt 
###############################################################################

# MaxEnt: splice acceptor and donor score fold-change
# SMN1 intron backbone
cor.test(data$delta_acc_score, data$dpsi_smn1)
gg <- data %>% 
  gather(key = 'splice_site', value = 'fold_change', 
         don_score_fold_change, acc_score_fold_change) %>% 
  filter(fold_change != 0) %>% 
  mutate(fold_change_bin = cut(fold_change, c(-184, -2, -1, 0, 1))) %>%
  filter(!is.na(fold_change_bin), seq_type == 'mut') %>% 
  # reorder binned intervals
  mutate(fold_change_bin = factor(fold_change_bin,
                                  levels = c('(-184,-2]', '(-2,-1]',
                                             '(-1,0]', '(0,1]'))) %>% 
  # select(acc_score_fold_change, don_score_fold_change, splice_site, 
  # fold_change_bin)
  ggplot(aes(fold_change_bin, dpsi_smn1)) + 
  geom_jitter(alpha = jitter_alpha, aes(color = nat_index_smn1)) + 
  geom_boxplot(alpha = 0) +
  facet_grid(~ splice_site,
             labeller = 
               as_labeller(c('acc_score_fold_change' = 'Splice Acceptor',
                             'don_score_fold_change' = 'Splice Donor'))) +
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = seq(0, 1, by = 0.25), 
                         colors = pal(321)) +
  scale_x_discrete(labels = c('<= -2', '(-2, -1]', '(-1, -0)', '(0, 1]')) +
  labs(x = 'MaxEnt score (fold change)', 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  theme(strip.text = element_text(size = 16),
        strip.background = element_rect(fill = "#E8E8E8", color = "white"),
        axis.title.x = element_text(size = 16, vjust = -2), 
        axis.title.y = element_text(size = 18), 
        axis.text = element_text(size = 12, color = "grey20"), 
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        text = element_text(size = 12),
        plot.margin = unit(c(0,0,3,0),"mm"),
        legend.position = 'none') 

ggggsave(paste0('../../figs/splicemod/smn1/splicemod_smn1_don_acc', 
              plot_format_main), gg, 
       width = 4.6, height = 3, dpi = hi_res, scale = 1.3)
ggsave(paste0('../../figs/splicemod/smn1/splicemod_smn1_don_acc', 
              plot_format), gg, 
       width = 4.6, height = 3, dpi = lo_res, scale = 1.3)

# MaxEnt: splice acceptor and donor score fold-change)
# DHFR intron backbone

gg <- data %>% 
  gather(key = 'splice_site', value = 'fold_change', 
         don_score_fold_change, acc_score_fold_change) %>% 
  filter(fold_change != 0) %>% 
  mutate(fold_change_bin = cut(fold_change, c(-184, -2, -1, 0, 1))) %>%
  filter(!is.na(fold_change_bin), seq_type == 'mut') %>% 
  # reorder binned intervals
  mutate(fold_change_bin = factor(fold_change_bin,
                                  levels = c('(-184,-2]', '(-2,-1]',
                                             '(-1,0]', '(0,1]'))) %>% 
  # select(acc_score_fold_change, don_score_fold_change, splice_site,
  # fold_change_bin)
  ggplot(aes(fold_change_bin, dpsi_dhfr)) + 
  geom_jitter(alpha = jitter_alpha, aes(color = nat_index_dhfr)) + 
  geom_boxplot(alpha = 0) +
  facet_grid(~ splice_site,
             labeller = 
               as_labeller(c('acc_score_fold_change' = 'Splice Acceptor',
                             'don_score_fold_change' = 'Splice Donor'))) +
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = seq(0, 1, by = 0.25), 
                         colors = pal(321)) +
  scale_x_discrete(labels = c('<= -2', '(-2, -1]', '(-1, -0)', '(0, 1]')) +
  labs(x = 'MaxEnt score (fold change)', 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  theme(strip.text = element_text(size = 16),
        strip.background = element_rect(fill = "#E8E8E8", color = "white"),
        axis.title.x = element_text(size = 16, vjust = -2), 
        axis.title.y = element_text(size = 18), 
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        axis.text = element_text(size = 12, color = "grey20"), 
        text = element_text(size = 12),
        plot.margin = unit(c(0,0,3,0),"mm"),
        legend.position = 'none') 

ggsave(paste0('../../figs/splicemod/dhfr/splicemod_dhfr_don_acc', 
              plot_format), 
       gg, width = 4.6, height = 3, dpi = lo_res, scale = 1.3)

###############################################################################
# Splicemod categories
###############################################################################
tmp <- data %>% 
    mutate(test_cat = case_when(.$category == 'rmv_Ke2011_ESE' ~ T,
                                .$category == 'clst_Ke2011_ESE' ~ T,
                                TRUE ~ F))
t.test(dpsi_smn1 ~ test_cat, tmp)

esr_categories <- c('rmv_Ke2011_ESE', 'clst_Ke2011_ESE', 'rmv_Ke2011_ESS', 
                    'clst_Ke2011_ESS')
esr_labels <- c('weaken ESEs', 'destroy strongest ESE', 'weaken ESSs', 
                'destroy strongest ESS')

splice_site_categories <- c( 'weak_spl_a', 'weak_spl_d', 'p_weak_spl', 
                             'no_spl_a', 'no_spl_d',
                             'same_splice_a', 'same_splice_d',
                             'rmv_me_splice_acceptor', 'rmv_me_splice_donor', 
                             'csplice_a', 'csplice_d')
splice_site_labels <- c('weaken acceptor', 'weaken donor', 
                        'weaken donor + acceptor', 
                        'destroy acceptor', 'destroy donor',
                        'same score acceptor', 'same score donor',
                        'weaken spurious acceptor', 'weaken spurious donor',
                        'destroy spurious acceptor', 'destroy spurious donor')

intron_categories <- c('clst_Vlkr07_AICS', 'clst_Vlkr07_DICS', 
                       'rmv_Vlkr07_AICS', 'rmv_Vlkr07_DICS')
intron_labels <- c('weaken intronic conserved (acceptor side)', 
                   'weaken intronic conserved (donor side)', 
                   'destroy intronic conserved (acceptor side)', 
                   'destroy intronic conserved (donor side)')

random_exon_categories <- c('rnd_exon_1nt', 'rnd_exon_2nt', 
                            'rnd_exon_3nt', 'rnd_exon_5nt', 
                            'aggr_exon')
random_exon_labels <- c('random 1nt exon', 'random 2nt exon', 
                        'random 3nt exon', 'random 5nt exon',
                        'aggressive exon (only syn. mut.)')

random_intron_categories <- c('rnd_intron_1nt', 'rnd_intron_2nt', 
                              'rnd_intron_3nt', 'rnd_intron_5nt', 
                              'aggr_intron', 'p_aggr_intr', 'aggr_both')
random_intron_labels <- c('random 1nt intron', 'random 2nt intron', 
                          'random 3nt intron', 'random 5nt intron',
                          'aggressive intron', 'aggr. + random intron', 
                          'aggr. intron + exon')

other_categories <- c('cnsrv_1nt', 'cnsrv_3nt', 
                      'RBPmats', 'variation')
other_labels <- c('conserved 1nt', 'conservered 3nt', 
                  'destroy RBP motifs', 'dbSNPs')

# SMN1 intron backbone
gg <- data %>% 
  filter(seq_type == 'mut') %>% 
  mutate(category_fctr = factor(category,
                                levels = c(splice_site_categories, 
                                           esr_categories, 
                                           random_exon_categories,
                                           intron_categories, 
                                           random_intron_categories, 
                                           other_categories))) %>% 
  ggplot(aes(x = category_fctr, y = dpsi_smn1)) + 
  geom_jitter(alpha = jitter_alpha, aes(color = nat_index_smn1)) + 
  geom_boxplot(alpha = 0) +
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = seq(0, 1, by = 0.25), colors = pal(321)) + 
  scale_x_discrete(labels = c(splice_site_labels, esr_labels, 
                              random_exon_labels, intron_labels, 
                              random_intron_labels, other_labels)) +
  theme(axis.title.y = element_text(size = 18, vjust = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "grey20"),
        axis.text.y = element_text(size = 10, color = "grey20"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        legend.position = 'none') +
  labs(x = '', y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "]))

gg_no_x_axis <- gg + theme(axis.text.x = element_blank())

ggsave(paste0('../../figs/splicemod/smn1/', 
              'splicemod_smn1_all_categories_no_x_axis', plot_format_main), 
       gg_no_x_axis, width = 12, height = 3, dpi = hi_res)
ggsave(paste0('../../figs/splicemod/smn1/', 
              'splicemod_smn1_all_categories_no_x_axis', plot_format), 
       gg_no_x_axis, width = 12, height = 3, dpi = lo_res)
ggsave(paste0('../../figs/splicemod/smn1/',
              'splicemod_smn1_all_categories', plot_format), 
       gg, width = 12, height = 5, dpi = lo_res)

# DHFR intron backbone
gg <- data %>% 
  filter(seq_type == 'mut') %>% 
  mutate(category_fctr = factor(category,
                                levels = c(splice_site_categories, 
                                           esr_categories, 
                                           random_exon_categories,
                                           intron_categories, 
                                           random_intron_categories, 
                                           other_categories))) %>% 
  ggplot(aes(x = category_fctr, y = dpsi_dhfr)) + 
  geom_jitter(alpha = jitter_alpha, aes(color = nat_index_dhfr)) + 
  geom_boxplot(alpha = 0) +
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = seq(0, 1, by = 0.25), colors = pal(321)) + 
  scale_x_discrete(labels = c(splice_site_labels, esr_labels, 
                              random_exon_labels, intron_labels, 
                              random_intron_labels, other_labels)) +
  theme(axis.title.y = element_text(size = 18, vjust = -1), 
        axis.text.x = element_text(angle = 45, hjust = 1, color = "grey20"), 
        axis.text.y = element_text(size = 10, color = "grey20"), 
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        legend.position = 'none') +
  labs(x = '', y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "]))

gg_no_x_axis <- gg + theme(axis.text.x = element_blank())

ggsave(paste0('../../figs/splicemod/dhfr/',
              'splicemod_dhfr_all_categories_no_x_axis', plot_format_main), 
       gg_no_x_axis, width = 12, height = 3, dpi = lo_res)
ggsave(paste0('../../figs/splicemod/dhfr/',
              'splicemod_dhfr_all_categories_no_x_axis', plot_format), 
       gg_no_x_axis, width = 12, height = 3, dpi = lo_res)
ggsave(paste0('../../figs/splicemod/dhfr/',
              'splicemod_dhfr_all_categories', plot_format), 
       gg, width = 12, height = 5, dpi = lo_res)

###############################################################################
# Exonic motifs
###############################################################################

# calculate change in average HAL score
data <- data %>% 
  rename('avg_HAL_score' = 'avg_exon_effect_score') %>% 
  group_by(ensembl_id) %>% 
  mutate(delta_avg_HAL_score = 
           avg_HAL_score - avg_HAL_score[sub_id == '000']) %>% 
  ungroup()

################################
# HAL (hexamer additive linear)
# PMID: 26496609
################################
data <- data %>%
  mutate(HAL_bin = case_when(.$delta_avg_HAL_score < 0 ~ 'down',
                             .$delta_avg_HAL_score > 0 ~ 'up',
                             .$delta_avg_HAL_score == 0 ~ 'same'))

# SMN1 intron backbone
gg <- data %>%
  filter(HAL_bin != 'same', seq_type == 'mut') %>%
  ggplot(aes(HAL_bin, dpsi_smn1)) + geom_jitter(alpha = jitter_alpha, 
                                                aes(color = nat_index_smn1)) + 
  geom_violin(alpha = 0, color = "grey35") +
  scale_colour_gradientn(limits = c(-0.005, 1), breaks = seq(0, 1, by = 0.25), 
                         colors = pal(321)) + 
  labs(x = expression(paste(Delta, ' avg. exonic hexamer score')), 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  ggsignif::geom_signif(comparisons = list(c('down', 'up')),
                        test = 't.test', map_signif_level = T, 
                        tip_length = 0) +
  stat_summary(fun.y = mean, geom = "point", size = 2, color = "black") +
  theme( axis.title.x = element_text(size = 16, vjust = -2, hjust = 1), 
         axis.title.y = element_text(size = 18, vjust = 40),
         axis.ticks.x = element_blank(),
         axis.ticks.y = element_line(color = "grey50"),
         axis.line.x = element_line(color = "grey50"),
         axis.line.y = element_line(color = "grey50"),
         axis.text.x = element_text(size = 18, color = "grey20"), 
         axis.text.y = element_text(color = "grey20"), 
         plot.margin = unit(c(1,1,1,1.5),"mm"),
         legend.position = 'none')

ggsave(paste0('../../figs/splicemod/smn1/splicemod_smn1_hal', 
              plot_format_main), gg, 
       width = 2.5, height = 3, dpi = hi_res, scale = 1.3)

# DHFR intron backbone
gg <- data %>%
  filter(HAL_bin != 'same', seq_type == 'mut') %>%
  ggplot(aes(HAL_bin, dpsi_dhfr)) + geom_jitter(alpha = jitter_alpha, 
                                                aes(color = nat_index_dhfr)) + 
  geom_violin(alpha = 0, color = "grey35") +
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = seq(0, 1, by = 0.25), colors = pal(321)) + 
  labs(x = expression(paste(Delta, ' avg. exonic hexamer score (HAL)')), 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  ggsignif::geom_signif(comparisons = list(c('down', 'up')),
                        test = 't.test', map_signif_level = T, 
                        tip_length = 0) +
  stat_summary(fun.y = mean, geom = "point", size = 2, color = "black") +
  theme(axis.title.x = element_text(size = 14, vjust = -2, hjust = 1), 
        axis.title.y = element_text(size = 18, vjust = 40),
        axis.text.x = element_text(size = 18, color = "grey20"),  
        axis.text.y = element_text(color = "grey20"), 
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        plot.margin = unit(c(1,1,1,1.5),"mm"),
        legend.position = 'none')

ggsave(paste0('../../figs/splicemod/dhfr/splicemod_dhfr_hal', 
              plot_format), gg, 
       width = 2.5, height = 3, dpi = hi_res, scale = 1.3)

########################
# Ke 2011 ESE/ESS 6-mer
# PMID: 21659425
########################
data <- data %>% 
  group_by(ensembl_id) %>% 
  mutate(delta_Ke2011_avg_score = 
           Ke2011_avg_score - Ke2011_avg_score[sub_id == '000'],
         delta_Ke2011_ESE_avg_score = 
           Ke2011_ESE_avg_score - Ke2011_ESE_avg_score[sub_id == '000'],
         delta_Ke2011_ESS_avg_score = 
           Ke2011_ESS_avg_score - Ke2011_ESS_avg_score[sub_id == '000']) %>% 
  ungroup()

data <- data %>%
  mutate(Ke_bin = case_when(.$delta_Ke2011_avg_score < 0 ~ 'down',
                            .$delta_Ke2011_avg_score > 0 ~ 'up',
                            .$delta_Ke2011_avg_score == 0 ~ 'same'),
         Ke_ESE_bin = case_when(.$delta_Ke2011_ESE_avg_score < 0 ~ 'down',
                                .$delta_Ke2011_ESE_avg_score > 0 ~ 'up',
                                .$delta_Ke2011_ESE_avg_score == 0 ~ 'same'),
         Ke_ESS_bin = case_when(.$delta_Ke2011_ESS_avg_score < 0 ~ 'down',
                                .$delta_Ke2011_ESS_avg_score > 0 ~ 'up',
                                .$delta_Ke2011_ESS_avg_score == 0 ~ 'same'))

# SMN1 intron backbone
gg <- data %>%
  filter(Ke_bin != 'same', seq_type == 'mut') %>%
  ggplot(aes(Ke_bin, dpsi_smn1)) + geom_jitter(alpha = jitter_alpha, 
                                               aes(color = nat_index_smn1)) + 
  geom_violin(alpha = 0, color = "grey35") +
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = seq(0, 1, by = 0.25), colors = pal(321)) +
  labs(x = expression(paste(Delta, ' avg. exonic hexamer score (Ke)')), 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  ggsignif::geom_signif(comparisons = list(c('down', 'up')),
                        test = 't.test', 
                        map_signif_level = T, tip_length = 0) +
  stat_summary(fun.y = mean, geom = "point", size = 2, color = "black") +
  theme(axis.title.x = element_text(size = 15, vjust = -2, hjust = 1), 
        axis.title.y = element_text(size = 18, vjust = 40),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        axis.text.x = element_text(size = 18, color = "grey20"),  
        axis.text.y = element_text(color = "grey20"),
        plot.margin = unit(c(1,1,1,1),"mm"),
        legend.position = 'none')

ggsave(paste0('../../figs/splicemod/smn1/splicemod_smn1_Ke11', plot_format),
       gg, width = 2.5, height = 3, dpi = lo_res, scale = 1.3)

# DHFR intron backbone
gg <- data %>%
  filter(Ke_bin != 'same', seq_type == 'mut') %>%
  ggplot(aes(Ke_bin, dpsi_dhfr)) + geom_jitter(alpha = jitter_alpha, 
                                               aes(color = nat_index_dhfr)) + 
  geom_violin(alpha = 0, color = "grey35") +
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = seq(0, 1, by = 0.25), colors = pal(321)) +
  labs(x = expression(paste(Delta, ' avg. exonic hexamer score (Ke)')), 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  ggsignif::geom_signif(comparisons = list(c('down', 'up')),
                        test = 't.test', 
                        map_signif_level = T, tip_length = 0) +
  stat_summary(fun.y = mean, geom = "point", size = 2, color = "black") +
  theme(axis.title.x = element_text(size = 15, vjust = -2, hjust = 1), 
        axis.title.y = element_text(size = 18, vjust = 40),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "grey50"),
        axis.line.x = element_line(color = "grey50"),
        axis.line.y = element_line(color = "grey50"),
        axis.text.x = element_text(size = 18, color = "grey20"),  
        axis.text.y = element_text(color = "grey20"),
        plot.margin = unit(c(1,1,1,1),"mm"),
        legend.position = 'none')

ggsave(paste0('../../figs/splicemod/dhfr/splicemod_dhfr_Ke11', plot_format),
       gg, width = 2.5, height = 3, dpi = lo_res, scale = 1.3)

########################
# Legend for Figure
########################
# grab legend
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

gg <- data %>%
  filter(acc_score_fold_change != 0) %>%
  mutate(acc_fold_change_bin = 
           cut(acc_score_fold_change, c(-184, -2, -1, 0, 1))) %>%
  filter(!is.na(acc_fold_change_bin)) %>%
  ggplot(aes(acc_fold_change_bin, dpsi_dhfr)) + 
  geom_jitter(alpha = jitter_alpha, aes(color = nat_index_dhfr)) + 
  scale_x_discrete(labels = c('<= -2', '(-2, -1]', '(-1, -0)', '(0, 1]')) +
  geom_boxplot(alpha = 0) + 
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = c(0.2, 0.4, 0.6, 0.8, 1), 
                         colors = pal(321)) +
  labs(x = '', y = '',
       color = expression(index["WT "])) +
  theme(axis.text = element_text(size = 12, color = "grey20"), 
        text = element_text(size = 12),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1.8, "cm"),
        legend.text = element_text(size = 24, color = "grey20"),
        legend.title = element_text(size = 32),
        plot.margin = unit(c(0,0,0,0),"mm"))

legend <- g_legend(gg)
tiff(paste0('../../figs/splicemod/both/legend', plot_format_main), 
     width = 40, height = 125, units = 'mm', res = hi_res)
grid.newpage()
grid.draw(legend)
dev.off()

legend <- g_legend(gg)
png(paste0('../../figs/splicemod/both/legend', plot_format), 
     width = 40, height = 125, units = 'mm', res = lo_res)
grid.newpage()
grid.draw(legend)
dev.off()

