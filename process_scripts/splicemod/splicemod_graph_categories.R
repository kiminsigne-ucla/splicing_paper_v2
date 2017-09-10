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
plot_format <- '.png'

axis_title_x <- 16
axis_title_y <- 16
axis_text <- 12
general_text <- 12
jitter_alpha <- 0.50

# custom color palette
source("../color_palette.R")
# Specify new color palette
steps <- c("blue2", "cyan", "white", "yellow", "red2")
pal <- color.palette(steps, c(160,1,1,160), space = "rgb")

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



# calculate change in donor/acceptor score between mutant and natural
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

###############################################################################
# MaxEnt 
###############################################################################
# MaxEnt splice donor score fold-change
# SMN1

gg <- data %>%
  filter(don_score_fold_change != 0) %>%
  mutate(don_fold_change_bin = cut(don_score_fold_change, c(-12, -2, -1, 0, 1))) %>%
  filter(!is.na(don_fold_change_bin)) %>%
  ggplot(aes(don_fold_change_bin, dpsi_smn1)) + 
  geom_jitter(alpha = jitter_alpha, aes(color = nat_index_smn1)) + 
  scale_x_discrete(labels = c('<= -2', '(-2, -1]', '(-1, -0)', '(0, 1]')) +
  geom_boxplot(alpha = 0) + 
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = seq(0, 1, by = 0.25), 
                         colors = pal(321)) +
  labs(x = 'fold change MaxEnt splice donor score', 
       y = expression(paste(Delta, ' inclusion index (SMN1)')),
       color = expression(index["WT "])) +
  theme(axis.text = element_text(size = axis_text), 
        text = element_text(size = general_text))
       
# grab legend
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
legend <- g_legend(gg)
tiff(paste0('../../figs/splicemod/nat_index_legend', plot_format), width = 350, 350)
grid.newpage()
grid.draw(legend)
dev.off()

gg_no_legend <- gg + theme(legend.position = 'none')

ggsave(paste0('../../figs/splicemod/smn1/splicemod_smn1_donor_fc', plot_format), 
       gg_no_legend, width = 2.4, height = 3, dpi = 300, scale = 1.3)

       
# MaxEnt splice donor score fold-change
# SMN1

gg <- data %>%
  filter(acc_score_fold_change != 0) %>%
  mutate(acc_fold_change_bin = cut(acc_score_fold_change, c(-184, -2, -1, 0, 1))) %>%
  filter(!is.na(acc_fold_change_bin)) %>%
  ggplot(aes(acc_fold_change_bin, dpsi_smn1)) + geom_jitter(alpha = jitter_alpha , aes(color = nat_index_smn1)) + 
  scale_x_discrete(labels = c('<= -2', '(-2, -1]', '(-1, -0)', '(0, 1]')) +
  geom_boxplot(alpha = 0) + 
  scale_colour_gradientn(limits = c(-0.005, 1), breaks = seq(0, 1, by = 0.25), colors = pal(321)) +
  labs(x = 'fold change\nMaxEnt splice acceptor score', 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  theme(axis.text = element_text(size = axis_text), text = element_text(size = general_text),
        legend.position = 'none')

ggsave(paste0('../../figs/splicemod/smn1/splicemod_smn1_acceptor_fc', plot_format), 
       gg, width = 2.4, height = 3, dpi = 300, scale = 1.3)


# MaxEnt both splice acceptor and donor score fold-change
# SMN1

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
  # select(acc_score_fold_change, don_score_fold_change, splice_site, fold_change_bin)
  ggplot(aes(fold_change_bin, dpsi_smn1)) + 
  geom_jitter(alpha = jitter_alpha, aes(color = nat_index_smn1)) + 
  geom_boxplot(alpha = 0) +
  facet_grid(~ splice_site,
             labeller = as_labeller(c('acc_score_fold_change' = 'Splice Acceptor',
                                      'don_score_fold_change' = 'Splice Donor'))) +
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = seq(0, 1, by = 0.25), 
                         colors = pal(321)) +
  scale_x_discrete(labels = c('<= -2', '(-2, -1]', '(-1, -0)', '(0, 1]')) +
  labs(x = 'MaxEnt score (fold change)', 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  theme(strip.background = element_rect(fill="lightgrey"),
        strip.text = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.title.x = element_text(size = 16), 
        axis.text = element_text(size = axis_text), 
        text = element_text(size = general_text),
        legend.position = 'none') 

ggsave(paste0('../../figs/splicemod/smn1/splicemod_smn1_don_acc', plot_format), 
       gg, width = 4.6, height = 3, dpi = 600, scale = 1.3)


# MaxEnt splice donor score fold-change
# DHFR

gg <- data %>%
  filter(don_score_fold_change != 0) %>%
  mutate(don_fold_change_bin = cut(don_score_fold_change, c(-12, -2, -1, 0, 1))) %>%
  filter(!is.na(don_fold_change_bin)) %>%
  ggplot(aes(don_fold_change_bin, dpsi_smn1)) + 
  geom_jitter(alpha = jitter_alpha , aes(color = nat_index_dhfr)) + 
  scale_x_discrete(labels = c('<= -2', '(-2, -1]', '(-1, -0)', '(0, 1]')) +
  geom_boxplot(alpha = 0) + 
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = seq(0, 1, by = 0.25), 
                         colors = pal(321)) +
  labs(x = 'fold change MaxEnt splice donor score', 
       y = expression(paste(Delta, ' inclusion index (SMN1)')),
       color = expression(index["WT "])) +
  theme(axis.text = element_text(size = axis_text), 
        text = element_text(size = general_text),
        legend.position = 'none')

ggsave(paste0('../../figs/splicemod/dhfr/splicemod_dhfr_donor_fc', plot_format), 
       gg, width = 2.4, height = 3, dpi = 300, scale = 1.3)

# Legend

gg <- data %>%
  filter(acc_score_fold_change != 0) %>%
  mutate(acc_fold_change_bin = cut(acc_score_fold_change, c(-184, -2, -1, 0, 1))) %>%
  filter(!is.na(acc_fold_change_bin)) %>%
  ggplot(aes(acc_fold_change_bin, dpsi_smn1)) + 
  geom_jitter(alpha = jitter_alpha, aes(color = nat_index_smn1)) + 
  scale_x_discrete(labels = c('<= -2', '(-2, -1]', '(-1, -0)', '(0, 1]')) +
  geom_boxplot(alpha = 0) + 
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = seq(0, 1, by = 0.25), 
                         colors = pal(321)) +
  labs(x = 'fold change MaxEnt splice acceptor score',
       y = expression(paste(Delta, ' inclusion index (SMN1)')),
       color = expression(index["WT "])) +
  theme(axis.text = element_text(size = axis_text), 
        text = element_text(size = general_text),
        legend.key.width=unit(1,"cm"),
        legend.key.height=unit(1.8,"cm"),
        legend.text=element_text(size=28),
        legend.title=element_text(size=32))

ggsave(paste0('../../figs/splicemod/legend', plot_format), 
       gg, width = 3, height = 6, dpi = 300, scale = 1.3)


# MaxEnt splice acceptor score fold-change
# DHFR

gg <- data %>%
  filter(acc_score_fold_change != 0) %>%
  mutate(acc_fold_change_bin = cut(acc_score_fold_change, c(-184, -2, -1, 0, 1))) %>%
  filter(!is.na(acc_fold_change_bin)) %>%
  ggplot(aes(acc_fold_change_bin, dpsi_dhfr)) + 
  geom_jitter(alpha = jitter_alpha, aes(color = nat_index_dhfr)) + 
  scale_x_discrete(labels = c('<= -2', '(-2, -1]', '(-1, -0)', '(0, 1]')) +
  geom_boxplot(alpha = 0) + 
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = seq(0, 1, by = 0.25), 
                         colors = pal(321)) +
  labs(x = 'fold change MaxEnt splice acceptor score',
       y = expression(paste(Delta, ' inclusion index (DHFR)')),
       color = expression(index["WT "])) +
  theme(axis.text = element_text(size = axis_text), 
        text = element_text(size = general_text),
        legend.position = 'none')

ggsave(paste0('../../figs/splicemod/smn1/splicemod_dhfr_acceptor_fc', plot_format), 
       gg, width = 3, height = 3, dpi = 300, scale = 1.3)


# MaxEnt both splice acceptor and donor score fold-change
# DHFR

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
  # select(acc_score_fold_change, don_score_fold_change, splice_site, fold_change_bin)
  ggplot(aes(fold_change_bin, dpsi_dhfr)) + 
  geom_jitter(alpha = jitter_alpha, aes(color = nat_index_dhfr)) + 
  geom_boxplot(alpha = 0) +
  facet_grid(~ splice_site,
             labeller = as_labeller(c('acc_score_fold_change' = 'Splice Acceptor',
                                      'don_score_fold_change' = 'Splice Donor'))) +
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = seq(0, 1, by = 0.25), 
                         colors = pal(321)) +
  scale_x_discrete(labels = c('<= -2', '(-2, -1]', '(-1, -0)', '(0, 1]')) +
  labs(x = 'MaxEnt score (fold change)', 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  theme(strip.text = element_text(size = 16),
        axis.title.y = element_text(size = 16), 
        axis.title.x = element_text(size = 16), 
        axis.text = element_text(size = axis_text), 
        text = element_text(size = general_text),
        legend.position = 'none') 

ggsave(paste0('../../figs/splicemod/dhfr/splicemod_dhfr_don_acc', plot_format), 
       gg, width = 4.6, height = 3, dpi = 300, scale = 1.3)
              

# facet (SMN1)
data %>% 
  gather(key = 'splice_site', value = 'fold_change', 
         don_score_fold_change, acc_score_fold_change) %>% 
  filter(fold_change != 0) %>% 
  mutate(fold_change_bin = cut(fold_change, c(-184, -2, -1, 0, 1))) %>%
  filter(!is.na(fold_change_bin), seq_type == 'mut') %>% 
  # reorder binned intervals
  mutate(fold_change_bin = factor(fold_change_bin,
                                  levels = c('(-184,-2]', '(-2,-1]',
                                             '(-1,0]', '(0,1]'))) %>% 
  # select(acc_score_fold_change, don_score_fold_change, splice_site, fold_change_bin)
  ggplot(aes(fold_change_bin, dpsi_smn1)) + 
  geom_jitter(alpha = jitter_alpha, aes(color = nat_index_smn1)) + 
  geom_boxplot(alpha = 0) +
  facet_grid(~ splice_site,
             labeller = as_labeller(c('acc_score_fold_change' = 'SA',
                                      'don_score_fold_change' = 'SD'))) +
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = seq(0, 1, by = 0.25), 
                         colors = pal(321)) +
  scale_x_discrete(labels = c('<= -2', '(-2, -1]', '(-1, -0)', '(0, 1]')) +
  labs(x = 'fold-change MaxEnt score', 
       y = expression(paste(Delta, ' inclusion index (SMN1)')),
       color = expression(index["WT "])) +
  theme(axis.text = element_text(size = axis_text), 
        text = element_text(size = general_text),
        legend.position = 'none')

# facet (DHFR)
data %>% 
  gather(key = 'splice_site', value = 'fold_change', 
         don_score_fold_change, acc_score_fold_change) %>% 
  filter(fold_change != 0) %>% 
  mutate(fold_change_bin = cut(fold_change, c(-184, -2, -1, 0, 1))) %>%
  filter(!is.na(fold_change_bin), seq_type == 'mut') %>% 
  # reorder binned intervals
  mutate(fold_change_bin = factor(fold_change_bin,
                                  levels = c('(-184,-2]', '(-2,-1]',
                                             '(-1,0]', '(0,1]'))) %>% 
  # select(acc_score_fold_change, don_score_fold_change, splice_site, fold_change_bin)
  ggplot(aes(fold_change_bin, dpsi_dhfr)) + 
  geom_jitter(alpha = jitter_alpha, aes(color = nat_index_dhfr)) + 
  geom_boxplot(alpha = 0) +
  facet_grid(~ splice_site,
             labeller = as_labeller(c('acc_score_fold_change' = 'SA',
                                      'don_score_fold_change' = 'SD'))) +
  scale_colour_gradientn(limits = c(-0.005, 1), 
                         breaks = seq(0, 1, by = 0.25), 
                         colors = pal(321)) +
  scale_x_discrete(labels = c('<= -2', '(-2, -1]', '(-1, -0)', '(0, 1]')) +
  labs(x = 'fold-change MaxEnt score', 
       y = expression(paste(Delta, ' inclusion index (DHFR)')),
       color = expression(index["WT "])) +
  theme(axis.text = element_text(size = axis_text), 
        text = element_text(size = general_text),
        legend.position = 'none')

###############################################################################
# Splicemod categories
###############################################################################
esr_categories <- c('rmv_Ke2011_ESE', 'clst_Ke2011_ESE', 'rmv_Ke2011_ESS', 'clst_Ke2011_ESS')
esr_labels <- c('weaken ESEs', 'destroy strongest ESE', 'weaken ESSs', 'destroy strongest ESS')

splice_site_categories <- c( 'weak_spl_a', 'weak_spl_d', 'p_weak_spl', 'no_spl_a', 'no_spl_d',
                             'same_splice_a', 'same_splice_d',
                             'rmv_me_splice_acceptor', 'rmv_me_splice_donor', 
                             'csplice_a', 'csplice_d')
splice_site_labels <- c('weaken acceptor', 'weaken donor', 'weaken donor + acceptor', 'destroy acceptor', 'destroy donor',
                        'same score acceptor', 'same score donor',
                        'weaken spurious acceptor', 'weaken spurious donor',
                        'destroy spurious acceptor', 'destroy spurious donor')

intron_categories <- c('clst_Vlkr07_AICS', 'clst_Vlkr07_DICS', 
                       'rmv_Vlkr07_AICS', 'rmv_Vlkr07_DICS')
intron_labels <- c('weaken intronic conserved (acceptor side)', 'weaken intronic conserved (donor side)', 
                   'destroy intronic conserved (acceptor side)', 'destroy intronic conserved (donor side)')

random_intron_categories <- c('rnd_intron_1nt', 'rnd_intron_2nt', 'rnd_intron_3nt', 'rnd_intron_5nt', 
                              'aggr_intron', 'p_aggr_intr', 'aggr_both')
random_intron_labels <- c('random 1nt intron', 'random 2nt intron', 'random 3nt intron', 'random 5nt intron',
                          'aggressive intron', 'aggr. + random intron', 'aggr. intron + exon')

random_exon_categories <- c('rnd_exon_1nt', 'rnd_exon_2nt', 'rnd_exon_3nt', 'rnd_exon_5nt', 
                            'aggr_exon')
random_exon_labels <- c('random 1nt exon', 'random 2nt exon', 'random 3nt exon', 'random 5nt exon',
                        'aggressive exon (only syn. mut.)')

other_categories <- c('cnsrv_1nt', 'cnsrv_3nt', 'RBPmats', 'variation')
other_labels <- c('conserved 1nt', 'conservered 3nt', 'destroy RBP motifs', 'dbSNPs')

# SMN1
gg <- data %>% 
  filter(seq_type == 'mut') %>% 
  mutate(category_fctr = factor(category,
                                levels = c(splice_site_categories, esr_categories, random_exon_categories,
                                           intron_categories, random_intron_categories, other_categories))) %>% 
  ggplot(aes(x = category_fctr, y = dpsi_smn1)) + 
  geom_jitter(alpha = jitter_alpha, aes(color = nat_index_smn1)) + 
  geom_boxplot(alpha = 0) +
  scale_colour_gradientn(limits = c(-0.005, 1), breaks = seq(0, 1, by = 0.25), colors = pal(321)) + 
  scale_x_discrete(labels = c(splice_site_labels, esr_labels, random_exon_labels, intron_labels, random_intron_labels, other_labels)) +
  theme(axis.title.y = element_text(size = axis_title_y), axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = 'none') +
  labs(x = '', y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "]))

gg_no_x_axis <- gg + theme(axis.text.x = element_blank())

ggsave(paste0('../../figs/splicemod/smn1/splicemod_smn1_all_categories_no_x_axis', plot_format), 
       gg_no_x_axis, width = 12, height = 3, dpi = 600)
ggsave(paste0('../../figs/splicemod/smn1/splicemod_smn1_all_categories', plot_format), 
       gg, width = 12, height = 5, dpi = 300)

# DHFR
gg <- data %>% 
  filter(seq_type == 'mut') %>% 
  mutate(category_fctr = factor(category,
                                levels = c(splice_site_categories, esr_categories, random_exon_categories,
                                           intron_categories, random_intron_categories, other_categories))) %>% 
  ggplot(aes(x = category_fctr, y = dpsi_dhfr)) + 
  geom_jitter(alpha = jitter_alpha, aes(color = nat_index_dhfr)) + 
  geom_boxplot(alpha = 0) +
  scale_colour_gradientn(limits = c(-0.005, 1), breaks = seq(0, 1, by = 0.25), colors = pal(321)) + 
  scale_x_discrete(labels = c(splice_site_labels, esr_labels, random_exon_labels, intron_labels, random_intron_labels, other_labels)) +
  theme(axis.title.y = element_text(size = axis_title_y), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none') +
  labs(x = '', y = expression(paste(Delta, ' inclusion index (DHFR)')),
       color = expression(index["WT "]))

gg_no_x_axis <- gg + theme(axis.text.x = element_blank())

ggsave(paste0('../../figs/splicemod/dhfr/splicemod_dhfr_all_categories_no_x_axis', plot_format), 
       gg_no_x_axis, width = 12, height = 3, dpi = 300)
ggsave(paste0('../../figs/splicemod/dhfr/splicemod_dhfr_all_categories', plot_format), 
       gg, width = 12, height = 5, dpi = 300)

# explore no splice vs weak splice
# score differential
data %>%
  filter(category == 'no_spl_a' | category == 'no_spl_d' | category == 'weak_spl_a' | category == 'weak_spl_d') %>%
  ggplot(aes(index_smn1, (new_score - orig_score), color = category)) + geom_point()

# number of changes
data %>% filter(category == 'no_spl_a' | category == 'no_spl_d' | category == 'weak_spl_a' | category == 'weak_spl_d') %>%
  ggplot(aes(num_changes, index_smn1, color = category)) + geom_jitter() + facet_grid (~ category)


###############################################################################
# Exonic motifs
###############################################################################

# calculate change in average HAL score
data <- data %>% 
  rename('avg_HAL_score' = 'avg_exon_effect_score') %>% 
  group_by(ensembl_id) %>% 
  mutate(delta_avg_HAL_score = avg_HAL_score - avg_HAL_score[sub_id == '000']) %>% 
  ungroup()

###########
# HAL
###########
data <- data %>%
  mutate(HAL_bin = case_when(.$delta_avg_HAL_score < 0 ~ 'down',
                             .$delta_avg_HAL_score > 0 ~ 'up',
                             .$delta_avg_HAL_score == 0 ~ 'same'))

# SMN1
gg <- data %>%
  filter(HAL_bin != 'same', seq_type == 'mut') %>%
  ggplot(aes(HAL_bin, dpsi_smn1)) + geom_jitter(alpha = jitter_alpha, aes(color = nat_index_smn1)) + 
  geom_violin(alpha = 0) +
  scale_colour_gradientn(limits = c(-0.005, 1), breaks = seq(0, 1, by = 0.25), colors = pal(321)) + 
  labs(x = expression(paste(Delta, ' average exonic hexamer score (HAL)')), 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  ggsignif::geom_signif(comparisons = list(c('down', 'up')),
                        test = 't.test', map_signif_level = T) +
  stat_summary(fun.y=mean, geom="point", size=2, color="black") +
  theme( axis.title.x = element_text(size = axis_title_x), 
         axis.title.y = element_text(size = axis_title_y),
         axis.text.x = element_text(size = 18), 
         legend.position = 'none')

ggsave(paste0('../../figs/splicemod/smn1/splicemod_smn1_hal', plot_format),
       gg, width = 2.5, height = 3, dpi = 300, scale = 1.3)

# DHFR
gg <- data %>%
  filter(HAL_bin != 'same', seq_type == 'mut') %>%
  ggplot(aes(HAL_bin, dpsi_dhfr)) + geom_jitter(alpha = jitter_alpha , aes(color = nat_index_dhfr)) + 
  geom_violin(alpha = 0) +
  scale_colour_gradientn(limits = c(-0.005, 1), breaks = seq(0, 1, by = 0.25), colors = pal(321)) + 
  labs(x = expression(paste(Delta, ' average exonic hexamer score (HAL)')), 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  ggsignif::geom_signif(comparisons = list(c('down', 'up')),
                        test = 't.test', map_signif_level = T) +
  stat_summary(fun.y=mean, geom="point", size=2, color="black") +
  theme(axis.title.x = element_text(size = axis_title_x), 
        axis.title.y = element_text(size = axis_title_y),
        axis.text.x = element_text(size = 18),  
        legend.position = 'none')

ggsave(paste0('../../figs/splicemod/dhfr/splicemod_dhfr_hal', plot_format),
       gg, width = 2.5, height = 3, dpi = 300, scale = 1.3)

# data %>% 
#     group_by(ensembl_id) %>% 
#     mutate(avg_HAL_score_nat = avg_HAL_score[sub_id == '000']) %>% 
#     ungroup() %>% 
#     mutate(avg_HAL_score_fc = (avg_HAL_score - avg_HAL_score_nat) / abs(correct_don_score_nat))


###########
# Ke 2011
###########
data <- data %>% 
  group_by(ensembl_id) %>% 
  mutate(delta_Ke2011_avg_score = Ke2011_avg_score - Ke2011_avg_score[sub_id == '000'],
         delta_Ke2011_ESE_avg_score = Ke2011_ESE_avg_score - Ke2011_ESE_avg_score[sub_id == '000'],
         delta_Ke2011_ESS_avg_score = Ke2011_ESS_avg_score - Ke2011_ESS_avg_score[sub_id == '000']) %>% 
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
# SMN1
gg <- data %>%
  filter(Ke_bin != 'same', seq_type == 'mut') %>%
  ggplot(aes(Ke_bin, dpsi_smn1)) + geom_jitter(alpha = jitter_alpha, aes(color = nat_index_smn1)) + 
  geom_violin(alpha=0) +
  scale_colour_gradientn(limits = c(-0.005, 1), breaks = seq(0, 1, by = 0.25), colors = pal(321)) +
  labs(x = expression(paste(Delta, ' Ke 2011 average exonic hexamer score')), 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  ggsignif::geom_signif(comparisons = list(c('down', 'up')),
                        test = 't.test', map_signif_level = T) +
  stat_summary(fun.y=mean, geom="point", size=2, color="black") +
  theme(axis.title = element_text(size = general_text), axis.text = element_text(size = axis_text), 
        legend.position = 'none')

ggsave(paste0('../../figs/splicemod/smn1/splicemod_smn1_Ke11', plot_format),
       gg, width = 3, height = 3, dpi = 300, scale = 1.3)

# DHFR
gg <- data %>%
  filter(Ke_bin != 'same', seq_type == 'mut') %>%
  ggplot(aes(Ke_bin, dpsi_dhfr)) + geom_jitter(alpha = jitter_alpha, aes(color = nat_index_dhfr)) + 
  geom_violin(alpha=0) +
  scale_colour_gradientn(limits = c(-0.005, 1), breaks = seq(0, 1, by = 0.25), colors = pal(321)) +
  labs(x = expression(paste(Delta, ' Ke 2011 average exonic hexamer score')), 
       y = expression(paste(Delta, ' inclusion index')),
       color = expression(index["WT "])) +
  ggsignif::geom_signif(comparisons = list(c('down', 'up')),
                        test = 't.test', map_signif_level = T) +
  stat_summary(fun.y=mean, geom="point", size=2, color="black") +
  theme(axis.title = element_text(size = general_text), axis.text = element_text(size = axis_text), 
        legend.position = 'none')

ggsave(paste0('../../figs/splicemod/dhfr/splicemod_dhfr_Ke11', plot_format),
       gg, width = 3, height = 3, dpi = 300, scale = 1.3)

# # splicemod exon mutation categories
# esr_labels <- c('weaken\nESEs', 'destroy\nstrongest\nESE', 'weaken\nESSs', 'destroy\nstrongest\nESS')
# data %>% 
#     filter(category %in% esr_categories) %>% 
#     ggplot(aes(x = category, y = dpsi_smn1)) + geom_jitter(alpha = 0.75, aes(color = nat_index_dhfr)) + 
#     geom_boxplot(alpha = 0) +
#     scale_colour_gradientn(limits = c(-0.005, 1), breaks = seq(0, 1, by = 0.25), colors = pal(321)) + 
#     scale_x_discrete(limits = esr_categories, labels = esr_labels) +
#     labs(x = '', y = expression(paste(Delta, ' exon inclusion index (SMN1)')),
#          color = 'exon\ninclusion\nindex (WT)')
