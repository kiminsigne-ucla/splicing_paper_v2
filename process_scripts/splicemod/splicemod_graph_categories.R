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

# custom color palette
color.palette <- function(steps, n.steps.between=NULL, ...){
    
    if (is.null(n.steps.between)) n.steps.between <- rep(0, (length(steps) - 1))
    if (length(n.steps.between) != length(steps) - 1) stop("Must have one less n.steps.between value than steps")
    
    fill.steps <- cumsum(rep(1, length(steps)) + c(0,n.steps.between))
    RGB <- matrix(NA, nrow = 3, ncol = fill.steps[length(fill.steps)])
    RGB[,fill.steps] <- col2rgb(steps)
    
    for (i in which(n.steps.between > 0)) {
        col.start = RGB[,fill.steps[i]]
        col.end = RGB[,fill.steps[i + 1]]
        for (j in seq(3)) {
            vals <- seq(col.start[j], col.end[j], length.out = n.steps.between[i] + 2)[2:(2 + n.steps.between[i] - 1)]  
            RGB[j,(fill.steps[i] + 1):(fill.steps[i + 1] - 1)] <- vals
        }
    }
    
    new.steps <- rgb(RGB[1,], RGB[2,], RGB[3,], maxColorValue = 255)
    pal <- colorRampPalette(new.steps, ...)
    return(pal)
}

# Specify new color palette
steps <- c("blue2", "cyan", "white", "yellow", "red2")
pal <- color.palette(steps, c(160,1,1,160), space = "rgb")



data <- read.table('../../processed_data/splicemod/splicemod_data_clean.txt', sep = '\t', header = T, 
                   colClasses = c('sub_id' = 'character')) %>% 
    filter(rep_quality == 'high')

# read in re-scored reference file
updated_ref <- read.csv('../../ref/splicemod/splicemod_ref_rescored.txt', header = T, sep = '\t',
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
           don_score_fold_change = (correct_don_score - correct_don_score_nat) / abs(correct_don_score_nat),
           acc_score_fold_change = (correct_acc_score - correct_acc_score_nat) / abs(correct_acc_score_nat)) %>%
    ungroup()

# calculate change in average HAL score
data <- data %>% 
    rename('avg_HAL_score' = 'avg_exon_effect_score') %>% 
    group_by(ensembl_id) %>% 
    mutate(delta_avg_HAL_score = avg_HAL_score - avg_HAL_score[sub_id == '000']) %>% 
    ungroup()

###############################################################################
# MaxEnt splice donor score fold-change
###############################################################################
data %>%
    filter(don_score_fold_change != 0) %>%
    mutate(don_fold_change_bin = cut(don_score_fold_change, c(-12, -2, -1, 0, 1))) %>%
    filter(!is.na(don_fold_change_bin)) %>%
    ggplot(aes(don_fold_change_bin, dpsi_smn1)) + 
    geom_jitter(alpha = 0.50, aes(color = nat_index_smn1)) + 
    scale_x_discrete(labels = c('<= -2', '(-2, -1]', '(-1, -0)', '(0, 1]')) +
    geom_boxplot(alpha = 0) + 
    scale_colour_gradientn(limits = c(-0.005, 1), breaks = seq(0, 1, by = 0.25), colors = pal(321)) +
    labs(x = 'fold change MaxEnt splice donor score', 
         y = expression(paste(Delta, ' inclusion index (SMN1)')),
         color = 'inclusion\nindex (WT)') +
    theme(axis.text = element_text(size = 15), text = element_text(size = 18))

ggsave('../../figs/splicemod/splicemod_donor_fc.tiff', width = 6, height = 4, dpi = 100)

###############################################################################
# MaxEnt splice acceptor score fold-change
###############################################################################
data %>%
    filter(acc_score_fold_change != 0) %>%
    mutate(acc_fold_change_bin = cut(acc_score_fold_change, c(-184, -2, -1, 0, 1))) %>%
    filter(!is.na(acc_fold_change_bin)) %>%
    ggplot(aes(acc_fold_change_bin, dpsi_smn1)) + geom_jitter(alpha = 0.50, aes(color = nat_index_smn1)) + 
    scale_x_discrete(labels = c('<= -2', '(-2, -1]', '(-1, -0)', '(0, 1]')) +
    geom_boxplot(alpha = 0) + 
    scale_colour_gradientn(limits = c(-0.005, 1), breaks = seq(0, 1, by = 0.25), colors = pal(321)) +
    labs(x = 'fold change MaxEnt splice acceptor score', y = expression(paste(Delta, ' inclusion index (SMN1)')),
         color = 'inclusion\nindex (WT)') +
    theme(axis.text = element_text(size = 15), text = element_text(size = 18))

ggsave('../../figs/splicemod/splicemod_acceptor_fc.tiff', width = 6, height = 4, dpi = 100)

###############################################################################
# All splicemod categories
###############################################################################
esr_categories <- c('rmv_Ke2011_ESE', 'clst_Ke2011_ESE', 'rmv_Ke2011_ESS', 'clst_Ke2011_ESS')
esr_labels <- c('weaken ESEs', 'destroy strongest ESE', 'weaken ESSs', 'destroy strongest ESS')

splice_site_categories <- c( 'weak_spl_a', 'weak_spl_d', 'p_weak_spl',
                             'same_splice_a', 'same_splice_d',
                             'rmv_me_splice_acceptor', 'rmv_me_splice_donor', 
                             'csplice_a', 'csplice_d',
                             'no_spl_a', 'no_spl_d')
splice_site_labels <- c('weaken acceptor', 'weaken donor', 'weaken donor + acceptor',
                        'same score acceptor', 'same score donor',
                        'weaken spurious acceptor', 'weaken spurious donor',
                        'destroy spurious acceptor', 'destroy spurious donor',
                        'destroy acceptor', 'destroy donor')

intron_categories <- c('clst_Vlkr07_AICS', 'clst_Vlkr07_DICS', 'rmv_Vlkr07_AICS', 'rmv_Vlkr07_DICS')
intron_labels <- c('weaken AICS', 'weaken DICS', 'destroy AICS', 'destroy DICS')

random_categories <- c('rnd_intron_1nt', 'rnd_intron_2nt', 'rnd_intron_3nt', 'rnd_intron_5nt', 
                       'aggr_intron', 'p_aggr_intr', 
                       'rnd_exon_1nt', 'rnd_exon_2nt', 'rnd_exon_3nt', 'rnd_exon_5nt', 'aggr_exon',
                       'aggr_both')
random_labels <- c('random 1nt intron', 'random 2nt intron', 'random 3nt intron', 'random 5nt intron',
                   'aggressive intron (all syn. mut.)', 'aggr. + random intron',
                   'random 1nt exon', 'random 2nt exon', 'random 3nt exon', 'random 5nt exon',
                   'aggressive exon (all syn. mut.)', 'aggr. intron + exon')

other_categories <- c('cnsrv_1nt', 'cnsrv_3nt', 'RBPmats', 'variation')
other_labels <- c('conserved 1nt', 'conservered 3nt', 'destroy RBP motifs', 'dbSNPs')

data %>% 
    filter(seq_type == 'mut') %>% 
    mutate(category_fctr = factor(category,
                                  levels = c(splice_site_categories, esr_categories,
                                             intron_categories, random_categories, other_categories))) %>% 
    ggplot(aes(x = category_fctr, y = dpsi_smn1)) + geom_jitter(alpha = 0.50, aes(color = nat_index_smn1)) + geom_boxplot(alpha=0) +
    scale_colour_gradientn(limits = c(-0.005, 1), breaks = seq(0, 1, by = 0.25), colors = pal(321)) + 
    scale_x_discrete(labels = c(splice_site_labels, esr_labels, intron_labels, random_labels, other_labels)) +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = expression(paste(Delta, ' inclusion index (SMN1)')),
         color = 'exon\ninclusion\nindex (WT)')

ggsave('../../figs/splicemod/splicemod_all_categories.tiff', width = 12, height = 5, dpi = 300)

###############################################################################
# Exonic changes
###############################################################################

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

# HAL
data <- data %>%
    mutate(HAL_bin = case_when(.$delta_avg_HAL_score < 0 ~ 'down',
                               .$delta_avg_HAL_score > 0 ~ 'up',
                               .$delta_avg_HAL_score == 0 ~ 'same'))

data %>%
    filter(HAL_bin != 'same', seq_type == 'mut') %>%
    ggplot(aes(HAL_bin, dpsi_smn1)) + geom_jitter(alpha = 0.25, aes(color = nat_index_smn1)) + 
    geom_violin(alpha=0) +
    scale_colour_gradientn(limits = c(-0.005, 1), breaks = seq(0, 1, by = 0.25), colors = pal(321)) + 
    labs(x = expression(paste(Delta, ' HAL mean exonic hexamer score')), 
         y = expression(paste(Delta, ' exon inclusion index (SMN1)')),
         color = 'exon\ninclusion\nindex (WT)') +
    ggsignif::geom_signif(comparisons = list(c('down', 'up')),
                          test = 't.test', map_signif_level = T) +
    stat_summary(fun.y=mean, geom="point", size=2, color="black") +
    theme(axis.title.y = element_text(size = 19), axis.text = element_text(size = 16), 
          axis.title.x = element_text(size = 16))

data %>% 
    group_by(ensembl_id) %>% 
    mutate(avg_HAL_score_nat = avg_HAL_score[sub_id == '000']) %>% 
    ungroup() %>% 
    mutate(avg_HAL_score_fc = (avg_HAL_score - avg_HAL_score_nat) / abs(correct_don_score_nat))

# Ke 2011
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

data %>%
    filter(Ke_bin != 'same', seq_type == 'mut') %>%
    ggplot(aes(Ke_bin, dpsi_smn1)) + geom_jitter(alpha = 0.25, aes(color = nat_index_smn1)) + geom_violin(alpha=0) +
    scale_colour_gradientn(limits = c(-0.005, 1), breaks = seq(0, 1, by = 0.25), colors = pal(321)) +
    labs(x = expression(paste(Delta, ' Ke 2011 mean exonic hexamer score')), 
         y = expression(paste(Delta, ' exon inclusion index (SMN1)')),
         color = 'exon\ninclusion\nindex (WT)') +
    ggsignif::geom_signif(comparisons = list(c('down', 'up')),
                          test = 't.test', map_signif_level = T) +
    stat_summary(fun.y=mean, geom="point", size=2, color="black") +
    theme(axis.title.y = element_text(size = 19), axis.text = element_text(size = 16), 
          axis.title.x = element_text(size = 16))
