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

exac_spanr <- read.table('../../processed_data/exac/exac_SPANR_scores.txt', 
                         sep='\t', header=T)
exac_exon_vars <- read.table('../../processed_data/exac/exac_HAL_scores.txt', 
                             sep='\t', header=T)
exac_intron_cons <- read.table('../../processed_data/exac/exac_data_intron_cons.txt',
                               sep = '\t', header = T)
exac_ref_rescored <- read.table('../../ref/exac/exac_ref_rescored.txt',
                                sep = '\t', header = T)

data_all <- exac_ref_rescored %>% 
    select(id, avg_exon_effect_score:correct_don_score, -correct_acc_seq, -correct_don_seq) %>% 
    full_join(select(exac_exon_vars, id, nat_v2_index, v2_dpsi, hal_dpsi = DPSI_pred), by = 'id') %>% 
    full_join(select(exac_spanr, id, spanr_dpsi = dpsi_spanr_capped), by = 'id') %>% 
    full_join(select(exac_intron_cons, id, upstr_intron_mean_cons:junc_avg_don), by = 'id')

corr <- cor(select(data_all, -id), use ='p')
model <- lm(v2_dpsi ~ correct_acc_score + correct_don_score, data_all)