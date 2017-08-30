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
    

# read in master data frames to get splicing index
exac_spanr <- exac_spanr %>%
    left_join(read.table('../../processed_data/exac/exac_data_clean.txt', 
                         sep='\t', header=T) %>% 
                  select(id, v2_dpsi, category, nat_v2_index, v2_index, strong_lof), 
              by='id') %>% 
    filter(!is.na(v2_dpsi))

exac_spanr <- exac_spanr %>%
    mutate(dpsi_max_tissue = dpsi_max_tissue/100, 
           dpsi_spanr_capped = ifelse(abs(dpsi_max_tissue) >= nat_v2_index,
                                      sign(dpsi_max_tissue)*nat_v2_index, 
                                      dpsi_max_tissue))


fit <- signif(summary(lm(dpsi_spanr_capped ~ v2_dpsi, exac_spanr))$r.squared, 3)
pearson <- signif(cor(exac_spanr$dpsi_spanr_capped, exac_spanr$v2_dpsi), 3)
spearman <- signif(cor(exac_spanr$dpsi_spanr_capped, exac_spanr$v2_dpsi, method = 'spearman'), 3)
ggplot(exac_spanr, aes(v2_dpsi, dpsi_spanr_capped)) + geom_point(alpha = 0.25) +
    geom_smooth(method = 'lm') + 
    labs(x = expression(paste(Delta, ' inclusion index')), 
         y = expression(paste('SPANR predicted ', Delta, 'PSI')),
         title = 'ExAC SNVs') +
    # annotate('text', parse=T, label=paste('R^2==', fit), x = 0.75, y = -0.70)
    annotate('text', label=paste0('R^2 = ', fit, 
                                  '\nr = ', pearson,
                                  '\nspearman = ', spearman), x = 0.75, y = -0.50)


exac_spanr <- exac_spanr %>%
    mutate(same_dir_change = ifelse(sign(dpsi_max_tissue) == sign(v2_dpsi), TRUE, FALSE))

num_success <- exac_spanr %>% filter(category == 'mutant', same_dir_change == TRUE) %>% nrow()
num_trials <- exac_spanr %>% filter(category == 'mutant') %>% nrow()

binom.test(num_success, num_trials)

dpsi_threshold <- -0.50
exac_spanr <- exac_spanr %>% 
    mutate(spanr_strong_lof = ifelse(dpsi_spanr_capped <= dpsi_threshold, T, F))
    

# true positive, both spanr and our calls agree
num_true_pos <- filter(exac_spanr, spanr_strong_lof == T, strong_lof == T) %>% nrow()
# all positive calls from spanr
num_spanr_pos <- filter(exac_spanr, spanr_strong_lof == T) %>% nrow()
# false positives, called as positive by spanr but not strong LoF in our calls
num_false_pos <- filter(exac_spanr, spanr_strong_lof == T, strong_lof == F) %>% nrow()
# true negative, both spanr and our calls agree
num_true_neg <- filter(exac_spanr, spanr_strong_lof == F, strong_lof == F) %>% nrow()
# false negative, called as negative by spanr but called strong LoF in our assay
num_false_neg <- filter(exac_spanr, spanr_strong_lof == F, strong_lof == T) %>% nrow()
# all negative spanr calls
num_spanr_neg <- filter(exac_spanr, spanr_strong_lof == F) %>% nrow()

# precision, how many selected items are relevant
precision <- (num_true_pos / (num_true_pos + num_false_pos)) * 100
# recall/sensitivity, how many relevant items are selected
recall <- (num_true_pos / (num_true_pos + num_false_neg)) * 100
# specificity, ability to correctly detect negatives 
specificity <- (num_true_neg / (num_true_neg + num_false_pos)) * 100

print(paste0('precision:', precision))
print(paste0('recall/sensitivity:', recall))
print(paste0('specificity:', specificity))

                 