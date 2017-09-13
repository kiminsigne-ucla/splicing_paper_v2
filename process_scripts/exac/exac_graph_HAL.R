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


exac_exon_vars <- read.table('../../processed_data/exac/exac_HAL_scores.txt', 
                             sep='\t', header=T) %>% 
    filter(!is.na(v2_dpsi))


fit <- signif(summary(lm(DPSI_pred ~ v2_dpsi, exac_exon_vars))$r.squared, 3)
pearson <- signif(cor(exac_exon_vars$DPSI_pred, exac_exon_vars$v2_dpsi), 3)
spearman <- signif(cor(exac_exon_vars$DPSI_pred, exac_exon_vars$v2_dpsi, method = 'spearman'), 3)
ggplot(exac_exon_vars, aes(v2_dpsi, DPSI_pred)) + geom_point(alpha = 0.25) +
    geom_smooth(method = 'lm') +
    labs(x = expression(paste(Delta, ' inclusion index ')),
         y = expression(paste('HAL predicted ', Delta, 'PSI')),
         title = 'ExAC exonic variants') +
    annotate('text', label=paste0('R^2 = ', fit, 
                                  '\nr = ', pearson,
                                  '\nspearman = ', spearman), x = 0.75, y = -0.50)


exac_exon_vars <- exac_exon_vars %>%
    mutate(same_dir_change = ifelse(sign(DPSI_pred) == sign(v2_dpsi), TRUE, FALSE))

num_success <- exac_exon_vars %>% filter(category == 'mutant', same_dir_change == TRUE) %>% nrow()
num_trials <- exac_exon_vars %>% filter(category == 'mutant') %>% nrow()

binom.test(num_success, num_trials)

dpsi_threshold <- -0.50

exac_exon_vars <- exac_exon_vars %>% 
    mutate(hal_strong_lof = ifelse(DPSI_pred <= dpsi_threshold, 'True', 'False'))



# true positive, both HAL and our calls agree
num_true_pos <- filter(exac_exon_vars, hal_strong_lof == 'True', strong_lof == 'True') %>% nrow()
# all positive calls from HAL
num_HAL_pos <- filter(exac_exon_vars, hal_strong_lof == 'True') %>% nrow()
# false positives, called as positive by HAL but not strong LoF in our calls
num_false_pos <- filter(exac_exon_vars, hal_strong_lof == 'True', strong_lof == 'False') %>% nrow()
# true negative, both HAL and our calls agree
num_true_neg <- filter(exac_exon_vars, hal_strong_lof == 'False', strong_lof == 'False') %>% nrow()
# false negative, called as negative by HAL but called strong LoF in our assay
num_false_neg <- filter(exac_exon_vars, hal_strong_lof == 'False', strong_lof == 'True') %>% nrow()
# all negative HAL calls
num_HAL_neg <- filter(exac_exon_vars, hal_strong_lof == 'False') %>% nrow()

# precision, how many selected items are relevant
precision <- (num_true_pos / (num_true_pos + num_false_pos)) * 100
# recall/sensitivity, how many relevant items are selected
recall <- (num_true_pos / (num_true_pos + num_false_neg)) * 100
# specificity, ability to correctly detect negatives 
specificity <- (num_true_neg / (num_true_neg + num_false_pos)) * 100

print(paste0('precision:', precision))
print(paste0('recall/sensitivity:', recall))
print(paste0('specificity:', specificity))

