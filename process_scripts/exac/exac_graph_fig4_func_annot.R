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
# custom color palette
source("../color_palette.R")
# Specify new color palette
steps <- c("blue2", "cyan", "white", "yellow", "red2")
pal <- color.palette(steps, c(160,1,1,160), space = "rgb")

data <- read.table('../../processed_data/exac/exac_func_annot.txt',
                   sep = '\t', header = T)

lof <- data %>% 
    filter(strong_lof == T, category == 'mutant')

# Figure 4A, conservation intron vs. exon
data_cons <- data %>% 
    mutate(cons_bin = cut(mean_cons_score, breaks = c(0, 0.5, 1.0), 
                          include.lowest = T,  labels = c('low', 'high'))) %>%
    filter(!is.na(cons_bin))

lof_cons <- lof %>% 
    mutate(cons_bin = cut(mean_cons_score, breaks = c(0, 0.5, 1.0), 
                          include.lowest = T, labels = c('low', 'high'))) %>%
    filter(!is.na(cons_bin))

data_cons$label_renamed <- factor(data_cons$label, 
                                  levels=c("upstr_intron", "exon", "downstr_intron"), 
                                  labels=c("Intron\nUpstr.", "Exon", "Intron\nDownstr."))

lof_cons$label_renamed <- factor(lof_cons$label, 
                                 levels=c("upstr_intron", "exon", "downstr_intron"), 
                                 labels=c("Intron\nUpstr.", "Exon", "Intron\nDownstr.") )

data_cons_count <- data_cons %>%
    group_by(label_renamed, cons_bin) %>%
    summarise(`total SNP count` = n())

lof_cons_count <- lof_cons %>%
    group_by(label_renamed, cons_bin) %>%
    summarise(`SNP count (loss-of-splicing)` = n())

cons_count <- full_join(data_cons_count, lof_cons_count, 
                        by = c("label_renamed", "cons_bin")) %>% 
    mutate(propFreq = `SNP count (loss-of-splicing)` / `total SNP count` * 100,
           label_cons = paste(label_renamed, "\n", cons_bin, "\nconserv.", sep = " "))

cons_count$cons_bin <- factor(cons_count$cons_bin, levels=c("low","high"), 
                              labels=c("Low Conservation", "High Conservation"))

fig4a <- cons_count %>% 
    ggplot(aes(label_renamed, propFreq)) +
    geom_histogram(stat = 'identity', width = 0.5) +
    ylab("% loss-of-splicing SNVs") +
    xlab("") +
    facet_wrap(~ cons_bin) +
    geom_hline(yintercept = 1050/29531*100, linetype = "dashed", color = "grey10") + 
    scale_y_continuous(expand = c(0,0)) +
    expand_limits(y = 25) +
    theme(strip.text = element_text(size = 16),
          strip.background = element_rect(fill="lightgrey"),
          axis.title.y = element_text(size = 20),
          axis.text = element_text(size = 16))

ggsave(paste0("../../figs/exac/exac_fig4A_phastCons_comparison_prop", plot_format), 
       width = 6, height = 4.5, units = 'in')

###############################################################################
# absolute counts for LoF SNPs, by conservation
###############################################################################
fig4b <- cons_count %>% 
    ggplot(aes(label_renamed, `SNP count (loss-of-splicing)`)) +
    geom_histogram(stat = 'identity', width = 0.5) +
    ylab("Number loss-of-splicing SNVs") +
    xlab("") +
    facet_wrap(~ cons_bin) +
    scale_y_continuous(expand = c(0,0)) +
    expand_limits(y = 425) +
    theme(strip.text = element_text(size = 16),
          strip.background = element_rect(fill="lightgrey"),
          axis.title.y = element_text(size = 18),
          axis.text = element_text(size = 16))


ggsave(paste0("../../figs/exac/exac_fig4B_phastCons_abs_num", plot_format), 
       width = 6, height = 4.5, units = 'in')

###############################################################################
# Figure 4C, Allele Frequency 
###############################################################################

data <- data %>% 
    mutate(AF_bin = cut(AF, 
                        breaks = c(0, 0.00001, 0.000025, 0.0001, 0.001, 1), 
                        include.lowest = T,
                        labels = c('Singleton', '2 or 3\nalleles', 
                                   '4-10\nalleles', '0.01%', '>0.1%' ))) 

fig4c <- data %>% 
    filter(!is.na(AF_bin)) %>% 
    ggplot(aes(AF_bin, v2_dpsi)) + 
    geom_jitter(alpha = 0.50, aes(color = nat_v2_index), width = 0.35) + 
    scale_colour_gradientn(limits = c(-0.005,1), 
                           breaks=seq(0,1, by = 0.25), 
                           colors=pal(321), 
                           expression(index["WT "])) +
    geom_boxplot(alpha = 0) +
    labs(x = " ", y = expression(paste(Delta, ' inclusion index'))) +
    theme(legend.position = 'none', 
          axis.title.y = element_text(size = 18),
          axis.text = element_text(size = 16)) 

ggsave(paste0("../../figs/exac/exac_fig4c_allele_frequency_binned", plot_format), 
       width = 5.5, height = 4, units = 'in')

###############################################################################
# Figure 4D, probability of gene being loss-of-function intolerant 
###############################################################################
dpsi_threshold <- -0.50

num_intolerant_lof <- data %>% 
    filter(category == 'mutant', pLI_high == TRUE, v2_dpsi <= dpsi_threshold ) %>% nrow()
num_intolerant_not_lof <- data %>% 
    filter(category == 'mutant', pLI_high == TRUE, v2_dpsi > dpsi_threshold ) %>% nrow()
num_tolerant_lof <- data %>%
    filter(category == 'mutant', pLI_high == FALSE, v2_dpsi <= dpsi_threshold ) %>% nrow()
num_tolerant_not_lof <- data %>% 
    filter(category == 'mutant', pLI_high == FALSE, v2_dpsi > dpsi_threshold ) %>% nrow()

df <- matrix(c(num_intolerant_lof, num_tolerant_lof,
                num_intolerant_not_lof, num_tolerant_not_lof), 
              nrow = 2, 
              dimnames = list(c('pLI >= 0.90', 'pLI < 0.90'),
                              c('dPSI <= -0.50', 'dPSI > -0.50')))

intolerant_ratio =  num_intolerant_lof / (num_intolerant_lof + num_intolerant_not_lof)
tolerant_ratio = num_tolerant_lof / (num_tolerant_not_lof + num_tolerant_lof)
ratio.df <-data.frame( fraction_of_strong_LoF_genes = c(intolerant_ratio, tolerant_ratio), tolerance = c("intolerant", "tolerant"))
# colnames(ratio.df) <- c("ratio")

ratio.df %>%
    ggplot(aes(tolerance, fraction_of_strong_LoF_genes)) + 
    geom_col(width = 0.5, position = position_dodge(width = 1.8)) + 
    ylab("Fraction loss-of-splicing SNVs") + xlab("pLI") + ylim(0,0.05) +
    geom_hline(yintercept = 1050/29531, linetype = "dashed", color = "grey10") +
    scale_y_continuous(expand = c(0, 0)) + 
    coord_equal(1/0.015)

fisher.test(df, alternative = 'less')

ggsave(paste0("../../figs/exac/exac_fig3D_pLI_enrichment", plot_format), 
       width = 2.8, height = 4, units = 'in')

