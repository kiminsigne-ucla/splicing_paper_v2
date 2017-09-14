load_pkgs <- function(pkgs){
    new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, 'Package'])]
    if(length(new_pkgs)) install.packages(new_pkgs)
    for(pkg in pkgs){
        suppressWarnings(suppressMessages(library(pkg, character.only = T)))
    }
}

pkgs <- c('dplyr', 'tidyr', 'ggplot2', 'cowplot', 'grid', 'gtable', 'ggsignif')
load_pkgs(pkgs)

options(stringsAsFactors = F, warn = -1, warnings = -1)

plot_format <- '.tiff'
hi_res <- 600

# custom color palette
source("../color_palette.R")
# Specify new color palette
steps <- c("blue2", "cyan", "white", "yellow", "red2")
pal <- color.palette(steps, c(160, 1, 1, 160), space = "rgb")

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
                                  labels=c("Intron\nupstr.", "Exon", "Intron\ndownstr."))

lof_cons$label_renamed <- factor(lof_cons$label, 
                                 levels=c("upstr_intron", "exon", "downstr_intron"), 
                                 labels=c("Intron\nupstr.", "Exon", "Intron\ndownstr.") )

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
                              labels=c("Low conservation", "High conservation"))

cons_count %>% 
    ggplot(aes(label_renamed, propFreq, fill = factor(cons_bin))) +
    geom_histogram(stat = 'identity', width = 0.75, position = "dodge") +
    ylab("Percentage of\nloss-of-splicing variants") +
    xlab("") +
    # facet_wrap(~ cons_bin) +
    geom_hline(yintercept = 3.6, linetype = "dashed", color = "grey20") + 
    scale_y_continuous(breaks = c(0, 3.6, 10, 15, 20, 25), expand = c(0,0)) +
    expand_limits(y = 26.5) +
    theme_bw() + 
    theme(strip.text = element_text(size = 18.5),
          strip.background = element_rect(fill = "#E0E0E0", color = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(1, "lines"),
          panel.border = element_rect(fill = NA, color = "grey50"),
          axis.title.y = element_text(size = 19, vjust = 2.75),
          axis.text.y = element_text(size = 14, color = "grey20"),
          axis.text.x = element_text(size = 18, color = "black"),
          axis.ticks.x = element_blank(),
          legend.title = element_blank(),
          legend.position = c(0.4,0.9),
          legend.text = element_text(size = 16),
          legend.key.height = unit(0.25, "inch"),
          legend.key.width = unit(0.4, "inch")
          ) +
    scale_fill_manual(values = c('#458B00','#2E9FFE'))

ggsave(paste0("../../figs/exac/exac_fig4A_phastCons_comparison_prop", plot_format),
        width = 4.5, height = 5, units = 'in', dpi = hi_res)

###############################################################################
# absolute counts for LoF SNPs, by conservation
###############################################################################

cons_count %>% 
    ggplot(aes(label_renamed, `SNP count (loss-of-splicing)`, fill = factor(cons_bin))) +
    geom_histogram(stat = 'identity', width = 0.75, position = "dodge") +
    ylab("Number of\nloss-of-splicing variants") +
    xlab("") +
    scale_y_continuous(expand = c(0,0)) +
    expand_limits(y = 450) +
    theme_bw() + 
    theme(strip.text = element_text(size = 19),
          strip.background = element_rect(fill = "#E0E0E0", color = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, colour = "grey50"),
          axis.title.y = element_text(size = 19, vjust = 2.75),
          axis.text.y = element_text(size = 14, color = "grey20"),
          axis.text.x = element_text(size = 18, color = "black"),
          axis.ticks.x = element_blank(),
          legend.position = c(0.4,0.9),
          legend.title = element_blank(),
          legend.text = element_text(size = 16),
          legend.key.height = unit(0.25, "inch"),
          legend.key.width = unit(0.4, "inch")
    ) +
  scale_fill_manual(values = c('#458B00',
                               '#2E9FFE'))
                               
ggsave(paste0("../../figs/exac/exac_fig4B_phastCons_abs_num", plot_format),
       width = 4.5, height = 5, units = 'in', dpi = hi_res)


##############################
# Figure 4C: pLI
##############################
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

intolerant_percent =  num_intolerant_lof / (num_intolerant_lof + num_intolerant_not_lof) * 100
tolerant_percent = num_tolerant_lof / (num_tolerant_not_lof + num_tolerant_lof) * 100
percent.df <-data.frame( fraction_of_strong_LoF_genes = c(intolerant_percent, tolerant_percent), tolerance = c('intolerant', 'tolerant'))

percent.df %>%
  ggplot(aes(tolerance, fraction_of_strong_LoF_genes)) + 
  geom_col(width = 0.425) + 
  ylab("% of loss-of-splicing variants") + xlab("pLI") + ylim(0,5) +
  geom_hline(yintercept = 1050/29531*100, linetype = "dashed", color = "grey40") +
  scale_y_continuous(expand = c(0, 0)) + 
  expand_limits(y = 4.8) +
  # coord_equal(1/1.125) +
  theme_bw() +
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color = "grey50"),
        axis.title.y = element_text(size = 14, color = "black", vjust = 2),
        axis.title.x = element_text(size = 20, color = "black", vjust = -0.5),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12, color = "grey20"),
        axis.text.x = element_text(size = 13, color = "black")) 

fisher.test(df, alternative = 'less')

ggsave(paste0("../../figs/exac/exac_fig4C_pLI_enrichment", plot_format), 
       width = 2.25, height = 3.5, units = 'in', dpi = hi_res)

###############################################################################
# Figure 4D, ExAC global allele frequency 
###############################################################################

data <- data %>% 
    mutate(AF_bin = cut(AF, 
                        breaks = c(0, 0.00001, 0.000025, 0.0001, 0.001, 1), 
                        include.lowest = T,
                        labels = c('Singleton', 'AC = 2-3', 
                                   'AC = 4-10', '0.01%', '>0.1%' ))) 

fig4d <- data %>% 
    filter(!is.na(AF_bin)) %>% 
    ggplot(aes(AF_bin, v2_dpsi)) + 
    geom_jitter(alpha = 0.25, aes(color = nat_v2_index), width = 0.35) + 
    scale_colour_gradientn(limits = c(-0.005,1), 
                           breaks=seq(0,1, by = 0.25), 
                           colors=pal(321), 
                           expression(index["WT "])) +
    geom_violin(alpha = 0, color = "grey10", size = 0.5) +
    stat_summary(fun.y = median, geom = "point", size = 1, color = "grey10") +
    labs(x = "ExAC global allele frequency", y = expression(paste(Delta, ' inclusion index'))) +
    theme_bw() + 
    theme(legend.position = 'none', 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, color = "grey50"),
          axis.title.y = element_text(size = 17, vjust = 12),
          axis.title.x = element_text(size = 17, vjust = -0.5),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 10, color = "grey20"),
          axis.text.x = element_text(size = 14, color = "grey10", angle = 45, vjust = 0.55)) 

ggsave(paste0("../../figs/exac/exac_fig4D_allele_frequency_binned", plot_format), 
       width = 4, height = 3.5, units = 'in', dpi = hi_res)





