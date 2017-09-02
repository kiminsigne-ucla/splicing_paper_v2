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

pkgs <- c('dplyr', 'tidyr', 'ggplot2', 'cowplot', 'ggExtra', 'grid', 'Unicode')
load_pkgs(pkgs)

options(stringsAsFactors = F, warn = -1, warnings = -1)
plot_format <- '.tiff'

data <- read.table('../../processed_data/splicemod/splicemod_data_clean.txt', sep = '\t', header = T, 
                   colClasses = c('sub_id' = 'character'))

###############################################################################

# Index across SMN1 replicates
corr <- signif(cor(data$index_R1_smn1, data$index_R2_smn1, use = 'p'), 3)
gg <- data %>% 
    ggplot(aes(index_R1_smn1, index_R2_smn1)) + geom_point(alpha = 0.10, aes(color = replicability_smn1)) +
    geom_smooth(method = 'lm', color = 'blue') +
    scale_color_manual(values = c('black', '#0033CC')) +
    labs(x = 'inclusion index (Rep. 1)', y = 'inclusion index (Rep. 2)') +
    theme(legend.position = 'None', text = element_text(size = 20), axis.text = element_text(size = 16)) +
    annotate('text', x = 0.05, y = 0.895, size = 7, label = paste("italic(r)"), parse=TRUE) +
    annotate('text', x = 0.22, y = 0.90, size = 7, label = paste("=", corr))
ggsave(paste0('../../figs/splicemod/smn1/splicemod_smn1_replicates', plot_format), gg, width = 5, height = 5, dpi = 300)

# Index across DHFR replicates
corr <- signif(cor(data$index_R1_dhfr, data$index_R2_dhfr, use = 'p'), 3)
gg <- data %>% 
    ggplot(aes(index_R1_dhfr, index_R2_dhfr)) + geom_point(alpha = 0.10, aes(color = replicability_dhfr)) +
    geom_smooth(method = 'lm', color = 'blue') +
    scale_color_manual(values = c('black', '#0033CC')) +
    labs(x = 'inclusion index (Rep. 1)', y = 'inclusion index (Rep. 2)') +
    theme(legend.position = 'None', text = element_text(size = 20), axis.text = element_text(size = 16)) +
    annotate('text', x = 0.05, y = 0.895, size = 7, label = paste("italic(r)"), parse=TRUE) +
    annotate('text', x = 0.22, y = 0.90, size = 7, label = paste("=", corr))
ggsave(paste0('../../figs/splicemod/dhfr/splicemod_dhfr_replicates', plot_format), gg, width = 5, height = 5, dpi = 300)

# Index across SMN1 and DHFR
corr <- signif(cor(data$index_smn1, data$index_dhfr, use = 'p'), 3)
gg <- ggplot(data %>% filter(replicability_dhfr == 'high', replicability_smn1 == 'high'), aes(index_smn1, index_dhfr)) + geom_point(alpha = 0.10) +
    # scale_color_manual(values = c('black', 'grey')) +
    geom_smooth(method = 'lm', color = 'blue') + 
    labs(x = 'inclusion index (SMN1)', y = 'inclusion index (DHFR)') +
    theme(text = element_text(size = 20), axis.text = element_text(size = 16),
          legend.position = 'none') +
    annotate('text', x = 0.05, y = 0.895, size = 7, label = paste("italic(r)"), parse=TRUE) +
    annotate('text', x = 0.22, y = 0.90, size = 7, label = paste("=", corr))
ggsave(paste0('../../figs/splicemod/splicemod_smn1_dhfr_replicates', plot_format), gg, width = 5, height = 5, dpi = 300)

###############################################################################
# Read distribution across bins
###############################################################################
bin_labels <- c('GFP+', 'GFPint', 'GFP-')   
# SMN1
read_dist <- data %>% 
    filter(rep_quality == 'high') %>%
    rowwise() %>% 
    select(id, DP_R1_norm_smn1, INT_R1_norm_smn1, SP_R1_norm_smn1) %>% 
    mutate(DP = DP_R1_norm_smn1 / (DP_R1_norm_smn1 + INT_R1_norm_smn1 + SP_R1_norm_smn1) * 100,
           INT = INT_R1_norm_smn1 / (DP_R1_norm_smn1 + INT_R1_norm_smn1 + SP_R1_norm_smn1) * 100,
           SP = SP_R1_norm_smn1 / (DP_R1_norm_smn1 + INT_R1_norm_smn1 + SP_R1_norm_smn1) * 100) %>% 
    select(-DP_R1_norm_smn1, -INT_R1_norm_smn1, -SP_R1_norm_smn1) %>% 
    gather(key = 'bin', value = 'bin_percent', DP:SP) %>% 
    # add index for graph order
    left_join(select(data, id, index_R1_smn1), by = 'id')

gg1 <- ggplot(read_dist, aes(x = reorder(id, index_R1_smn1), y = bin)) + 
    geom_tile(aes(fill = bin_percent)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          legend.position = 'none',
          plot.margin = unit(c(0,0,0,0), "in")) +
    labs(x = 'Construct #', y = '', fill = '% contigs') +
    scale_y_discrete(labels = bin_labels) +
    viridis::scale_fill_viridis(option = 'plasma')

# filled
gg2 <- data %>%
    filter(abs(index_R1_smn1 - index_R2_smn1) <= 0.30) %>%
    ggplot(aes(reorder(id, index_R1_smn1), index_R1_smn1)) +
    geom_bar(stat = 'identity', color = 'darkgrey') +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = unit(c(0,0,0, 0.5), "in")) +
    labs(y = 'index')

png(paste0('../../figs/splicemod/smn1/splicemod_smn1_bin_dist.png'),
    width = 13, height = 4, units = 'in', res = 300)
g <- rbind(ggplotGrob(gg2), ggplotGrob(gg1), size = 'last')
id <- g$layout$t[g$layout$name == "panel"]
g$heights[id] <- unit(c(1, 2.5), 'null')
grid.newpage()
grid.draw(g)
dev.off()
# ggsave(paste0('../../figs/splicemod/smn1/splicemod_smn1_bin_dist', plot_format),
#        width = 13, height = 4)

# DHFR
read_dist <- data %>%
    filter(rep_quality == 'high') %>%
    rowwise() %>%
    select(id, DP_R1_norm_dhfr, INT_R1_norm_dhfr, SP_R1_norm_dhfr) %>%
    mutate(DP = DP_R1_norm_dhfr / (DP_R1_norm_dhfr + INT_R1_norm_dhfr + SP_R1_norm_dhfr) * 100,
           INT = INT_R1_norm_dhfr / (DP_R1_norm_dhfr + INT_R1_norm_dhfr + SP_R1_norm_dhfr) * 100,
           SP = SP_R1_norm_dhfr / (DP_R1_norm_dhfr + INT_R1_norm_dhfr + SP_R1_norm_dhfr) * 100) %>%
    select(-DP_R1_norm_dhfr, -INT_R1_norm_dhfr, -SP_R1_norm_dhfr) %>%
    gather(key = 'bin', value = 'bin_percent', DP:SP) %>%
    # add index for graph order
    left_join(select(data, id, index_R1_dhfr), by = 'id')

gg1 <- ggplot(read_dist, aes(x = reorder(id, index_R1_dhfr), y = bin)) +
    geom_tile(aes(fill = bin_percent)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          legend.position = 'none',
          plot.margin = unit(c(0,0,0,0), "in")) +
    labs(x = 'Construct #', y = '', fill = '% contigs') +
    scale_y_discrete(labels = bin_labels) +
    viridis::scale_fill_viridis(option = 'plasma')

gg2 <- data %>% 
    filter(abs(index_R1_dhfr - index_R2_dhfr) <= 0.30) %>%
    ggplot(aes(reorder(id, index_R1_dhfr), index_R1_dhfr)) + 
    geom_bar(stat = 'identity', color = 'darkgrey') + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = unit(c(0,0,0, 0.5), "in")) +
    labs(y = 'index')

png(paste0('../../figs/splicemod/dhfr/splicemod_dhfr_bin_dist.png'),
    width = 13, height = 4, units = 'in', res = 300)
g <- rbind(ggplotGrob(gg2), ggplotGrob(gg1), size = 'last')
id <- g$layout$t[g$layout$name == "panel"]
g$heights[id] <- unit(c(1, 2.5), 'null')
grid.newpage()
grid.draw(g)
dev.off()


###############################################################################
# Histogram of inclusion index (naturals)
###############################################################################

# SMN1

smn1_natCount <- data %>% filter (sub_id == '000', !is.na(index_smn1), !is.na(chr), replicability_smn1 == 'high') %>% nrow()
dhfr_natCount <- data %>% filter (sub_id == '000', !is.na(index_dhfr), !is.na(chr), replicability_dhfr == 'high') %>% nrow()
data %>% filter(sub_id == '000', index_smn1 >= 0.80) %>% nrow()
data %>% filter(sub_id == '000', index_dhfr >= 0.80) %>% nrow()

gg <- data %>% 
  filter (sub_id == '000', !is.na(index_smn1), replicability_smn1 == 'high') %>%
  ggplot(aes(index_smn1)) + 
  geom_histogram(binwidth = 0.02) +
  geom_vline(xintercept = 0.79, color = 'red') +
  xlab('inclusion index') +
  annotate('text', x = 0.115, y = 558, size = 9, label = paste("italic(n)"), parse=TRUE) +
  annotate('text', x = 0.27, y = 560, size = 9, label = paste("=", smn1_natCount)) +
  scale_y_continuous(expand = c(0,0)) +
  expand_limits(y = 650) +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 16))
ggsave(paste0('../../figs/splicemod/smn1/splicemod_smn1_index_hist', plot_format), gg, width = 5, height = 5, dpi = 100)

gg <- data %>%
  filter (sub_id == '000', !is.na(index_dhfr), replicability_dhfr == 'high') %>%
  ggplot(aes(index_dhfr)) + 
  geom_histogram(binwidth = 0.02) +
  # geom_vline(xintercept = 0.79, color = 'red') +
  # geom_rect(aes(xmin = 0.8, xmax = 1.02, ymin = 0, ymax = 675), alpha = 0.1, fill = "lightgray") +
  xlab('inclusion index') +
  annotate('text', x = 0.115, y = 588, size = 9, label = paste("italic(n)"), parse=TRUE) +
  annotate('text', x = 0.33, y = 590, size = 9, label = paste("=", dhfr_natCount)) +
  scale_y_continuous(expand = c(0,0)) +
  expand_limits(y = 700) +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 18))
ggsave(paste0('../../figs/splicemod/dhfr/splicemod_dhfr_index_hist', plot_format), gg, width = 5, height = 5, dpi = 100)


# DHFR