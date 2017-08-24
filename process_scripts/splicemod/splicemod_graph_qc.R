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

pkgs <- c('dplyr', 'tidyr', 'ggplot2', 'cowplot', 'ggExtra', 'grid')
load_pkgs(pkgs)

options(stringsAsFactors = F, warn = -1, warnings = -1)
plot_format <- '.png'

data <- read.table('../../processed_data/splicemod/splicemod_data_clean.txt', sep = '\t', header = T, 
                   colClasses = c('sub_id' = 'character'))

###############################################################################

# Index across SMN1 replicates
corr <- signif(cor(data$index_R1_smn1, data$index_R2_smn1, use = 'p'), 3)
data %>% 
    ggplot(aes(index_R1_smn1, index_R2_smn1)) + geom_point(alpha = 0.10) +
    geom_smooth(method = 'lm', color = 'red') +
    scale_color_manual(values = c('black', 'red')) +
    labs(x = 'inclusion index (Rep. 1)', y = 'inclusion index (Rep. 2)') +
    theme(legend.position = 'None', text = element_text(size = 20), axis.text = element_text(size = 16)) +
    annotate('text', x = 0.85, y = 0.10, size = 6,  label = paste('r = ', corr))
ggsave(paste0('../../figs/splicemod/smn1/splicemod_smn1_replicates', plot_format), width = 5, height = 5, dpi = 100)

# Index across DHFR replicates
corr <- signif(cor(data$index_R1_dhfr, data$index_R2_dhfr, use = 'p'), 3)
data %>% 
    ggplot(aes(index_R1_dhfr, index_R2_dhfr)) + geom_point(alpha = 0.10) +
    geom_smooth(method = 'lm', color = 'red') +
    scale_color_manual(values = c('black', 'red')) +
    labs(x = 'inclusion index (Rep. 1)', y = 'inclusion index (Rep. 2)') +
    theme(legend.position = 'None', text = element_text(size = 20), axis.text = element_text(size = 16)) +
    annotate('text', x = 0.85, y = 0.10, size = 6,  label = paste('r = ', corr))
ggsave(paste0('../../figs/splicemod/dhfr/splicemod_dhfr_replicates', plot_format), width = 5, height = 5, dpi = 100)

# Index across SMN1 and DHFR
corr <- signif(cor(data$index_smn1, data$index_dhfr, use = 'p'), 3)
ggplot(data, aes(index_smn1, index_dhfr)) + geom_point(alpha = 0.10) +
    # scale_color_manual(values = c('black', 'grey')) +
    geom_abline(linetype = 'dotted', slope = 1, intercept = 0.30) +
    geom_abline(linetype = 'dotted', slope = 1, intercept = -0.30) +
    geom_smooth(method = 'lm', color = 'red') + 
    labs(x = 'inclusion index (SMN1)', y = 'inclusion index (DHFR)') +
    theme(text = element_text(size = 20), axis.text = element_text(size = 16),
          legend.position = 'none') +
    annotate('text', x = 0.85, y = 0.10, size = 6,  label = paste('r = ', corr))
ggsave(paste0('../../figs/splicemod/splicemod_smn1_dhfr_replicates', plot_format), width = 5, height = 5, dpi = 100)

###############################################################################
# Read distribution across bins
###############################################################################
bin_labels <- c('GFP+', 'GFPint', 'GFP-')   
# SMN1
read_dist <- data %>% 
    filter(abs(index_R1_smn1 - index_R2_smn1) <= 0.30) %>%
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
    geom_bar(stat = 'identity', color = 'lightgrey') +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = unit(c(0,0,0, 0.5), "in")) +
    labs(y = 'index')

# # line
# gg2 <- data %>% 
#     filter(abs(index_R1_smn1 - index_R2_smn1) <= 0.30) %>%
#     ggplot(aes(x = reorder(id, index_R1_smn1), y = index_R1_smn1, group = 1)) + 
#     geom_line(color = 'black') + 
#     theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
#           axis.title.x = element_blank(),
#           plot.margin = unit(c(0,0,0, 0.5), "in")) +
#     labs(y = 'index')

png(paste0('../../figs/splicemod/smn1/splicemod_smn1_bin_dist', plot_format),
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
    filter(abs(index_R1_dhfr - index_R2_dhfr) <= 0.30) %>%
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

png(paste0('../../figs/splicemod/dhfr/splicemod_dhfr_bin_dist', plot_format),
    width = 13, height = 4, units = 'in', res = 300)
g <- rbind(ggplotGrob(gg2), ggplotGrob(gg1), size = 'last')
id <- g$layout$t[g$layout$name == "panel"]
g$heights[id] <- unit(c(1, 2.5), 'null')
grid.newpage()
grid.draw(g)
dev.off()