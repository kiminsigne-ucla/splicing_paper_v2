### Replicate graphs and bin distribution ### 

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

pkgs <- c('dplyr', 'tidyr', 'ggplot2', 'cowplot', 'ggExtra', 'grid', 'Unicode', 'weights')
load_pkgs(pkgs)

options(stringsAsFactors = F, warn = -1, warnings = -1)
plot_format_main <- '.tiff'
plot_format <- '.png'
hi_res <- 600
lo_res <- 300

data <- read.table('../../processed_data/splicemod/splicemod_data_clean.txt', 
                   sep = '\t', header = T, 
                   colClasses = c('sub_id' = 'character'))

###############################################################################

# Index across SMN1 biological replicates
corr <- signif(cor(data$index_R1_smn1, data$index_R2_smn1, use = 'p'), 2)
gg <- data %>% 
    ggplot(aes(index_R1_smn1, index_R2_smn1)) + 
    geom_point(alpha = 0.10, aes(color = replicability_smn1)) +
    scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_color_manual(values = c('black', '#0033CC')) +
    labs(x = 'inclusion index (Rep. 1)', y = 'inclusion index (Rep. 2)') +
    theme(legend.position = 'none',
          axis.title.x = element_text(size = 20, vjust = -2), 
          axis.title.y = element_text(size = 20, vjust = +4),
          axis.text.x = element_text(size = 14, color = 'grey20'),
          axis.text.y = element_text(size = 14, color = 'grey20'),
          axis.ticks.x = element_line(color = 'grey50'),
          axis.ticks.y = element_line(color = 'grey50'),
          axis.line.x = element_line(color = 'grey50'),
          axis.line.y = element_line(color = 'grey50'),
          plot.margin = unit(c(2,2,3,3),"mm")) +
    annotate('text', x = 0.15, y = 0.895, size = 7, parse = T,
             label = paste("italic(r)==", signif(corr, 3))) +
    annotate('text', x = 0.165, y = 0.80, size = 7, parse = T,
             label = paste('italic(p) < 10^-16'))
ggsave(paste0('../../figs/splicemod/smn1/splicemod_smn1_replicates', 
              plot_format_main), gg, width = 5, height = 5, dpi = hi_res)

# Index across DHFR biological replicates
corr <- signif(cor(data$index_R1_dhfr, data$index_R2_dhfr, use = 'p'), 2)
gg <- data %>% 
    ggplot(aes(index_R1_dhfr, index_R2_dhfr)) + 
    geom_point(alpha = 0.10, aes(color = replicability_dhfr)) +
    scale_color_manual(values = c('black', '#0033CC')) +
    scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    labs(x = 'inclusion index (Rep. 1)', y = 'inclusion index (Rep. 2)') +
    theme(legend.position = 'none',
          axis.title.x = element_text(size = 20, vjust = -2), 
          axis.title.y = element_text(size = 20, vjust = +4),
          axis.text.x = element_text(size = 14, color = 'grey20'),
          axis.text.y = element_text(size = 14, color = 'grey20'),
          axis.ticks.x = element_line(color = 'grey50'),
          axis.ticks.y = element_line(color = 'grey50'),
          axis.line.x = element_line(color = 'grey50'),
          axis.line.y = element_line(color = 'grey50'),
          plot.margin = unit(c(2,2,3,3),"mm")) +
    annotate('text', x = 0.15, y = 0.895, size = 7, parse = T,
             label = paste("italic(r)==", signif(corr, 3))) +
    annotate('text', x = 0.165, y = 0.80, size = 7, parse = T,
             label = paste('italic(p) < 10^-16'))
ggsave(paste0('../../figs/splicemod/dhfr/splicemod_dhfr_replicates', 
              plot_format_main), gg, width = 5, height = 5, dpi = hi_res)

# Index across intron backbones
corr <- signif(cor(data$index_smn1, data$index_dhfr, use = 'p'), 2)
gg <- ggplot(data %>% 
               filter(replicability_dhfr == 'high', 
                      replicability_smn1 == 'high'), 
               aes(index_smn1, index_dhfr)
             ) + geom_point(alpha = 0.10) + 
    scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    labs(x = 'inclusion index (SMN1)', y = 'inclusion index (DHFR)') +
    theme(legend.position = 'none',
          axis.title.x = element_text(size = 20, vjust = -2), 
          axis.title.y = element_text(size = 20, vjust = +4),
          axis.text.x = element_text(size = 14, color = 'grey20'),
          axis.text.y = element_text(size = 14, color = 'grey20'),
          axis.ticks.x = element_line(color = 'grey50'),
          axis.ticks.y = element_line(color = 'grey50'),
          axis.line.x = element_line(color = 'grey50'),
          axis.line.y = element_line(color = 'grey50'),
          plot.margin = unit(c(2,2,3,3),"mm")) + 
    annotate('text', x = 0.15, y = 0.895, size = 7, parse = T,
             label = paste("italic(r)==", signif(corr, 3))) +
    annotate('text', x = 0.165, y = 0.80, size = 7, parse = T,
             label = paste('italic(p) < 10^-16'))
ggsave(paste0('../../figs/splicemod/both/splicemod_smn1_dhfr_replicates', 
              plot_format_main), gg, width = 5, height = 5, dpi = hi_res)
ggsave(paste0('../../figs/splicemod/both/splicemod_smn1_dhfr_replicates', 
              plot_format), gg, width = 5, height = 5, dpi = lo_res)

###############################################################################
# Read distribution across GFP bins
###############################################################################
bin_labels <- c('GFP+', 'GFPint', 'GFP-')   

# SMN1 intron backbone
read_dist <- data %>% 
    filter(rep_quality == 'high') %>%
    rowwise() %>% 
    select(id, DP_R1_norm_smn1, INT_R1_norm_smn1, SP_R1_norm_smn1) %>% 
    mutate(DP = DP_R1_norm_smn1 / 
             (DP_R1_norm_smn1 + INT_R1_norm_smn1 + SP_R1_norm_smn1) * 100,
           INT = INT_R1_norm_smn1 / 
             (DP_R1_norm_smn1 + INT_R1_norm_smn1 + SP_R1_norm_smn1) * 100,
           SP = SP_R1_norm_smn1 / 
             (DP_R1_norm_smn1 + INT_R1_norm_smn1 + SP_R1_norm_smn1) * 100) %>% 
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
          plot.margin = unit(c(1,0,0,0.5), "in")) +
    labs(y = 'index')

tiff(paste0('../../figs/splicemod/smn1/splicemod_smn1_bin_dist.tiff'),
    width = 13, height = 4, units = 'in', res = hi_res)
tiff(paste0('../../figs/splicemod/smn1/splicemod_smn1_bin_dist.tiff'),
     width = 13, height = 4, units = 'in', res = hi_res)
g <- rbind(ggplotGrob(gg2), ggplotGrob(gg1), size = 'last')
id <- g$layout$t[g$layout$name == "panel"]
g$heights[id] <- unit(c(1, 2.5), 'null')
grid.newpage()
grid.draw(g)
dev.off()

# DHFR intron backbone
read_dist <- data %>%
    filter(rep_quality == 'high') %>%
    rowwise() %>%
    select(id, DP_R1_norm_dhfr, INT_R1_norm_dhfr, SP_R1_norm_dhfr) %>%
    mutate(DP = DP_R1_norm_dhfr / 
             (DP_R1_norm_dhfr + INT_R1_norm_dhfr + SP_R1_norm_dhfr) * 100,
           INT = INT_R1_norm_dhfr / 
             (DP_R1_norm_dhfr + INT_R1_norm_dhfr + SP_R1_norm_dhfr) * 100,
           SP = SP_R1_norm_dhfr / 
             (DP_R1_norm_dhfr + INT_R1_norm_dhfr + SP_R1_norm_dhfr) * 100) %>%
    select(-DP_R1_norm_dhfr, -INT_R1_norm_dhfr, -SP_R1_norm_dhfr) %>%
    gather(key = 'bin', value = 'bin_percent', DP:SP) %>%
    # add index for graph order
    left_join(select(data, id, index_R1_dhfr), by = 'id')

gg1 <- ggplot(read_dist, aes(x = reorder(id, index_R1_dhfr), y = bin)) +
    geom_tile(aes(fill = bin_percent)) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          legend.position = 'none',
          plot.margin = unit(c(1,0,0,0.5), "in")) +
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
    width = 13, height = 4, units = 'in', res = lo_res)
g <- rbind(ggplotGrob(gg2), ggplotGrob(gg1), size = 'last')
id <- g$layout$t[g$layout$name == "panel"]
g$heights[id] <- unit(c(1, 2.5), 'null')
grid.newpage()
grid.draw(g)
dev.off()

###############################################################################
# Histogram of inclusion index (wild-type human exons)
###############################################################################

# SMN1 intron backbone

smn1_natCount <- data %>% filter(sub_id == '000', 
                                 !is.na(index_smn1), 
                                 !is.na(chr), 
                                 replicability_smn1 == 'high') %>% nrow()
dhfr_natCount <- data %>% filter(sub_id == '000', 
                                 !is.na(index_dhfr), 
                                 !is.na(chr), 
                                 replicability_dhfr == 'high') %>% nrow()

gg <- data %>% 
  filter (sub_id == '000', 
          !is.na(index_smn1), 
          replicability_smn1 == 'high') %>%
  ggplot(aes(index_smn1)) + 
  geom_histogram(binwidth = 0.02, fill = "grey", color = 'black', size = 0.1) +
  geom_histogram(data = subset(data %>% filter (sub_id == '000', 
                                     !is.na(index_smn1), 
                                     replicability_smn1 == 'high'),
                             index_smn1 >= 0.79 & index_smn1 <= 1.0), 
                 binwidth = 0.02, 
                 color = 'black', 
                 fill = 'grey',
                 size = 0.1) +
  geom_vline(xintercept = 0.79, color = 'red') +
  labs(x = 'inclusion index', title = 'Wild-type human exons') +
  annotate('text', x = 0.215, y = 558, size = 5, 
           label = paste("italic(n)"), parse = TRUE) +
  annotate('text', x = 0.425, y = 560, size = 5, 
           label = paste("=", smn1_natCount)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  expand_limits(y = 650) +
  theme(axis.title.y = element_text(size = 18, vjust = 30),
        axis.title.x = element_text(size = 18, vjust = -30),
        axis.text.y = element_text(size = 12, color = 'grey20'),
        axis.text.x = element_text(size = 12, color = 'grey20'),
        plot.title = element_text(size = 16, hjust = 0.55, vjust = 2),
        axis.ticks.x = element_line(color = 'grey50'),
        axis.ticks.y = element_line(color = 'grey50'),
        plot.margin = unit(c(2,2,2,2),"mm"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = 'grey50'))
ggsave(paste0('../../figs/splicemod/smn1/splicemod_smn1_index_hist', 
              plot_format), 
       gg, width = 3, height = 4, dpi = lo_res)

# DHFR intron backbone

gg <- data %>%
  filter (sub_id == '000', 
          !is.na(index_dhfr), 
          replicability_dhfr == 'high') %>%
  ggplot(aes(index_dhfr)) + 
  geom_vline(xintercept = 0.79, color = 'red') +
  geom_histogram(binwidth = 0.02, fill = "grey", color = 'black', size = 0.1) +
  geom_histogram(data = subset(data %>%
                               filter (sub_id == '000', 
                                       !is.na(index_dhfr), 
                                       replicability_dhfr == 'high'),
                               index_dhfr > 0.79 & index_dhfr <= 1.0), 
                 binwidth = 0.02, 
                 color = 'black', 
                 fill = 'grey',
                 size = 0.1) +
  # geom_vline(xintercept = 0.79, color = 'red') +
  # geom_rect(aes(xmin = 0.8, xmax = 1.02, ymin = 0, 
  # ymax = 675), alpha = 0.1, fill = "lightgray") +
  labs(x = 'inclusion index', title = 'Wild-type human exons') +
  annotate('text', x = 0.135, y = 628, size = 5, 
           label = paste("italic(n)"), parse = TRUE) +
  annotate('text', x = 0.345, y = 630, size = 5, 
           label = paste("=", dhfr_natCount)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  expand_limits(y = 700) +
  theme(axis.title.x = element_text(size = 18, vjust = -30),
        axis.title.y = element_text(size = 18, vjust = 30),
        axis.text.y = element_text(size = 12, color = 'grey20'),
        axis.text.x = element_text(size = 12, color = 'grey20'),
        plot.title = element_text(size = 16, hjust = 0.55, vjust = 2, face = "plain"),
        axis.ticks.x = element_line(color = 'grey50'),
        axis.ticks.y = element_line(color = 'grey50'),
        plot.margin = unit(c(2,4,2,2),"mm"),
        panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, color = 'grey50'))
                                 
ggsave(paste0('../../figs/splicemod/dhfr/splicemod_dhfr_index_hist', 
              plot_format_main), 
       gg, width = 3, height = 4, dpi = hi_res)
