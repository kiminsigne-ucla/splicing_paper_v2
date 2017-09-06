load_pkgs <- function(pkgs){
    new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, 'Package'])]
    if(length(new_pkgs)) install.packages(new_pkgs)
    for(pkg in pkgs){
        suppressWarnings(suppressMessages(library(pkg, character.only = T)))
    }
}

pkgs <- c('dplyr', 'tidyr', 'ggplot2', 'cowplot', 'forcats', 'gridExtra', 'grid')
load_pkgs(pkgs)

options(stringsAsFactors = F, warn = -1, warnings = -1)

plot_format <- '.png'

# custom color palette
source("../color_palette.R")
# Specify new color palette
steps <- c("blue2", "cyan", "white", "yellow", "red2")
pal <- color.palette(steps, c(160,1,1,160), space = "rgb")


data <- read.table('../../processed_data/exac/exac_data_clean.txt',
                   sep = '\t', header = T)
data_snps <- data %>% 
    filter(category == 'mutant')

dpsi_threshold <- -0.50

###############################################################################
# Figure 3A, rank ordered splicing
###############################################################################

# calculate number of mutants per exon background
data_num_mut_per_exon <- data_snps %>% 
    count(ensembl_id) %>% 
    rename(num_muts_per_exon = n) %>% 
    mutate(ensembl_id = fct_reorder(ensembl_id, num_muts_per_exon, .desc = T)) %>% 
    arrange(desc(num_muts_per_exon))

# top panel, number of mutants per exon background
gg1 <-
    data_num_mut_per_exon %>% 
    ggplot(aes(x = ensembl_id, y = num_muts_per_exon, group = 1)) + geom_line() + 
    labs(x = '', y = 'SNV count') +
    theme_classic() + 
    theme(
        plot.margin = unit(c(1,0,0,0), units = "lines"),
        panel.grid = element_blank(),
        plot.title = element_text(size=16, face="bold", margin = margin(50, 0, 0, 50)),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=20),
        axis.line = element_line(colour = "grey20", size = 0.25),
        axis.text.y = element_text(colour="grey20",size=18,angle=0,hjust=.5,vjust=.5) ,
        axis.text.x = element_blank(),
        axis.ticks.x= element_blank(),
        legend.background = element_rect(),
        legend.title = element_text(size=18, face="bold"),
        legend.key.size = unit(8, "mm"),
        legend.key.width = unit(4, "mm"),
        legend.text = element_text(size=18),
        legend.position = 'none')


# jitter formatting function
delta.psi.jitter <- function(data, x, y) {
    ggplot(data, aes_string(x = x, y = y), na.rm = TRUE) + 
        geom_jitter(aes(color = nat_v2_index), alpha = 0.10) + 
        scale_colour_gradientn(limits = c(-0.005,1), 
                               breaks = seq(0, 1, by = 0.25), 
                               colors = pal(321),
                               expression(index["WT "])) +
        scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.5)) + 
        xlab("") + 
        ylab(expression(paste(Delta, ' inclusion index'))) +
        theme_classic() + 
        theme(
            plot.margin = unit(c(0.0,0.5,0.2,0.2), units = "lines"),
            panel.grid = element_blank(),
            plot.title = element_text(size = 16, margin = margin(0, 0, 10, 10)),
            axis.title.y = element_text(size = 18),
            axis.title.x = element_text(size = 20),
            axis.text.y = element_text(colour = "grey20",size = 18,angle = 0,
                                       hjust = 0.5,vjust = 0.5 ),
            axis.line = element_line(colour = "grey20", size = 0.25),
            axis.text.x = element_text(colour = "grey20", size = 24, angle = 0,
                                       hjust = 0.5, margin = margin(5, 0, 0, 0)),
            axis.ticks.x = element_blank()
        )  
}

# Figure 3A, bottom panel, change in exon inclusion index for mutant library,
# sorted by decreasing number of mutants per exon background
gg2 <- data_snps %>% 
    inner_join(data_num_mut_per_exon, by = 'ensembl_id') %>% 
    mutate(ensembl_id = fct_reorder(ensembl_id, num_muts_per_exon, .desc = T)) %>%
    delta.psi.jitter("ensembl_id", "v2_dpsi") + 
    theme(axis.text.x = element_blank() ) + 
    xlab("Wild-type exon ID") + 
    theme(legend.position = "none") +
    geom_hline(yintercept = dpsi_threshold, linetype = "dashed", color = "grey10") 

# combine into one graph
g <- rbind(ggplotGrob(gg1), ggplotGrob(gg2), size = 'last')
id <- g$layout$t[g$layout$name == "panel"]
g$heights[id] <- unit(c(0.6, 2.4), 'null')
grid.newpage()
grid.draw(g)

ggsave( 
    paste0('../../figs/exac/exac_rank_order_splicing_3A', plot_format), 
    plot = g,
    width = 4, height = 6, units = 'in'
)

###############################################################################
# Figure 3B, index and conservation
###############################################################################
### index, binned by relative position for boxplot ###
data$label_renamed <- factor(data$label, 
                             levels=c("upstr_intron", "exon", "downstr_intron"), 
                             labels=c("5' intron", "exon", "3' intron"))

data <- data %>%
    mutate(rel_pos_binned = cut(rel_position_scaled, 
                                breaks = seq(-0.80, 1.80, 0.01)))

group.colors <- c("5' intron" = "#b90c0d", # orange 
                  "exon" = "black",#838383",   # belize
                  "3' intron" = "#b90c0d") #b90c0d") # orangefac364

index_boxplot <- data %>% 
    filter(category == "mutant") %>% 
    filter(!is.na(rel_pos_binned)) %>% 
    ggplot(aes(rel_pos_binned, v2_dpsi)) + 
    geom_boxplot(aes(color = label_renamed), outlier.size = 0.10, 
                 outlier.colour = "black", notch = FALSE) +
    theme(legend.position = "none",
          axis.text.x = element_blank(), 
          axis.title.y = element_text(margin = margin(0, 0, -65, -65), size = 20),
          axis.title.x =  element_text(margin = margin(0, 0, -45, -65), size = 20),
          axis.text.y = element_text(colour = "grey20", size = 16),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(x = '', 
         y = expression(paste(Delta, ' inclusion index')),
         color = '') + 
    scale_color_manual(values = group.colors) + ylim(-1, 1)

ggsave(paste0('../../figs/exac/exac_index_boxplot', plot_format),
       width = 11, height = 4, units = 'in')

#### index, average index at each binned relative position, heatmap tile ###
index_tile_with_legend <- data %>%
    filter(category == "mutant") %>%
    filter(!is.na(rel_pos_binned)) %>% 
    group_by(rel_pos_binned) %>% 
    summarise(mean_dpsi_per_rel_pos = mean(v2_dpsi, na.rm = T)) %>% 
    ggplot(aes(x = rel_pos_binned, y = 0.5)) + 
    geom_tile(aes(fill = mean_dpsi_per_rel_pos)) +
    viridis::scale_fill_viridis(option = "viridis", direction = -1, 
                                breaks = seq(-1, 0.1, 0.25), 
                                limits = c(-1, 0.10)) +
    labs(x = '', y = '', fill = '') +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.line = element_blank(),
          axis.title.y = element_text(margin = margin(0,0,-65,-65), size = 20),
          axis.title.x =  element_text(margin = margin(0,0,-45,-65), size = 20),
          plot.margin = unit(c(0, 0, 0, 0), "in"))

# save legend separately
g_legend <- function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}
legend <- g_legend(index_tile_with_legend)
tiff('../../figs/exac/fig3b_legend.tiff', width = 350)
grid.newpage()
grid.draw(legend)
dev.off()

index_tile <- index_tile_with_legend + theme(legend.position = 'none')

ggsave(paste0('../../figs/exac/exac_index_tile', plot_format),
    width = 11, height = 0.8, units = 'in')



### phastCons scores for all natural sequences in library
create_positions <- function(ensembl_id, chr, start, end, strand) {
    all_positions = data.frame()
    if(strand == '-'){ counter = 170}
    else { counter = 1}
    for(i in seq(start, end + 1)) {
        if(counter > 170 | counter <= 0) { next }
        id = paste0(ensembl_id, '_pos', counter)
        if(strand == '-') {counter = counter - 1}
        else { counter = counter + 1}
        all_positions <- bind_rows(all_positions,
                                   data.frame(chr = chr, start = i, 
                                              end = i + 1, id = id))
    }
    return(all_positions)
}

# for each natural exon background, generate data frame with each genomic position
# of the intron-exon-intron construct
if(!file.exists('../../processed_data/exac/exac_nat_cons_scaled.rds')) {
    data %>% 
        filter(category == 'natural') %>% 
        group_by(ensembl_id) %>% 
        do(data.frame(create_positions(.$ensembl_id, .$chr, .$start_hg38, 
                                       .$end_hg38, .$strand))) %>% 
        ungroup() %>% 
        select(-ensembl_id) %>% 
        write.table('../../processed_data/exac/exac_nat_positions.bed',
                    sep = '\t', quote = F, row.names = F)
    
    # grab conservation at each base
    system(paste('bash',
                 '../run_phastCons.sh',
                 '../../processed_data/exac/exac_nat_positions.bed', 
                 '../../processed_data/exac/exac_nat_cons_scores_all.bed'))
    
    # read in data
    exac_nat_cons <- read.table('../../processed_data/exac/exac_nat_cons_scores_all.bed', 
                                sep = '\t', header = F,
                                col.names = c('name', 'size', 'bases_covered', 
                                              'snp_sum', 'mean0', 'mean_cons_score')) %>%
        filter(bases_covered != 0) %>% 
        separate(name, into = c('ensembl_id', 'position'), sep = '_', remove = F) %>%
        mutate(rel_position = as.numeric(gsub('pos', '', position))) %>% 
        select(-(position:mean0), id = name) %>% 
        # join other information necessary to calculate relative position
        left_join(select(data, ensembl_id, intron1_len, 
                         exon_len, intron2_len, strand) %>% 
                      distinct(), 
                  by = 'ensembl_id')
    
    rel_position <- function(intron1_len, exon_len, intron2_len, strand, rel_position) {
        in_interval <- function(start, end, x) {
            if (start <= x & x <= end) { return(TRUE) }
            else {return(FALSE)}
        }
        if (strand == '-' | strand == -1) {
            # set appropriate lengths
            upstr_intron_len <- intron2_len
            downstr_intron_len <- intron1_len
        }
        else{
            upstr_intron_len <- intron1_len
            downstr_intron_len <- intron2_len
        }
        
        regions <- data.frame(label = c('upstr_intron', 'exon', 'downstr_intron'),
                              start = c(1, upstr_intron_len + 1, upstr_intron_len + exon_len + 1),
                              end = c(upstr_intron_len, upstr_intron_len + exon_len, 
                                      upstr_intron_len + exon_len + downstr_intron_len))
        # select region the SNP falls in
        label <- regions %>% 
            rowwise() %>% 
            filter(in_interval(start, end, x = rel_position))
        
        start <- label$start[1]
        end <- label$end[1]
        label <- label$label[1]
        
        # get distance from end of upstream intron/exon boundary
        boundary <- upstr_intron_len
        distance <- rel_position - boundary
        
        # normalize to feature length
        if (label == 'downstr_intron') {
            scaled_distance <- 1 + (distance - exon_len) / (end - start)
        }
        else {
            scaled_distance <- distance / (end - start + 1)
        }
        
        # additionally, find distance from respective intron/exon boundary
        # left side of boundary is negative, right side is positive
        if (label == 'upstr_intron') {
            rel_pos_feature <- rel_position - end - 1
        }
        if (label == 'exon') {
            # closer to downstream intron
            if ( rel_position - start + 1 >= end - rel_position + 1) {
                rel_pos_feature <- rel_position - end - 1 # negative position, left side of boundary
            }
            else {
                # closer to upstream intron, right side of boundary, positive
                rel_pos_feature <- rel_position - start + 1
            }
        }
        if (label == 'downstr_intron') {
            rel_pos_feature <- rel_position - start + 1
        }
        return(paste(label, rel_pos_feature, scaled_distance, sep = ':'))
    }
    
    exac_nat_cons <- exac_nat_cons %>% 
        select(intron1_len, exon_len, intron2_len, strand, rel_position, id) %>%
        na.omit() %>%
        rowwise() %>%
        mutate(rel_position_info = rel_position(intron1_len, exon_len, 
                                                intron2_len, strand, rel_position)) %>%
        separate(rel_position_info, 
                 c('label', 'rel_position_feature', 'rel_position_scaled'), 
                 sep = ':', convert = T) %>%
        select(id, label, rel_position_feature, rel_position_scaled) %>%
        left_join(exac_nat_cons, ., by = 'id') %>% distinct()
    
    exac_nat_cons <- exac_nat_cons %>% 
        mutate(rel_pos_binned = cut(rel_position_scaled, 
                                    breaks = seq(-.80, 1.80, 0.01)))
    # save as RDS so factor is retained
    saveRDS(exac_nat_cons, '../../processed_data/exac/exac_nat_cons_scaled.rds')
} else {
    exac_nat_cons <- readRDS('../../processed_data/exac/exac_nat_cons_scaled.rds')
}

nat_cons_tile <- exac_nat_cons %>%
    filter(!is.na(rel_pos_binned)) %>% 
    group_by(rel_pos_binned) %>%
    summarise(mean_cons = mean(mean_cons_score)) %>%
    ggplot(aes(x = rel_pos_binned, y = 0.5)) + geom_tile(aes(fill = mean_cons)) +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          legend.position = 'none',
          axis.line = element_blank(),
          axis.title.y = element_text(margin = margin(0, 0, -65, -65), size = 20),
          axis.title.x =  element_text(margin = margin(0, 0, -45, -65), size = 20),
          plot.margin = unit(c(0, 0, 0, 0), "in")) +
    labs(y = '', x = '', fill = '') +
    viridis::scale_fill_viridis()

ggsave(paste0("../../figs/exac/exac_phastCons_nat_tile", plot_format), 
       width = 11, height = 0.8, units = 'in')


### SNP density tile ###
ref <- read.table('../../ref/exac/exac_ref_formatted_converted.txt',
                  sep = '\t', header = T)

total_snps <- ref %>%
    filter(sub_id != '000', sub_id != 'BRK') %>%
    nrow()

snp_density_tile <- ref %>%
    mutate(rel_pos_binned = cut(rel_position_scaled, breaks = seq(-.80, 1.80, 0.01))) %>%
    filter(!is.na(rel_pos_binned), sub_id != '000', sub_id != 'BRK') %>%
    group_by(rel_pos_binned) %>%
    summarise(snp_density = n() / total_snps) %>%
    ggplot(aes(rel_pos_binned, 0.5)) + geom_tile(aes(fill = snp_density)) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), 
          legend.position = 'none',
          axis.line = element_blank(),
          # axis.title.y = element_blank(),
          # axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(0,0,-65,-65), size = 20),
          axis.title.x =  element_text(margin = margin(0,0,-45,-65), size = 20),
          plot.margin = unit(c(0,0,0,0),'in')) +
    viridis::scale_fill_viridis() +
    labs(x = '', y = '', fill = '')

ggsave(paste0('../../figs/exac/exac_snv_density', plot_format),
    width = 11, height = 0.8, units = 'in')




