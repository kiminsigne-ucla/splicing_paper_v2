load(exac_intron_cons.RData)

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

data <- read.table('../../processed_data/exac/exac_data_clean.txt', 
                   sep = '\t', header = T)

###############################################################################
# calculate genomic coordinates for the upstream and downstream introns for
# conservation
###############################################################################
intron_coords <- data %>% 
    filter(category == 'mutant') %>% 
    distinct(ensembl_id, .keep_all = T) %>% 
    mutate(intron1_start = start_hg38_0based, 
           intron1_end = start_hg38_0based + intron1_len - 1,
           intron2_start = start_hg38_0based + intron1_len + exon_len, 
           intron2_end = end_hg38_0based) %>% 
    select(id, ensembl_id, chr, strand, intron1_len, intron2_len, intron1_start:intron2_end)

intron_coords %>% 
    mutate(upstr_intron = ifelse(strand == '+', 
                                 paste(intron1_start, intron1_end, sep = '-'),
                                 paste(intron2_start, intron2_end, sep = '-')),
           downstr_intron = ifelse(strand == '+', 
                                   paste(intron2_start, intron2_end, sep = '-'),
                                   paste(intron1_start, intron1_end, sep = '-'))) %>% 
    select(chr, upstr_intron, downstr_intron, ensembl_id) %>% 
    gather('feature_type', 'region', upstr_intron:downstr_intron) %>% 
    separate(region, into = c('start', 'end'), sep = '-') %>% 
    mutate(id = paste(ensembl_id, feature_type, sep = '-')) %>% 
    select(chr, start, end, id) %>% 
    write.table(file = '../../processed_data/exac/exac_intron_coords.bed', 
                sep = '\t', quote = F, row.names = F, col.names = F)

system(paste('bash',
             '../run_phastCons.sh',
             '../../processed_data/exac/exac_intron_coords.bed',
             '../../processed_data/exac/exac_intron_cons_scores_all.bed'))

intron_cons <- read.table('../../processed_data/exac/exac_intron_cons_scores_all.bed', 
                          sep = '\t', header = F,
                          col.names = c('name', 'size', 'bases_covered', 
                                        'snp_sum', 'mean0', 'mean_cons_score')) %>% 
    filter(bases_covered != 0) %>% 
    separate(name, into = c('ensembl_id', 'feature_type'), sep = '-') %>% 
    mutate(feature_type = paste0(feature_type, '_mean_cons')) %>% 
    spread(feature_type, mean_cons_score)

data <- data %>% 
    left_join(select(intron_cons, ensembl_id, upstr_intron_mean_cons) %>% 
                  na.omit(),
              by = 'ensembl_id') %>% 
    left_join(select(intron_cons, ensembl_id, downstr_intron_mean_cons) %>% 
                  na.omit(), 
              by = 'ensembl_id')

###############################################################################
# Let's get coordinates for 100bp of the flanking intron instead of shorter 
# version we synthesized.
###############################################################################
intron_coords %>% 
    mutate(intron1_start = intron1_start - (100 - intron1_len + 1),
           intron2_end = intron2_end + (100 - intron2_len + 1),
           upstr_intron_100 = ifelse(strand == '+', 
                                     paste(intron1_start, intron1_end, sep = '-'),
                                     paste(intron2_start, intron2_end, sep = '-')),
           downstr_intron_100 = ifelse(strand == '+', 
                                       paste(intron2_start, intron2_end, sep = '-'),
                                       paste(intron1_start, intron1_end, sep = '-'))) %>% 
    select(chr, upstr_intron_100, downstr_intron_100, ensembl_id) %>% 
    gather('feature_type', 'region', upstr_intron_100:downstr_intron_100) %>% 
    separate(region, into = c('start', 'end'), sep = '-') %>% 
    mutate(id = paste(ensembl_id, feature_type, sep = '-')) %>% 
    select(chr, start, end, id) %>% 
    write.table(file = '../../processed_data/exac/exac_intron_coords_100.bed', 
                sep = '\t', quote = F, row.names = F, col.names = F)

system(paste('bash',
             '../run_phastCons.sh',
             '../../processed_data/exac/exac_intron_coords_100.bed', 
             '../../processed_data/exac/exac_intron_cons_scores_100_all.bed'))

intron_cons_100 <- read.table('../../processed_data/exac/exac_intron_cons_scores_100_all.bed', 
                              sep = '\t', header = F,
                              col.names = c('name', 'size', 'bases_covered', 
                                            'snp_sum', 'mean0', 'mean_cons_score')) %>% 
    filter(bases_covered != 0) %>% 
    separate(name, into = c('ensembl_id', 'feature_type'), sep = '-') %>% 
    mutate(feature_type = paste0(feature_type, '_mean_cons')) %>% 
    spread(feature_type, mean_cons_score)

data <- data %>% 
    left_join(select(intron_cons_100, ensembl_id, upstr_intron_100_mean_cons) %>% 
                  na.omit(), 
              by = 'ensembl_id') %>% 
    left_join(select(intron_cons_100, ensembl_id, downstr_intron_100_mean_cons) %>% 
                  na.omit(), 
              by = 'ensembl_id')

data <- data %>% 
    mutate(upstr_intron_len = ifelse(strand == '+', intron1_len, intron2_len),
           downstr_intron_len = ifelse(strand == '+', intron2_len, intron1_len))

# supplement, synthetic short intron conservation vs. 100 bp conservation
# upstream
up <- ggplot(data, aes(upstr_intron_mean_cons, upstr_intron_100_mean_cons)) + 
    geom_point(alpha = 0.1) +
    geom_abline(intercept = 0, slope = 1, type = 'dashed') +
    labs(x = 'average phastCons score\nupstream intron (40-81bp)', 
         y = 'average phastCons score\nupstream intron (100 bp)') +
    theme(legend.position = 'none',
        axis.title.x = element_text(size = 20, vjust = -2), 
        axis.title.y = element_text(size = 20, vjust = +4),
        axis.text.x = element_text(size = 14, color = 'grey20'),
        axis.text.y = element_text(size = 14, color = 'grey20'),
        axis.ticks.x = element_line(color = 'grey50'),
        axis.ticks.y = element_line(color = 'grey50'),
        axis.line.x = element_line(color = 'grey50'),
        axis.line.y = element_line(color = 'grey50'),
        plot.margin = unit(c(2,2,3,3),"mm")) 
up

ggsave(paste0('../../figs/supplement/exac_upstr_intron_cons_short_vs_full', plot_format),
       height = 6, width = 6, unit = 'in')

ggplot(data, aes(upstr_intron_mean_cons - upstr_intron_100_mean_cons)) +
    geom_density() + 
    labs(x = 'mean phastCons upstream intron (40-81bp) -
         mean phastCons upstream intron (100 bp)')

# downstream
down <- ggplot(data, aes(downstr_intron_mean_cons, downstr_intron_100_mean_cons)) + 
    geom_point(alpha = 0.1) + 
    geom_abline(intercept = 0, slope = 1, type = 'dashed') +
    labs(x = 'average phastCons score\ndownstream intron (30-71bp)', 
         y = 'average phastCons score\ndownstream intron (100 bp)') +
    theme(legend.position = 'none',
        axis.title.x = element_text(size = 20, vjust = -2), 
        axis.title.y = element_text(size = 20, vjust = +4),
        axis.text.x = element_text(size = 14, color = 'grey20'),
        axis.text.y = element_text(size = 14, color = 'grey20'),
        axis.ticks.x = element_line(color = 'grey50'),
        axis.ticks.y = element_line(color = 'grey50'),
        axis.line.x = element_line(color = 'grey50'),
        axis.line.y = element_line(color = 'grey50'),
        plot.margin = unit(c(2,2,3,3),"mm")) 
down

ggsave(paste0('../../figs/supplement/exac_downstr_intron_cons_short_vs_full', plot_format),
       height = 6, width = 6, unit = 'in')

plot_grid(up, down, labels = "AUTO")

ggsave(paste0('../../figs/supplement/exac_upstr_downstr_intron_cons_short_vs_full', plot_format),
       height = 6, width = 12, unit = 'in')

# combined
data %>% 
    mutate(upstr_cons_diff = upstr_intron_mean_cons - upstr_intron_100_mean_cons,
           downstr_cons_diff = downstr_intron_mean_cons - downstr_intron_100_mean_cons) %>% 
    gather(key = 'intron_type', value = 'mean_cons_diff', upstr_cons_diff:downstr_cons_diff) %>% 
    ggplot(aes(mean_cons_diff)) + geom_density(aes(color = intron_type)) +
    scale_color_discrete(labels = c('downstream intron', 'upstream intron')) +
    labs(color = '',
         x = 'intron (30-81bp) mean phastCons - \nintron (100bp) mean phastCons') +
    theme(legend.position = c(0.675, 0.675),
          axis.title.x = element_text(size = 14),
          legend.text = element_text(size = 10))
    
ggsave(paste0('../../figs/supplement/exac_upstr_downstr_cons', plot_format),
       height = 3, width = 5, unit = 'in')

###############################################################################
# Let's get conservation for the intronic portions of the splice donor and 
# acceptor for each sequence. First, let's grab the locations.
###############################################################################
data %>% 
    filter(category == 'natural') %>% 
    distinct(ensembl_id, .keep_all = T) %>% 
    mutate(exon_start_hg38 = start_hg38_0based + intron1_len,
           exon_end_hg38 = exon_start_hg38 + exon_len - 1,
           acc_intron_coord = ifelse(strand == '+',
                                     paste(exon_start_hg38 - 21, exon_start_hg38 - 1, sep = '-'),
                                     paste(exon_end_hg38 + 1, exon_end_hg38 + 21, sep = '-')),
           don_intron_coord = ifelse(strand == '+',
                                     paste(exon_end_hg38 + 1, exon_end_hg38 + 7, sep = '-'),
                                     paste(exon_start_hg38 - 7, exon_start_hg38 - 1, sep = '-'))) %>% 
    select(chr, ensembl_id, acc_intron_coord, don_intron_coord) %>% 
    gather('feature_type', 'region', acc_intron_coord:don_intron_coord) %>% 
    separate(region, into = c('start', 'end'), sep = '-') %>% 
    mutate(id = paste(ensembl_id, feature_type, sep = '-')) %>% 
    select(chr, start, end, id) %>% 
    write.table(sep = '\t', row.names = F, col.names = F,
            file = '../../processed_data/exac/splice_site_intron_coords.bed', quote = F)

system(paste('bash',
             '../run_phastCons.sh',
             '../../processed_data/exac/splice_site_intron_coords.bed', 
             '../../processed_data/exac/splice_site_intron_cons.bed'))

splice_site_cons <- read.table('../../processed_data/exac/splice_site_intron_cons.bed', 
                               sep = '\t', header = F,
                               col.names = c('name', 'size', 'bases_covered', 
                                             'snp_sum', 'mean0', 'mean_cons_score')) %>% 
    filter(bases_covered != 0) %>% 
    separate(name, into = c('ensembl_id', 'feature_type'), sep = '-') %>%
    mutate(feature_type = gsub('_coord', '', feature_type),
           feature_type = paste0(feature_type, '_ss_mean_cons')) %>% 
    spread(feature_type, mean_cons_score)

data <- data %>% 
    left_join(select(splice_site_cons, ensembl_id, acc_intron_ss_mean_cons) %>% 
                  na.omit(), 
              by = 'ensembl_id') %>% 
    left_join(select(splice_site_cons, ensembl_id, don_intron_ss_mean_cons) %>% 
                  na.omit(), 
              by = 'ensembl_id')

data <- data %>% 
    mutate(junc_avg_acc = acc_intron_ss_mean_cons / upstr_intron_100_mean_cons,
           junc_avg_don = don_intron_ss_mean_cons / downstr_intron_100_mean_cons)

write.table(data, '../../processed_data/exac/exac_data_intron_cons.txt',
            sep = '\t', row.names = F, quote = F)

###############################################################################
# Conservation of in-frame vs. out-of-frame exons (whole genome)
###############################################################################
# release 89 5/17 based on hg38
ensembl <- biomaRt::useMart('ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl')
attributes <- c('ensembl_gene_id', 'description', 'chromosome_name', 
                'start_position', 'end_position', 'strand', 'ensembl_transcript_id', 
                'transcript_start', 'transcript_end', 'ensembl_exon_id', 
                'exon_chrom_start', 'exon_chrom_end', 'is_constitutive', 'rank', 
                'phase', 'end_phase')

all_exon_ids <- biomaRt::getBM(attributes = attributes, mart = ensembl)

# calculate various coordinates necessary to determine relative scaled position
# of SNVs within intron/exon
all_exon_ids <- all_exon_ids %>% 
    mutate(intron1_end = exon_chrom_start - 1,
           intron1_start = intron1_end - 99,
           intron2_start = exon_chrom_end + 1,
           intron2_end = intron2_start + 99,
           upstr_intron_start = ifelse(strand == 1, intron1_start, intron2_start),
           upstr_intron_end = ifelse(strand == 1, intron1_end, intron2_end),
           downstr_intron_start = ifelse(strand == 1, intron2_start, intron1_start),
           downstr_intron_end = ifelse(strand == 1, intron2_end, intron1_end),
           upstr_intron_len = upstr_intron_end - upstr_intron_start + 1,
           downstr_intron_len = downstr_intron_end - downstr_intron_start + 1,
           exon_len = exon_chrom_end - exon_chrom_start + 1)

inframe_outframe_exons <- all_exon_ids %>% 
    mutate(exon_type = case_when(.$phase == 0 & .$end_phase == 0 ~ 'inframe',
                                 TRUE ~ 'outframe'))

# # get all exons that start and end on phase 0
# inframe_outframe_exons <- all_exon_ids %>% 
#     filter(phase == 0, end_phase == 0) %>% 
#     mutate(exon_type = 'inframe')
# 
# # get out-frame exons and label them
# inframe_outframe_exons <- all_exon_ids %>% 
#     filter(phase != -1, end_phase != -1) %>% 
#     filter(phase != 0 & end_phase != 0) %>% 
#     mutate(exon_type = 'outframe') %>% 
#     bind_rows(inframe_outframe_exons)

inframe_outframe_exons <- inframe_outframe_exons %>% 
    distinct(ensembl_exon_id, .keep_all = T)

# distribution of exon numbers across phases
inframe_outframe_exons %>% 
    filter(phase != -1, end_phase != -1) %>% 
    arrange(phase, end_phase) %>% 
    mutate(frame = paste(phase, end_phase, sep=',')) %>% 
    ggplot(aes(frame)) + geom_bar() +
    labs(x = 'exon start phase, end phase') +
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust=1))

ggsave(paste0('../../figs/supplement/genome_exon_phase_dist', plot_format),
       height = 3, width = 4, units = 'in')

# let's get upstream intron coordinates first. We will generate a bed file
# containing position of each nucleotide in range so we can get conservation per
# nucleotide
inframe_outframe_exons %>% 
    filter(chromosome_name %in% c(1:22, 'X', 'Y')) %>%
    mutate(chromosome_name = paste0('chr', chromosome_name)) %>% 
    select(chromosome_name, upstr_intron_start, upstr_intron_end, ensembl_exon_id) %>% 
    group_by(ensembl_exon_id) %>% 
    mutate(position = list(seq(upstr_intron_start, upstr_intron_end, by = 1))) %>% 
    unnest(position) %>%
    ungroup() %>% 
    mutate(id = paste(ensembl_exon_id, position, sep = '_'),
           start = position,
           end = position + 1) %>% 
    select(chromosome_name, start, end, id) %>% 
    write.table(file = '../../processed_data/exac/nat_upstr_intron_positions.bed', 
                sep = '\t', col.names = F, row.names = F, quote = F)

# system(paste('bash',
#              '../run_phastCons.sh',
#              '../../processed_data/exac/nat_upstr_intron_positions.bed', 
#              '../../processed_data/exac/nat_upstr_intron_cons_scores_all.bed'))

# downstream intron
inframe_outframe_exons %>% 
    filter(chromosome_name %in% c(1:22, 'X', 'Y')) %>%
    mutate(chromosome_name = paste0('chr', chromosome_name)) %>% 
    select(chromosome_name, downstr_intron_start, downstr_intron_end, ensembl_exon_id) %>% 
    group_by(ensembl_exon_id) %>% 
    mutate(position = list(seq(downstr_intron_start,downstr_intron_end, by = 1))) %>% 
    unnest(position) %>%
    ungroup() %>% 
    mutate(id = paste(ensembl_exon_id, position, sep = '_'),
           start = position,
           end = position + 1) %>% 
    select(chromosome_name, start, end, id) %>% 
    write.table(file = '../../processed_data/exac/nat_downstr_intron_positions.bed', 
                sep = '\t', col.names = F, row.names = F, quote = F)

# system(paste('bash',
#              '../run_phastCons.sh',
#              '../../processed_data/exac/nat_downstr_intron_positions.bed', 
#              '../../processed_data/exac/nat_downstr_intron_cons_scores_all.bed'))

# exon
inframe_outframe_exons %>%                             
    filter(chromosome_name %in% c(1:22, 'X', 'Y')) %>%
    mutate(chromosome_name = paste0('chr', chromosome_name)) %>% 
    select(chromosome_name, exon_chrom_start, exon_chrom_end, ensembl_exon_id) %>% 
    group_by(ensembl_exon_id) %>% 
    mutate(position = list(seq(exon_chrom_start,exon_chrom_end, by = 1))) %>% 
    unnest(position) %>%
    ungroup() %>% 
    mutate(id = paste(ensembl_exon_id, position, sep = '_'),
           start = position,
           end = position + 1) %>% 
    select(chromosome_name, start, end, id) %>% 
    write.table(file = '../../processed_data/exac/nat_exon_positions.bed', 
                sep = '\t', col.names = F, row.names = F, quote = F)

# system(paste('bash',
#              '../run_phastCons.sh',
#              '../../processed_data/exac/nat_exon_positions.bed', 
#              '../../processed_data/exac/nat_exon_cons_scores_all.bed'))

# it takes awhile to read in the conservation score files and calculate summary,
# so we'll only do this if they don't already exist
if(!file.exists('../../processed_data/exac/nat_upstr_cons_summary.rds')){
    nat_upstr_cons <- read.table('../../processed_data/exac/nat_upstr_intron_cons_scores_all.bed', 
                                 sep = '\t', header = F,
                                 col.names = c('name', 'size', 'bases_covered', 
                                               'snp_sum', 'mean0', 'mean_cons_score')) %>%
        filter(bases_covered != 0) %>% 
        tidyr::separate(name, into = c('ensembl_id', 'position'), sep = '_', remove = T) %>% 
        select(-(size:mean0)) %>% 
        left_join(select(inframe_outframe_exons, ensembl_exon_id, upstr_intron_len, 
                         upstr_intron_start, upstr_intron_end, strand, exon_type),
                  by = c('ensembl_id' = 'ensembl_exon_id')) %>%
        arrange(ensembl_id, position) %>% 
        mutate(position = as.numeric(position),
               rel_position = ifelse(strand == 1, 
                                     upstr_intron_end - position,
                                     position - upstr_intron_start),
               # upstream intron, keep them all negative,
               rel_position = -1 * rel_position,
               rel_position_scaled = rel_position / upstr_intron_len,
               rel_pos_binned = cut(rel_position_scaled, breaks = seq(-1, 0, 0.01))) %>% 
        group_by(rel_pos_binned, exon_type) %>%
        summarise(mean_cons_per_rel_pos = mean(mean_cons_score, na.rm = T))
    # save as RDS so levels in factor variable are saved correctly
    saveRDS(nat_upstr_cons, '../../processed_data/exac/nat_upstr_cons_summary.rds')
} else{
    nat_upstr_cons <- readRDS('../../processed_data/exac/nat_upstr_cons_summary.rds')
}

if(!file.exists('../../processed_data/exac/nat_downstr_cons_summary.rds')) {
    nat_downstr_cons <- read.table('../../processed_data/exac/nat_downstr_intron_cons_scores_all.bed', 
                                   sep = '\t', header = F,
                                   col.names = c('name', 'size', 'bases_covered', 
                                                 'snp_sum', 'mean0', 'mean_cons_score')) %>%
        filter(bases_covered != 0) %>%
        tidyr::separate(name, into = c('ensembl_id', 'position'), sep = '_', remove = T) %>% 
        select(-(size:mean0)) %>% 
        left_join(select(inframe_outframe_exons, ensembl_exon_id, downstr_intron_len, 
                         downstr_intron_start, downstr_intron_end, strand, exon_type),
                  by = c('ensembl_id' = 'ensembl_exon_id')) %>%
        arrange(ensembl_id, position) %>% 
        mutate(position = as.numeric(position),
               rel_position = ifelse(strand == 1,
                                     position - downstr_intron_start,
                                     downstr_intron_end - position),
               rel_position_scaled = 1 + (rel_position / downstr_intron_len),
               rel_pos_binned = cut(rel_position_scaled, breaks = seq(1, 2, 0.01), 
                                    include.lowest = T)) %>% 
        group_by(rel_pos_binned, exon_type) %>%
        summarise(mean_cons_per_rel_pos = mean(mean_cons_score, na.rm = T))
    saveRDS(nat_downstr_cons, '../../processed_data/exac/nat_downstr_cons_summary.rds')
} else {
    nat_downstr_cons <- readRDS('../../processed_data/exac/nat_downstr_cons_summary.rds')
}

if(!file.exists('../../processed_data/exac/nat_exon_cons_summary.rds')) {
    nat_exon_cons <- read.table('../../processed_data/exac/nat_exon_cons_scores_all.bed', 
                                sep = '\t', header = F,
                                col.names = c('name', 'size', 'bases_covered', 
                                              'snp_sum', 'mean0', 'mean_cons_score')) %>%
        filter(bases_covered != 0) %>%
        tidyr::separate(name, into = c('ensembl_id', 'position'), sep = '_', remove = T) %>% 
        select(-(size:mean0)) %>% 
        left_join(select(inframe_outframe_exons, ensembl_exon_id, exon_len, 
                         exon_chrom_start, exon_chrom_end, strand, exon_type),
                  by = c('ensembl_id' = 'ensembl_exon_id')) %>%
        arrange(ensembl_id, position) %>% 
        mutate(position = as.numeric(position),
               rel_position = ifelse(strand == 1,
                                     position - exon_chrom_start,
                                     exon_chrom_end - position),
               rel_position_scaled = rel_position / exon_len,
               rel_pos_binned = cut(rel_position_scaled, breaks = seq(0, 1, 0.01), 
                                    include.lowest = T)) %>% 
        group_by(rel_pos_binned, exon_type) %>%
        summarise(mean_cons_per_rel_pos = mean(mean_cons_score, na.rm = T))
    saveRDS(nat_exon_cons, '../../processed_data/exac/nat_exon_cons_summary.rds')
} else {
    nat_exon_cons <- readRDS('../../processed_data/exac/nat_exon_cons_summary.rds')
}

nat_cons <- bind_rows(nat_upstr_cons, nat_exon_cons, nat_downstr_cons)
nat_cons$rel_pos_binned <- factor(nat_cons$rel_pos_binned, levels = nat_cons$rel_pos_binned)
nat_cons$exon_type <- factor(nat_cons$exon_type)
levels(nat_cons$exon_type) <- c('phase (0-0)', 'other phases')

# nat_cons %>%
#     ggplot(aes(rel_pos_binned, 0.5)) + 
#     geom_tile(aes(fill = mean_cons_per_rel_pos)) + 
#     facet_grid(exon_type ~ .) + 
#     theme(axis.text = element_blank(), axis.ticks = element_blank()) +
#     viridis::scale_fill_viridis(limits = c(0, 1)) +
#     labs(x = '', y = '', fill = 'phastCons\nscore')

# natural conservation between in-frame and out-of-frame exons, genome-wide
ggplot(nat_cons, aes(rel_pos_binned, mean_cons_per_rel_pos, color = exon_type)) +
    geom_point() + scale_color_manual(values = c('black', 'red')) +
    ylim(c(0, 1)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    theme(legend.position = c(0.85, 0.80)) +
    labs(x = 'relative scaled position', 
         y = 'average phastCons score', color = 'exon phase')

ggsave(paste0('../../figs/supplement/genome_exon_cons_inframe_vs_outframe', plot_format),
       height = 4, width = 5, unit = 'in')

t.test(mean_cons_per_rel_pos ~ exon_type, nat_cons)

save.image("../../processed_data/exac/exac_intron_cons.RData")
