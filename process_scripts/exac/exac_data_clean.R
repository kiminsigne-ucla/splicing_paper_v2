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

###############################################################################
# Read in data
###############################################################################

# ExAC first sequencing run, three bins
# select sort 2 only
exac_v1 <- read.csv('../../processed_data/exac/exac_v1_all_alignments.csv') %>% 
    select(-ends_with('S1')) %>% 
    rename(DP_R1 = DP1S2, DP_R2 = DP2S2, INT_R1 = INT1S2, INT_R2 = INT2S2,
           PS_R1 = R1.PS, PS_R2 = R2.PS, SP_R1 = SP1S2, SP_R2 = SP2S2)

exac_v1 <- exac_v1 %>%
    select(-id) %>% 
    mutate_all(funs(norm = . / (sum(.) / 1000000))) %>% 
    bind_cols(select(exac_v1, id), .)

# proportion of cells in each bin (relative to pre-sort)
# DP1S2, DP2S2, INT1S2, INT2S2, PS_R1, PS_R2, SP1S2, SP2S2
bin_prop <- c(0.092, 0.091, 0.127, 0.122, 1, 1, 0.596, 0.603)

# multiply each bin count by bin proportion
exac_v1 <- bind_cols(select(exac_v1, header = id, DP_R1:SP_R2), 
                     data.frame(mapply(`*`, select(exac_v1, DP_R1_norm:SP_R2_norm), 
                                       bin_prop, SIMPLIFY = FALSE)))


# ExAC second sequencing run, four bins
exac_v2 <- read.csv('../../processed_data/exac/exac_v2_all_alignments.csv')
exac_v2 <- exac_v2 %>% 
    select(Hi.R1:Lo.R2) %>% 
    mutate_all(funs(norm = . / (sum(., na.rm = T) / 1000000))) %>% 
    bind_cols(select(exac_v2, id), .)

# Hi.R1, Hi.R2, IntHi.R1, IntHi.R2, IntLo.R1, IntLo.R2, Lo.R1, Lo.R2
bin_prop <- c(0.070, 0.075, 0.029, 0.032, 0.032, 0.033, 0.227, 0.252)

# multiply each bin count by bin proportion
exac_v2 <- bind_cols(select(exac_v2, header = id, Hi.R1:Lo.R2), 
                     data.frame(mapply(`*`, select(exac_v2, Hi.R1_norm:Lo.R2_norm), 
                                       bin_prop, SIMPLIFY = FALSE)))

print(paste("Initial number of sequences (v1, v2):", nrow(exac_v1), nrow(exac_v2)))

###############################################################################
# Filtering
###############################################################################

hi_read_threshold <- 10
# nat_read_threshold <- 500

rep_agreement <- 0.20

# read filter
exac_v1 <- exac_v1 %>%
    mutate(v1_R1_sum = DP_R1 + INT_R1 + SP_R1,
           v1_R2_sum = DP_R2 + INT_R2 + SP_R2) %>%
    filter(v1_R1_sum >= hi_read_threshold , v1_R2_sum >= hi_read_threshold)

exac_v2 <- exac_v2 %>%
    mutate(v2_R1_sum = Hi.R1 + IntHi.R1 + IntLo.R1 + Lo.R1,
           v2_R2_sum = Hi.R2 + IntHi.R2 + IntLo.R2 + Lo.R2) %>% 
    # read filter
    filter(v2_R1_sum >= hi_read_threshold , v2_R2_sum >= hi_read_threshold)

print(paste("Number of sequences after read filter (v1, v2):", nrow(exac_v1), nrow(exac_v2)))

# index agreement between replicates filter
exac_v1 <- exac_v1 %>% 
    # calculate index
    mutate(v1_index_R1 = (DP_R1_norm * 0 + INT_R1_norm * 0.85 + SP_R1_norm * 1) / 
               (DP_R1_norm + INT_R1_norm + SP_R1_norm),
           v1_index_R2 = (DP_R2_norm * 0 + INT_R2_norm * 0.85 + SP_R2_norm * 1) / 
               (DP_R2_norm + INT_R2_norm + SP_R2_norm)) %>%
    # rep agreement
    filter(abs( v1_index_R1 - v1_index_R2 ) <= rep_agreement)

exac_v2 <- exac_v2 %>% 
    mutate(v2_index_R1 = (Hi.R1_norm * 0 + IntHi.R1_norm * 0.80 + IntLo.R1_norm * 0.95 + Lo.R1_norm * 1) / 
               (Hi.R1_norm + IntHi.R1_norm + IntLo.R1_norm + Lo.R1_norm), 
           v2_index_R2 = (Hi.R2_norm * 0 + IntHi.R2_norm * 0.80 + IntLo.R2_norm * 0.95 + Lo.R2_norm * 1) / 
               (Hi.R2_norm + IntHi.R2_norm + IntLo.R2_norm + Lo.R2_norm),
           v2_R1_norm_sum = Hi.R1_norm + IntHi.R1_norm + IntLo.R1_norm + Lo.R1_norm,
           v2_R2_norm_sum = Hi.R2_norm + IntHi.R2_norm + IntLo.R2_norm + Lo.R2_norm) %>%
    # rep agreement
    filter(abs(v2_index_R1 - v2_index_R2) <= rep_agreement)

print(paste("Number of sequences after index filter (v1, v2):", nrow(exac_v1), nrow(exac_v2)))

###############################################################################
# join data
###############################################################################
data_all <- full_join(exac_v1, exac_v2, by = 'header') %>% 
    # small substitutions so separate will work easier
    mutate(header = gsub('strand= ', 'strand=', header),
           header = gsub('>', '', header)) %>%
    separate(header, into = c('id', 'chr', 'strand', 'length'), sep = ' ') %>% 
    select(-chr, -strand, -length)

data_all$v1_index <- rowMeans(select(data_all, v1_index_R1, v1_index_R2))
data_all$v2_index <- rowMeans(select(data_all, v2_index_R1, v2_index_R2))

# Correlation between v1 and v2
fit <- summary(lm(v2_index ~ v1_index, data_all))$r.squared
gg <- ggplot(data_all, aes(v1_index, v2_index)) + geom_point(alpha = 0.25) +
    geom_smooth(method = 'lm', color = 'red') +
    labs(x = 'ExAC V1', y = 'ExAC V2') +
    annotate('text', x = 0.95, y = 0.10, parse = T, 
             label = paste0('R^2==', signif(fit, 3)), size = 5)

ggsave(paste0('../../figs/exac/exac_v1_v2_replicates', plot_format), gg,
       width = 6, height = 6, dpi = 300)

# read in updated ref
ref <- read.table('../../ref/exac/exac_ref_formatted_converted_flipped.txt', 
                  sep='\t', header=T)

# combine with data
data_all <- left_join(data_all, ref, by = 'id') %>% 
    arrange(ensembl_id, sub_id) %>% 
    # update the category to either control, natural, or mutant
    # ifelse() structure: ifelse(condition, action if true, action if false)
    mutate(category = ifelse(endsWith(id, '000'), 'natural', 'mutant'), 
           category = ifelse(endsWith(id, 'BRK'), 'control', category),
           category = ifelse(endsWith(id, 'SKP'), 'skipped_exon', category),
           category = ifelse(startsWith(ensembl_id, 'RANDOM-EXON'), 'random_exon', category)) %>%
    # get rid of misc. data columns
    select(id, ensembl_id:category, v1_R1_sum:v1_index_R2, v1_index, 
           v2_R1_sum:v2_index_R2, v2_index)

# replace NaN with NA
data_all[data_all == 'NaN'] <- NA

###############################################################################
# filtering on natural exon inclusion 
###############################################################################

# keep SKP and RANDOM-EXON categories
data_other <- data_all %>% 
    group_by(ensembl_id) %>% 
    filter(any(sub_id == 'SKP') | (any(ensembl_id == 'RANDOM-EXON')) ) %>%
    ungroup()

data <- data_all %>% 
    group_by(ensembl_id) %>% 
    # mutant must have corresponding natural that passed previous filters
    filter(any(sub_id == '000')) %>% 
    # calculate dPSI
    mutate(nat_v1_index_R1 = v1_index_R1[sub_id == '000'],
           nat_v1_index_R2 = v1_index_R2[sub_id == '000'],
           nat_v2_index_R1 = v2_index_R1[sub_id == '000'], 
           nat_v2_index_R2 = v2_index_R2[sub_id == '000'],
           v1_dpsi_R1 = v1_index_R1 - nat_v1_index_R1,
           v1_dpsi_R2 = v1_index_R2 - nat_v1_index_R2,
           v2_dpsi_R1 = v2_index_R1 - nat_v2_index_R1,
           v2_dpsi_R2 = v2_index_R2 - nat_v2_index_R2,
           nat_v1_index = (nat_v1_index_R1 + nat_v1_index_R2) / 2, 
           nat_v2_index = (nat_v2_index_R1 + nat_v2_index_R2) / 2,
           nat_seq = original_seq[sub_id == '000']) %>% 
    filter(abs(nat_v2_index_R1 - nat_v2_index_R2) <= rep_agreement) %>%
    ungroup()

# get averages between replicates
data <- data %>%
    rowwise() %>%
    mutate(v1_dpsi = mean(c(v1_dpsi_R1, v1_dpsi_R2)),
           v2_dpsi = mean(c(v2_dpsi_R1, v2_dpsi_R2)),
           delta_dpsi = abs(v1_dpsi - v2_dpsi) ) %>%
    ungroup()

dpsi_threshold <- -0.50
dpsi_threshold_stringent <- -0.70

data <- data %>% 
    mutate(strong_lof = ifelse(v2_dpsi <= dpsi_threshold, TRUE, FALSE),
           stringent_lof = ifelse(v2_dpsi <= dpsi_threshold_stringent, TRUE, FALSE))

write.table(data, '../../processed_data/exac/exac_data_clean.txt', sep = '\t',
            row.names = F, quote = F)
