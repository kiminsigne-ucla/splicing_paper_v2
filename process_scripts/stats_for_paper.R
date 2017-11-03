library(dplyr)
library(tidyr)

data <- read.table('../processed_data/splicemod/splicemod_data_clean.txt',
                   sep = '\t', header = T) %>% 
    filter(rep_quality == 'high')

### statistical tests ###
# compare ESE changes to random
data %>% 
    mutate(category = ifelse(grepl('ESE', category), 'perturb_ESE', category),
           category = ifelse(grepl('rnd_exon', category), 'rnd_exon', category)) %>% 
    filter(category == 'perturb_ESE' | category == 'rnd_exon') %>% 
    wilcox.test(dpsi_smn1 ~ category, data = ., alternative = 'less')

# compare ESS changes to random
data %>% 
    mutate(category = ifelse(grepl('ESS', category), 'perturb_ESS', category),
           category = ifelse(grepl('rnd_exon', category), 'rnd_exon', category)) %>% 
    filter(category == 'perturb_ESS' | category == 'rnd_exon') %>% 
    wilcox.test(dpsi_smn1 ~ category, data = ., alternative = 'greater')

# strongest ESE
data %>% 
    mutate(category = ifelse(grepl('rnd_exon', category), 'rnd_exon', category)) %>% 
    filter(category == 'clst_Ke2011_ESE' | category == 'rnd_exon') %>% 
    wilcox.test(dpsi_smn1 ~ category, data = ., alternative = 'less')

# strongest ESS
data %>% 
    mutate(category = ifelse(grepl('rnd_exon', category), 'rnd_exon', category)) %>% 
    filter(category == 'clst_Ke2011_ESS' | category == 'rnd_exon') %>% 
    wilcox.test(dpsi_smn1 ~ category, data = ., alternative = 'greater')

# intron changes
data %>% 
    mutate(category = ifelse(grepl('ICS', category), 'intron', category),
           category = ifelse(grepl('rnd_intron', category), 'rnd_intron', category)) %>% 
    filter(category == 'intron' | category == 'rnd_intron') %>% 
    wilcox.test(dpsi_smn1 ~ category, data = .)

# label variants as exonic or intronic
data_labelled <- data %>% 
    filter(seq_type == 'mut') %>% 
    separate(loc, into = c('loc_start', 'loc_end'), sep = ':', remove = F) %>% 
    mutate(exon_start = intron1_len,
           exon_end = intron1_len + exon_len,
           loc_start = as.numeric(loc_start),
           loc_end = as.numeric(loc_end),
           loc_start_stranded = ifelse(strand == -1, 170 - loc_end + 1, loc_start + 1),
           loc_end_stranded = ifelse(strand == -1, 170 - loc_start - 1, loc_end - 1)) %>% 
    select(id, category, strand, dpsi_smn1, exon_start, exon_end, loc_start, loc_end,
           loc_start_stranded, loc_end_stranded) %>% 
    na.omit() %>% 
    rowwise() %>% 
    mutate(exon_overlap_len = length(intersect(seq(loc_start_stranded, loc_end_stranded), 
                                               seq(exon_start, exon_end))),
           overlaps_exon = ifelse(exon_overlap_len > 0, T, F))

# how many RBP variants are completely intronic?
data_labelled %>% filter(category == 'RBPmats') %>% count(overlaps_exon)
# compare intronic RBP variants to random
data_labelled %>% 
    mutate(category = ifelse(grepl('rnd_intron', category), 'rnd_intron', category)) %>% 
    filter(overlaps_exon == F) %>%
    filter(category == 'RBPmats' | category == 'rnd_intron') %>% 
    wilcox.test(dpsi_smn1 ~ category, data = ., alternative = 'less')

data_labelled %>% 
    mutate(category = ifelse(grepl('rnd_exon', category), 'rnd_exon', category)) %>% 
    filter(overlaps_exon == T) %>%
    filter(category == 'RBPmats' | category == 'rnd_exon') %>% 
    wilcox.test(dpsi_smn1 ~ category, data = ., alternative = 'less')

data_labelled %>% 
    mutate(category = ifelse(grepl('rnd_intron', category), 'rnd_intron', category)) %>% 
    filter(overlaps_exon == F) %>% 
    filter(category == 'variation' | category == 'rnd_intron') %>% 
    wilcox.test(dpsi_smn1 ~ category, data = .)

data_labelled %>% 
    mutate(category = ifelse(grepl('rnd_exon', category), 'rnd_exon', category)) %>% 
    filter(overlaps_exon == T) %>%
    filter(category == 'variation' | category == 'rnd_exon') %>% 
    wilcox.test(dpsi_smn1 ~ category, data = ., alternative = 'less')

data_labelled %>% 
    mutate(strong_lof = ifelse(dpsi_smn1 <= -0.50, T, F)) %>% 
    count(strong_lof, overlaps_exon)

data_labelled %>% 
    mutate(strong_lof = ifelse(dpsi_smn1 <= -0.50, T, F)) %>% 
    filter(category == 'RBPmats') %>% 
    count(strong_lof)

data %>% 
    filter(seq_type == 'mut', HAL_bin != 'same') %>% 
    wilcox.test(dpsi_smn1 ~ HAL_bin, data = .)

snv_data <- read.table('../processed_data/exac/exac_func_annot.txt',
                       sep='\t', header = T)

snv_data <- snv_data %>% 
    mutate(cons_bin = cut(mean_cons_score, breaks = c(0, 0.5, 1.0), 
                          include.lowest = T,  labels = c('low', 'high')))

# are intronic SDVs enriched for high conservation
snv_data %>% 
    filter(grepl('intron', label), !is.na(cons_bin), !is.na(strong_lof)) %>% 
    mutate(cons_bin = factor(cons_bin, levels = c('high', 'low')),
           strong_lof = factor(strong_lof, levels = c(T, F))) %>%
    select(strong_lof, cons_bin) %>% 
    table() %>% 
    fisher.test()

# are exonic SDVs enriched for high conservation
snv_data %>% 
    filter(label == 'exon', !is.na(cons_bin), !is.na(strong_lof)) %>% 
    mutate(cons_bin = factor(cons_bin, levels = c('high', 'low')),
           strong_lof = factor(strong_lof, levels = c(T, F))) %>%
    select(strong_lof, cons_bin) %>% 
    table() %>% 
    fisher.test()

# is conservation enriched in introns for all SDVs
snv_data %>% 
    filter(strong_lof == T) %>% 
    mutate(label = ifelse(grepl('intron', label), 'intron', label),
           label = factor(label, levels = c('intron', 'exon')),
           cons_bin = factor(cons_bin, levels = c('high', 'low'))) %>% 
    select(cons_bin, label) %>% 
    table() %>% 
    fisher.test(alternative = 'greater')

# is conservation enriched for Lof vs. not LoF for upstream introns
snv_data %>% 
    filter(label == 'downstr_intron') %>% 
    mutate(cons_bin = factor(cons_bin, levels = c('high', 'low')),
           strong_lof = factor(strong_lof, levels = c(T, F))) %>% 
    select(cons_bin, strong_lof) %>% 
    table() %>% 
    fisher.test(alternative = 'greater')

# is conservation enriched in exons
snv_data %>% 
    mutate(label = ifelse(grepl('intron', label), 'intron', label),
           cons_bin = factor(cons_bin, levels = c('high', 'low'))) %>% 
    select(cons_bin, label) %>% 
    table() %>% 
    fisher.test(alternative = 'greater')

# difference in AF
snv_data <- snv_data %>% 
    mutate(AF_bin = cut(AF, 
                        breaks = c(0, 0.00001, 0.000025, 0.0001, 0.001, 1), 
                        include.lowest = T,
                        labels = c('Singleton', 'AC = 2-3', 
                                   'AC = 4-10', '0.01%', '>0.1%' ))) 

kruskal.test(v2_dpsi ~ AF_bin, snv_data)
