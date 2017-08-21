load_pkgs <- function(pkgs){
    new_pkgs <- pkgs[!(pkgs %in% installed.packages()[, 'Package'])]
    if(length(new_pkgs)) install.packages(new_pkgs)
    for(pkg in pkgs){
        suppressWarnings(suppressMessages(library(pkg, character.only = T)))
    }
}

pkgs <- c('dplyr', 'tidyr')
load_pkgs(pkgs)

options(stringsAsFactors = F, warn = -1, warnings = -1)

ref <- read.table('../../ref/splicemod/splicemod_ref_formatted_converted.txt', sep = '\t',
                  header = T, colClasses = c('sub_id' = 'character'))

smn1 <- read.csv('../../processed_data/splicemod/smn1/smn1_all_alignments.csv') %>% 
    rename(DP_R1 = R1.DP, INT_R1 = R1.INT, PS_R1 = R1.PRESORT, SP_R1 = R1.SP, 
           DP_R2 = R2.DP, INT_R2 = R2.INT, PS_R2 = R2.PRESORT, SP_R2 = R2.SP)

dhfr <- read.csv('../../processed_data/splicemod/dhfr/dhfr_all_alignments.csv') %>% 
    rename(DP_R1 = R1.DP, INT_R1 = R1.INT, PS_R1 = R1.PRESORT, SP_R1 = R1.SP, DN_R1 = R1.DN,
           DP_R2 = R2.DP, INT_R2 = R2.INT, PS_R2 = R2.PRESORT, SP_R2 = R2.SP, DN_R2 = R2.DN)

# normalize to number of million reads in each sample
smn1 <- smn1 %>%
    select(-id) %>% 
    mutate_all(funs(norm = . / (sum(.) / 1000000))) %>% 
    bind_cols(select(smn1, id), .)

dhfr <- dhfr %>% 
    select(-id) %>% 
    mutate_all(funs(norm = . / (sum(.) / 1000000))) %>% 
    bind_cols(select(dhfr, id), .)

# proportion of cells sorted into each bin
smn1_bin_prop <- c(0.320, 0.057, 1, 0.443, 0.318, 0.057, 1, 0.430)
dhfr_bin_prop <- c(0.024, 0.346, 0.077, 1, 0.363, 0.027, 0.365, 0.078, 1, 0.388)

# multiply each bin count by bin proportion
smn1 <- bind_cols(select(smn1, header = id, DP_R1:SP_R2), 
                  data.frame(mapply(`*`, select(smn1, DP_R1_norm:SP_R2_norm), 
                                    smn1_bin_prop, SIMPLIFY = FALSE)))

dhfr <- bind_cols(select(dhfr, header = id, DN_R1:SP_R2),
                  data.frame(mapply(`*`, select(dhfr, DN_R1_norm:SP_R2_norm),
                                    dhfr_bin_prop, SIMPLIFY = FALSE)))

# read sum across bins
smn1 <- smn1 %>% 
    mutate(R1_sum = DP_R1 + INT_R1 + SP_R1,
           R2_sum = DP_R2 + INT_R2 + SP_R2)

dhfr <- dhfr %>% 
    mutate(R1_sum = DP_R1 + INT_R1 + SP_R1,
           R2_sum = DP_R2 + INT_R2 + SP_R2)

# filter for low read count
low_read_threshold <- 5
smn1 <- smn1 %>% 
    filter(R1_sum >= low_read_threshold, R2_sum >= low_read_threshold)

dhfr <- dhfr %>% 
    filter(R1_sum >= low_read_threshold, R2_sum >= low_read_threshold)

# combine
data <- full_join(smn1, dhfr, by = 'header', suffix = c('_smn1', '_dhfr'))

# join reference 
data <- data %>% 
    rowwise %>% 
    mutate(id = unlist(strsplit(header, split = ' '))[1]) %>% 
    select(-header) %>% 
    left_join(ref, by = 'id') %>% 
    # reorder
    select(id:end_hg38, DP_R1_smn1:SP_R2_norm_dhfr) %>% 
    ungroup() %>% 
    mutate(id = ifelse(grepl('_', id), id, paste(ensembl_id, sub_id, sep = '_'))) %>% 
    arrange(ensembl_id, sub_id)

# convert blanks to NA
data[data == ''] <- NA

write.table(data, '../../processed_data/splicemod/splicemod_data_clean.txt', sep ='\t', row.names = F)
