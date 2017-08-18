# Format splicemod reference into distinct fields and convert genomic coordinates
# to lastest version

suppressMessages(library(dplyr))
suppressWarnings(suppressMessages(library(tidyr)))

options(stringsAsFactors = F, warn = -1, warnings = -1)

############################################################################### 
# File formatting 
###############################################################################
# read in ref, format ID field
ref <- read.table('../../ref/splicemod/splicemod_ref.txt', sep = '\t', 
                  col.names = c('header', 'seq')) %>% 
    rowwise() %>% 
    mutate(id = unlist(strsplit(header, split = ' '))[[1]]) %>% 
    distinct(id, .keep_all = T) %>% 
    # small substitutions to make separating easier
    mutate(header = gsub('cat= ', 'cat=', header))

# format nat and mut ref separately
nat_ref <- ref %>% 
    filter(!grepl('mut', header)) %>% 
    select(-id) %>% 
    separate(header, into = c('id', 'chr', 'strand', 'length', 'ccds'), sep = ' ') %>% 
    mutate(seq_type = 'nat')

mut_ref <- ref %>% 
    filter(grepl('mut', header)) %>% 
    select(-id) %>% 
    # don't modify equal sign for control sequences because they don't have any information for those fields
    mutate(header = ifelse(grepl('CTRL', header), header, gsub('= ', '=', header)),
           header = gsub(', ', ',', header),
           header = gsub(';', '', header)) %>% 
    separate(header, c('id', 'seq_type', 'chr', 'strand', 'length', 'ccds', 'category', 
                       'loc', 'str', 'mutations', 'num_changes', 'orig_score',
                       'new_score', 'criteria_met', 'max_iter_reached', 
                       'no_new_mutants'), sep = ' ')

# more reference formatting 
ref <- bind_rows(nat_ref, mut_ref) %>% 
    separate(chr, c('chr', 'region'), sep = ':') %>%
    separate(region, c('start', 'end'), sep = '-', fill = 'right', convert = T) %>%
    separate(id, c('ensembl_id', 'sub_id'), sep = '_', remove = F) %>% 
    # if natural sequence and no sub id, add 000
    mutate(sub_id = ifelse(is.na(sub_id), '000', sub_id)) %>% 
    # get rid of leftover field identifiers
    mutate(strand = gsub('strand=', '', strand),
           length = gsub('len=', '', length),
           chr = gsub('chr', '', chr),
           ccds = gsub('ccds=', '', ccds),
           category = gsub('cat=', '', category),
           loc = gsub('loc=', '', loc),
           mutations = gsub('set=', '', mutations),
           num_changes = as.numeric(gsub('changes=', '', num_changes)),
           orig_score = as.numeric(gsub('orig_score=', '', orig_score)),
           new_score = as.numeric(gsub('new_score=', '', new_score)),
           criteria_met = as.logical(gsub('criteria_met=', '', criteria_met)),
           max_iter_reached = as.logical(gsub('MaxIterReached=', '', max_iter_reached)),
           no_new_mutants = as.logical(gsub('NoNewMutants=', '', no_new_mutants))) %>% 
    # separate mutant range into separate columns
    separate(loc, c('mut_start', 'mut_end'), sep = ':', remove = F, convert = T) %>% 
    # update seq type for controls
    mutate(seq_type = ifelse(grepl('CTRL', ensembl_id), 'control', seq_type),
           seq_type = ifelse(grepl('RAND', ensembl_id), 'control', seq_type)) %>% 
    # finally, separate length into intron-exon-intron lengths
    extract(length, c("intron1_len","exon_len","intron2_len"),
            "([[:alnum:]]+).([[:alnum:]]+).([[:alnum:]]+)", convert = T) %>% 
    arrange(ensembl_id, sub_id)

write.table(ref, '../../ref/splicemod/splicemod_ref_formatted.txt', sep = '\t', 
row.names = F, col.names = F, quote = F)

###############################################################################
# Convert genomic coordinates
###############################################################################

system('bash ./splicemod_ref_liftover.sh')

###############################################################################
# Join converted coordinates
###############################################################################
ref <- ref %>% 
    left_join(read.table('../../processed_data/splicemod/splicemod_ref_liftover.bed', sep = '\t',
      header = F, col.names = c('chr', 'start_hg38', 'end_hg38', 'id'))  %>%
                  select(-chr), by = 'id')

write.table(ref, '../../ref/splicemod/splicemod_ref_formatted_converted.txt', 
  sep = '\t', row.names = F, col.names = T, quote = F)