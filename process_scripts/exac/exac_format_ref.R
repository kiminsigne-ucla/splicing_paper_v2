# Format ExAC reference into distinct fields and convert genomic coordinates
# to lastest version

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

###############################################################################
# File formatting
###############################################################################
ref <- read.table('../../ref/exac/exac_ref_all.txt', sep = '\t', 
                  col.names = c('header', 'sequence')) %>% 
    # remove spaces after equal signs and colons
    mutate(header = gsub(': ', ':', gsub('= ', '=', header))) %>% 
    # separate fields
    separate(header, c('id', 'region', 'strand', 'length', 'ref_allele', 'alt_allele', 
                       'snp_position', 'vcf_id', 'rel_position'), sep = ' ') %>% 
    separate(region, c('chr', 'region'), sep = ':') %>%
    separate(region, c('start', 'end'), sep = '-', fill = 'right', convert = T) %>%
    separate(id, c('ensembl_id', 'sub_id'), sep = '_', remove = F) %>%
    # get rid of leftover field identifiers
    mutate(strand = gsub('strand=', '', strand),
           length = gsub('len=', '', length),
           ref_allele = gsub('ref:', '', gsub('ref=', '', ref_allele)),
           alt_allele = gsub('alt=', '', alt_allele),
           # make sure this is numeric
           start = as.numeric(start)) %>%
    # finally, separate length into intron-exon-intron lengths
    extract(length, c("intron1_len","exon_len","intron2_len"),
            "([[:alnum:]]+).([[:alnum:]]+).([[:alnum:]]+)", convert = T)

ref <- ref %>% 
    # the BRK controls have different formatting for reference alleles, 
    # separate these and deal with separately
    filter(endsWith(id, 'BRK')) %>% 
    # controls do not have SNP positions or VCF ids but have an additional set of reference and 
    # alternate alleles which occupy these columns. Combine this information into the reference 
    # and alternate allele columns
    mutate(snp_position = gsub(';', '', gsub('ref=', '', snp_position)),
           ref_allele = paste0(ref_allele, snp_position),
           snp_position = NA,
           vcf_id = gsub('alt=', '', vcf_id),
           alt_allele = paste0(alt_allele, ';', vcf_id),
           vcf_id = NA) %>% 
    # combine back with reference
    bind_rows(filter(ref, !endsWith(id, 'BRK'))) %>% 
    # now re-format and convert the rest of the necessary columns
    mutate(snp_position = as.numeric(gsub('pos=', '', snp_position)),
           rel_position = as.numeric(gsub('rel_pos=', '', rel_position)),
           vcf_id = as.numeric(gsub('vcf-id=', '', vcf_id)))

write.table(ref, '../../ref/exac/exac_ref_formatted.txt', sep = '\t',
            row.names = F, col.names = F, quote = F)

###############################################################################
# convert genome coordinates
###############################################################################
system('bash ./exac_ref_liftover.sh')

# join
ref <- ref %>% 
    left_join(read.table('../../processed_data/exac/exac_ref_liftover.bed', sep = '\t',
                         col.names = c('chr', 'start_hg38', 'end_hg38', 'id')) %>% select(-chr), by = 'id') %>% 
    left_join(read.table('../../processed_data/exac/exac_snp_liftover.bed', sep = '\t',
                         col.names = c('chr', 'snp_position_hg38', 'snp_position_end_hg38', 'id')) %>% select(-chr), by = 'id')

write.table(ref, '../../ref/exac/exac_ref_formatted_converted.txt', sep = '\t',
            quote = F, row.names = F)
