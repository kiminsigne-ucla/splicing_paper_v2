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

# read in data
data <- read.table('../../processed_data/exac/exac_data_clean.txt',
                   sep = '\t', header = T)

###############################################################################
# extract conservation for each SNP
###############################################################################
# output bed file of SNP positions to get conservation
write.table(file = '../../processed_data/exac/snp_positions.bed', 
            x = data %>% 
                select(chr, snp_position_hg38, snp_position_end_hg38, id) %>% 
                na.omit(), 
            sep = '\t', col.names = F, row.names = F, quote = F)

system(paste('bash',
             '../run_phastCons.sh',
             '../../processed_data/exac/snp_positions.bed',
             '../../processed_data/exac/snp_position_cons_scores.bed'))

# read in phastCons
phastCons <- read.table('../../processed_data/exac/snp_position_cons_scores.bed', 
                        sep = '\t', header = F,
                        col.names = c('name', 'size', 'bases_covered', 
                                      'snp_sum', 'mean0', 'mean_cons_score')) %>% 
    filter(bases_covered != 0)

data <- data %>% 
    left_join(select(phastCons, id = name, mean_cons_score), by = 'id')

###############################################################################
# grab ExAC annotation for all sequences
###############################################################################
ref <- read.table('../../ref/exac/exac_ref_formatted_converted.txt',
                  sep = '\t', header = T)
# write BED file of SNP positions
write.table(ref %>%
                mutate(chr = gsub('chr', '', chr),
                       # ExAC built on hg19, supply region so tabix works properly
                       snp_region = paste0(chr, ':', snp_position, '-',  snp_position + 1)) %>% 
                select(snp_region) ,
            file = '../../processed_data/exac/tabix_input_snp_regions.txt',
            quote = F, row.names = F, col.names = F, sep = '\t')

# system(paste('while read line; do tabix',
#              '../../ref/exac/ExAC.r0.3.1.sites.vep.vcf.gz',
#              '$line >> ../../ref/exac/snp_exac_annot_all.txt;',
#              'done < ../../processed_data/exac/tabix_input_snp_regions.txt'))

# read in annotation and parse
exac_annot <- read.table('../../ref/exac/snp_exac_annot_all.txt', 
                         sep = '\t', header = F)
colnames(exac_annot) <- c('chr', 'position', 'snp_id', 'ref_allele', 
                          'alt_allele', 'quality', 'filter', 'info')

# grab allele frequency column separately
exac_annot <- exac_annot %>% 
    separate(info, into = c('other', 'AF'), sep = 'AF=', remove = F) %>% 
    select(-other) %>% 
    separate(AF, into = c('AF', 'other'), sep = ';', convert = T) %>% 
    select(-other)

exac_annot <- exac_annot %>% 
    # only care about CSQ field, split by it and get rid of everything before
    separate(info, into = c('other', 'CSQ'), sep = 'CSQ=') %>% 
    select(-other) %>% 
    # get rid of all fields after CSQ
    separate(CSQ, into = c('CSQ', 'other'), sep = ';') %>% 
    select(-other)

# parse CSQ column, consequence annotations from Ensembl VEP (Variant Effect 
# Predictor). Description of fields taken from header of vcf file
csq_fields <- 'Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|ASN_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF_info|LoF_flags|LoF_filter|LoF|context|ancestral'
csq_fields <- unlist(strsplit(csq_fields, split = '[|]'))
exac_annot <- separate(exac_annot, CSQ, into = csq_fields, sep = '[|]')

exac_annot <- exac_annot %>% 
    # give each alternate allele its own entry
    mutate(alt_allele = strsplit(alt_allele, split = ','),
           AF = strsplit(AF, split = ',')) %>%
    rename(annot = Consequence) %>% 
    distinct(chr, position, ref_allele, alt_allele, annot, .keep_all = T) %>%
    # make sure alt allele matches allele the annotation is describing
    filter(alt_allele == Allele) %>% 
    mutate(AF = as.numeric(AF)) %>%
    # rearrange
    select(chr:Allele, annot, IMPACT:AF)


# join exac data
data <- data %>%
    mutate(chr = gsub('chr', '', chr)) %>% 
    left_join(select(exac_annot, -ref_allele, -alt_allele, alt_allele = Allele),
              by = c('chr', 'snp_position' = 'position', 'alt_allele')) %>% 
    distinct(id, annot, .keep_all = T)

###############################################################################
# re-score variants using Ensembl VEP
###############################################################################
# write VCF file for Ensembl VEP
data %>% 
    filter(category == 'mutant') %>% 
    dplyr::select(chr, snp_position_hg38, id, ref_allele, alt_allele) %>% 
    write.table(file = '../../ref/exac/exac.vcf', 
                sep='\t', quote = F, row.names = F, col.names = F)

# system(paste('../ensembl-vep/vep -i ../../ref/exac/exac.vcf',
#              '--cache --dir_cache ../../ref/vep/'))
# system('mv variant_effect_output.txt ../../processed_data/exac/exac_vep.txt')
# system('rm -f variant_effect_output.txt_summary.html')

# join updated Ensembl exon IDs so we can join the Ensembl VEP annotation
updated_ids <- read.table('../../ref/exon_ids_updated.txt', 
                          sep = '\t', header = T)

# exons can fall in more than one transcript, annotations are dependent on both
# transcript and exon ID
exon_data <- data %>%
    select(id, chr, snp_position_hg38, ref_allele, alt_allele, 
           exon_id_old = ensembl_id) %>% 
    left_join(select(updated_ids, exon_id_old, exon_id_new = new_exon_id,
                     ensembl_gene_id, ensembl_transcript_id, is_constitutive, 
                     rank), by = c('exon_id_old'))

vep <- read.table('../../processed_data/exac/exac_vep.txt', sep = '\t', comment.char = '#')
colnames(vep) <- c('id', 'location', 'alt_allele', 'ensembl_gene_id', 
                   'ensembl_transcript_id', 'feature_type', 'consequence', 
                   'cDNA_position', 'CDS_position', 'protein_position', 
                   'amino_acids', 'codons', 'existing_variation', 'extra')

exon_data <- exon_data %>% 
    left_join(select(vep, id, alt_allele, ensembl_transcript_id, feature_type, 
                     consequence, extra), 
              by = c('id', 'alt_allele', 'ensembl_transcript_id')) %>%
    distinct()

# some entries have multiple consequences, create new row for each
exon_data <- exon_data %>% 
    mutate(consequence = strsplit(consequence, ',')) %>% 
    unnest(consequence)

# let's try to filter it so there is one transcript per exon. For now, let's 
# just choose the transcript where the consequence is the most severe
# convert to factor, order levels according to severity can be found here:
# (http://www.ensembl.org/info/genome/variation/predicted_data.html)
consequence_levels <- c('transcript_ablation', 'splice_acceptor_variant', 
                        'splice_donor_variant', 'stop_gained', 'frameshift_variant', 
                        'stop_lost', 'start_lost', 'transcript_amplification', 
                        'inframe_insertion', 'inframe_deletion', 'missense_variant', 
                        'protein_altering_variant', 'splice_region_variant',
                        'incomplete_terminal_codon_variant', 'stop_retained_variant',
                        'synonymous_variant', 'coding_sequence_variant', 
                        'mature_miRNA_variant', '5_prime_UTR_variant', 
                        '3_prime_UTR_variant', 'non_coding_transcript_exon_variant', 
                        'intron_variant', 'NMD_transcript_variant',
                        'non_coding_transcript_variant', 'upstream_gene_variant', 
                        'downstream_gene_variant', 'TFBS_ablation', 
                        'TFBS_amplification', 'TF_binding_site_variant',
                        'regulatory_region_ablation', 'feature_elongation', 
                        'regulatory_region_variant', 'feature_truncation', 
                        'intergenic_variant')
exon_data$consequence <- factor(exon_data$consequence, levels = consequence_levels)

exon_data_filtered <- exon_data %>% 
    group_by(id) %>% 
    arrange(consequence) %>%
    filter(row_number() == 1) %>% 
    distinct(id, consequence, .keep_all = T) %>% 
    ungroup() %>% 
    arrange(id)

# join to original data
data <- data %>%
    left_join(select(exon_data_filtered, -(chr:ref_allele), -exon_id_old), 
              by = c('id', 'alt_allele')) 

###############################################################################
# CADD annotation
###############################################################################
# submitted VCF file (exac.vcf) at http://cadd.gs.washington.edu/score to 
# score variants

cadd_annot <- read.table('../../processed_data/exac/snp_cadd_annot.txt', 
                         sep = '\t', header = F, 
                         col.names = c('chr', 'snp_position', 'ref_allele', 
                                       'alt_allele', 'cadd_score', 'phred')) %>%
    distinct()

data <- data %>%
    left_join(dplyr::select(cadd_annot, -ref_allele), 
              by = c('chr', 'snp_position', 'alt_allele'))
