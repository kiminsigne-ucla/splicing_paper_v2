# get all SNVs from ExAC annotation, this takes awhile
zgrep 'SNV' ../../ref/exac/ExAC.r0.3.1.sites.vep.vcf.gz > ../../ref/exac/ExAC_SNVs.vcf.gz

# reformat to bed format, add ID and strand column. Easier to just do this in Python
python exac_snv_vcf2bed.py ../../ref/exac/ExAC_SNVs.vcf.gz ../../ref/exac/ExAC_SNVs.bed.gz

# intersect with bed file containing exonic and intronic regions to get location
# of SNVs, print the original feature a and original feature b
intersectBed -a ../../ref/exac/ExAC_SNVS.bed.gz -b ../../ref/exon_intron_all.bed.gz -wa -wb | 
gzip > ../../ref/exac/ExAC_SNVs_genomic_overlap.bed.gz