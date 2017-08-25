# Runs UCSC liftOver command to convert coordinates from hg19 to hg38
# called in exac_format_ref.R and depends on output file

# start and end coordinates
liftOver <(cat ../../ref/exac/exac_ref_formatted.txt | \
	awk '{OFS="\t"; print $4,$5,$6,$1}' | grep -v 'NA' | sort | uniq) \
  ../../ref/hg19ToHg38.over.chain \
  ../../processed_data/exac/exac_ref_liftover.bed \
  ../../processed_data/exac/exac_unlifted.bed

# SNP position
liftOver <(cat ../../ref/exac/exac_ref_formatted.txt | \
	awk '{OFS="\t"; print $4,$13,$13+1,$1}' | grep -v 'NA' | sort | uniq) \
  ../../ref/hg19ToHg38.over.chain \
  ../../processed_data/exac/exac_snp_liftover.bed \
  ../../processed_data/exac/exac_snp_unlifted.bed