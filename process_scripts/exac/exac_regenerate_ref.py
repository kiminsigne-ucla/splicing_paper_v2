import pandas as pd
from cogent.db.ensembl import HostAccount, Genome
from numpy import isnan


# def pycog_grab_seq(chrom, start, end, strand, genome):
# 	# remove 'chr' if exists, convert to int
# 	chrom = chrom.replace('chr', '')
# 	# convert to int, may be double
# 	start = int(start)
# 	end = int(end)

# 	# assume start and end are 1-based, pycog is 0-based
# 	cog_reg = genome.getRegion(CoordName=chrom, Start=start - 1, End=end)
# 	cog_seq = cog_reg.Seq

# 	if strand == '-':
# 		cog_seq = cog_seq.rc()
# 	seq = cog_seq._seq

# 	return seq

def pycog_grab_seq_by_group(df_group, genome):
	# for use with pandas groupby
	# whole group has some chromosome, start, end, strand
	# remove 'chr' if exists, convert to int
	chrom = df_group.chr.iloc[0].replace('chr', '')
	# convert to int, may be double
	start = int(df_group.start.iloc[0])
	end = int(df_group.end.iloc[0])
	strand = df_group.strand.iloc[0]

	# assume start and end are 1-based, pycog is 0-based
	cog_reg = genome.getRegion(CoordName=chrom, Start=start - 1, End=end)
	cog_seq = cog_reg.Seq

	if strand == '-':
		cog_seq = cog_seq.rc()
	seq = cog_seq._seq

	return seq


def sub_alt_allele(seq, alt_allele, start, end, snp_position, strand):
	# assume sequence is strand-aware (negative strand already reverse complemented)
	# get relative position of SNV
	if isnan(snp_position):
		return seq
	if len(alt_allele) != 1: # controls have more than one
		return seq
	seq_list = list(seq)
	if strand == '+':
		rel_position = snp_position - start
	elif strand == '-':
		rel_position = end - snp_position

	seq_list[rel_position] = alt_allele
	return ''.join(seq_list)	


# use pycogent to get sequences, set up the connection
pycog = HostAccount('ensembldb.ensembl.org', # host
                    'anonymous', # user
                    '', # password
                    3306) # port

hs37 = Genome('human', Release=78, account=pycog)

# read in formatted reference file
ref = pd.read_table('../../ref/exac/exac_ref_formatted_converted.txt',
	sep='\t', header=0)
# convert to int
ref[['start', 'end', 'snp_position']] = ref[['start', 'end', 'snp_position']].fillna(0.0).astype(int)
# fill NaN in alt_allele with blank string
ref['alt_allele'] = ref.alt_allele.fillna('')

# ref['nat_seq'] = ref.apply(lambda x: pycog_grab_seq(x['chr'], x['start'],
# 	x['end'], x['strand'], hs37), axis=1)
nat_seqs = ref.groupby('ensembl_id').apply(pycog_grab_seq_by_group, genome=hs37)
# returns named Series, convert to data frame
nat_seqs = pd.DataFrame({'ensembl_id' : nat_seqs.index, 
	'nat_seq' : nat_seqs.values})
# merge to ref
ref = ref.merge(nat_seqs, on='ensembl_id', how='left')

ref['original_seq'] = ref.apply(lambda x: sub_alt_allele(x['nat_seq'],
	x['alt_allele'], x['start'], x['end'], x['snp_position'], x['strand']), axis=1)

ref.to_csv('../../ref/exac/exac_ref_formatted_converted_regen_seq.txt',
	sep='\t', index=False)



