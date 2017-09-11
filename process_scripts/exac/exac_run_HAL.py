import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import pickle
from helpful_utils import reverse_complement
import re
import scipy.io as sio
import scipy.stats
from math import log, exp, isnan
from copy import copy
import csv
from collections import Counter
from interval import interval
import warnings
warnings.filterwarnings('ignore')
import argparse

# Let's load in the Shendure model (Hexamer Additive Linear model, HAL), 
# determined an effect size for each exonic hexamer, as well as a model for the 
# exonic portion of the splice acceptor. The scoring functions and helper 
# functions are adapted from Rosenberg's github code used for the paper.

# DNA helper functions
bases = ['A','T','C','G']
def add_base(li):
    """Used in make_mer_list to add one more base to list"""
    new_li = []
    for s in li:
        for b in bases:
            new_li.append(s+b)
    return new_li
   

def make_mer_list(mer_len):
    """Makes a list of all n-mers"""
    li = bases
    for i in range(mer_len-1):
        li = add_base(li)
    return li


def logit(x):
    # refer to Wiki page on logit for asymptotic behavior at bounds
    if x == 0:
        return float('-Inf')
    elif x == 1:
        return float('Inf')
    else:
        return log(x) - log(1-x)


expit = lambda x: 1./(1.+exp(-x))


def extract_exon_seq(strand, intron1_len, exon_len, intron2_len, seq, extend=0, 
	rc=False):
    # extract exon sequence based on strand and intron/exon lengths. Grab n 
    # extra bp of downstream intron so hexamers at exon/intron border can be 
    # properly scored
    if not isinstance(seq, str):
        return ''
    
    seq = seq.upper()
    
    if strand not in ['+', '-', 1, -1, '1', '-1']:
        return ''
    if strand in ['-', -1, '-1']:
        if rc:
            seq = reverse_complement(seq)
        exon_seq = seq[intron2_len : intron2_len + exon_len + extend]
    else:
        exon_seq = seq[intron1_len : intron1_len + exon_len + extend]
    
    return exon_seq


def score_exon_seq(seq, mer_scores, exonic_acceptor_scores, sd_scores, 
	mult_factor=1):
    if seq == '':
        return float('NaN')
    if 'N' in seq:
        return float('NaN')
    # need extra 5bp of right intron to properly score hexamers at exon/intron 
    # junction
    score = 0.
    # score first 3 bp of exon as part of splice acceptor
    score += exonic_acceptor_scores.ix[seq[:3]]*mult_factor
    # score rest of exon up until donor with exon hexamer scores
    for b in range(len(seq)-5-6-3): # don't score last 3bp
        score += mer_scores[seq[b:b+6]]*mult_factor
    # Score the hexamers overlapping with the exonic portion of the splice donor 
    # at -3, -2, -1:
    for b in range(3):
        score += sd_scores.ix[seq[len(seq)-8+b:len(seq)-8+6+b],b]*mult_factor
    return score


def make_exon_skipping_predictions(exon_var_group, df, y_name, mult_factor=None):
    # y_name is column name of splicing index to use for predictions
    if mult_factor==None:
        mult_factor = np.ones(len(exon_var_group))*2.

    # check if there's a natural sequence, it may have been filtered out for 
    # low quality
    num_ref = ((df.ensembl_id == exon_var_group.ensembl_id.iloc[0]) & 
    	(df.sub_id == '000')).tolist().count(True)
    if num_ref == 0:
        exon_var_group['PSI_pred'] = float('NaN') * len(exon_var_group)
        exon_var_group['DPSI_pred'] = float('NaN') * len(exon_var_group)
        return exon_var_group
    # grab natural sequence for mutant
    else:
        ref = df[(df.ensembl_id == exon_var_group.ensembl_id.iloc[0]) & 
        (df.sub_id == '000')]
        # only 1 row
        ref_seq = ref.exon_seq.iloc[0]
        ref_score = score_exon_seq(ref_seq, exonic_mer6_scores, 
        	exonic_acceptor_scores, sd_scores)
    
    psi_pred_list = []
    
    for i in range(len(exon_var_group)):
        mut_seq = exon_var_group.exon_seq.iloc[i]
        mut_score = score_exon_seq(mut_seq, exonic_mer6_scores, 
        	exonic_acceptor_scores, sd_scores, mult_factor[i])
        psi_pred_list.append(expit(logit(ref[y_name].iloc[0]) + mut_score - ref_score))
    
    exon_var_group['PSI_pred'] = psi_pred_list
    exon_var_group['DPSI_pred'] = exon_var_group.PSI_pred - ref[y_name].iloc[0]
    
    return exon_var_group


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('infile', help='Name of input file')
	parser.add_argument('outfile', help='Name of output file')

	args = parser.parse_args()
	# read in data
	exac_data = pd.read_table(args.infile, sep = '\t')
	# reformat so intron/exon length columns are int and not float, convert NA to 0 so we can use int
	exac_data[['start', 'end', 'intron1_len', 'exon_len', 'intron2_len']] = \
		exac_data[['start', 'end', 'intron1_len', 'exon_len', 
		'intron2_len']].fillna(0.0).astype(int)


	# load splice site models
	exonic_mer6_scores = pd.read_pickle('../../ref/Rosenberg_2015/exonic_mer6_scores.series')
	exonic_acceptor_scores = pd.read_pickle('../../ref/Rosenberg_2015/exonic_acceptor_scores.series')
	model = sio.loadmat('../../ref/Rosenberg_2015/model_full_data.mat')

	# grab model information for positions overlapping the splice donor (-3, -2, -1, +1)
	sd_scores = pd.DataFrame(index=make_mer_list(6),data=model['Mer_scores'][:4**6*8].reshape(4**6,8)[:,2:6])

	exac_data['exon_seq'] = exac_data.apply(lambda x : 
		extract_exon_seq(x['strand'], x['intron1_len'], 
			x['exon_len'], x['intron2_len'], x['unflipped_seq'], 
			extend=5, rc=False), axis=1)

	# HAL only valid for exonic mutations, so let's subset our data and make 
	# our predictions.

	exac_exon_vars = exac_data[exac_data.label == 'exon']

	exac_exon_vars = exac_exon_vars.groupby('ensembl_id').apply(make_exon_skipping_predictions, 
		df=exac_data, y_name='v2_index')

	exac_exon_vars.to_csv(args.outfile, sep='\t', na_rep='NA', index=False, 
		quoting = csv.QUOTE_NONNUMERIC)
	