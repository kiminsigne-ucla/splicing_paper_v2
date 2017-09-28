'''
@author: Kimberly Insigne
kiminsigne@gmail.com

This script demultiplexes two FASTQ file corresponding to paired reads into 
separate FASTQ files based on their index. If an index read is 1bp away from an 
index in the given reference file, it will be included in that sample FASTQ file.

Input:
reads_file_1	: FASTQ file for read 1
reads_file_2	: FASTQ file for read 2
index_reads_file  : FASTQ file of index read, where index is first N bp of read
index_file	: tab-separated text file, no headers, first column is sample name, 
second column is index
index_length	: integer, length of index 

Optional:
output	: name of output directory, default is current directory
rev		: if enabled, sequences in index file are reverse complement 
relative to index reads file

Output:
Separate FASTQ files, 2 for each sample, one for each paired read
'''

import argparse
from itertools import islice
from helpful_utils import reverse_complement, mutate_1_bp


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('read1_file', help='.fastq file of R1')
	parser.add_argument('read2_file', help='.fastq file of R2')
	parser.add_argument('index_reads_file', 
		help='.fastq file, where the index is the first 6bp of the read')
	parser.add_argument('index_file', 
		help='tab separated text file where first column is sample name \
												and second column is index')
	parser.add_argument('index_length', type=int, 
		help="length of index, assuming it's the first n bases")
	# parser.add_argument('output', default='./',
	# 	help='Name of output directory, default is current directory')
	parser.add_argument('-rev', action='store_true', 
		help='If sequences in index_file are reverse complement relative to \
		index reads file''')

	args = parser.parse_args()		

	# read in index file
	index_file = open(args.index_file, 'rU')
	indices = {}
	for line in index_file:
		name, index = line.strip().split('\t')
		if args.rev:
			index = reverse_complement(index)
		indices[index] = name
		# generate all 1bp mutants and add to look-up so we don't have to 
		# check each index
		for x in mutate_1_bp(index):
			indices[x] = name

	# because we added 1bp mutant indices there will be lots of repeated 
	# sample entries in the dictionary, so make this a set
	samples = set(indices.values())
	index_to_handle = {}
	for sample in samples:
		index_to_handle[sample+'_read1'] = open(sample+'_read1.fastq', 'w')
		index_to_handle[sample+'_read2'] = open(sample+'_read2.fastq', 'w')

	
	bad_index_1 = open('bad_index_1.fastq', 'w')
	bad_index_3 = open('bad_index_3.fastq', 'w')

	read1_file = open(args.read1_file)
	read2_file = open(args.read2_file)
	index_reads_file = open(args.index_reads_file)

	count = 0

	while True:
		read1_info = list(islice(read1_file, 4))
		read2_info = list(islice(read2_file, 4))
		index_info = list(islice(index_reads_file, 4))
		if not read1_info:
			break

		if count % 10000000 == 0:
			print count, '...'

		index = index_info[1].strip()[:args.index_length]

		if index in indices:
			filehandle1 = index_to_handle[indices[index] + '_read1']
			filehandle2 = index_to_handle[indices[index] + '_read2']
		else:
		 	filehandle1 = bad_index_1
			filehandle2 = bad_index_3
			
		filehandle1.writelines(read1_info)
		filehandle2.writelines(read2_info)
		count += 1

	# close everything
	for x in index_to_handle:
		index_to_handle[x].close()

	read1_file.close()
	read2_file.close()
	index_reads_file.close()
	bad_index_1.close()
	bad_index_3.close()


