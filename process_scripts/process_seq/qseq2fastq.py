'''
This script converts qseq files to FASTQ format. 

Input: text file of filenames, one per line, of files to be concatenated and
converted.
Output: one converted FASTQ filter
'''

import sys
import argparse
import gzip

if __name__ == '__main__':

	parser = argparse.ArgumentParser('converts qseq files to FASTQ format')
	parser.add_argument('qseqs', help='.txt file of filenames, one per line,\
		of files to be concatenated together and converted')
	parser.add_argument('output', help='name of output file')

	args = parser.parse_args()

	outfile = open(args.output, 'w')

	with open(args.qseqs) as infile:
		for filename in infile:
			filename = filename.strip()
			if filename.endswith('.gz'):
				qseq_file = gzip.open(filename, 'rb')
			else:
				qseq_file = open(filename)
			for line in qseq_file:
				fields = line.strip().split('\t')
				# remove reads that don't pass filter which == 0
				if fields[10] == '1': 
					header = '@' + ':'.join(fields[0:6]) + '#' + \
							fields[6] + '/' + fields[7]

					# replace uncalled bases with N
					seq = fields[8].replace('.', 'N')

					outfile.write('\n'.join([header, seq, '+', fields[9], '']))

