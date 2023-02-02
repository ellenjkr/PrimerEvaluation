import os
import pandas as pd

from Bio import SeqIO
from primers_analysis_primer3 import PrimersAnalyzer


def add_barcodes2primers(primer_pair, barcodes):
	barcoded_primers_f = [barcode + primer_pair['forward_primer'] for barcode in barcodes]
	barcoded_primers_r = [barcode + primer_pair['reverse_primer'] for barcode in barcodes]

	return(barcoded_primers_f, barcoded_primers_r)


def create_primers_files(primer_pair, barcoded_primers_f, barcoded_primers_r, fasta_directory):
	primer_f_id = primer_pair['forward_id']
	primer_r_id = primer_pair['reverse_id']

	os.mkdir(f'{fasta_directory}/{primer_f_id}-{primer_r_id}')
	with open(f'{fasta_directory}/{primer_f_id}-{primer_r_id}/barcoded_{primer_f_id}_F.fasta', 'w') as f:
		for pos, primer in enumerate(barcoded_primers_f):
			f.write(f'>{pos}\n')
			f.write(f'{primer}\n')

	with open(f'{fasta_directory}/{primer_f_id}-{primer_r_id}/barcoded_{primer_r_id}_R.fasta', 'w') as f:
		for pos, primer in enumerate(barcoded_primers_r):
			f.write(f'>{pos}\n')
			f.write(f'{primer}\n')


def join_barcodes2primers(primer_pair, barcodes, files_directory, fasta_directory):
	barcoded_primers_f, barcoded_primers_r = add_barcodes2primers(primer_pair, barcodes)
	create_primers_files(primer_pair, barcoded_primers_f, barcoded_primers_r, fasta_directory)


def get_barcodes(directory, barcodes_file):
	if 'txt' in barcodes_file:
		with open(f"{directory}/{barcodes_file}", 'r') as f:
			barcodes = f.read().split('\n')
			if '' in barcodes:
				barcodes.remove('')
	elif barcodes_file.split('.')[1] in ['fa', 'fna', 'fasta', 'ffn', 'faa', 'frn']:
		barcodes = [str(record.seq) for record in SeqIO.parse(f"{directory}/{barcodes_file}", "fasta")]

	return barcodes


def get_primers(directory, primers_file):
	primers = pd.read_csv(f'{directory}/{primers_file}', sep=';')

	return primers


if __name__ == '__main__':
	ADD_BARCODES2PRIMERS = True
	PRIMERS_AND_BARCODES_DIR = 'Primers And Barcodes Files'
	BARCODES_FILE = 'barcodes_l5_n16000000_filtered.txt'
	PRIMERS_FILE = 'primer_dgHCO_dgLCO.csv'
	FASTA_FILES_DIR = 'Fasta Files - Primers'
	RESULTS_DIR = 'Results'
	GC_PERC = (35, 65)

	if ADD_BARCODES2PRIMERS:  # Read primers from fasta file, barcodes already added
		barcodes = get_barcodes(PRIMERS_AND_BARCODES_DIR, BARCODES_FILE)
		primers = get_primers(PRIMERS_AND_BARCODES_DIR, PRIMERS_FILE)
		for index, primer_pair in primers.iterrows():
			join_barcodes2primers(primer_pair, barcodes, PRIMERS_AND_BARCODES_DIR, FASTA_FILES_DIR)
			
	for folder in os.listdir(FASTA_FILES_DIR):
		files = os.listdir(f'{FASTA_FILES_DIR}/{folder}')
		files_path = [f'{FASTA_FILES_DIR}/{folder}/{file}' for file in files]
		pair_name = folder
		primers_analyzer = PrimersAnalyzer(files_path[0], files_path[1], pair_name, GC_PERC, RESULTS_DIR)
		primers_analyzer.run()
