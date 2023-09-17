import os
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from primers_analysis_primer3 import PrimersAnalyzer


def add_barcodes2primers(primer_pair, forward_barcodes, reverse_barcodes, REVERSE_COMPLEMENT):
	if REVERSE_COMPLEMENT:
		barcoded_primers_f = [barcode + primer_pair['forward_primer_reverse_complement'] for barcode in forward_barcodes]
		barcoded_primers_r = [barcode + primer_pair['reverse_primer_reverse_complement'] for barcode in reverse_barcodes]
	else:
		barcoded_primers_f = [barcode + primer_pair['forward_primer'] for barcode in forward_barcodes]
		barcoded_primers_r = [barcode + primer_pair['reverse_primer'] for barcode in reverse_barcodes]

	return(barcoded_primers_f, barcoded_primers_r)


def create_primers_files(primer_pair, barcoded_primers_f, barcoded_primers_r, fasta_directory):
	primer_f_id = primer_pair['forward_id']
	primer_r_id = primer_pair['reverse_id']

	if os.path.exists(f'{fasta_directory}/{primer_f_id}-{primer_r_id}') is False:
		os.mkdir(f'{fasta_directory}/{primer_f_id}-{primer_r_id}')
	with open(f'{fasta_directory}/{primer_f_id}-{primer_r_id}/barcoded_{primer_f_id}_F.fasta', 'w') as f:
		for pos, primer in enumerate(barcoded_primers_f):
			f.write(f'>{pos}\n')
			f.write(f'{primer}\n')

	with open(f'{fasta_directory}/{primer_f_id}-{primer_r_id}/barcoded_{primer_r_id}_R.fasta', 'w') as f:
		for pos, primer in enumerate(barcoded_primers_r):
			f.write(f'>{pos}\n')
			f.write(f'{primer}\n')


def join_barcodes2primers(primer_pair, forward_barcodes, reverse_barcodes, files_directory, fasta_directory, REVERSE_COMPLEMENT):
	barcoded_primers_f, barcoded_primers_r = add_barcodes2primers(primer_pair, forward_barcodes, reverse_barcodes, REVERSE_COMPLEMENT)
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


def discard_baseruns(barcodes, primer_pair, REVERSE_COMPLEMENT):
	if REVERSE_COMPLEMENT:
		forward_primer = primer_pair['forward_primer_reverse_complement']
		reverse_primer = primer_pair['reverse_primer_reverse_complement']
	else:
		forward_primer = primer_pair['forward_primer']
		reverse_primer = primer_pair['reverse_primer']

	forward_filtered_barcodes = []
	reverse_filtered_barcodes = []
	for barcode in barcodes:

		if barcode[-1] != forward_primer[0]:
			forward_filtered_barcodes.append(barcode)
		if barcode[-1] != reverse_primer[0]:
			reverse_filtered_barcodes.append(barcode)


	return forward_filtered_barcodes, reverse_filtered_barcodes


def get_primers(directory, primers_file, REVERSE_COMPLEMENT):
	primers = pd.read_csv(f'{directory}/{primers_file}', sep=';')

	if REVERSE_COMPLEMENT:
		forward_reverse_complement = primers['forward_primer'].apply(lambda x: str(Seq(x).reverse_complement()))
		primers.insert (2, 'forward_primer_reverse_complement', forward_reverse_complement)
		reverse_reverse_complement = primers['reverse_primer'].apply(lambda x: str(Seq(x).reverse_complement()))
		primers.insert (5, 'reverse_primer_reverse_complement', reverse_reverse_complement)


	return primers


if __name__ == '__main__':
	ADD_BARCODES2PRIMERS = True
	REVERSE_COMPLEMENT = False
	PRIMERS_AND_BARCODES_DIR = 'Primers And Barcodes Files'
	BARCODES_FILE = 'barcodes_l5_n16000000_nobaseruns_filtered.fa'
	PRIMERS_FILE = 'primers.csv'
	FASTA_FILES_DIR = 'Fasta Files - Primers'
	OUTPUT_FILE_PATTERN = 'barcoded_l5_primer'
	RESULTS_DIR = 'Results'
	DISCARD_BASERUNS = True  # Discard baseruns between barcodes and primers
	GC_PERC = (35, 65)


	primers = get_primers(PRIMERS_AND_BARCODES_DIR, PRIMERS_FILE, REVERSE_COMPLEMENT)
	if ADD_BARCODES2PRIMERS:  # Read primers from fasta file, barcodes already added
		barcodes = get_barcodes(PRIMERS_AND_BARCODES_DIR, BARCODES_FILE)
	
	for index, primer_pair in primers.iterrows():
		forward_primer = primer_pair['forward_id']
		reverse_primer = primer_pair['reverse_id']
		folder = f'{FASTA_FILES_DIR}/{forward_primer}-{reverse_primer}'
		if ADD_BARCODES2PRIMERS:
			if DISCARD_BASERUNS:
				forward_filtered_barcodes, reverse_filtered_barcodes = discard_baseruns(barcodes, primer_pair, REVERSE_COMPLEMENT)
				join_barcodes2primers(primer_pair, forward_filtered_barcodes, reverse_filtered_barcodes, PRIMERS_AND_BARCODES_DIR, FASTA_FILES_DIR, REVERSE_COMPLEMENT)
			else:
				join_barcodes2primers(primer_pair, barcodes, barcodes, PRIMERS_AND_BARCODES_DIR, FASTA_FILES_DIR, REVERSE_COMPLEMENT)

		files = os.listdir(folder)
		for file in files:

			if file.split('.')[0][-1] == 'F':
				f_file = f'{folder}/{file}'
			if file.split('.')[0][-1] == 'R':
				r_file = f'{folder}/{file}' 


		OUTPUT_FILE = OUTPUT_FILE_PATTERN.replace('primer', f'{forward_primer}-{reverse_primer}')
		primers_analyzer = PrimersAnalyzer(f_file, r_file, primer_pair, GC_PERC, RESULTS_DIR, REVERSE_COMPLEMENT, OUTPUT_FILE)
		primers_analyzer.run()
