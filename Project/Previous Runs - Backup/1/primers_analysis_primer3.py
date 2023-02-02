import os
import pandas as pd
import primer3

from Bio import SeqIO
from Bio.SeqUtils import GC


def get_primers_info(file):
	primers_dict = {
		'ID': [],
		'Seq': [],
		'Tm': [],
		'GC %': [],
	}

	for record in SeqIO.parse(file, "fasta"):
		primers_dict['ID'].append(record.id)

		seq = str(record.seq)
		
		primers_dict['Seq'].append(seq)
		primers_dict['Tm'].append(primer3.calcTm(seq))
		primers_dict['GC %'].append(round(GC(seq), 4))

		hairpin = primer3.calcHairpin(seq).todict()
		hairpin.pop('ascii_structure')
		primers_dict.update({f"Hairpin-{key}": [] for key in hairpin.keys() if f"Hairpin-{key}" not in primers_dict.keys()})
		for key in hairpin.keys():
			primers_dict[f'Hairpin-{key}'].append(hairpin[key])

		homodimer = primer3.calcHomodimer(seq).todict()
		homodimer.pop('ascii_structure')
		primers_dict.update({f"Homodimer-{key}": [] for key in homodimer.keys() if f"Homodimer-{key}" not in primers_dict.keys()})
		for key in homodimer.keys():
			primers_dict[f'Homodimer-{key}'].append(homodimer[key])


	return primers_dict


def build_excel_file(file, primers_dict):
	df = pd.DataFrame(primers_dict)

	writer = pd.ExcelWriter(file.replace('.fasta', '.xlsx'), engine='xlsxwriter')
	df.to_excel(writer, sheet_name='Sheet1', startrow=1, header=False, index=False)

	worksheet = writer.sheets['Sheet1']

	(max_row, max_col) = df.shape
	column_settings = [{'header': column} for column in df.columns]
	worksheet.add_table(0, 0, max_row, max_col - 1, {'columns': column_settings})
	worksheet.set_column(0, max_col - 1, 15)

	writer.save()


for file in os.listdir(os.getcwd()):
	if '.fasta' in file:
		primers_dict = get_primers_info(file)

		build_excel_file(file, primers_dict)
