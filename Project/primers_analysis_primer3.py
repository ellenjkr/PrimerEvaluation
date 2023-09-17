import os
import pandas as pd
import primer3

from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp as mt
from statistics import mean


class PrimersAnalyzer():
	def __init__(self, input_file_f, input_file_r, primer_pair, GC_PERC, RESULTS_DIR, REVERSE_COMPLEMENT, OUTPUT_FILE):
		super(PrimersAnalyzer, self).__init__()
		self.input_file_f = input_file_f
		self.input_file_r = input_file_r
		self.primer_pair = primer_pair
		self.GC_PERC = GC_PERC
		self.RESULTS_DIR = RESULTS_DIR
		self.REVERSE_COMPLEMENT = REVERSE_COMPLEMENT
		self.OUTPUT_FILE = OUTPUT_FILE


	def get_primers_from_fasta(self, file):
		primers_dict = {
			'ID': [],
			'Seq': [],
		}

		for record in SeqIO.parse(file, "fasta"):
			primers_dict['ID'].append(record.id)
			primers_dict['Seq'].append(str(record.seq))

		return primers_dict


	def get_tm_avg(self, seq):
		tm_methods = ['Tm Wallace', 'Tm NN', ]
		tm_methods.extend([f'Tm NN DNA_NN{i} Table' for i in range(1, 5)])
		tm_methods.extend([f'Tm GC valueset {i}' for i in range(1, 9)])

		results = []
		for method in tm_methods:
			if method == 'Tm Wallace':
				tm = round(mt.Tm_Wallace(seq, strict=False))
			elif method == 'Tm NN':
				tm = round(mt.Tm_NN(seq, strict=False), 2)
			elif 'DNA_NN' in method:
				nn_table = getattr(mt, method.split(' ')[2])
				tm = round(mt.Tm_NN(seq, nn_table=nn_table, strict=False), 2)
			elif 'GC' in method:
				valueset = int(method[-1])
				tm = round(mt.Tm_GC(seq, valueset=valueset, strict=False), 2)

			results.append(tm)

		tm = mean(results)

		return tm


	def get_primers_info(self, primers_dict, filter_gc=True):
		filtered_primers_dict = {
			'ID': [],
			'Seq': [],
			'Tm Primer3': [],
			'Tm Médio': [],
			'GC %': []
		}

		for pos, seq in enumerate(primers_dict['Seq']):
			
			gc_perc = round(GC(seq), 4)
			
			if gc_perc >= self.GC_PERC[0] and gc_perc <= self.GC_PERC[1] or filter_gc is False:
				filtered_primers_dict['GC %'].append(gc_perc)
				filtered_primers_dict['ID'].append(primers_dict['ID'][pos])
				filtered_primers_dict['Seq'].append(seq)

				filtered_primers_dict['Tm Primer3'].append(primer3.calcTm(seq))

				Tm = self.get_tm_avg(seq)
				# Tm = mt.Tm_GC(seq, valueset=2, strict=False)
				filtered_primers_dict['Tm Médio'].append(Tm)

				hairpin = primer3.calcHairpin(seq).todict()
				hairpin.pop('ascii_structure')
				filtered_primers_dict.update({f"Hairpin-{key}": [] for key in hairpin.keys() if f"Hairpin-{key}" not in filtered_primers_dict.keys()})
				for key in hairpin.keys():
					filtered_primers_dict[f'Hairpin-{key}'].append(hairpin[key])

				homodimer = primer3.calcHomodimer(seq).todict()
				homodimer.pop('ascii_structure')
				filtered_primers_dict.update({f"Homodimer-{key}": [] for key in homodimer.keys() if f"Homodimer-{key}" not in filtered_primers_dict.keys()})
				for key in homodimer.keys():
					filtered_primers_dict[f'Homodimer-{key}'].append(homodimer[key])

		return filtered_primers_dict

	def get_primer_pair_desc(self):
		primers_dict = {
			'ID': [],
			'Seq': [],
		}

		if self.REVERSE_COMPLEMENT:
			primers_dict['ID'].append(self.primer_pair['forward_id'] + ' Reverse Complement')
			primers_dict['ID'].append(self.primer_pair['reverse_id'] + ' Reverse Complement')
			primers_dict['Seq'].append(self.primer_pair['forward_primer_reverse_complement'])
			primers_dict['Seq'].append(self.primer_pair['reverse_primer_reverse_complement'])
		else:
			primers_dict['ID'].append(self.primer_pair['forward_id'])
			primers_dict['ID'].append(self.primer_pair['reverse_id'])
			primers_dict['Seq'].append(self.primer_pair['forward_primer'])
			primers_dict['Seq'].append(self.primer_pair['reverse_primer'])


		primers_dict = self.get_primers_info(primers_dict, filter_gc=False)

		df = pd.DataFrame(primers_dict)

		return df


	def build_primer_desc_sheet(self, writer):  # Build primer description sheet
		self.primer_pair.to_excel(writer, sheet_name='Primer Pair Description', startrow=0, header=False, index=True)
		worksheet = writer.sheets['Primer Pair Description']
		primer_pair_desc = self.get_primer_pair_desc()

		max_row = self.primer_pair.size

		primer_pair_desc.to_excel(writer, sheet_name='Primer Pair Description', startrow=max_row + 2, header=False, index=False)
		
		(max_row_table, max_col_table) = primer_pair_desc.shape
		column_settings = [{'header': column} for column in primer_pair_desc.columns]
		worksheet.add_table(max_row + 1, 0, max_row_table + max_row + 1, max_col_table - 1, {'columns': column_settings})
		worksheet.set_column(2, max_col_table - 1, 15)
		worksheet.set_column(0, 1, 40)

		return (writer, worksheet)

	def build_excel_sheet(self, writer, sheet_name, df):
		df.to_excel(writer, sheet_name=sheet_name, startrow=1, header=False, index=False)

		worksheet = writer.sheets[sheet_name]

		(max_row, max_col) = df.shape
		column_settings = [{'header': column} for column in df.columns]
		worksheet.add_table(0, 0, max_row, max_col - 1, {'columns': column_settings})
		worksheet.set_column(0, max_col - 1, 15)

		return (writer, worksheet)

	def build_tm_dg_chart(self, writer, worksheet, sheet_name, df, dg_column, dg_column_name, column_to_place_chart):
		workbook = writer.book
		chart = workbook.add_chart({'type': 'scatter'})
		(max_row, max_col) = df.shape

		chart.add_series({'name': f'Tm Médio X {dg_column_name}-DG', 'categories': f'={sheet_name}!${dg_column}$2:${dg_column}${max_row + 1}', 'values': f'={sheet_name}!$D$2:$D${max_row + 1}'})
		chart.set_title({'name': f'Tm Médio X {dg_column_name}-DG'})
		chart.set_x_axis({'name': f'{dg_column_name}-DG'})
		chart.set_y_axis({'name': 'Tm Médio'})
		chart.set_size({'width': 1000, 'height': 576})
		worksheet.insert_chart(f'{column_to_place_chart}{max_row + 3}', chart)

		return (writer, worksheet)

	def run(self):
		forward_id = self.primer_pair['forward_id']
		reverse_id = self.primer_pair['reverse_id']
		writer = pd.ExcelWriter(f'{self.RESULTS_DIR}/{self.OUTPUT_FILE}.xlsx', engine='xlsxwriter')
		self.build_primer_desc_sheet(writer)
		
		primers_dict_f = self.get_primers_from_fasta(self.input_file_f)
		primers_dict_f = self.get_primers_info(primers_dict_f)
		df = pd.DataFrame(primers_dict_f)
		writer, worksheet = self.build_excel_sheet(writer, 'Forward', df)
		writer, worksheet = self.build_tm_dg_chart(writer, worksheet, 'Forward', df, 'H', 'Hairpin', 'B')
		writer, worksheet = self.build_tm_dg_chart(writer, worksheet, 'Forward', df, 'M', 'Homodimer', 'L')

		primers_dict_r = self.get_primers_from_fasta(self.input_file_r)

		primers_dict_r = self.get_primers_info(primers_dict_r)
		df = pd.DataFrame(primers_dict_r)
		writer, worksheet = self.build_excel_sheet(writer, 'Reverse', df)
		writer, worksheet = self.build_tm_dg_chart(writer, worksheet, 'Reverse', df, 'H', 'Hairpin', 'B')
		writer, worksheet = self.build_tm_dg_chart(writer, worksheet, 'Reverse', df, 'M', 'Homodimer', 'L')
		writer.close()
