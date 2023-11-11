import pandas as pd
import numpy as np
import itertools
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.SeqUtils import IUPACData
import primer3


def merge_additional_info(all_data, forward_data, reverse_data):
    all_data = all_data.merge(forward_data, left_on="Forward ID", right_on="ID", how="left").drop(columns=["ID", "Seq"])
    all_data = all_data.merge(reverse_data, left_on="Reverse ID", right_on="ID", how="left").drop(columns=["ID", "Seq"])
    all_data = all_data.rename(columns={"Forward Tm Primer_x": "Forward Tm Primer", "Forward Homodimer Tm_x": "Forward Homodimer Tm", "Forward Hairpin Tm_x": "Forward Hairpin Tm", "Reverse Tm Primer_y": "Reverse Tm Primer", "Reverse Homodimer Tm_y": "Reverse Homodimer Tm", "Reverse Hairpin Tm_y": "Reverse Hairpin Tm"})
    return all_data


def degenerate_seq_variants(seq):
    variants = ['']
    for letter in seq:
        if letter in IUPACData.ambiguous_dna_values:
            new_variants = []
            for base in IUPACData.ambiguous_dna_values[letter]:
                for variant in variants:
                    new_variants.append(variant + base)
            variants = new_variants
        else:
            variants = [variant + letter for variant in variants]
    return variants


def calculate_heterodimer_properties(Forward_sequence, Reverse_sequence):
    # obtem variantes para os primers degenerados
    forward_variants = degenerate_seq_variants(Forward_sequence)
    reverse_variants = degenerate_seq_variants(Reverse_sequence)

    # calcula as propriedades para cada combinação degenerada do par
    dg_values = []
    tm_values = []
    for f_seq in forward_variants:
        for r_seq in reverse_variants:
            heterodimer = primer3.calcHeterodimer(f_seq, r_seq)
            dg_values.append(heterodimer.dg)
            tm_values.append(heterodimer.tm)

    # faz a media entre as combinações degeneradas
    heterodimer_dg = sum(dg_values) / len(dg_values)
    heterodimer_tm = sum(tm_values) / len(tm_values)
    return heterodimer_dg, heterodimer_tm


def analyze_primers(data, direction):
    tm_values = []
    homodimer_tm_values = []
    hairpin_tm_values = []

    for _, row in data.iterrows():
        sequence = row["Seq"]

        tm = mt.Tm_NN(Seq(sequence))
        hairpin = primer3.calcHairpin(sequence)
        homodimer = primer3.calcHomodimer(sequence)

        tm_values.append(tm)
        homodimer_tm_values.append(homodimer.tm)
        hairpin_tm_values.append(hairpin.tm)

    data[f"{direction} Tm Primer"] = tm_values
    data[f"{direction} Homodimer Tm"] = homodimer_tm_values
    data[f"{direction} Hairpin Tm"] = hairpin_tm_values

    return data


def analyze_excel(file_path):
    df = pd.read_excel(file_path, sheet_name=["Forward", "Reverse"])
    forward_data = df["Forward"][["ID", "Seq"]]  # primer-barcode id e primer-barcode sequence
    reverse_data = df["Reverse"][["ID", "Seq"]]

    forward_data = analyze_primers(forward_data, "Forward")  # calcula as propriedades fisico quimicas de cada combinação primer-barcode
    reverse_data = analyze_primers(reverse_data, "Reverse")

    combinations = list(itertools.product(forward_data.iterrows(), reverse_data.iterrows()))  # combinações possíveis entre primers-barcodes forward e reverse
    results = []
    for (f_index, f_row), (r_index, r_row) in combinations:
        heterodimer_dg, heterodimer_tm = calculate_heterodimer_properties(f_row["Seq"], r_row["Seq"])  # analisa as propriedades de heterodimer do par
        results.append({"Forward ID": f_row["ID"], "Reverse ID": r_row["ID"], "Heterodimer-dG": heterodimer_dg, "Heterodimer-Tm": heterodimer_tm})


    pair_results_df = pd.DataFrame(results)
    all_data = pair_results_df.copy()
    all_data.to_csv('temp.tsv', sep='\t', index=False)
    all_data = merge_additional_info(all_data, forward_data, reverse_data)

    return forward_data, reverse_data, pair_results_df, all_data


def part2(all_data):
    # Adicionar uma nova coluna 'Pair-ID' e numerar as linhas que contêm informações
    all_data.insert(0, 'Pair-ID', np.arange(1, len(all_data) + 1))

    # Calcular o quantil inferior de 50% para a coluna 'Heterodimer-Tm'
    q50_inf = all_data['Heterodimer-Tm'].quantile(0.5)

    # Selecionar os registros que ocorrem na faixa de 50% quantil inferior
    selected_records = all_data[all_data['Heterodimer-Tm'] <= q50_inf]

    return all_data, selected_records


def part3(inf_het_quantil):
    # Inicializar variáveis
    interval = 2.5
    step = 0.1
    results = []

    # Obter os valores mínimo e máximo dos conjuntos de dados
    min_value = min(inf_het_quantil['Forward Tm Primer'].min(), inf_het_quantil['Reverse Tm Primer'].min())
    max_value = max(inf_het_quantil['Forward Tm Primer'].max(), inf_het_quantil['Reverse Tm Primer'].max())

    # Avaliar o intervalo para a maior co-ocorrência de valores
    while (min_value + interval) <= max_value:
        forward_count = inf_het_quantil[(inf_het_quantil['Forward Tm Primer'] >= min_value) & (inf_het_quantil['Forward Tm Primer'] <= min_value + interval)].count()[0]
        reverse_count = inf_het_quantil[(inf_het_quantil['Reverse Tm Primer'] >= min_value) & (inf_het_quantil['Reverse Tm Primer'] <= min_value + interval)].count()[0]
        co_occurrence = forward_count + reverse_count
        results.append({'co-ocorr_min': min_value, 'co-ocorr_max': min_value + interval, 'co-ocorrencia-contag': co_occurrence})
        min_value += step

    # Criar DataFrame com os resultados
    results_df = pd.DataFrame(results)

    # Selecionar os registros de "Pair-ID" que estão dentro do intervalo de co-ocorrência
    max_co_occurrence = results_df['co-ocorrencia-contag'].max()
    max_interval = results_df[results_df['co-ocorrencia-contag'] == max_co_occurrence].iloc[0]
    coocurring_tm_selected = inf_het_quantil[((inf_het_quantil['Forward Tm Primer'] >= max_interval['co-ocorr_min']) & (inf_het_quantil['Forward Tm Primer'] <= max_interval['co-ocorr_max'])) & ((inf_het_quantil['Reverse Tm Primer'] >= max_interval['co-ocorr_min']) & (inf_het_quantil['Reverse Tm Primer'] <= max_interval['co-ocorr_max']))]

    return results_df, coocurring_tm_selected


def part4(coocurring_tm_selected):
    # Determinar o quantil de 50% inferior para as variáveis especificadas
    q50 = {
        'Forward Homodimer Tm': coocurring_tm_selected['Forward Homodimer Tm'].quantile(0.5),
        'Forward Hairpin Tm': coocurring_tm_selected['Forward Hairpin Tm'].quantile(0.5),
        'Reverse Homodimer Tm': coocurring_tm_selected['Reverse Homodimer Tm'].quantile(0.5),
        'Reverse Hairpin Tm': coocurring_tm_selected['Reverse Hairpin Tm'].quantile(0.5),
    }
    # Selecionar os "Pair-ID" que possuam valores dentro do quantil de 50% inferior para todas as variáveis simultaneamente

    inf_homohair_quantil = coocurring_tm_selected[
        (coocurring_tm_selected['Forward Homodimer Tm'] <= q50['Forward Homodimer Tm']) &
        (coocurring_tm_selected['Forward Hairpin Tm'] <= q50['Forward Hairpin Tm']) &
        (coocurring_tm_selected['Reverse Homodimer Tm'] <= q50['Reverse Homodimer Tm']) &
        (coocurring_tm_selected['Reverse Hairpin Tm'] <= q50['Reverse Hairpin Tm'])
    ]    
    
    return inf_homohair_quantil


def part5(inf_homohair_quantil):
    # Criar a lista ranqueada e contagem de ocorrências para "Forward ID"
    forward_id_count = inf_homohair_quantil['Forward ID'].value_counts().reset_index()
    forward_id_count.columns = ['Forward ID', 'Count']

    # Criar a lista ranqueada e contagem de ocorrências para "Reverse ID"
    reverse_id_count = inf_homohair_quantil['Reverse ID'].value_counts().reset_index()
    reverse_id_count.columns = ['Reverse ID', 'Count']

    # Adicionar a coluna "List-Reverse ID" na forward_id_count
    forward_id_count['List-Reverse ID'] = forward_id_count['Forward ID'].apply(lambda x: inf_homohair_quantil[inf_homohair_quantil['Forward ID'] == x]['Reverse ID'].unique())

    # Adicionar a coluna "List-Forward ID" na reverse_id_count
    reverse_id_count['List-Forward ID'] = reverse_id_count['Reverse ID'].apply(lambda x: inf_homohair_quantil[inf_homohair_quantil['Reverse ID'] == x]['Forward ID'].unique())

    return forward_id_count, reverse_id_count


def part6(forward_id_count, reverse_id_count, inf_homohair_quantil):
    # Selecionar os primeiros 20 registros de "Forward ID" e "Reverse ID"
    top_forward = forward_id_count.head(20)
    top_reverse = reverse_id_count.head(20)

    # Determinar todas as combinações possíveis entre um registro "Forward ID" e um registro "Reverse ID"
    combinations = list(itertools.product(top_forward['Forward ID'], top_reverse['Reverse ID']))

    # Criar um DataFrame com as combinações possíveis
    combina_possib_df = pd.DataFrame(combinations, columns=['Forward ID', 'Reverse ID'])

    # Selecionar os registros de "Pair-ID" que forem idênticos quanto ao "Forward ID" e "Reverse ID" das combinações consideradas
    selected_pairs = inf_homohair_quantil[inf_homohair_quantil.apply(lambda row: (row['Forward ID'], row['Reverse ID']) in combinations, axis=1)]

    selected_pairs = selected_pairs.reset_index(drop=True)
    return combina_possib_df, selected_pairs


if __name__ == '__main__':
    
    INPUT_FILE = 'Results/barcoded_l5_dgHCO-2198-dgLCO-1490.xlsx'

    # =======================================================================================
    # PARTE 1 - FAZ CÁLCULOS
    # forward_data, reverse_data, pair_results_df, all_data = analyze_excel(INPUT_FILE)
    forward_data = pd.read_excel('result_1.xlsx', engine='openpyxl', sheet_name='Forward')
    reverse_data = pd.read_excel('result_1.xlsx', engine='openpyxl', sheet_name='Reverse')
    pair_results_df = pd.read_excel('result_1.xlsx', engine='openpyxl', sheet_name='Results')
    all_data = pd.read_excel('result_1.xlsx', engine='openpyxl', sheet_name='All')
    # =======================================================================================
    # PARTE 2 - FAZER A SELEÇÃO - QUANTIL INFERIOR DE HETERODIMERO
    all_data, inf_het_quantil = part2(all_data)

    # PARTE 3 - FAZER A SELEÇÃO PARA TM CO-OCORRENTE
    interval_count, coocurring_tm_selected = part3(inf_het_quantil)

    # PARTE 4 - FAZER A SELEÇÃO - QUANTIL INFERIOR HOMODIMERO E HAIRPIN
    inf_homohair_quantil = part4(coocurring_tm_selected)

    # PARTE 5 - FAZER RANK E SELEÇÃO DOS PRIMERS MAIS OCORRENTES
    forward_id_count, reverse_id_count = part5(inf_homohair_quantil)

    # PARTE 6 - SELEÇÃO DE PRIMERS MAIS OCORRENTES, COMBINAÇÕES E SELEÇÃO
    combina_possib_df, selected_pairs = part6(forward_id_count, reverse_id_count, inf_homohair_quantil)

    # PARTE 7 - RECUPERAR AS INFORMAÇÕES DOS PRIMERS SELECIONADOS
    selected_pairs = selected_pairs.merge(forward_data[['ID', 'Seq']], left_on='Forward ID', right_on='ID', how='inner')
    selected_pairs = selected_pairs.drop(columns='ID')
    selected_pairs = selected_pairs.rename(columns={'Seq': 'Seq Forward'})

    selected_pairs = selected_pairs.merge(reverse_data[['ID', 'Seq']], left_on='Reverse ID', right_on='ID', how='inner')
    selected_pairs = selected_pairs.drop(columns='ID')
    selected_pairs = selected_pairs.rename(columns={'Seq': 'Seq Reverse'})

    # PARTE 8 - CRIAR ARQUIVO EXCEL

    output_path = "result_ALL.xlsx"
    
    with pd.ExcelWriter(output_path) as writer:
        forward_data.to_excel(writer, sheet_name="Forward", index=False)
        reverse_data.to_excel(writer, sheet_name="Reverse", index=False)
        pair_results_df.to_excel(writer, sheet_name="Results", index=False)
        all_data.to_excel(writer, sheet_name="All", index=False)
        all_data.to_excel(writer, sheet_name="All", index=False)
        selected_pairs.to_excel(writer, sheet_name="Selected Pairs", index=False)
