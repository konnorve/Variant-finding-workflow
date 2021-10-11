from pathlib import Path
import pandas as pd
import numpy as np

import gffpandas.gffpandas as gffpd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def convert_list_column_to_columns(df, col, output_names, deliminator):
    col_list = df[col].to_numpy()
    col_arr = np.empty((len(col_list), len(output_names)), dtype=object)
    r = range(len(output_names))
    for i in range(len(col_list)):
        s = col_list[i]
        if s is None:
            elements = [None] * len(output_names)
        else:
            elements = s.split(deliminator)
        for j in r:
            col_arr[i][j] = elements[j]
    return col_arr

def convert_dict_column_to_columns(df, col, output_names, element_deliminator, key_value_deliminator, DEBUG=False):
    col_list = df[col].to_numpy()
    col_arr = np.empty((len(col_list), len(output_names)), dtype=object)
    for i in range(len(col_list)):
        row_dict = dict(ele.split("=") for ele in col_list[i].split(element_deliminator) if len(ele.split(key_value_deliminator))==2)
        for j, col in enumerate(output_names):
            col_arr[i][j] = row_dict.get(col)
        if DEBUG and i % 10000 == 0: print(i, "\t", row_dict, "\t", col_arr[i])
    return col_arr


def get_unique_INFO_elements(position_call_df, DEBUG=False):
    info_list = position_call_df['INFO'].to_numpy()
    
    info_columns = []
    for i in range(info_list.size):
        s = info_list[i]
        if DEBUG and i % 10000 == 0: print(i, "\t", info_columns)
        for ele in s.split(';'):
            col = ele[:ele.find('=')]
            if col not in info_columns:
                info_columns.append(col)

    return info_columns


def get_parsed_position_call_df(position_call_df, info_columns=None, DEBUG=False):

    if info_columns is None:
        info_columns = get_unique_INFO_elements(position_call_df)

    if DEBUG: print('current columns:', position_call_df.columns)
    if DEBUG: print('info columns:', info_columns)

    info_arr = convert_dict_column_to_columns(position_call_df, 'INFO', info_columns, ';', '=', DEBUG=DEBUG)
    position_call_df[info_columns] = info_arr

    # must convert DP4 and PV4 to columns
    dp4_columns = ['ref_forward_reads', 'ref_reverse_reads', 'alt_forward_reads', 'alt_reverse_reads']
    dp4_arr = convert_list_column_to_columns(position_call_df, 'DP4', dp4_columns, ',')
    position_call_df[dp4_columns] = dp4_arr

    pv4_columns = ['strand_bias_pVal', 'baseQ_bias_pVal', 'mapQ_bias_pVal', 'tail_dist_bias_pVal']
    pv4_arr = convert_list_column_to_columns(position_call_df, 'PV4', pv4_columns, ',')
    position_call_df[pv4_columns] = pv4_arr

    position_call_parsed_df = change_datatypes(position_call_df)

    position_call_parsed_df['fwd_strand_coverage'] = position_call_parsed_df['ref_forward_reads'] + position_call_parsed_df['alt_forward_reads']
    position_call_parsed_df['rev_strand_coverage'] = position_call_parsed_df['ref_reverse_reads'] + position_call_parsed_df['alt_reverse_reads']

    if DEBUG: print(position_call_df.head())
    if DEBUG: print(position_call_df.columns)

    print(position_call_parsed_df.dtypes)
    for col in position_call_parsed_df.columns:
        print(col, '\t', position_call_parsed_df[col].unique()[:5])

    return position_call_parsed_df


def change_datatypes(df):
    all_columns_datatype_dict = {'POS':'uint64', 'ID':'str_', 'REF':'str_', 'ALT':'str_', 'QUAL':"float64", 
                            'FILTER':"str_", 'INFO':"str_", 'FORMAT':"str_", 'DP':"uint64", 'FS':"float64", 
                            'MQ0F':"float64", 'AF1':"float64", 'AC1':"float64", 'DP4':"str_", 'MQ':"float64", 
                            'FQ':"float64", 'SGB':"float64", 'RPBZ':"float64", 'BQBZ':"float64", 'PV4':"str_",
                            'SCBZ':"float64", 'VDB':"float64", 'INDE':"str_", 'IDV':"float64", 'IMF':"float64",
                            'ref_forward_reads':'uint64', 'ref_reverse_reads':'uint64', 'alt_forward_reads':'uint64', 'alt_reverse_reads':'uint64',
                            'strand_bias_pVal':"float64", 'baseQ_bias_pVal':"float64", 'mapQ_bias_pVal':"float64", 'tail_dist_bias_pVal':"float64",
                            'rev_strand_coverage':'uint64', 'fwd_strand_coverage':'uint64'}

    column_datatype_dict = {key:value for key, value in zip(all_columns_datatype_dict.keys(), all_columns_datatype_dict.values()) if key in list(df.columns)}
    return df.astype(column_datatype_dict)


def get_position_call_df(all_calls_path):
    raw_df = pd.read_csv(all_calls_path, sep='\t')
    narrow_df = raw_df[['POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']]
    return narrow_df


def get_gff_df(gff_path):
    annotation = gffpd.read_gff3(gff_path)
    attributes_df = annotation.attributes_to_columns()
    return attributes_df


def get_position_call_df_all_samples(all_calls_dir):
    dfs2concat = []
    for sample_path in all_calls_dir.iterdir():
        sample_name = sample_path.stem.split('_')[0]
        sample_df = get_position_call_df(sample_path)
        col_order = list(sample_df.columns)
        sample_df['genome_name'] = sample_name
        col_order.insert(0, 'genome_name')
        sample_df = sample_df[col_order]
        print(sample_df.columns)
        dfs2concat.append(sample_df)
    return pd.concat(dfs2concat, axis=0, ignore_index=True)


def filterVariants(all_position_calls_df):
    return all_position_calls_df[all_position_calls_df['ALT'] != '.']


def filterOutNonSigVariants(variant_df):
    # each strand must have at least 2x coverage
    filtered_variant_df = variant_df[variant_df['DP'] > 5]
    filtered_variant_df = filtered_variant_df[filtered_variant_df['FQ'] >= 30]
    filtered_variant_df = filtered_variant_df[(filtered_variant_df['fwd_strand_coverage'] >= 2) & (filtered_variant_df['rev_strand_coverage'] >= 2)]
    return filtered_variant_df


def indicate_nearby_genes(filtered_variant_df, gff_df):

    nearby_genes_df_list = []

    for i, row in filtered_variant_df.iterrows():
        loc = filtered_variant_df.loc[i, 'POS']
        nearby_genes_df = gff_df[((gff_df['end'] > loc - 1000) & (gff_df['end'] < loc + 1000)) | ((gff_df['start'] > loc - 1000) & (gff_df['start'] < loc + 1000))]
        nearby_genes_df['genome_name'] = filtered_variant_df.loc[i, 'genome_name']
        nearby_genes_df['POS'] = loc
        nearby_genes_df['REF'] = filtered_variant_df.loc[i, 'REF']
        nearby_genes_df['ALT'] = filtered_variant_df.loc[i, 'ALT']

        nearby_genes_df_list.append(nearby_genes_df)

    nearby_genes_df = pd.concat(nearby_genes_df_list, axis=0)

    filtered_variant_df_w_genes =  filtered_variant_df.merge(nearby_genes_df)

    return filtered_variant_df_w_genes[['genome_name','POS','REF','ALT','QUAL','start','end','strand','Name','locus_tag','product','translation','transl_table']]


def saveVariantDF(df, outpath):
    df.to_csv(outpath, sep='\t', index=False)


def readVariantDF(inpath):
    variant_df = pd.read_csv(inpath, sep='\t')
    variant_df = change_datatypes(variant_df)
    return variant_df


def df_creation(all_calls_dir, variant_df_outpath):
    called_df = get_position_call_df_all_samples(all_calls_dir)
    full_called_df = get_parsed_position_call_df(called_df, DEBUG=True)
    full_variant_df = filterVariants(full_called_df)
    saveVariantDF(full_variant_df, variant_df_outpath)


def analyse_variant_effect(variants_w_nearby_genes_df, ref_seq):
    variants_w_nearby_genes_df['CDS loc'] = "Intergenic"
    variants_w_nearby_genes_df['Effect'] = None
    variants_w_nearby_genes_df['Genbank nuc equals ref nuc'] = None
    variants_w_nearby_genes_df['Genbank aa equals ref nuc'] = None
    variants_w_nearby_genes_df['Ref Translation'] = None
    variants_w_nearby_genes_df['Alt Translation'] = None


    ['genome_name','POS','REF','ALT','QUAL','start','end','Name','locus_tag','product','translation', 'transl_table']

    for i, variant_gene in variants_w_nearby_genes_df.iterrows():
        var_pos = variant_gene['POS']
        gene_start = variant_gene['start']
        gene_end = variant_gene['end']

        if (gene_start <= var_pos <= gene_end):
            variants_w_nearby_genes_df.loc[i, 'CDS loc'] = "Intragenic"

            gene_CDS = ref_seq[gene_start-1:gene_end]
            var_gene_pos = var_pos - gene_start

            strand = variant_gene['strand']
            transl_table = variant_gene['transl_table']

            var_ref_seq = Seq(variant_gene['REF'])
            var_alt_seq = Seq(variant_gene['ALT'])
            
            start_seq = gene_CDS[:var_gene_pos]
            end_seq = gene_CDS[var_gene_pos+len(var_ref_seq):]

            nuc_seq_ref = start_seq + var_ref_seq + end_seq
            nuc_seq_alt = start_seq + var_alt_seq + end_seq

            variants_w_nearby_genes_df.loc[i, 'Genbank nuc equals ref nuc'] = nuc_seq_ref == gene_CDS

            if strand == '-':
                nuc_seq_ref = nuc_seq_ref.reverse_complement()
                nuc_seq_alt = nuc_seq_alt.reverse_complement()

            aa_seq_ref = nuc_seq_ref.translate(table=transl_table, stop_symbol="")
            aa_seq_alt = nuc_seq_alt.translate(table=transl_table, stop_symbol="")
            previously_translated_seq = variant_gene['translation']
            
            if aa_seq_ref == aa_seq_alt:
                variants_w_nearby_genes_df.loc[i, 'Effect'] = 'Synonymous'
            else:
                variants_w_nearby_genes_df.loc[i, 'Effect'] = 'Non-Synonymous'

            variants_w_nearby_genes_df.loc[i, 'Genbank aa equals ref nuc'] = previously_translated_seq == aa_seq_ref

            # if previously_translated_seq != aa_seq_ref:
            #     for tab in [1,2,3,4,5,6,9,10,11,12,13,14,16,21,22,23,24,25,26,27,28,29,30,31,33]:
            #         aa_seq_ref = nuc_seq_ref.translate(table=tab, stop_symbol="")
            #         print(tab, ":\t", previously_translated_seq == aa_seq_ref)

            variants_w_nearby_genes_df.loc[i, 'Alt Translation'] = str(aa_seq_alt)
            variants_w_nearby_genes_df.loc[i, 'Ref Translation'] = str(aa_seq_ref)
    
    return variants_w_nearby_genes_df


def df_analysis(full_variant_df_path, gff_path, ref_path, filtered_variant_df_outpath, variant_effect_df_outpath):
    full_variant_df = readVariantDF(full_variant_df_path)
    filtered_variant_df = filterOutNonSigVariants(full_variant_df)
    saveVariantDF(filtered_variant_df, filtered_variant_df_outpath)

    gff_df = get_gff_df(gff_path)
    gff_cds_df = gff_df[gff_df['type']=='CDS']
    ref_genome_record = SeqIO.read(ref_path, 'fasta')
    ref_genome_seq = ref_genome_record.seq

    print(gff_cds_df.columns)

    variants_w_nearby_genes_df = indicate_nearby_genes(filtered_variant_df, gff_cds_df)
    variant_effect_df = analyse_variant_effect(variants_w_nearby_genes_df, ref_genome_seq)
    saveVariantDF(variant_effect_df, variant_effect_df_outpath)
    