import pandas as pd
from pathlib import Path
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt
from venn import venn

def get_unique_column_elements(df, dict_col, element_deliminator, key_value_deliminator):
    info_columns = set()
    for i, s in enumerate(df[dict_col].to_list()):
        try:
            for ele in s.split(element_deliminator):
                info_columns.add(ele[:ele.find(key_value_deliminator)])
        except:
            print(i, s, df.iloc[i])
            raise
    return [f"{dict_col}:{split_col}" for split_col in info_columns]

def convert_dict_column_to_columns(df, dict_col, element_deliminator, key_value_deliminator):
    out_col_names = get_unique_column_elements(df, dict_col, element_deliminator, key_value_deliminator)
    col_arr = np.empty((len(df), len(out_col_names)), dtype=object)
    for i, s in enumerate(df[dict_col].to_list()):
        row_dict = dict(ele.split(key_value_deliminator,2) for ele in s.split(element_deliminator) if len(ele.split(key_value_deliminator,2))==2)
        for j, col in enumerate(out_col_names):
            col_arr[i][j] = row_dict.get(col)

    df[out_col_names] = col_arr
    return df

def convert_list_column_to_columns(df, col, out_col_names, deliminator):
    col_arr = np.empty((len(df), len(out_col_names)), dtype=object)
    for i, s in enumerate(df[col].to_list()):
        if s is None:
            col_arr[i][:] = None
        else:
            for j, ele in enumerate(s.split(deliminator)):
                col_arr[i][j] = ele

    df[out_col_names] = col_arr
    return df


def change_datatypes(df):
    all_columns_datatype_dict = {'POS':'uint64', 'ID':'str_', 'REF':'str_', 'ALT':'str_', 'QUAL':"float64", 
                            'FILTER':"str_", 'INFO':"str_", 'FORMAT':"str_"}
    df = df.replace('.', None)
    for k,v in all_columns_datatype_dict.items():
        if k in list(df.columns):
            df = df.astype({k:v})
    return df


def venn_set(include_sets, exclude_sets, set_dict):
    included = set.intersection(*[set_dict[i] for i in include_sets])
    excluded = set.union(*[set_dict[e] for e in exclude_sets])
    return included.difference(excluded)

def venn_name(include_sets, exclude_sets, sep=', '):
    return sep.join([f"+{i}" for i in include_sets]+[f"-{e}" for e in exclude_sets])

def sets_to_venn(set_dict):
    venn_dict = dict()
    set_keys = set(set_dict.keys())
    for i in range(1,len(set_keys)+1):
        for group in combinations(set_keys, i):
            include = set(group)
            exclude = set_keys-set(group)
            if len(exclude) > 0:
                venn_dict[venn_name(include, exclude)] = venn_set(include, exclude, set_dict)
            else:
                name = ", ".join([f"+{i}" for i in include])
                venn_dict[name] = set.intersection(*[set_dict[i] for i in include])
    return venn_dict

all_vcfs = []

for vcf_path in Path("/nfs/chisholmlab001/kve/2022_SNPs_Dark_Adapted_Genomes/results/raw_data/vcfs").glob("N2_A1A_afterT1_A/*.vcf"):
    print(vcf_path)
    metadata = []
    with open(vcf_path) as f:
        for l in f:
            if l.startswith("##"):
                metadata.append(l[2:])
            if l.startswith("#"):
                header = l[1:-1].split('\t')
                
    temp_df = pd.read_table(vcf_path, names=header, comment='#')
    temp_df = convert_dict_column_to_columns(temp_df, 'INFO', ";", "=")

    if 'INFO:DP4' in temp_df.columns:
        out_columns = [f"INFO:DP4:{s}" for s in ["REF_FWD","REF_REV","ALT_FWD","ALT_REV"]]
        temp_df = convert_list_column_to_columns(temp_df, 'INFO:DP4', out_columns, ",")
    if 'INFO:PV4' in temp_df.columns:
        out_columns = [f"INFO:PV4:{s}" for s in ["STRAND_BIAS_P","BASQ_BIAS_P","MAPQ_BIAS_P","TAIL_DISTANCE_BIAS_P"]]
        temp_df = convert_list_column_to_columns(temp_df, 'INFO:PV4', out_columns, ",")

    temp_df.loc[:,'SAMPLE'] = vcf_path.name.split('.')[0]
    temp_df.loc[:,'METHOD'] = vcf_path.name.split('.')[1]

    all_vcfs.append(temp_df)

vcfs = pd.concat(all_vcfs, ignore_index=True)

vcfs = vcfs[[col for col in vcfs.columns if '.bam' not in col]]

vcfs.loc[:,'ID'] = vcfs.apply(lambda r: "|".join([str(r[i]) for i in ['POS', 'REF', 'ALT']]), axis=1)

vcfs = change_datatypes(vcfs)

for sample in vcfs['SAMPLE'].unique():
    print(sample)
    sample_df = vcfs[vcfs['SAMPLE']==sample]

    #filtering
    sample_df = sample_df[sample_df['QUAL']>20]

    set_dict = {method : set(sample_df[sample_df['METHOD']==method]['ID'].values) for method in sample_df['METHOD'].unique()}
    venn_dict = sets_to_venn(set_dict)

    for k,s in set_dict.items():
        print(k, len(s), sep='\t')
    venn_sizes = {k:len(s) for k,s in venn_dict.items()}
    venn_sizes = dict(sorted(venn_sizes.items(), key=lambda i: i[1]))
    for k,l in venn_sizes.items():
        if l > 0:
            print(k, l, sep='\t')

    fig = plt.figure()
    venn(set_dict, ax=fig.gca())
    plt.savefig(f"{sample}_venn.png")
