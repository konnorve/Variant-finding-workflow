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

def gen_venn_set(included_vector, set_dict):
    included = set.intersection(*[set_dict[k] for k,i in zip(set_dict.keys(),included_vector) if i])
    if all(included_vector):
        return included
    else:
        excluded = set.union(*[set_dict[k] for k,i in zip(set_dict.keys(),included_vector) if not i])
        return included.difference(excluded)

def sets_to_venn(set_dict):
    venn_dict = dict()
    venn_len_arr = []
    set_keys = list(set_dict.keys())
    for i in range(1,len(set_keys)+1):
        for group in combinations(set_keys, i):
            included_vector = [k in set(group) for k in set_keys]
            venn_name = ", ".join([f"+{k}" if k in set(group) else f"-{k}" for k in set_keys])
            venn_set = gen_venn_set(included_vector, set_dict)
            venn_dict[venn_name] = venn_set
            venn_len_arr.append(included_vector + [len(venn_set)])
    return venn_dict, pd.DataFrame(venn_len_arr, columns=set_keys+['length'])

all_vcfs = []

Path(snakemake.output[0]).mkdir(parents=True, exist_ok=True)
for vcf_path in [Path(f) for f in snakemake.input]:
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

vcf_length_dfs = []

for sample in vcfs['SAMPLE'].unique():
    
    sample_df = vcfs[vcfs['SAMPLE']==sample]

    #filtering
    sample_df = sample_df[sample_df['QUAL']>20]

    set_dict = {method : set(sample_df[sample_df['METHOD']==method]['ID'].values) for method in sample_df['METHOD'].unique()}
    venn_dict, venn_df = sets_to_venn(set_dict)

    venn_df['sample'] = sample

    vcf_length_dfs.append(venn_df)

    fig = plt.figure(figsize=(10,10))
    venn(set_dict, ax=fig.gca())
    plt.savefig(f"{snakemake.output[0]}/{sample}_venn.png")
    plt.close()

pd.concat(vcf_length_dfs, ignore_index=True).to_csv(f'{snakemake.output[0]}/vcf_set_lengths.tsv', sep='\t', index=False)

