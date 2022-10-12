from pathlib import Path
import pandas as pd
import numpy as np
import gffpandas.gffpandas as gffpd
from Bio import SeqIO, pairwise2, SeqUtils

def get_unique_column_elements(df, dict_col, element_deliminator, key_value_deliminator):
    info_columns = set()
    for i, s in enumerate(df[dict_col].to_list()):
        try:
            for ele in s.split(element_deliminator):
                info_columns.add(ele[:ele.find(key_value_deliminator)])
        except:
            print("erroneous string", i, s)
            print(df.iloc[i])
            print(df)
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


def parse_vcf(vcf_path):
    with open(vcf_path) as f:
        for l in f:
            if l.startswith("#"):
                header = l[1:].split('\t')
    vcf_df = pd.read_table(vcf_path, names=header, comment='#')

    print(vcf_df)

    vcf_df = convert_dict_column_to_columns(vcf_df, 'INFO', ";", "=")

    if 'INFO:DP4' in vcf_df.columns:
        out_columns = [f"INFO:DP4:{s}" for s in ["REF_FWD","REF_REV","ALT_FWD","ALT_REV"]]
        vcf_df = convert_list_column_to_columns(vcf_df, 'INFO:DP4', out_columns, ",")
    if 'INFO:PV4' in vcf_df.columns:
        out_columns = [f"INFO:PV4:{s}" for s in ["STRAND_BIAS_P","BASQ_BIAS_P","MAPQ_BIAS_P","TAIL_DISTANCE_BIAS_P"]]
        vcf_df = convert_list_column_to_columns(vcf_df, 'INFO:PV4', out_columns, ",")

    vcf_df.loc[:,'SAMPLE'] = vcf_path.name.split('.')[0]
    vcf_df.loc[:,'METHOD'] = vcf_path.name.split('.')[1]

    return vcf_df

def get_annotation_df(gff_path):
    annotation = gffpd.read_gff3(gff_path)
    annotation_df = annotation.attributes_to_columns()
    return annotation_df

def get_ref_seq(ref_path):
    with open(ref_path) as in_seq_handle:
        return SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))

def determineEffect(aa_ref, aa_alt):
    """
    Returns effect type
    examples of effects: 
        Silent -- aa sequences are identical
        Misense -- aa sequences differ by a single aa. There can be multiple Misense mutations.
        Addition -- addition of an aa. There can be multiple additions.
        Deletion -- deletion of an aa. There can be multiple deletions.
        
        This method can not call frameshift or nonsense mutations well. That is best done
        by comparing an offset in nucleotide sequences. 
    """
    if aa_ref == aa_alt:
        return "Silent"

    # align pairs  
    alignment = pairwise2.align.globalms(aa_ref, aa_alt, 2, 0, -1, 0)[0]

    aligned_ref = alignment.seqA
    aligned_alt = alignment.seqB

    mismatches = []
    insertions = []
    deletions = []

    for i in range(len(aligned_ref)):
        if aligned_ref[i] != aligned_alt[i]:
            if aligned_ref[i] == "-":
                insertions.append(i)
            elif aligned_alt[i] == "-":
                deletions.append(i)
            else:
                mismatches.append(i)

    i = 0
    insertion_locs = []
    while i < len(insertions):
        start = insertions[i]
        while i + 1 < len(insertions) and insertions[i] + 1 == insertions[i+1]:
            i += 1
        end = insertions[i]
        i += 1
        insertion_locs.append((start, end))

    i = 0
    deletion_locs = []
    while i < len(deletions):
        start = deletions[i]
        while i + 1 < len(deletions) and deletions[i] + 1 == deletions[i+1]:
            i += 1
        end = deletions[i]
        i += 1
        deletion_locs.append((start, end))

    mismatch_str = None
    if len(mismatches) != 0:
        mismatch_strs = [f"{SeqUtils.seq3(aligned_ref[n])}->{SeqUtils.seq3(aligned_alt[n])} at {n}" for n in mismatches]
        mismatch_str = "mismatches at " + ", ".join(mismatch_strs)

    insertion_str = None
    if len(insertion_locs) != 0:
        insertion_strs = [f"{SeqUtils.seq3(aligned_alt[s:e+1])} from {s} to {e}" for s, e in insertion_locs]
        insertion_str = "insertion of " + ", ".join(insertion_strs)
    
    deletion_str = None
    if len(deletion_locs) != 0:
        deletion_strs = [f"{SeqUtils.seq3(aligned_ref[s:e+1])} from {s} to {e}" for s, e in deletion_locs]
        deletion_str = "deletion of " + ", ".join(deletion_strs)

    effect = "; ".join([s for s in [mismatch_str, insertion_str, deletion_str] if s])

    return effect

def filter_overlapping_variants(vcf_df):
    var_groups = {}
    for ix, var in vcf_df.iterrows():
        if var['CHROM'] not in var_groups.keys():
            var_groups[var['CHROM']] = []
        var_set = set(range(var['POS'],var['POS']+len(var['REF'])))

        no_overlap = True
        for g in var_groups[var['CHROM']]:
            if g['set'].intersection(var_set):
                g['set'] |= var_set
                if var['QUAL'] > g['qual']:
                    g['qual'] = var['QUAL']
                    g['ix'] = ix
                no_overlap = False
                break
        if no_overlap:
            var_groups[var['CHROM']].append({
                'set': var_set,
                'ix': ix,
                'qual': var['QUAL'],
            })
    
    indicies2keep = [   
        group['ix']
        for chrom in var_groups.keys()
        for group in var_groups[chrom]
    ]

    return vcf_df.loc[indicies2keep]
                


def applyVariants2Seq(variant_subset, nuc_seq_ref):
    nuc_seq_alt = nuc_seq_ref
    arr = variant_subset[['POS', 'REF', 'ALT']].to_numpy()
    delta_sum = 0
    for i in range(len(arr)):
        p, r, a = arr[i]
        lr = len(r)
        la = len(a)
        delta = la - lr
        delta_sum += delta
        nuc_seq_alt = nuc_seq_alt[:p] + a + nuc_seq_alt[p+lr:]
        arr[i:, 0] = arr[i:, 0] + delta
        assert len(nuc_seq_ref) == len(nuc_seq_alt) - delta_sum
    return nuc_seq_alt
        

# read in variant DF
def main(vcf_path, ref_path, gff_path, gene_table_output):
    vcf_df = parse_vcf(Path(vcf_path))
    vcf_df = filter_overlapping_variants(vcf_df)
    annotation_df = get_annotation_df(gff_path)
    ref_seq_dict = get_ref_seq(ref_path)

    acceptible_gene_types = ['CDS']
    filtered_annotation_df = annotation_df[annotation_df['type'].isin(acceptible_gene_types)]
    varied_genes=[]
    for ix, gene in filtered_annotation_df.iterrows():

        gene_variants = vcf_df[(vcf_df['CHROM']==gene['seq_id']) & (vcf_df['POS']>=gene['start']) & (vcf_df['POS']<=gene['end'])]
        gene_variants.loc[:,'POS'] = gene_variants['POS'] - gene['start']

        if len(gene_variants) > 0:
            varied_genes.append(ix)

            nuc_seq_ref = ref_seq_dict[gene['seq_id']].seq[gene['start']-1:gene['end']]
            nuc_seq_alt = applyVariants2Seq(gene_variants, nuc_seq_ref)
            
            if gene['strand'] == '-':
                nuc_seq_ref = nuc_seq_ref.reverse_complement()
                nuc_seq_alt = nuc_seq_alt.reverse_complement()

            transl_table=11
            aa_seq_ref = nuc_seq_ref.translate(table=transl_table, stop_symbol="")
                    
            excess_nucs = len(nuc_seq_alt)%3
            try:
                if excess_nucs != 0:
                    aa_seq_alt = nuc_seq_alt[:-excess_nucs].translate(table=transl_table, stop_symbol="")
                else:
                    aa_seq_alt = nuc_seq_alt.translate(table=transl_table, stop_symbol="")
            except TranslationError as e:
                logging.error(f"nuc_seq_alt\n{nuc_seq_alt}")
                raise e

            annotation_df.loc[ix, 'synonymous'] = aa_seq_alt == aa_seq_ref

            effect = None
            if excess_nucs != 0:
                k = 0
                while k < len(aa_seq_alt) - 1 and k < len(aa_seq_ref) - 1:
                    if (aa_seq_alt[k] != aa_seq_ref[k]) and (aa_seq_alt[k+1] != aa_seq_ref[k+1]):
                        break
                    k += 1
                
                effect = f"frameshift starting at {k}. alt translation {len(aa_seq_alt)} aas, ref translation {len(aa_seq_ref)} aas"
            else:
                effect = determineEffect(aa_seq_ref, aa_seq_alt)

            annotation_df.loc[ix, 'Effect'] = effect

    annotation_df.loc[varied_genes].to_csv(gene_table_output, sep='\t', index=False)


main(snakemake.input["vcf"], snakemake.input["ref"], snakemake.input["gff"], snakemake.output[0])