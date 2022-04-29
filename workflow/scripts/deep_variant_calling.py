from os import name
from pathlib import Path
import pandas as pd
import numpy as np

import gffpandas.gffpandas as gffpd
from Bio import SeqIO, pairwise2
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import seq3
from BCBio import GFF
from Bio.Seq import MutableSeq, Seq
from dna_features_viewer import BiopythonTranslator
from Bio.SeqFeature import SeqFeature, FeatureLocation

import matplotlib.pyplot as plt
from textwrap import wrap


def get_unique_INFO_elements(vcf_df):
    info_list = vcf_df['INFO'].to_numpy()
    
    info_columns = []
    for i in range(info_list.size):
        s = info_list[i]
        for ele in s.split(';'):
            col = ele[:ele.find('=')]
            if col not in info_columns:
                info_columns.append(col)

    return info_columns


def convert_dict_column_to_columns(df, col, element_deliminator, key_value_deliminator):
    
    out_col_names = get_unique_INFO_elements(vcf_df)

    col_list = df[col].to_numpy()
    col_arr = np.empty((len(col_list), len(out_col_names)), dtype=object)
    for i in range(len(col_list)):
        row_dict = dict(ele.split("=") for ele in col_list[i].split(element_deliminator) if len(ele.split(key_value_deliminator))==2)
        for j, col in enumerate(out_col_names):
            col_arr[i][j] = row_dict.get(col)

    df[out_col_names] = col_arr

    return df


def convert_list_column_to_columns(df, col, out_col_names, deliminator):
    col_list = df[col].to_numpy()
    col_arr = np.empty((len(col_list), len(out_col_names)), dtype=object)
    r = range(len(out_col_names))
    for i in range(len(col_list)):
        s = col_list[i]
        if s is None:
            elements = [None] * len(out_col_names)
        else:
            elements = s.split(deliminator)
        for j in r:
            col_arr[i][j] = elements[j]

    df[out_col_names] = col_arr

    return df


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


def get_parsed_vcf_df(vcf_df):

    vcf_df = convert_dict_column_to_columns(vcf_df, 'INFO', ';', '=')

    # must convert DP4 and PV4 to columns
    dp4_columns = ['ref_forward_reads', 'ref_reverse_reads', 'alt_forward_reads', 'alt_reverse_reads']
    vcf_df = convert_list_column_to_columns(vcf_df, 'DP4', dp4_columns, ',')

    pv4_columns = ['strand_bias_pVal', 'baseQ_bias_pVal', 'mapQ_bias_pVal', 'tail_dist_bias_pVal']
    vcf_df = convert_list_column_to_columns(vcf_df, 'PV4', pv4_columns, ',')

    vcf_parsed_df = change_datatypes(vcf_df)

    vcf_parsed_df['fwd_strand_coverage'] = vcf_parsed_df['ref_forward_reads'] + vcf_parsed_df['alt_forward_reads']
    vcf_parsed_df['rev_strand_coverage'] = vcf_parsed_df['ref_reverse_reads'] + vcf_parsed_df['alt_reverse_reads']

    print(vcf_parsed_df.dtypes)
    for col in vcf_parsed_df.columns:
        print(col, '\t', vcf_parsed_df[col].unique()[:5])

    return vcf_parsed_df


def get_vcf_df(all_calls_path):
    raw_df = pd.read_csv(all_calls_path, sep='\t')
    narrow_df = raw_df[['POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']]
    return narrow_df


def get_gff_df(gff_path):
    annotation = gffpd.read_gff3(gff_path)
    attributes_df = annotation.attributes_to_columns()
    return attributes_df


def agg_vcf_df(files):
    dfs2concat = []
    for sample_path in files:
        sample_name = Path(sample_path).name.split('.')[0]
        sample_df = get_vcf_df(sample_path)
        col_order = list(sample_df.columns)
        sample_df['genome_name'] = sample_name
        col_order.insert(0, 'genome_name')
        sample_df = sample_df[col_order]
        print(sample_df.columns)
        dfs2concat.append(sample_df)
    return pd.concat(dfs2concat, axis=0, ignore_index=True)


def removeNonVariants(all_vcfs_df):
    return all_vcfs_df[all_vcfs_df['ALT'] != '.']


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


def df_creation(vcf_files, variant_df_outpath):
    called_df = agg_vcf_df(vcf_files)
    full_called_df = get_parsed_vcf_df(called_df, DEBUG=True)
    full_variant_df = removeNonVariants(full_called_df)
    saveVariantDF(full_variant_df, variant_df_outpath)


def applyVariants2Genome(variant_df, ref_seq, gff_df):
    gff_alt_df = gff_df.copy()
    alt_seq = MutableSeq(ref_seq)
    assert alt_seq is not ref_seq

    start_pos = gff_alt_df['start'].to_numpy()
    end_pos = gff_alt_df['end'].to_numpy()
    num_genes = len(gff_alt_df)

    offset_sum = 0

    for i, variant in variant_df.iterrows():
        var_pos = variant['POS'] + offset_sum - 1
        var_ref_seq = variant['REF']
        var_alt_seq = variant['ALT']

        start_seq = alt_seq[:var_pos]
        end_seq = alt_seq[var_pos+len(var_ref_seq):]

        alt_seq = start_seq + var_alt_seq + end_seq

        offset = len(var_alt_seq) - len(var_ref_seq)

        offset_sum += offset

        if offset != 0:
            for j in range(num_genes):
                if start_pos[j] > var_pos:
                    start_pos[j] += offset
                if end_pos[j] > var_pos:
                    end_pos[j] += offset

    gff_alt_df['start'] = start_pos
    gff_alt_df['end'] = end_pos

    assert len(ref_seq) == len(alt_seq) - offset_sum

    return gff_alt_df, Seq(alt_seq)


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
        mismatch_strs = [f"{seq3(aligned_ref[n])}->{seq3(aligned_alt[n])} at {n}" for n in mismatches]
        mismatch_str = "mismatches at " + ", ".join(mismatch_strs)

    insertion_str = None
    if len(insertion_locs) != 0:
        insertion_strs = [f"{seq3(aligned_alt[s:e+1])} from {s} to {e}" for s, e in insertion_locs]
        insertion_str = "insertion of " + ", ".join(insertion_strs)
    
    deletion_str = None
    if len(deletion_locs) != 0:
        deletion_strs = [f"{seq3(aligned_ref[s:e+1])} from {s} to {e}" for s, e in deletion_locs]
        deletion_str = "deletion of " + ", ".join(deletion_strs)

    effect = "; ".join([s for s in [mismatch_str, insertion_str, deletion_str] if s])

    return effect

def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)

def df2gff3(annotation_df, outpath, genome_record):
    gff_df = annotation_df[['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']]
    gff_df['seq_id'] = genome_record.id
    gff_df.to_csv(outpath, sep='\t', index=False, header=False)
    line_prepender(outpath, f"##sequence-region {genome_record.id} 1 {len(genome_record.seq)}")
    line_prepender(outpath, "##gff-version 3")

def analyse_variant_effects(filtered_variant_df, ref_record, ref_gff_df, alt_genome_dir=None, plotting_dir=None):
    ref_seq = ref_record.seq
    
    unique_genomes = filtered_variant_df['genome_name'].unique()

    num_genes = len(ref_gff_df)

    dfs_to_concat = []

    graphic_record = BiopythonTranslator().translate_record(ref_record)

    for unique_genome in unique_genomes:
        print(unique_genome)

        genome_variants = filtered_variant_df[filtered_variant_df['genome_name'] == unique_genome]

        genome_df = ref_gff_df.copy()
        genome_df['genome_name'] = unique_genome
        
        # genome_df['Nearby Variant Ref Positions'] = ""
        genome_df['Intragenic Variant Ref Positions'] = ""
        genome_df['Genbank aa equals ref aa'] = None
        genome_df['ref nuc equals alt nuc'] = None
        genome_df['Ref Translation'] = None
        genome_df['Alt Translation'] = None
        genome_df['ref aa equals alt aa'] = None
        genome_df['Effect'] = None

        alt_gff_df, alt_seq = applyVariants2Genome(genome_variants, ref_seq, ref_gff_df)

        if alt_genome_dir:
            alt_seq_id = f"{ref_record.id}_{unique_genome}_alt"
            alt_seq_record = SeqRecord(alt_seq, id=alt_seq_id)

            gff_out_path = alt_genome_dir / f"{unique_genome}_alt.gff3"
            alt_ref_path = alt_genome_dir / f"{unique_genome}_alt.fasta"
            df2gff3(alt_gff_df, gff_out_path, alt_seq_record)
            
            SeqIO.write(alt_seq_record, alt_ref_path, "fasta")

        assert len(ref_gff_df) == len(alt_gff_df) == len(genome_df)

        for i in range(num_genes):
            ref_start = int(genome_df.loc[i, 'start'])
            alt_start = int(alt_gff_df.loc[i, 'start'])
            ref_end = int(genome_df.loc[i, 'end'])
            alt_end = int(alt_gff_df.loc[i, 'end'])

            intragenic = False
            intragenic_vars = []
            for j, variant in genome_variants.iterrows():
                var_pos = variant['POS'] 

                if ref_start <= var_pos <= ref_end:
                    intragenic = True
                    genome_df.loc[i, 'Intragenic Variant Ref Positions'] = str(var_pos) + ", " + str(genome_df.loc[i, 'Intragenic Variant Ref Positions'])
                    intragenic_vars.append(var_pos)

                # elif ((ref_end > var_pos - 1000) & (ref_end < var_pos + 1000)) | ((ref_start > var_pos - 1000) & (ref_start < var_pos + 1000)):
                #     genome_df.loc[j, 'Nearby Variant Ref Positions'] = str(var_pos) + ", " + genome_df.loc[j, 'Nearby Variant Ref Positions']

            if intragenic:
                print(i)

                strand = genome_df.loc[i, 'strand']
                transl_table = genome_df.loc[i, 'transl_table']
                aa_seq_genbank = genome_df.loc[i, 'translation']

                nuc_seq_ref = ref_seq[ref_start-1:ref_end]
                nuc_seq_alt = alt_seq[alt_start-1:alt_end]

                offset = len(nuc_seq_alt) - len(nuc_seq_ref)

                genome_df.loc[i, 'ref nuc equals alt nuc'] = nuc_seq_ref == nuc_seq_alt

                if strand == '-':
                    nuc_seq_ref = nuc_seq_ref.reverse_complement()
                    nuc_seq_alt = nuc_seq_alt.reverse_complement()

                aa_seq_ref = nuc_seq_ref.translate(table=transl_table, stop_symbol="")
                aa_seq_alt = nuc_seq_alt.translate(table=transl_table, stop_symbol="")

                genome_df.loc[i, 'Genbank aa equals ref aa'] = aa_seq_genbank[1:] == aa_seq_ref[1:]
                genome_df.loc[i, 'Ref Translation'] = str(aa_seq_ref)
                genome_df.loc[i, 'Alt Translation'] = str(aa_seq_alt)
                genome_df.loc[i, 'ref aa equals alt aa'] = aa_seq_ref == aa_seq_alt

                effect = None
                if offset % 3 != 0:
                    k = 0
                    while k < len(aa_seq_alt) and k < len(aa_seq_ref):
                        if aa_seq_alt[k] != aa_seq_ref[k]:
                            break
                        k += 1
                    effect = "frameshift starting at {} in translated sequence".format(k)
                else:
                    effect = determineEffect(aa_seq_ref, aa_seq_alt)

                genome_df.loc[i, 'Effect'] = effect
                
                print(effect)

                if plotting_dir is not None and effect != "Silent":
                    plot_intragenic_variant(graphic_record, unique_genome, genome_df.loc[i], effect, 
                                            intragenic_vars, plotting_dir, WINDOW_SIZE=5000)

                print()


        # filter genome_df for just nearby or intergenic rows
        # genome_df = genome_df[(genome_df['Nearby Variant Ref Positions'] != "") | (genome_df['Intragenic Variant Ref Positions'] != "")]
        genome_df = genome_df[genome_df['Intragenic Variant Ref Positions'] != ""]
        dfs_to_concat.append(genome_df)

    # concatenate genome_dfs
    concatenated_genomes_df = pd.concat(dfs_to_concat, axis=0)

    # rename columns
    columns2keep = ['genome_name','seq_id','type','start','end','strand','Intragenic Variant Ref Positions','Genbank aa equals ref aa',
                    'ref nuc equals alt nuc','ref aa equals alt aa','Effect','Ref Translation','Alt Translation',
                    'Name','locus_tag','note','product','translation']
    sliced_concat_genomes_df = concatenated_genomes_df[columns2keep]

    return sliced_concat_genomes_df


def plot_intergenic_variants(filtered_variant_df, ref_gff_df, ref_record, figure_out_dir, WINDOW_SIZE = 5000):

    unique_genomes = filtered_variant_df['genome_name'].unique()

    graphic_record = BiopythonTranslator().translate_record(ref_record)

    all_intergenic_variant_dfs = []

    for unique_genome in unique_genomes:
        print(unique_genome)

        genome_variants = filtered_variant_df[filtered_variant_df['genome_name'] == unique_genome]

        for j, variant in genome_variants.iterrows():
            var_pos = variant['POS'] 
            
            genome_df = ref_gff_df.copy()
            genome_df['var_right_of_start'] = genome_df['start'] <= var_pos
            genome_df['var_left_of_end'] = var_pos <=genome_df['end']
            genome_df = genome_df[genome_df['var_right_of_start'] & genome_df['var_left_of_end']]

            if len(genome_df) == 0:
                var_effect = f"{variant['REF']} to {variant['ALT']}"
                plot_intergenic_variant(graphic_record, unique_genome, var_pos, var_effect, figure_out_dir, WINDOW_SIZE)
            
        all_intergenic_variant_dfs.append(genome_variants)
    
    return pd.concat(all_intergenic_variant_dfs)

def plot_intergenic_variant(graphic_record, unique_genome, var_pos, var_effect, figure_out_dir, WINDOW_SIZE=5000):

    out_path = figure_out_dir / f"{unique_genome}_{var_pos:08}_{WINDOW_SIZE:06}.png"
    
    if WINDOW_SIZE >= 2000:
        height = 3*WINDOW_SIZE/5000 
    else:
        height = 3

    fig, ax = plt.subplots(1, 1, figsize=(10, height))
    fig.suptitle(f"{unique_genome} variant {var_pos} from {var_effect}", size="xx-large")

    cropped_record = graphic_record.crop((var_pos-WINDOW_SIZE, var_pos+WINDOW_SIZE))

    ax.axvline(var_pos, c='r')
    cropped_record.plot(ax=ax, figure_width=10, strand_in_label_threshold=7)

    plt.savefig(out_path)
    plt.close()

def plot_intragenic_variant(graphic_record, unique_genome, gene_series, var_effect, var_pos_list, figure_out_dir, WINDOW_SIZE=5000):

    out_path = figure_out_dir / f"{unique_genome}_{gene_series['locus_tag']}_{WINDOW_SIZE:06}.png"
    
    if WINDOW_SIZE >= 2000:
        height = 3*WINDOW_SIZE/5000 
    else:
        height = 3

    fig, ax = plt.subplots(1, 1, figsize=(10, height))

    long_title = f"{unique_genome} gene {gene_series['Name']} with {var_effect}"
    title = "\n".join(wrap(long_title, 90))

    fig.suptitle(title, size="x-large")

    cropped_record = graphic_record.crop((gene_series['start']-WINDOW_SIZE, gene_series['end']+WINDOW_SIZE))

    for var_pos in var_pos_list:
        ax.axvline(var_pos, c='r')
    cropped_record.plot(ax=ax, figure_width=10, strand_in_label_threshold=7)

    plt.savefig(out_path)
    plt.close()

def analyze_chi_sites(filtered_variant_df, chi_gff_df):
    
    chi_gff_df["midpoint"] = (chi_gff_df["start"] + chi_gff_df["end"]) / 2
    chi_sites = chi_gff_df["midpoint"].to_numpy()

    unique_genomes = filtered_variant_df['genome_name'].unique()

    var2chi_dfs = []

    for unique_genome in unique_genomes:
        print(unique_genome)

        genome_variants = filtered_variant_df[filtered_variant_df['genome_name'] == unique_genome]

        var_pos = genome_variants['POS'].to_numpy()

        var_dist2chi = np.zeros(shape=(len(var_pos)))

        for i in range(len(var_pos)):
            var_dist2chi[i] = np.min(np.abs(chi_sites - var_pos[i]))

        genome_variants['dist_to_chi_site'] = var_dist2chi

        var2chi_dfs.append(genome_variants)
    
    return pd.concat(var2chi_dfs)


def chi_summary_statistics(filtered_variant_df_chi):

    unique_genomes = filtered_variant_df_chi['genome_name'].unique()

    dist_dict = {}

    for unique_genome in unique_genomes:

        genome_variants = filtered_variant_df_chi[filtered_variant_df_chi['genome_name'] == unique_genome]

        chi_distances = genome_variants['dist_to_chi_site'].to_numpy()

        dist_dict[unique_genome] = np.median(chi_distances)

    return dist_dict

    
        
def get_chi_distribution(chi_gff_df, genome_len):
    # calculate distribution for genome

    chi_gff_df["midpoint"] = (chi_gff_df["start"] + chi_gff_df["end"]) / 2
    chi_sites = chi_gff_df["midpoint"].to_numpy()

    dist2chi = np.zeros(shape=(genome_len))

    for i in range(genome_len):
        dist2chi[i] = int(np.min(np.abs(chi_sites - i + 1)))
        if i % 100000 == 0: print(i,'\t',dist2chi[i])

    print(f"median dist to chi whole genome:\t{np.sort(dist2chi)[int(len(dist2chi)/2)]}")

    dist2chi_distribution = np.zeros(shape=(int(max(dist2chi)+1)))

    unique, counts = np.unique(dist2chi, return_counts=True)
    distribution_dictionary = dict(zip(unique, counts))

    for i in range(len(dist2chi_distribution)):
        if i in distribution_dictionary.keys():
            dist2chi_distribution[i] = distribution_dictionary[i]

    return dist2chi_distribution

def plot_chi_distribution(figure_outpath, chi_dist_distribution, average_feature_dist_dict=None):

    fig, ax = plt.subplots(1, 1, figsize=(30, 5))

    ax.fill_between(range(len(chi_dist_distribution)), 0, chi_dist_distribution, label="distance to chi site (nt)")
    
    ax.set_xlabel("distance from chi site")

    ax.set_ylabel("count of nucleotides with this distance")

    sum_distribution = np.sum(chi_dist_distribution*np.arange(len(chi_dist_distribution)))
    count_distribution = np.sum(chi_dist_distribution)
    avg_distribution = sum_distribution / count_distribution

    ax.axvline(avg_distribution, label="average dist to chi across genome", c='k')
    print(f"average dist to chi whole genome:\t{avg_distribution}")

    if average_feature_dist_dict is not None:
        for key, value in average_feature_dist_dict.items():
            color = np.random.rand(3) / 2 + 0.25
            ax.axvline(value, label=key, c=color)
            print(f"median dist var to chi {key}:\t{value}")
        ax.legend()

    plt.savefig(figure_outpath)
    plt.close()

def get_seq_record(fasta_path, gff_path):
    
    limit_info = dict(gff_type=["CDS"])

    in_seq_handle = open(fasta_path)
    seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
    in_seq_handle.close()

    record = None
    in_handle = open(gff_path)
    for rec in GFF.parse(in_handle, base_dict=seq_dict, limit_info=limit_info):
        record = rec
    in_handle.close()

    return record


def df_filtered_analysis(full_variant_df_path, gff_path, ref_path, filtered_variant_df_outpath, variant_effect_df_outpath, alt_genome_dir=None):
    full_variant_df = readVariantDF(full_variant_df_path)
    filtered_variant_df = filterOutNonSigVariants(full_variant_df)
    saveVariantDF(filtered_variant_df, filtered_variant_df_outpath)
    df_analysis(filtered_variant_df, gff_path, ref_path, variant_effect_df_outpath, alt_genome_dir)


def df_standard_analysis(filtered_variant_df_path, gff_path, ref_path, variant_effect_df_outpath=None, alt_genome_dir=None, intergenic_figure_dir=None, intragenic_figure_dir=None, intergenic_df_outpath=None, chi_distribution_figure_outpath=None):
    filtered_variant_df = readVariantDF(filtered_variant_df_path)
    df_analysis(filtered_variant_df, gff_path, ref_path, variant_effect_df_outpath, alt_genome_dir, intergenic_figure_dir, intragenic_figure_dir, intergenic_df_outpath, chi_distribution_figure_outpath)


def df_analysis(filtered_variant_df, gff_path, ref_path, variant_effect_df_outpath=None, alt_genome_dir=None, intergenic_figure_dir=None, intragenic_figure_dir=None, intergenic_df_outpath=None, chi_distribution_figure_outpath=None):
    gff_df = get_gff_df(gff_path)

    gff_cds_df = gff_df[gff_df['type']=='CDS']
    gff_cds_df = gff_cds_df.reset_index()

    gff_chi_df = gff_df[gff_df['type']=='CHI_Pro']
    gff_chi_df = gff_chi_df.reset_index()

    ref_genome_record = get_seq_record(ref_path, gff_path)

    if chi_distribution_figure_outpath is not None:

        var2chi_df = analyze_chi_sites(filtered_variant_df, gff_chi_df)

        average_feature_dist_dict = chi_summary_statistics(var2chi_df)

        dist2chi_distribution = get_chi_distribution(gff_chi_df, len(ref_genome_record.seq))

        plot_chi_distribution(chi_distribution_figure_outpath, dist2chi_distribution, average_feature_dist_dict=average_feature_dist_dict)

    if intergenic_figure_dir is not None:
        for ws in [5000, 10000]:
            intergenic_variants_df = plot_intergenic_variants(filtered_variant_df, gff_cds_df, ref_genome_record, intergenic_figure_dir, WINDOW_SIZE = ws)

            if intergenic_df_outpath is not None:
                saveVariantDF(intergenic_variants_df, intergenic_df_outpath)

    if variant_effect_df_outpath is not None and intragenic_figure_dir is not None:
        variant_effect_df = analyse_variant_effects(filtered_variant_df, ref_genome_record, gff_cds_df, alt_genome_dir, intragenic_figure_dir)
        saveVariantDF(variant_effect_df, variant_effect_df_outpath)
    