import pandas as pd

def check_tensor(df):
    for chrom in df.index.get_level_values('chrom').unique():
        temp_df = df.loc[chrom]
        idf = temp_df.index.to_frame()
        idf['in_index'] = True
        occurances = idf.pivot(index='pos',columns='sample',values='in_index')
        if not occurances.to_numpy().all():
            raise ValueError("not all samples are present at a variant position")

def decorate_mpileup(row):
    # catches exception where row has no data because there are no reads mapped to that position
    if row['cov'] == 0:
        return {}
    else:
        # creates temporary dataframe with each read as a row

        bases = list(row['bases'].replace('.',row['ref'].upper()).replace(',',row['ref'].lower()))
        quals = [ord(c) - 33 for c in list(row['quals'])]
        mapqs = [ord(c) - 33 for c in list(row['MAPQ'])]
        dists = [int(x) for x in row['tail_dist'].split(",")]

        tdf = pd.DataFrame(dict(
            base = bases, q = quals, m = mapqs, t = dists
        ))

        # filters based on tail distance
        tdf = tdf[
            (tdf['t'] < snakemake.config['candidate mutation table params']['maximum tail distance'])
            &
            (tdf['t'] > snakemake.config['candidate mutation table params']['minimum tail distance'])
        ]
        tdf = tdf[tdf['base']!='*']

        # creates a dataframe of means for quality score, mapq, and tail distances
        mean_df = tdf.groupby('base').mean()
        mean_df = mean_df.melt(ignore_index=False).set_index('variable',append=True)
        mean_df.index = mean_df.index.map(''.join)
        mean_dict = mean_df['value'].to_dict()

        # adds count information for each base in string
        info = tdf['base'].value_counts().to_dict()
        info.update(mean_dict)
        return info

df = pd.read_table(snakemake.input[0], index_col=['chrom','pos','sample'])

check_tensor(df)

info_df = pd.DataFrame.from_dict(
    {
        ix: decorate_mpileup(row)
        for ix, row in df.iterrows()
    },
    orient='index'
)

# sorts columns kinda okay
info_df = info_df.reindex(sorted(info_df.columns, key=lambda s: sum([ord(c) for c in list(s)])), axis=1)

# adds index names
info_df.index = info_df.index.set_names(['chrom','pos','sample'])

# joins candidate mutation table with info
df = df.join(info_df)

df.to_csv(snakemake.output[0], sep='\t')

# [A T C G a t c g Aq ... gq Am .... gm  At .... gt Ps Pb Pm Pftd Prtd E I D]
# List of statistics by index:
# DONE [1-4] A is the number of forward reads supporting A
# DONE [5-8] a is the number of reverse reads supporting A
# DONE [9-16] Aq is the average phred qualities of all A's
# DONE [17-24] Am is the average mapping qualities of all A's
# DONE [25-32] At is the average tail distance of all A's
# TO ADD [33] Ps is the p value for strand bias (fishers test)
# TO ADD [34] Pb is the p value for the base qualities being the same for the two different types of calls (1st major, 2nd major nt, either strand) (ttest)
# TO ADD [35] Pm is the p value for the mapping qualities being the same for the two different types of calls (ttest)
# TO ADD [36] Pftd is the p value for the tail distantces on the forward strand being the same for the two different types of calls (ttest)
# TO ADD [37] Pftd is the p value for the tail distantces on the reverse strand being the same for the two  different types of calls (ttest)
# TO ADD [38] E is number of calls at ends of a read
# TO ADD [39] I is number of reads supporting indels in the +/- (indelregion) bp region
# TO ADD [40] D is number of reads supporting deletions in the +/- (indelregion) bp region