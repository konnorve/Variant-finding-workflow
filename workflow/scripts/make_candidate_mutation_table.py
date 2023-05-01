import pandas as pd
from pathlib import Path
import logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

# loads variant positions
var_pos_idx = pd.MultiIndex.from_frame(pd.read_table(snakemake.input['var_pos']))
logging.info(f"loaded var_pos_idx")

pileups = []
for pileup in snakemake.input['pileups']:
    # loads pileup table
    tdf = pd.read_table(pileup)
    logging.info(f"loaded {Path(pileup).stem}")

    # adds sample to table
    tdf['sample'] = Path(pileup).stem

    # sets index
    tdf = tdf.set_index(['chrom','pos','sample'])

    # filters index
    tdf = tdf[tdf.index.droplevel('sample').isin(var_pos_idx)]
    logging.info(f"filtered {pileup} by variant positions")

    pileups.append(tdf)

df = pd.concat(pileups)
logging.info(f"concatenated mpileups")

# indexes dataframe on scaffold, position, and sample
df = df.sort_index()

# saves positions to table
logging.info(f"saving final table")
df.to_csv(snakemake.output['table'], sep='\t')