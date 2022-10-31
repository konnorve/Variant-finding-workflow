import pandas as pd
import plotly.graph_objects as go
import numpy as np

html_outpath = snakemake.output['html']
png_outpath = snakemake.output['png']
phenotype_variants = pd.read_table(snakemake.input[0], index_col=[0,1])

heatmap_data = pd.get_dummies(phenotype_variants, prefix='', prefix_sep='', columns=['control']).groupby(['pheno']).sum()
heatmap_data = heatmap_data.sort_index(ascending=False)

fig = go.Figure()
fig.add_trace(
    go.Heatmap(
        z=np.log1p(heatmap_data.to_numpy()),
        zsmooth = 'best',
        x=heatmap_data.columns,
        y=heatmap_data.index,
        colorscale = 'tempo'
    )
)

fig.update_layout(
    xaxis_title='Proportion of control lineages with mutated gene',
    yaxis_title='Proportion of phenotype lineages with mutated gene'
)

fig.write_html(html_outpath)
fig.write_image(png_outpath, width=1500, height=1500)