import pandas as pd
import plotly.graph_objects as go

image_outpath = snakemake.output[0]
phenotype_variants = pd.read_csv(snakemake.input[0])

heatmap_data = pd.get_dummies(phenotype_variants, prefix='', prefix_sep='', columns=['control']).groupby(['pheno']).sum()
heatmap_data = heatmap_data.sort_index(ascending=False)

fig = go.Figure()
fig.add_trace(
    go.Heatmap(
        z=heatmap_data.to_numpy(),
        zsmooth = 'best',
        x=heatmap_data.columns,
        y=heatmap_data.index,
    )
)
fig.update_layout(
    xaxis_title='Proportion of control lineages with mutated gene',
    yaxis_title='Proportion of phenotype lineages with mutated gene'
)
fig.write_image(image_outpath)