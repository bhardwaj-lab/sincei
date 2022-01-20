# for bokeh plotting
from bokeh.models import ColumnDataSource, Div, Select, Slider, TextInput, Legend
from bokeh.io import curdoc
from bokeh.embed import components
from bokeh.plotting import figure, output_file, show
from bokeh.palettes import Category20, Blues
from bokeh.transform import factor_cmap, linear_cmap
import pandas as pd
import subprocess

# execute a command using args passed from WTforms
def execute_command(command, args):
    final_arg = [command]
    for field in args:
        if str(field.id) == "submit":
            continue
        final_arg.append("--"+str(field.id)+" "+str(field.data))
    final_arg=" ".join(final_arg)
    try:
        subprocess.Popen(final_arg, shell=True, stdout=subprocess.PIPE).stdout.read()
    except:
        pass
    return final_arg

def fetch_results_scFilterStats():
    df = pd.read_csv("./example_data/scFilterStats.txt", sep="\t", index_col=0)
    df['color']="#000000"
    pretty_labels={}
    for x in df.columns:
        pretty_labels[" ".join(x.split("_"))] = x
    source = ColumnDataSource(df)
    xlabel="Total sampled"
    ylabel="Wrong motif"
    fig = figure(plot_height=500, plot_width=500, tooltips=[("ID", "@index")])# level_0 refers to "index"
    fig.circle(x=pretty_labels[xlabel], y=pretty_labels[ylabel], source=source, size=8, color="color", line_color=None)
    fig.xaxis.axis_label = xlabel
    fig.yaxis.axis_label = ylabel
    script, div = components(fig)
    return script, div


def fetch_results_UMAP(gene=None):
    df = pd.read_csv("./example_data/normalized_counts.tsv", sep="\t", index_col=0)
    pretty_labels={}
    for x in df.columns:
        pretty_labels[" ".join(x.split("_"))] = x
    source = ColumnDataSource(df)
    ## map celltypes to a colormap
    xlabel="UMAP1"
    ylabel="UMAP2"
    fig = figure(plot_height=500, plot_width=700, tooltips=[("ID", "@cell")])# level_0 refers to "index"
    fig.add_layout(Legend(), 'right')
    if gene:
        fig.circle(x=pretty_labels[xlabel], y=pretty_labels[ylabel], source=source, size=8,
                   fill_color=linear_cmap(gene, palette=Blues[256][::-1], low=min(df[gene]), high=max(df[gene])),
                   line_color=None)
    else:
        fig.circle(x=pretty_labels[xlabel], y=pretty_labels[ylabel], source=source, size=8,
                   fill_color=factor_cmap('celltype', palette=Category20[20], factors=df.celltype.unique().tolist()),
                   legend_group='celltype', line_color=None)
    fig.xaxis.axis_label = xlabel
    fig.yaxis.axis_label = ylabel
    script, div = components(fig)
    return script, div
