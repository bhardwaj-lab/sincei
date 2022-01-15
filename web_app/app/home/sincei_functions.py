# for bokeh plotting
from bokeh.models import ColumnDataSource, Div, Select, Slider, TextInput, Legend
from bokeh.io import curdoc
from bokeh.embed import components
from bokeh.plotting import figure, output_file, show
from bokeh.palettes import Category20
from bokeh.transform import factor_cmap
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
    df = pd.read_csv("/Users/vivek/programs/sincei/web_app/example_data/scFilterStats.txt", sep="\t", index_col=0)
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


def fetch_results_UMAP():
    df = pd.read_csv("/Users/vivek/programs/sincei/web_app/example_data/umap_annotation.tsv", sep="\t",
                     index_col=0,
                     names=['cell', 'celltype', 'UMAP1', 'UMAP2'])
    pretty_labels={}
    for x in df.columns:
        pretty_labels[" ".join(x.split("_"))] = x
    source = ColumnDataSource(df)
    ## map celltypes to a colormap
    xlabel="UMAP1"
    ylabel="UMAP2"
    fig = figure(plot_height=500, plot_width=800, tooltips=[("ID", "@cell")])# level_0 refers to "index"
    fig.add_layout(Legend(), 'right')
    fig.circle(x=pretty_labels[xlabel], y=pretty_labels[ylabel], source=source, size=8,
               fill_color=factor_cmap('celltype', palette=Category20[20], factors=df.celltype.unique().tolist()),
               legend_group='celltype',
               line_color=None)
    fig.xaxis.axis_label = xlabel
    fig.yaxis.axis_label = ylabel
    script, div = components(fig)
    return script, div
