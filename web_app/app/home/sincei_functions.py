# for bokeh plotting
from bokeh.models import ColumnDataSource, Div, Select, Slider, TextInput
from bokeh.io import curdoc
from bokeh.embed import components
from bokeh.plotting import figure, output_file, show
import pandas as pd

# execute a command using args passed from WTforms
def execute_command(command, args):
    if command == "scFilterStats":
        result = True
    return result

def fetch_results_scFilterStats():
    def getFilterStats(txtpath="/Users/vivek/programs/sincei/web_app/example_data/scFilterStats.txt"):
        res = pd.read_csv(txtpath, sep="\t", index_col=0)
        return res

    df = getFilterStats()
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
