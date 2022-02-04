# for bokeh plotting
from bokeh.models import ColumnDataSource, Div, Select, Slider, TextInput, Legend, ColorBar, LinearColorMapper
from bokeh.io import curdoc
from bokeh.embed import components
from bokeh.plotting import figure, output_file, show
from bokeh.palettes import Category20, Blues
from bokeh.transform import factor_cmap, linear_cmap
import pandas as pd
import numpy as np
import subprocess
import scanpy as sc
from deeptoolsintervals import GTF
import sys
import os
sys.path.append("../sincei")
from RegionQuery import *
WORKDIR="example_data"
import warnings


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

# load anndata object from default output dir
def load_anndata(name):
    """
    Load filtered/unfiltered anndata object
    """
    print("LOADING ANNDATA")
    ad = sc.read_loom(os.path.join(WORKDIR, "output", "loomfiles", name+'.loom'))
    ad.obs.set_index('obs_names', inplace=True)
    ad.var.set_index('var_names', inplace=True)
    ## if data is not binarized, normalize counts per cell
    if (np.sum(ad.X > 1) != 0):
        sc.pp.normalize_total(ad, target_sum=100000, exclude_highly_expressed=True)
    df = pd.DataFrame(ad.obsm['X_umap'])
    df.index = ad.obs_names
    df.columns = ['UMAP1', 'UMAP2']
    df['Cluster'] = ad.obs['cluster_lsi']#WIP:expeted colname for clusters
    df['celltype'] = ad.obs['hour']#WIP: expected colname for celltypes
    return ad, df

def get_gtf_olaps(ad, gtf_prefix="GRCz11"):
    """
    Check for overlap of ad.var with gtf regions
    """
    if all([len(x.split('_')) == 3 for x in ad.var_names]):
        gtf_file = os.path.join(WORKDIR, "annotations", gtf_prefix+".gtf")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gtf = GTF(gtf_file, transcript_id_designator='gene_name',keepExons=False)
        olap = get_gtf_adata_olaps(ad, gtf)
        return olap
    else:
        return None

def fetch_results_scFilterStats():
    """
    Load and plot the output of scFilterStats
    """
    df = pd.read_csv(os.path.join(WORKDIR, "output", "textfiles", "scFilterStats.txt"), sep="\t", index_col=0)
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


def fetch_results_UMAP(ad, df, gene=None, gene_bin_dict=None):
    """
    Plot UMAP with annotation, and if asked, gene/region signal
    """
    pretty_labels={}
    for x in df.columns:
        pretty_labels[" ".join(x.split("_"))] = x
    xlabel="UMAP1"
    ylabel="UMAP2"
    ## map celltypes to a colormap
    if gene:
        try:
            gene_idx = ad.var.index.get_loc(gene)
        except KeyError:# gene name is not in index
            if gene_bin_dict:
                bins=get_bins_by_gene(gene_bin_dict, gene, firstBin=False)
                gene_idx = [ad.var.index.get_loc(x) for x in bins]
            else:
                print("Gene not found in anndata.var. Gene->bin mapping unavailable")
        out_df = pd.DataFrame(np.sum(ad.X[:, gene_idx].todense(), axis=1))
        out_df.index = df.index
        out_df.rename({0:gene}, axis=1, inplace=True)
        # scale counts to range 0:100 across cells
        out_df=out_df.apply(lambda a: (100*(a - np.min(a))/np.ptp(a)), axis=0)

        df = df.join(out_df)
        source = ColumnDataSource(df)
        fig = figure(title="Activity of Feature: {}".format(gene),
                    plot_height=500, plot_width=500, tooltips=[("ID", "@celltype"),
                                                               ("Intensity", "@"+gene)])# level_0 refers to "index"
        fig.circle(x=pretty_labels[xlabel], y=pretty_labels[ylabel], source=source, size=8,
                   fill_color=linear_cmap(gene, palette=Blues[256][::-1], low=min(df[gene]), high=max(df[gene])),
                   line_color=None)
        ## add color bar
        color_mapper=LinearColorMapper(palette=Blues[256][::-1], low=0, high=100)
        color_bar = ColorBar(color_mapper=color_mapper, label_standoff=12)
        fig.add_layout(color_bar, 'right')
    else:
        source = ColumnDataSource(df)
        fig = figure(title="Cell types", plot_height=500, plot_width=700, tooltips=[("ID", "@celltype")])# level_0 refers to "index"
        fig.add_layout(Legend(), 'right')
        fig.circle(x=pretty_labels[xlabel], y=pretty_labels[ylabel], source=source, size=8,
                   fill_color=factor_cmap('celltype', palette=Category20[20],
                                          factors=df.celltype.unique().tolist()),
                   legend_group='celltype', line_color=None)
    fig.xaxis.axis_label = xlabel
    fig.yaxis.axis_label = ylabel
    script, div = components(fig)
    return script, div
