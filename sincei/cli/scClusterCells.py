from __future__ import annotations

import logging
import sys
import warnings
from pathlib import Path
from typing import Annotated

import anndata as ad
import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import typer

from sincei.tools.ExponentialFamily import GLMFamily
from sincei.tools.GLMPCA import GLMPCA
from sincei.tools.TopicModels import LDA, LSA
from sincei.utils import DEPTH_CORRELATION_CUTOFF, depth_components

from ._common_args import (
    CM_PER_INCH,
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    PLOT_OPTS,
    DimRed,
    PlotFileFormat,
    configure_logging,
    log_parameters,
    preprocess_args,
)
from ._parsers import validate_anndata

DESCRIPTION = (
    "Cluster cells from a cell-by-feature matrix.\n\n"
    "``scClusterCells`` clusters cells based on the input count matrix (output of "
    "``scCountReads``) and performs dimensionality reduction, community detection and "
    "2D projection (UMAP) of the cells. The result is an updated h5ad object, and "
    "(optionally) a plot file and a .tsv file with UMAP coordinates and the "
    "corresponding cluster id for each barcode.\n\n"
    'Each dimensionality reduction is stored in ``obsm["X_<method>"]``, and the '
    'region loadings of LSA and LDA in ``varm["X_<method>"]``, with the method name '
    'in lowercase (for example ``obsm["X_glmpca"]``). If the input AnnData object '
    "already contains a reduction of the requested method, the existing reduction"
    "is used instead of computing a new one, unless ``--recomputeDimRed`` is given. "
    "With ``--logNormalize``, the counts of every cell are scaled to a total of 10,000 "
    "and log1p-transformed before the reduction. The counts in the output file are not "
    "changed.\n\n"
    "``scClusterCells`` provides the following dimensionality reduction methods:\n"
    "* glmPCA: generalized PCA, with an exponential family distribution such as "
    "Poisson, Bernoulli, etc.\n"
    "* PCA: Principal Component Analysis.\n"
    "* LSA: Latent Semantic Analysis. One extra component is fitted, and the "
    "component that correlates the most with the cell depth is dropped if its "
    f"|Spearman rho| is above {DEPTH_CORRELATION_CUTOFF}, else the extra component.\n"
    "* LDA: Latent Dirichlet Allocation."
)


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help=DESCRIPTION,
    context_settings={"help_option_names": []},
)

_REDUCTION = "Dimensionality reduction options"
_LDA = "LDA options"
_GLMPCA = "glmPCA options"
_CLUSTERING = "Clustering options"


@app.callback(invoke_without_command=True)
def main(
    # Input / Output options
    input: Annotated[str, INPUT_OUTPUT_OPTS["h5ad_file"]],
    out_file: Annotated[str, INPUT_OUTPUT_OPTS["out_file"]],
    # Dimensionality reduction options
    method: Annotated[
        DimRed,
        typer.Option(
            "-m",
            "--method",
            metavar="METHOD",
            rich_help_panel=_REDUCTION,
            help=(
                "The dimensionality reduction method to use before clustering cells."
                "\n\n"
                "One of: "
                "[bold yellow]glmPCA[/bold yellow], "
                "[bold yellow]PCA[/bold yellow], "
                "[bold yellow]LSA[/bold yellow], "
                "[bold yellow]LDA[/bold yellow]. "
            ),
        ),
    ] = DimRed.LSA,
    n_prin_comps: Annotated[
        int,
        typer.Option(
            "-n",
            "--nPrinComps",
            rich_help_panel=_REDUCTION,
            help=(
                "Number of principal components or topics to reduce the dimensionality "
                "to. Use a higher number for samples with more expected heterogeneity."
            ),
        ),
    ] = 20,
    n_neighbors: Annotated[
        int,
        typer.Option(
            "-nk",
            "--nNeighbors",
            rich_help_panel=_REDUCTION,
            help=(
                "Number of nearest neighbours to consider for clustering and UMAP. "
                "Choose this considering the total number of cells and the expected "
                "number of clusters; smaller numbers lead to more fragmented clusters."
            ),
        ),
    ] = 30,
    binarize: Annotated[
        bool,
        typer.Option(
            "--binarize",
            rich_help_panel=_REDUCTION,
            help=(
                "Binarize the counts per region before dimensionality reduction. Only "
                "used with [bold yellow]LSA[/bold yellow] and "
                "[bold yellow]LDA[/bold yellow]."
            ),
        ),
    ] = False,
    log_normalize: Annotated[
        bool,
        typer.Option(
            "--logNormalize",
            rich_help_panel=_REDUCTION,
            help=(
                "Depth-normalize the counts of every cell to a total of 10,000 and "
                "then log1p transform before dimensionality reduction. Ignored with "
                "[bold yellow]LSA[/bold yellow] and [bold yellow]LDA[/bold yellow], "
                "which take the counts."
            ),
        ),
    ] = False,
    recompute_dim_red: Annotated[
        bool,
        typer.Option(
            "--recomputeDimRed",
            rich_help_panel=_REDUCTION,
            help=(
                "Recompute the dimensionality reduction even if a precomputed version "
                "exists."
            ),
        ),
    ] = False,
    # glmPCA options
    glmpca_family: Annotated[
        GLMFamily,
        typer.Option(
            "-gf",
            "--glmPCAfamily",
            metavar="FAMILY",
            rich_help_panel=_GLMPCA,
            help=(
                "The choice of exponential family distribution to use for the glmPCA "
                "method. Only used with [bold yellow]glmPCA[/bold yellow].\n\n"
                "One of: [bold yellow]gaussian[/bold yellow], "
                "[bold yellow]poisson[/bold yellow], "
                "[bold yellow]bernoulli[/bold yellow], "
                "[bold yellow]beta[/bold yellow], "
                "[bold yellow]gamma[/bold yellow], "
                "[bold yellow]lognormal[/bold yellow], "
                "[bold yellow]sigmoid_beta[/bold yellow]."
            ),
        ),
    ] = GLMFamily.poisson,
    # LDA options
    n_iterations: Annotated[
        int,
        typer.Option(
            "--nIterations",
            rich_help_panel=_LDA,
            help=(
                "Number of Gibbs sampling iterations over all the counts. Only used "
                "with [bold yellow]LDA[/bold yellow]."
            ),
        ),
    ] = 1000,
    alpha: Annotated[
        float | None,
        typer.Option(
            "--alpha",
            rich_help_panel=_LDA,
            help=(
                "Prior of the topic proportions of a cell, the same for every topic. "
                "If not given, 1/nPrinComps is used. Only used with "
                "[bold yellow]LDA[/bold yellow]."
            ),
        ),
    ] = None,
    eta: Annotated[
        float,
        typer.Option(
            "--eta",
            rich_help_panel=_LDA,
            help=(
                "Prior of the region weights of a topic. Only used with "
                "[bold yellow]LDA[/bold yellow]."
            ),
        ),
    ] = 0.1,
    # Clustering options
    out_file_umap: Annotated[
        str | None,
        typer.Option(
            "-op",
            "--outFileUMAP",
            rich_help_panel=_CLUSTERING,
            help=(
                "The output plot file (for UMAP). If specified, a 4-column .tsv file "
                "with the same prefix is also created with the cell IDs, raw UMAP "
                "coordinates (UMAP1 and UMAP2) and Leiden cluster number."
            ),
        ),
    ] = None,
    cluster_resolution: Annotated[
        float,
        typer.Option(
            "-cr",
            "--clusterResolution",
            rich_help_panel=_CLUSTERING,
            help=(
                "Resolution parameter for Leiden clustering. Values lower than 1.0 "
                "result in fewer clusters, while higher values lead to splitting of "
                "clusters. In most cases the optimum is between 0.8 and 1.2."
            ),
        ),
    ] = 1.0,
    # Plot options
    plot_width: Annotated[float, PLOT_OPTS["plot_width"]] = 25.0,
    plot_height: Annotated[float, PLOT_OPTS["plot_height"]] = 25.0,
    plot_file_format: Annotated[
        PlotFileFormat, PLOT_OPTS["plot_file_format"]
    ] = PlotFileFormat.png,
    dpi: Annotated[int, PLOT_OPTS["dpi"]] = 300,
    # Other options
    verbose: Annotated[bool, OTHER_OPTS["verbose"]] = False,
    help: Annotated[bool, OTHER_OPTS["help"]] = False,
) -> int:
    if verbose:
        log_parameters(
            input=input,
            out_file=out_file,
            method=method,
            n_prin_comps=n_prin_comps,
            n_neighbors=n_neighbors,
            binarize=binarize,
            log_normalize=log_normalize,
            recompute_dim_red=recompute_dim_red,
            n_iterations=n_iterations,
            alpha=alpha,
            eta=eta,
            glmpca_family=glmpca_family,
            out_file_umap=out_file_umap,
            cluster_resolution=cluster_resolution,
            plot_width=plot_width,
            plot_height=plot_height,
            plot_file_format=plot_file_format,
            dpi=dpi,
        )
    else:
        warnings.filterwarnings("ignore")

    try:
        method = DimRed(method)
    except ValueError:
        choices = ", ".join(member.value for member in DimRed)
        sys.stderr.write(f"{method!r} is not a valid method. Use one of {choices}.\n")
        raise typer.Exit(code=1) from None
    if log_normalize and method in {DimRed.LSA, DimRed.LDA}:
        logging.warning(
            "--logNormalize is ignored with %s, which takes raw counts.", method.value
        )
        log_normalize = False

    adata = validate_anndata(ad.read_h5ad(input), input)
    reduction = f"X_{method.value.lower()}"

    if reduction in adata.obsm and not recompute_dim_red:
        logging.info("Using the existing reduction in obsm[%r].", reduction)
    else:
        matrix = adata
        if log_normalize:
            matrix = adata.copy()
            sc.pp.normalize_total(matrix, target_sum=1e4)
            sc.pp.log1p(matrix)

        if method is DimRed.PCA:
            sc.pp.pca(matrix, n_comps=n_prin_comps)
            adata.obsm[reduction] = matrix.obsm["X_pca"]
            adata.varm["PCs"] = matrix.varm["PCs"]
            adata.uns["pca"] = matrix.uns["pca"]
        elif method is DimRed.glmPCA:
            glmpca = GLMPCA(n_pc=n_prin_comps, family=glmpca_family)
            glmpca.fit(matrix)
            assert glmpca.saturated_loadings_ is not None
            adata.obsm[reduction] = glmpca.saturated_loadings_.detach().numpy()
        elif method is DimRed.LSA:
            lsa = LSA(n_pc=n_prin_comps + 1, binarize=binarize)
            scores = lsa.fit_transform(matrix)
            assert lsa.components_ is not None
            assert lsa.singular_values_ is not None
            depth = np.asarray(adata.X.sum(axis=1), dtype=np.float32).ravel()
            keep = depth_components(scores, depth)
            adata.obsm[reduction] = scores[:, keep]
            adata.varm[reduction] = lsa.components_[:, keep]
            adata.var["lsa_idf"] = lsa.idf_
            adata.uns["lsa"] = {
                "pivot": lsa.pivot_,
                "singular_values": lsa.singular_values_[keep],
                "binarize": binarize,
            }
        elif method is DimRed.LDA:
            lda = LDA(
                n_topics=n_prin_comps,
                alpha=alpha,
                eta=eta,
                iterations=n_iterations,
                binarize=binarize,
            )
            adata.obsm[reduction] = lda.fit_transform(matrix)
            assert lda.topics_ is not None
            adata.varm[reduction] = lda.topics_.T
            adata.uns["lda"] = {
                "alpha": 1.0 / n_prin_comps if alpha is None else alpha,
                "eta": eta,
                "binarize": binarize,
            }

    sc.pp.neighbors(adata, use_rep=reduction, n_neighbors=n_neighbors)
    sc.tl.leiden(adata, resolution=cluster_resolution)
    sc.tl.paga(adata)
    sc.pl.paga(adata, plot=False, threshold=0.1)
    sc.tl.umap(adata, min_dist=0.1, spread=5, init_pos="paga")

    adata.write_h5ad(out_file)

    if out_file_umap is not None:
        figure, axes = plt.subplots(
            figsize=(plot_width / CM_PER_INCH, plot_height / CM_PER_INCH)
        )
        sc.pl.umap(adata, color="leiden", legend_loc="on data", ax=axes, show=False)
        figure.savefig(out_file_umap, dpi=dpi, format=plot_file_format.value)
        plt.close(figure)

        umap = pd.DataFrame(
            adata.obsm["X_umap"],
            index=adata.obs_names,
            columns=pd.Index(["UMAP1", "UMAP2"]),
        )
        umap["cluster"] = adata.obs["leiden"]
        umap.to_csv(
            Path(out_file_umap).with_suffix(".tsv"), sep="\t", index_label="Cell_ID"
        )

    return 0


def cli() -> None:
    configure_logging()
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
