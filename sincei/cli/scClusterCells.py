from __future__ import annotations

import logging
import warnings
from pathlib import Path
from typing import Annotated

import anndata as ad
import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import typer

from sincei.tools.ExponentialFamily import GLMFamily
from sincei.tools.GLMPCA import GLMPCA
from sincei.tools.TopicModels import TOPICMODEL

from ._common_args import (
    AVAILABLE_PROCESSORS,
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
    'Each dimensionality reduction is stored in ``obsm["X_<method>"]``. If the input '
    "AnnData object already contains a reduction of the requested method, the "
    "existing reduction is used instead of computing a new one, unless "
    "``--recomputeDimRed`` is given.\n\n"
    "``scClusterCells`` provides the following dimensionality reduction methods:\n"
    "* glmPCA: generalized PCA, with an exponential family distribution such as "
    "Poisson, Bernoulli, etc.\n"
    "* logPCA: Principal Component Analysis preceded by a logarithm transform.\n"
    "* LSA: Latent Semantic Analysis.\n"
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
                "[bold yellow]glmPCA[/bold yellow],"
                "[bold yellow]logPCA[/bold yellow], "
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
    out_file_trained_model: Annotated[
        str | None,
        typer.Option(
            "-om",
            "--outFileTrainedModel",
            rich_help_panel=_REDUCTION,
            help=(
                "The output file for the trained model. The saved model can be used "
                "later to embed/compare new cells to the existing cluster of cells. "
                "Only used with [bold yellow]LSA[/bold yellow] and "
                "[bold yellow]LDA[/bold yellow]."
            ),
        ),
    ] = None,
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
    n_passes: Annotated[
        int,
        typer.Option(
            "--nPasses",
            rich_help_panel=_LDA,
            help=(
                "Number of passes through the corpus for LDA model fitting. Only used "
                "with [bold yellow]LDA[/bold yellow]."
            ),
        ),
    ] = 5,
    n_iterations: Annotated[
        int,
        typer.Option(
            "--nIterations",
            rich_help_panel=_LDA,
            help=(
                "Maximum number of iterations through the corpus when inferring the "
                "topic distribution. Only used with [bold yellow]LDA[/bold yellow]."
            ),
        ),
    ] = 500,
    alpha: Annotated[
        float,
        typer.Option(
            "--alpha",
            rich_help_panel=_LDA,
            help=(
                "Prior to initialize cell-topic vectors. Only used with "
                "[bold yellow]LDA[/bold yellow]."
            ),
        ),
    ] = 50.0,
    eta: Annotated[
        float,
        typer.Option(
            "--eta",
            rich_help_panel=_LDA,
            help=(
                "Prior to initialize feature-topic vectors. Only used with "
                "[bold yellow]LDA[/bold yellow]."
            ),
        ),
    ] = 0.1,
    gamma_threshold: Annotated[
        float,
        typer.Option(
            "--gammaThreshold",
            rich_help_panel=_LDA,
            help=(
                "Minimum change in the value of the gamma parameters to continue "
                "iterating. Only used with [bold yellow]LDA[/bold yellow]."
            ),
        ),
    ] = 0.001,
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
    number_of_processors: Annotated[
        int, OTHER_OPTS["number_of_processors"]
    ] = AVAILABLE_PROCESSORS,
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
            recompute_dim_red=recompute_dim_red,
            out_file_trained_model=out_file_trained_model,
            n_passes=n_passes,
            n_iterations=n_iterations,
            alpha=alpha,
            eta=eta,
            gamma_threshold=gamma_threshold,
            glmpca_family=glmpca_family,
            out_file_umap=out_file_umap,
            cluster_resolution=cluster_resolution,
            plot_width=plot_width,
            plot_height=plot_height,
            plot_file_format=plot_file_format,
            dpi=dpi,
            number_of_processors=number_of_processors,
        )
    else:
        warnings.filterwarnings("ignore")

    adata = validate_anndata(ad.read_h5ad(input), input)
    reduction = f"X_{method.value}"

    model = None
    if reduction in adata.obsm and not recompute_dim_red:
        logging.info("Using the existing reduction in obsm[%r].", reduction)
    elif method is DimRed.logPCA:
        normalized = adata.copy()
        sc.pp.normalize_total(normalized, target_sum=1e4)
        sc.pp.log1p(normalized)
        sc.pp.pca(normalized, n_comps=n_prin_comps)
        adata.obsm[reduction] = normalized.obsm["X_pca"]
    elif method is DimRed.glmPCA:
        glmpca = GLMPCA(n_pc=n_prin_comps, family=glmpca_family)
        glmpca.fit(adata)
        assert glmpca.saturated_loadings_ is not None
        adata.obsm[reduction] = glmpca.saturated_loadings_.detach().numpy()
    else:
        topics = TOPICMODEL(
            adata,
            n_topics=n_prin_comps,
            binarize=binarize,
            n_passes=n_passes,
            n_workers=number_of_processors,
        )
        if method is DimRed.LSA:
            topics.runLSA()
            model = topics.lsi_model
        else:
            topics.runLDA(
                iterations=n_iterations,
                alpha=alpha,
                eta=eta,
                gamma_threshold=gamma_threshold,
            )
            model = topics.lda_model
        # The first component is dropped.
        adata.obsm[reduction] = topics.get_cell_topic().to_numpy()[:, 1:n_prin_comps]

    sc.pp.neighbors(adata, use_rep=reduction, n_neighbors=n_neighbors)
    sc.tl.leiden(adata, resolution=cluster_resolution)
    sc.tl.paga(adata)
    sc.pl.paga(adata, plot=False, threshold=0.1)
    sc.tl.umap(adata, min_dist=0.1, spread=5, init_pos="paga")

    adata.write_h5ad(out_file)

    if out_file_trained_model is not None:
        if model is not None:
            model.save(out_file_trained_model)
        elif method in {DimRed.LSA, DimRed.LDA}:
            logging.warning(
                "No model was trained, because obsm[%r] was reused. Give "
                "--recomputeDimRed to train and save a model.",
                reduction,
            )

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
