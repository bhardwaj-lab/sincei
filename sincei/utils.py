from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Iterable

    from scipy import sparse

    # The matrix layouts `anndata.AnnData.X` holds for a count matrix.
    CountMatrix = (
        np.ndarray
        | sparse.csr_matrix
        | sparse.csc_matrix
        | sparse.csr_array
        | sparse.csc_array
    )


def gini(i: int, X: CountMatrix) -> float:
    r"""Computes the Gini coefficient for each row of a sparse matrix (Obs*Var).

    Parameters
    ----------
    i : int
        row index
    X : numpy array
        matrix

    Returns
    -------
    float
        Gini coefficient for the given row

    Examples
    --------
    >>> X = np.matrix([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
    >>> gini(0, X)
    0.0
    >>> gini(1, X)
    0.0
    >>> gini(2, X)
    0.0
    """

    # based on bottom eq:
    # http://www.statsdirect.com/help/generatedimages/equations/equation154.svg
    # from:
    # http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
    # All values are treated equally, arrays must be 1d:
    row = X[i, :]  # get all bins from i'th cell
    to_dense = getattr(row, "toarray", None)
    array = (row if to_dense is None else to_dense()).flatten()
    array = array[array.nonzero()]

    if array.shape[0] <= 1:
        return np.nan
    else:
        array = np.sort(array)
        # Index per array element:
        index = np.arange(1, array.shape[0] + 1)
        # Number of array elements:
        n = array.shape[0]
        # Gini coefficient:
        return (np.sum((2 * index - n - 1) * array)) / (n * np.sum(array))


def chromosome_to_numeric(chroms: Iterable[str]) -> dict[str, int]:
    """Map chromosome labels to a numeric sort order.

    Only the unique labels in ``chroms`` are considered. The ordering is:
    numbered (autosomal) chromosomes first, keeping their own number
    regardless of a "chr"/"CHR" prefix (e.g. "chr7" -> 7); then the sex
    chromosomes (X, Y) immediately after the highest-numbered autosome;
    then the mitochondrial chromosome (MT/M); and finally any remaining
    contigs/scaffolds, which get increasing values in sorted order.

    Returns a dict mapping each original label to its numeric order.
    """
    autosomes: dict[str, int] = {}  # label -> chromosome number
    sex: dict[str, str] = {}  # label -> "X"/"Y"
    mito: list[str] = []  # mitochondrial labels
    other: list[str] = []  # contigs / scaffolds

    for chrom in set(chroms):
        # Drop a leading "chr"/"CHR" prefix, case-insensitively.
        name = chrom.upper()
        if name.startswith("CHR"):
            name = name[3:]
        try:
            autosomes[chrom] = int(name)
        except ValueError:
            if name in {"X", "Y"}:
                sex[chrom] = name
            elif name in {"MT", "M"}:
                mito.append(chrom)
            else:
                other.append(chrom)

    mapping = dict(autosomes)
    counter = max(autosomes.values(), default=0)

    # Sex chromosomes ("X", "Y")
    for target in ("X", "Y"):
        for chrom in sorted(label for label, s in sex.items() if s == target):
            counter += 1
            mapping[chrom] = counter
    # Mitochondrial chromosomes ("MT", "M")
    for chrom in sorted(mito):
        counter += 1
        mapping[chrom] = counter
    # Remaining contigs/scaffolds get increasing values
    for chrom in sorted(other):
        counter += 1
        mapping[chrom] = counter

    return mapping
