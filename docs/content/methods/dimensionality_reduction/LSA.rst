Latent Semantic Analysis (LSA)
==============================

Sincei's LSA applies the Term Frequency-Inverse Document Frequency (TF-IDF) transform
on a matrix followed by truncated Singular Value Decomposition (SVD) to obtain a
low-dimensional representation of the data.

The TF-IDF used is the ``lfu`` formulation. In
`SMART <https://en.wikipedia.org/wiki/SMART_Information_Retrieval_System>`_ notation,
``l`` is a logarithmic term frequency, ``f`` is the inverse document frequency, and
``u`` is the pivoted unique normalization.
Here, a cell is a document and a feature is a term.

Let :math:`X \in \mathbb{R}^{N \times P}` be the input matrix, :math:`x_{ij}` be the
count of feature :math:`j` in cell :math:`i`, :math:`N` the total number of cells,
:math:`P` the total number of features, and :math:`\mathrm{n}_j` the number of cells with
:math:`x_{ij} > 0` for feature :math:`j`.

The term frequency (``l``) is:

.. math::

   \mathrm{tf}_{ij} = 1 + \log_2 x_{ij} \quad \text{for } x_{ij} > 0

The inverse document frequency (``f``) is:

.. math::

   \mathrm{idf}_j = \log_2 \frac{N}{\mathrm{n}_j}

A feature in all cells has :math:`\mathrm{idf}_j = 0`.

And the pivoted unique normalization (``u``) divides the weights of a cell by:

.. math::

   d_i = (1 - \beta)\,\bar{\mathrm{n}} + \beta\,|\mathcal{F}_i|

with the slope :math:`\beta = 0.25`, the pivot
:math:`\bar{\mathrm{n}} = \frac{1}{N}\sum_j \mathrm{n}_j = \mathrm{nnz} / N`
(the mean number of features with counts in a cell), and
:math:`\mathcal{F}_i = \{\, j : x_{ij} > 0,\; \mathrm{idf}_j > 0 \,\}` the features of
cell :math:`i` with a count above zero and an idf above zero.

The TF-IDF weight matrix :math:`W \in \mathbb{R}^{N \times P}` is then:

.. math::

   w_{ij} =
   \begin{cases}
   \dfrac{\mathrm{tf}_{ij}\,\mathrm{idf}_j}{d_i} & \text{if } x_{ij} > 0 \\
   0 & \text{if } x_{ij} = 0
   \end{cases}

The zeros of the input matrix stay at zero, so :math:`W` stays sparse.

LSA then computes the truncated SVD :math:`W \approx U \Sigma V^\top` of the weight
matrix, with the top :math:`K` components in :math:`V \in \mathbb{R}^{P \times K}`
(features by components).

The scores of the cells are :math:`U \Sigma = W V \in \mathbb{R}^{N \times K}`, and
the score of cell :math:`i` is :math:`w_i V`, with :math:`w_i` the row of the cell
in :math:`W`.

Effect of cell depth
--------------------

The first component of LSA tends to capture cell depth (the total number of sequencing
reads in each cell). To deal with this, we compute one extra component. If the component
that correlates the most with the cell depth has an absolute Spearman correlation above
a cutoff (0.75), LSA drops that component. Otherwise, the extra component is dropped.

SVD solver
----------

Truncated SVD computes the top :math:`K` singular vectors of :math:`W`, with :math:`K`
much smaller than :math:`N` and :math:`P`.

Sincei's LSA calls ``scipy.sparse.linalg.svds``, which uses the ARPACK library
:footcite:p:`lehoucq1998arpack`. It uses :math:`W` only through sparse matrix products,
so :math:`W` never becomes dense. It finds the top eigenvectors of the Gram matrix
of :math:`W` (:math:`W^\top W` or :math:`W W^\top`, whichever is smaller) with the
implicitly restarted Lanczos method:

1. Build a Krylov subspace :math:`\{y, Hy, H^2y, \dots\}` of the Gram matrix :math:`H`.
   Each step needs one product with :math:`W` and another with :math:`W^\top`.
2. Keep a small set of Lanczos vectors (:math:`\max(2K + 1, 20)`), and restart to
   remove the directions of unwanted eigenvalues.
3. Stop when the top :math:`K` eigenpairs converge to machine precision.

The singular values are the square roots of the eigenvalues, and the singular vectors
come from one more product with :math:`W`.

ARPACK gives the optimal rank-:math:`K` SVD approximation, up to numerical precision.
The number of Lanczos steps increases automatically when singular values are close
to each other.

Space and time complexity
-------------------------

With :math:`N` cells, :math:`P` features, :math:`\mathrm{nnz}` the number of non-zero
entries in the matrix and :math:`K` components:

- **TF-IDF**

  - Time: :math:`O(\mathrm{nnz})`.
  - Space: :math:`O(\mathrm{nnz})`. :math:`W` has the same sparsity as the input matrix.

- **ARPACK**

  - Time: :math:`O\big(t\,(\mathrm{nnz} + \min(N, P)\,K) + t_r\,\min(N, P)\,K^2\big)`,
    with :math:`t` Lanczos steps and :math:`t_r` restarts. The number of steps
    depends on the gaps between the singular values.
  - Space: :math:`O(\mathrm{nnz} + (N + P)\,K)`.

- **Scores and depth check**

  - Time: :math:`O(\mathrm{nnz}\,K)` for :math:`W V`, and :math:`O(K\,N \log N)` for
    the Spearman correlations.
  - Space: :math:`O(N\,K)` for the scores.

References
----------

.. footbibliography::
