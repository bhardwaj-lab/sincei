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

Let :math:`X \in \mathbb{R}^{N \times G}` be the input matrix, :math:`x_{ci}` be the
count of feature :math:`i` in cell :math:`c`, :math:`N` the total number of cells,
:math:`G` the total number of features, and :math:`\mathrm{n}_i` the number of cells with
:math:`x_{ci} > 0` for feature :math:`i`.

The term frequency (``l``) is:

.. math::

   \mathrm{tf}_{ci} = 1 + \log_2 x_{ci} \quad \text{for } x_{ci} > 0

The inverse document frequency (``f``) is:

.. math::

   \mathrm{idf}_i = \log_2 \frac{N}{\mathrm{n}_i}

A feature in all cells has :math:`\mathrm{idf}_i = 0`.

And the pivoted unique normalization (``u``) divides the weights of a cell by:

.. math::

   d_c = (1 - \beta)\,\bar{\mathrm{n}} + \beta\,|\mathcal{F}_c|

with the slope :math:`\beta = 0.25`, the pivot
:math:`\bar{\mathrm{n}} = \frac{1}{N}\sum_i \mathrm{n}_i = \mathrm{nnz} / N`
(the mean number of features with counts in a cell), and
:math:`\mathcal{F}_c = \{\, i : x_{ci} > 0,\; \mathrm{idf}_i > 0 \,\}` the features of
cell :math:`c` with a count above zero and an idf above zero.

The TF-IDF weight matrix :math:`W \in \mathbb{R}^{N \times G}` is then:

.. math::

   w_{ci} =
   \begin{cases}
   \dfrac{\mathrm{tf}_{ci}\,\mathrm{idf}_i}{d_c} & \text{if } x_{ci} > 0 \\
   0 & \text{if } x_{ci} = 0
   \end{cases}

The zeros of the input matrix stay at zero, so :math:`W` stays sparse.

LSA then computes the truncated SVD :math:`W \approx U \Sigma V^\top` of the weight
matrix, with the top :math:`K` components in :math:`V \in \mathbb{R}^{G \times K}`
(features by components).

The scores of the cells are :math:`U \Sigma = W V \in \mathbb{R}^{N \times K}`, and
the score of cell :math:`c` is :math:`w_c V`, with :math:`w_c` the row of the cell
in :math:`W`.

Effect of cell depth
--------------------

The first component of LSA tends to capture cell depth (the total number of sequencing
reads in each cell). To deal with this, we compute one extra component. If the component
that correlates the most with the cell depth has an absolute Spearman correlation above
a cutoff (0.75), LSA drops that component. Otherwise, the extra component is dropped.

SVD solvers
-----------

Truncated SVD computes the top :math:`K` singular vectors of :math:`W`, with :math:`K`
much smaller than :math:`N` and :math:`G`.

Sincei's LSA offers two solvers: ARPACK (default) and randomized SVD using torch. Both
solvers use :math:`W` only through sparse matrix products, so :math:`W` never becomes
dense.

Use ``"arpack"`` for exact results and ``"torch"`` for very large matrices, or when
a GPU is available.

ARPACK
~~~~~~

``solver="arpack"`` calls ``scipy.sparse.linalg.svds``, which uses the ARPACK library
:footcite:p:`lehoucq1998arpack`. It finds the top eigenvectors of the Gram matrix
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
The number of Lanczos steps increases auttomatically when singular values are close
to each other.

Randomized SVD in torch
~~~~~~~~~~~~~~~~~~~~~~~

``solver="torch"`` computes a randomized SVD :footcite:p:`halko2009finding` using
torch :footcite:p:`paszke2019pytorch`, on the CPU or on a GPU. With
:math:`\ell = K + o` columns where :math:`o` is the oversample parameter (equal to
:math:`K` by default) and :math:`q` is the number of power iterations (8 by default):

1. Sample a random Gaussian matrix :math:`\Omega \in \mathbb{R}^{G \times \ell}`.
2. Compute an orthonormal basis :math:`P \in \mathbb{R}^{N \times \ell}` of the range
   of :math:`W`: :math:`P = \operatorname{orth}(W \Omega)`.
3. Repeat :math:`q` times (subspace iteration):

   .. math::

      P' = \operatorname{orth}(W^\top P), \qquad P = \operatorname{orth}(W P')

   Each iteration multiplies the weight of the singular value :math:`\sigma_k` by
   :math:`\sigma_k^2`, so the top components separate from the rest. The
   orthonormalization after each product keeps the basis numerically stable (see
   below).
4. Compute the exact SVD of the small matrix
   :math:`B = P^\top W = \tilde U \Sigma V^\top`.
5. Keep the top :math:`K` right singular vectors :math:`V`, and compute the scores
   :math:`W V`.

The result is an approximation. Its error decreases quickly with :math:`q` and
:math:`o`, and is small when the singular values decrease quickly. The work is in
large sparse-dense matrix products and in QR steps, which are fast on a GPU. A
fixed seed makes the result reproducible on the same device.

The orthonormalization step
^^^^^^^^^^^^^^^^^^^^^^^^^^^

:math:`\operatorname{orth}(A)` uses the QR decomposition :math:`A = QR`, with
:math:`Q` a matrix with orthonormal columns (the same size as :math:`A`) and
:math:`R` an upper triangular :math:`\ell \times \ell` matrix. It keeps :math:`Q`
and discards :math:`R`. :math:`R` is invertible when :math:`A` has full rank, so
:math:`Q` and :math:`A` span the same space.

Without this step, the products with :math:`W` and :math:`W^\top` make all the columns
of the basis converge to the top singular vector, because :math:`\sigma_1` grows the
fastest. This keeps each column pointing in a different direction and the values in
a stable range.

Space and time complexity
-------------------------

With :math:`N` cells, :math:`G` features, :math:`\mathrm{nnz}` the number of non-zero
entries in the matrix and :math:`K` components:

- **TF-IDF**

  - Time: :math:`O(\mathrm{nnz})`.
  - Space: :math:`O(\mathrm{nnz})`. :math:`W` has the same sparsity as the input matrix.

- **ARPACK**

  - Time: :math:`O\big(t\,(\mathrm{nnz} + \min(N, G)\,K) + t_r\,\min(N, G)\,K^2\big)`,
    with :math:`t` Lanczos steps and :math:`t_r` restarts. The number of steps
    depends on the gaps between the singular values.
  - Space: :math:`O(\mathrm{nnz} + (N + G)\,K)`.

- **Randomized SVD** (see the breakdown below)

  - Time: :math:`O\big((q + 1)\,(\mathrm{nnz}\,\ell + (N + G)\,\ell^2)\big)`, with
    :math:`\ell = K + o` columns and :math:`q` power iterations.
  - Space: :math:`O(\mathrm{nnz} + (N + G)\,\ell)`, with two copies of :math:`W`.

- **Scores and depth check**

  - Time: :math:`O(\mathrm{nnz}\,K)` for :math:`W V`, and :math:`O(K\,N \log N)` for
    the Spearman correlations.
  - Space: :math:`O(N\,K)` for the scores.

Randomized SVD in detail
~~~~~~~~~~~~~~~~~~~~~~~~

The randomized SVD works with blocks of :math:`\ell = K + o` columns. With the
defaults (:math:`o = K`, :math:`q = 8`), :math:`\ell = 2K`.

A product of the sparse :math:`W` (or :math:`W^\top`) with a dense matrix of
:math:`\ell` columns costs :math:`O(\mathrm{nnz}\,\ell)`: each non-zero of :math:`W`
is used once for each column.

A QR step of a tall matrix with :math:`r` rows and :math:`\ell` columns costs 
:math:`O(r\,\ell^2)`.

The table shows the cost of each step of the algorithm:

.. list-table::
   :header-rows: 1
   :widths: 40 15 45

   * - Step
     - Count
     - Cost of one
   * - :math:`W \Omega`, then QR (:math:`N \times \ell`)
     - 1
     - :math:`O(\mathrm{nnz}\,\ell + N\,\ell^2)`
   * - :math:`W^\top P`, then QR (:math:`G \times \ell`)
     - :math:`q`
     - :math:`O(\mathrm{nnz}\,\ell + G\,\ell^2)`
   * - :math:`W P'`, then QR (:math:`N \times \ell`)
     - :math:`q`
     - :math:`O(\mathrm{nnz}\,\ell + N\,\ell^2)`
   * - :math:`B = (W^\top P)^\top`, then its SVD (:math:`\ell \times G`)
     - 1
     - :math:`O(\mathrm{nnz}\,\ell + G\,\ell^2)`
   * - Scores :math:`W V` (:math:`K` columns)
     - 1
     - :math:`O(\mathrm{nnz}\,K)`

In total, there are :math:`2q + 2` sparse products with :math:`\ell` columns, and
:math:`2q + 1` QR steps. The time is thus:

.. math::

   O\big((q + 1)\,(\mathrm{nnz}\,\ell + (N + G)\,\ell^2)\big)

The sparse products usually take most of the time, because :math:`\mathrm{nnz}` is
much larger than :math:`(N + G)\,\ell`.

The memory is :math:`O(\mathrm{nnz} + (N + G)\,\ell)`:

- :math:`O(\mathrm{nnz})` for :math:`W` in two copies: one for :math:`W` and one for
  :math:`W^\top`. Both are in compressed sparse row format, so both products read
  the rows of a matrix, which is fast.
- :math:`O(N\,\ell)` for :math:`P`, and :math:`O(G\,\ell)` for :math:`\Omega`,
  :math:`P'` and :math:`B`.

On a GPU, all of these matrices must fit in the memory of the GPU.

Compared with ARPACK, the randomized SVD does a fixed number of products, and each
product uses :math:`\ell` columns at the same time. ARPACK computes more products,
each with one vector. Thus, the randomized SVD makes better use of the hardware,
especially on a GPU.

References
----------

.. footbibliography::
