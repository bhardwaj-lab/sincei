Generalized Principal Component Analysis (GLMPCA)
=================================================

Principal Component Analysis (PCA) finds the linear subspace that gives the smallest
squared reconstruction error for a given matrix. This is the maximum likelihood
solution for Gaussian data. However, single-cell data is usually not Gaussian.
For example, it can be:

- Counts (gene expression, open chromatin, histone modifications, etc.)
- Proportions (methylation ratios, allele fractions, etc.)
- Binary or categorical (mutations, presence/absence, etc.)

Generalized PCA (GLMPCA) replaces the Gaussian likelihood of PCA with the likelihood
of an exponential family distribution that fits the data.

sincei uses the saturated parameter formulation of GLMPCA
:footcite:p:`landgraf2020binary,mourragui2023percolate`. In this formulation, GLMPCA
projects the saturated parameters of the data onto a linear subspace.

Exponential family distributions
--------------------------------

Let :math:`X \in \mathbb{R}^{N \times G}` be the data matrix, with :math:`N` cells
in rows and :math:`G` features (bins, genes or other genomic regions) in columns. GLMPCA
models every entry :math:`x_{ci}` with an exponential family distribution:

.. math::

   p(x_{ci} \mid \theta_{ci}) = h(x_{ci})
   \exp\big( T(x_{ci})\,\eta(\theta_{ci}) - A(\theta_{ci}) \big)

with the parameter :math:`\theta_{ci}`, the sufficient statistic :math:`T`, the
natural parametrization :math:`\eta`, the log-partition function :math:`A` and the
base measure :math:`h`. The cost of a model :math:`\hat\Theta` is the negative
log-likelihood of the data:

.. math::

   \mathcal{L}(\hat\Theta; X) = - \sum_{c=1}^{N} \sum_{i=1}^{G}
   \log p(x_{ci} \mid \hat\theta_{ci})

For the Gaussian distribution with unit variance, :math:`\mathcal{L}` is the squared
error, and so Gaussian GLMPCA is equal to PCA.

Saturated parameters
--------------------

The saturated parameter :math:`\tilde\theta_{ci}` of an entry is the parameter that
gives the highest likelihood to that entry alone:

.. math::

   \tilde\theta_{ci} = g^{-1}(x_{ci})
   = \operatorname*{arg\,max}_{\theta} \; p(x_{ci} \mid \theta)

The saturated model :math:`\tilde\Theta = g^{-1}(X)` has one parameter for each
entry, so it fits the data perfectly. :math:`g^{-1}` is applied to every entry
separately.

The low-rank model
------------------

GLMPCA looks for a rank-:math:`K` approximation of the saturated parameters. It
centers them with a feature intercept :math:`\mu \in \mathbb{R}^G`, and projects
them onto the span of an orthonormal matrix of loadings
:math:`V \in \mathbb{R}^{G \times K}`:

.. math::

   \hat\Theta = \mathbf{1}\mu^\top + (\tilde\Theta - \mathbf{1}\mu^\top)\,V V^\top,
   \qquad V^\top V = I_K

The scores of the cells are the coordinates of this projection:

.. math::

   U = (\tilde\Theta - \mathbf{1}\mu^\top)\,V \in \mathbb{R}^{N \times K}

GLMPCA finds the :math:`V` and :math:`\mu` that minimize the cost:

.. math::

   \min_{V,\,\mu} \; \mathcal{L}\big(\mathbf{1}\mu^\top +
   (\tilde\Theta - \mathbf{1}\mu^\top)\,V V^\top;\; X\big)
   \quad \text{subject to} \quad V^\top V = I_K

The loadings :math:`V` are on the Stiefel manifold, the set of :math:`G \times K`
matrices with orthonormal columns.

This formulation has these advantages over the factorization
:math:`\hat\Theta = \mathbf{1}\mu^\top + U V^\top` with free :math:`U`:

- **Few parameters**: the model learns only :math:`V` and :math:`\mu`
  (:math:`G\,K + G` values). The number of parameters does not increase with the
  number of cells.
- **Exact out-of-sample projection**: the scores of a new cell are
  :math:`u = V^\top(g^{-1}(x) - \mu)`. This is a linear projection, and it does not
  need a new fit.
- **Stochastic optimization**: the cost is a sum over cells, and every cell depends
  only on :math:`V` and :math:`\mu`. Thus, mini-batches of cells give unbiased
  gradients, and GLMPCA scales to large datasets.
- **No derivations**: a new distribution needs only its likelihood and
  :math:`g^{-1}`. Automatic differentiation computes the gradients.

Truncation of the saturated parameters
--------------------------------------

For some distributions, the saturated parameter is infinite. For example, the
Poisson parameter :math:`\theta = \log \lambda` gives
:math:`\tilde\theta = \log 0 = -\infty` for a zero count. sincei clips the data
before it applies :math:`g^{-1}`, as in :footcite:t:`landgraf2020binary`. The clipping
value is a regularization parameter: a larger :math:`|\tilde\theta|` gives more
weight to the zeros.

The table shows the saturated parameter of each family, with the default clipping
values:

.. list-table::
   :header-rows: 1
   :widths: 18 32 50

   * - Family
     - Parameter :math:`\theta`
     - Saturated parameter :math:`\tilde\theta = g^{-1}(x)`
   * - ``gaussian``
     - mean, unit variance
     - :math:`x`
   * - ``poisson``
     - :math:`\log` rate
     - :math:`\log \max(x, 0.02)`
   * - ``negative_binomial``
     - :math:`\log` mean, inverse dispersion :math:`\kappa_i`
     - :math:`\log \max(x, 0.02)`
   * - ``bernoulli``
     - logit of the probability
     - :math:`\operatorname{clip}(\operatorname{logit} x, -4, 4)`
   * - ``gamma``
     - shape :math:`\theta + 1`, rate :math:`\nu_i`
     - :math:`\psi^{-1}\big(\log(\nu_i x)\big) - 1`
   * - ``lognormal``
     - mean of :math:`\log x`, standard deviation :math:`\nu_i`
     - :math:`\log \max(x, 10^{-5})`
   * - ``beta``
     - mean :math:`\pi \in (0, 1)`, concentration :math:`\nu_i`
     - the root :math:`\pi` of
       :math:`\operatorname{logit} x + \psi\big((1-\pi)\nu_i\big) - \psi(\pi\nu_i) = 0`
   * - ``sigmoid_beta``
     - logit of the mean, concentration :math:`\nu_i`
     - :math:`\operatorname{logit}` of the ``beta`` root

:math:`\psi` is the digamma function. ``sigmoid_beta`` is a ``beta`` distribution
with the mean :math:`\sigma(\theta)`, so :math:`\theta` has no bounds. With ``beta``,
the projection can move the mean out of :math:`(0, 1)`.

Nuisance parameters
~~~~~~~~~~~~~~~~~~~

Some families have one more parameter for each feature, shown with the index
:math:`i` in the table. GLMPCA estimates it once from the column statistics of the
data, before the fit, and keeps it fixed:

- ``negative_binomial``: method of moments,
  :math:`\kappa_i = \bar x_i^2 / (\sigma_i^2 - \bar x_i)` with the mean
  :math:`\bar x_i` and the variance :math:`\sigma_i^2` of the feature, clipped to
  :math:`[0.1, 500]`.
- ``beta`` and ``sigmoid_beta``: maximum likelihood with the fixed-point iteration of
  Minka, from a method of moments start.
- ``gamma``: maximum likelihood of the shape :math:`a_i` with Newton's method, and the
  rate :math:`\nu_i = a_i / \bar x_i`.
- ``lognormal``: the standard deviation of :math:`\log x`.

Sparse data
-----------

:math:`g^{-1}` is applied to every entry separately. Thus, all zeros of feature
:math:`i` have the same saturated parameter :math:`\theta^0_i = g^{-1}(0)`. sincei
writes the saturated parameters of sparse data as:

.. math::

   \tilde\Theta = \mathbf{1}\,(\theta^0)^\top + \Delta

where :math:`\Delta` is sparse and has the non-zeros of :math:`X`. The scores are then:

.. math::

   U = \mathbf{1}\,(\theta^0 - \mu)^\top V + \Delta V

This calculation is exact for all families. It costs :math:`O(\mathrm{nnz}\,K)`
and not :math:`O(N\,G\,K)`. GLMPCA never makes the full matrix dense: it reads
cells in chunks.

Sequencing depth correction
---------------------------

We can include a size factor to correct the effects of sequencing depth per cell.
For cell :math:`c`, the size factor is :math:`s_c = L_c / \operatorname{median}(L)`,
with :math:`L_c = \sum_i x_{ci}` the total counts of the cell. It can be applied in one
of two ways:

- ``"normalize"`` fits GLMPCA on :math:`x_{ci} / s_c`. It is available for
  ``poisson``, ``negative_binomial``, ``gamma`` and ``lognormal``.
- ``"offset"`` keeps the counts in the likelihood and adds :math:`\log s_c` to the
  parameters of the cell:

  .. math::

     \hat\theta_{ci} = \log s_c + \mu_i + u_c^\top v_i

  with :math:`u_c` the row of cell :math:`c` in :math:`U` and :math:`v_i` the row of
  feature :math:`i` in :math:`V`. The non-zero saturated parameters of the cell move by
  :math:`-\log s_c`, which is :math:`g^{-1}(x_{ci} / s_c)` for the log-link families. It is available for
  ``poisson``, ``negative_binomial`` and ``lognormal``.

Optimization
------------

**Initialization.** The intercept starts at the saturated parameter of the column
means, :math:`\mu_0 = g^{-1}(\bar{x})`, with :math:`\bar x_i` the mean of feature :math:`i`. This is the maximum likelihood intercept
when :math:`T(x) = x`. The loadings start with a randomized SVD of
:math:`\tilde\Theta - \mathbf{1}\mu_0^\top` on a random set of up to 8,192 cells
(``init="spectral"``).

**Riemannian Adam (default).** The loadings :math:`V` are a parameter on the Stiefel
manifold, and the optimizer is Riemannian Adam (geoopt). Every step uses a
mini-batch of cells (1,024 by default), and a retraction keeps :math:`V` orthonormal.
The learning rate of :math:`V` is :math:`0.5 / \sqrt{G}` (default) because the entries
of an orthonormal :math:`V` are of order :math:`1/\sqrt{G}`. When the mean cost per
entry does not decrease for some epochs, the learning rate decreases by ``gamma`` (0.5
by default). The fit stops at the fourth plateau or after 100 epochs. If the cost
becomes infinite or NaN, the fit starts again with a smaller learning rate.

**L-BFGS.** With ``optimizer="lbfgs"``, GLMPCA minimizes the full cost with L-BFGS
and a strong Wolfe line search. It writes :math:`V = \operatorname{qf}(Y)`, the QR
retraction of a free matrix :math:`Y`. The gradients are exact, and there is no
learning rate. Use it for up to approximately :math:`10^5` cells, or with
``fit_mode="subsample"``.

**Fit modes.** ``fit_mode`` selects the cells of each step: ``"minibatch"``
(default), ``"exact"`` (all cells in each step), or ``"subsample"`` (a fit on a
random set of cells, then a projection of all cells).

Whitening of the components
---------------------------

The cost depends only on the subspace :math:`V V^\top`. Thus, a rotation
:math:`V P` with an orthogonal :math:`P \in \mathbb{R}^{K \times K}` gives the same
cost. After the fit, sincei
selects the rotation that makes the scores uncorrelated, as in PCA. It computes the
eigendecomposition of the Gram matrix of the scores:

.. math::

   U^\top U = P\,\Lambda\,P^\top

and replaces :math:`V` with :math:`V P`, with the eigenvalues in decreasing order.
The first component then has the largest variance, and the singular values are
:math:`\sqrt{\Lambda}`.

Space and time complexity
-------------------------

With :math:`N` cells, :math:`G` features, :math:`\mathrm{nnz}` non-zero entries,
:math:`K` components, mini-batches of :math:`b` cells and :math:`E` epochs:

- **Scores**: :math:`O(\mathrm{nnz}\,K)` for sparse data. :math:`O(N\,G\,K)` for
  dense data.
- **Cost and gradients**: the likelihood includes every entry, zeros included. Thus,
  every epoch must compute :math:`\hat\Theta = \mathbf{1}\mu^\top + U V^\top` for all
  :math:`N\,G` entries: :math:`O(N\,G\,K)` time. With ``feature_subsample``, each step
  uses only a random fraction :math:`f` of the features: :math:`O(N\,f\,G\,K)`.
- **Retraction**: each step keeps :math:`V` on the Stiefel manifold with an SVD or
  QR step of :math:`O(G\,K^2)`. There are :math:`N / b` steps per epoch.
- **Total fit time**: :math:`O\big(E\,(N\,G\,K + (N/b)\,G\,K^2)\big)` with
  Riemannian Adam. The initialization adds a randomized SVD on up to 8,192 cells, and
  the whitening adds :math:`O(\mathrm{nnz}\,K + N\,K^2 + K^3)`.
- **Memory**: :math:`O(G\,K)` for :math:`V`, :math:`\mu` and the optimizer state,
  plus one dense chunk of cells (at most :math:`5 \times 10^7` values, approximately
  200 MB in float32). The saturated parameters are cached: :math:`O(\mathrm{nnz})`
  for sparse data, or :math:`O(N\,G)` for dense data if they fit in 4 GB. The full
  matrix never becomes dense.
- **Projection of new cells**: :math:`O(\mathrm{nnz}\,K)`, with no new fit.

The fit time increases with :math:`N\,G` and not with :math:`\mathrm{nnz}`. Thus,
GLMPCA is slower than LSA on sparse data, but it gives the exact likelihood of each
entry. The chunks of a pass can run in parallel on ``n_jobs`` threads, or on a GPU.

References
----------

.. footbibliography::
