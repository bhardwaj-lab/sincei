Latent Dirichlet Allocation (LDA)
=================================

Latent Dirichlet Allocation (LDA) is a probabilistic topic model
:footcite:p:`blei2003latent` originally proposed for natural language processing.
In single-cell genomics, a cell is considered the analogue of a document and a
genomic feature the analogue of a word. LDA describes a cell as a mixture of
:math:`K` topics, and every topic as a distribution over the features. The
representation of a cell is its vector of topic proportions.

sincei fits LDA by collapsed Gibbs sampling :footcite:p:`griffiths2004finding` using
the `tomotopy <https://bab2min.github.io/tomotopy/>`_ library.

The generative model
--------------------

Let :math:`N` be the number of cells, :math:`P` the number of features, and :math:`K`
the number of topics to fit. LDA assumes that the data comes from this
process:

1. For each topic :math:`k`, draw a distribution over the features:
   :math:`\phi_k \sim \operatorname{Dirichlet}(\eta)`.
2. For each cell :math:`c`, draw the topic proportions:
   :math:`\vartheta_c \sim \operatorname{Dirichlet}(\alpha)`.
3. For each token :math:`t` of cell :math:`c`, draw a topic
   :math:`z_t \sim \operatorname{Categorical}(\vartheta_c)`, then draw a feature
   :math:`w_t \sim \operatorname{Categorical}(\phi_{z_t})`.

:math:`\alpha` controls the topic proportions of a cell. A small :math:`\alpha`
results in cell representations with few topics.

:math:`\eta` controls the features of a topic. A small :math:`\eta` will give topics
with few features.

From counts to tokens
---------------------

LDA requires integer counts. A count :math:`x_{ci}` of feature :math:`i` in cell
:math:`c` gives :math:`x_{ci}` tokens of :math:`i`. The total number of tokens in
the dataset is the sum of all counts :math:`T = \sum_{c,i} x_{ci}`.

If the data is binarized before LDA, the number of tokens is equivalent to the number
of non-zero entries of the matrix :math:`T = \mathrm{nnz}`.

To binarize or not to binarize
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Many scATAC-seq methods binarize counts, because a genomic region of a diploid cell has
at most two accessible copies. :footcite:t:`martens2024modeling` showed that this
may remove useful information. In four scATAC-seq datasets they found that:

- More than 65% of the non-zero counts were above 1.
- Fragment counts (not read counts) fit a Poisson distribution: their variance
  was approximately equal to their mean. Read counts of paired-end data are mostly
  even, because both reads of a fragment usually fall in the same peak, and their
  variance was approximately twice their mean.
- High counts were more frequent at promoters of highly expressed genes, and low
  counts at distal enhancers.
- A model of the fragment counts (a Poisson variational autoencoder) gave a better
  reconstruction of the data than the same model trained on binarized counts, and in
  most datasets it also resulted in better batch integration. The advantage increased
  with the number of fragments per cell.

This may or may not hold for other data modalities such as histone marks, but it
does point out that binarization of genomic counts may remove useful information.

One advantage of binarized counts is that they are computationally cheaper to model.
Collapsed Gibbs sampling gives each token its own topic assignment, so the time and
memory of the fit increase with the number of tokens :math:`T`
(see `Space and time complexity`_). With binarized counts, :math:`T` decreases to
:math:`\mathrm{nnz}` from the total number of counts.

Collapsed Gibbs sampling
------------------------

Let :math:`\mathbf{w} = (w_1, \dots, w_T)` be the features of all the tokens. The
posterior :math:`p(z, \vartheta, \phi \mid \mathbf{w})` has no closed form. Collapsed
Gibbs sampling integrates out :math:`\vartheta` and :math:`\phi`, and samples only
the topic assignment :math:`z_t` of each token. Let:

- :math:`n_{ck}`: the number of tokens of cell :math:`c` with topic :math:`k`,
- :math:`n_{ki}`: the number of tokens of feature :math:`i` with topic :math:`k`,
- :math:`n_k = \sum_i n_{ki}`: the number of tokens with topic :math:`k`,

and let the superscript :math:`-t` mean "without token :math:`t`". Each iteration
samples the topic of every token again from its conditional distribution:

.. math::

   p(z_t = k \mid z^{-t}, \mathbf{w}) \;\propto\;
   \big(n_{ck}^{-t} + \alpha_k\big)\,
   \frac{n_{k w_t}^{-t} + \eta}{n_k^{-t} + P\eta}

with :math:`c` the cell of token :math:`t`.

The first term favors the frequent topics in the cell, and the second term favors
the topics that already contain the feature.

After each sample, the counts :math:`n_{ck}`, :math:`n_{ki}` and :math:`n_k` change
to include the new topic of the token.

Estimates
~~~~~~~~~

Results come from the last sample in a run. The topic proportions of cell :math:`c` are:

.. math::

   \hat\vartheta_{ck} = \frac{n_{ck} + \alpha_k}{n_c + \sum_{k'} \alpha_{k'}}

with :math:`n_c` the number of tokens of the cell.

The feature distribution of topic :math:`k` is:

.. math::

   \hat\phi_{ki} = \frac{n_{ki} + \eta}{n_k + P\eta}

over all features, also the features without tokens.

Optimization of alpha
~~~~~~~~~~~~~~~~~~~~~

tomotopy optimizes :math:`\alpha` during the sampling. By default, every 10 iterations
from the start (``optim_interval=10``, ``burn_in=0``), it estimates one :math:`\alpha_k`
for each topic with Minka's fixed-point iteration :footcite:p:`minka2000estimating`.
Thus, the given :math:`\alpha` (:math:`1/K` by default) is only the start value of a
symmetric prior, and the final prior is asymmetric: frequent topics get a larger
:math:`\alpha_k`. :math:`\eta` stays fixed.

Parallel sampling
~~~~~~~~~~~~~~~~~

tomotopy samples the tokens on all available processors. The threads sample different
parts of the data at the same time, and then combine their counts. This is an
approximation of the sequential sampler. This means that for a fixed random seed, the
result depends on the number of threads used in the computation.

Space and time complexity
-------------------------

With :math:`N` cells, :math:`P` features, :math:`T` tokens, :math:`K` topics and
:math:`S` iterations:

- **Time**: :math:`O(S\,T\,K)`. Each iteration samples every token from :math:`K`
  topics. Runtime is linear with the number of tokens. Binarized counts reduce
  :math:`T` to :math:`\mathrm{nnz}`.
- **Memory**: :math:`O(T + (N + P)\,K)`, for one topic assignment per token and the
  tables :math:`n_{ck}` and :math:`n_{ki}`. sincei passes tokens to tomotopy as lists
  of strings, which needs an additional :math:`O(T)` memory while the model is built.

References
----------

.. footbibliography::
