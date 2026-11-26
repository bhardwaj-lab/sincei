Variable Chromatin Region Finder
================================

:ref:`scFindVCRs` segments the chromosomes in an input AnnData into variable chromatin
regions (VCRs): genomic regions with signal in the same cells. In other words,
in aggregates co-occurring genomic bins into larger related regions.

It takes a count matrix :math:`X \in \mathbb{R}^{N \times P}` of :math:`N` cells
and :math:`P` genomic bins of equal size (from :ref:`scCountReads bins <scCountReads>`)
as input, and processes it in two stages: a prefilter that splits each chromosome into
signal blocks, and segmentation of each block into regions.

The prefilter keeps bins which have nonzero counts in at least a percentage of the total
cells. This reduces the noise in the input count matrix and roughly separates the
genomic bins into signal blocks and background.

The segmentation minimizes a cost function within each signal block to aggregate
neighboring bins into regions of high co-occurrence. The cost is based on the
normalized co-occurrence of the binarized input count matrix. It is the cost of kernel
change point detection with a cosine kernel, and the segmentation is calculated with
`ruptures <https://centre-borelli.github.io/ruptures-docs/>`_.

The result of the algorithm is a list of regions in the genome that co-occur for the
genomic signal from the input count matrix.

Notation
--------

In addition to the :doc:`common notation <notation>`, this page uses:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Symbol
     - Meaning
   * - :math:`B \in \{0, 1\}^{N \times P}`
     - Binarized counts: :math:`B_{ij} = 1` if :math:`x_{ij} > 0`, 0 otherwise
   * - :math:`f_j = \sum_{i=1}^{N} B_{ij}`
     - Number of cells with signal in bin :math:`j`
   * - :math:`x_\mathrm{pct}`
     - Minimum percentage of cells with signal for a bin to pass the prefilter
       (``--minCellPct``)
   * - :math:`k`
     - Slack parameter: the number of neighbor bins kept by the prefilter (``--slack``)
   * - :math:`\mathcal{B}`
     - Set of bins that pass the prefilter
   * - :math:`\widetilde{B}`
     - :math:`B` after the prefilter
   * - :math:`\mathcal{C}`, :math:`B_\mathcal{C}`, :math:`w`
     - A signal block, the columns of :math:`\widetilde{B}` in it, and its number of
       bins
   * - :math:`O = B_\mathcal{C}^\top B_\mathcal{C}`
     - Co-occurrence: :math:`O_{jl}` is the number of cells with signal in bins
       :math:`j` and :math:`l`, with :math:`O_{jj} = f_j`
   * - :math:`W = D^{-1/2} O D^{-1/2}`
     - Normalized co-occurrence of the block, with :math:`D = \operatorname{diag}(O)`
   * - :math:`[a, b)`, :math:`m`
     - A segment made up of the bins :math:`a, \dots, b-1` of a block, and its number
       of bins :math:`m = b - a`
   * - :math:`C(a, b)`
     - Cost of a segment
   * - :math:`S`
     - Number of segments in a block
   * - :math:`\lambda`
     - Penalty for each segment, in units of bins

Algorithm
---------

**1. Binarize**

Set the nonzero entries of :math:`X` to 1.

.. math::

  B_{ij} = \begin{cases}
  1, & \text{if } x_{ij} > 0, \\
  0, & \text{otherwise.}
  \end{cases}

**2. Prefilter**

We denote :math:`f_j` as the number of cells with counts in bin :math:`j`:

.. math::

  f_j = \sum_{i=1}^N B_{ij}

A bin passes the prefilter if it has nonzero counts in at least :math:`x_\mathrm{pct}`
percentage of cells:

.. math::

   \mathcal{B} = \left\{ j \mid f_j \geq N \cdot \frac{x_\mathrm{pct}}{100} \right\}

The bins in :math:`\mathcal{B}`, and the bins at most :math:`k` bins from one of them,
are kept. The entries of all other bins are set to 0:

.. math::

   \widetilde{B}_{ij} = \begin{cases}
   B_{ij}, & \text{if } \exists\, b \in \mathcal{B}: |j - b| \leq k, \\
   0, & \text{otherwise.}
   \end{cases}

**3. Signal blocks**

:math:`\widetilde{B}` is split into blocks :math:`\mathcal{C}` of adjacent bins with
signal. A bin without signal in :math:`\widetilde{B}` always splits two blocks. Thus,
all bins in a block have :math:`f_j \geq 1`.

**4. Normalized co-occurrence**

We compute the cosine similarity :math:`W_{jl}` of the columns :math:`B_{\cdot j}` and
:math:`B_{\cdot l}` for each block:

.. math::

  W_{jl} = \frac{O_{jl}}{\sqrt{f_j f_l}}

This is the number of cells with signal in both bins, divided by the geometric mean of
the total number of cells with signal in each bin. Its values are in :math:`[0, 1]`,
and all the diagonal elements are :math:`W_{jj} = 1`.

We can write this in matrix notation as:

.. math::

  W = D^{-1/2}\, O\, D^{-1/2}, \quad D = \operatorname{diag}(O)

**5. Segment cost**

The cost of a segment :math:`[a, b)` of :math:`m` bins is

.. math::

   C(a, b) = \sum_{j=a}^{b-1} W_{jj} \;-\; \frac{1}{m} \sum_{j=a}^{b-1} \sum_{l=a}^{b-1} W_{jl}

This is the cost of kernel change point detection with a cosine kernel
:math:`\kappa(B_{\cdot j}, B_{\cdot l}) = W_{jl}`.

All :math:`W_{jj}` are 1, so we can also write the cost as

.. math::

   C(a, b) = m \big(1 - \overline{W}_{[a, b)}\big)

with :math:`\overline{W}_{[a, b)}` the mean of all :math:`m^2` entries of :math:`W` in
the segment, diagonal included. Thus, the cost is the length of the segment multiplied
by its mean cosine dissimilarity. It follows that:

- A segment of one bin costs 0.
- A segment of bins that share all their cells (:math:`\overline{W}_{[a, b)} = 1`)
  also costs 0.
- A segment of bins that share no cells (:math:`\overline{W}_{[a, b)} = 1/m`) costs
  :math:`m - 1`.

**6. Segmentation**

The segments of a block are determined by the minimization of

.. math::

   \sum_{s=1}^{S} C(a_{s-1}, a_s) \;+\; \lambda S,
   \quad 0 = a_0 < a_1 < \dots < a_S = w

with :math:`a_1, \dots, a_{S-1}` the ends of the segments. The minimum is calculated
with the PELT algorithm (`Killick et al., 2012
<https://doi.org/10.1080/01621459.2012.737745>`_), as implemented in
``ruptures.KernelCPD`` with the cosine kernel and a minimum segment length of one bin.
PELT gives the exact minimum, with no limit on the length of a segment.

If two segmentations have the same minimum, the result depends on floating-point
rounding. Ties occur at :math:`\lambda = 1`. For example, :math:`n` bins that share no
cells cost :math:`(n - 1) + \lambda` as one segment, and :math:`n \lambda` as :math:`n`
segments of one bin. At :math:`\lambda = 1`, these two costs are the same.

Choosing a penalty
~~~~~~~~~~~~~~~~~~

Each bin in a segment of more than one bin adds about :math:`1 - \bar{W}` to the cost,
where :math:`\bar{W}` is the mean off-diagonal value of :math:`W` in the segment. Adding
a segment adds :math:`\lambda` to the cost.

For sparse single-cell data the cosines are small, so the number of regions changes
quickly for penalties near 1: below about 0.9 almost every bin is a separate region,
and above about 1.1 almost every signal block is one region.

:ref:`scFindVCRs` uses :math:`\lambda = 1`, our recommended default. To use a different
penalty, call :func:`sincei.tools.VCRfinder.VCRfinder` changing the ``penalty`` argument.

Output
------

A 6-column BED file (``chrom``, ``start``, ``end``, ``name``, ``score``, ``strand``),
with one row per region. The score is the segment cost :math:`C(a, b)` (0 for a region
of one bin). See :ref:`the file formats <fmt-bed>` for details.

Space and time complexity
-------------------------

For a block of :math:`w` bins:

.. list-table::
   :header-rows: 1
   :widths: 40 35 25

   * - Step
     - Time
     - Space
   * - Binarize, prefilter and blocks
     - :math:`O(\mathrm{nnz} + P)`
     - :math:`O(\mathrm{nnz} + P)`
   * - Segmentation of one block (PELT)
     - :math:`O(w^2 N)` at most
     - :math:`O(w N)`
   * - Segment costs of one block
     - :math:`O(w N)`
     - :math:`O(w N)`

ruptures needs the block as a dense matrix of :math:`w \times N` values, but it does not
store :math:`W`: it calculates each kernel value when it needs it. The pruning of PELT
removes most of the :math:`w^2` pairs of bins when the segments are short, so the time
is usually much smaller than :math:`O(w^2 N)`.

The blocks are segmented in parallel by ``--numberOfProcessors`` processes, so the
working memory scales with :math:`O(n_\mathrm{proc}\, w_\mathrm{max} N)`.
