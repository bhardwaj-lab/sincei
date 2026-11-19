Variable Chromatin Region Finder
================================

:ref:`scFindVCRs` segments the chromosomes in an input AnnData into variable chromatin
regions (VCRs): bins with signal in the same cells. In other words, co-occurring bins.

It takes a count matrix :math:`X \in \mathbb{R}^{N \times P}` of :math:`N` cells
and :math:`P` genomic bins of equal size (from :ref:`scCountReads` ``bins``) as input,
and processes it in two stages: a prefilter that splits each chromosome into signal
blocks, and a segmentation of each block into regions.

The prefilter keeps bins which have nonzero counts in at least a percentage of the total
cells. This reduces the noise in the input count matrix and roughly separates the
genomic bins into signal blocks and background.

The segmentation minimizes a cost function within each signal block to aggregate
neighboring bins into regions of high co-occurrence. The cost is based on the
normalized co-occurrence of the binarized input count matrix.

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
     - Slack parameter, in bins (``--slack``)
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
     - Normalized co-occurrence of the block
   * - :math:`[a, b)`, :math:`m`
     - A segment made up of the bins :math:`a, \dots, b-1` of a block, and its number
       of bins :math:`m = b - a`
   * - :math:`C(a, b)`
     - Cost of a segment
   * - :math:`S`
     - Number of segments in a block
   * - :math:`\lambda`
     - Penalty for each segment, in units of bins (``--penalty``)
   * - :math:`L`
     - Maximum segment length, in bins (``--maxBinsRegion``)

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

:math:`\widetilde{B}` is split into blocks :math:`\mathcal{C}` of
contiguous bins with signal, with at most :math:`k` bins without signal between two bins
with signal. Bins outside all blocks get no region. A block of one bin is one region.

**4. Normalized co-occurrence**

We compute the cosine similarity :math:`W_{jl}` of the columns :math:`B_{\cdot j}` and
:math:`B_{\cdot l}` for each block:

.. math::

  W_{jl} = \frac{O_{jl}}{\sqrt{f_j f_l}}

This is the number of cells with signal in both bins, divided by the geometric mean of
the total number of cells with signal in each bin. It is in :math:`[0, 1]`, and
:math:`W_{jj} = 1` for every bin with signal.

We can write this in matrix notation as:

.. math::

  W = D^{-1/2}\, O\, D^{-1/2}, \qquad D_{jj} = \max(O_{jj},\, 1)

With :math:`k > 0`, a block can contain bins without signal (:math:`O_{jj} = 0`). For
these bins :math:`D_{jj}` is set to 1 to prevent a division by 0, so their row and
column of :math:`W` are 0.

**5. Segment cost**

The cost of a segment :math:`[a, b)` of :math:`m` bins is

.. math::

   C(a, b) = \sum_{j=a}^{b-1} W_{jj} \;-\; \frac{1}{m} \sum_{j=a}^{b-1} \sum_{l=a}^{b-1} W_{jl}

In other terms, if all bins in the segment have signal its cost is

.. math::

   C(a, b) = m \big(1 - \overline{W}_{[a, b)}\big)

with :math:`\overline{W}_{[a, b)}` the mean of all :math:`m^2` entries of :math:`W` in
the segment, diagonal included. Thus, the cost is the length of the segment multiplied
by its mean dissimilarity. A segment of one bin costs 0. A segment of bins that share
all their cells (:math:`\overline{W}_{[a, b)} = 1`) also costs 0. A segment of bins
that share no cells (:math:`\overline{W}_{[a, b)} = 1/m`) costs :math:`m - 1`.

**6. Segmentation**

The segments of a block are determined by the minimization of

.. math::

   \sum_{s=1}^{S} C(a_{s-1}, a_s) \;+\; \lambda S,
   \qquad 0 = a_0 < a_1 < \dots < a_S = w, \qquad a_s - a_{s-1} \le L

with :math:`a_1, \dots, a_{S-1}` the ends of the segments. The minimum is calculated
by dynamic programming. :math:`F(b)` is the minimum cost of the bins :math:`0, \dots,
b-1`, and the last segment of this minimum is :math:`[a, b)`:

.. math::

   F(0) = 0, \qquad F(b) = \min_{\max(0,\, b - L) \le a < b}
   \big[ F(a) + C(a, b) + \lambda \big]

Tracing the best :math:`a` back from :math:`b = w` gives the segment ends. If two values
of :math:`a` have the same minimum, the smallest :math:`a` is chosen, and thus the
longest last segment :math:`[a, b)`.

Ties occur at :math:`\lambda = 1`. For example, :math:`n` bins that share no cells cost
:math:`(n - 1) + \lambda` as one segment, and :math:`n \lambda` as :math:`n` segments of
one bin. At :math:`\lambda = 1`, these two costs are the same.

Maximum segment length
~~~~~~~~~~~~~~~~~~~~~~

The memory required to compute the full normalized co-occurrence matrix :math:`W` of a
block scales with :math:`O(w^2)`. To avoid a memory blow-up, we limit segments to a
maximum length of :math:`L` bins (``--maxBinsRegion``, 2,000 bins by default). This allows us to compute
only the band :math:`W_{jl}` with :math:`0 \le l - j < L` for the segmentation. For a
block of length :math:`L` or shorter, the segmentation is exact. For larger blocks, the
result is the same as the exact segmentation if no exact segment is longer than
:math:`L` bins.

In practice, such long blocks of signal rarely occur in the genome, so this heuristic
only safeguards a hypothetical the worst case scenario.

Choosing a penalty
~~~~~~~~~~~~~~~~~~

Each bin in a segment of more than one bin adds about :math:`1 - \bar{W}` to the cost,
where :math:`\bar{W}` is the mean off-diagonal value of :math:`W` in the segment. Adding
a segment adds :math:`\lambda` to the cost.

For sparse single-cell data the cosines are small, so the number of regions changes
quickly for penalties near 1: below about 0.9 almost every bin is a separate region,
and above about 1.1 almost every signal block is one region.

Output
------

A 6-column BED file (``chrom``, ``start``, ``end``, ``name``, ``score``, ``strand``),
with one row per region. The score is the segment cost :math:`C(a, b)` (0 for a region
of one bin). See :ref:`the file formats
<fmt-bed>` for details.

Space and time complexity
-------------------------

For a block of :math:`w` bins, :math:`r_i` is the number of its bins with signal in cell
:math:`i`.

.. list-table::
   :header-rows: 1
   :widths: 40 35 25

   * - Step
     - Time
     - Space
   * - Binarize, prefilter and blocks
     - :math:`O(\mathrm{nnz} + P)`
     - :math:`O(\mathrm{nnz} + P)`
   * - Band of :math:`W` and DP of one block
     - :math:`O\big(\sum_i r_i \min(r_i, L) + w L\big)`
     - :math:`O(w L)`

The blocks are segmented one at a time, so the working memory scales with
:math:`O(w_\mathrm{max} L)`.  Larger values of the slack parameter :math:`k` give wider
blocks, so the :math:`w L` terms increase.
