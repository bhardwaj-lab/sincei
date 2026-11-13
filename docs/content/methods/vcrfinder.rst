Variable Chromatin Region Finder
================================

:ref:`scFindVCRs` segments the chromosomes of an input AnnData into variable chromatin
regions (VCRs): runs of bins whose signal is in the same cells.

It starts from a count matrix :math:`X \in \mathbb{R}^{N \times P}` of :math:`N` cells
and :math:`P` genomic bins of equal size (from :ref:`scCountReads` ``bins``), and
processes each chromosome separately. It has two parts: a prefilter that splits the
chromosome into signal blocks, and an exact segmentation of each block with a cost
based on normalized co-occurrence.

Notation
--------

In addition to the :doc:`common notation <notation>`, this page uses:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Symbol
     - Meaning
   * - :math:`B \in \{0, 1\}^{N \times P}`
     - Binarized counts: :math:`B_{ci} = 1` if :math:`x_{ci} > 0`, 0 otherwise
   * - :math:`f_i = \sum_{c=1}^{N} B_{ci}`
     - Number of cells with signal in bin :math:`i`
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
     - Co-occurrence: :math:`O_{ij}` is the number of cells with signal in bins
       :math:`i` and :math:`j`, with :math:`O_{ii} = f_i`
   * - :math:`D = \operatorname{diag}(O)`
     - Degree matrix of the block (1 for a bin without signal)
   * - :math:`W = D^{-1/2} O D^{-1/2}`
     - Normalized co-occurrence of the block
   * - :math:`[a, b)`, :math:`m`
     - A segment of the bins :math:`a, \dots, b-1` of a block, and its number of bins
       :math:`m = b - a`
   * - :math:`C(a, b)`
     - Cost of a segment
   * - :math:`S`
     - Number of segments in a block
   * - :math:`\lambda`
     - Penalty for each segment, in units of bins (``--penalty``)

Algorithm
---------

**1. Binarize**

Set the nonzero entries of :math:`X` to 1.

.. math::

  {B}_{ci} = \begin{cases}
  1, & \text{if } X_{cf} > 0, \\
  0, & \text{otherwise.}
  \end{cases}

**2. Prefilter**

We denote :math:`f_i` as the fraction of cells with counts in bin :math:`i`:

.. math::

  f_i = \sum_c^N B_{ci}

A bin passes the prefilter if it has nonzero counts in more than the threshold fraction:

.. math::

   \mathcal{B} = \left\{ i \mid f_i \geq N \cdot \frac{x_\mathrm{pct}}{100} \right\}

The bins in :math:`\mathcal{B}`, and the bins at most :math:`k` bins from one of them,
are kept. The entries of all other bins are set to 0:

.. math::

   \widetilde{B}_{ci} = \begin{cases}
   B_{ci}, & \text{if } \exists\, b \in \mathcal{B}: |i - b| \leq k, \\
   0, & \text{otherwise.}
   \end{cases}

**3. Signal blocks**

:math:`\widetilde{B}` is split into blocks :math:`\mathcal{C}` of
contiguous bins with signal, with at most :math:`k` bins without signal between two bins
with signal. Bins outside all blocks get no region. A block of one bin is one region.

**4. Normalized co-occurrence**

For each block:

.. math::

   W = D^{-1/2}\, O\, D^{-1/2},
   \qquad
   W_{ij} = \frac{O_{ij}}{\sqrt{f_i f_j}}

:math:`W_{ij}` is the cosine similarity of the columns :math:`B_{\cdot i}` and
:math:`B_{\cdot j}`: the number of cells with signal in both bins, divided by the
geometric mean of the number of cells in each bin. It is in :math:`[0, 1]`, and
:math:`W_{ii} = 1` for every bin with signal. With :math:`k > 0`, a block can contain bins
without signal (:math:`O_{ii} = 0`); for these bins :math:`D_{ii}` is set to 1 to prevent
a division by 0, so their row, column and diagonal of :math:`W` are 0.

**5. Segment cost**

The cost of a segment :math:`[a, b)` of :math:`m` bins is

.. math::

   C(a, b) = \sum_{i=a}^{b-1} W_{ii} \;-\; \frac{1}{m} \sum_{i=a}^{b-1} \sum_{j=a}^{b-1} W_{ij}

This is the kernel cost of change-point detection with a cosine kernel. The first term,
summed over a block, is the same for every segmentation. Thus a low total cost means a
high co-occurrence inside the segments, divided by their length. With the 2D prefix sums
of :math:`W`, each :math:`C(a, b)` costs :math:`O(1)`.

**6. Segmentation**

The segments of a block are the exact minimum of

.. math::

   \sum_{s=1}^{S} C(a_{s-1}, a_s) \;+\; \lambda S,
   \qquad 0 = a_0 < a_1 < \dots < a_S = w

calculated by dynamic programming over all start positions:

.. math::

   F(0) = 0, \qquad F(b) = \min_{0 \le a < b} \big[ F(a) + C(a, b) + \lambda \big]

Tracing the best :math:`a` back from :math:`b = w` gives the segment ends.

Choosing a penalty
~~~~~~~~~~~~~~~~~~

Each bin in a segment of more than one bin adds about :math:`1 - \bar{W}` to the cost,
where :math:`\bar{W}` is the mean off-diagonal value of :math:`W` in the segment, and each
segment adds :math:`\lambda`. A penalty of 1 is therefore about the cost of one bin. On
sparse single-cell data the cosines are small, so the number of regions changes quickly
for penalties near 1: below about 0.9 almost every bin is a separate region, and above
about 1.1 almost every signal block is one region.

Output
------

A 6-column BED file (``chrom``, ``start``, ``end``, ``name``, ``score``, ``strand``),
with one row per region. The score is the segment cost :math:`C(a, b)` (0 for a region
of one bin). See :ref:`the file formats
<fmt-bed>` for details.

Space and time complexity
-------------------------

For a block of :math:`w` bins, :math:`r_c` is the number of its bins with signal in cell
:math:`c`.

.. list-table::
   :header-rows: 1
   :widths: 50 25 25

   * - Step
     - Time
     - Space
   * - Binarize, prefilter and blocks
     - :math:`O(\mathrm{nnz} + P)`
     - :math:`O(\mathrm{nnz} + P)`
   * - :math:`O` of one block
     - :math:`O\big(\sum_c r_c^2 + w^2\big)`
     - :math:`O(w^2)`
   * - :math:`W`, prefix sums and dynamic programming of one block
     - :math:`O(w^2)`
     - :math:`O(w^2)`

The blocks are segmented one at a time (in parallel with ``-p``), so the working memory
scales with :math:`O(w_\mathrm{max}^2)`. Runtime does not depend on the penalty or on
the number of regions. A larger :math:`k` gives wider blocks, so :math:`w^2` terms
increase.
