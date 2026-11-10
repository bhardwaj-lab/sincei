Notation
========

The method pages all use the following notation:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Symbol
     - Meaning
   * - :math:`N`, :math:`c`
     - Number of cells, and index of a cell
   * - :math:`P`, :math:`i`
     - Number of features (bins, genes, peaks, etc.), and index of a feature
   * - :math:`K`, :math:`k`
     - Number of components or topics, and index of a component or topic
   * - :math:`X \in \mathbb{R}^{N \times P}`, :math:`x_{ci}`
     - Data matrix (cells in rows, features in columns), and data entry for cell
       :math:`c` and feature :math:`i`
   * - :math:`\mathrm{nnz}`
     - Number of non-zero entries of :math:`X`
   * - :math:`U \in \mathbb{R}^{N \times K}`, :math:`u_c`
     - Scores of the cells (dimensionality reduction), and the row of cell :math:`c`
   * - :math:`V \in \mathbb{R}^{P \times K}`, :math:`v_i`
     - Loadings of the features, and the row of feature :math:`i`
   * - :math:`\mathbf{1}`
     - Vector of ones
   * - :math:`I_K`
     - Identity matrix of size :math:`K \times K`

Symbols that only one method uses are defined on the page of that method.
