from __future__ import annotations

import logging
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Any

import anndata as ad

from sincei import _sincei as internal
from sincei.cli import _parsers as backend
from sincei.cli._common_args import DuplicateFilter, FilterRNAStrand

if TYPE_CHECKING:
    from collections.abc import Sequence

log = logging.getLogger(__name__)


class CountReadsPerBin:
    r"""Count reads per cell barcode in genomic bins or features.

    The counting runs in the Rust backend of ``scCountReads``. The result goes
    to a temporary ``.h5ad`` file, which is read back and then deleted.

    Parameters
    ----------
    bamFilesList : list of str
        Indexed BAM files.
    binLength : int
        Bin size in bp. Ignored when ``bedFile`` is given.
    barcodes : list of str or str
        Cell barcodes to count, or the path to a file with one barcode per line.
    stepSize : int, optional
        Distance between bin starts. Defaults to ``binLength`` (contiguous bins).
    bedFile : str, optional
        BED/GTF/GFF file. If given, reads are counted on its features instead of
        bins.
    cellTag : str
        BAM tag that holds the cell barcode.
    umiTag : str
        BAM tag that holds the UMI. Used only by UMI duplicate filters.
    groupTag : str, optional
        BAM tag that holds the sample group, for a merged BAM. Needs exactly one
        BAM file; ``groupLabels`` is then ignored.
    groupLabels : list of str, optional
        One sample label per BAM file. Defaults to the file names.
    valueTag : str, optional
        BAM tag whose value is added instead of 1 per read.
    region : str, optional
        Region to count, as ``chrom[:start-end]``.
    blackListFileName : str, optional
        BED file with regions to exclude.
    chrsToSkip : list of str, optional
        Chromosomes to exclude.
    minMappingQuality : int, optional
        Minimum mapping quality of a read.
    samFlag_include, samFlag_exclude : int, optional
        SAM flag bits a read must have / must not have.
    minFragmentLength, maxFragmentLength : int
        Fragment length limits. 0 means no limit.
    filterRNAstrand : {"forward", "reverse"}, optional
        Keep only reads from this RNA strand.
    extendReads : int, optional
        Extend reads to this fragment length.
    center_read : bool
        Count only the centre of each fragment.
    duplicateFilter : DuplicateFilter or str, optional
        Remove duplicates with this method: ``"start_bc"``, ``"start_bc_umi"``,
        ``"start_end_bc"`` or ``"start_end_bc_umi"``.
    motifFilter : list of str, optional
        ``"read_motif,ref_motif"`` pairs. Needs ``genome2bit``.
    genome2bit : str, optional
        Reference genome in 2bit format.
    GCcontentFilter : str, optional
        ``"<low>,<high>"`` GC content limits of a read.
    minAlignedFraction : float, optional
        Minimum aligned fraction of a read.
    featureId, exonId : list of str, optional
        GTF/GFF feature and exon types to count. Ignored for BED files.
    featureIdTag : str, optional
        GTF/GFF attribute that names a feature. Ignored for BED files.
    metagene : bool
        Count GTF/GFF features over their exons only.
    genomeChunkSize : int, optional
        Genome chunk size given to each thread.
    numberOfProcessors : int
        Number of threads. 0 uses all cores.
    tmpDir : str, optional
        Directory for the temporary ``.h5ad`` file. Defaults to the system
        temporary directory.
    verbose : bool
        Log the backend parameters.

    Examples
    --------
    >>> counter = CountReadsPerBin(["cells.bam"], binLength=5000, barcodes="bc.txt")
    >>> adata = counter.run()
    """

    def __init__(
        self,
        bamFilesList: Sequence[str],
        binLength: int = 10_000,
        barcodes: Sequence[str] | str | Path | None = None,
        stepSize: int | None = None,
        bedFile: str | None = None,
        cellTag: str = "BC",
        umiTag: str = "RX",
        groupTag: str | None = None,
        groupLabels: Sequence[str] | None = None,
        valueTag: str | None = None,
        region: str | None = None,
        blackListFileName: str | None = None,
        chrsToSkip: Sequence[str] | None = None,
        minMappingQuality: int | None = None,
        samFlag_include: int | None = None,
        samFlag_exclude: int | None = None,
        minFragmentLength: int = 0,
        maxFragmentLength: int = 0,
        filterRNAstrand: FilterRNAStrand | str | None = None,
        extendReads: int | None = None,
        center_read: bool = False,
        duplicateFilter: DuplicateFilter | str | None = None,
        motifFilter: Sequence[str] | None = None,
        genome2bit: str | None = None,
        GCcontentFilter: str | None = None,
        minAlignedFraction: float | None = None,
        featureId: Sequence[str] | None = None,
        exonId: Sequence[str] | None = None,
        featureIdTag: str | None = None,
        metagene: bool = False,
        genomeChunkSize: int | None = None,
        numberOfProcessors: int = 1,
        tmpDir: str | None = None,
        verbose: bool = False,
    ) -> None:
        if barcodes is None:
            msg = "barcodes is required: give a list of barcodes or a barcode file."
            raise ValueError(msg)

        self.bamFilesList = list(bamFilesList)
        self.binLength = binLength
        self.barcodes = (
            backend.read_barcodes(str(barcodes))
            if isinstance(barcodes, str | Path)
            else list(barcodes)
        )
        self.stepSize = stepSize if stepSize is not None else binLength
        self.bedFile = bedFile
        self.cellTag = cellTag
        self.umiTag = umiTag
        self.groupTag = groupTag
        self.groupLabels = list(groupLabels) if groupLabels else []
        self.valueTag = valueTag
        self.region = region
        self.blackListFileName = blackListFileName
        self.chrsToSkip = list(chrsToSkip) if chrsToSkip else []
        self.minMappingQuality = minMappingQuality
        self.samFlag_include = samFlag_include
        self.samFlag_exclude = samFlag_exclude
        self.minFragmentLength = minFragmentLength
        self.maxFragmentLength = maxFragmentLength
        self.filterRNAstrand = (
            FilterRNAStrand(filterRNAstrand) if filterRNAstrand is not None else None
        )
        self.extendReads = extendReads
        self.center_read = center_read
        self.duplicateFilter = (
            DuplicateFilter(duplicateFilter) if duplicateFilter is not None else None
        )
        self.motifFilter = list(motifFilter) if motifFilter else None
        self.genome2bit = genome2bit
        self.GCcontentFilter = GCcontentFilter
        self.minAlignedFraction = minAlignedFraction
        self.featureId = list(featureId) if featureId else None
        self.exonId = list(exonId) if exonId else None
        self.featureIdTag = featureIdTag
        self.metagene = metagene
        self.genomeChunkSize = genomeChunkSize
        self.numberOfProcessors = numberOfProcessors
        self.tmpDir = tmpDir
        self.verbose = verbose

    def _backend_kwargs(self, output_path: str) -> dict[str, Any]:
        backend.require_single_bam_for_group_tag(self.bamFilesList, self.groupTag)
        min_gc, max_gc = backend.parse_gc_content(self.GCcontentFilter)

        kwargs: dict[str, Any] = {
            "barcodes": self.barcodes,
            "labels": [] if self.groupTag is not None else self.groupLabels,
            "output_path": output_path,
            "bc_tag": self.cellTag,
            "umi_tag": backend.umi_tag_if_used(self.umiTag, self.duplicateFilter),
            "count_tag": self.valueTag,
            "group_tag": self.groupTag,
            "min_mapq": self.minMappingQuality,
            "sam_flag_include": self.samFlag_include,
            "sam_flag_exclude": self.samFlag_exclude,
            "chr_to_skip": self.chrsToSkip,
            "region": self.region,
            "blacklist_path": self.blackListFileName,
            "extend_reads": self.extendReads,
            "center_reads": self.center_read,
            "dup_method": backend.dup_method(self.duplicateFilter),
            "filter_rna_strand": (
                self.filterRNAstrand.value if self.filterRNAstrand else None
            ),
            "genome_2bit": self.genome2bit,
            "motif_filter": backend.parse_motif_filter(self.motifFilter),
            "min_gc": min_gc,
            "max_gc": max_gc,
            "min_aligned_fraction": self.minAlignedFraction,
            "min_fragment_length": backend.optional_length(self.minFragmentLength),
            "max_fragment_length": backend.optional_length(self.maxFragmentLength),
            "num_threads": self.numberOfProcessors,
        }
        if self.genomeChunkSize:
            kwargs["chunk_size"] = self.genomeChunkSize
        return kwargs

    def run(self) -> ad.AnnData:
        r"""Count the reads and return an AnnData object with the results.

        Returns
        -------
        AnnData
            Counts in ``X`` (cells x bins or features). ``obs`` has the
            ``sample`` and ``barcode`` columns; ``var`` has ``chrom``, ``start``,
            ``end`` and ``name``.
        """
        with tempfile.TemporaryDirectory(dir=self.tmpDir) as tmp:
            output_path = str(Path(tmp) / "counts.h5ad")
            kwargs = self._backend_kwargs(output_path)
            if self.verbose:
                log.info("Counting %s with %s", self.bamFilesList, kwargs)

            if self.bedFile is None:
                internal.count_bins(
                    self.bamFilesList,
                    bin_size=self.binLength,
                    step_size=self.stepSize,
                    **kwargs,
                )
            else:
                internal.count_features(
                    self.bamFilesList,
                    self.bedFile,
                    feature_type=self.featureId,
                    exon_type=self.exonId,
                    name_attr=self.featureIdTag,
                    metagene=self.metagene,
                    **kwargs,
                )

            return ad.read_h5ad(output_path)
