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


def _as_list(value: Sequence[str] | str | Path | None) -> list[str]:
    """Wrap a single string or path in a list; copy any other sequence."""
    if value is None:
        return []
    if isinstance(value, str | Path):
        return [str(value)]
    return [str(v) for v in value]


def count_reads(
    bamFiles: Sequence[str] | str | Path,
    barcodes: Sequence[str] | str | Path | None = None,
    binLength: int = 10_000,
    stepSize: int | None = None,
    bedFile: str | None = None,
    cellTag: str = "BC",
    umiTag: str = "RX",
    groupTag: str | None = None,
    groupLabels: Sequence[str] | str | None = None,
    valueTag: str | None = None,
    region: str | None = None,
    blackListFileName: str | None = None,
    chrsToSkip: Sequence[str] | str | None = None,
    minMappingQuality: int | None = None,
    samFlag_include: int | None = None,
    samFlag_exclude: int | None = None,
    minFragmentLength: int = 0,
    maxFragmentLength: int = 0,
    filterRNAstrand: FilterRNAStrand | str | None = None,
    extendReads: int | None = None,
    center_read: bool = False,
    duplicateFilter: DuplicateFilter | str | None = None,
    motifFilter: Sequence[str] | str | None = None,
    genome2bit: str | None = None,
    GCcontentFilter: str | None = None,
    minAlignedFraction: float | None = None,
    featureId: Sequence[str] | str | None = None,
    exonId: Sequence[str] | str | None = None,
    featureIdTag: str | None = None,
    metagene: bool = False,
    genomeChunkSize: int | None = None,
    numberOfProcessors: int = 1,
    tmpDir: str | None = None,
    verbose: bool = False,
) -> ad.AnnData:
    r"""Count reads per cell barcode in genomic bins or features.

    The counting runs in the Rust backend of ``scCountReads``. The result goes
    to a temporary ``.h5ad`` file, which is read back and then deleted.

    Parameters
    ----------
    bamFiles : list of str or str
        Indexed BAM files, or one BAM file.
    barcodes : list of str or str, optional
        Cell barcodes to count, or the path to a file with one barcode per line.
        If not given, every barcode found in the ``cellTag`` of the BAM files is
        counted, and only the (sample, barcode) pairs with counts become cells.
    binLength : int
        Bin size in bp. Ignored when ``bedFile`` is given.
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
    groupLabels : list of str or str, optional
        One sample label per BAM file. Defaults to the file names.
    valueTag : str, optional
        BAM tag whose value is added instead of 1 per read.
    region : str, optional
        Region to count, as ``chrom[:start-end]``.
    blackListFileName : str, optional
        BED file with regions to exclude.
    chrsToSkip : list of str or str, optional
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
    motifFilter : list of str or str, optional
        ``"read_motif,ref_motif"`` pairs. Needs ``genome2bit``.
    genome2bit : str, optional
        Reference genome in 2bit format.
    GCcontentFilter : str, optional
        ``"<low>,<high>"`` GC content limits of a read.
    minAlignedFraction : float, optional
        Minimum aligned fraction of a read.
    featureId, exonId : list of str or str, optional
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

    Returns
    -------
    AnnData
        Counts in ``X`` (cells x bins or features). ``obs`` has the ``sample``
        and ``barcode`` columns; ``var`` has ``chrom``, ``start``, ``end`` and
        ``name``.

    Examples
    --------
    >>> adata = count_reads("cells.bam", "bc.txt", binLength=5000)
    """
    bam_files = _as_list(bamFiles)
    if isinstance(barcodes, str | Path):
        barcodes = backend.read_barcodes(str(barcodes))
    backend.require_single_bam_for_group_tag(bam_files, groupTag)

    dup_filter = DuplicateFilter(duplicateFilter) if duplicateFilter else None
    rna_strand = FilterRNAStrand(filterRNAstrand) if filterRNAstrand else None
    min_gc, max_gc = backend.parse_gc_content(GCcontentFilter)

    kwargs: dict[str, Any] = {
        "barcodes": list(barcodes) if barcodes is not None else None,
        "labels": [] if groupTag is not None else _as_list(groupLabels),
        "bc_tag": cellTag,
        "umi_tag": backend.umi_tag_if_used(umiTag, dup_filter),
        "count_tag": valueTag,
        "group_tag": groupTag,
        "min_mapq": minMappingQuality,
        "sam_flag_include": samFlag_include,
        "sam_flag_exclude": samFlag_exclude,
        "chr_to_skip": _as_list(chrsToSkip),
        "region": region,
        "blacklist_path": blackListFileName,
        "extend_reads": extendReads,
        "center_reads": center_read,
        "dup_method": backend.dup_method(dup_filter),
        "filter_rna_strand": rna_strand.value if rna_strand else None,
        "genome_2bit": genome2bit,
        "motif_filter": backend.parse_motif_filter(_as_list(motifFilter)),
        "min_gc": min_gc,
        "max_gc": max_gc,
        "min_aligned_fraction": minAlignedFraction,
        "min_fragment_length": backend.optional_length(minFragmentLength),
        "max_fragment_length": backend.optional_length(maxFragmentLength),
        "num_threads": numberOfProcessors,
    }
    if genomeChunkSize:
        kwargs["chunk_size"] = genomeChunkSize
    if verbose:
        log.info("Counting %s with %s", bam_files, kwargs)

    with tempfile.TemporaryDirectory(dir=tmpDir) as tmp:
        output_path = str(Path(tmp) / "counts.h5ad")
        if bedFile is None:
            internal.count_bins(
                bam_files,
                output_path=output_path,
                bin_size=binLength,
                step_size=stepSize if stepSize is not None else binLength,
                **kwargs,
            )
        else:
            internal.count_features(
                bam_files,
                bedFile,
                output_path=output_path,
                feature_type=_as_list(featureId) or None,
                exon_type=_as_list(exonId) or None,
                name_attr=featureIdTag,
                metagene=metagene,
                **kwargs,
            )

        return ad.read_h5ad(output_path)
