//! Python access to the annotation parsers and their overlap index.

use std::path::PathBuf;

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::annotation::{Feature, GenomeIndex, parse_annotation_files};

/// A parsed genome annotation: its features, and an index to find them by
/// position.
#[pyclass(frozen, module = "sincei._sincei")]
pub struct GenomeAnnotation {
    index: GenomeIndex,
    features: Vec<Feature>,
}

impl GenomeAnnotation {
    fn overlapping(&self, chrom: &str, start: usize, end: usize) -> Vec<usize> {
        let Some(chrom_index) = self.index.get(chrom) else {
            return Vec::new();
        };
        let mut hits: Vec<usize> = chrom_index.find(start, end).map(|iv| iv.var_idx).collect();
        hits.sort_unstable();
        hits.dedup();
        hits
    }
}

#[pymethods]
impl GenomeAnnotation {
    /// Indices of the features that overlap the half-open region `[start, end)`.
    ///
    /// Index `i` is row `i` of `features()`. A chromosome that is not in the
    /// annotation gives an empty list.
    fn find_overlaps(&self, chrom: &str, start: usize, end: usize) -> Vec<usize> {
        self.overlapping(chrom, start, end)
    }

    /// The features as columns (`chrom`, `start`, `end`, `name`, `strand`), in
    /// file order. Coordinates are 0-based half-open; a missing name is `None`.
    fn features<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let features = &self.features;
        let columns = PyDict::new(py);
        columns.set_item(
            "chrom",
            features
                .iter()
                .map(|f| f.chrom.as_str())
                .collect::<Vec<_>>(),
        )?;
        columns.set_item(
            "start",
            features.iter().map(|f| f.start).collect::<Vec<_>>(),
        )?;
        columns.set_item("end", features.iter().map(|f| f.end).collect::<Vec<_>>())?;
        columns.set_item(
            "name",
            features
                .iter()
                .map(|f| f.name.as_deref())
                .collect::<Vec<_>>(),
        )?;
        columns.set_item(
            "strand",
            features.iter().map(|f| f.strand).collect::<Vec<_>>(),
        )?;
        Ok(columns)
    }

    /// Names of the chromosomes in the annotation, sorted.
    #[getter]
    fn chroms(&self) -> Vec<String> {
        let mut chroms: Vec<String> = self.index.keys().cloned().collect();
        chroms.sort_unstable();
        chroms
    }

    fn __len__(&self) -> usize {
        self.features.len()
    }
}

/// Parse BED / GTF / GFF3 files into one searchable annotation.
///
/// The files are merged in the given order. `feature_types`, `exon_types`,
/// `name_attr` and `metagene` are those of `count_features`, and only affect
/// GTF / GFF3 files. By default a GTF / GFF3 gives one feature per gene.
#[pyfunction(signature = (
    paths,
    feature_types = None,
    exon_types = None,
    name_attr = None,
    metagene = false,
))]
pub fn parse_annotation(
    py: Python<'_>,
    paths: Vec<PathBuf>,
    feature_types: Option<Vec<String>>,
    exon_types: Option<Vec<String>>,
    name_attr: Option<String>,
    metagene: bool,
) -> PyResult<GenomeAnnotation> {
    py.detach(|| {
        let feature_types: Option<Vec<&str>> = feature_types
            .as_ref()
            .map(|t| t.iter().map(String::as_str).collect());
        let exon_types: Option<Vec<&str>> = exon_types
            .as_ref()
            .map(|t| t.iter().map(String::as_str).collect());
        parse_annotation_files(
            &paths,
            feature_types.as_deref(),
            exon_types.as_deref(),
            name_attr.as_deref(),
            metagene,
        )
    })
    .map(|(index, features)| GenomeAnnotation { index, features })
    .map_err(|e| PyRuntimeError::new_err(format!("{e:#}")))
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use super::*;

    fn annotation(suffix: &str, text: &str, metagene: bool) -> GenomeAnnotation {
        let mut file = tempfile::Builder::new().suffix(suffix).tempfile().unwrap();
        file.write_all(text.as_bytes()).unwrap();
        let (index, features) =
            parse_annotation_files([file.path()], None, None, None, metagene).unwrap();
        GenomeAnnotation { index, features }
    }

    const BED: &str = "chr1\t100\t200\ta\t0\t+\n\
                       chr1\t150\t300\tb\t0\t-\n\
                       chr2\t0\t50\tc\t0\t.\n";

    const GTF: &str = "chr1\tsrc\tgene\t101\t400\t.\t+\t.\tgene_id \"g1\";\n\
                       chr1\tsrc\texon\t101\t200\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n\
                       chr1\tsrc\texon\t301\t400\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n";

    #[test]
    fn a_bed_file_keeps_its_features_in_file_order() {
        let ann = annotation(".bed", BED, false);

        let names: Vec<_> = ann.features.iter().map(|f| f.name.as_deref()).collect();
        let strands: Vec<_> = ann.features.iter().map(|f| f.strand).collect();
        assert_eq!(names, [Some("a"), Some("b"), Some("c")]);
        assert_eq!(strands, ['+', '-', '*']);
        assert_eq!(ann.chroms(), ["chr1", "chr2"]);
    }

    #[test]
    fn overlaps_are_half_open() {
        let ann = annotation(".bed", BED, false);

        assert_eq!(ann.overlapping("chr1", 199, 200), [0, 1]);
        assert_eq!(ann.overlapping("chr1", 200, 250), [1]);
        assert!(ann.overlapping("chr1", 50, 100).is_empty());
        assert!(ann.overlapping("chr1", 300, 400).is_empty());
    }

    #[test]
    fn an_unknown_chromosome_has_no_overlaps() {
        let ann = annotation(".bed", BED, false);

        assert!(ann.overlapping("1", 0, 1_000).is_empty());
    }

    #[test]
    fn a_gtf_gives_one_feature_per_gene() {
        let ann = annotation(".gtf", GTF, false);

        assert_eq!(ann.features.len(), 1);
        assert_eq!(ann.features[0].name.as_deref(), Some("g1"));
        assert_eq!((ann.features[0].start, ann.features[0].end), (100, 400));
        assert_eq!(ann.overlapping("chr1", 250, 260), [0]);
    }

    #[test]
    fn metagene_exons_of_one_gene_are_reported_once() {
        let ann = annotation(".gtf", GTF, true);

        assert_eq!(ann.overlapping("chr1", 150, 350), [0]);
        assert!(ann.overlapping("chr1", 250, 260).is_empty());
    }
}
