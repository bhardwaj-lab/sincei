# Benchmarks

Profiling / benchmarking drivers for the sincei Rust backend. These are
standalone scripts — they are not part of the published package.

## Benchmarking `scCountReads`

**Always build in release mode first** — a debug build is many times slower
and distorts the relative breakdown:

```bash
maturin develop --release
```

Then:

```bash
python benchmarks/profile_count_reads.py \
    --bam tests/count_testdata/test-1.bam \
    --barcodes tests/count_testdata/maya_384NLA.bc \
    --bin-size 10000 --threads 0 --compression none --repeats 3
```

Run `python benchmarks/profile_count_reads.py --help` for all options
(thread count, chunk size, compression type/level, bin/step size, repeats).
