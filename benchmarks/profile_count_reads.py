#!/usr/bin/env python
"""Profile / benchmark the scCountReads backend (`_sincei.count_bins`).

Runs the Rust counting backend on a BAM file and reports a per-phase
breakdown (setup / parallel counting / merge / CSR build / AnnData write)
plus total wall-clock time.

The per-phase numbers come from the library's opt-in instrumentation: this
script sets ``SINCEI_PROFILE=1`` and the backend prints
``[profile] <phase> <seconds>`` lines to stderr.

Build the backend in release mode first, otherwise the numbers are
meaningless:

    maturin develop --release

Usage:

    python benchmarks/profile_count_reads.py \
        --bam tests/count_testdata/test-1.bam \
        --barcodes tests/count_testdata/maya_384NLA.bc \
        --bin-size 10000 --threads 0 --compression none --repeats 1
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
import time


def read_barcodes(path: str) -> list[str]:
    with open(path) as handle:
        return [line.strip() for line in handle if line.strip()]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", default="tests/count_testdata/test-1.bam")
    parser.add_argument("--barcodes", default="tests/count_testdata/maya_384NLA.bc")
    parser.add_argument("--bin-size", type=int, default=10_000)
    parser.add_argument("--step-size", type=int, default=None, help="defaults to --bin-size")
    parser.add_argument("--threads", type=int, default=0, help="0 = all cores")
    parser.add_argument("--chunk-size", type=int, default=1_000_000)
    parser.add_argument("--compression", choices=["none", "gzip"], default="none")
    parser.add_argument("--compression-level", type=int, default=4)
    parser.add_argument("--cell-tag", default="BC")
    parser.add_argument("--repeats", type=int, default=1)
    args = parser.parse_args()

    # Enable the backend's per-phase timing (prints to stderr).
    os.environ["SINCEI_PROFILE"] = "1"

    from sincei import _sincei as internal

    barcodes = read_barcodes(args.barcodes)
    step_size = args.step_size or args.bin_size

    print(
        f"bam={args.bam}  barcodes={len(barcodes)}  bin_size={args.bin_size}  "
        f"threads={args.threads or 'all'}  compression={args.compression}",
        file=sys.stderr,
    )

    totals: list[float] = []
    for i in range(args.repeats):
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "counts.h5ad")
            t0 = time.perf_counter()
            internal.count_bins(
                [args.bam],
                barcodes,
                out,
                bin_size=args.bin_size,
                step_size=step_size,
                bc_tag=args.cell_tag,
                compression=args.compression,
                compression_level=args.compression_level,
                num_threads=args.threads,
                chunk_size=args.chunk_size,
            )
            dt = time.perf_counter() - t0
            size_mb = os.path.getsize(out) / 1e6
        totals.append(dt)
        print(f"[total]   run {i + 1}/{args.repeats}      {dt:8.3f} s   (output {size_mb:.1f} MB)", file=sys.stderr)

    if args.repeats > 1:
        best = min(totals)
        avg = sum(totals) / len(totals)
        print(f"[total]   best {best:.3f} s | avg {avg:.3f} s", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
