"""Benchmark GPLANE factor construction/reuse against a real PLPROC executable.

Run from the repository root with the newly compiled pestutils library loaded:
python -m benchmarks.benchmark_gplane --plproc-exe /path/to/plproc --output results.json
"""

import argparse
import hashlib
import json
import os
import platform
import statistics
import subprocess
import tempfile
import time
from pathlib import Path

# Keep BLAS parallelism reproducible; PLPROC and the native kernel are serial.
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import numpy as np

from pypestutils.gplane import GPlane
from pypestutils.pestutilslib import PestUtilsLib
from tests.gplane_reference import prepare_plproc, run_plproc


def median_time(function, repeats):
    """Return median wall time and last result, with one untimed warm-up."""
    result = function()
    times = []
    for _ in range(repeats):
        start = time.perf_counter()
        result = function()
        times.append(time.perf_counter() - start)
    return statistics.median(times), result


def main():
    """Check executable parity, time each workload and save reproducible results."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plproc-exe", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sizes", type=int, nargs="+", default=[1000, 10000, 100000])
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    if args.repeats < 1 or any(n < 1 for n in args.sizes):
        parser.error("sizes and repeats must be positive")
    exe = args.plproc_exe.resolve(strict=True)
    lib = PestUtilsLib()
    library = Path(lib.pestutils._name)
    results = dict(
        platform=platform.platform(),
        processor=platform.processor(),
        python=platform.python_version(),
        numpy=np.__version__,
        repeats=args.repeats,
        plproc_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),
        library_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),
        timing_scope={
            "ppu_factor_apply_s": (
                "In-process factor construction and one application; excludes file I/O."
            ),
            "ppu_cached_apply_s": (
                "In-process application of precomputed factors; excludes file I/O."
            ),
            "plproc_process_s": (
                "Subprocess startup, input reads, geometry/interpolation and output "
                "write; inputs prepared before timing."
            ),
            "ppu_build_write_s": (
                "Python plane construction, native factors, application and output "
                "write; existing input arrays in memory."
            ),
        },
        cases=[],
    )
    rng = np.random.default_rng(20261109)
    with tempfile.TemporaryDirectory(prefix="ppu-gplane-benchmark-") as tmp:
        for bounded in [False, True]:
            boundary = (
                [[0, 0], [700, 0], [700, 400], [400, 400], [400, 1000], [0, 1000]]
                if bounded
                else None
            )
            geometry = dict(
                points=[[0, 0, 20], [1000, 0, -80], [0, 1000, 20]],
                grid=[[0, 0], [1000, 0], [0, 1000]],
                shape=(10, 10),
                spacing=(100.0, 100.0),
                boundary=boundary,
            )
            plane = GPlane(**geometry)
            values = np.exp(rng.normal(2, 0.5, (10, 10)))
            settings = dict(conwidth=2.0, a=30.0, anis_ratio=2.0, anis_bearing=30.0)
            apply_settings = dict(
                action="blend", interp_transform="log", transform="log"
            )
            for n in args.sizes:
                work = Path(tmp) / f"{bounded}-{n}"
                uv = rng.uniform([-200, -200], [1200, 1200], (n, 2))
                targets = (
                    plane.to_global(uv) + rng.uniform(-80, 80, (n, 1)) * plane.normal
                )
                background = np.full(n, 5.0)
                prepare_plproc(
                    work,
                    plane,
                    targets,
                    values,
                    background,
                    **settings,
                    **apply_settings,
                )
                factors = plane.calc_factors(targets, lib=lib, **settings)
                actual = factors.apply(values, background, **apply_settings)
                expected = run_plproc(exe, work)
                np.testing.assert_allclose(actual, expected, rtol=1.0e-6, atol=5.0e-7)

                def factor_apply():
                    return plane.calc_factors(targets, lib=lib, **settings).apply(
                        values, background, **apply_settings
                    )

                def build_write():
                    p = GPlane(**geometry)
                    out = p.calc_factors(targets, lib=lib, **settings).apply(
                        values, background, **apply_settings
                    )
                    np.savetxt(work / "ppu_result.dat", out, fmt="%.7g")

                def plproc_process():
                    result = subprocess.run(
                        [str(exe), "run.plproc"],
                        cwd=work,
                        capture_output=True,
                        text=True,
                        errors="replace",
                        timeout=180,
                        check=True,
                    )
                    if "no more commands to process" not in result.stdout.lower():
                        raise RuntimeError(result.stdout + result.stderr)

                build_s, _ = median_time(factor_apply, args.repeats)
                cached_s, _ = median_time(
                    lambda: factors.apply(values, background, **apply_settings),
                    args.repeats,
                )
                total_s, _ = median_time(build_write, args.repeats)
                reference_s, _ = median_time(plproc_process, args.repeats)
                case = dict(
                    targets=n,
                    sources=len(plane.source_uv),
                    bounded=bounded,
                    max_absolute_error=float(np.max(np.abs(actual - expected))),
                    ppu_factor_apply_s=build_s,
                    ppu_cached_apply_s=cached_s,
                    ppu_build_write_s=total_s,
                    plproc_process_s=reference_s,
                    factor_memory_bytes=factors.weights.nbytes
                    + factors.influence.nbytes,
                )
                results["cases"].append(case)
                print(json.dumps(case), flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
