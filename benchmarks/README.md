# GPLANE validation and benchmark

Measured on 2026-11-09, Windows 11 x64 (build 26200), Intel Family 6 Model 183
Stepping 1, Python 3.11.15, NumPy 2.4.6, GNU Fortran 16.2 with Meson
optimization level 2. Numerical thread counts were fixed at one. The reference
is the real Watermark **PLPROC 5.0** executable. Executable and native-library
SHA-256 hashes and unrounded measurements are in
[results/gplane_windows.json](results/gplane_windows.json).

## Validation

The core suite, excluding notebook execution, passed **155 tests**, with eight
existing skips. Of these, **36 tests invoke PLPROC** independently. They cover
horizontal, dipping, and vertical planes; all three reference coordinate pairs;
bounded/unbounded features; concave limits; anisotropy; varying IDW powers;
four actions; separate log/linear interpolation and influence transforms;
inactive targets; exact hits; moved planes; and left-handed grids. Other tests
check analytic projection/decay, source ordering, translated map coordinates,
chunk consistency, saved factors/spec round trips, and invalid inputs.

Each executable comparison writes geometry, targets, source values and a PLPROC
command file, runs PLPROC, requires its successful-completion marker, then reads
its property output. There is no substitute Python oracle when PLPROC is
missing: these tests skip explicitly. PLPROC's column writer emits about seven
significant figures, so comparisons use `rtol=1e-6, atol=5e-7`.

The six larger benchmark cases also pass that tolerance before timing. Their
maximum absolute differences are 3.90e-6–6.23e-6 for the generated property
values, consistent with PLPROC's text-output precision. This establishes
agreement at that precision, rather than bitwise agreement inside PLPROC.

## Results

Each case uses a dipping 10-by-10 plane grid (100 source values), anisotropy
ratio 2, bearing 30 degrees, `conwidth=2`, `a=30`, log interpolation and log
blending. Bounded cases have a concave polygon. Targets span both sides of the
plane and the polygon limits. Data use a fixed random seed. Times are median
wall seconds over three repeats, following one warm-up per operation.

| Targets | Polygon | PPU factors + apply | PPU cached apply | PPU build + write | PLPROC process | Factors MB |
| ---: | :--- | ---: | ---: | ---: | ---: | ---: |
| 1,000 | No | 0.00210 | 0.000087 | 0.00445 | 0.02970 | 0.808 |
| 10,000 | No | 0.02206 | 0.000598 | 0.03700 | 0.15010 | 8.08 |
| 100,000 | No | 0.19533 | 0.008572 | 0.32703 | 1.22486 | 80.8 |
| 1,000 | Yes | 0.00238 | 0.000082 | 0.00557 | 0.03037 | 0.808 |
| 10,000 | Yes | 0.02074 | 0.000607 | 0.03624 | 0.14533 | 8.08 |
| 100,000 | Yes | 0.21016 | 0.008901 | 0.34506 | 1.14423 | 80.8 |

Timing scopes differ:

- **PPU factors + apply:** in-process native factor construction and one NumPy
  application; excludes input/output files and plane construction.
- **PPU cached apply:** applies already computed factors; excludes files and
  factor loading. This represents repeated value updates for fixed geometry.
- **PPU build + write:** Python plane construction, native factors, application,
  and writing results at seven significant figures. Input arrays are in memory.
- **PLPROC process:** executable startup, reading prewritten inputs, geometry,
  interpolation, influence application, and writing results. Input preparation
  and parsing results in Python are outside the timer.

At 100,000 targets, the PPU build-and-write path is approximately 3.3–3.7 times
faster than the measured PLPROC process. This is a workflow comparison with
different input paths, not an isolated interpolation-kernel speedup. Results
apply to this Windows machine and these data sizes; no Linux/macOS performance
claim is made. Storage grows with targets times sources; chunk large target
arrays as described in [the API guide](../docs/gplane.md).

## Reproduce

Build and install this branch with the tools in [BUILD.md](../BUILD.md). Run
from the repository root with the newly compiled library discoverable. Set
`PLPROC_EXE` to a locally available executable; it is not redistributed here.

```sh
export PLPROC_EXE=/path/to/plproc
python -m pytest tests --ignore=tests/test_notebooks.py -q
python -m benchmarks.benchmark_gplane \
    --plproc-exe "$PLPROC_EXE" \
    --sizes 1000 10000 100000 --repeats 3 \
    --output benchmarks/results/gplane_local.json
```

PowerShell uses `$env:PLPROC_EXE = 'C:\path\to\plproc.exe'` and
`--plproc-exe $env:PLPROC_EXE`. Activate the Python/Conda environment so its
numerical runtime libraries are available. `PESTUTILS_LIBRARY` may be set to
the directory containing this branch's newly built shared library. PLPROC
comparisons are optional for contributors without the executable; analytic
regression tests still run. Notebook execution is excluded from the reported
local validation.
