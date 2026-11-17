# Changes in the pypestutils fork

## ppuplus integration branch

`ppuplus` is the combined branch for all fork features, with package version
`0.3.3+usg.1`. It merges `codex/gplane-3d` and `feat/fieldgen3d-sva2`, including
the already integrated `feat/disu-grids`, `feat/mod2obs-disu-layers`,
`feat/mfusg-grids`, and `perf/interp-factor-ring-search` changes. The existing
feature branches are preserved.

In addition to all capabilities below, `PestUtilsLib.fieldgen3d_sva2` computes
3D stochastic fields from user-supplied independent normal variates and source
volumes. Source points can differ from targets; anisotropy, rotations,
natural/log properties, per-target statistics, and an optional summation
cutoff are supported. Tests compare four fields with saved FIELDGEN3D_SVA2
utility outputs and verify variance, mean and variate scaling.

The combined release is `v0.3.3-usg.1`, including compiled wheels for all
features above. See the [README](README.md) for installation. The older
`v0.3.2-usg.1` wheels predate this integration and lack `fieldgen3d_sva2`.

Integration validation on Windows x64 / Python 3.11 (2026-11-17): a fresh
Fortran build and installed wheel passed the GPLANE/SVA2 wheel check and
164 non-notebook tests, including 36 live PLPROC comparisons and all four
SVA2 reference fields. Eight existing driver cases were skipped for missing
full input datasets. Notebook tests were excluded from this local run.
Python lint and formatting checks passed for the changed files.

## GPLANE release additions

`codex/gplane-3d` builds on the USG additions below and adds grid-independent
3D GPLANE overlays: arbitrary plane orientation, optional polygon limits,
anisotropic IDW, Gaussian influence, add/blend/clamped actions, independent
log transforms, target masks, and reusable factors. A native Fortran kernel
is exposed through `pypestutils.gplane.GPlane` and
`PestUtilsLib.calc_gplane_factors`. Geometry specifications can be exchanged
with PLPROC.

See [the API guide](docs/gplane.md) and
[PLPROC validation and benchmarks](benchmarks/README.md). GPLANE is included in
`v0.3.2-usg.1` (Python version `0.3.2+usg.1`). This covers the core overlay
workflow, not PLPROC's entire script language or source-selection machinery.

## Fork wheels and earlier USG additions

Prebuilt wheels and installation instructions are available in the
[v0.3.2-usg.1 release](https://github.com/oscarfasanchez/pypestutils/releases/tag/v0.3.2-usg.1).
Use [environment-usg.yml](environment-usg.yml) to create a Conda environment
that automatically downloads the appropriate wheel.

This note summarizes the additions on the `feat/mfusg-grids` branch compared
with upstream `develop` at `6721975` (`[pre-commit.ci] pre-commit autoupdate
(#65)`). It covers the fork work in this branch; it is not a list of changes
from every branch in the fork.

## MODFLOW 6 DISU grids

- Read unstructured DISU grids from MODFLOW 6 binary grid (`.grb`) files and
  install them for use by pypestutils.
- Associate nodes with model layers and provide DISU node centres.
- Interpolate points layer by layer, and let `mod2obs_mf6` use `node_layers`
  (model type 33) so an observation is evaluated against nodes in its layer.

## MODFLOW-USG grids

- Read legacy MODFLOW-USG DISU and GSF grid files, including `OPEN/CLOSE`
  paths that contain Windows backslashes or are quoted.
- Install a USG grid from DISU and GSF arrays, deriving layer assignments from
  `NODLAY`.
- Read single- and double-precision unstructured dependent-variable files and
  interpolate their values at points.
- Add `mod2obs_mfusg` and grid-reading helpers corresponding to the MODFLOW 6
  functionality.

## Layer-aware 2D helpers

The 2D grid information, realization generation, and pilot-point interpolation
helpers can work on one layer of a DISU grid using `node_layers`. Grid
information can also be read for one GSF layer, and input tables with a layer
column are filtered to the requested layer.

## Interpolation performance

`calc_mf6_interp_factors` now searches nearby cell-centre rings in increasing
size, from three centres through seven. The resulting interpolation factors
are unchanged; the changelog reports roughly a 30x speedup for Voronoi-cell
DISV grids.

## Python and Fortran interfaces and coverage

The exported Fortran interfaces and Python `ctypes` declarations were extended
for the new grid readers, installers, interpolation routines, and observation
mapping. Regression tests and small grid and head-file fixtures cover DISU,
USG, layer handling, and interpolation behavior.

See `CHANGELOG.md` for the public API details by version.
