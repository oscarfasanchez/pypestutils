# Changes in the pypestutils fork

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

See `CHANGELOG.md` under **Unreleased** for the public API details.
