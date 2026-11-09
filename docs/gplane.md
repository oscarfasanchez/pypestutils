# 3D GPLANE structural overlays

The `codex/gplane-3d` development branch adds PLPROC-style GPLANE overlays to
the fork. Its Python version is `0.3.2.dev0+gplane.1`. The existing
`v0.3.1-usg.1` release wheels do **not** contain this feature. Build this branch's
native library using [BUILD.md](../BUILD.md) before using the API; installing
the Python source alone does not compile Fortran.

A GPLANE defines a planar structural feature in global XYZ coordinates. Values
on its regular local grid are interpolated onto the plane, then applied to
arbitrary target XYZ coordinates with distance-dependent influence. Targets can
be cell centres from DIS, DISV, DISU, MODFLOW-USG, or a point cloud. No model
connectivity is required.

## Example

```python
import numpy as np
from pypestutils.gplane import GPlane, GPlaneFactors

plane = GPlane(
    points=[[0, 0, 0], [1000, 0, -200], [0, 1000, 0]],
    grid=[[0, 0], [1000, 0], [0, 1000]],
    shape=(10, 10),              # nx, ny
    spacing=(100, 100),          # local dx, dy; model length units
    reference="xy",
    boundary=[[0, 0], [1000, 0], [1000, 1000], [0, 1000]],
)

# Replace these with the model's cell-centre XYZ and background property.
xyz = np.array([[200, 300, -40], [200, 300, -65], [1200, 300, -240]])
background_k = np.full(len(xyz), 10.0)
plane_k = np.full((10, 10), 0.1)  # shape (ny, nx); x changes fastest

factors = plane.calc_factors(
    xyz, conwidth=4, a=20,
    inverse_power=2, anis_ratio=2, anis_bearing=30,
    active=np.ones(len(xyz)),
)
k = factors.apply(
    plane_k, background_k, action="blend_never_higher",
    interp_transform="log", transform="log",
)

# Reuse for new property values while geometry/settings remain fixed.
factors.save("fault_factors.npz")
loaded = GPlaneFactors.load("fault_factors.npz")
k_next = loaded.apply(0.2, background_k, transform="log")

# Exchange geometry with PLPROC. Property values are stored separately.
plane.write_spec("fault.spec")
same_plane = GPlane.from_spec("fault.spec")
```

`PestUtilsLib.calc_gplane_factors(plane, xyz, **settings)` is also available.
Projection and factor construction use the exported native
`calc_gplane_factors` routine; property application uses NumPy.

## Geometry and influence

Three noncollinear XYZ `points` define the plane. `reference` chooses which two
global coordinates specify the `grid` and polygon vertices: `xy`, `xz`, or `yz`.
Use `xz` or `yz` when `xy` cannot resolve a vertical plane. The three `grid`
points define the origin, positive local x direction, and approximate positive
local y direction. The latter is orthogonalized. Left-handed grids are allowed.
Source locations are cell centres, starting at `(dx/2, dy/2)`, with x varying
fastest. Source values may be scalar, a flat vector, or a `(ny, nx)` array.

Targets project orthogonally onto the plane. An optional simple polygon of
3–10 vertices limits the feature; omit the repeated closing vertex. Projections
outside the polygon snap to its nearest edge for interpolation. The influence
for distance `d` is:

```text
f(d) = 1                                      if d <= conwidth / 2
       exp(-((d - conwidth / 2) / a)**2)       otherwise
```

`conwidth` is the **full** constant-width core. Both `conwidth` and positive `a`
use model length units. Outside polygon limits, the perpendicular and lateral
distance influences multiply. Very small tails (`(d-conwidth/2)/a >= 27`) are
treated as zero. Zero entries in the target `active` mask preserve background
values.

Anisotropy bearing is in degrees clockwise from local +y toward local +x.
The ratio scales distances across this direction. For PLPROC 5 compatibility,
the kernel uses its executable-observed IDW rule:

```text
weight ~ (along**4 + across**4)**(-inverse_power / 2)
```

Here `across` includes anisotropy scaling. This differs from Euclidean IDW.
The convention was established by perturbing individual source values in the
PLPROC executable and is covered by a fixed regression and executable tests.
Coincident source points receive unit weight within the numerical tolerance.

## Actions and transforms

For background `b`, interpolated plane value `p`, and influence `f`:

| Action | Result |
| --- | --- |
| `add` | `b + f*p` |
| `blend` | `b + f*(p-b)` |
| `blend_never_higher` | `min(b, b + f*(p-b))` |
| `blend_never_lower` | `max(b, b + f*(p-b))` |

`interp_transform="log"` interpolates base-10 logs of source values.
Independently, `transform="log"` performs the action in base-10 log space and
back-transforms the result. Log operations require positive participating
values. In log space, `add` therefore multiplies the background by `p**f`.
The default for both transforms is `"none"`.

## PEST / pyEMU use

For fixed geometry, save factors during interface setup. A worker's forward
script can read a templated plane-value file and apply the cached factors:

```python
def apply_fault():
    import numpy as np
    from pypestutils.gplane import GPlaneFactors

    values = np.loadtxt("fault_values.dat")
    background = np.loadtxt("background_k.dat")
    factors = GPlaneFactors.load("fault_factors.npz")
    k = factors.apply(values, background, transform="log")
    np.savetxt("model_k.dat", k)
```

Template `fault_values.dat` with `PstFrom.add_parameters(..., par_type="grid",
par_style="direct", transform="log")`, using that file's array dimensions;
copy the factor/background files into the worker template. Register
`apply_fault()` with
`pf.add_py_function("fault_preprocess.py", "apply_fault()", is_pre_cmd=True)`
before the model runs. Save the function in `fault_preprocess.py` for that call.
`model_k.dat` must use the array or node ordering expected by the model.

If plane position, dip, polygon, grid, target coordinates, mask, width, decay,
or IDW settings are parameters, rebuild the `GPlane` and recalculate factors
inside each forward run. Saved factors do not track geometry movement. Multiple
planes can be applied sequentially; their order can affect the result.

## Scope and scaling

This implementation covers core 3D GPLANE projection, polygon limits, regular
source grids, anisotropic IDW, four actions, two independent transforms, target
masks, and factor caching. Spec readers/writers exchange geometry; optional
PLIST/SLIST sections are ignored on import. It does not implement PLPROC's
command language, source-selection expressions, VTK export, or report writers.

Factors are dense: approximately `8 * ntarget * (nsource + 1)` bytes. For
100,000 targets and 100 sources this is 80.8 MB, excluding temporary arrays.
For larger problems, calculate/apply target chunks and assemble the resulting
property vector. Each chunk produces the same result as the corresponding
part of a full calculation. Rebuild the native library together with the
Python wrapper; older released libraries do not export the new symbol.

See [benchmark and validation results](../benchmarks/README.md) for the
independent PLPROC comparisons, reproducible commands, and measured timings.
The reference behavior comes from the Watermark PLPROC manual, section 2.5,
`NEW_GPLANE_CLIST`, and `INFORM_FROM_GPLANE`, checked against PLPROC 5.0.
