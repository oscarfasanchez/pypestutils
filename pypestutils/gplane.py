"""Reusable, grid-independent 3D GPLANE overlays compatible with PLPROC 5.

Plane geometry is expressed in model length units. Target XYZ coordinates may
come from any grid. This module does not modify model files or model topology.
"""

from ctypes import byref, c_double, c_int
from dataclasses import dataclass
from pathlib import Path

import numpy as np


def _array(value, shape, name):
    a = np.array(value, dtype=np.float64, copy=True, order="C")
    if a.shape != shape or not np.isfinite(a).all():
        raise ValueError(f"{name} must be a finite array with shape {shape}")
    if np.any(np.abs(a) > 1.0e100):
        raise ValueError(f"{name} coordinates exceed the safe numerical range")
    return a


def _scalar(value, name, positive=True):
    value = float(value)
    if not np.isfinite(value) or (value <= 0 if positive else value < 0):
        raise ValueError(
            f"{name} must be finite and {'positive' if positive else 'nonnegative'}"
        )
    if value > 1.0e100:
        raise ValueError(f"{name} exceeds the safe numerical range")
    return value


def _validate_polygon(p, tol):
    """Reject repeated vertices, zero area and intersecting polygon edges."""
    n = len(p)
    if not 3 <= n <= 10:
        raise ValueError(
            "boundary must have 3 to 10 vertices, without a closing repeat"
        )
    for i in range(n):
        if np.any(np.linalg.norm(p[i + 1 :] - p[i], axis=1) <= tol):
            raise ValueError("boundary has repeated vertices")
    if (
        abs(np.sum(p[:, 0] * np.roll(p[:, 1], -1) - p[:, 1] * np.roll(p[:, 0], -1)))
        <= tol**2
    ):
        raise ValueError("boundary has zero area")

    def cross(a, b, c):
        d, e = b - a, c - a
        return d[0] * e[1] - d[1] * e[0]

    def on_segment(a, b, c):
        return (
            abs(cross(a, b, c)) <= tol * np.linalg.norm(b - a)
            and np.all(c >= np.minimum(a, b) - tol)
            and np.all(c <= np.maximum(a, b) + tol)
        )

    for i in range(n):
        a, b = p[i], p[(i + 1) % n]
        for j in range(i + 1, n):
            if j == (i + 1) % n or (j + 1) % n == i:
                continue
            c, d = p[j], p[(j + 1) % n]
            if (
                cross(a, b, c) * cross(a, b, d) < 0
                and cross(c, d, a) * cross(c, d, b) < 0
            ) or any(
                (
                    on_segment(a, b, c),
                    on_segment(a, b, d),
                    on_segment(c, d, a),
                    on_segment(c, d, b),
                )
            ):
                raise ValueError("boundary has intersecting edges")


@dataclass(frozen=True)
class GPlaneFactors:
    """Cached source weights and target influence for one fixed geometry.

    ``weights`` has shape (ntarget, nsource); ``influence`` has shape (ntarget,).
    Dense storage uses about ``8 * ntarget * (nsource + 1)`` bytes. Calculate
    factors in target chunks when this is too large. Geometry, masks, widths,
    decay, or interpolation-setting changes require recalculation.
    """

    weights: np.ndarray
    influence: np.ndarray
    source_shape: tuple | None = None

    def __post_init__(self):
        """Protect cached factors against accidental mutation."""
        self.weights.setflags(write=False)
        self.influence.setflags(write=False)

    def apply(
        self,
        source_values,
        background,
        *,
        action="blend",
        interp_transform="none",
        transform="none",
    ):
        """Apply new source values to a background without mutating either.

        Source values are a flat PLPROC-ordered vector or an array of shape
        ``(ny, nx)`` (x changes fastest). Actions are ``add``, ``blend``,
        ``blend_never_higher`` and ``blend_never_lower``. ``interp_transform``
        controls in-plane IDW; ``transform`` independently controls blending.
        ``log`` requires positive values and uses base-10 logarithms.
        """
        if action not in ("add", "blend", "blend_never_higher", "blend_never_lower"):
            raise ValueError("unknown GPLANE action")
        if interp_transform not in ("none", "log") or transform not in ("none", "log"):
            raise ValueError("transforms must be 'none' or 'log'")
        source = np.asarray(source_values, dtype=float)
        if source.ndim > 2 or (
            source.ndim == 2
            and self.source_shape is not None
            and source.shape != self.source_shape
        ):
            raise ValueError(
                f"source_values must be flat or have shape {self.source_shape}"
            )
        source = source.ravel()
        if source.size == 1:
            source = np.full(self.weights.shape[1], source.item())
        if source.shape != (self.weights.shape[1],) or not np.isfinite(source).all():
            raise ValueError("source_values must match the source grid and be finite")
        base = np.broadcast_to(
            np.asarray(background, dtype=float), self.influence.shape
        ).copy()
        if not np.isfinite(base).all():
            raise ValueError("background must be finite")
        use = self.influence > 0
        if not use.any():
            return base
        if interp_transform == "log":
            if np.any(source <= 0):
                raise ValueError("log interpolation requires positive source values")
            source = np.log10(source)
        feature = (self.weights @ source)[use]
        if interp_transform == "log":
            feature = np.power(10.0, feature)
        old = base[use].copy()
        if transform == "log":
            if np.any(old <= 0) or np.any(feature <= 0):
                raise ValueError(
                    "log blending requires positive source and background values"
                )
            old, feature = np.log10(old), np.log10(feature)
        change = feature if action == "add" else feature - old
        new = old + self.influence[use] * change
        if action == "blend_never_higher":
            new = np.minimum(new, old)
        elif action == "blend_never_lower":
            new = np.maximum(new, old)
        base[use] = np.power(10.0, new) if transform == "log" else new
        if not np.isfinite(base).all():
            raise ValueError("GPLANE application produced nonfinite values")
        return base

    def save(self, filename):
        """Save numeric factors to an NPZ file, with no Python object payload."""
        shape = () if self.source_shape is None else self.source_shape
        np.savez_compressed(
            filename,
            weights=self.weights,
            influence=self.influence,
            source_shape=np.array(shape, dtype=np.int64),
        )

    @classmethod
    def load(cls, filename):
        """Load and validate a factor file written by :meth:`save`."""
        with np.load(filename, allow_pickle=False) as data:
            if not {"weights", "influence", "source_shape"} <= set(data.files):
                raise ValueError("incomplete GPLANE factor file")
            w = np.array(data["weights"], dtype=float, order="C")
            f = np.array(data["influence"], dtype=float)
            shape = tuple(data["source_shape"].tolist())
        if (
            w.ndim != 2
            or min(w.shape) == 0
            or f.shape != (len(w),)
            or not np.isfinite(w).all()
            or not np.isfinite(f).all()
            or np.any(w < 0)
            or np.any((f < 0) | (f > 1))
            or not np.allclose(w[f > 0].sum(axis=1), 1)
        ):
            raise ValueError("invalid GPLANE factor file")
        if shape and (
            len(shape) != 2 or min(shape) < 1 or np.prod(shape) != w.shape[1]
        ):
            raise ValueError("invalid GPLANE source shape")
        return cls(w, f, shape or None)


class GPlane:
    """A PLPROC-style plane, regular source grid, and optional polygon limits.

    Parameters
    ----------
    points : array_like, shape (3, 3)
        Noncollinear global XYZ points defining the plane.
    grid : array_like, shape (3, 2)
        Global reference-coordinate pairs defining the grid origin, positive
        local x direction, and approximate positive local y direction.
    shape : (int, int)
        Number of source cells (nx, ny). Values flatten with x changing fastest.
    spacing : (float, float)
        Positive local cell widths (dx, dy), in model length units.
    reference : {'xy', 'xz', 'yz'}
        Which two global coordinates define grid and boundary points. Choose
        xz or yz for vertical planes.
    boundary : array_like, shape (n, 2), optional
        Simple polygon of 3 to 10 reference-coordinate vertices. No closing
        repeat. Outside projections snap to the nearest edge; lateral and
        perpendicular Gaussian influence factors multiply.
    tolerance : float
        Positive coincidence tolerance in model length units.
    """

    def __init__(
        self,
        points,
        grid,
        shape,
        spacing,
        *,
        reference="xy",
        boundary=None,
        tolerance=1.0e-9,
    ):
        self.points = _array(points, (3, 3), "points")
        self.grid = _array(grid, (3, 2), "grid")
        self.tolerance = _scalar(tolerance, "tolerance")
        if reference not in ("xy", "xz", "yz"):
            raise ValueError("reference must be 'xy', 'xz', or 'yz'")
        self.reference = reference
        if len(shape) != 2 or any(int(n) != n or n < 1 for n in shape):
            raise ValueError("shape must contain two positive integers (nx, ny)")
        self.shape = tuple(int(n) for n in shape)
        if np.prod(self.shape, dtype=object) > np.iinfo(np.int32).max:
            raise ValueError("source grid exceeds the native integer capacity")
        if len(spacing) != 2:
            raise ValueError("spacing must contain dx and dy")
        self.spacing = tuple(_scalar(s, "spacing") for s in spacing)
        ab, ac = self.points[1:] - self.points[0]
        norm_ab = np.linalg.norm(ab / np.max(np.abs(ab))) if np.any(ab) else 0
        norm_ac = np.linalg.norm(ac / np.max(np.abs(ac))) if np.any(ac) else 0
        if not norm_ab or not norm_ac:
            raise ValueError("plane points are coincident")
        ab_unit = (ab / np.max(np.abs(ab))) / norm_ab
        ac_unit = (ac / np.max(np.abs(ac))) / norm_ac
        normal = np.cross(ab_unit, ac_unit)
        length = np.linalg.norm(normal)
        if (
            np.linalg.norm(ab) <= self.tolerance
            or np.linalg.norm(ac) <= self.tolerance
            or length <= 1.0e-12
        ):
            raise ValueError("plane points are coincident or nearly collinear")
        self.normal = np.ascontiguousarray(normal / length)
        g = self._lift(self.grid)
        self.origin = np.ascontiguousarray(g[0])
        ex = g[1] - g[0]
        exlen = np.linalg.norm(ex)
        if exlen <= self.tolerance:
            raise ValueError("grid x direction is degenerate")
        self.xaxis = np.ascontiguousarray(ex / exlen)
        ey = g[2] - g[0]
        ey -= np.dot(ey, self.xaxis) * self.xaxis
        if np.linalg.norm(ey) <= self.tolerance:
            raise ValueError("grid y direction is degenerate")
        self.yaxis = np.ascontiguousarray(ey / np.linalg.norm(ey))
        self.boundary = None
        self.boundary_uv = np.empty((0, 2))
        if boundary is not None:
            p = np.asarray(boundary, dtype=float)
            if p.ndim != 2 or p.shape[1] != 2 or not np.isfinite(p).all():
                raise ValueError("boundary must be a finite (n, 2) array")
            self.boundary = p.copy()
            if np.any(np.abs(p) > 1.0e100):
                raise ValueError("boundary coordinates exceed the safe numerical range")
            self.boundary_uv = self.to_local(self._lift(p))
            _validate_polygon(self.boundary_uv, self.tolerance)
        nx, ny = self.shape
        u, v = np.meshgrid(
            (np.arange(nx) + 0.5) * self.spacing[0],
            (np.arange(ny) + 0.5) * self.spacing[1],
        )
        self.source_uv = np.column_stack((u.ravel(), v.ravel()))
        for array in (
            self.points,
            self.grid,
            self.normal,
            self.origin,
            self.xaxis,
            self.yaxis,
            self.boundary_uv,
            self.source_uv,
            self.boundary,
        ):
            if array is not None:
                array.setflags(write=False)

    def _lift(self, pairs):
        axes = ["xyz".index(c) for c in self.reference]
        missing = next(i for i in range(3) if i not in axes)
        if abs(self.normal[missing]) < 1.0e-12:
            raise ValueError("reference coordinates cannot resolve this plane")
        xyz = np.zeros((len(pairs), 3))
        xyz[:, axes] = pairs
        # Centre differences before solving to preserve precision at map offsets.
        xyz[:, missing] = (
            self.points[0, missing]
            - ((pairs - self.points[0, axes]) @ self.normal[axes])
            / self.normal[missing]
        )
        return xyz

    def to_local(self, xyz):
        """Orthogonally project XYZ points to local (u, v) coordinates."""
        delta = np.asarray(xyz, dtype=float) - self.origin
        return np.column_stack((delta @ self.xaxis, delta @ self.yaxis))

    def to_global(self, uv):
        """Convert local (u, v) coordinates to global XYZ points on the plane."""
        uv = np.asarray(uv, dtype=float)
        return self.origin + uv[:, :1] * self.xaxis + uv[:, 1:] * self.yaxis

    def calc_factors(
        self,
        targets,
        *,
        conwidth,
        a,
        inverse_power=2.0,
        anis_ratio=1.0,
        anis_bearing=0.0,
        active=1,
        lib=None,
    ):
        """Calculate reusable factors using the compiled GPLANE kernel.

        ``targets`` has shape (ntarget, 3). ``conwidth`` is the full constant
        core width; ``a`` is the Gaussian decay length. Bearing is in degrees
        clockwise from local +v toward +u. PLPROC 5's empirically verified IDW
        convention is ``(du**4 + dv**4)**(-inverse_power/2)`` after anisotropic
        rotation/scaling. This deliberately preserves PLPROC's actual weights.
        ``active`` is a scalar or target-length mask; zero means unchanged.
        """
        from .pestutilslib import PestUtilsLib, PestUtilsLibError

        targets = np.asarray(targets, dtype=float)
        if (
            targets.ndim != 2
            or targets.shape[1] != 3
            or len(targets) == 0
            or not np.isfinite(targets).all()
        ):
            raise ValueError("targets must be a finite nonempty (n, 3) array")
        if len(targets) > np.iinfo(np.int32).max or np.any(np.abs(targets) > 1.0e100):
            raise ValueError("targets exceed native dimensions or coordinate range")
        settings = [
            _scalar(conwidth, "conwidth", False),
            _scalar(a, "a"),
            _scalar(inverse_power, "inverse_power"),
            _scalar(anis_ratio, "anis_ratio"),
        ]
        if not np.isfinite(anis_bearing):
            raise ValueError("anis_bearing must be finite")
        settings += [float(anis_bearing) % 360.0, self.tolerance]
        active = np.broadcast_to(np.asarray(active, dtype=float), len(targets))
        if not np.isfinite(active).all():
            raise ValueError("active mask must be finite")
        mask = np.ascontiguousarray(active != 0, dtype=np.int32)
        vectors = [np.ascontiguousarray(self.source_uv[:, i]) for i in range(2)]
        vectors += [np.ascontiguousarray(targets[:, i]) for i in range(3)]
        bounds = [np.ascontiguousarray(self.boundary_uv[:, i]) for i in range(2)]
        weights = np.zeros((len(targets), len(self.source_uv)))
        influence = np.zeros(len(targets))
        lib = PestUtilsLib() if lib is None else lib
        result = lib.pestutils.calc_gplane_factors(
            byref(c_int(len(self.source_uv))),
            *vectors[:2],
            byref(c_int(len(targets))),
            *vectors[2:],
            self.origin,
            self.xaxis,
            self.yaxis,
            self.normal,
            byref(c_int(len(self.boundary_uv))),
            *bounds,
            mask,
            *(byref(c_double(s)) for s in settings),
            weights,
            influence,
        )
        if result:
            raise PestUtilsLibError(lib.retrieve_error_message())
        return GPlaneFactors(weights, influence, self.shape[::-1])

    def write_spec(self, filename):
        """Write geometry as a PLPROC GPLANE specification (no PLIST values)."""

        def fmt(rows):
            return "\n".join(" ".join(f"{v:.17g}" for v in row) for row in rows)

        boundary = (
            "0"
            if self.boundary is None
            else f"{len(self.boundary)}\n{fmt(self.boundary)}"
        )
        text = (
            f"* gplane definition\n{fmt(self.points)}\n* gplane reference coordinates\n"
            f"{' '.join(self.reference)}\n* gplane limits\n{boundary}\n"
            f"* gplane grid\n{fmt(self.grid)}\n"
            f"{self.shape[0]} {self.spacing[0]:.17g}\n"
            f"{self.shape[1]} {self.spacing[1]:.17g}\n"
        )
        Path(filename).write_text(text, encoding="ascii")

    @classmethod
    def from_spec(cls, filename, **kwargs):
        """Read PLPROC geometry sections; optional PLIST/SLIST data are ignored."""
        sections = {}
        current = None
        for line in Path(filename).read_text().splitlines():
            line = line.split("#", 1)[0].strip()
            if line.startswith("*"):
                current = line[1:].strip().lower()
                sections[current] = []
            elif line and current:
                sections[current].append(line.replace(",", " ").split())
        try:
            points = np.array(sections["gplane definition"], dtype=float)
            reference = "".join(sections["gplane reference coordinates"][0]).lower()
            limits = sections["gplane limits"]
            n = int(limits[0][0])
            if len(limits) != n + 1:
                raise ValueError("boundary vertex count does not match specification")
            boundary = np.array(limits[1:], dtype=float) if n else None
            grid = sections["gplane grid"]
            if len(grid) != 5:
                raise ValueError("gplane grid section must have five lines")
            return cls(
                points,
                np.array(grid[:3], dtype=float),
                (int(grid[3][0]), int(grid[4][0])),
                (float(grid[3][1]), float(grid[4][1])),
                reference=reference,
                boundary=boundary,
                **kwargs,
            )
        except (KeyError, IndexError) as error:
            raise ValueError("incomplete GPLANE specification") from error
