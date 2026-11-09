"""Analytic, validation, round-trip and real-executable GPLANE regression tests."""

import os

import numpy as np
import pytest

from pypestutils.gplane import GPlane, GPlaneFactors
from pypestutils.pestutilslib import PestUtilsLib

from .gplane_reference import prepare_plproc, run_plproc


def make_plane(kind="horizontal", bounded=False):
    if kind == "vertical":
        points = [[0, 0, 0], [0, 4, 0], [0, 0, 3]]
        reference = "yz"
    elif kind == "xz":
        points = [[0, 0, 0], [4, 0, 0], [0, 0, 3]]
        reference = "xz"
    else:
        points = [[0, 0, 0], [4, 0, -2 if kind == "dipping" else 0], [0, 3, 0]]
        reference = "xy"
    boundary = [[0, 0], [3, 0], [3, 1], [1, 1], [1, 2], [0, 2]] if bounded else None
    return GPlane(
        points,
        [[0, 0], [4, 0], [0, 3]],
        (4, 3),
        (1.0, 1.0),
        reference=reference,
        boundary=boundary,
    )


def test_projection_source_order_and_exact_hits():
    p = make_plane("dipping")
    assert np.allclose(
        p.source_uv[:5], [[0.5, 0.5], [1.5, 0.5], [2.5, 0.5], [3.5, 0.5], [0.5, 1.5]]
    )
    targets = p.to_global(p.source_uv)
    factors = PestUtilsLib().calc_gplane_factors(p, targets, conwidth=0, a=1)
    np.testing.assert_allclose(factors.weights, np.eye(12), atol=1.0e-12)
    np.testing.assert_allclose(factors.apply(np.arange(1, 13), 5), np.arange(1, 13))
    np.testing.assert_allclose(p.to_local(targets), p.source_uv)


@pytest.mark.parametrize(
    "action", ["add", "blend", "blend_never_higher", "blend_never_lower"]
)
@pytest.mark.parametrize("transform", ["none", "log"])
def test_decay_actions_and_masks(action, transform):
    p = make_plane()
    targets = [[1, 1, 0], [1, 1, 1], [1, 1, 3], [1, 1, 500], [1, 1, 0]]
    f = p.calc_factors(targets, conwidth=2, a=2, active=[1, 1, 1, 1, 0])
    np.testing.assert_allclose(f.influence, [1, 1, np.exp(-1), 0, 0])
    base, value = 5.0, 2.0
    old, feature = (
        (np.log10(base), np.log10(value)) if transform == "log" else (base, value)
    )
    expected = old + f.influence[:4] * (feature if action == "add" else feature - old)
    if action == "blend_never_higher":
        expected = np.minimum(expected, old)
    if action == "blend_never_lower":
        expected = np.maximum(expected, old)
    if transform == "log":
        expected = 10**expected
    np.testing.assert_allclose(
        f.apply(value, base, action=action, transform=transform), np.r_[expected, base]
    )


def test_boundary_snap_and_multiplicative_decay():
    p = GPlane(
        [[0, 0, 0], [4, 0, 0], [0, 4, 0]],
        [[0, 0], [4, 0], [0, 4]],
        (4, 4),
        (1, 1),
        boundary=[[0, 0], [2, 0], [2, 2], [0, 2]],
    )
    targets = [[2, 1, 2], [4, 1, 2], [2, 1, 0], [4, 1, 0]]
    f = p.calc_factors(targets, conwidth=0, a=2)
    np.testing.assert_allclose(f.influence, [np.exp(-1), np.exp(-2), 1, np.exp(-1)])
    np.testing.assert_allclose(f.weights[0], f.weights[1])
    np.testing.assert_allclose(f.weights[0], f.weights[2])


def test_plproc_component_distance_convention():
    p = GPlane(
        [[0, 0, 0], [2, 0, 0], [0, 3, 0]], [[0, 0], [2, 0], [0, 3]], (2, 3), (1, 1)
    )
    f = p.calc_factors([[0, 0, 0]], conwidth=0, a=1)
    weights = 1 / np.sum(p.source_uv**4, axis=1)
    np.testing.assert_allclose(f.weights[0], weights / weights.sum())
    # Measured independently from PLPROC 5 by perturbing each source in turn.
    np.testing.assert_allclose(
        f.weights[0],
        [0.9370719, 0.02285541, 0.02285541, 0.01156879, 0.00299384, 0.00265459],
        atol=5.0e-8,
    )


def test_log_interpolation_and_roundtrips(tmp_path):
    p = make_plane("vertical", True)
    p.write_spec(tmp_path / "plane.spec")
    other = GPlane.from_spec(tmp_path / "plane.spec")
    np.testing.assert_allclose(other.source_uv, p.source_uv)
    np.testing.assert_allclose(other.boundary_uv, p.boundary_uv)
    targets = p.to_global([[1, 1], [2, 1], [3, 2]])
    f = other.calc_factors(targets, conwidth=0, a=1)
    f.save(tmp_path / "factors.npz")
    loaded = GPlaneFactors.load(tmp_path / "factors.npz")
    source = np.arange(1, 13)
    expected = 10 ** (f.weights @ np.log10(source))
    expected = 5 + f.influence * (expected - 5)
    np.testing.assert_allclose(
        loaded.apply(source, 5, interp_transform="log"), expected
    )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"points": [[0, 0, 0]] * 3},
        {"points": [[0, 0, 0], [1, 1, 1], [2, 2, 2]]},
        {"grid": [[0, 0], [0, 0], [0, 1]]},
        {"grid": [[0, 0], [1, 0], [2, 0]]},
        {"shape": (0, 2)},
        {"shape": (2.5, 2)},
        {"spacing": (-1, 1)},
        {"boundary": [[0, 0], [1, 0], [0, 1], [0, 0]]},
        {"boundary": [[0, 0], [3, 0], [0, 2], [3, 3]]},
        {"boundary": [[0, 0], [1, 0], [2, 0]]},
    ],
)
def test_invalid_geometry(kwargs):
    arguments = dict(
        points=[[0, 0, 0], [4, 0, 0], [0, 3, 0]],
        grid=[[0, 0], [4, 0], [0, 3]],
        shape=(4, 3),
        spacing=(1, 1),
    )
    arguments.update(kwargs)
    with pytest.raises(ValueError):
        GPlane(**arguments)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"a": 0},
        {"conwidth": -1},
        {"inverse_power": 0},
        {"anis_ratio": 0},
        {"anis_bearing": np.nan},
    ],
)
def test_invalid_settings(kwargs):
    settings = dict(conwidth=0, a=1)
    settings.update(kwargs)
    with pytest.raises(ValueError):
        make_plane().calc_factors([[0, 0, 0]], **settings)


def test_invalid_values():
    f = make_plane().calc_factors([[0, 0, 0]], conwidth=0, a=1)
    for kwargs in [
        {"action": "unknown"},
        {"transform": "bad"},
        {"interp_transform": "bad"},
    ]:
        with pytest.raises(ValueError):
            f.apply(1, 1, **kwargs)
    with pytest.raises(ValueError):
        f.apply(-1, 1, interp_transform="log")
    with pytest.raises(ValueError):
        f.apply(1, 0, transform="log")


def test_geometry_movement_offsets_and_chunking():
    p = make_plane("dipping", True)
    targets = p.to_global([[1, 1], [2, 2], [4, 4]]) + p.normal * 2
    f = p.calc_factors(targets, conwidth=1, a=2)
    offset = np.array([500000.0, 6500000.0, 100.0])
    # Rebuild the plane when geometry moves; fixed factors cannot track movement.
    moved = GPlane(
        p.points + offset,
        p.grid + offset[:2],
        p.shape,
        p.spacing,
        boundary=p.boundary + offset[:2],
    )
    shifted = moved.calc_factors(targets + offset, conwidth=1, a=2)
    np.testing.assert_allclose(f.weights, shifted.weights, atol=1.0e-8)
    np.testing.assert_allclose(f.influence, shifted.influence, atol=1.0e-8)
    chunks = [p.calc_factors(t[None], conwidth=1, a=2) for t in targets]
    np.testing.assert_allclose(np.vstack([c.weights for c in chunks]), f.weights)
    np.testing.assert_allclose(
        np.concatenate([c.influence for c in chunks]), f.influence
    )
    with pytest.raises(ValueError):
        p.points[0, 0] = 100


def test_left_handed_grid_and_reversed_boundary():
    p = GPlane(
        [[0, 0, 0], [4, 0, 0], [0, 3, 0]],
        [[0, 0], [4, 0], [0, -3]],
        (4, 3),
        (1, 1),
        boundary=[[0, 0], [4, 0], [4, -3], [0, -3]],
    )
    other = GPlane(p.points, p.grid, p.shape, p.spacing, boundary=p.boundary[::-1])
    xyz = p.to_global([[1, 1], [4, 5]])
    np.testing.assert_allclose(p.to_local(xyz), [[1, 1], [4, 5]])
    f = p.calc_factors(xyz, conwidth=0, a=2)
    g = other.calc_factors(xyz, conwidth=0, a=2)
    np.testing.assert_allclose(f.influence, g.influence)
    np.testing.assert_allclose(f.weights, g.weights)


def test_extreme_decay_and_invalid_mask_and_shape(tmp_path):
    p = make_plane()
    f = p.calc_factors([[1, 1, 1.0e6]], conwidth=0, a=np.finfo(float).tiny)
    assert f.influence[0] == 0
    with pytest.raises(ValueError):
        p.calc_factors([[0, 0, 0]], conwidth=0, a=1, active=np.nan)
    with pytest.raises(ValueError):
        p.calc_factors([[1.0e200, 0, 0]], conwidth=0, a=1)
    with pytest.raises(ValueError):
        f.apply(np.ones((4, 3)), 1)  # Source array must be (ny, nx).
    with pytest.raises(ValueError):
        f.weights[0, 0] = 5
    np.savez(tmp_path / "bad.npz", weights=np.ones((2, 3)))
    with pytest.raises(ValueError):
        GPlaneFactors.load(tmp_path / "bad.npz")


@pytest.mark.parametrize("kind", ["horizontal", "dipping", "vertical", "xz"])
@pytest.mark.parametrize("bounded", [False, True])
@pytest.mark.parametrize(
    "action,interp_transform,transform,power",
    [
        ("blend", "none", "none", 2.0),
        ("add", "log", "log", 1.0),
        ("blend_never_higher", "log", "none", 2.5),
        ("blend_never_lower", "none", "log", 2.0),
    ],
)
def test_against_plproc(
    tmp_path, kind, bounded, action, interp_transform, transform, power
):
    exe = os.environ.get("PLPROC_EXE")
    if not exe:
        pytest.skip("set PLPROC_EXE to run the real-executable parity tests")
    p = make_plane(kind, bounded)
    rng = np.random.default_rng(1979)
    uv = rng.uniform([-2, -2], [6, 5], (31, 2))
    targets = p.to_global(uv) + rng.uniform(-4, 4, (31, 1)) * p.normal
    values = rng.uniform(0.2, 3.0, (3, 4))
    background = rng.uniform(0.2, 5.0, 31)
    active = np.ones(31, int)
    active[::7] = 0
    settings = dict(
        conwidth=0.5,
        a=1.5,
        inverse_power=power,
        anis_ratio=2.0,
        anis_bearing=30.0,
        active=active,
    )
    f = p.calc_factors(targets, **settings)
    actual = f.apply(
        values,
        background,
        action=action,
        interp_transform=interp_transform,
        transform=transform,
    )
    prepare_plproc(
        tmp_path,
        p,
        targets,
        values,
        background,
        **settings,
        action=action,
        interp_transform=interp_transform,
        transform=transform,
    )
    expected = run_plproc(exe, tmp_path)
    # PLPROC's column writer emits seven significant figures by default.
    np.testing.assert_allclose(actual, expected, rtol=1.0e-6, atol=5.0e-7)


@pytest.mark.parametrize(
    "moved,left_handed,ratio",
    [
        (False, False, 1),
        (True, False, 1),
        (False, True, 0.5),
        (True, True, 0.5),
    ],
)
def test_moved_and_left_handed_against_plproc(tmp_path, moved, left_handed, ratio):
    exe = os.environ.get("PLPROC_EXE")
    if not exe:
        pytest.skip("set PLPROC_EXE to run the real-executable parity tests")
    p = make_plane("dipping", True)
    xyz = p.to_global([[0.5, 0.5], [1.5, 0.5], [2, 1], [-2, 1], [3, 4]])
    xyz += np.arange(5)[:, None] * p.normal
    values = np.arange(1, 13).reshape(3, 4)
    if moved:
        # Move the plane relative to fixed targets, then recompute its factors.
        points = p.points + [0, 0, 1.25]
        p = GPlane(points, p.grid, p.shape, p.spacing, boundary=p.boundary)
    if left_handed:
        p = GPlane(
            p.points, [[0, 0], [4, 0], [0, -3]], p.shape, p.spacing, boundary=p.boundary
        )
    settings = dict(conwidth=0.2, a=1.5, anis_ratio=ratio)
    actual = p.calc_factors(xyz, **settings).apply(values, 5)
    prepare_plproc(tmp_path, p, xyz, values, 5, **settings)
    expected = run_plproc(exe, tmp_path)
    np.testing.assert_allclose(actual, expected, rtol=1.0e-6, atol=5.0e-7)
