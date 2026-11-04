"""MODFLOW-USG unstructured grids: reading, install, interpolation and mod2obs.

vorflow_usg is a layered model built with vorflow and run with MODFLOW-USG 1.5
(mfusg, single precision) and mfusgdbl (double precision): three layers of 55
cells, the second split into two sublayers inside a pit, so four layers of 55,
55, 8 and 55 nodes. Two output times. vorflow_usg_points.csv holds heads at eight
points, from a weighted sum of the heads read by flopy's HeadUFile.
"""

import shutil
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pypestutils.helpers import (
    MFUSG_GRID_ARGS,
    get_2d_grid_info_from_file,
    get_grid_info_from_mfusg_files,
    mod2obs_mfusg,
    read_mfusg_disu,
)
from pypestutils.pestutilslib import PestUtilsLib, PestUtilsLibError

data_dir = Path(__file__).parent / "data"
DISU = data_dir / "vorflow_usg.disu"
GSF = data_dir / "vorflow_usg.gsf"
NODELAY = [55, 55, 8, 55]


@pytest.fixture
def lib_grid():
    lib = PestUtilsLib()
    grid = get_grid_info_from_mfusg_files(DISU, GSF)
    lib.install_mfusg_grid("usg", **{k: grid[k] for k in MFUSG_GRID_ARGS})
    yield lib, grid
    lib.uninstall_mf6_grid("usg")


def test_read_mfusg_files():
    grid = get_grid_info_from_mfusg_files(DISU, GSF)
    assert grid["ncells"] == sum(NODELAY)
    assert grid["nodelay"].tolist() == NODELAY
    np.testing.assert_array_equal(grid["layer"], np.repeat(np.arange(1, 5), NODELAY))
    assert grid["ia"][0] == 1 and grid["ia"][-1] == len(grid["ja"]) + 1
    # each node is listed first among its own connections
    np.testing.assert_array_equal(grid["ja"][grid["ia"][:-1] - 1], np.arange(1, 174))
    assert (grid["top"] > grid["bot"]).all()
    assert (grid["area"] > 0).all()
    assert (np.diff(grid["iavert"]) >= 3).all()


def test_read_mfusg_disu_array_controls(tmp_path):
    (tmp_path / "iac.txt").write_text("2 2\n2 2\n")
    disu = tmp_path / "t.disu"
    disu.write_text(
        "# two layers of two nodes\n"
        "4 2 8 0 1 4 2 0\n"
        "0 0\n"
        "CONSTANT 2\n"
        "CONSTANT 10.0\nINTERNAL 1.0 (FREE) -1\n5.0 4.0\n"
        "CONSTANT 5.0\nCONSTANT 0.0\n"
        "CONSTANT 1.0\nCONSTANT 1.0\n"
        "OPEN/CLOSE iac.txt 1 (FREE) -1\n"
        "INTERNAL 1 (FREE) -1\n1 2 2 1 3 4 4 3\n"
    )
    res = read_mfusg_disu(str(disu))
    assert res["nodelay"].tolist() == [2, 2]
    np.testing.assert_allclose(res["top"], [10.0, 10.0, 5.0, 4.0])
    np.testing.assert_allclose(res["bot"], [5.0, 5.0, 0.0, 0.0])
    assert res["iac"].tolist() == [2, 2, 2, 2]
    assert res["ja"].tolist() == [1, 2, 2, 1, 3, 4, 4, 3]


def test_install_mfusg_grid_cell_centres(lib_grid):
    lib, grid = lib_grid
    x, y, z = lib.get_cell_centres_mf6("usg", grid["ncells"])
    np.testing.assert_allclose(x, grid["x"])
    np.testing.assert_allclose(y, grid["y"])
    np.testing.assert_allclose(z, 0.5 * (grid["top"] + grid["bot"]))


def test_install_mfusg_grid_checks_its_input():
    lib = PestUtilsLib()
    grid = get_grid_info_from_mfusg_files(DISU, GSF)
    args = {k: grid[k] for k in MFUSG_GRID_ARGS}
    bad = dict(args, ia=args["ia"] + 1)
    with pytest.raises(PestUtilsLibError, match="IA array"):
        lib.install_mfusg_grid("usg_bad", **bad)
    bad = dict(args, layer=np.zeros_like(args["layer"]))
    with pytest.raises(PestUtilsLibError, match="LAYER"):
        lib.install_mfusg_grid("usg_bad", **bad)
    # nothing was left installed
    lib.install_mfusg_grid("usg_bad", **args)
    lib.uninstall_mf6_grid("usg_bad")


@pytest.mark.parametrize("precision, iprec", [("single", 1), ("double", 2)])
def test_interp_from_mfusg_depvar_file(tmp_path, lib_grid, precision, iprec):
    lib, _ = lib_grid
    pts = pd.read_csv(data_dir / "vorflow_usg_points.csv")
    factorfile = tmp_path / "usg.fac"
    success = lib.calc_mf6_interp_factors(
        "usg", pts.x, pts.y, pts.layer, factorfile, "text", tmp_path / "usg.bln"
    )
    assert success.all()
    hds = data_dir / f"vorflow_usg_{precision}.hds"
    info = lib.inquire_modflow_binary_file_specs(hds, "", 22, 1)
    assert info["iprec"] == iprec and info["ntime"] == 2
    res = lib.interp_from_mfusg_depvar_file(
        hds, factorfile, "text", 2, "head", iprec, 1e20, 0, -999.0, len(pts)
    )
    assert res["nproctime"] == 2
    np.testing.assert_allclose(res["simtime"], [10.0, 20.0])
    expected = pts[[f"{precision}_t1", f"{precision}_t2"]].to_numpy().T
    np.testing.assert_allclose(res["simstate"], expected, atol=1e-8)


def test_interp_from_mfusg_depvar_file_wrong_precision(tmp_path, lib_grid):
    lib, _ = lib_grid
    pts = pd.read_csv(data_dir / "vorflow_usg_points.csv")
    factorfile = tmp_path / "usg.fac"
    lib.calc_mf6_interp_factors(
        "usg", pts.x, pts.y, pts.layer, factorfile, "text", tmp_path / "usg.bln"
    )
    with pytest.raises(PestUtilsLibError):
        lib.interp_from_mfusg_depvar_file(
            data_dir / "vorflow_usg_double.hds",
            factorfile,
            "text",
            2,
            "head",
            1,
            1e20,
            0,
            -999.0,
            len(pts),
        )


def test_mod2obs_mfusg(tmp_path):
    for f in [DISU, GSF, data_dir / "vorflow_usg_double.hds"]:
        shutil.copy(f, tmp_path)
    pts = pd.read_csv(data_dir / "vorflow_usg_points.csv")
    obs = pts[["name", "x", "y", "layer"]].rename(columns={"name": "site"})
    obs["datetime"] = datetime(2000, 1, 21)  # day 20, the last output time
    res = mod2obs_mfusg(
        str(tmp_path / DISU.name),
        str(tmp_path / GSF.name),
        str(tmp_path / "vorflow_usg_double.hds"),
        obs,
        pd.Timestamp(datetime(2000, 1, 1)),
    )
    sim = res["interpolated_results"].set_index("site").simulated
    np.testing.assert_allclose(sim[pts.name], pts.double_t2, atol=1e-8)
    assert res["all_results"].shape == (2, len(pts))


def test_2d_grid_info_of_a_mfusg_layer():
    grid = get_grid_info_from_mfusg_files(DISU, GSF)
    info = get_2d_grid_info_from_file(str(GSF), layer=3)
    assert info["nnodes"] == 8
    np.testing.assert_array_equal(info["nodes"], np.flatnonzero(grid["layer"] == 3))
    np.testing.assert_allclose(info["x"], grid["x"][grid["layer"] == 3])
    np.testing.assert_allclose(info["area"], grid["area"][grid["layer"] == 3])
