"""PLPROC executable oracle. Never substitute a Python result for this oracle."""

import subprocess
from pathlib import Path

import numpy as np


def prepare_plproc(
    work,
    plane,
    targets,
    values,
    background,
    *,
    conwidth,
    a,
    inverse_power=2.0,
    anis_ratio=1.0,
    anis_bearing=0.0,
    action="blend",
    interp_transform="none",
    transform="none",
    active=1,
):
    """Write independent PLPROC inputs for the same geometry and property data."""
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    plane.write_spec(work / "plane.spec")
    n = len(targets)
    mask = np.broadcast_to(active, n) != 0
    table = np.column_stack(
        (np.arange(1, n + 1), targets, np.broadcast_to(background, n), mask.astype(int))
    )
    np.savetxt(work / "targets.dat", table, fmt="%.17g")
    values = np.asarray(values)
    if values.size == 1:
        values = np.full(len(plane.source_uv), values.item())
    else:
        values = values.reshape(plane.shape[1], plane.shape[0]).ravel()
    np.savetxt(work / "values.dat", values, fmt="%.17g")
    (work / "run.plproc").write_text(
        "t=read_list_file(file='targets.dat',dimensions=3,id_type='indexed',"
        "plist=b;column=5,slist=mask;column=6)\n"
        "g=new_gplane_clist(gplane='gp',gplane_file='plane.spec',"
        f"interp_transform='{interp_transform}',grid_interp='ivd';inverse_power={inverse_power:.17g},"
        f"interp_anis_ratio={anis_ratio:.17g},interp_anis_bearing={anis_bearing:.17g})\n"
        "v=new_plist(reference_clist=g,value=1)\n"
        "read_column_data_file(file='values.dat',plist=v;column=1)\n"
        f"b(select=(mask==1))=v.inform_from_gplane(action='{action}';"
        f"transform='{transform}',conwidth={conwidth:.17g},a={a:.17g})\n"
        "write_column_data_file(file='result.dat';delim='space',plist=b,header='no')\n",
        encoding="ascii",
    )


def run_plproc(exe, work):
    """Run PLPROC and require its success marker, not merely exit code zero."""
    work = Path(work)
    output = work / "result.dat"
    output.unlink(missing_ok=True)
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
    return np.loadtxt(output, ndmin=1)
