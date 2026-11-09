"""Execute worked notebooks without modifying repository inputs."""

import os
import shutil
import sys
from pathlib import Path
from subprocess import run

import pytest

nbformat = pytest.importorskip("nbformat")
examples_dir = Path(__file__).parent.parent / "examples"


@pytest.fixture(scope="session")
def example_binary_cache(tmp_path_factory):
    return os.environ.get(
        "PYPESTUTILS_EXAMPLE_CACHE", str(tmp_path_factory.mktemp("example-executables"))
    )


@pytest.mark.parametrize(
    "nb_file",
    [
        pth.name
        for pth in sorted(examples_dir.glob("*.ipynb"))
        if "understanding_variograms" not in pth.name
    ],
)
def test_notebooks(nb_file, tmp_path, example_binary_cache):
    nb = nbformat.read(examples_dir / nb_file, as_version=4)
    for cell in nb["cells"]:
        for line in cell["source"].splitlines():
            if line.startswith("import ") and "ppu_helpers" not in line:
                pytest.importorskip(line.split()[1])

    # Copy authoritative inputs, never previous generated workspaces or binaries.
    for name in ("freyberg_aux_files", "freyberg_monthly", "freyberg_quadtree"):
        shutil.copytree(examples_dir / name, tmp_path / name)
    for pattern in ("*.py", "*.json"):
        for path in examples_dir.glob(pattern):
            shutil.copy2(path, tmp_path / path.name)
    nbformat.write(nb, tmp_path / nb_file)
    env = dict(
        os.environ,
        PYPESTUTILS_EXAMPLE_CACHE=str(example_binary_cache),
        MPLBACKEND="Agg",
    )
    run(
        [
            sys.executable,
            "-m",
            "jupyter",
            "nbconvert",
            "--to",
            "notebook",
            "--execute",
            "--ExecutePreprocessor.timeout=300",
            "--inplace",
            nb_file,
        ],
        cwd=tmp_path,
        env=env,
        check=True,
    )
