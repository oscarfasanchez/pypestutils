# GPLANE, MODFLOW-USG and DISU fork release

Tag **`v0.3.2-usg.1`**, Python package version **`0.3.2+usg.1`**.
The distribution and import name remain `pypestutils`. This independent fork
release is hosted on GitHub.

This release adds native 3D GPLANE structural overlays: arbitrary XYZ targets,
polygon limits, anisotropic PLPROC-compatible IDW, Gaussian influence, four
actions, independent log transforms, masks, reusable factors and PLPROC
geometry-specification exchange. It also includes the previous MODFLOW-USG,
MODFLOW 6 DISU, layer-aware helper and interpolation-performance changes.

[Changes](https://github.com/oscarfasanchez/pypestutils/blob/v0.3.2-usg.1/FORK_CHANGES.md)
· [GPLANE guide](https://github.com/oscarfasanchez/pypestutils/blob/v0.3.2-usg.1/docs/gplane.md)
· [PLPROC benchmarks](https://github.com/oscarfasanchez/pypestutils/blob/v0.3.2-usg.1/benchmarks/README.md)

## Install with pip

Run this on a supported platform; pip downloads the matching compiled wheel:

```sh
python -m pip install -r https://github.com/oscarfasanchez/pypestutils/releases/download/v0.3.2-usg.1/requirements-usg.txt
```

No repository checkout, Fortran compiler, or manual wheel download is required.
This also upgrades an earlier fork or upstream installation in that environment.
An ordinary `pip install pypestutils` continues to select upstream from PyPI.

For a direct Windows x64 wheel installation:

```sh
python -m pip install "pypestutils @ https://github.com/oscarfasanchez/pypestutils/releases/download/v0.3.2-usg.1/pypestutils-0.3.2+usg.1-py3-none-win_amd64.whl"
```

## Install with Conda / environment.yml

Download `environment-usg.yml` from the release assets, then run:

```sh
conda env create -f environment-usg.yml
conda activate ppu-usg
```

Or add this to your own `environment.yml`:

```yaml
dependencies:
  - python=3.12
  - pip
  - pip:
      - "-r https://github.com/oscarfasanchez/pypestutils/releases/download/v0.3.2-usg.1/requirements-usg.txt"
```

This installs the GitHub wheel through pip inside the Conda environment. Remove
an existing `pypestutils` Conda/PyPI dependency from that file when using the
fork entry. Upstream and fork installations share the import name; use separate
environments to keep both.

## Platforms and verification

Python 3.10 or newer is required. The automatic requirements and environment
file support Windows x64, Linux glibc 2.27+ on x64 or ARM64, and macOS 15+ on
Intel or Apple Silicon. Alpine/musl Linux wheels are separate release assets;
install the matching musllinux wheel URL directly. Windows ARM64 and 32-bit
systems are not covered.

After installation:

```sh
python -c "from pypestutils.version import __version__; from pypestutils.pestutilslib import PestUtilsLib; print(__version__); print(PestUtilsLib().pestutils.calc_gplane_factors)"
```

The GPLANE implementation passed 155 local tests (eight existing skips),
including 36 real PLPROC 5.0 comparisons. The feature commit also passed all
six Linux/macOS/Windows CI jobs with Python 3.10 and 3.14. The wheel pipeline
checks that Python imports the installed wheel and loads its bundled native
library, then runs the non-notebook regression suite on every wheel platform,
including Linux ARM64 and musllinux. PLPROC-specific tests skip in CI because
the executable is not redistributed. Optional-dependency tests can also skip.
Wheel validation uses Python 3.10; the source CI covers Python 3.10 and 3.14.

`SHA256SUMS` records the checksums of all other release assets.
