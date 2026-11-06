# MODFLOW-USG and DISU fork release

This release contains the changes in `feat/mfusg-grids`, including MODFLOW 6
DISU support, MODFLOW-USG grid and dependent-variable readers, layer-aware 2D
helpers, and the interpolation ring-search speedup. See
[FORK_CHANGES.md](https://github.com/oscarfasanchez/pypestutils/blob/v0.3.1-usg.1/FORK_CHANGES.md)
for the scope and upstream baseline.

The Git tag is `v0.3.1-usg.1`; the Python distribution version is
`0.3.1+usg.1`. The distribution and import name remain `pypestutils`.
This is an independent fork release hosted on GitHub.

## Install with Conda

Download `environment-usg.yml` from the release assets, or copy its `pip:`
entries into your environment file, then run:

```sh
conda env create -f environment-usg.yml
conda activate ppu-usg
```

Pip automatically downloads the matching compiled wheel. No local Fortran
compiler is required. The environment file supports Windows x86_64, Linux
glibc 2.27+ on x86_64 or ARM64, and macOS 15+ on Intel or Apple Silicon.
Alpine/musl Linux wheels are also provided as separate release assets.

## Install with pip

For Windows x86_64:

```sh
python -m pip install "pypestutils @ https://github.com/oscarfasanchez/pypestutils/releases/download/v0.3.1-usg.1/pypestutils-0.3.1%2Busg.1-py3-none-win_amd64.whl"
```

For other supported platforms, use the matching wheel URL from
`environment-usg.yml`, or download the appropriate wheel and run
`python -m pip install <wheel-file>`. Python 3.10 or newer is required.

Use a separate environment to retain an upstream installation: the fork and
upstream both provide `import pypestutils`.

## Build and validation

The tag triggers native-library compilation and wheel builds in GitHub
Actions. Windows, Linux x86_64, and both macOS wheels use the existing
wheel test suite. Linux ARM64 and musllinux wheel tests are skipped by the existing
build configuration. `SHA256SUMS` records the release asset checksums.
