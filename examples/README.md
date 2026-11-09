# Running the worked notebooks

Use an environment containing the fork's compiled library, FloPy, pyEMU,
Jupyter and Matplotlib. From this directory, open the notebooks with
`python -m jupyter lab`. The two PstFrom notebooks run real MODFLOW 6 and PEST++
executables.

`example_utils.prepare_executables()` prepares **MF6 6.7.0** and
**PEST++ 5.2.25**, copying only `mf6` and `pestpp-ies` into each notebook's bin
directory. Versions, platform assets, and publisher SHA-256 checksums are
pinned in [example_executables.json](example_executables.json). Mac Intel and
Apple Silicon use their respective native binaries.

Archives come from the official
[MODFLOW 6 release](https://github.com/MODFLOW-ORG/modflow6/releases/tag/6.7.0)
and [PEST++ release](https://github.com/pestpp/pestpp/releases/tag/5.2.25).
Download failures fall back to GitHub's asset API and a fresh public download
route. Every archive must match its pinned checksum before extraction. Only
the requested executable is extracted. GitHub authentication is used only for
API requests and is removed on redirects.

Verified binaries are cached under `.cache/example-executables` at the
repository root. Set `PYPESTUTILS_EXAMPLE_CACHE` to use another directory.
The cache is keyed by platform and manifest contents; damaged cached binaries
are fetched again. Preload it from the repository root with:

```sh
python examples/example_utils.py
```

The notebook tests copy their input models into temporary directories and
execute there, preserving tracked notebooks and model inputs. Run them from
the repository root with:

```sh
python -m pytest tests/test_notebooks.py -v
```

CI prepares the binaries once per job and caches them across runs. Notebook
execution still runs the complete model and preprocessing workflow; download
or model failures fail the checks. The source and wheel test runners propagate
test failures as nonzero exit codes.
