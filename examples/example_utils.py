"""Fetch and cache the verified executables used by the worked examples."""

import argparse
import hashlib
import json
import os
import platform
import shutil
import tarfile
import tempfile
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path


class _PublicRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        redirected = super().redirect_request(req, fp, code, msg, headers, newurl)
        if redirected is not None:
            redirected.remove_header("Authorization")
        return redirected


def _sha256(path):
    with Path(path).open("rb") as stream:
        digest = hashlib.sha256()
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _platform():
    system, machine = platform.system(), platform.machine().lower()
    machine = {"amd64": "x86_64", "aarch64": "arm64"}.get(machine, machine)
    return f"{system.lower()}-{machine}"


def _download(spec, destination):
    """Try the public asset and API routes, validating publisher checksums."""
    repo, version, asset = spec["repository"], spec["version"], spec["asset"]
    browser = f"https://github.com/{repo}/releases/download/{version}/{asset}"
    api = f"https://api.github.com/repos/{repo}/releases/assets/{spec['asset_id']}"
    opener = urllib.request.build_opener(_PublicRedirect())
    last_error = None
    for attempt, url in enumerate((browser, api, browser + "?download=1")):
        headers = {"User-Agent": "pypestutils-examples"}
        if url == api:
            headers["Accept"] = "application/octet-stream"
            if token := os.environ.get("GITHUB_TOKEN"):
                headers["Authorization"] = f"Bearer {token}"
        try:
            request = urllib.request.Request(url, headers=headers)
            with opener.open(request, timeout=45) as response:
                with Path(destination).open("wb") as stream:
                    shutil.copyfileobj(response, stream)
            if _sha256(destination) != spec["sha256"]:
                raise ValueError(f"checksum mismatch for {asset}")
            return
        except (OSError, ValueError) as error:
            last_error = error
            if attempt < 2:
                print(f"Retrying {asset} using another download route", flush=True)
                time.sleep(attempt + 1)
    raise RuntimeError(
        f"Could not download verified executable archive {asset}"
    ) from last_error


def _extract(archive, executable, destination):
    """Extract only the named regular file, without using archive paths."""
    if zipfile.is_zipfile(archive):
        with zipfile.ZipFile(archive) as bundle:
            matches = [n for n in bundle.namelist() if Path(n).name == executable]
            if len(matches) != 1:
                raise ValueError(f"expected exactly one {executable} in archive")
            with (
                bundle.open(matches[0]) as source,
                Path(destination).open("wb") as target,
            ):
                shutil.copyfileobj(source, target)
    else:
        with tarfile.open(archive) as bundle:
            matches = [
                m
                for m in bundle.getmembers()
                if m.isfile() and Path(m.name).name == executable
            ]
            if len(matches) != 1:
                raise ValueError(f"expected exactly one {executable} in archive")
            with (
                bundle.extractfile(matches[0]) as source,
                Path(destination).open("wb") as target,
            ):
                shutil.copyfileobj(source, target)
    Path(destination).chmod(0o755)


def prepare_executables(destination=None, *, cache_dir=None):
    """Return pinned MF6/PEST++ binaries, optionally copying them to a model bin."""
    manifest_path = Path(__file__).with_name("example_executables.json")
    manifest = json.loads(manifest_path.read_text())
    platform_key = _platform()
    if platform_key not in manifest:
        raise ValueError(f"example executables are unavailable for {platform_key}")
    specs = manifest[platform_key]
    key = hashlib.sha256(json.dumps(specs, sort_keys=True).encode()).hexdigest()[:16]
    if cache_dir is None:
        cache_dir = os.environ.get(
            "PYPESTUTILS_EXAMPLE_CACHE",
            Path(__file__).parent.parent / ".cache" / "example-executables",
        )
    cache = Path(cache_dir).resolve() / platform_key / key
    cache.mkdir(parents=True, exist_ok=True)
    receipt_path = cache / "checksums.json"
    try:
        receipt = json.loads(receipt_path.read_text())
    except (OSError, ValueError):
        receipt = {}
    if not isinstance(receipt, dict):
        receipt = {}
    for program, spec in specs.items():
        executable = program + (".exe" if platform_key.startswith("windows-") else "")
        target = cache / executable
        valid = target.is_file() and receipt.get(executable) == {
            "archive": spec["sha256"],
            "binary": _sha256(target),
        }
        if not valid:
            print(
                f"Preparing {program} {spec['version']} for {platform_key}", flush=True
            )
            with tempfile.TemporaryDirectory(prefix="ppu-example-bin-") as temp:
                archive = Path(temp) / spec["asset"]
                binary = Path(temp) / executable
                _download(spec, archive)
                _extract(archive, executable, binary)
                shutil.copy2(binary, target)
            receipt[executable] = {"archive": spec["sha256"], "binary": _sha256(target)}
            receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
        target.chmod(0o755)
    if destination is None:
        return cache
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    for program in specs:
        executable = program + (".exe" if platform_key.startswith("windows-") else "")
        if destination != cache:
            shutil.copy2(cache / executable, destination / executable)
    return destination


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path)
    options = parser.parse_args()
    print(prepare_executables(cache_dir=options.cache_dir))
