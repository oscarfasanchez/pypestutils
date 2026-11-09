"""Regression coverage for example downloads, archive handling and caching."""

import hashlib
import importlib.util
import io
import json
import shutil
import tarfile
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

import pytest

module_path = Path(__file__).parent.parent / "examples" / "example_utils.py"
spec = importlib.util.spec_from_file_location("ppu_example_utils", module_path)
example_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(example_utils)


def test_download_fallback_validates_checksum_and_limits_auth(tmp_path, monkeypatch):
    payload = b"verified archive content"
    calls = []

    class Opener:
        def open(self, request, timeout):
            calls.append(request)
            if len(calls) == 1:
                raise urllib.error.URLError("temporary outage")
            return io.BytesIO(payload)

    monkeypatch.setattr(
        example_utils.urllib.request, "build_opener", lambda _: Opener()
    )
    monkeypatch.setattr(example_utils.time, "sleep", lambda _: None)
    monkeypatch.setenv("GITHUB_TOKEN", "test-token")
    spec = dict(
        repository="publisher/tool",
        version="1",
        asset="tool.zip",
        asset_id=42,
        sha256=hashlib.sha256(payload).hexdigest(),
    )
    target = tmp_path / "archive.zip"
    example_utils._download(spec, target)
    assert target.read_bytes() == payload
    assert "Authorization" not in calls[0].headers
    assert calls[1].get_header("Authorization") == "Bearer test-token"
    redirected = example_utils._PublicRedirect().redirect_request(
        calls[1],
        None,
        302,
        "Found",
        {},
        "https://release-assets.githubusercontent.com/file",
    )
    assert redirected.get_header("Authorization") is None


def test_bad_download_never_passes_checksum(tmp_path, monkeypatch):
    class Opener:
        def open(self, request, timeout):
            return io.BytesIO(b"incorrect archive")

    monkeypatch.setattr(
        example_utils.urllib.request, "build_opener", lambda _: Opener()
    )
    monkeypatch.setattr(example_utils.time, "sleep", lambda _: None)
    spec = dict(
        repository="publisher/tool",
        version="1",
        asset="tool.zip",
        asset_id=42,
        sha256="0" * 64,
    )
    with pytest.raises(RuntimeError, match="verified executable archive"):
        example_utils._download(spec, tmp_path / "bad.zip")


def test_binary_cache_reuse_and_corruption_recovery(tmp_path, monkeypatch):
    archive = tmp_path / "bundle.zip"
    with zipfile.ZipFile(archive, "w") as bundle:
        bundle.writestr("bin/mf6.exe", b"mf6 executable")
        bundle.writestr("bin/pestpp-ies.exe", b"pestpp executable")
    entry = dict(
        repository="publisher/tool",
        version="1",
        asset="bundle.zip",
        asset_id=42,
        sha256=example_utils._sha256(archive),
    )
    manifest = {"windows-x86_64": {"mf6": entry, "pestpp-ies": entry}}
    (tmp_path / "example_executables.json").write_text(json.dumps(manifest))
    monkeypatch.setattr(example_utils, "__file__", str(tmp_path / "example_utils.py"))
    monkeypatch.setattr(example_utils, "_platform", lambda: "windows-x86_64")
    calls = []

    def download(spec, target):
        calls.append(spec)
        shutil.copyfile(archive, target)

    monkeypatch.setattr(example_utils, "_download", download)
    cache = example_utils.prepare_executables(
        tmp_path / "first", cache_dir=tmp_path / "cache"
    )
    assert (cache / "mf6.exe").read_bytes() == b"mf6 executable"
    assert len(calls) == 2
    stored = example_utils.prepare_executables(cache_dir=tmp_path / "cache")
    (stored / "mf6.exe").write_bytes(b"corrupt")
    example_utils.prepare_executables(tmp_path / "second", cache_dir=tmp_path / "cache")
    assert len(calls) == 3
    assert (tmp_path / "second" / "mf6.exe").read_bytes() == b"mf6 executable"
    example_utils.prepare_executables(cache_dir=tmp_path / "cache")
    assert len(calls) == 3


def test_extract_uses_only_requested_file_and_rejects_ambiguity(tmp_path):
    archive = tmp_path / "bundle.tar.gz"
    with tarfile.open(archive, "w:gz") as bundle:
        member = tarfile.TarInfo("../../mf6")
        member.size = 3
        bundle.addfile(member, io.BytesIO(b"exe"))
    target = tmp_path / "safe" / "mf6"
    target.parent.mkdir()
    example_utils._extract(archive, "mf6", target)
    assert target.read_bytes() == b"exe"
    assert not (tmp_path / "mf6").exists()
    ambiguous = tmp_path / "ambiguous.zip"
    with zipfile.ZipFile(ambiguous, "w") as bundle:
        bundle.writestr("one/mf6", b"one")
        bundle.writestr("two/mf6", b"two")
    with pytest.raises(ValueError, match="exactly one"):
        example_utils._extract(ambiguous, "mf6", target)
