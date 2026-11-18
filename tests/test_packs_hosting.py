"""Checks for the pack hosting added for instruction 99 step 1."""
import json, os, pytest
from parascape import packs


def test_default_host_is_hugging_face():
    assert packs.packs_host().startswith("https://huggingface.co/")


def test_host_setting_is_honoured(monkeypatch):
    monkeypatch.setenv(packs.ENV_PACKS_HOST, "https://huggingface.co/datasets/x/y/resolve/main/p/")
    assert packs.packs_host() == "https://huggingface.co/datasets/x/y/resolve/main/p"


@pytest.mark.parametrize("url", ["http://huggingface.co/a.zip", "https://example.com/a.zip",
                                 "ftp://huggingface.co/a.zip"])
def test_other_hosts_and_plain_http_are_refused(url):
    with pytest.raises(ValueError):
        packs._checked_pack_url(url)


def test_a_known_host_is_allowed():
    packs._checked_pack_url("https://huggingface.co/datasets/a/b/resolve/main/packs/x.zip")


def test_index_is_empty_or_format_one():
    assert isinstance(packs.pack_index(), dict)


def test_unknown_organism_is_refused_by_name(monkeypatch):
    monkeypatch.setattr(packs, "pack_index", lambda: {})
    with pytest.raises(LookupError) as e:
        packs.pack_release("Tg")
    assert "Tg" in str(e.value)


def test_newest_version_is_chosen_and_the_checksum_comes_from_the_index(monkeypatch):
    monkeypatch.setattr(packs, "pack_index", lambda: {
        "Tg": {"0.1.0": {"sha256": "a" * 64, "size": 1, "path": "Tg/0.1.0/Tg-0.1.0.zip"},
               "0.2.0": {"sha256": "b" * 64, "size": 2, "path": "Tg/0.2.0/Tg-0.2.0.zip"}}})
    url, sha, version = packs.pack_release("Tg")
    assert version == "0.2.0" and sha == "b" * 64 and url.endswith("Tg/0.2.0/Tg-0.2.0.zip")
    assert url.startswith("https://huggingface.co/")


def test_download_pack_never_takes_a_checksum_from_the_caller(monkeypatch):
    monkeypatch.setattr(packs, "pack_index", lambda: {
        "Tg": {"0.1.0": {"sha256": "c" * 64, "size": 1, "path": "Tg/0.1.0/Tg-0.1.0.zip"}}})
    seen = {}
    monkeypatch.setattr(packs, "download", lambda url, sha, organism=None: seen.update(
        url=url, sha=sha, organism=organism) or "installed")
    assert packs.download_pack("Tg") == "installed"
    assert seen["sha"] == "c" * 64 and seen["organism"] == "Tg"


def test_url_override_still_must_be_a_known_host(monkeypatch):
    monkeypatch.setattr(packs, "pack_index", lambda: {
        "Tg": {"0.1.0": {"sha256": "d" * 64, "size": 1, "path": "Tg/0.1.0/Tg-0.1.0.zip"}}})
    monkeypatch.setenv(packs.ENV_PACKS_URL, "https://example.com/evil.zip")
    with pytest.raises(ValueError):
        packs.download_pack("Tg")


def test_the_index_ships_with_the_wheel():
    """The checksum must travel inside the application, or `download_pack` has nothing to trust."""
    import tomllib
    from pathlib import Path
    root = Path(__file__).resolve().parent.parent
    cfg = tomllib.loads((root / "pyproject.toml").read_text())
    listed = cfg["tool"]["setuptools"]["package-data"]["parascape"]
    assert packs.PACK_INDEX in listed, f"{packs.PACK_INDEX} must be in package-data"
    assert (root / "parascape" / packs.PACK_INDEX).exists()
