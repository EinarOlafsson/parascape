"""The ParaScape relay client (instructions 94, 97): what a report holds, and that nothing personal
gets into it. No network: the relay is never called here."""
import sys

import pytest
from PySide6 import QtCore

from parascape import relay as R


@pytest.fixture(autouse=True)
def settings(tmp_path, monkeypatch):
    from parascape.backdrops import _compat as AC
    path = str(tmp_path / "relay.ini")
    monkeypatch.setattr(AC, "_settings",
                        lambda: QtCore.QSettings(path, QtCore.QSettings.Format.IniFormat))


def _raise(message):
    try:
        raise ValueError(message)
    except ValueError:
        return sys.exc_info()


def test_paths_emails_tokens_and_the_user_name_are_removed(monkeypatch):
    monkeypatch.setenv("USER", "alice")
    text = R.scrub("cannot read /home/alice/data/hits.csv for alice@lab.org from 10.1.2.3 "
                   "token abcdefghijklmnopqrstuvwxyz0123456789")
    assert "/home" not in text and "alice" not in text
    assert "hits.csv" in text, "the file name is kept, the path is not"
    assert "<email>" in text and "<ip>" in text and "<token>" in text


def test_a_report_holds_file_names_not_paths_and_names_the_parascape_module():
    report = R.report_from_exception(*_raise("bad /home/x/y/z.csv"))
    assert report["exception_type"] == "ValueError"
    assert all("/" not in f["file"] for f in report["frames"])
    assert report["frames"][-1]["file"] == "test_relay.py"
    assert "/home" not in report["message"]
    assert report["schema_version"] == 1 and report["os"] in ("linux", "windows", "macos")


def test_the_install_id_is_made_once_and_kept():
    first = R.install_id()
    assert first == R.install_id()


def test_consent_defaults_to_asking():
    assert R.consent() == "ask"
    R.set_consent("always")
    assert R.consent() == "always"
    R.set_consent("nonsense")
    assert R.consent() == "ask"


def test_a_hand_written_report_is_a_user_report_with_no_traceback():
    report = R.problem_report("The map went blank after I picked a colour")
    assert report["kind"] == "user_error" and report["frames"] == []
    assert report["exception_type"] == "ProblemReport"


def test_the_preview_shows_what_is_sent_but_not_the_install_id():
    report = R.report_from_exception(*_raise("x"))
    text = R.preview(report)
    assert "ValueError" in text and report["install_id"] not in text


def test_the_same_error_has_the_same_signature():
    a = R.report_from_exception(*_raise("one"))
    b = R.report_from_exception(*_raise("two"))
    assert R.signature(a) == R.signature(b)


def test_only_https_or_this_computer_can_override_the_relay(monkeypatch):
    monkeypatch.setenv(R.RELAY_ENV, "http://evil.example.com")
    assert R.relay_url() == R.RELAY_URL
    monkeypatch.setenv(R.RELAY_ENV, "http://localhost:8787")
    assert R.relay_url() == "http://localhost:8787"
