"""Talking to the ParaScape relay: anonymised error reports, and the free ParaScape AI.

The relay is the Cloudflare Worker in ``services/parascape-relay`` (instructions 94 and 97). It
files an error report as an issue on ``EinarOlafsson/parascape`` -- a repeat of the same error adds
a comment to the open issue instead -- and answers questions from a free model, with no GitHub or
AI account on the user's side.

Nothing is sent without the user's say-so: `report_from_exception` only BUILDS a report, the
window shows exactly what it holds (`preview`), and `send_report` runs when the user clicks Send.
What a report holds is deliberately thin -- the exception type, its message with paths, e-mail
addresses and long tokens removed, file names (never paths) and function names from the
traceback, the ParaScape / Python / Qt versions and the operating system. The relay checks the
same rules again and refuses anything that still looks personal.

No window here, so every rule is testable.
"""
from __future__ import annotations

import json
import os
import platform
import re
import sys
import traceback
import urllib.error
import urllib.request
import uuid

#: The deployed relay. ``PARASCAPE_RELAY`` overrides it (an ``https://`` address, or ``http://``
#: on this computer for a local ``wrangler dev``).
RELAY_URL = "https://parascape-relay.einar-olafsson.workers.dev"
RELAY_ENV = "PARASCAPE_RELAY"
SCHEMA_VERSION = 1
TIMEOUT_S = 60

#: The relay's own limits (`services/parascape-relay/src/validate.js`), applied here first so a
#: report is never refused for being long.
MAX_MESSAGE = 300
MAX_NOTE = 1000
MAX_QUESTION = 2000
MAX_FRAMES = 40

_SETTINGS_KEY = "relay/install_id"
_CONSENT_KEY = "relay/send_errors"          # "ask" (default), "always", "never"


def relay_url() -> str:
    """The relay address: `RELAY_URL`, or a safe ``PARASCAPE_RELAY`` override."""
    override = os.environ.get(RELAY_ENV, "").strip().rstrip("/")
    if override.startswith("https://") or re.match(r"^http://(localhost|127\.0\.0\.1)(:\d+)?$",
                                                   override):
        return override
    return RELAY_URL


# ---------------------------------------------------------------------------- settings
def _settings():
    from .backdrops import _compat as AC
    return AC._settings()


def install_id() -> str:
    """A random id for this installation, made once and kept: the relay counts quotas by it.
    It identifies nothing about the person or the computer."""
    s = _settings()
    value = str(s.value(_SETTINGS_KEY, "") or "")
    try:
        uuid.UUID(value)
    except ValueError:
        value = str(uuid.uuid4())
        s.setValue(_SETTINGS_KEY, value)
    return value


def consent() -> str:
    """'ask' (show the report and ask), 'always' (send without asking) or 'never'."""
    value = str(_settings().value(_CONSENT_KEY, "ask") or "ask")
    return value if value in ("ask", "always", "never") else "ask"


def set_consent(value: str) -> None:
    """Remember 'ask', 'always' or 'never' (anything else means 'ask')."""
    _settings().setValue(_CONSENT_KEY, value if value in ("ask", "always", "never") else "ask")


# ---------------------------------------------------------------------------- anonymising
_PATH_RE = re.compile(r"(?:[A-Za-z]:)?(?:[\\/][^\s\\/'\"<>|:]+){2,}[\\/]?")
_EMAIL_RE = re.compile(r"[\w.+-]+@[\w-]+(?:\.[\w-]+)+")
_TOKEN_RE = re.compile(r"[A-Za-z0-9_+=-]{32,}")
_IP_RE = re.compile(r"(?<![\w.])(?:\d{1,3}\.){3}\d{1,3}(?![\w.])")
_HOME_RE = re.compile(r"~[\\/]\S*")


def scrub(text: str, limit: int | None = None) -> str:
    """Remove what could identify someone: paths (kept as their last name), home-relative paths,
    e-mail addresses, IP addresses and long tokens. Then one line, cut to `limit`."""
    text = str(text or "")
    text = _HOME_RE.sub("<path>", text)
    text = _PATH_RE.sub(lambda m: "<path>/" + re.split(r"[\\/]", m.group(0).rstrip("\\/"))[-1],
                        text)
    text = _EMAIL_RE.sub("<email>", text)
    text = _IP_RE.sub("<ip>", text)
    text = _TOKEN_RE.sub("<token>", text)
    user = os.environ.get("USER") or os.environ.get("USERNAME") or ""
    if len(user) >= 3:
        text = re.sub(re.escape(user), "<user>", text, flags=re.IGNORECASE)
    text = " ".join(text.split())
    return text[:limit].rstrip() if limit else text


def _safe(name: str, pattern: str, fallback: str) -> str:
    name = re.sub(r"[^A-Za-z0-9_.<> -]", "_", str(name or ""))[:100]
    return name if re.match(pattern, name) else fallback


# ---------------------------------------------------------------------------- environment
def _version(dist: str) -> str:
    try:
        from importlib.metadata import version
        return version(dist)
    except Exception:
        return ""


def parascape_version() -> str:
    """The installed ParaScape version, or 0.0.0 when it cannot be read."""
    for dist in ("parascape", "starplast"):
        v = _version(dist)
        if v:
            return v
    return "0.0.0"


def environment() -> dict:
    """The run, as the relay asks for it: versions and operating system, nothing else."""
    system = {"Linux": "linux", "Windows": "windows", "Darwin": "macos"}.get(platform.system(),
                                                                              "linux")
    try:
        from PySide6 import __version__ as qt
    except Exception:
        qt = None
    torch = _version("torch") or "none"
    gpu = "cpu"
    if "torch" in sys.modules:
        try:
            gpu = "cuda" if sys.modules["torch"].cuda.is_available() else "cpu"
        except Exception:
            gpu = "cpu"
    return {"schema_version": SCHEMA_VERSION, "install_id": install_id(),
            "parascape_version": parascape_version(), "os": system,
            "python": platform.python_version(), "qt": qt, "torch": torch, "gpu": gpu}


# ---------------------------------------------------------------------------- reports
def _module_of(frames, exc_type) -> str:
    """The ParaScape module the error surfaced in: the last parascape frame, else the
    exception's own module."""
    for fs in reversed(frames):
        path = fs.filename.replace("\\", "/")
        if "/parascape/" in path:
            rel = path.split("/parascape/", 1)[1].rsplit(".py", 1)[0].replace("/", ".")
            return _safe("parascape." + rel, r"^[A-Za-z0-9_][A-Za-z0-9_.:-]{0,79}$", "parascape")
    return _safe(getattr(exc_type, "__module__", "") or "parascape",
                 r"^[A-Za-z0-9_][A-Za-z0-9_.:-]{0,79}$", "parascape")


def report_from_exception(exc_type, exc, tb, note: str = "", screen: str = "") -> dict:
    """An anonymised report of one exception -- built, not sent."""
    frames = traceback.extract_tb(tb) if tb is not None else []
    rows = []
    for fs in list(frames)[-MAX_FRAMES:]:
        rows.append({"file": _safe(os.path.basename(fs.filename) or "<unknown>",
                                   r"^(?:<[A-Za-z0-9_ -]{1,40}>|[A-Za-z0-9_][A-Za-z0-9_.-]{0,99})$",
                                   "<file>"),
                     "function": _safe(fs.name, r"^(?:<[A-Za-z0-9_ -]{1,40}>|[A-Za-z_][A-Za-z0-9_.]{0,119})$",
                                       "<function>"),
                     "line": int(fs.lineno or 0)})
    report = environment()
    report.update({
        "module": _module_of(frames, exc_type),
        "exception_type": _safe(getattr(exc_type, "__name__", "Exception"),
                                r"^[A-Za-z_][A-Za-z0-9_.]{0,119}$", "Exception"),
        "message": scrub(str(exc), MAX_MESSAGE),
        "frames": rows,
    })
    if note.strip():
        report["user_note"] = scrub(note, MAX_NOTE)
    if screen:
        report["screen"] = _safe(screen, r"^[A-Za-z0-9][A-Za-z0-9 &()+.,'-]{0,59}$", "ParaScape")
    return report


def problem_report(note: str, screen: str = "") -> dict:
    """A report the user writes by hand (Help ▸ Report a problem), with no exception."""
    report = environment()
    report.update({"module": "parascape.app", "exception_type": "ProblemReport",
                   "kind": "user_error", "message": scrub(note, 800) or "(no description)",
                   "frames": [], "user_note": scrub(note, MAX_NOTE)})
    if screen:
        report["screen"] = _safe(screen, r"^[A-Za-z0-9][A-Za-z0-9 &()+.,'-]{0,59}$", "ParaScape")
    return report


def signature(report: dict) -> str:
    """What makes two reports the same error, for asking once per session."""
    frames = ";".join(f"{f['file']}:{f['function']}" for f in report.get("frames", [])[-5:])
    return f"{report.get('exception_type')}|{report.get('module')}|{frames}"


def preview(report: dict) -> str:
    """Exactly what would be sent, readable: the dialog shows this."""
    return json.dumps({k: v for k, v in report.items() if k != "install_id"}, indent=2)


# ---------------------------------------------------------------------------- sending
class RelayError(RuntimeError):
    """The relay refused or could not be reached; the message says which, in plain words."""


def _post(path: str, body: dict, timeout: float = TIMEOUT_S) -> dict:
    data = json.dumps(body).encode("utf-8")
    request = urllib.request.Request(relay_url() + path, data=data, method="POST",
                                     headers={"content-type": "application/json",
                                              "user-agent": f"parascape/{parascape_version()}"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as error:
        try:
            payload = json.loads(error.read().decode("utf-8"))
        except Exception:
            payload = {}
        if error.code == 429:
            raise RelayError("The daily limit for this has been reached; it resets at midnight UTC.")
        raise RelayError(payload.get("error") or f"the relay answered {error.code}")
    except (urllib.error.URLError, TimeoutError, OSError) as error:
        raise RelayError(f"could not reach the ParaScape relay ({type(error).__name__})")


def send_report(report: dict) -> dict:
    """File the report. Returns the relay answer (``issue_url``, ``deduplicated``...)."""
    return _post("/report", report)


def ask(question: str, context: str = "") -> str:
    """The free ParaScape AI's answer to a question (no exception attached)."""
    body = environment()
    text = scrub(question, MAX_QUESTION)
    if context:
        room = MAX_QUESTION - len(text) - 40
        if room > 200:
            text = f"{text}\n\nWhat is on screen: {scrub(context, room)}"
    body.update({"module": "", "exception_type": "", "message": "", "frames": [],
                 "question": text})
    answer = _post("/diagnose", body)
    return str(answer.get("advice") or "").strip() or "(no answer)"
