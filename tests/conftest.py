"""Pin the Qt binding before anything imports pyqtgraph.

parascape imports PySide6 directly, and pyqtgraph binds itself to whichever binding it finds in
sys.modules first. When the dependency was PyQt6 this was a live hazard: pytest-qt imported PySide6 as its
plugin loaded, pyqtgraph bound to that, the app kept using PyQt6, and the two halves of the GL widget came
from different bindings -- so the app failed to import under pytest while working perfectly when launched
normally. On PySide6 both halves agree by default, but pinning is kept: the failure was silent, load-order
dependent, and cost an afternoon, and nothing says a future dependency cannot drag the other binding in.

Setting PYQTGRAPH_QT_LIB makes pyqtgraph's choice explicit instead of load-order dependent.
"""
import os

os.environ.setdefault("PYQTGRAPH_QT_LIB", "PySide6")
os.environ.setdefault("PYTEST_QT_API", "PySide6")
# The GL widget needs a platform plugin; offscreen keeps the suite headless.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
# And a SOFTWARE GL context, which the suite has always used by accident and now uses on purpose.
# pyqtgraph 0.14 parses the driver's version string as a float, so an NVIDIA context reporting
# `4.6.0 NVIDIA 580.173.02` raises "Requires >= OpenGL 2.1" on a driver that plainly exceeds it --
# and every GL test then fails on a machine that has a card, while passing on one that does not.
# Software GL also makes the render comparisons reproducible, which is what they are for.
os.environ.setdefault("LIBGL_ALWAYS_SOFTWARE", "1")
# Mesa reads the variable above; a machine whose GLX vendor library is NVIDIA's does not, so Qt is
# told directly as well. Both, because either one alone leaves the suite passing or failing
# according to which graphics stack the machine happens to have.
os.environ.setdefault("QT_OPENGL", "software")
os.environ.setdefault("__GLX_VENDOR_LIBRARY_NAME", "mesa")

# Cap the threads of OpenMP and BLAS, because without a cap one test can take over ten minutes. Measured on
# 2026-11-12 on a 24-CPU machine at a load average of 70: `test_boosting_answers_the_same_question_template...`
# (a 400-row histogram gradient boosting over 24 classes) took 7 seconds with four threads and was
# still inside the first fit after ten minutes with the default of 24. OpenMP threads wait for each other
# at a barrier, so when the machine has more runnable threads than CPUs every barrier waits for a thread
# that is not running, and a call that is short on a quiet machine becomes unbounded on a busy one. A
# test of wiring gains nothing from 24 threads. `setdefault`, so a person who wants more sets it.
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")


# Preferences persist through QSettings, and a test run must not write into the settings file of
# whoever is running it -- nor read one, since a developer with logging enabled would otherwise have
# every Window in the suite writing to their real log.
#
# The path is set for BOTH formats, which is the part that matters and that a first attempt got
# wrong. On Unix NativeFormat IS the ini format, so a QSettings built from (organisation,
# application) resolves as NativeFormat whatever `setDefaultFormat` says, and setting the path for
# IniFormat alone left every test writing to ~/.config. It did: a test's temporary directory ended
# up in the real settings file and the application then started with a log directory that no longer
# existed. `test_settings_are_isolated_from_the_user_running_the_suite` is the check that this holds.
import tempfile  # noqa: E402

from PySide6 import QtCore  # noqa: E402

SETTINGS_DIR = tempfile.mkdtemp(prefix="parascape-test-settings-")
QtCore.QSettings.setDefaultFormat(QtCore.QSettings.Format.IniFormat)
for _fmt in (QtCore.QSettings.Format.NativeFormat, QtCore.QSettings.Format.IniFormat):
    QtCore.QSettings.setPath(_fmt, QtCore.QSettings.Scope.UserScope, SETTINGS_DIR)
    QtCore.QSettings.setPath(_fmt, QtCore.QSettings.Scope.SystemScope, SETTINGS_DIR)

# The same rule for the user's own state, and for the same reason. Every `Window` the suite builds
# keeps its clusterings in `paths.user_cache_dir()/runs` and its annotations beside them, so the
# tests were writing into the saved runs of whoever ran them: a real store on this machine had
# `test_run_a`, `test_two_b` and `before` sitting among a user's own work, and any test that keeps a
# run left another one behind. Saved embeddings go the same way, into the package's own directory.
#
# `PARASCAPE_STATE` points all three somewhere temporary for the duration of the run.
# `test_the_user_s_saved_runs_are_isolated_from_the_suite` is the check that this holds.
STATE_DIR = tempfile.mkdtemp(prefix="parascape-test-state-")
os.environ["PARASCAPE_STATE"] = STATE_DIR


# ONE application for the whole session. Twenty-two test modules create their own with
# `QApplication.instance() or QApplication([])`; when a module's last reference went, PyQt destroyed
# that application and the next module built a second one. Qt keeps style state in globals that
# outlive an application, so a style installed on the first (the glass style) left the second
# pointing at freed memory, and `setStyleSheet` segfaulted -- only in full runs, where garbage
# collection happened to fall between the two. The desktop program never builds two applications;
# the suite now does not either.
import pytest  # noqa: E402


@pytest.fixture(scope="session", autouse=True)
def _one_qt_application():
    from PySide6 import QtWidgets
    if QtWidgets.QApplication.instance() is None:
        # What `parascape.app` does at import, and it must precede the application: a shared GL
        # context group and a surface format pyqtgraph accepts.
        from parascape import sprite
        sprite.ensure_gl_format()
        QtCore.QCoreApplication.setAttribute(QtCore.Qt.ApplicationAttribute.AA_ShareOpenGLContexts)
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    yield app


@pytest.fixture(scope="module", autouse=True)
def _close_windows_after_each_module(_one_qt_application):
    """Delete every top-level window a test module leaves behind, once the module is done.

    Destroying the application used to do this as a side effect. With one application for the
    session, leftover windows would otherwise accumulate, and every stylesheet change re-polishes
    all of them: the suite slowed several-fold before this existed.
    """
    yield
    from PySide6 import QtCore, QtWidgets
    app = _one_qt_application
    for widget in app.topLevelWidgets():
        widget.close()
        widget.deleteLater()
    QtCore.QCoreApplication.sendPostedEvents(None, QtCore.QEvent.Type.DeferredDelete)
    app.processEvents()
