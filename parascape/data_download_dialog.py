"""Offer the versioned data download before creating a data-dependent main window."""
from __future__ import annotations

from PyQt6 import QtCore, QtWidgets

from . import data_manager as D


class _Download(QtCore.QThread):
    progress = QtCore.pyqtSignal(str, int)
    completed = QtCore.pyqtSignal(bool, str)

    def run(self):
        try:
            D.download(progress=lambda phase, done, total: self.progress.emit(phase, int(100 * done / max(total, 1))),
                       cancelled=self.isInterruptionRequested)
            self.completed.emit(True, 'Data downloaded and verified.')
        except Exception as exc:
            self.completed.emit(False, str(exc))


class DataDownloadDialog(QtWidgets.QDialog):
    """A cancellable download with visible destination, progress, errors and retry."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Download ParaScape data')
        self.resize(520, 260)
        self.worker = None
        self._cancelling = False
        layout = QtWidgets.QVBoxLayout(self)
        state = D.status()
        label = QtWidgets.QLabel(f'ParaScape needs its data to open. Download the data for version {state["version"]} '
                                 f'from GitHub ({state["download_bytes"] / 1024**2:.0f} MiB)?\n\n'
                                 f'Data folder: {state["directory"]}\n\n'
                                 'Once installed, browsing works offline. Existing data is preserved if the download fails.')
        label.setWordWrap(True)
        layout.addWidget(label)
        self.message = QtWidgets.QLabel('Ready to download.')
        self.message.setWordWrap(True)
        layout.addWidget(self.message)
        self.bar = QtWidgets.QProgressBar()
        self.bar.setValue(0)
        layout.addWidget(self.bar)
        buttons = QtWidgets.QHBoxLayout()
        self.start = QtWidgets.QPushButton('Download data')
        self.start.setToolTip('Download the versioned public GitHub data asset and verify its checksums before opening ParaScape.')
        self.start.clicked.connect(self._start)
        buttons.addWidget(self.start)
        self.cancel = QtWidgets.QPushButton('Cancel')
        self.cancel.setToolTip('Cancel the download and close ParaScape; a partial download can resume next time.')
        self.cancel.clicked.connect(self.reject)
        buttons.addWidget(self.cancel)
        layout.addLayout(buttons)

    def _start(self):
        self.start.setEnabled(False)
        self.worker = _Download(self)
        self.worker.progress.connect(self._progress)
        self.worker.completed.connect(self._completed)
        self.worker.start()

    def _progress(self, phase, percent):
        self.message.setText(phase)
        self.bar.setValue(percent)

    def _completed(self, success, message):
        self.worker.wait()
        self.message.setText(message)
        if success:
            self.accept()
        else:
            self.start.setText('Retry download')
            self.start.setEnabled(True)

    def reject(self):
        """Wait for an interrupted download to leave its staging files safely."""
        if self.worker is not None and self.worker.isRunning():
            if self._cancelling:
                return
            self._cancelling = True
            self.worker.requestInterruption()
            self.message.setText('Cancelling…')
            self.cancel.setEnabled(False)
            self.worker.completed.disconnect(self._completed)
            self.worker.finished.connect(super().reject)
            return
        super().reject()


def ensure_data() -> bool:
    """Return immediately when data is present; otherwise ask before downloading."""
    if D.status()['ready']:
        return True
    return DataDownloadDialog().exec() == QtWidgets.QDialog.DialogCode.Accepted
