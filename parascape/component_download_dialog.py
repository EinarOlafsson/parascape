"""Offer the missing components of one dataset while keeping downloads cancellable."""
from PyQt6 import QtCore, QtWidgets

from . import data_components as C


class _Download(QtCore.QThread):
    progress = QtCore.pyqtSignal(str, int)
    completed = QtCore.pyqtSignal(bool, str)

    def __init__(self, components, parent=None):
        super().__init__(parent)
        self.components = components

    def run(self):
        try:
            for number, component in enumerate(self.components, 1):
                if self.isInterruptionRequested():
                    raise InterruptedError('Dataset download cancelled; verified files retained')
                caption = f'Data component {number} of {len(self.components)}'
                C.download(component, cancelled=self.isInterruptionRequested,
                           progress=lambda phase, done, total: self.progress.emit(
                               caption + ': ' + phase, int(100 * done / max(total, 1))))
            self.completed.emit(True, 'Dataset files downloaded and verified.')
        except Exception as error:
            self.completed.emit(False, str(error))


class ComponentDownloadDialog(QtWidgets.QDialog):
    """Show size, destination, progress, cancellation and retry for optional data."""

    def __init__(self, components, parent=None):
        super().__init__(parent)
        self.components = list(components)
        for component in self.components:
            C.validate(component)
        self.worker = None
        self._cancelling = False
        self.setWindowTitle('Download selected dataset')
        self.resize(560, 270)
        layout = QtWidgets.QVBoxLayout(self)
        size = sum(c['archive_size'] for c in self.components)
        description = QtWidgets.QLabel(
            f'Download or verify the selected dataset from GitHub?\n\n'
            f'Up to {size / 1024**2:.1f} MiB across {len(self.components)} components. '
            'Verified shared files are reused. After installation, this dataset works offline.\n\n'
            f'Data folder: {C.storage_root()}')
        description.setWordWrap(True)
        layout.addWidget(description)
        self.message = QtWidgets.QLabel('Ready to download.')
        self.message.setWordWrap(True)
        layout.addWidget(self.message)
        self.bar = QtWidgets.QProgressBar()
        self.bar.setValue(0)
        layout.addWidget(self.bar)
        buttons = QtWidgets.QHBoxLayout()
        self.start = QtWidgets.QPushButton('Download / verify dataset')
        self.start.setToolTip('Download only these source files and shared metadata, verifying all checksums before use.')
        self.start.clicked.connect(self._start)
        buttons.addWidget(self.start)
        self.cancel = QtWidgets.QPushButton('Cancel')
        self.cancel.setToolTip('Cancel safely; verified files and resumable parts are kept for a later attempt.')
        self.cancel.clicked.connect(self.reject)
        buttons.addWidget(self.cancel)
        layout.addLayout(buttons)

    def _start(self):
        if self.worker is not None and self.worker.isRunning():
            return
        self.start.setEnabled(False)
        self.worker = _Download(self.components, self)
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
        if self.worker is not None and self.worker.isRunning():
            if self._cancelling:
                return
            self._cancelling = True
            self.worker.requestInterruption()
            self.message.setText('Cancelling…')
            self.start.setEnabled(False)
            self.cancel.setEnabled(False)
            self.worker.completed.disconnect(self._completed)
            self.worker.finished.connect(super().reject)
            return
        super().reject()
