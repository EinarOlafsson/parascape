"""The AI symbol, drawn and behaving as spaCR draws it.

Einar: "the ai symbol should work in the same way". In spaCR the symbol is not a picture and not a
button with a word on it. It is a ring with a filled dot inside, sitting at the right-hand end of
the search field, and it is a switch:

* **Click** turns the assistant on or off. The dot is only there while it is on, so the state is
  readable at a glance without a label.
* **The colour of the dot is the provider**: white for the free ParaScape AI, orange for Claude,
  green for Codex, blue for Gemini.
* **Typing makes the dot pulse**, larger for faster typing, and it settles back within a second.
  Nothing is sent by typing; the pulse only says that the field is listening.
* **While an answer is being written** the dot is circled by a marker (the free provider) or a
  glow (the others), so waiting does not look like nothing happening.
* **Right-click chooses the provider.** The first click, with no provider chosen yet, asks the
  same question instead of guessing.
* The motion stops when the symbol is hidden or Reduce motion is on, and it never runs while the
  assistant is off.

This module holds the drawing and the choice of provider. What the page does with an answer lives
in :mod:`parascape.navigator`, because that is where the answer is shown.
"""
from __future__ import annotations

import math

from PySide6 import QtCore, QtGui, QtWidgets

#: The colour of the dot for each provider of :data:`parascape.chat.PROVIDERS`. spaCR's values.
PROVIDER_COLORS = {"parascape": "#FFFFFF", "claude": "#DE7356", "codex": "#A8E6CF",
                   "gemini": "#4285F4"}

#: How long the dot keeps pulsing after the last key, in milliseconds.
PULSE_MS = 850

#: Edge of the square the symbol is drawn in. spaCR's size.
SIZE = 36


def _reduced_motion() -> bool:
    """Whether the person asked for no animation, so the symbol sits still.

    Read from the same preference the cards use, so one switch quiets both.
    """
    try:
        from .backdrops import _compat as AC
        return bool(AC.get_reduce_motion())
    except Exception:
        return False


class AiButton(QtWidgets.QPushButton):
    """The ring-and-dot switch at the end of the Navigator field.

    Checkable, keyboard-focusable, and painted by hand because the whole point is that it is not
    a labelled button. It owns only its own drawing; the page connects ``clicked`` to the turning
    on and off, and ``customContextMenuRequested`` to the choice of provider.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("NavigatorAiButton")
        self.setAccessibleName("Navigator AI")
        self.setCheckable(True)
        self.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        self.setStyleSheet("background: transparent; border: none; padding: 0; "
                           f"min-width: {SIZE}px; min-height: {SIZE}px;")
        self.setFixedSize(SIZE, SIZE)
        self.provider = "parascape"
        self.busy = False
        self._pulse = 0.0
        self._last_key = -10_000
        self._phase = 0.0
        self._clock = QtCore.QElapsedTimer()
        self._clock.start()
        self._motion = QtCore.QTimer(self)
        self._motion.setInterval(33)
        self._motion.timeout.connect(self._tick)
        self.toggled.connect(self._tick)

    def typed(self, _text: str = "") -> None:
        """Pulse the dot for a key just typed, larger the faster the typing.

        Connected to the field, so the symbol shows that the field is listening. It never sends
        anything.

        :param _text: the text of the field, ignored.
        """
        now = self._clock.elapsed()
        speed = max(0.0, 1.0 - (now - self._last_key) / 500.0)
        self._last_key = now
        self._pulse = 0.35 + 0.65 * speed
        self._tick()

    def set_busy(self, busy: bool) -> None:
        """Say whether an answer is being written, which makes the waiting motion run.

        :param busy: True from the moment a question is sent until the answer is complete or
            cancelled.
        """
        self.busy = bool(busy)
        self._tick()

    def set_provider(self, name: str) -> None:
        """Change the colour of the dot to that of a provider.

        :param name: a provider name from :data:`PROVIDER_COLORS`; an unknown one draws white.
        """
        self.provider = name
        self.update()

    def color(self) -> QtGui.QColor:
        """The colour the dot is filled with now: that of the chosen provider."""
        return QtGui.QColor(PROVIDER_COLORS.get(self.provider, PROVIDER_COLORS["parascape"]))

    def _tick(self, *_args) -> None:
        """Repaint, and run the timer only while there is something visible to animate.

        Idle symbols cost nothing: the timer is stopped when the assistant is off, when the symbol
        is hidden, when motion is reduced, and when the pulse has died away.
        """
        now = self._clock.elapsed()
        self.update()
        wanted = (self.isVisible() and self.isChecked() and not _reduced_motion()
                  and (self.busy or now - self._last_key < PULSE_MS))
        if wanted:
            if not self._motion.isActive():
                self._motion.start()
        else:
            self._motion.stop()

    def showEvent(self, event):                            # noqa: N802
        """Resume the motion only when the symbol becomes visible."""
        super().showEvent(event)
        self._tick()

    def hideEvent(self, event):                            # noqa: N802
        """Release the timer while another page is on screen."""
        self._motion.stop()
        super().hideEvent(event)

    def paintEvent(self, event):                           # noqa: N802
        """Draw the ring, and inside it the dot when on, with the waiting motion when busy."""
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        painter.translate(self.rect().center())
        color = self.color()
        reduced = _reduced_motion()
        hosted = self.provider == "parascape"
        if self.isChecked() and self.busy and not hosted and not reduced:
            painter.save()
            painter.rotate((self._clock.elapsed() * 0.18) % 360)
            painter.setBrush(QtCore.Qt.BrushStyle.NoBrush)
            for width, alpha in ((4.0, 18), (2.5, 35), (1.0, 90)):
                glow = QtGui.QColor(color)
                glow.setAlpha(alpha)
                painter.setPen(QtGui.QPen(glow, width, QtCore.Qt.PenStyle.SolidLine,
                                          QtCore.Qt.PenCapStyle.RoundCap))
                painter.drawArc(QtCore.QRectF(-12, -12, 24, 24), 20 * 16, 100 * 16)
            painter.restore()
        painter.setPen(QtGui.QPen(QtGui.QColor("white"), 2.5 if self.hasFocus() else 1.5))
        painter.setBrush(QtCore.Qt.BrushStyle.NoBrush)
        painter.drawEllipse(QtCore.QRectF(-14, -14, 28, 28))
        if self.isChecked():
            elapsed = max(0, self._clock.elapsed() - self._last_key)
            pulse = 0.0 if reduced else self._pulse * math.exp(-elapsed / 220.0)
            radius = (8.0 + 5.0 * pulse) * (1.0125 if hosted else 1.0)
            painter.setPen(QtCore.Qt.PenStyle.NoPen)
            painter.setBrush(color)
            painter.drawEllipse(QtCore.QRectF(-radius, -radius, radius * 2, radius * 2))
            if hosted and self.busy and not reduced:
                angle = self._clock.elapsed() * math.tau / 1800.0
                size = radius * 0.2
                x, y = 0.5 * math.cos(angle) * radius, 0.5 * math.sin(angle) * radius
                painter.setBrush(QtGui.QColor(0, 0, 0))
                painter.drawEllipse(QtCore.QRectF(x - size, y - size, size * 2, size * 2))
        painter.end()


class ProviderDialog(QtWidgets.QDialog):
    """Choose which assistant answers, saying plainly what is shared before it is chosen.

    spaCR asks this on the first click and on a right-click. The wording says what leaves the
    computer, because choosing a provider is the moment a person agrees to that.
    """

    def __init__(self, providers, current: str = "", parent=None):
        super().__init__(parent)
        self.setObjectName("NavigatorAiProviderDialog")
        self.setWindowTitle("Navigator AI")
        lay = QtWidgets.QVBoxLayout(self)
        explanation = QtWidgets.QLabel(
            "The assistant answers questions about ParaScape and about what is on screen. In the "
            "Navigator, a search that finds nothing goes to the assistant you choose here, and "
            "the answer appears below the field. It does not run an analysis from the field.\n\n"
            "The free ParaScape AI uses the hosted service and needs no account. Claude, Codex "
            "and Gemini use the command line tool you have installed and signed in to. What you "
            "type is shared with the assistant you choose. Check what it says before using it.")
        explanation.setWordWrap(True)
        explanation.setMinimumWidth(480)
        lay.addWidget(explanation)
        self.choice = QtWidgets.QComboBox()
        self.choice.setObjectName("NavigatorAiProviderChoice")
        self.choice.setAccessibleName("Choose the assistant")
        self.choice.setToolTip("Which assistant answers. One that is not installed says how to "
                               "install it when you ask it a question.")
        for p in providers:
            self.choice.addItem(p.label if p.available() else f"{p.label} (not installed)", p.name)
        at = self.choice.findData(current)
        if at >= 0:
            self.choice.setCurrentIndex(at)
        lay.addWidget(self.choice)
        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        buttons.button(QtWidgets.QDialogButtonBox.StandardButton.Ok).setText("Use AI")
        buttons.button(QtWidgets.QDialogButtonBox.StandardButton.Ok).setToolTip(
            "Turn the assistant on with this provider.")
        buttons.button(QtWidgets.QDialogButtonBox.StandardButton.Cancel).setToolTip(
            "Leave the assistant as it was.")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        lay.addWidget(buttons)

    def selected(self) -> str:
        """The name of the provider chosen in the list."""
        return str(self.choice.currentData())
