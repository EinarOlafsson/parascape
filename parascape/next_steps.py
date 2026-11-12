"""What to do next with the thing in your hand, offered the way spaCR offers it.

PORTED IN SHAPE FROM spaCR
==========================
spaCR puts a strip above the Run button of every pipeline screen
(``spacr/qt/chaining.py``, ``spacr.chaining.next_steps``): when a run finishes it lists
the modules that can work on what the run just produced, each **pre-filled** with it, each
run through ``spacr.ports.check_ready`` first, and **an offer that cannot work is still
shown -- greyed, with the reason it cannot and what to do about it**. Instruction 98 asks
for exactly that in ParaScape, and Einar asked for it in his own words: "when i click the
gene card once i should see the moduals downstream of one gene".

    Copyright (c) 2025-2026, Einar Birnir Olafsson and the spaCR contributors
    BSD 3-Clause with additional noncommercial terms
    SPDX-License-Identifier: LicenseRef-spaCR-BSD-3-Clause-Additional-Terms
    Full text: LICENSES/spaCR-BSD-3-Clause-Additional-Terms.txt

**What changed in porting.** spaCR pre-fills a FOLDER PATH and checks it against an
artifact registry on disk. ParaScape pre-fills a SUBJECT -- one gene, a list, a trait --
and checks it against :func:`parascape.ports.check_ready`, which needs nothing on disk.
There is no staleness row here, because a gene cannot go out of date the way a folder of
crops can. What is kept is the part Einar is looking at: ready offers first, blocked ones
last and greyed with their reason, and every offer carrying the subject so nothing has to
be typed twice.

A blocked offer is never hidden. Dropping it would leave somebody with no idea the module
exists, which is the failure spaCR records and the reason its own list keeps them.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

from PySide6 import QtCore, QtWidgets

from . import modules as MOD
from . import ports as PORTS
from . import theme as TH

#: What a blocked offer adds to its name, so the greying is also readable as words.
NOT_READY = " — not ready"


@dataclass(frozen=True)
class Offer:
    """One module offered as the next step, already checked.

    :param module: the module key offered.
    :param name: what the catalogue calls it.
    :param note: the one-line plain-words description from the catalogue.
    :param subject: the subject being handed on, e.g. ``"gene"``.
    :param value: what is in it -- the gene identifier -- so the offer pre-fills.
    :param shown: how to SAY what is in it, when the value itself is too long to read. A gene
        is its own name; a list of forty genes is "40 genes", and printing the forty
        accessions in a tooltip would be a wall of text where a count was wanted. Empty
        means the value reads well enough on its own, which is the case for one gene.
    :param readiness: the verdict of :func:`parascape.ports.check_ready`.
    :param reason: why it cannot be opened in this window at all, from
        :func:`parascape.modules.unavailable_reason`. Separate from `readiness`
        because they are different failures: "you have no gene" is yours to fix,
        "this view has not been built" is not.
    """

    module: str
    name: str
    note: str
    subject: str
    value: str
    readiness: PORTS.Readiness
    reason: str = ""
    shown: str = ""

    @property
    def words(self) -> str:
        """What is in hand, in words fit to read: :attr:`shown` when there is one."""
        return self.shown or self.value

    @property
    def ok(self) -> bool:
        """Whether this offer can be taken right now."""
        return bool(self.readiness.ok) and not self.reason

    @property
    def blocked(self) -> str:
        """Why it cannot be taken, in one line, or ``""`` when it can."""
        if self.reason:
            return self.reason
        if self.readiness.ok:
            return ""
        errors = self.readiness.errors
        if not errors:
            return "it cannot work on this"
        more = len(errors) - 1
        return errors[0].message + (f" (and {more} more)" if more else "")

    @property
    def fix(self) -> str:
        """What to do about the blockage, or ``""`` when there is none."""
        if self.reason:
            return ("Nothing to do here yet: this is a view ParaScape has not built, and it "
                    "is shown so you know it is coming rather than hidden.")
        errors = self.readiness.errors
        return errors[0].fix if errors else ""

    @property
    def label(self) -> str:
        """What the button says."""
        return self.name if self.ok else self.name + NOT_READY

    @property
    def tip(self) -> str:
        """The tooltip: what it would do, or why it cannot and what to do instead."""
        if self.ok:
            what = PORTS.describe_subject(self.subject)
            with_it = f" with {what} {self.words} already filled in" if self.value else ""
            # Only the first letter: `str.capitalize` lowercases the rest, which turned
            # "its RNA and protein" into "its rna and protein" on five tooltips.
            note = self.note[:1].upper() + self.note[1:]
            return f"{note}. Opens {self.name}{with_it}."
        return f"{self.blocked}\n\n{self.fix}".strip()


def offers(subject: str, value: str = "", *, without: str = "",
           available=None, shown: str = "") -> Tuple[Offer, ...]:
    """Every module that can take `subject`, checked, ready ones first.

    The list comes from :func:`parascape.modules.downstream_of_subject`, which reads the
    declared subjects of :mod:`parascape.ports` -- so it is derived from what the modules
    say they read, never from a list anybody keeps by hand.

    :param subject: the subject in hand, e.g. :data:`parascape.ports.GENE`.
    :param value: what is in it. Empty is allowed and is a different answer: the offer is
        shown with "none is chosen" rather than hidden, because a person looking at an
        empty gene card should still see what a gene would let them do.
    :param without: a module key to leave out, normally the one making the offer.
    :param available: ``{key: reason}`` for the modules this window cannot open, from
        `landing_install._reasons`. Omitted, the catalogue answers on its own.
    :param shown: how to say what is in hand when the value is too long to read -- "40
        genes" for a list. See :attr:`Offer.shown`.
    :returns: one :class:`Offer` per module, ready ones first.
    """
    reasons = dict(available or {})
    out = []
    for key in MOD.downstream_of_subject(subject, without=without):
        module = MOD.get(key)
        held = PORTS.Subject(subject, value=value or None)
        out.append(Offer(module=key, name=module.name, note=module.note, subject=subject,
                         value=str(value or ""), shown=str(shown or ""),
                         readiness=PORTS.check_ready(key, held),
                         reason=reasons.get(key, MOD.unavailable_reason(module))))
    return tuple(sorted(out, key=lambda o: not o.ok))


class NextSteps(QtWidgets.QWidget):
    """The row of next steps: spaCR's chaining strip, for a ParaScape subject.

    One line saying what is in hand, then a button per module that can take it. A blocked
    offer keeps its place in the row, greyed and reading "not ready", and says why in its
    tooltip. Clicking a ready one emits :attr:`chosen` with the module key and the subject,
    which is how the gene gets carried across without being typed again.
    """

    #: ``(module key, subject kind, value)`` when a ready offer is clicked.
    chosen = QtCore.Signal(str, str, str)

    def __init__(self, subject: str, parent=None, lead: str = "Next, with this gene:"):
        super().__init__(parent)
        self.subject = subject
        self._offers: Tuple[Offer, ...] = ()
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(8, 6, 8, 8)
        layout.setSpacing(6)
        self.lead = QtWidgets.QLabel(lead)
        self.lead.setObjectName("ShellMuted")
        self.lead.setWordWrap(True)
        self.lead.setToolTip(TH.tip(
            "Every module that reads " + PORTS.describe_subject(subject) + ". Each opens with "
            "what you are looking at already filled in. One that cannot work is still listed, "
            "greyed, and its tooltip says why and what would fix it."))
        layout.addWidget(self.lead)
        self._row = QtWidgets.QWidget(self)
        self._row_layout = QtWidgets.QHBoxLayout(self._row)
        self._row_layout.setContentsMargins(0, 0, 0, 0)
        self._row_layout.setSpacing(8)
        self._row_layout.addStretch(1)
        layout.addWidget(self._row)
        self.buttons: dict = {}

    def set_subject(self, value: str, *, available=None, without: str = "",
                    shown: str = "") -> Tuple[Offer, ...]:
        """Redraw the offers for one value of the subject; returns what was offered.

        :param value: the gene (or other subject) in hand; empty means none chosen yet.
        :param available: ``{key: reason}`` for modules this window cannot open.
        :param without: a module key to leave out, normally the one offering.
        :param shown: how to say what is in hand when the value is too long to read -- "40
            genes" rather than forty accessions. See :attr:`Offer.shown`.
        """
        self._offers = offers(self.subject, value, without=without, available=available,
                              shown=shown)
        while self._row_layout.count() > 1:
            item = self._row_layout.takeAt(0)
            if item.widget() is not None:
                item.widget().setParent(None)
                item.widget().deleteLater()
        self.buttons = {}
        for offer in self._offers:
            button = QtWidgets.QPushButton(offer.label, self._row)
            button.setEnabled(offer.ok)
            button.setToolTip(TH.tip(offer.tip))
            button.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
            button.clicked.connect(
                lambda _checked=False, o=offer: self.chosen.emit(o.module, o.subject, o.value))
            self._row_layout.insertWidget(self._row_layout.count() - 1, button)
            self.buttons[offer.module] = button
        ready = sum(1 for o in self._offers if o.ok)
        what = PORTS.describe_subject(self.subject)
        if value:
            self.lead.setText(f"With {what} <b>{shown or value}</b> in hand, {ready} of "
                              f"{len(self._offers)} of these can run now:")
        else:
            # Some of these read a gene they can do without, so "none of them can run" would
            # be false. What is true is that none of them will be filled in for you.
            self.lead.setText(f"Choose {what} first and these open with it already filled in. "
                              f"{ready} of {len(self._offers)} can open without one:")
        return self._offers

    @property
    def offers(self) -> Tuple[Offer, ...]:
        """What is currently on offer, in the order drawn."""
        return self._offers
