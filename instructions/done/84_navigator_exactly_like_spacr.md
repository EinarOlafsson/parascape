# 84 · The Navigator, exactly as spaCR's

**Status: built 2026-11-06, waiting on the gate (Einar compares it with spaCR).** All five plan steps are
in: View ▸ Navigator and Help ▸ Navigator (checkable, remembered in QSettings `landing/navigator`); search
on Return only; 224 px square cards with a 96 px spaCR icon and an 18 px weight-300 name; details panel
below the row on hover; grow-in animation with faded row edges. Tests: `tests/test_navigator.py`,
`tests/test_modules_in_window.py::test_the_navigator_can_be_turned_off_from_view_and_help`.

Requested 2026-11-06. Einar, after using it: "I couldn't turn off the navigator in View
or Tools like I can in spaCR ... when I type things in the navigator they instantly appear, they should wait
until I hit enter, and the cards that appear do not look at all like those in spaCR, and I want them to look
exactly the same (I mean be in exactly the same style)."

## What spaCR does, measured from `spacr/qt/widgets/home.py` and `spacr/qt/app.py`

| | spaCR | ParaScape now |
|---|---|---|
| Turning it off | a checkable **Navigator** action (`EasyModeAction`); off, the home is the tile page | no switch at all |
| When suggestions appear | after an input attempt, not per keystroke | on every keystroke (180 ms debounce) |
| Card shape | **square, 224 px**, 20 px gap, row centred | 220 × 104 px rectangles |
| Card content | module **icon** centred on top, **name centred, 18 px, weight 300**, no description | coloured KIND label, left-aligned title, two-line subtitle |
| Description | in the details panel **below the row**, shown on hover (`explain_card`) | on the card, and in a pane |
| Arrival | cards **grow in** (geometry animation, masked to the row, edges faded) | appear at once |
| Hover colour | by module maturity (`_stage_button_rules`) | accent rim only |
| Card QSS | surface fill, 1 px border, radius 16, padding 12 | similar values, different layout |

## The plan
1. **A Navigator switch in View** (and in Help, where spaCR keeps it), checkable, remembered. Off: the landing
   opens on the module tiles and the Navigator page and its switch disappear. On: as now.
2. **Search on Return.** Typing changes nothing below the field; Return shows the cards; Return again, or a
   click, opens the highlighted card. The feedback line says "Press Return to search" while there is text
   that has not been searched.
3. **Port `_EasySuggestionCard`'s layout and style**: 224 px squares, icon above an 18 px light centred name,
   no description on the card; the details panel below the row shows the hovered or highlighted card's
   description; the row centred with spaCR's 20 px gap; spaCR's card QSS.
4. **Icons for every kind**, since a gene, a dataset and a question need a mark where spaCR has module icons.
5. **The grow-in animation and the faded row edges**, if they port without spaCR's ambient-field coupling;
   otherwise recorded here as not done and why.

## Gate
Einar looks at it beside spaCR and says the cards match.
