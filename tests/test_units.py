"""Units of numeric columns (instruction 64.11): curated with a basis, patterns said so, never guessed."""
from parascape import organisms as O, units as U


def test_every_declared_value_target_has_a_curated_unit_and_basis():
    for code in ("Tg", "Pf"):
        for column in O.get(code).numbers:
            unit, basis = U.unit_of(column)
            assert unit != U.UNRECORDED and basis, column


def test_patterns_say_they_are_patterns_and_unknown_columns_are_not_guessed():
    unit, basis = U.unit_of("brain168465_14d_lfc")
    assert unit == "log2 fold change" and "not curated" in basis
    assert U.unit_of("some_score")[0] == U.UNRECORDED
    assert "not recorded" in U.sentence("some_score")
    assert "Sidik" in U.sentence("fit_invitro_hff")
