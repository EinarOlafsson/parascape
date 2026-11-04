"""The information map stays aligned, and a new column is one registration plus its evidence rows."""
import pandas as pd
import pytest

from parascape import information_map as im

GENES = ['G1', 'G2', 'G3']


def rows(*items):
    return [dict(zip(im.EVIDENCE_COLUMNS, item)) for item in items]


def test_three_tables_share_shape_and_empty_cells():
    ev = rows(('G1', 'host_facing', 'yes', 0.7, 'observed-proximity', 'r1', 'd'),
              ('G1', 'host_facing', 'yes', 0.95, 'observed-direct', 'r2', 'd'),
              ('G2', 'compartment_hyperlopit', 'micronemes', 0.9, 'observed-curated', 'r3', 'd'))
    values, probability, source = im.build(ev, GENES)
    assert values.shape == probability.shape == source.shape
    assert (values.isna() == probability.isna()).all().all() and (values.isna() == source.isna()).all().all()
    assert values.loc['G1', 'host_facing'] == True and probability.loc['G1', 'host_facing'] == 0.95  # noqa: E712
    assert source.loc['G1', 'host_facing'] == 'observed-direct'
    assert pd.isna(values.loc['G3', 'host_facing'])                 # unlabeled stays empty, never False


def test_adding_a_column_leaves_the_others_unchanged():
    base = rows(('G1', 'host_facing', 'yes', 0.95, 'observed-direct', 'r', 'd'),
                ('G2', 'function_go_bp', 'GO:1 x', 0.95, 'observed-direct', 'r', 'd'))
    before = im.build(base, GENES)
    im.register(im.Trait('test_essential', 'measure', 'test-only column'))
    try:
        after = im.build(base + rows(('G3', 'test_essential', 'False', 0.8, 'observed-screen', 'r', 'd')), GENES)
        for old, new in zip(before, after):
            pd.testing.assert_frame_equal(old, new[old.columns])
        assert after[0].loc['G3', 'test_essential'] == False        # noqa: E712  a measured negative is kept
    finally:
        im.TRAITS.pop('test_essential')


def test_terms_show_the_best_supported_term_and_count():
    ev = rows(('G1', 'function_go_mf', 'GO:2 b', 0.95, 'observed-direct', 'p1', 'd'),
              ('G1', 'function_go_mf', 'GO:1 a', 0.95, 'observed-direct', 'p2', 'd'),
              ('G1', 'function_go_mf', 'GO:1 a', 0.95, 'observed-direct', 'p3', 'd'))
    values, _, _ = im.build(ev, GENES)
    assert values.loc['G1', 'function_go_mf'] == 'GO:1 a' and values.loc['G1', 'function_go_mf_n'] == 2


def test_unregistered_traits_and_bad_probabilities_are_refused():
    with pytest.raises(ValueError):
        im.build(rows(('G1', 'no_such_trait', 'x', 0.5, 'observed-direct', 'r', 'd')), GENES)
    with pytest.raises(ValueError):
        im.build(rows(('G1', 'host_facing', 'yes', 1.5, 'observed-direct', 'r', 'd')), GENES)


def test_partner_summary_counts_both_directions():
    pairs = pd.DataFrame([('G1', 'G2', 'parasite', 0.97, 'observed-direct', 'x'),
                          ('G1', 'STAT6', 'host', 0.9, 'observed-curated', 'y')], columns=im.PAIR_COLUMNS)
    values, _, _ = im.build([], GENES, pairs=pairs)
    assert values.loc['G1', 'partner_n_ge50'] == 2 and values.loc['G2', 'partner_best'] == 'G1'
    assert pd.isna(values.loc['G3', 'partner_best'])


def test_the_fill_bar_is_ninety_percent():
    assert im.FILL_PRECISION == 0.90
