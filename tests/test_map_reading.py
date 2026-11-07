"""Reading a map: neighbours and how far to trust them (the Map module, instruction 96)."""
import numpy as np

from parascape import map_reading as MR


def _two_clusters():
    rng = np.random.default_rng(0)
    a = rng.normal(0, 0.1, (30, 3))
    b = rng.normal(5, 0.1, (30, 3))
    xyz = np.vstack([a, b])
    labels = ["dense granule"] * 30 + ["rhoptry"] * 30
    return xyz, labels


def test_a_gene_is_read_by_its_own_cluster():
    xyz, labels = _two_clusters()
    hood = MR.neighbours(xyz, labels, [3])
    assert hood.members == [3] and 3 not in hood.neighbours
    assert hood.top == ("dense granule", 10, 10)


def test_a_list_is_read_by_what_surrounds_it_not_by_itself():
    xyz, labels = _two_clusters()
    hood = MR.neighbours(xyz, labels, list(range(5)))
    assert not set(hood.neighbours) & set(range(5))
    assert hood.top[0] == "dense granule"


def test_unlabelled_neighbours_are_counted_apart():
    xyz, labels = _two_clusters()
    labels = ["" if i % 2 else v for i, v in enumerate(labels)]
    hood = MR.neighbours(xyz, labels, [0])
    assert hood.unlabelled > 0
    assert all(label for label, _ in hood.counts)


def test_unplaced_genes_are_neither_selected_nor_neighbours():
    xyz, labels = _two_clusters()
    placed = np.ones(60, bool)
    placed[1:10] = False
    hood = MR.neighbours(xyz, labels, [0, 1], placed=placed)
    assert hood.members == [0]
    assert not set(hood.neighbours) & set(range(1, 10))


def test_agreement_is_reported_against_chance():
    xyz, labels = _two_clusters()
    result = MR.agreement(xyz, labels)
    assert result.accuracy == 1.0 and result.chance == 0.5 and result.n == 60
    assert result.for_label("rhoptry") == (1.0, 30)


def test_a_map_that_mixes_labels_scores_near_chance():
    rng = np.random.default_rng(1)
    xyz = rng.normal(0, 1, (200, 3))
    labels = list(rng.choice(["a", "b"], 200))
    acc, chance = MR.label_agreement(xyz, labels)
    assert abs(acc - 0.5) < 0.15 and 0.45 < chance < 0.6


def test_too_few_labelled_genes_gives_no_number():
    xyz, labels = _two_clusters()
    assert MR.agreement(xyz[:5], labels[:5]) is None


# ------------------------------------------------------------------ groups
def _groups():
    # group 0: 8 dense granule + 2 rhoptry + 3 unknown; group 1: 10 rhoptry; noise: 5 dense granule
    labels = (["dense granule"] * 8 + ["rhoptry"] * 2 + [""] * 3
              + ["rhoptry"] * 10 + ["dense granule"] * 5)
    groups = [0] * 13 + [1] * 10 + [-1] * 5
    return groups, labels


def test_a_group_is_read_by_precision_and_recall():
    groups, labels = _groups()
    rows = {r.group: r for r in MR.group_rows(groups, labels)}
    g0 = rows[0]
    assert (g0.label, g0.hits, g0.labelled, g0.unknown) == ("dense granule", 8, 10, 3)
    assert g0.precision == 0.8 and abs(g0.recall - 8 / 13) < 1e-9
    g1 = rows[1]
    assert g1.label == "rhoptry" and g1.precision == 1.0 and abs(g1.recall - 10 / 12) < 1e-9


def test_each_label_names_its_best_group():
    groups, labels = _groups()
    rows = {r.label: r for r in MR.label_rows(groups, labels)}
    assert rows["rhoptry"].group == 1 and rows["dense granule"].group == 0


def test_learning_leaves_each_gene_out_and_skips_the_unclustered():
    groups, labels = _groups()
    learned = MR.learn_from_groups(groups, labels)
    # group 0: 8 dense granules called right, 2 rhoptry called wrong; group 1: 10 right
    assert learned.called == 20 and abs(learned.accuracy - 18 / 20) < 1e-9
    assert learned.labelled == 25 and learned.chance == 13 / 25
    assert learned.candidates == 3, "group 0 is 80% pure with 10 known, so its 3 unknowns count"


def test_a_group_of_unknowns_is_kept_but_not_scored():
    rows = MR.group_rows([0, 0, 1, 1], ["", "", "a", "a"])
    unknown = next(r for r in rows if r.group == 0)
    assert unknown.labelled == 0 and unknown.unknown == 2 and unknown.label == ""


def test_a_label_split_over_groups_is_reported_as_fragmented():
    """64.14: groups needed for 80 % of the label's grouped genes, and its unclustered share."""
    import numpy as np
    from parascape import map_reading as MR
    groups = np.array([0, 0, 1, 1, 2, 2, -1, 3, 3, 3])
    labels = np.array(["a", "a", "a", "a", "a", "a", "a", "b", "b", "b"], dtype=object)
    rows = {r.label: r for r in MR.label_rows(groups, labels)}
    assert rows["a"].groups_for_80 == 3 and abs(rows["a"].unclustered - 1 / 7) < 1e-9
    assert rows["b"].groups_for_80 == 1 and rows["b"].unclustered == 0.0
