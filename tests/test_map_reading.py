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
