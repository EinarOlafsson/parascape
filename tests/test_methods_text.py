"""Methods paragraph, references and BibTeX from the dataset registry (instruction 60, change 2)."""
from __future__ import annotations

from parascape import datasets as D
from parascape import methods_text as M


def _confirmed_column():
    d = next(d for d in D.REGISTRY if d.citation and d.columns and d.organism == "Tg" and not d.derived_from)
    return d, d.columns[0]


def _unconfirmed_column():
    d = next(d for d in D.REGISTRY if not d.citation and d.columns and d.organism == "Tg" and not d.derived_from)
    return d, d.columns[0]


def test_every_citation_comes_from_the_registry_entry_of_a_column_used():
    d, column = _confirmed_column()
    src = M.sources([column], "Tg")
    assert [s.dataset.key for s in src.found] == [d.key]
    refs = M.reference_list(src)
    assert refs == [f"[1] " + "; ".join(x for x in (d.citation, f"PMID {d.pmid}" if d.pmid else "",
                                                     d.accession or "") if x)]


def test_an_unconfirmed_source_is_flagged_and_never_cited():
    d, column = _unconfirmed_column()
    src = M.sources([column], "Tg")
    assert src.numbered() == {}
    assert "confirm before citing" in M.reference_list(src)[0]
    assert "@misc" not in M.bibtex(src)
    assert "not yet confirmed" in M.methods_paragraph([column], "Tg")


def test_a_column_with_no_source_is_named_not_hidden():
    src = M.sources(["no_such_column_anywhere"], "Tg")
    assert src.unregistered == ["no_such_column_anywhere"]
    assert "1 further column has no recorded source" in M.methods_paragraph(["no_such_column_anywhere"], "Tg")


def test_the_bibtex_holds_only_fields_the_registry_records():
    d, column = _confirmed_column()
    bib = M.bibtex(M.sources([column], "Tg"))
    fields = {line.split("=")[0].strip() for line in bib.splitlines() if "=" in line}
    assert fields <= {"note", "pmid", "howpublished", "url"}
    assert "title" not in bib and "doi" not in bib and "journal" not in bib


def test_a_derived_layer_cites_the_measurements_it_was_computed_from():
    derived = next((d for d in D.REGISTRY if d.derived_from and d.columns and d.organism == "Tg"), None)
    if derived is None:
        return
    src = M.sources([derived.columns[0]], "Tg")
    parents = {D.provenance(c, "Tg").key for c in derived.derived_from if D.provenance(c, "Tg") is not None}
    assert parents <= {s.dataset.key for s in src.found}


def test_the_paragraph_states_strategy_settings_and_seed():
    _d, column = _confirmed_column()
    text = M.methods_paragraph([column], "Tg", strategy="feature_knn", settings={"k": 15}, seed=3)
    assert "feature_knn with k = 15, random seed 3." in text
