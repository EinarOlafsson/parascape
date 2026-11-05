"""A dropped file of genes: every format, and nothing silently dropped.

The expensive failure here is not a crash, it is a list that quietly shrinks: a user drops 300
identifiers, gets 250 genes and never learns which 50 went or why. So the tests that matter are that
an unresolved identifier is REPORTED, that an ambiguous one is not guessed at, and that an
unsupported extension says what is accepted instead of failing inside a parser.
"""
import pandas as pd
import pytest

from parascape import gene_drop as GD
from parascape import importer as IMP

KNOWN = ["TGME49_208830", "TGME49_233460", "TGME49_214320"]


@pytest.fixture
def table():
    return pd.DataFrame({"gene": ["TGME49_208830", "TgME49.233460", "not a gene",
                                  "TGME49_999999"],
                         "FPKM": [1200.0, 3.4, 0.0, 55.0]})


def _write(table, path, sheet="hits"):
    if str(path).endswith((".xlsx", ".xls")):
        with pd.ExcelWriter(path) as writer:
            table.to_excel(writer, sheet_name=sheet, index=False)
    elif str(path).endswith(".tsv"):
        table.to_csv(path, sep="\t", index=False)
    else:
        table.to_csv(path, index=False)
    return str(path)


@pytest.mark.parametrize("name", ["hits.csv", "hits.tsv", "hits.txt", "hits.xlsx"])
def test_a_dropped_file_of_each_format_reaches_the_node(tmp_path, table, name):
    drop = GD.read_drop(_write(table, tmp_path / name), known_ids=KNOWN)
    assert drop.name == name
    assert drop.gene_column == "gene"
    assert drop.genes == ["TGME49_208830"] or "TGME49_233460" in drop.genes
    assert "identifier column 'gene'" in drop.headline()


def test_an_excel_drop_reports_its_sheet_and_offers_the_others(tmp_path, table):
    path = tmp_path / "supplement.xlsx"
    with pd.ExcelWriter(path) as writer:
        pd.DataFrame({"legend": ["see methods"]}).to_excel(writer, sheet_name="legend", index=False)
        table.to_excel(writer, sheet_name="Table S3", index=False)
    first = GD.read_drop(path, known_ids=KNOWN)
    assert first.sheets == ("legend", "Table S3")
    assert first.sheet == "legend"                     # the published supplement's usual trap
    chosen = GD.read_drop(path, sheet="Table S3", known_ids=KNOWN)
    assert chosen.sheet == "Table S3"
    assert chosen.gene_column == "gene"
    assert "sheet Table S3" in chosen.headline()


def test_a_tab_separated_txt_is_read_as_tsv_and_the_separator_reported(tmp_path, table):
    path = tmp_path / "really_a_tsv.txt"
    table.to_csv(path, sep="\t", index=False)
    drop = GD.read_drop(path, known_ids=KNOWN)
    assert drop.separator == "\t"
    assert "separator tab" in drop.headline()
    assert list(drop.table.columns) == ["gene", "FPKM"]


def test_a_semicolon_csv_from_a_european_excel_is_not_read_as_one_column(tmp_path, table):
    path = tmp_path / "euro.csv"
    table.to_csv(path, sep=";", index=False)
    drop = GD.read_drop(path, known_ids=KNOWN)
    assert drop.separator == ";"
    assert len(drop.table.columns) == 2


def test_an_unresolved_identifier_is_reported_and_not_discarded(tmp_path, table):
    drop = GD.read_drop(_write(table, tmp_path / "hits.csv"), known_ids=KNOWN)
    report = drop.report
    assert len(report.rows) == len(table), "every input row must come back"
    raws = {r.raw: r for r in report.unresolved}
    assert "not a gene" in raws
    assert "no gene identifier" in raws["not a gene"].reason
    assert "TGME49_999999" in raws
    assert "does not hold" in raws["TGME49_999999"].reason
    assert "unresolved and listed below, not discarded" in report.summary()


def test_a_malformed_identifier_is_repaired_rather_than_lost():
    report = GD.resolution_report(["TgME49.233460", "TgME49-208830"], known_ids=KNOWN)
    assert report.genes == ["TGME49_233460", "TGME49_208830"]
    assert report.repair
    assert "reshaped" in report.summary()


def test_an_identifier_naming_two_genes_is_ambiguous_rather_than_guessed():
    report = GD.resolution_report(["TGME49_214320 / TGME49_208830"], known_ids=KNOWN)
    assert len(report.ambiguous) == 1
    assert report.ambiguous[0].candidates == ("TGME49_214320", "TGME49_208830")
    assert report.genes == []
    assert "ambiguous" in report.summary()


def test_repeated_identifiers_are_counted_rather_than_hidden():
    report = GD.resolution_report(["TGME49_208830", "TGME49_208830"], known_ids=KNOWN)
    assert report.genes == ["TGME49_208830"]
    assert report.duplicates == 1
    assert "repeated identifiers" in report.summary()


def test_without_a_node_table_only_the_identifier_form_is_checked_and_it_says_so():
    report = GD.resolution_report(["TGME49_999999"])
    assert report.genes == ["TGME49_999999"]
    assert "no node table to check against" in report.summary()


def test_an_unsupported_extension_gives_the_readable_list(tmp_path):
    path = tmp_path / "hits.docx"
    path.write_text("not a table")
    assert not GD.accepts(path)
    with pytest.raises(GD.UnsupportedFormat) as raised:
        GD.read_drop(path)
    message = str(raised.value)
    assert ".docx" in message
    for extension in IMP.READABLE:
        assert extension in message
    assert raised.value.accepted == tuple(IMP.READABLE)


def test_extra_columns_are_judged_by_their_range_and_never_by_their_name(tmp_path):
    """importer's own rule: a column reaching 16,520 has not been logged whatever it is called."""
    table = pd.DataFrame({"gene": KNOWN,
                          "log2_looking_name": [1200.0, 9000.0, 16520.0],
                          "fpkm_looking_name": [1.2, 4.0, 9.7]})
    drop = GD.read_drop(_write(table, tmp_path / "ranges.csv"), known_ids=KNOWN)
    by_column = {e.column: e for e in drop.extras}
    assert by_column["log2_looking_name"].treated_as == "fpkm"
    assert "16,520" in by_column["log2_looking_name"].why
    assert by_column["fpkm_looking_name"].treated_as == "log_intensity"
    assert "already logged" in by_column["fpkm_looking_name"].why
    assert all(e.line() for e in drop.extras)


def test_a_list_with_no_data_columns_is_still_a_list(tmp_path):
    path = tmp_path / "ids_only.txt"
    path.write_text("gene\n" + "\n".join(KNOWN) + "\n")
    drop = GD.read_drop(path, known_ids=KNOWN)
    assert drop.genes == KNOWN
    assert drop.extras == []


def test_an_empty_file_is_refused_with_a_reason(tmp_path):
    path = tmp_path / "empty.csv"
    path.write_text("gene\n")
    with pytest.raises(ValueError):
        GD.read_drop(path, known_ids=KNOWN)


def test_annotating_keeps_the_unresolved_rows_visible(tmp_path, table):
    drop = GD.read_drop(_write(table, tmp_path / "hits.csv"), known_ids=KNOWN)
    nodes = pd.DataFrame({"gene_id": KNOWN, "product": ["a", "b", "c"]})
    frame = GD.annotated(drop, nodes, columns=["product"])
    assert len(frame) == len(table)
    assert set(frame["state"]) == {"resolved", "unresolved"}
    assert frame.loc[frame["state"] == "unresolved", "reason"].notna().all()
