"""The ParaScape manuscript: figures and legends assembled from the results, not retyped.

Instruction 101, item A8. Every number in every legend is read from the `summary.json` the figure
script wrote, so a legend cannot drift from the figure beside it. Each figure lands in a table of
one column and two rows, the image on top and its legend beneath, which is the layout Einar asked
for on 2026-11-17 and which `brady_metabolic_classes/report/build_manuscript_docx.py` implements
for the bradyzoite paper; that module is imported rather than copied.

Run:  python report/build_parascape_manuscript.py
"""
from __future__ import annotations

import glob
import json
import os
import shutil
import subprocess
import xml.etree.ElementTree as ET
import zipfile
from datetime import date

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: WordprocessingML namespace.
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def newest(pattern: str) -> str:
    found = sorted(glob.glob(os.path.join(ROOT, "results", pattern)))
    if not found:
        raise SystemExit(f"missing result: {pattern}")
    return found[-1]


def summary(pattern: str) -> dict:
    return json.load(open(newest(pattern)))


def panel(pattern: str) -> str:
    """The drawn panel for a figure, or "" when none has been drawn yet.

    A figure whose result is measured but whose panel is not drawn must say so. Pointing the
    manuscript at the result TABLE instead, which the first version did for seven figures, puts a
    .tsv where an image belongs and pandoc quietly replaces it with its own description.
    """
    found = sorted(glob.glob(os.path.join(ROOT, "results", pattern)))
    return os.path.relpath(found[-1], ROOT) if found else ""


def figure(number: str, title: str, image: str, legend: str) -> str:
    """One figure: a title, the panel if it exists, and the legend, each its own paragraph.

    The blank lines around the image matter: without them pandoc folds the title, the image and
    the legend into ONE paragraph and the two-cell table cannot be built.
    """
    body = (f"![]({image}){{width=6.5in}}\n\n" if image
            else "*[panel not yet drawn; the result is measured and the numbers below are from "
                 "it]*\n\n")
    return (f"\n**Figure {number}. {title}**\n\n{body}{legend.strip()}\n\n")


def build() -> str:
    f1, f1b = summary("figure1_schematic_*/summary.json"), summary("corpus_inventory_*/summary.json")
    f2 = summary("control_comparability_*/summary.json")
    f3 = summary("coexpression_modules_*/summary.json")
    f4 = summary("what_structures_contrasts_*/summary.json")
    f5 = summary("bradyzoite_reach_*/summary.json")
    f6b = summary("what_sets_protein_abundance_*/summary.json")
    f6c = summary("which_data_explains_which_*/summary.json")
    f7 = summary("network_layers_*/summary.json")
    f8 = summary("strain_variation_by_class_*/summary.json")
    f9 = summary("species_transfer_*/summary.json")
    f10 = summary("host_response_shared_*/summary.json")
    f10sc = summary("dual_perturb_single_cell_*/BUT290A27*_summary.json")
    f11 = summary("filling_the_map_*/summary.json")
    f12 = summary("measured_absence_*/summary.json")

    cover = {c["category"]: c for c in f1b["coverage"]}
    worst = sorted(f12["by_category"], key=lambda r: -r["genes_never_measured"])[0]
    best11 = sorted(f11["best"], key=lambda r: -r["lift_over_baseline"])[0]
    place = next(r for r in f9["layers"] if r.get("measurable") and r["layer"].startswith("where"))
    fitness = next(r for r in f9["layers"] if r.get("measurable")
                   and r["layer"].startswith("whether"))
    control9 = next(r for r in f9["layers"] if r.get("measurable") and r["layer"].startswith("how long"))
    classes = {c["class"]: c for c in f8["classes"]}
    union = next(r for r in f7["layers"] if r["layer"] == "ALL LAYERS AT ONCE")
    single = max((r for r in f7["layers"] if r.get("enough") and "ALL LAYERS" not in r["layer"]),
                 key=lambda r: r["accuracy"])

    out = [f"% ParaScape: one gene space for the apicomplexan literature",
           f"% {date.today():%-d %B %Y}", ""]
    out.append("## Figures\n")

    out.append(figure(
        "1", "What ParaScape is, and how much of the genome it reaches",
        panel("figure1_schematic_*/panel_schematic.pdf"),
        f"**(a)** {f1['datasets']} datasets imported into one gene space of 8,140 genes and 362 "
        f"slots, made comparable across studies, and read through {f1['modules']} modules. "
        f"Counts are read from the dataset catalogue at draw time. "
        f"**(b)** The share of genes anything has ever measured, per category of knowledge. "
        f"Transcript reaches {cover['transcript']['share_measured']:.1%} and chromatin "
        f"{cover['chromatin']['share_measured']:.1%}, while function reaches "
        f"{cover['function']['share_measured']:.1%}, phenotype "
        f"{cover['phenotype']['share_measured']:.1%} and host response "
        f"{cover['host_response']['share_measured']:.1%}. "
        f"{f1b['columns_set_aside']} columns are set aside as derivations, flags or counts in "
        f"which absence was written as zero."))

    out.append(figure(
        "2", "Which studies can be compared at all",
        panel("control_comparability_*/panel_*.pdf"),
        f"Agreement between the untreated arms of {f2['studies']} independent studies, over "
        f"{f2['genes']:,} genes. Removing the cell-cycle axis raises the median agreement from "
        f"{f2['median_agreement_before']:.3f} to {f2['median_agreement_after']:.3f}. A random axis "
        f"of the same size and shape gains {f2['null_gain_mean']:+.3f} "
        f"(sd {f2['null_gain_sd']:.3f}, {f2['n_null']} draws), so the cell-cycle excess is "
        f"{f2['gain'] - f2['null_gain_mean']:+.3f} at p = {f2['p']:.4f}. A study whose control "
        f"profile carries one value is refused rather than scored."))

    out.append(figure(
        "3", "The co-expression modules, and what each one is",
        "",
        f"{f3['genes']:,} genes cut into {f3['modules']} modules on co-expression components from "
        f"57 independent studies. {f3['named']} modules hold a compartment above their own matched "
        f"baseline; the rest hold none, which is reported rather than forced. The baseline is a "
        f"random gene set matched to each module in size and in how much is already known about "
        f"its members ({f3['n_null']} draws), never the prior size of the compartment."))

    out.append(figure(
        "4", "The laboratory structures the corpus, not the perturbation",
        "",
        f"{f4['pairs']:,} pairs of stored contrasts. Two from the same study agree at a median "
        f"absolute Spearman of {f4['median_abs_rho_same_study']:.3f}, two from different studies at "
        f"{f4['median_abs_rho_different_study']:.3f} -- {f4['ratio']:.1f} times. Of "
        f"{f4['parsed_contrasts']} parsed contrasts over {f4['studies']} studies, "
        f"{len(f4['perturbations_repeated_across_studies'])} perturbation is stated by two "
        f"independent laboratories, so cross-laboratory replication cannot be asked of this corpus."))

    out.append(figure(
        "5", "The bradyzoite subtype nothing published reaches",
        panel("bradyzoite_reach_*/panel_reach_by_subtype.pdf"),
        f"Five subtypes of in vivo bradyzoite tested against all {f5['contrasts']} stored "
        f"within-study contrasts, {f5['n_null']:,} expression-matched random gene sets behind each. "
        f"{f5['reaching_any']} contrasts reach at least one subtype. Subtype "
        f"{f5['least_reached']} is reached by {f5['by_subtype'][f5['least_reached']]}, against "
        f"{max(f5['by_subtype'].values())} for the most reached. Each subtype is marked by "
        f"{list(f5['markers_per_subtype'].values())[0]} genes."))

    out.append(figure(
        "6", "Which measurement explains which other measurement",
        panel("coexpression_components_localization_*/panel_which_evidence_places_a_protein.pdf"),
        f"**(a)** The 13 coarse hyperLOPIT compartments recovered from RNA co-expression alone "
        f"against all allowed evidence, as cross-validated AUC. Five are placed as well or better "
        f"by transcription alone, and each is built in a wave during division. "
        f"**(b)** What sets protein abundance in *Plasmodium*: transcript abundance explains "
        f"{f6b['r2_by_layer']['transcript abundance']:.1%} of the variance on the "
        f"{f6b['genes_used']} genes carrying every layer, and every layer together "
        f"{list(f6b['r2_by_layer'].values())[-1]:.1%}. On the wider {f6b['wide_genes']}-gene set "
        f"transcript explains {f6b['wide_transcript_only_r2']:.1%}. Comparable with the human "
        f"literature, one averaged mRNA measure against protein gives r2 = "
        f"{f6b['single_mrna_vs_protein_pearson_r2']:.3f} over "
        f"{f6b['single_mrna_vs_protein_n']:,} genes. "
        f"**(c)** Every pair of categories: {f6c['pairs_above_5_points']} of {f6c['pairs']} "
        f"recover more than five points above their own shuffled baseline. Pairs in which the "
        f"predictors determine the target by arithmetic are excluded by name "
        f"({', '.join(f6c['structural_pairs_excluded'])})."))

    out.append(figure(
        "7", "The whole network against its best single layer",
        "",
        f"Compartment recovered from a gene's neighbours, layer by layer, over "
        f"{f7['genes_with_label']:,} labelled genes against a commonest-answer baseline of "
        f"{f7['baseline']:.1%}. The best single layer is {single['layer']} at "
        f"{single['accuracy']:.1%} on {single['answered']:,} genes; every layer at once is "
        f"{union['accuracy']:.1%} on {union['answered']:,} genes. Weighting each layer by what it "
        f"is worth, fitted on half the labelled genes and scored on the other half, gives "
        f"{f7['weighted_on_held_out']['accuracy']:.1%} against "
        f"{f7['unweighted_same_genes']['accuracy']:.1%} unweighted. The `{f7['excluded_layer']}` "
        f"layer is excluded because it joins genes that share the label being recovered."))

    out.append(figure(
        "8", "What differs between strains is what goes into the host",
        "",
        f"Nonsynonymous over synonymous SNPs across the sequenced strains, by class of gene, over "
        f"{f8['genes_with_a_ratio']:,} genes carrying a ratio. Dense granules reach "
        f"{classes['dense granules (hyperLOPIT)']['ratio']:.2f} and apical secretory "
        f"{classes['apical secretory (rhoptries, micronemes)']['ratio']:.2f}, against a "
        f"length-matched baseline near "
        f"{classes['dense granules (hyperLOPIT)']['length_matched_baseline']:.2f}; nucleus "
        f"{classes['nucleus']['ratio']:.2f} and cytosol {classes['cytosol']['ratio']:.2f} fall "
        f"below it. The null is matched gene by gene on coding length ({f8['n_null']:,} draws), "
        f"because a longer gene carries more SNPs of both kinds. The ratio is an index and not "
        f"dN/dS: the counts are not normalised by the number of sites at which each change could "
        f"occur."))

    out.append(figure(
        "9", "Which layers transfer between Toxoplasma and Plasmodium",
        "",
        f"{f9['pairs']['one_to_one']:,} one-to-one orthologue pairs; "
        f"{f9['pairs']['many_to_many_set_aside']} orthogroups are set aside as many-to-many rather "
        f"than averaged. Where the protein is transfers at {place['transfer']:.3f} against a "
        f"{place['baseline']:.3f} baseline; whether losing it costs fitness at "
        f"{fitness['transfer']:.3f}. Protein length transfers at {control9['transfer']:.3f}, which "
        f"is the control on the pairing: a sequence property must transfer between true "
        f"orthologues, so a layer that does not is reporting biology rather than a mapping fault."))

    out.append(figure(
        "10", "The host response is mostly private to one effector",
        "",
        f"**(a)** {f10['host_genes']:,} host genes moved by at least one of {f10['effectors']} "
        f"effectors: {f10['private_host_genes']:,} are moved by exactly one and "
        f"{f10['shared_host_genes']} by {f10['shared_cut']} or more. The translocon is the shared "
        f"part, overlapping itself at {f10['overlap_translocon_within']:.3f} against "
        f"{f10['overlap_cargo_within']:.3f} for cargo with cargo, which is the control on the "
        f"import. **(b)** In the single-cell data, {f10sc['cells_assigned']:,} of "
        f"{f10sc['cells']:,} cells carry one guide; of {f10sc['targets_tested']} targets with "
        f"enough cells, {f10sc['targets_moving_the_parasite']} changes the parasite total "
        f"transcriptional output. The assumption every pseudobulk host analysis rests on -- that "
        f"the parasite is otherwise unchanged -- therefore largely holds."))

    out.append(figure(
        "11", "What the map can fill, and what it cannot",
        "",
        f"{f11['held_out_predictions']:,} held-out predictions over {f11['targets']} targets and "
        f"{f11['strategies']} strategies, under frozen partitions. "
        f"{f11['targets_beating_baseline']} of {f11['targets_summarised']} targets are called "
        f"better than by guessing the commonest answer; the best is {best11['target']} at "
        f"{best11['right_when_answered']:.1%} against {best11['commonest_answer_baseline']:.1%}, "
        f"answering {best11['answer_rate']:.1%} of the time. Accuracy, answer rate and baseline "
        f"are reported together because any one of them alone can be made to look good."))

    out.append(figure(
        "12", "What nothing measures, and which experiment would close the most",
        panel("corpus_inventory_*/panel_coverage.pdf"),
        f"Across {f12['genes']:,} genes and {f12['categories']} categories of knowledge, built "
        f"from {f12['columns']} measurement columns. The largest gap is {worst['category']}: "
        f"{worst['genes_never_measured']:,} genes "
        f"({1 - worst['share_measured']:.1%}) have never been measured that way, so one assay of "
        f"that kind would give them a first measurement. A column with no missing value at all is "
        f"not a measurement of every gene but a count or a flag in which absence was written as "
        f"zero, and such columns are excluded."))

    return "".join(out)


def w(tag: str) -> str:
    return W + tag


def figures_into_tables(path: str) -> None:
    """Put every figure in a table of one column and two rows: image on top, legend beneath.

    The same layout the bradyzoite paper uses. It is implemented here rather than imported from
    that project, because importing that module RUNS it -- the first version of this script
    rebuilt the bradyzoite manuscript as a side effect of an import.

    The body is rebuilt in one linear pass: inserting tables while walking the same list shifts
    the indices underneath you and merges two figures into one table.
    """
    zin = zipfile.ZipFile(path)
    doc = ET.fromstring(zin.read("word/document.xml"))
    body = doc.find(w("body"))
    kids = list(body)

    def text_of(el):
        return "".join(t.text or "" for t in el.iter(w("t"))).strip()

    def is_para(el):
        return el.tag == w("p")

    def has_image(el):
        return is_para(el) and el.find(".//" + w("drawing")) is not None

    def is_title(el):
        return is_para(el) and text_of(el).startswith("Figure ")

    def make_table(top, bottom):
        tbl = ET.Element(w("tbl"))
        pr = ET.SubElement(tbl, w("tblPr"))
        tw = ET.SubElement(pr, w("tblW")); tw.set(w("w"), "0"); tw.set(w("type"), "auto")
        bd = ET.SubElement(pr, w("tblBorders"))
        for side in ("top", "left", "bottom", "right", "insideH", "insideV"):
            e = ET.SubElement(bd, w(side)); e.set(w("val"), "none"); e.set(w("sz"), "0")
        grid = ET.SubElement(tbl, w("tblGrid"))
        ET.SubElement(grid, w("gridCol")).set(w("w"), "9360")
        for paragraphs in (top, bottom):
            tr = ET.SubElement(tbl, w("tr"))
            tc = ET.SubElement(tr, w("tc"))
            tc_pr = ET.SubElement(tc, w("tcPr"))
            cw = ET.SubElement(tc_pr, w("tcW")); cw.set(w("w"), "9360"); cw.set(w("type"), "dxa")
            for paragraph in paragraphs:
                tc.append(paragraph)
            if not paragraphs:
                ET.SubElement(tc, w("p"))
        return tbl

    out, i, made = [], 0, 0
    while i < len(kids):
        el = kids[i]
        if is_title(el) and i + 1 < len(kids) and has_image(kids[i + 1]):
            title, image, j, legend = el, kids[i + 1], i + 2, []
            while (j < len(kids) and is_para(kids[j]) and not is_title(kids[j])
                   and not has_image(kids[j]) and text_of(kids[j])):
                legend.append(kids[j]); j += 1
            out.append(make_table([image], [title] + legend))
            made += 1
            i = j
        else:
            out.append(el); i += 1
    for child in list(body):
        body.remove(child)
    for child in out:
        body.append(child)
    xml = ET.tostring(doc, xml_declaration=True, encoding="UTF-8")
    tmp = path + ".tmp"
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as outz:
        for name in zin.namelist():
            outz.writestr(name, xml if name == "word/document.xml" else zin.read(name))
    zin.close()
    shutil.move(tmp, path)
    print(f"figures wrapped in a two-row table: {made}")


def main() -> None:
    text = build()
    out_dir = os.path.join(ROOT, "report")
    md = os.path.join(out_dir, "parascape_manuscript.md")
    open(md, "w").write(text)
    print("wrote", md, f"({len(text.splitlines())} lines)")

    docx = os.path.join(out_dir, "parascape_manuscript.docx")
    try:
        subprocess.run(["pandoc", md, "-o", docx, f"--resource-path={ROOT}"], check=True)
    except (subprocess.CalledProcessError, FileNotFoundError) as error:
        print("pandoc did not run:", error)
        return
    figures_into_tables(docx)
    print("wrote", docx)


if __name__ == "__main__":
    main()
