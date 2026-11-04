"""Turn Einar's review of the 140 candidate claims into the per-gene host-facing truth set (instruction 72).

Einar reviewed every claim on the review page (artifact TbAj2n9cnBPxZ7Leb5duCh) on 2026-11-04: 107 rows
host-facing, 33 not. Several genes carry more than one published claim, so the rows are collapsed per gene.

Rules he set the same day:
* **Strongest evidence wins** -- a direct experimental claim outranks a proximity hit, whichever way it points.
* GRA60, GRA44 and GRA47, where he marked the proximity row yes and the more specific row no, are **not**
  host-facing: they follow the 22 other proximity-only dense granule proteins he rejected. Recorded as explicit
  overrides rather than derived, because for GRA44 the rule alone would say otherwise.

The decision itself is expert curation on published evidence, so a gene decided by a direct experimental claim
keeps `observed-direct` (0.95) and one decided by a proximity or inferred claim becomes `observed-curated`
(0.90): Einar's judgement is worth more than the raw screen call it replaces, and less than a published
localisation.
"""
import argparse
import glob
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/host_facing_truth_2026_11_04'
ARTIFACT = 'https://claude.ai/artifact/TbAj2n9cnBPxZ7Leb5duCh'
#: Lower rank wins when a gene's claims disagree.
STRENGTH_RANK = {'direct': 0, 'proximity': 1, 'inferred': 2}
#: Einar 2026-11-04, explicit: proximity yes but a more specific no -> not host-facing.
OVERRIDE = {'TGME49_204270': 'no', 'TGME49_228170': 'no', 'TGME49_254000': 'no'}
TIER = {'direct': ('observed-direct', 0.95), 'proximity': ('observed-curated', 0.90),
        'inferred': ('observed-curated', 0.90)}


def decisions(folder: Path) -> pd.DataFrame:
    rows = []
    for path in sorted(glob.glob(str(folder / '*.json'))):
        raw = json.loads(Path(path).read_text())
        value = raw.get('data', raw)
        rows.append(dict(key=raw.get('id', Path(path).stem),
                         **{k: value.get(k) for k in ('decision', 'note', 'gene', 'protein', 'place',
                                                      'pmid', 'strength')}))
    frame = pd.DataFrame(rows)
    frame['gene'] = frame.gene.fillna('')
    return frame


def per_gene(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame[frame.gene != ''].copy()
    frame['rank'] = frame.strength.map(STRENGTH_RANK).fillna(9)
    out = []
    for gene, rows in frame.groupby('gene'):
        best = rows.sort_values(['rank', 'decision']).iloc[0]
        decision, basis = best.decision, 'strongest evidence'
        if gene in OVERRIDE:
            decision, basis = OVERRIDE[gene], 'Einar explicit override'
        tier, probability = TIER.get(str(best.strength), ('observed-curated', 0.90))
        out.append(dict(gene_id=gene, protein=best.protein, decision=decision, deciding_strength=best.strength,
                        deciding_place=best.place, deciding_pmid=best.pmid, tier=tier, probability=probability,
                        claims=len(rows), disagreed=rows.decision.nunique() > 1, basis=basis,
                        notes='; '.join(sorted({str(n) for n in rows.note.dropna() if str(n).strip()}))))
    return pd.DataFrame(out).sort_values(['decision', 'gene_id'])


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--decisions', type=Path, required=True, help='folder of decision JSON files')
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    raw = decisions(args.decisions)
    raw.to_csv(OUT / 'reviewed_claims.tsv', sep='\t', index=False)
    truth = per_gene(raw)
    truth.to_csv(OUT / 'host_facing_truth.tsv', sep='\t', index=False)
    evidence = pd.DataFrame([dict(gene_id=r.gene_id, trait='host_facing',
                                  value='yes' if r.decision == 'yes' else 'no', probability=r.probability,
                                  source_tier=r.tier,
                                  reference=f'PMID:{r.deciding_pmid}' if str(r.deciding_pmid) else '',
                                  dataset=f'Einar review 2026-11-04 of {r.protein} ({r.deciding_place}, '
                                          f'{r.deciding_strength}; {r.basis}); {ARTIFACT}')
                             for r in truth.itertuples()])
    evidence.to_parquet(OUT / 'evidence.parquet', index=False)
    counts = truth.decision.value_counts().to_dict()
    (OUT / 'README.md').write_text(f"""# Host-facing truth set, reviewed by Einar (2026-11-04)

Built by `scripts/build_host_facing_truth.py` from his decisions on every one of the 140 candidate claims
({ARTIFACT}). {len(raw)} claims over {truth.gene_id.nunique()} genes.

* **{counts.get('yes', 0)} genes host-facing**, **{counts.get('no', 0)} not**.
* Rule: strongest evidence wins (direct > proximity > inferred). {int(truth.disagreed.sum())} genes had
  disagreeing claims; GRA60, GRA44 and GRA47 are explicit overrides to "not host-facing".
* Probability: a gene decided by a direct experimental claim is `observed-direct` 0.95, one decided by a
  proximity or inferred claim is `observed-curated` 0.90.
* What makes this set valuable: {counts.get('no', 0)} negatives are **vacuole-resident** proteins, the hard
  kind. The earlier map had 9 such negatives and several hundred markers of compartments inside the parasite,
  which answer an easier question. Einar rejected 22 dense granule proteins that the PVM proximity screen had
  called positive, which is the distinction the first model failed to learn.

## Files
* `host_facing_truth.tsv` -- one row per gene: decision, the claim that decided it, tier, probability.
* `reviewed_claims.tsv` -- all {len(raw)} reviewed claims as recorded.
* `evidence.parquet` -- evidence rows for the information map.
""")
    print(json.dumps(dict(claims=len(raw), genes=int(truth.gene_id.nunique()), **counts,
                          disagreed=int(truth.disagreed.sum()), out=str(OUT))))


if __name__ == '__main__':
    main()
