#!/usr/bin/env python3
"""Build the contrast ladder: every comparison a study supports, each against its own baseline.

The keyword rule this replaces asked one question of each sample group -- "does the label look like a
control?" -- and so it read "WT" in `T. gondii VEG WT, Switch Conditions` and called the alkaline-stressed
arm a baseline. The mistake was reading a level on the *strain* axis as an answer about the *treatment*
axis, which is why a patched keyword list would not have helped.

The rule here is different, and it is the one Einar stated: a study with several levels of baseline gets
several levels of contrast. Untreated knockouts are compared to untreated wild type; treated samples are
compared to the untreated samples of their own genotype. So:

  1. Parse each group label into factors -- treatment, genotype, strain, timepoint, tissue.
  2. Keep the axes that actually vary inside that study.
  3. Pair groups that differ on EXACTLY ONE axis. That pair is a contrast on that axis.
  4. Orient it: the group sitting at the axis reference level (untreated, wild type, earliest time) is the
     baseline.
  5. Depth = how many of the OTHER axes are off their own reference level. Depth 0 is the clean effect on an
     otherwise-baseline background; depth 1 is the same effect inside a mutant or a stressed background.

Nothing is averaged across depths. A depth-1 treatment effect inside a knockout is a different measurement
from the depth-0 effect in wild type, and the ladder keeps them apart so the map can say which it holds.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# Axis names as the depositors write them, mapped onto the axis we reason about. `condition` is folded into
# treatment because GEO submitters use the two interchangeably for the same thing.
AXIS_ALIAS = {
    'treatment': 'treatment', 'condition': 'treatment', 'agent': 'treatment', 'drug': 'treatment',
    'infection': 'treatment', 'stress': 'treatment', 'stimulus': 'treatment',
    'genotype': 'genotype', 'genotype/variation': 'genotype', 'cell line': 'genotype', 'line': 'genotype',
    'strain': 'strain', 'isolate': 'strain',
    'timepoint': 'time', 'time': 'time', 'time point': 'time', 'age': 'time', 'hours post infection': 'time',
    'tissue': 'tissue', 'cell type': 'tissue', 'source': 'tissue',
    'developmental stage': 'stage', 'stage': 'stage', 'life stage': 'stage',
    'sex': 'sex', 'replicate': None, 'biological replicate': None, 'library': None,
}

# A level is the reference on the treatment axis when it says so. Ordered longest-first on purpose: "no atc"
# must win over "atc".
TREATMENT_REFERENCE = (
    'control media', 'control condition', 'standard condition', 'normal condition', 'complete media',
    'no treatment', 'not treated', 'non-treated', 'nontreated', 'untreated', 'unstimulated', 'uninfected',
    'no atc', 'minus atc', '-atc', 'without atc', 'no iaa', '-iaa', 'no rapamycin', 'no shield',
    'vehicle', 'mock', 'dmso', 'control', 'ctrl', 'none', 'normal', 'unstressed', 'naive', 'baseline',
    'wild type', 'wild-type', 'parental',
)
# ...and on the genotype axis.
GENOTYPE_REFERENCE = (
    'wild type', 'wild-type', 'wildtype', 'wt', 'parental', 'parent', 'empty vector', 'scrambled',
    'non-targeting', 'nontargeting', 'control', 'cas9', 'luciferase', 'untagged', 'tir1',
)
# Reference levels we must never infer. A bare strain name is a choice of organism, not a control, and
# picking one silently is how a strain difference gets reported as a treatment effect.
NO_AUTO_REFERENCE = ('strain', 'tissue', 'stage', 'sex')

TIME_RE = re.compile(r'(-?\d+(?:\.\d+)?)\s*(h|hr|hrs|hour|hours|d|day|days|min|minutes|m|wk|week|weeks)\b')
# `TgROCY1 KO`, `BFD1 KO`, `cKD`, `BFD1-OE`. The symbol is whatever sits in front of the operation word.
EDIT_RE = re.compile(
    r'\b([A-Za-z][A-Za-z0-9_.-]{1,24})[\s_-]*'
    r'(ko|k\.o\.|knockout|knock-out|kd|knockdown|knock-down|ckd|cko|ikd|ki|oe|overexpress\w*|'
    r'depleted|deletion|null|mutant|tagged|ha|myc|aid)\b', re.I)
WT_RE = re.compile(r'\b(wt|wild[\s-]?type|parental)\b', re.I)


def _clean(s: str) -> str:
    return re.sub(r'\s+', ' ', s.replace('_', ' ')).strip(' ,;:.').lower()


def parse_factors(label: str) -> tuple[dict, str]:
    """Split one group label into {axis: level}. Returns the factors and the free-text tail.

    GEO writes the structured part as `strain: VEG; treatment: high pH, low C02, low serum conditions`, so
    the split is on `;` and never on `,` -- a treatment level legitimately contains commas.
    """
    head, _, tail = label.partition('|')
    factors: dict[str, str] = {}
    for part in head.split(';'):
        if ':' not in part:
            continue
        key, _, value = part.partition(':')
        axis = AXIS_ALIAS.get(_clean(key))
        if axis is None:
            continue
        value = _clean(value)
        if value and value not in ('na', 'n/a', 'none given', 'not specified', 'unknown'):
            # Keep the first mention of an axis; GEO sometimes repeats one.
            factors.setdefault(axis, value)
    return factors, _clean(tail)


def genotype_from_tail(tail: str) -> str | None:
    """Recover a genotype that the depositor only wrote into the free text.

    GSE218921 is the case that matters: its structured head is identical for WT, ROCY1 KO and BFD1 KO --
    all three read `strain: VEG; treatment: ...; tissue: ...`. The genotype exists only in the tail, so
    without this the three strains collapse into one group and no genotype axis is seen at all.
    """
    edit = EDIT_RE.search(tail)
    if edit:
        symbol, op = edit.group(1), edit.group(2).lower()
        if _clean(symbol) in ('t', 'tg', 'gondii', 'the', 'a', 'in', 'of', 'and'):
            return None
        return f'{symbol.lower()} {op}'
    if WT_RE.search(tail):
        return 'wild type'
    return None


# Tokens that mark a replicate rather than a condition. They may differ between two groups without making
# the comparison unsafe; anything else that differs is an axis nobody parsed.
REPLICATE_RE = re.compile(r'^(array|rep|replicate|r|run|batch|exp|biorep|br|sample|s|lane|set)[\s_-]?\d*$')
NOISE_TOKENS = frozenset(('t', 'tg', 'gondii', 'toxoplasma', 'in', 'vitro', 'vivo', 'the', 'of', 'and',
                          'conditions', 'condition', 'cells', 'cell', 'parasites', 'parasite', 'strain'))


def residual_tokens(label: str, factors: dict) -> frozenset:
    """Everything in a label that no parsed factor explains.

    This is the safety rule, and it is what makes the ladder trustworthy rather than merely tidy. Two groups
    may only be compared when the text left over after removing every parsed factor is the same. If it is
    not, the study varies on an axis that was never parsed, and pairing across it compares two things at
    once -- `WT_stress_24h` against `KO_unstress_24h` in GSE132238, where the head declares only genotype
    and timepoint while `stress`/`unstress` hides in the free text. Refusing such a pair costs a contrast;
    forming it silently reports a stress response as a genotype effect.
    """
    head, _, _ = label.partition('|')
    text = _clean(label.replace('|', ' ').replace('(', ' ').replace(')', ' '))
    # Strip only the axis names this label actually declares as keys. Stripping every known alias would
    # delete `stress` from `WT_stress_24h`, where it is a treatment LEVEL hiding in free text and the very
    # thing that has to be noticed.
    for part in head.split(';'):
        key = _clean(part.partition(':')[0])
        if key:
            text = re.sub(r'\b' + re.escape(key) + r'\b', ' ', text)
    for axis, level in factors.items():
        text = text.replace(level, ' ')
        # A genotype recovered from free text reads `bfd1 ko` there but `BFD1 KO` or `{delta}BFD1` in the
        # label, so strip its words individually too.
        for word in level.split():
            if len(word) > 1:
                text = re.sub(r'\b' + re.escape(word) + r'\b', ' ', text)
    tokens = set()
    for raw in re.split(r'[^a-z0-9+]+', text):
        token = raw.strip('+')
        if not token or token in NOISE_TOKENS or REPLICATE_RE.match(token) or token.isdigit():
            continue
        tokens.add(token)
    return frozenset(tokens)


def time_from_tail(tail: str) -> str | None:
    """Recover a timepoint the depositor only wrote into the free text (GSE140197: `4hr`, `8hr`, `12hr`)."""
    m = TIME_RE.search(tail)
    return m.group(0) if m else None


# Tokens that look like a plate, array or batch identifier rather than a condition: `1a`, `2b`, a lone
# letter. They are not treated as replicates automatically -- guessing wrong would merge two real
# conditions -- but refusals caused only by these are marked so the review list can be read at a glance.
PROBABLE_ID_RE = re.compile(r'^(\d+[a-z]|[a-z]|[a-z]\d+)$')


def reference_rank(axis: str, level: str) -> tuple[int, str]:
    """How strongly `level` reads as the reference on `axis`. Lower rank wins; 99 means no claim.

    Returns the rank and the evidence, so every orientation in the output can be traced to the word that
    decided it rather than to the order rows happened to arrive in.
    """
    if axis in NO_AUTO_REFERENCE:
        return 99, ''
    if axis == 'time':
        m = TIME_RE.search(level)
        if m:
            value = float(m.group(1))
            unit = m.group(2).lower()
            mult = {'min': 1 / 60, 'm': 1 / 60, 'h': 1, 'hr': 1, 'hrs': 1, 'hour': 1, 'hours': 1,
                    'd': 24, 'day': 24, 'days': 24, 'wk': 168, 'week': 168, 'weeks': 168}[unit]
            # Earliest time is the reference. Offset so an absolute ordering survives the sort.
            return int(round(value * mult)) + 1, f'earliest time ({m.group(0)})'
        return 99, ''
    vocabulary = GENOTYPE_REFERENCE if axis == 'genotype' else TREATMENT_REFERENCE
    for i, word in enumerate(vocabulary):
        if re.search(r'(?<![a-z])' + re.escape(word) + r'(?![a-z])', level):
            # An edited genotype that also says "control" is not the reference; the edit wins.
            if axis == 'genotype' and EDIT_RE.search(level) and word in ('control', 'cas9', 'tir1'):
                continue
            return i, f'"{word}"'
    return 99, ''


def build(rows: list[dict], group_key: str, study_key: str) -> tuple[list[dict], list[dict]]:
    groups: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in rows:
        groups[(row[study_key], row[group_key])].append(row)

    parsed: dict[str, dict[str, dict]] = defaultdict(dict)
    for (study, label), members in groups.items():
        factors, tail = parse_factors(label)
        if 'genotype' not in factors:
            found = genotype_from_tail(tail) or genotype_from_tail(_clean(label.split('|')[0]))
            if found:
                factors['genotype'] = found
        if 'time' not in factors:
            found = time_from_tail(tail)
            if found:
                factors['time'] = found
        parsed[study][label] = {'factors': factors, 'tail': tail, 'n': len(members),
                                'residual': residual_tokens(label, factors)}

    # An axis some groups declare and others omit is the subtlest version of the original mistake: the
    # comparison simply skips it. GSE136123 has two different "untreated" baselines -- `Pru MORC KD_UT` and
    # `Pru_US` -- and because only the first parses a genotype, both read as the same baseline and the
    # knockdown is compared against the parental line as if it were the same strain. An explicit `unstated`
    # level makes them differ, and since `unstated` matches no reference vocabulary the resulting genotype
    # contrast goes to review instead of being oriented on a guess.
    for study, labels in parsed.items():
        declared = Counter(a for info in labels.values() for a in info['factors'])
        for axis, count in declared.items():
            if count == len(labels):
                continue
            for info in labels.values():
                info['factors'].setdefault(axis, 'unstated')
        for label, info in labels.items():
            info['residual'] = residual_tokens(label, info['factors'])

    # Two groups with the same factors and the same leftover text are replicates of one condition, not two
    # conditions. GSE132238 deposits `WT_stress_24h_array1` and `..._array2` as separate sample groups;
    # left apart they emit the same contrast twice and each carries half the replicates it really has.
    merged: dict[str, dict[str, dict]] = {}
    for study, labels in parsed.items():
        buckets: dict[tuple, list] = defaultdict(list)
        for label, info in labels.items():
            key = (tuple(sorted(info['factors'].items())), info['residual'])
            buckets[key].append((label, info))
        out: dict[str, dict] = {}
        for members_ in buckets.values():
            label, info = sorted(members_)[0]
            out[label] = {**info, 'n': sum(i['n'] for _, i in members_),
                          'n_groups_merged': len(members_)}
        merged[study] = out
    parsed = merged

    contrasts, unresolved, refused = [], [], []
    for study, labels in sorted(parsed.items()):
        # Only axes that take more than one level inside this study can carry a contrast.
        levels: dict[str, set] = defaultdict(set)
        for info in labels.values():
            for axis, level in info['factors'].items():
                levels[axis].add(level)
        varying = {a for a, v in levels.items() if len(v) > 1}
        if not varying:
            continue

        # Attribute each leftover token to an axis, by asking the study rather than a word list: a token
        # belongs to axis A when its presence is constant inside every level of A. `Switch Conditions` in
        # the GSE218921 tail is then recognised as a restatement of `treatment: high pH`, while `stress` in
        # GSE132238 is constant inside no axis and so stands revealed as an axis nobody parsed. Only those
        # unattributable tokens can refuse a pair.
        # The test runs inside strata, because tissue, strain and stage separate what are really different
        # experiments sharing one accession. GSE218921 is why: its tachyzoite arm words the stress as
        # `Switch Conditions` and its sporozoite arm as `bradyzoite induction conditions`, so study-wide
        # `switch` looks unexplained, and the alkaline contrast -- the one this whole item exists to
        # recover -- would be refused a second time for a second wrong reason.
        def stratum(info: dict) -> tuple:
            return tuple(sorted((a, l) for a, l in info['factors'].items() if a in NO_AUTO_REFERENCE))

        all_tokens = set().union(*(i['residual'] for i in labels.values())) if labels else set()
        unattributed = set()
        for token in all_tokens:
            carriers = {l for l, i in labels.items() if token in i['residual']}
            if len(carriers) == len(labels):
                continue  # present everywhere, explains nothing
            for axis in varying:
                cells = defaultdict(set)
                for label, info in labels.items():
                    if axis in info['factors']:
                        cells[(stratum(info), info['factors'][axis])].add(label in carriers)
                if cells and all(len(v) == 1 for v in cells.values()):
                    break  # determined by this axis within every stratum
            else:
                unattributed.add(token)
        if unattributed:
            unresolved.append({
                'study': study, 'axis': 'UNPARSED (free text)', 'n_levels': len(unattributed),
                'levels': ' | '.join(sorted(unattributed)[:8]),
                'why': ('these words look like plate or batch identifiers; confirm they are not conditions'
                        if all(PROBABLE_ID_RE.match(t) for t in unattributed)
                        else 'these words vary between groups and no declared axis explains them'),
            })

        # The reference level of each varying axis, used to measure depth.
        axis_reference = {}
        for axis in varying:
            ranked = sorted((reference_rank(axis, l)[0], l) for l in levels[axis])
            axis_reference[axis] = ranked[0][1] if ranked[0][0] < 99 else None
            if axis_reference[axis] is None:
                unresolved.append({
                    'study': study, 'axis': axis, 'n_levels': len(levels[axis]),
                    'levels': ' | '.join(sorted(levels[axis])[:6]),
                    'why': 'axis never carries an automatic reference' if axis in NO_AUTO_REFERENCE
                           else 'no level matches the reference vocabulary',
                })

        for (la, ia), (lb, ib) in combinations(sorted(labels.items()), 2):
            fa, fb = ia['factors'], ib['factors']
            # Compare only on axes both groups declare, so a group missing one annotation is not treated as
            # differing on it.
            shared = set(fa) & set(fb)
            differ = [a for a in shared if fa[a] != fb[a]]
            if len(differ) != 1:
                continue
            axis = differ[0]
            if axis not in varying:
                continue
            hidden = (ia['residual'] ^ ib['residual']) & unattributed
            if hidden:
                # An axis nobody parsed separates these two groups. Refuse rather than guess.
                refused.append({
                    'study': study, 'declared_axis': axis,
                    'unparsed_difference': ' / '.join(sorted(hidden)),
                    'looks_like': 'plate or batch identifier' if all(PROBABLE_ID_RE.match(t) for t in hidden)
                                  else 'a real condition axis',
                    'group_a': la, 'group_b': lb,
                })
                continue
            ra, ea = reference_rank(axis, fa[axis])
            rb, eb = reference_rank(axis, fb[axis])
            if ra == rb == 99:
                continue  # neither side claims to be the baseline; left for review
            if ra <= rb:
                base, treat, base_label, treat_label, evidence = fa, fb, la, lb, ea
            else:
                base, treat, base_label, treat_label, evidence = fb, fa, lb, la, eb

            # Depth: the other axes, held equal in this pair, that sit off their own reference level.
            off = sorted(a for a in shared
                         if a != axis and axis_reference.get(a) and base[a] != axis_reference[a])
            # The background has to name every other varying axis this pair holds fixed, including the ones
            # that carry no reference and so add no depth. Without strain and tissue in here, the two
            # depth-0 alkaline contrasts of GSE218921 -- one in replicating tachyzoites, one in freshly
            # excysted sporozoites -- print as identical rows and would be merged as one measurement.
            held = sorted(a for a in shared
                          if a != axis and a in varying and a not in off)
            # Unparsed words the pair holds in common still describe the background it was measured on. In
            # GSE132238 the genotype effect under stress and the genotype effect unstressed are two
            # different measurements, and `stress` is the only thing that says which is which -- without it
            # here they print as one row and would be averaged together.
            shared_unparsed = sorted((ia['residual'] & ib['residual']) & unattributed)
            background = '; '.join([f'{a}={base[a]}' for a in off]
                                   + [f'{a}={base[a]}' for a in held]
                                   + ([f'unparsed={",".join(shared_unparsed)}'] if shared_unparsed else []))
            contrasts.append({
                'study': study, 'axis': axis,
                'baseline_level': base[axis], 'treated_level': treat[axis],
                'depth': len(off),
                'background': background or 'reference background',
                'oriented_by': evidence,
                'baseline_group': base_label, 'treated_group': treat_label,
                'n_baseline': (ia if base is fa else ib)['n'],
                'n_treated': (ib if base is fa else ia)['n'],
            })
    return contrasts, unresolved, refused


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--samples', default='results/rna_sample_quality_2026_11_04/sample_quality.tsv')
    ap.add_argument('--out', default='results/contrast_ladder_2026_11_05')
    ap.add_argument('--group-column', default='condition_group')
    ap.add_argument('--study-column', default='study')
    args = ap.parse_args()

    src = (REPO / args.samples) if not Path(args.samples).is_absolute() else Path(args.samples)
    with src.open() as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    contrasts, unresolved, refused = build(rows, args.group_column, args.study_column)

    out = (REPO / args.out) if not Path(args.out).is_absolute() else Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    fields = ['study', 'axis', 'baseline_level', 'treated_level', 'depth', 'background', 'oriented_by',
              'n_baseline', 'n_treated', 'baseline_group', 'treated_group']
    contrasts.sort(key=lambda c: (c['study'], c['axis'], c['depth'], c['treated_level']))
    with (out / 'contrasts.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(contrasts)
    with (out / 'reference_needs_review.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['study', 'axis', 'n_levels', 'levels', 'why'], delimiter='\t')
        w.writeheader()
        w.writerows(unresolved)
    refused.sort(key=lambda r: (r['study'], r['unparsed_difference']))
    with (out / 'refused_unparsed_axis.tsv').open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['study', 'declared_axis', 'unparsed_difference', 'looks_like',
                                           'group_a', 'group_b'], delimiter='\t')
        w.writeheader()
        w.writerows(refused)

    by_axis = Counter(c['axis'] for c in contrasts)
    by_depth = Counter(c['depth'] for c in contrasts)
    signature = Counter((c['study'], c['axis'], c['baseline_level'], c['treated_level'], c['background'])
                        for c in contrasts)
    summary = {
        'studies_in': len({r[args.study_column] for r in rows}),
        'studies_with_contrasts': len({c['study'] for c in contrasts}),
        'contrasts': len(contrasts),
        'by_axis': dict(by_axis.most_common()),
        'by_depth': {str(k): v for k, v in sorted(by_depth.items())},
        'axes_needing_a_named_reference': len(unresolved),
        'pairs_refused_for_an_unparsed_axis': len(refused),
        'studies_with_a_refusal': len({r['study'] for r in refused}),
        'refusals_on_a_real_axis': sum(1 for r in refused if r['looks_like'] == 'a real condition axis'),
        'refusals_on_a_plate_label': sum(1 for r in refused
                                         if r['looks_like'] == 'plate or batch identifier'),
        # Must stay 0. A repeated signature means two different measurements print as one row.
        'indistinguishable_rows': sum(v for v in signature.values() if v > 1),
    }
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
