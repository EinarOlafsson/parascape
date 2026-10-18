#!/usr/bin/env python3
"""Stream public RNA-seq into ParaScape: download each run, quantify it every documented way, keep the
results on the NAS and delete the reads (instruction 69).

    python scripts/omics_quant.py references            # fetch + index the three references
    python scripts/omics_quant.py runs GSE60305 ...      # resolve studies to runs (ENA), write the plan
    python scripts/omics_quant.py quantify               # work through the plan, resumably

Conventions kept side by side, so a user can choose:
  * quantifier: salmon (selective alignment, genome decoys, GC + sequence bias correction) and kallisto
    (pseudoalignment);
  * reference: ME49 for every study (shared gene identifiers), and for strain studies also the strain's
    own genome (GT1 = Type I, VEG = Type III) -- mapping a Type I read to ME49 penalises SNP-dense genes;
  * per run: estimated counts, TPM, effective lengths -- the normalisations (CPM, TMM, median-of-ratios,
    upper quartile, within-study control-relative) are computed from these at contrast time.

Reads come from ENA as gzipped FASTQ (no fasterq-dump expansion) and are deleted after quantification;
they remain re-downloadable. Every output records the tool version, reference checksum and run metadata.
Each step runs with the binaries of the separate `parascape-quant` conda environment.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
import urllib.parse
import urllib.request

NAS = os.environ.get("PARASCAPE_OMICS", "/nas_mnt/carruthers/Einar/parascape_omics")
SCRATCH = os.environ.get("PARASCAPE_SCRATCH", "/media/carruthers/mnt3/claude/parascape_omics_scratch")
BIN = os.path.expanduser("~/anaconda3/envs/parascape-quant/bin")
THREADS = int(os.environ.get("PARASCAPE_THREADS", "4"))

ENSEMBL = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-63/fasta"
#: Resolved from the Ensembl Genomes release-63 FTP listing on 2026-10-17 (instruction 69).
REFERENCES = {
    "ME49": (f"{ENSEMBL}/toxoplasma_gondii/cdna/Toxoplasma_gondii.GCA000006565v2.cdna.all.fa.gz",
             f"{ENSEMBL}/toxoplasma_gondii/dna/Toxoplasma_gondii.GCA000006565v2.dna.toplevel.fa.gz"),
    "GT1": (f"{ENSEMBL}/protists_alveolata1_collection/toxoplasma_gondii_gt1_gca_000149715/cdna/"
            "Toxoplasma_gondii_gt1_gca_000149715.TGGT1.cdna.all.fa.gz",
            f"{ENSEMBL}/protists_alveolata1_collection/toxoplasma_gondii_gt1_gca_000149715/dna/"
            "Toxoplasma_gondii_gt1_gca_000149715.TGGT1.dna.toplevel.fa.gz"),
    "VEG": (f"{ENSEMBL}/protists_alveolata1_collection/toxoplasma_gondii_veg_gca_000150015/cdna/"
            "Toxoplasma_gondii_veg_gca_000150015.TGVEG.cdna.all.fa.gz",
            f"{ENSEMBL}/protists_alveolata1_collection/toxoplasma_gondii_veg_gca_000150015/dna/"
            "Toxoplasma_gondii_veg_gca_000150015.TGVEG.dna.toplevel.fa.gz"),
}

#: Host transcriptomes (Ensembl release 116, resolved from its FTP listing on 2026-10-17), for the
#: competitive "ME49+host" reference. In dual RNA-seq of intracellular parasites most reads are the host
#: cell's (the first Trx1 run mapped 33% to Toxoplasma); against a parasite-only index, host reads from
#: conserved genes -- actin, tubulin, ribosomal proteins -- can be counted as parasite. With the host's
#: transcripts present they compete and lose. rRNA (ncrna) is included for the same reason.
ENSEMBL_MAIN = "https://ftp.ensembl.org/pub/release-116/fasta"
HOST_TRANSCRIPTS = (
    f"{ENSEMBL_MAIN}/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
    f"{ENSEMBL_MAIN}/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz",
    f"{ENSEMBL_MAIN}/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz",
    f"{ENSEMBL_MAIN}/mus_musculus/ncrna/Mus_musculus.GRCm39.ncrna.fa.gz",
)
#: Gene identifiers given to host transcripts, so a matrix can drop them (and report the host share).
HOST_PREFIX = "HOST:"

#: Free space a download may leave on scratch, at least.
SCRATCH_MARGIN_GB = 25


def log(msg: str) -> None:
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    print(line, flush=True)
    os.makedirs(os.path.join(NAS, "logs"), exist_ok=True)
    with open(os.path.join(NAS, "logs", "omics_quant.log"), "a") as fh:
        fh.write(line + "\n")


def fetch(url: str, dest: str, tries: int = 5) -> str:
    """Download to dest (atomically), with retries; returns the file's sha256."""
    tmp = dest + ".part"
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as fh:
                shutil.copyfileobj(r, fh, length=1 << 22)
            os.replace(tmp, dest)
            return sha256(dest)
        except Exception as e:                       # network: retry with backoff
            log(f"  fetch failed ({i + 1}/{tries}) {url}: {e}")
            time.sleep(10 * (i + 1))
    raise RuntimeError(f"could not download {url}")


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def run(cmd: list, **kw) -> None:
    log("  $ " + " ".join(map(str, cmd)))
    subprocess.run(cmd, check=True, **kw)


def version(tool: str) -> str:
    out = subprocess.run([os.path.join(BIN, tool), "--version" if tool == "salmon" else "version"],
                         capture_output=True, text=True)
    return (out.stdout or out.stderr).strip()


# --------------------------------------------------------------------------- references
def references() -> None:
    for name, (cdna_url, dna_url) in REFERENCES.items():
        d = os.path.join(NAS, "references", name)
        os.makedirs(d, exist_ok=True)
        meta_path = os.path.join(d, "reference.json")
        if os.path.exists(os.path.join(d, "salmon_index", "versionInfo.json")) and \
                os.path.exists(os.path.join(d, "kallisto.idx")):
            log(f"{name}: indexes present")
            continue
        cdna = os.path.join(d, os.path.basename(cdna_url))
        dna = os.path.join(d, os.path.basename(dna_url))
        meta = {"name": name, "cdna_url": cdna_url, "dna_url": dna_url}
        meta["cdna_sha256"] = fetch(cdna_url, cdna) if not os.path.exists(cdna) else sha256(cdna)
        meta["dna_sha256"] = fetch(dna_url, dna) if not os.path.exists(dna) else sha256(dna)
        # Salmon with the genome as decoys (Srivastava et al. 2020): reads from unannotated sequence are
        # not forced onto a transcript. The gentrome is transcripts first, then genome.
        decoys = os.path.join(d, "decoys.txt")
        with gzip.open(dna, "rt") as fh, open(decoys, "w") as out:
            for line in fh:
                if line.startswith(">"):
                    out.write(line[1:].split()[0] + "\n")
        gentrome = os.path.join(SCRATCH, f"{name}_gentrome.fa.gz")
        os.makedirs(SCRATCH, exist_ok=True)
        with open(gentrome, "wb") as out:
            for part in (cdna, dna):
                with open(part, "rb") as fh:
                    shutil.copyfileobj(fh, out)
        run([os.path.join(BIN, "salmon"), "index", "-t", gentrome, "-d", decoys, "-p", str(THREADS),
             "-i", os.path.join(d, "salmon_index"), "-k", "31"])
        os.remove(gentrome)
        run([os.path.join(BIN, "kallisto"), "index", "-i", os.path.join(d, "kallisto.idx"), cdna])
        # Transcript -> gene, from the cDNA headers ("gene:XXXX").
        t2g = os.path.join(d, "tx2gene.tsv")
        with gzip.open(cdna, "rt") as fh, open(t2g, "w") as out:
            for line in fh:
                if line.startswith(">"):
                    tx = line[1:].split()[0]
                    gene = next((f.split(":", 1)[1] for f in line.split() if f.startswith("gene:")), tx)
                    out.write(f"{tx}\t{gene}\n")
        meta["salmon"] = version("salmon")
        meta["kallisto"] = version("kallisto")
        meta["built"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        with open(meta_path, "w") as fh:
            json.dump(meta, fh, indent=2)
        log(f"{name}: indexed")



def host_reference() -> None:
    """ME49 + human + mouse transcripts (with ncRNA), ME49 genome as decoys: the competitive reference."""
    d = os.path.join(NAS, "references", "ME49_host")
    if os.path.exists(os.path.join(d, "salmon_index", "versionInfo.json")) and \
            os.path.exists(os.path.join(d, "kallisto.idx")):
        log("ME49_host: indexes present")
        return
    os.makedirs(d, exist_ok=True)
    me49 = os.path.join(NAS, "references", "ME49")
    cdna_url, dna_url = REFERENCES["ME49"]
    parts = [os.path.join(me49, os.path.basename(cdna_url))]
    meta = {"name": "ME49_host", "parasite": cdna_url, "decoys": dna_url, "host": list(HOST_TRANSCRIPTS),
            "sha256": {}}
    for url in HOST_TRANSCRIPTS:
        dest = os.path.join(d, os.path.basename(url))
        meta["sha256"][os.path.basename(url)] = fetch(url, dest) if not os.path.exists(dest) else sha256(dest)
        parts.append(dest)
    # tx2gene: parasite genes as in ME49; every host transcript to HOST:<gene>.
    t2g = os.path.join(d, "tx2gene.tsv")
    with open(t2g, "w") as out:
        for i, part in enumerate(parts):
            with gzip.open(part, "rt") as fh:
                for line in fh:
                    if line.startswith(">"):
                        tx = line[1:].split()[0]
                        gene = next((f.split(":", 1)[1] for f in line.split() if f.startswith("gene:")), tx)
                        out.write(f"{tx}\t{gene if i == 0 else HOST_PREFIX + gene}\n")
    targets = os.path.join(SCRATCH, "ME49_host_targets.fa.gz")
    os.makedirs(SCRATCH, exist_ok=True)
    with open(targets, "wb") as out:
        for part in parts:
            with open(part, "rb") as fh:
                shutil.copyfileobj(fh, out)
    gentrome = os.path.join(SCRATCH, "ME49_host_gentrome.fa.gz")
    with open(gentrome, "wb") as out:
        for part in (targets, os.path.join(me49, os.path.basename(dna_url))):
            with open(part, "rb") as fh:
                shutil.copyfileobj(fh, out)
    run([os.path.join(BIN, "salmon"), "index", "-t", gentrome, "-d", os.path.join(me49, "decoys.txt"),
         "-p", str(THREADS), "-i", os.path.join(d, "salmon_index"), "-k", "31"])
    run([os.path.join(BIN, "kallisto"), "index", "-i", os.path.join(d, "kallisto.idx"), "-t",
         str(THREADS), targets])
    os.remove(gentrome)
    os.remove(targets)
    meta.update({"salmon": version("salmon"), "kallisto": version("kallisto"),
                 "built": time.strftime("%Y-%m-%dT%H:%M:%S")})
    with open(os.path.join(d, "reference.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    log("ME49_host: indexed")

# --------------------------------------------------------------------------- the plan
def ena_runs(accession: str) -> list:
    """Run accessions with their FASTQ URLs, sizes and checksums, from the ENA file report."""
    fields = ("run_accession,study_accession,secondary_study_accession,sample_accession,"
              "sample_alias,sample_title,library_layout,library_strategy,scientific_name,"
              "fastq_ftp,fastq_md5,fastq_bytes,read_count,instrument_model")
    url = ("https://www.ebi.ac.uk/ena/portal/api/filereport?" +
           urllib.parse.urlencode({"accession": accession, "result": "read_run", "fields": fields,
                                   "format": "tsv"}))
    with urllib.request.urlopen(url, timeout=120) as r:
        lines = r.read().decode().strip().split("\n")
    head = lines[0].split("\t")
    return [dict(zip(head, ln.split("\t"))) for ln in lines[1:] if ln.strip()]


def geo_project(gse: str) -> str | None:
    """A GEO series' BioProject, through E-utilities (GEO accessions are not ENA accessions)."""
    E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    q = f"{E}/esearch.fcgi?" + urllib.parse.urlencode({"db": "gds", "term": f"{gse}[ACCN] AND gse[ETYP]",
                                                        "retmode": "json"})
    with urllib.request.urlopen(q, timeout=60) as r:
        ids = json.loads(r.read())["esearchresult"]["idlist"]
    if not ids:
        return None
    time.sleep(0.4)
    q = f"{E}/esummary.fcgi?" + urllib.parse.urlencode({"db": "gds", "id": ids[0], "retmode": "json"})
    with urllib.request.urlopen(q, timeout=60) as r:
        doc = json.loads(r.read())["result"][ids[0]]
    for rel in doc.get("extrelations", []) or []:
        if str(rel.get("relationtype", "")).upper() == "SRA":
            return rel.get("targetobject")
    bp = doc.get("bioproject")
    return f"PRJNA{bp}" if bp and str(bp).isdigit() else bp


def plan(studies: list, strain_studies: set) -> None:
    path = os.path.join(NAS, "catalogue", "plan.jsonl")
    done = set()
    if os.path.exists(path):
        with open(path) as fh:
            done = {json.loads(ln)["run_accession"] for ln in fh if ln.strip()}
    with open(path, "a") as out:
        for study in studies:
            acc = geo_project(study) if study.startswith("GSE") else study
            if not acc:
                log(f"{study}: no SRA project found")
                continue
            time.sleep(0.4)
            runs = ena_runs(acc)
            refs = ["ME49", "ME49_host"] + (["GT1", "VEG"] if study in strain_studies else [])
            n = 0
            for r in runs:
                if r["run_accession"] in done or "Toxoplasma" not in r.get("scientific_name", "Toxoplasma"):
                    continue
                if r.get("library_strategy", "RNA-Seq") not in ("RNA-Seq", "OTHER", "Ribo-Seq", ""):
                    continue
                r.update({"study": study, "project": acc, "references": refs})
                out.write(json.dumps(r) + "\n")
                n += 1
            log(f"{study} ({acc}): {n} runs planned")


# --------------------------------------------------------------------------- quantification
def _free_gb(path: str) -> float:
    return shutil.disk_usage(path).free / 1e9


def quantify(limit: int | None = None) -> None:
    path = os.path.join(NAS, "catalogue", "plan.jsonl")
    with open(path) as fh:
        todo = [json.loads(ln) for ln in fh if ln.strip()]
    os.makedirs(SCRATCH, exist_ok=True)
    count = 0
    for r in todo:
        out_dir = os.path.join(NAS, "rnaseq", r["study"], r["run_accession"])
        if os.path.exists(os.path.join(out_dir, "done.json")):
            continue
        local = []
        _clear_scratch(r["run_accession"])          # whatever an interrupted attempt left behind
        try:
            local = _reads(r)
            if not local:
                log(f"{r['run_accession']}: no reads available from AWS or ENA; skipped")
                continue
            os.makedirs(out_dir, exist_ok=True)
            paired = len(local) == 2
            for ref in r["references"]:
                idx = os.path.join(NAS, "references", ref)
                s_out = os.path.join(out_dir, f"salmon_{ref}")
                reads = (["-1", local[0], "-2", local[1]] if paired else ["-r", local[0]])
                run([os.path.join(BIN, "salmon"), "quant", "-i", os.path.join(idx, "salmon_index"),
                     "-l", "A", *reads, "-p", str(THREADS), "--validateMappings", "--gcBias",
                     "--seqBias", "-o", s_out])
                _gzip(os.path.join(s_out, "quant.sf"))
                k_out = os.path.join(out_dir, f"kallisto_{ref}")
                k_reads = local if paired else ["--single", "-l", "200", "-s", "30", local[0]]
                run([os.path.join(BIN, "kallisto"), "quant", "-i", os.path.join(idx, "kallisto.idx"),
                     "-o", k_out, "-t", str(THREADS), *k_reads])
                _gzip(os.path.join(k_out, "abundance.tsv"))
                h5 = os.path.join(k_out, "abundance.h5")
                if os.path.exists(h5):
                    os.remove(h5)                      # the TSV holds the same numbers
            with open(os.path.join(out_dir, "done.json"), "w") as fh:
                json.dump({**r, "salmon": version("salmon"), "kallisto": version("kallisto"),
                           "finished": time.strftime("%Y-%m-%dT%H:%M:%S")}, fh, indent=2)
            log(f"{r['run_accession']} ({r['study']}): quantified against {', '.join(r['references'])}")
        except Exception as e:
            log(f"{r['run_accession']} ({r['study']}): FAILED {e}")
        finally:
            # Everything this run left on scratch -- reads, the .sra, fasterq-dump temp files -- whether
            # it finished or not: scratch is small and shared, and a failed run must not fill it.
            _clear_scratch(r["run_accession"])
        count += 1
        if limit and count >= limit:
            break


#: NCBI's copy of SRA on AWS (Open Data Program): measured 43 MB/s from here on 2026-10-17, against
#: 1.3 MB/s from ENA over HTTPS -- the difference between days and a month for the full set.
SRA_ODP = "https://sra-pub-run-odp.s3.amazonaws.com/sra/{run}/{run}"


def _head_size(url: str) -> int | None:
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=60) as r:
            return int(r.headers.get("Content-Length") or 0) or None
    except Exception:
        return None


def _reads(r: dict) -> list:
    """Local FASTQ files for one run: from AWS SRA (validated, then fasterq-dump), else ENA FASTQ."""
    run_acc = r["run_accession"]
    url = SRA_ODP.format(run=run_acc)
    size = _head_size(url)
    if size:
        # fasterq-dump writes uncompressed FASTQ, ~4-6x the .sra file, plus its temp files.
        need = size * 8 / 1e9 + SCRATCH_MARGIN_GB
        while _free_gb(SCRATCH) < need:
            log(f"{run_acc}: waiting for scratch space ({_free_gb(SCRATCH):.0f} GB free, {need:.0f} needed)")
            time.sleep(600)
        sra = os.path.join(SCRATCH, f"{run_acc}.sra")
        fetch(url, sra)
        run([os.path.join(BIN, "vdb-validate"), sra], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        run([os.path.join(BIN, "fasterq-dump"), "--split-3", "-e", str(THREADS), "-t", SCRATCH,
             "-O", SCRATCH, sra], stdout=subprocess.DEVNULL)
        os.remove(sra)
        pair = [os.path.join(SCRATCH, f"{run_acc}_{k}.fastq") for k in (1, 2)]
        single = os.path.join(SCRATCH, f"{run_acc}.fastq")
        if all(os.path.exists(f) for f in pair):
            if os.path.exists(single):
                os.remove(single)            # unpaired leftovers of a paired run
            return pair
        return [single] if os.path.exists(single) else []
    # Fallback: ENA's gzipped FASTQ, md5-checked.
    urls = [u for u in r.get("fastq_ftp", "").split(";") if u]
    md5s = [m for m in r.get("fastq_md5", "").split(";") if m]
    sizes = [int(b) for b in r.get("fastq_bytes", "").split(";") if b.isdigit()]
    if not urls:
        return []
    need = sum(sizes) / 1e9 + SCRATCH_MARGIN_GB
    while _free_gb(SCRATCH) < need:
        log(f"{run_acc}: waiting for scratch space ({_free_gb(SCRATCH):.0f} GB free, {need:.0f} needed)")
        time.sleep(600)
    pairs = [u for u in urls if u.endswith(("_1.fastq.gz", "_2.fastq.gz"))] or urls
    local = []
    for u in pairs:
        dest = os.path.join(SCRATCH, os.path.basename(u))
        fetch("https://" + u if not u.startswith("http") else u, dest)
        i = urls.index(u)
        if i < len(md5s):
            md5 = hashlib.md5()
            with open(dest, "rb") as fh:
                for block in iter(lambda: fh.read(1 << 22), b""):
                    md5.update(block)
            if md5.hexdigest() != md5s[i]:
                raise RuntimeError(f"md5 mismatch for {u}")
        local.append(dest)
    return local

def _clear_scratch(run_accession: str) -> None:
    """Remove one run's files (reads, .sra, fasterq-dump temp) from scratch."""
    for f in os.listdir(SCRATCH):
        if f.startswith(run_accession) or f.startswith("fasterq.tmp"):
            path = os.path.join(SCRATCH, f)
            shutil.rmtree(path, ignore_errors=True) if os.path.isdir(path) else os.remove(path)

def _gzip(path: str) -> None:
    if os.path.exists(path):
        with open(path, "rb") as fh, gzip.open(path + ".gz", "wb") as out:
            shutil.copyfileobj(fh, out)
        os.remove(path)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("references")
    p = sub.add_parser("runs")
    p.add_argument("studies", nargs="+")
    p.add_argument("--strain", nargs="*", default=[], help="studies to also quantify per strain")
    q = sub.add_parser("quantify")
    q.add_argument("--limit", type=int, default=None)
    args = ap.parse_args(argv)
    if args.cmd == "references":
        references()
        host_reference()
    elif args.cmd == "runs":
        plan(args.studies, set(args.strain))
    else:
        quantify(args.limit)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
