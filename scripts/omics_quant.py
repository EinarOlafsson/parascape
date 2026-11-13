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
import glob
import gzip
import hashlib
import json
import os
import re
from pathlib import Path
import shlex
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
RUNTIME_SOURCE_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
DOWNLOAD_CACHE = os.environ.get("PARASCAPE_DOWNLOAD_CACHE", "")
DOWNLOAD_RATE = int(os.environ.get("PARASCAPE_DOWNLOAD_BYTES_PER_SECOND", "0"))
PLAN = os.environ.get("PARASCAPE_PLAN", os.path.join(NAS, "catalogue", "plan.jsonl"))

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
#: Further single-host competitive references (resolved from the same Ensembl release-116 FTP listing on
#: 2026-11-01). ME49_host stays human + mouse so existing results keep their meaning; dual RNA-seq from
#: another host needs that host's transcripts present, or its conserved genes can be counted as parasite.
HOST_REFERENCES = {
    "ME49_host": HOST_TRANSCRIPTS,
    "ME49_pig": (f"{ENSEMBL_MAIN}/sus_scrofa/cdna/Sus_scrofa.Sscrofa11.1.cdna.all.fa.gz",
                 f"{ENSEMBL_MAIN}/sus_scrofa/ncrna/Sus_scrofa.Sscrofa11.1.ncrna.fa.gz"),
    "ME49_cow": (f"{ENSEMBL_MAIN}/bos_taurus/cdna/Bos_taurus.ARS-UCD2.0.cdna.all.fa.gz",
                 f"{ENSEMBL_MAIN}/bos_taurus/ncrna/Bos_taurus.ARS-UCD2.0.ncrna.fa.gz"),
    "ME49_cat": (f"{ENSEMBL_MAIN}/felis_catus/cdna/Felis_catus.F.catus_Fca126_mat1.0.cdna.all.fa.gz",
                 f"{ENSEMBL_MAIN}/felis_catus/ncrna/Felis_catus.F.catus_Fca126_mat1.0.ncrna.fa.gz"),
    "ME49_hamster": (f"{ENSEMBL_MAIN}/mesocricetus_auratus/cdna/Mesocricetus_auratus.MesAur1.0.cdna.all.fa.gz",
                     f"{ENSEMBL_MAIN}/mesocricetus_auratus/ncrna/Mesocricetus_auratus.MesAur1.0.ncrna.fa.gz"),
    "ME49_rat": (f"{ENSEMBL_MAIN}/rattus_norvegicus/cdna/Rattus_norvegicus.GRCr8.cdna.all.fa.gz",
                 f"{ENSEMBL_MAIN}/rattus_norvegicus/ncrna/Rattus_norvegicus.GRCr8.ncrna.fa.gz"),
}
#: ENA scientific names served by each competitive host reference.
HOST_SPECIES = {"Homo sapiens": "ME49_host", "Mus musculus": "ME49_host", "Sus scrofa": "ME49_pig",
                "Bos taurus": "ME49_cow", "Felis catus": "ME49_cat", "Mesocricetus auratus": "ME49_hamster",
                "Rattus norvegicus": "ME49_rat"}
#: Gene identifiers given to host transcripts, so a matrix can drop them (and report the host share).
HOST_PREFIX = "HOST:"

#: Free space a download may leave on scratch, at least.
SCRATCH_MARGIN_GB = float(os.environ.get("PARASCAPE_SCRATCH_MARGIN_GB", "25"))  # free space kept on the scratch disk


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
                while True:
                    started = time.monotonic()
                    block = r.read(1 << 20)
                    if not block:
                        break
                    fh.write(block)
                    if DOWNLOAD_RATE > 0:
                        time.sleep(max(0, len(block) / DOWNLOAD_RATE - (time.monotonic() - started)))
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
    out = subprocess.run([os.path.join(BIN, tool), "--version" if tool in ("salmon", "fastp") else "version"],
                         capture_output=True, text=True)
    return (out.stdout or out.stderr).strip()


# --------------------------------------------------------------------------- references
def _salmon_index_present(directory: str) -> bool:
    """Salmon 1.x wrote versionInfo.json; salmon 2.x (piscem-based) writes info.json instead."""
    return any(os.path.exists(os.path.join(directory, "salmon_index", name))
               for name in ("versionInfo.json", "info.json"))


def references() -> None:
    for name, (cdna_url, dna_url) in REFERENCES.items():
        d = os.path.join(NAS, "references", name)
        os.makedirs(d, exist_ok=True)
        meta_path = os.path.join(d, "reference.json")
        if _salmon_index_present(d) and os.path.exists(os.path.join(d, "kallisto.idx")):
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



def host_reference(name: str = "ME49_host") -> None:
    """ME49 + one host set's transcripts (with ncRNA), ME49 genome as decoys: a competitive reference."""
    host_urls = HOST_REFERENCES[name]
    d = os.path.join(NAS, "references", name)
    if _salmon_index_present(d) and os.path.exists(os.path.join(d, "kallisto.idx")):
        log(f"{name}: indexes present")
        return
    os.makedirs(d, exist_ok=True)
    me49 = os.path.join(NAS, "references", "ME49")
    cdna_url, dna_url = REFERENCES["ME49"]
    parts = [os.path.join(me49, os.path.basename(cdna_url))]
    meta = {"name": name, "parasite": cdna_url, "decoys": dna_url, "host": list(host_urls),
            "sha256": {}}
    for url in host_urls:
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
    targets = os.path.join(SCRATCH, f"{name}_targets.fa.gz")
    os.makedirs(SCRATCH, exist_ok=True)
    with open(targets, "wb") as out:
        for part in parts:
            with open(part, "rb") as fh:
                shutil.copyfileobj(fh, out)
    gentrome = os.path.join(SCRATCH, f"{name}_gentrome.fa.gz")
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
    log(f"{name}: indexed")

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
#: Optional second scratch on the NAS (Einar 2026-11-02: use NAS space when local space is the bottleneck).
#: Used only for a run that does not fit locally, and only while the NAS keeps this much free afterwards.
NAS_SCRATCH = os.environ.get("PARASCAPE_NAS_SCRATCH", "")
NAS_SCRATCH_RESERVE_GB = float(os.environ.get("PARASCAPE_NAS_SCRATCH_RESERVE_GB", "150"))
#: Runs held back by hand (comma-separated accessions), e.g. a giant run whose only feasible read source is
#: stalled; they stay pending in the plan and are reported, so the rest of the queue is never blocked.
DEFER_RUNS = {r.strip() for r in os.environ.get("PARASCAPE_DEFER_RUNS", "").split(",") if r.strip()}
#: Gzipped FASTQ made by convert_sra_to_gzip.py for runs too large to extract; used only with its done.json.
PRECONVERTED = os.environ.get("PARASCAPE_PRECONVERTED", "/nas_mnt/carruthers/Einar/parascape_omics/preconverted_fastq")
#: Extra salmon outputs (Einar 2026-11-03), so runs can be revisited without their raw reads: equivalence
#: classes with per-transcript weights (re-quantification against an updated annotation) and 30 Gibbs
#: posterior samples (per-transcript uncertainty). Runs finished before this are topped up later.
SALMON_EXTRA = ["--dumpEqWeights", "--numGibbsSamples", "30"]
#: Parallel quantifier workers in one service (Einar 2026-11-03: two). Each worker has its own scratch
#: folder, active file and THREADS // WORKERS threads; runs are claimed atomically so none runs twice.
WORKERS = max(1, int(os.environ.get("PARASCAPE_WORKERS", "1")))
#: Workers beyond WORKERS that claim runs only while the resource guard reports spare CPUs ('wide' mode:
#: Einar 2026-11-03, use the cores autotrade is not using). They get the same thread count as the others.
EXTRA_WORKERS = max(0, int(os.environ.get("PARASCAPE_EXTRA_WORKERS", "0")))
RESOURCE_STATUS = os.environ.get("PARASCAPE_RESOURCE_STATUS",
                                 os.path.expanduser("~/.cache/parascape-downloads/resource-status.json"))
EXTRA_RANK = 0          # set in an extra worker: 1 for the first extra, 2 for the second, ...


#: Salmon against these references runs in a later pass (Einar 2026-11-13, "Host salmon later"). Salmon
#: against the combined ME49_host index (1.14 M references) maps some projects at about 8,300 fragments a
#: second (GSE210671, GSE223819, GSE223621), so one deep run held a worker for more than 7 hours while the
#: cheap stages of hundreds of other runs waited. Now every run first gets its cheap stages -- archive check,
#: reads, fastp, salmon on the parasite-only indexes, kallisto on every index -- and the host salmon follows
#: once no run is left without them. Settings and outputs are unchanged; only the order is. The reads are not
#: kept between the passes (all queued runs together hold ~3.5 TB of reads); the later pass re-derives them
#: from the archive and checks that fastp keeps the same number of reads. Empty: one pass, as before.
HOST_SALMON_LATER = {ref.strip() for ref in os.environ.get("PARASCAPE_HOST_SALMON_LATER", "ME49_host").split(",")
                     if ref.strip()}
#: Studies whose host salmon goes to the later pass (comma-separated); empty means every study (the decision of
#: Einar). Re-deriving the reads costs more than host salmon itself for most projects (216 against 347 M
#: reads an hour per worker, 2026-11-06 to 11-13), so this can narrow the change to the slow projects.
HOST_SALMON_LATER_STUDIES = {s.strip() for s in os.environ.get("PARASCAPE_HOST_SALMON_LATER_STUDIES", "").split(",")
                             if s.strip()}
#: Written when a run's cheap stages are done and its host salmon is still to run. done.json is written only
#: when every stage is done, so everything that counts or reads done.json (queue verification, processing
#: ledger, app snapshots, re-quantification) still sees an unfinished run as unfinished.
STAGES_FILE = "stages.json"
CHEAP, HOST, ALL = "cheap stages", "host salmon", "all stages"
#: Finish the current run, then stop taking new ones while this file exists -- so a code change can reach the
#: workers between runs instead of killing work in flight. `quantify` removes it when it starts again.
DRAIN_FILE = os.environ.get("PARASCAPE_DRAIN_FILE", os.path.join(SCRATCH, "DRAIN"))


def _draining() -> bool:
    return bool(DRAIN_FILE) and os.path.exists(DRAIN_FILE)


def _active_elsewhere(run_acc: str) -> bool:
    """True when another live quantifier process (another service, or a bridge started during a code change)
    names this run in its active file; claims only coordinate the workers of one process."""
    if not DOWNLOAD_CACHE:
        return False
    for path in glob.glob(os.path.join(DOWNLOAD_CACHE, "active*.json")):
        if os.path.basename(path) == ACTIVE_NAME:
            continue
        try:
            with open(path) as fh:
                active = json.load(fh)
            if active.get("run") != run_acc or int(active.get("pid", 0)) == os.getpid():
                continue
            with open(f"/proc/{int(active['pid'])}/cmdline", "rb") as fh:
                if b"omics_quant.py" in fh.read():
                    return True
        except (OSError, ValueError, KeyError, TypeError):
            continue
    return False


_PLANNED = {"epoch": 0.0, "runs": None, "mtime": None}


def _still_planned(run_acc: str) -> bool:
    """Workers hold the queue in memory; a run removed from the plan since then is skipped (checked every 10 min)."""
    if time.time() - _PLANNED["epoch"] > 600:
        try:
            mtime = os.stat(PLAN).st_mtime_ns
            if _PLANNED["mtime"] is None:
                _PLANNED["mtime"] = mtime            # the plan this worker started from: nothing to skip
            elif mtime != _PLANNED["mtime"]:
                with open(PLAN) as fh:
                    _PLANNED["runs"] = {json.loads(line)["run_accession"] for line in fh if line.strip()}
                _PLANNED["mtime"] = mtime
        except (OSError, ValueError):
            pass                                     # keep the previous set if the plan is mid-rewrite
        _PLANNED["epoch"] = time.time()
    return _PLANNED["runs"] is None or run_acc in _PLANNED["runs"]


def _spare_workers() -> int:
    try:
        with open(RESOURCE_STATUS) as fh:
            return int(json.load(fh).get("extra_workers", 0))
    except (OSError, ValueError):
        return 0
ACTIVE_NAME = os.environ.get("PARASCAPE_ACTIVE_NAME", "active.json")
CLAIMS = ""                                  # claims folder, set in multi-worker mode only


class ScratchDeferred(Exception):
    """A run needs more scratch than is free now; quantify others and retry it later."""

    def __init__(self, message: str, need: float = 0.0):
        super().__init__(message)
        self.need = need


def _reserved_by_others(run_acc: str) -> float:
    """Scratch GB other workers have reserved on the shared disk (0 with a single worker)."""
    total = 0.0
    if CLAIMS and os.path.isdir(CLAIMS):
        for name in os.listdir(CLAIMS):
            if name != f"{run_acc}.claim":
                try:
                    total += float(open(os.path.join(CLAIMS, name)).read() or 0)
                except (OSError, ValueError):
                    continue
    return total


def _need_scratch(run_acc: str, need: float) -> None:
    # Waiting here blocked the whole queue behind one large run (2026-11-01: 328 GB needed, 107 free).
    free = _free_gb(SCRATCH) - (_reserved_by_others(run_acc) if SCRATCH != NAS_SCRATCH else 0.0)
    if free < need:
        raise ScratchDeferred(f"{run_acc}: deferred, scratch {free:.0f} GB free after other workers, {need:.0f} needed", need)
    if CLAIMS and SCRATCH != NAS_SCRATCH:
        with open(os.path.join(CLAIMS, f"{run_acc}.claim"), "w") as fh:   # reserve without the shared margin
            fh.write(f"{max(0.0, need - SCRATCH_MARGIN_GB):.1f}")


def _claim(run_acc: str) -> bool:
    """Atomically take a run for this worker; False if another worker has it."""
    if not CLAIMS:
        return True
    try:
        os.close(os.open(os.path.join(CLAIMS, f"{run_acc}.claim"), os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        return True
    except FileExistsError:
        return False


def _release(run_acc: str) -> None:
    if CLAIMS:
        Path(CLAIMS, f"{run_acc}.claim").unlink(missing_ok=True)


def _nas_scratch_fits(need: float) -> bool:
    """True when the run's own scratch use (need without the local margin) leaves the NAS reserve free."""
    if not NAS_SCRATCH:
        return False
    os.makedirs(NAS_SCRATCH, exist_ok=True)
    return _free_gb(NAS_SCRATCH) - (need - SCRATCH_MARGIN_GB) >= NAS_SCRATCH_RESERVE_GB


def _free_gb(path: str) -> float:
    return shutil.disk_usage(path).free / 1e9


# --------------------------------------------------------------------------- re-quantification
#: Re-quantification of runs verified before trimming and the salmon extras (instruction 81) is written
#: beside the verified result, never over it: rnaseq/<study>/<run>/requant/. Set by `requantify`.
REQUANT = {"subdir": ""}
REQUANT_DIR = "requant"
#: A re-quantified run whose transcripts agree with the verified ones at least this well (Spearman of TPM over
#: transcripts expressed in either) agrees; below it, the difference is reported, not overwritten.
REQUANT_AGREE = 0.9


def _tpm(path: str):
    import pandas as pd
    for name in ("quant.sf.gz", "quant.sf"):
        f = os.path.join(path, name)
        if os.path.exists(f):
            return pd.read_csv(f, sep="\t", usecols=["Name", "TPM"]).set_index("Name")["TPM"]
    return None


def compare_requant(run_dir: str, new_dir: str, references: list) -> dict:
    """Compare a re-quantification with the verified result and write requant/compare.json."""
    from scipy.stats import spearmanr
    out = {"compared": time.strftime("%Y-%m-%dT%H:%M:%S"), "references": {}}
    worst = None
    for ref in references:
        old = _tpm(os.path.join(run_dir, f"salmon_{ref}"))
        new = _tpm(os.path.join(new_dir, f"salmon_{ref}"))
        if old is None or new is None:
            continue
        both = old.to_frame("old").join(new.to_frame("new"), how="inner")
        expressed = both[(both["old"] > 1) | (both["new"] > 1)]
        old_mapped = salmon_mapped(os.path.join(run_dir, f"salmon_{ref}"))
        new_mapped = salmon_mapped(os.path.join(new_dir, f"salmon_{ref}"))
        rho = float(spearmanr(expressed["old"], expressed["new"]).correlation) if len(expressed) > 10 else None
        out["references"][ref] = {"spearman": rho, "transcripts": int(len(expressed)),
                                  "mapped_before": old_mapped, "mapped_now": new_mapped}
        if rho is not None and (worst is None or rho < worst):
            worst = rho
    zero_before = any((v.get("mapped_before") or 0) == 0 and (v.get("mapped_now") or 0) > 0
                      for v in out["references"].values())
    out["spearman"] = None if worst is None else round(worst, 4)
    out["verdict"] = ("recovered: mapped nothing before trimming, maps now" if zero_before
                      else "agrees" if worst is not None and worst >= REQUANT_AGREE
                      else "DIFFERS -- review before use" if worst is not None else "not comparable")
    _cache_module().save(Path(new_dir) / "compare.json", out)
    return out


def requant_todo() -> list:
    """Runs verified before trimming (no `trimming` in done.json) and not yet re-quantified, the runs
    that mapped nothing first -- they contribute nothing until this runs."""
    with open(PLAN) as fh:
        plan_rows = [json.loads(ln) for ln in fh if ln.strip()]
    first, rest = [], []
    for r in plan_rows:
        run_dir = os.path.join(NAS, "rnaseq", r["study"], r["run_accession"])
        done = os.path.join(run_dir, "done.json")
        if not os.path.exists(done) or os.path.exists(os.path.join(run_dir, REQUANT_DIR, "done.json")):
            continue
        try:
            with open(done) as fh:
                if "trimming" in json.load(fh):
                    continue
        except (OSError, ValueError):
            continue
        mapped = salmon_mapped(os.path.join(run_dir, f"salmon_{r['references'][0]}"))
        (first if mapped == 0 else rest).append(r)
    return first + rest


def pending_new_runs() -> int:
    """Planned runs not yet quantified (nor recorded as colour space): the new downloads."""
    with open(PLAN) as fh:
        plan_rows = [json.loads(ln) for ln in fh if ln.strip()]
    n = 0
    for r in plan_rows:
        run_dir = os.path.join(NAS, "rnaseq", r["study"], r["run_accession"])
        if not (os.path.exists(os.path.join(run_dir, "done.json"))
                or os.path.exists(os.path.join(run_dir, COLOUR_SPACE_FILE))):
            n += 1
    return n


def requantify(limit: int | None = None, wait: bool = True) -> None:
    """Re-quantify the earlier runs (instruction 81) -- only once the new downloads are done, as Einar
    asked: until then this waits, checking hourly."""
    while wait and pending_new_runs() > 0:
        log(f"re-quantification waits: {pending_new_runs()} new runs still to quantify")
        time.sleep(3600)
    REQUANT["subdir"] = REQUANT_DIR
    todo = requant_todo()
    log(f"re-quantifying {len(todo)} runs verified before trimming")
    os.makedirs(SCRATCH, exist_ok=True)
    count = 0
    while todo:
        deferred = _quantify_pass(todo, limit, count)
        if deferred is None or not deferred:
            return
        time.sleep(600)
        todo = deferred


def quantify(limit: int | None = None) -> None:
    path = PLAN
    with open(path) as fh:
        todo = [json.loads(ln) for ln in fh if ln.strip()]
    os.makedirs(SCRATCH, exist_ok=True)
    if DOWNLOAD_CACHE:
        os.makedirs(DOWNLOAD_CACHE, exist_ok=True)
    if _draining():
        os.remove(DRAIN_FILE)                       # a new start is the deployment the drain waited for
        log(f"drain marker {DRAIN_FILE} removed at start")
    if WORKERS > 1:
        return _quantify_workers(todo, limit)
    _run_phases(todo, limit)


def _run_phases(todo: list, limit: int | None) -> None:
    """Cheap stages of every queued run first, then the host salmon (Einar 2026-11-13); one pass when
    HOST_SALMON_LATER is empty. Runs deferred for scratch space are retried every 10 minutes; while only
    such runs are left in the cheap pass, one host-salmon run is done between retries."""
    phase, queue = (CHEAP if HOST_SALMON_LATER else ALL), todo
    while True:
        deferred = _quantify_pass(queue, limit, 0, phase)
        if deferred is None or _draining():
            return
        if not deferred:
            if phase != CHEAP:
                return
            phase, queue = HOST, todo               # every run has its cheap stages: the later pass
            log("cheap stages done for every run this worker can take; starting the host salmon pass")
            continue
        log(f"{len(deferred)} runs deferred for scratch space ({phase}); retrying in 10 minutes")
        if phase == CHEAP and _quantify_pass(todo, 1, 0, HOST) is None:
            queue = deferred                         # one host run done meanwhile; retry at once
            continue
        time.sleep(600)
        queue = deferred


def _quantify_workers(todo: list, limit: int | None) -> None:
    """Run WORKERS forked workers over one shared queue; the service fails if any worker fails."""
    import multiprocessing
    base = SCRATCH
    claims = os.path.join(base, "claims")
    shutil.rmtree(claims, ignore_errors=True)       # claims of a previous service instance are stale
    os.makedirs(claims)
    context = multiprocessing.get_context("fork")
    workers = [context.Process(target=_worker, args=(k, todo, limit, base, claims))
               for k in range(WORKERS + EXTRA_WORKERS)]
    for w in workers:
        w.start()
    for w in workers:
        w.join()
    if any(w.exitcode for w in workers):
        raise SystemExit(f"quantifier worker failed: {[w.exitcode for w in workers]}")


def _worker(k: int, todo: list, limit: int | None, base: str, claims: str) -> None:
    global SCRATCH, THREADS, ACTIVE_NAME, CLAIMS, EXTRA_RANK
    SCRATCH, ACTIVE_NAME, CLAIMS = os.path.join(base, f"w{k}"), f"active-w{k}.json", claims
    THREADS = max(1, THREADS // WORKERS)
    EXTRA_RANK = max(0, k - WORKERS + 1)
    os.makedirs(SCRATCH, exist_ok=True)
    log(f"worker {k}: scratch {SCRATCH}, {THREADS} threads" + (" (only while CPUs are spare)" if EXTRA_RANK else "")
        + (f"; salmon against {', '.join(sorted(HOST_SALMON_LATER))} in a later pass" if HOST_SALMON_LATER else ""))
    _run_phases(todo, limit)


def _trim(run_acc: str, local: list) -> tuple:
    """Cut adapters before quantification (Einar 2026-11-04).

    Some deposits are short-insert libraries left untrimmed: GSE243206/243207/245775 carry adapter in 99.5 % of
    reads, their 150 nt reads hold a ~44 nt insert, and salmon mapped 0.0 % of them while kallisto still found
    10-32 %. fastp detects the adapter from the data; a normal RNA-seq library loses almost nothing.
    """
    report = os.path.join(SCRATCH, f"{run_acc}.fastp.json")
    out = [os.path.join(SCRATCH, f"{run_acc}_trimmed_{i + 1}.fastq.gz") for i in range(len(local))]
    command = [os.path.join(BIN, "fastp"), "-w", str(min(4, THREADS)), "-l", "20", "-j", report, "-h", "/dev/null"]
    if len(local) == 2:
        command += ["-i", local[0], "-I", local[1], "-o", out[0], "-O", out[1], "--detect_adapter_for_pe"]
    else:
        command += ["-i", local[0], "-o", out[0]]
    run(command)
    with open(report) as fh:
        summary = json.load(fh)
    before, after = summary["summary"]["before_filtering"], summary["summary"]["after_filtering"]
    kept = dict(fastp=version("fastp"), reads_before=before["total_reads"], reads_after=after["total_reads"],
                read1_mean_length_before=before.get("read1_mean_length"),
                read1_mean_length_after=after.get("read1_mean_length"),
                adapter_trimmed_reads=summary.get("adapter_cutting", {}).get("adapter_trimmed_reads"))
    share = (kept["adapter_trimmed_reads"] or 0) / max(1, before["total_reads"])
    log(f"{run_acc}: trimmed adapters from {share:.1%} of reads; mean read length "
        f"{kept['read1_mean_length_before']} -> {kept['read1_mean_length_after']}")
    for path in local:                                  # the untrimmed copies are not needed again
        try:
            os.remove(path)
        except OSError:
            pass
    return out, kept


#: Instruments whose reads are in colour space (two-base encoding: "T0123..."), which no k-mer
#: quantifier can read. GSE60305 (AB SOLiD 4) failed kallisto 49 times in a day and "succeeded" in
#: salmon with 0 of 47 million reads mapped before this was seen (instruction 89).
COLOUR_SPACE_INSTRUMENTS = ("SOLiD",)
COLOUR_SPACE_FILE = "colour_space.json"


def colour_space(r: dict) -> bool:
    """Whether a planned run's reads are colour-space SOLiD reads."""
    model = str(r.get("instrument_model", "") or "")
    return any(name.lower() in model.lower() for name in COLOUR_SPACE_INSTRUMENTS)


def salmon_mapped(s_out: str) -> int | None:
    """Reads salmon mapped, from its meta_info.json (None when it is not there)."""
    path = os.path.join(s_out, "aux_info", "meta_info.json")
    try:
        with open(path) as fh:
            return int(json.load(fh).get("num_mapped", 0))
    except (OSError, ValueError, TypeError):
        return None


class NothingMapped(RuntimeError):
    """A quantifier finished but mapped no read: a failure, whatever its exit status says."""


def _salmon(ref: str, local: list, s_out: str) -> int | None:
    """Salmon against one reference; returns the reads it mapped (a run that maps none has failed)."""
    reads = (["-1", local[0], "-2", local[1]] if len(local) == 2 else ["-r", local[0]])
    # An interrupted salmon leaves its mapping file (GBs); this worker holds the run's claim.
    for stale in glob.glob(os.path.join(s_out, "intermediate_mappings.rad.partial-*")):
        os.remove(stale)
    run([os.path.join(BIN, "salmon"), "quant", "-i", os.path.join(NAS, "references", ref, "salmon_index"),
         "-l", "A", *reads, "-p", str(THREADS), "--validateMappings", "--gcBias",
         "--seqBias", *SALMON_EXTRA, "-o", s_out])
    _gzip(os.path.join(s_out, "quant.sf"))
    mapped = salmon_mapped(s_out)
    if mapped == 0:
        raise NothingMapped(f"salmon mapped 0 reads against {ref}: the reads are not "
                            f"this organism, or not in base space")
    return mapped


def _kallisto(ref: str, local: list, k_out: str) -> str | None:
    """Kallisto against one reference; returns the failure to record, or None."""
    k_reads = local if len(local) == 2 else ["--single", "-l", "200", "-s", "30", local[0]]
    try:
        run([os.path.join(BIN, "kallisto"), "quant", "-i", os.path.join(NAS, "references", ref, "kallisto.idx"),
             "-o", k_out, "-t", str(THREADS), *k_reads])
    except subprocess.CalledProcessError as error:
        return f"exit {error.returncode}"
    _gzip(os.path.join(k_out, "abundance.tsv"))
    h5 = os.path.join(k_out, "abundance.h5")
    if os.path.exists(h5):
        os.remove(h5)                      # the TSV holds the same numbers
    return None


def _skip(phase: str, run_dir: str, out_dir: str) -> bool:
    """Whether this pass leaves a run alone: finished, colour space, or not at this pass's stage."""
    if os.path.exists(os.path.join(out_dir, "done.json")) or os.path.exists(os.path.join(run_dir, COLOUR_SPACE_FILE)):
        return True
    staged = os.path.exists(os.path.join(out_dir, STAGES_FILE))
    return staged if phase == CHEAP else (not staged) if phase == HOST else False


def _quantify_pass(todo: list, limit: int | None, count: int, phase: str = ALL):
    """One pass over the queue; returns runs deferred for scratch space, or None when the limit is hit
    or a drain was asked for. `phase`: CHEAP (everything but the salmon in HOST_SALMON_LATER), HOST (that
    salmon, for runs whose cheap stages are done) or ALL (every stage, as before)."""
    global SCRATCH
    local_scratch = SCRATCH
    deferred = []
    for r in todo:
        run_dir = os.path.join(NAS, "rnaseq", r["study"], r["run_accession"])
        out_dir = os.path.join(run_dir, REQUANT["subdir"]) if REQUANT["subdir"] else run_dir
        if _skip(phase, run_dir, out_dir):
            continue
        if _draining():
            log(f"drain marker {DRAIN_FILE} present: no new run taken")
            return None
        if colour_space(r):
            os.makedirs(out_dir, exist_ok=True)
            _cache_module().save(Path(out_dir) / COLOUR_SPACE_FILE, {
                **r, "status": "not quantified",
                "reason": "colour-space reads (AB SOLiD); kallisto and salmon cannot read them -- "
                          "needs a colour-space aligner or the study's processed table",
                "recorded": time.strftime("%Y-%m-%dT%H:%M:%S")})
            log(f"{r['run_accession']} ({r['study']}): colour-space SOLiD reads; recorded, not retried")
            continue
        if not _still_planned(r["run_accession"]):
            continue                                 # removed from the plan after this worker started
        if r["run_accession"] in DEFER_RUNS:
            log(f"{r['run_accession']} ({r['study']}): held back by PARASCAPE_DEFER_RUNS; stays pending")
            continue
        while EXTRA_RANK and _spare_workers() < EXTRA_RANK:
            time.sleep(60)                           # wait for spare CPUs; others may finish this run meanwhile
            if os.path.exists(os.path.join(out_dir, "done.json")) or (
                    CLAIMS and os.path.exists(os.path.join(CLAIMS, f"{r['run_accession']}.claim"))):
                break
        if _active_elsewhere(r["run_accession"]):
            continue                                 # another quantifier process has this run
        if not _claim(r["run_accession"]):
            continue                                 # another worker has this run
        if _skip(phase, run_dir, out_dir):
            _release(r["run_accession"])             # moved on by another worker after the check above
            continue
        local = []
        if DOWNLOAD_CACHE:
            _cache_module().save(Path(DOWNLOAD_CACHE) / ACTIVE_NAME,
                                 {'run': r['run_accession'], 'study': r['study'], 'pid': os.getpid()})
        _clear_scratch(None)                        # whatever interrupted attempts left behind
        try:
            try:
                local = _reads(r)
            except ScratchDeferred as e:
                if not _nas_scratch_fits(e.need):
                    raise
                log(f"{r['run_accession']}: local scratch too small; extracting on NAS scratch {NAS_SCRATCH}")
                SCRATCH = NAS_SCRATCH
                _clear_scratch(r["run_accession"])
                local = _reads(r)
            if not local:
                log(f"{r['run_accession']}: no reads available from AWS or ENA; skipped")
                continue
            os.makedirs(out_dir, exist_ok=True)
            local, trim = _trim(r["run_accession"], local)
            if phase == HOST:
                _host_salmon(r, out_dir, local, trim)
            else:
                _cheap_or_all(r, out_dir, local, trim, phase)
            if REQUANT["subdir"] and os.path.exists(os.path.join(out_dir, "done.json")):
                verdict = compare_requant(run_dir, out_dir, r["references"])
                log(f"{r['run_accession']}: re-quantified -- {verdict['verdict']} "
                    f"(Spearman {verdict.get('spearman')})")
        except ScratchDeferred as e:
            log(str(e))
            deferred.append(r)
            continue
        except Exception as e:
            log(f"{r['run_accession']} ({r['study']}): FAILED {e}")
        finally:
            # Everything this run left on scratch -- reads, the .sra, fasterq-dump temp files -- whether
            # it finished or not: scratch is small and shared, and a failed run must not fill it.
            _clear_scratch(r["run_accession"])
            if NAS_SCRATCH:
                shutil.rmtree(Path(NAS_SCRATCH) / r["run_accession"], ignore_errors=True)
            SCRATCH = local_scratch                     # a NAS-scratch run never changes the next run's scratch
            if DOWNLOAD_CACHE:
                (Path(DOWNLOAD_CACHE) / ACTIVE_NAME).unlink(missing_ok=True)
            _release(r["run_accession"])
        count += 1
        if limit and count >= limit:
            return None
    return deferred


def _cheap_or_all(r: dict, out_dir: str, local: list, trim: dict, phase: str) -> None:
    """Every stage (ALL), or every stage but the later-pass salmon (CHEAP; then stages.json, not done.json)."""
    deferred_study = not HOST_SALMON_LATER_STUDIES or r["study"] in HOST_SALMON_LATER_STUDIES
    later = [ref for ref in r["references"] if ref in HOST_SALMON_LATER] if phase == CHEAP and deferred_study else []
    kallisto_failed, stages = {}, ["archive", "reads", "fastp"]
    for ref in r["references"]:
        mapped = None
        if ref not in later:
            mapped = _salmon(ref, local, os.path.join(out_dir, f"salmon_{ref}"))
            stages.append(f"salmon_{ref}")
        failed = _kallisto(ref, local, os.path.join(out_dir, f"kallisto_{ref}"))
        if failed:
            # Salmon mapped reads, so the run is kept on its numbers and kallisto's failure
            # recorded: one tool must not lose a whole study again (instruction 89).
            kallisto_failed[ref] = failed
            log(f"{r['run_accession']}: kallisto failed against {ref} ({failed}); kept on salmon"
                + (f" ({mapped:,} reads mapped)" if mapped is not None else ""))
        else:
            stages.append(f"kallisto_{ref}")
    record = {**r, "trimming": trim, "salmon": version("salmon"), "kallisto": version("kallisto"),
              "salmon_extra": SALMON_EXTRA, "kallisto_failed": kallisto_failed,
              "configured_threads": THREADS, "runtime_source_sha256": RUNTIME_SOURCE_SHA256}
    now = time.strftime("%Y-%m-%dT%H:%M:%S")
    if later:
        _cache_module().save(Path(out_dir) / STAGES_FILE, {
            **record, "stages_done": stages, "stages_pending": [f"salmon_{ref}" for ref in later],
            "status": "cheap stages done; host salmon pending (Einar 2026-11-13: host salmon later)",
            "cheap_stages_finished": now})
        log(f"{r['run_accession']} ({r['study']}): cheap stages done; salmon against {', '.join(later)} "
            f"waits for the later pass")
        return
    _cache_module().save(Path(out_dir) / "done.json",
                         {**record, "stages_done": stages, "stages_pending": [], "finished": now})
    log(f"{r['run_accession']} ({r['study']}): quantified against {', '.join(r['references'])}")


def _host_salmon(r: dict, out_dir: str, local: list, trim: dict) -> None:
    """The later pass: salmon against the deferred references, on reads re-derived from the archive."""
    with open(os.path.join(out_dir, STAGES_FILE)) as fh:
        staged = json.load(fh)
    before = staged["trimming"]
    if (trim.get("reads_before"), trim.get("reads_after")) != (before.get("reads_before"), before.get("reads_after")):
        raise RuntimeError(f"re-derived reads differ from the first pass (fastp kept {trim.get('reads_after')} "
                           f"of {trim.get('reads_before')}, before {before.get('reads_after')} of "
                           f"{before.get('reads_before')}); host salmon not run")
    stages = list(staged["stages_done"])
    for stage in staged["stages_pending"]:
        ref = stage[len("salmon_"):]
        _salmon(ref, local, os.path.join(out_dir, stage))
        stages.append(stage)
    now = time.strftime("%Y-%m-%dT%H:%M:%S")
    done = {key: value for key, value in staged.items() if key not in ("status", "stages_pending")}
    done.update({"stages_done": stages, "stages_pending": [], "host_salmon_later": True,
                 "host_salmon_finished": now, "host_salmon_threads": THREADS,
                 "host_salmon_runtime_source_sha256": RUNTIME_SOURCE_SHA256, "finished": now})
    _cache_module().save(Path(out_dir) / "done.json", done)
    Path(out_dir, STAGES_FILE).unlink(missing_ok=True)
    log(f"{r['run_accession']} ({r['study']}): host salmon done; quantified against {', '.join(r['references'])}")


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
    marker = Path(PRECONVERTED) / run_acc / "done.json"
    record = json.loads(marker.read_text()) if PRECONVERTED and marker.exists() else None
    files = [marker.parent / name for name in sorted(record["files"])] if record else []
    if record and all(f.exists() for f in files):
        if record.get("run") != run_acc or any(f.stat().st_size != record["files"][f.name]["bytes"] for f in files):
            raise RuntimeError(f"{run_acc}: preconverted FASTQ differs from its done.json")
        log(f"{run_acc}: using preconverted gzipped FASTQ ({record['tool']})")
        return [str(f) for f in files]
    if record:                      # used (and removed) by the first pass; the archive serves the later pass
        log(f"{run_acc}: preconverted FASTQ already used; reading the archive again")
    url = SRA_ODP.format(run=run_acc)
    size = _head_size(url)
    ena_bytes = sum(int(b) for b in r.get("fastq_bytes", "").split(";") if b.isdigit())
    ena_need = ena_bytes / 1e9 + SCRATCH_MARGIN_GB if r.get("fastq_ftp") else None
    if size:
        # Reads are streamed out of the archive straight into gzip (fastq-dump | split | pigz): about 1.2-1.5x
        # the archive on local scratch instead of fasterq-dump's ~8x uncompressed FASTQ plus temporary files,
        # and a 0.8 GB run takes about a minute (2026-11-03). Runs too large even for this are converted on
        # NAS scratch; ENA's gzip copy is used only for runs without an SRA archive (it stalled from here).
        need = size * 1.6 / 1e9 + SCRATCH_MARGIN_GB
        if NAS_SCRATCH and SCRATCH == NAS_SCRATCH:
            return _gzip_on_nas(r, url, size)
        local_free = _free_gb(SCRATCH) - _reserved_by_others(run_acc)
        if local_free < need and NAS_SCRATCH and _nas_scratch_fits(size * 3 / 1e9 + SCRATCH_MARGIN_GB):
            return _gzip_on_nas(r, url, size)
        _need_scratch(run_acc, need)
        cached = (_cache_module().cached_archive(Path(DOWNLOAD_CACHE), run_acc, url, size)
                  if DOWNLOAD_CACHE else None)
        sra = str(cached) if cached else os.path.join(SCRATCH, f"{run_acc}.sra")
        if cached:
            log(f"{run_acc}: using verified prefetched archive")
        else:
            fetch(url, sra)
        run([os.path.join(BIN, "vdb-validate"), sra], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        files = _stream_gzip(run_acc, sra)
        if files is None:                      # paired run with unpaired reads: keep them apart, as --split-3 does
            if not NAS_SCRATCH:
                raise RuntimeError(f"{run_acc}: unpaired reads in a paired run and no NAS scratch for --split-3")
            return _gzip_on_nas(r, url, size)
        os.remove(sra)
        if cached:
            cached.with_suffix('.json').unlink(missing_ok=True)
        return files
    # Fallback: ENA's gzipped FASTQ, md5-checked.
    urls = [u for u in r.get("fastq_ftp", "").split(";") if u]
    md5s = [m for m in r.get("fastq_md5", "").split(";") if m]
    sizes = [int(b) for b in r.get("fastq_bytes", "").split(";") if b.isdigit()]
    if not urls:
        return []
    need = sum(sizes) / 1e9 + SCRATCH_MARGIN_GB
    _need_scratch(run_acc, need)
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


def _stream_gzip(run_acc: str, sra: str):
    """Gzipped FASTQ on local scratch, streamed from the archive; None when R1/R2 counts differ."""
    r1, r2 = (os.path.join(SCRATCH, f"{run_acc}_{k}.fastq.gz") for k in (1, 2))
    counts = os.path.join(SCRATCH, f"{run_acc}.counts")
    pigz, threads = os.path.join(BIN, "pigz"), str(max(2, THREADS // 2))
    # Each read's header ends in /<read index> ($ri), which routes it; awk counts both sides.
    split = ('NR % 4 == 1 { n = split($0, a, "/"); side = (a[n] == "1") ? 1 : 2; seen[side]++ } '
             '{ print | (side == 1 ? c1 : c2) } END { print seen[1] + 0, seen[2] + 0 > counts }')
    pipeline = (f"set -o pipefail; {os.path.join(BIN, 'fastq-dump')} --split-spot --skip-technical "
                f"--defline-seq '@$ac.$si/$ri' --defline-qual '+' -Z {shlex.quote(sra)} 2>/dev/null | "
                f"awk -v c1={shlex.quote(f'{pigz} -p {threads} > {r1}')} -v c2={shlex.quote(f'{pigz} -p {threads} > {r2}')} "
                f"-v counts={shlex.quote(counts)} {shlex.quote(split)}")
    log(f"{run_acc}: streaming reads from the archive into gzip on local scratch")
    subprocess.run(["bash", "-c", pipeline], check=True)
    one, two = (int(x) for x in open(counts).read().split())
    os.remove(counts)
    if one == 0:
        raise RuntimeError(f"{run_acc}: no reads written")
    if two == 0:
        Path(r2).unlink(missing_ok=True)
        return [r1]
    if one != two:
        for f in (r1, r2):
            Path(f).unlink(missing_ok=True)
        log(f"{run_acc}: {one} read-1 vs {two} read-2 records; unpaired reads present")
        return None
    return [r1, r2]


def _gzip_on_nas(r: dict, url: str, size: int) -> list:
    """Gzipped FASTQ for one run, converted from its (cached) SRA archive into its own NAS scratch folder."""
    run_acc = r["run_accession"]
    folder = Path(NAS_SCRATCH) / run_acc
    shutil.rmtree(folder, ignore_errors=True)
    folder.mkdir(parents=True)
    cached = (_cache_module().cached_archive(Path(DOWNLOAD_CACHE), run_acc, url, size)
              if DOWNLOAD_CACHE else None)
    sra = str(cached) if cached else str(folder / f"{run_acc}.sra")
    if not cached:
        fetch(url, sra)
    log(f"{run_acc}: local scratch too small; converting the SRA archive to gzipped FASTQ on NAS scratch")
    run([os.path.join(BIN, "vdb-validate"), sra], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    run([os.path.join(BIN, "fastq-dump"), "--split-3", "--gzip", "-O", str(folder), sra], stdout=subprocess.DEVNULL)
    if not cached:
        os.remove(sra)
    pair = [folder / f"{run_acc}_{k}.fastq.gz" for k in (1, 2)]
    single = folder / f"{run_acc}.fastq.gz"
    if all(p.exists() for p in pair):
        single.unlink(missing_ok=True)       # unpaired leftovers of a paired run
        return [str(p) for p in pair]
    return [str(single)] if single.exists() else []


def _cache_module():
    """Import the optional transport helper for both direct and module execution."""
    try:
        from . import omics_download_cache
    except ImportError:
        import omics_download_cache
    return omics_download_cache

def _clear_scratch(run_accession: str | None) -> None:
    """Remove one run's files (reads, .sra, fasterq-dump temp) from scratch; None clears every run's.

    This is the only quantifier using SCRATCH, so files of other runs are leftovers of interrupted
    attempts (89 GB had accumulated by 2026-11-01).
    """
    for f in os.listdir(SCRATCH):
        stale_run = run_accession is None and re.match(r"[SED]RR\d+", f)
        if stale_run or (run_accession and f.startswith(run_accession)) or f.startswith("fasterq.tmp"):
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
    h = sub.add_parser("host-references", help="build further single-host competitive references")
    h.add_argument("names", nargs="*", help="default: every name in HOST_REFERENCES except ME49_host")
    q = sub.add_parser("quantify")
    q.add_argument("--limit", type=int, default=None)
    rq = sub.add_parser("requantify", help="re-quantify runs verified before trimming (instruction 81)")
    rq.add_argument("--limit", type=int, default=None)
    rq.add_argument("--now", action="store_true", help="do not wait for the new downloads")
    args = ap.parse_args(argv)
    if args.cmd == "references":
        references()
        host_reference()
    elif args.cmd == "host-references":
        for name in args.names or [n for n in HOST_REFERENCES if n != "ME49_host"]:
            host_reference(name)
    elif args.cmd == "runs":
        plan(args.studies, set(args.strain))
    elif args.cmd == "requantify":
        requantify(args.limit, wait=not args.now)
    else:
        quantify(args.limit)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
