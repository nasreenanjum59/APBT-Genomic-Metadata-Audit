#!/usr/bin/env python3
"""
03_APBT_CrossRepository_Linkage_Enrichment.py
=============================================

APBT Study 1: Cross-repository linkage and dataset enrichment

Public reproducibility version of the cross-repository enrichment stage used
in the APBT metadata audit.

Input
-----
A CSV produced by Program 01 containing a ``biosample_accession`` column.

Cross-repository checks enriched by this program
------------------------------------------------
P2  BioProject association
P3  Nucleotide/GenBank linkage
P4  SRA linkage
A1  Sequencing platform or technology
A3  Coverage or sequencing depth
A4  Assembly method

The program preserves the original BioSample cohort and follows public NCBI
Entrez relationships to BioProject, SRA, and Nucleotide/GenBank resources.
Information present directly in BioSample is retained as available. Information
recovered from linked NCBI records is recorded separately. When information
cannot be recovered from the queried resources, the status is recorded as
``not_found_in_queried_resources``.

That status refers only to the resources queried by this program and must not
be interpreted as proof that the information is absent from all NCBI or INSDC
resources.

For A3 and A4, this implementation examines BioSample-local values and linked
Nucleotide/GenBank metadata, including Assembly-Data structured comments when
available. It does not require a separate NCBI Assembly-database ELink pass.

The program is read-only with respect to NCBI and uses a SQLite cache so large
runs can be resumed without repeating completed requests.

Usage
-----
Pilot run:

    python 03_APBT_CrossRepository_Linkage_Enrichment.py \
        --input data/SARS_CoV2_BioSample_metadata_100000.csv \
        --email your_email@example.com \
        --limit 100

Full run:

    python 03_APBT_CrossRepository_Linkage_Enrichment.py \
        --input data/SARS_CoV2_BioSample_metadata_100000.csv \
        --email your_email@example.com

Requirements
------------
pip install biopython pandas
"""

import argparse
import json
import os
import re
import sqlite3
import time
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import Entrez


TOOL = "APBT_CrossRepository_Audit"

# ---------------------------------------------------------------------------
# Portable public configuration
# ---------------------------------------------------------------------------

DEFAULT_ENTREZ_EMAIL = os.environ.get("NCBI_EMAIL", "")
DEFAULT_RESULTS_DIR = Path.cwd() / "results" / "cross_repository"

# 0 means analyse every record in the input CSV.
DEFAULT_RECORD_LIMIT = 0

RESOLVE_BATCH = 40
ELINK_BATCH = 25
SUMMARY_BATCH = 100
NUCCORE_FETCH_BATCH = 50

BIOPROJECT_RE = re.compile(r"\bPRJ(?:NA|EB|DB)\d+\b", re.I)
SRA_RUN_RE = re.compile(r"\b(?:SRR|ERR|DRR)\d+\b", re.I)
SRA_EXP_RE = re.compile(r"\b(?:SRX|ERX|DRX)\d+\b", re.I)
SRA_STUDY_RE = re.compile(r"\b(?:SRP|ERP|DRP)\d+\b", re.I)
ASSEMBLY_RE = re.compile(r"\bGC[AF]_\d+\.\d+\b", re.I)
PLACEHOLDERS = {"", "na", "n/a", "none", "null", "unknown",
                "not provided", "not available", "not applicable", "missing"}


def clean(x):
    if x is None:
        return ""
    s = str(x).strip()
    return "" if s.lower() in PLACEHOLDERS else s


def unique(values):
    out, seen = [], set()
    for v in values:
        v = clean(v)
        if v and v not in seen:
            seen.add(v)
            out.append(v)
    return out


def joined(values):
    return ";".join(unique(values))


def batches(items, size):
    for i in range(0, len(items), size):
        yield items[i:i + size]


def all_text(obj):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield str(k)
            yield from all_text(v)
    elif isinstance(obj, (list, tuple)):
        for v in obj:
            yield from all_text(v)
    elif obj is not None:
        yield str(obj)


def read_entrez(handle):
    try:
        return Entrez.read(handle)
    finally:
        handle.close()


def retry(fn, label, attempts=5):
    last = None
    for n in range(1, attempts + 1):
        try:
            return fn()
        except Exception as e:
            last = e
            wait = 3 * n
            print(f"    {label} failed: {e}; retrying in {wait}s")
            time.sleep(wait)
    raise RuntimeError(f"{label} failed after {attempts} attempts: {last}")


class Cache:
    def __init__(self, filename):
        self.db = sqlite3.connect(filename)
        self.db.execute(
            "CREATE TABLE IF NOT EXISTS uid "
            "(accession TEXT PRIMARY KEY, uid TEXT)"
        )
        self.db.execute(
            "CREATE TABLE IF NOT EXISTS links "
            "(accession TEXT, dest TEXT, value TEXT, "
            " PRIMARY KEY(accession,dest))"
        )
        self.db.execute(
            "CREATE TABLE IF NOT EXISTS meta "
            "(db TEXT, uid TEXT, value TEXT, PRIMARY KEY(db,uid))"
        )
        self.db.commit()

    def get_uid(self, acc):
        r = self.db.execute(
            "SELECT uid FROM uid WHERE accession=?", (acc,)
        ).fetchone()
        return r[0] if r else None

    def put_uid(self, acc, uid):
        self.db.execute(
            "INSERT OR REPLACE INTO uid VALUES(?,?)", (acc, uid)
        )
        self.db.commit()

    def get_links(self, acc, dest):
        r = self.db.execute(
            "SELECT value FROM links WHERE accession=? AND dest=?",
            (acc, dest)
        ).fetchone()
        return json.loads(r[0]) if r else None

    def put_links(self, acc, dest, value):
        self.db.execute(
            "INSERT OR REPLACE INTO links VALUES(?,?,?)",
            (acc, dest, json.dumps(value))
        )
        self.db.commit()

    def get_meta(self, db, uid):
        r = self.db.execute(
            "SELECT value FROM meta WHERE db=? AND uid=?", (db, uid)
        ).fetchone()
        return json.loads(r[0]) if r else None

    def put_meta(self, db, uid, value):
        self.db.execute(
            "INSERT OR REPLACE INTO meta VALUES(?,?,?)",
            (db, uid, json.dumps(value))
        )
        self.db.commit()

    def close(self):
        self.db.close()


def resolve_biosample_uids(accessions, cache):
    """Resolve the exact BioSample accessions already present in the study CSV."""
    resolved, todo = {}, []
    for acc in accessions:
        uid = cache.get_uid(acc)
        if uid:
            resolved[acc] = uid
        else:
            todo.append(acc)

    print(f"[1/5] Resolving {len(todo):,} uncached BioSample accessions")

    for i, batch in enumerate(batches(todo, RESOLVE_BATCH), 1):
        term = " OR ".join(f'"{a}"[All Fields]' for a in batch)

        result = retry(
            lambda: read_entrez(Entrez.esearch(
                db="biosample", term=term, retmax=max(100, len(batch) * 3)
            )),
            "BioSample ESearch"
        )
        uids = [str(x) for x in result.get("IdList", [])]

        # Fetch BioSample XML because the BioSample element normally exposes
        # both accession and numeric Entrez ID.
        mapping = {}
        if uids:
            def _fetch():
                h = Entrez.efetch(
                    db="biosample", id=uids, rettype="full", retmode="xml"
                )
                try:
                    raw = h.read()
                finally:
                    h.close()
                if isinstance(raw, bytes):
                    raw = raw.decode("utf-8", errors="replace")
                return raw

            raw = retry(_fetch, "BioSample EFetch")
            try:
                root = ET.fromstring(raw)
            except ET.ParseError:
                root = ET.fromstring("<root>" + raw + "</root>")

            nodes = root.findall(".//BioSample")
            if root.tag == "BioSample":
                nodes = [root]

            for node in nodes:
                acc = node.get("accession", "")
                uid = node.get("id", "") or node.get("Id", "")
                if acc and uid:
                    mapping[acc] = uid

        # Fallback for any accession not mapped from XML.
        for acc in batch:
            if acc in mapping:
                continue

            one = retry(
                lambda a=acc: read_entrez(Entrez.esearch(
                    db="biosample",
                    term=f'"{a}"[All Fields]',
                    retmax=5
                )),
                f"BioSample ESearch {acc}"
            )
            ids = [str(x) for x in one.get("IdList", [])]
            if ids:
                mapping[acc] = ids[0]

        for acc, uid in mapping.items():
            resolved[acc] = uid
            cache.put_uid(acc, uid)

        if i == 1 or i % 25 == 0:
            print(f"    processed ~{min(i*RESOLVE_BATCH, len(todo)):,}/{len(todo):,}")

    return resolved


def parse_elink(result):
    out = {}
    for linkset in result:
        source = [str(x) for x in linkset.get("IdList", [])]
        if not source:
            continue
        dest = []
        for block in linkset.get("LinkSetDb", []):
            for link in block.get("Link", []):
                if link.get("Id"):
                    dest.append(str(link["Id"]))
        out[source[0]] = unique(dest)
    return out


def get_links(accessions, uid_map, dest, cache):
    """
    Retrieve BioSample cross-database links robustly.

    Large ELink responses can occasionally be truncated by the remote server.
    Instead of repeatedly retrying the same large request, this implementation
    automatically splits a failing batch into smaller pieces. Successful
    sub-batches are written to the SQLite cache immediately, so a stopped or
    interrupted run can safely resume without repeating completed work.
    """
    out, todo = {}, []

    for acc in accessions:
        old = cache.get_links(acc, dest)
        if old is not None:
            out[acc] = old
        elif acc in uid_map:
            todo.append(acc)
        else:
            out[acc] = []

    print(f"    BioSample -> {dest}: {len(todo):,} uncached records")

    if not todo:
        return out

    completed = 0
    next_report = 1000

    def process_batch(batch):
        nonlocal completed, next_report

        if not batch:
            return

        ids = [uid_map[a] for a in batch]

        def request():
            return read_entrez(Entrez.elink(
                dbfrom="biosample",
                db=dest,
                id=ids
            ))

        try:
            # Two attempts at the current batch size. If both fail, reduce
            # response size rather than repeatedly hammering the same request.
            result = retry(
                request,
                f"ELink BioSample->{dest} batch={len(batch)}",
                attempts=2
            )

        except RuntimeError as exc:
            if len(batch) > 1:
                mid = len(batch) // 2
                left = batch[:mid]
                right = batch[mid:]
                print(
                    f"    ELink {dest}: batch of {len(batch)} still failed; "
                    f"splitting into {len(left)} + {len(right)}"
                )
                process_batch(left)
                process_batch(right)
                return

            # A single BioSample should produce a very small response. Give
            # it additional retries. If it still fails, stop rather than
            # incorrectly recording a network failure as 'no link'.
            print(
                f"    ELink {dest}: single-record request failed. "
                f"Applying extended retries for {batch[0]}."
            )
            result = retry(
                request,
                f"ELink BioSample->{dest} single record {batch[0]}",
                attempts=8
            )

        by_uid = parse_elink(result)

        # Cache this successful batch immediately.
        for acc in batch:
            value = by_uid.get(uid_map[acc], [])
            out[acc] = value
            cache.put_links(acc, dest, value)

        completed += len(batch)

        if completed >= next_report or completed == len(todo):
            print(
                f"      {dest}: processed "
                f"{completed:,}/{len(todo):,} uncached BioSamples"
            )
            while next_report <= completed:
                next_report += 1000

    for batch in batches(todo, ELINK_BATCH):
        process_batch(batch)

    return out


def document_uid(doc):
    for key in ("Id", "uid", "UID"):
        if key in doc and clean(doc[key]):
            return clean(doc[key])
    attrs = getattr(doc, "attributes", None)
    if attrs:
        for key in ("uid", "Id"):
            if key in attrs:
                return str(attrs[key])
    return ""


def esummary(db, ids):
    result = retry(
        lambda: read_entrez(Entrez.esummary(db=db, id=ids, retmode="xml")),
        f"{db} ESummary"
    )
    if isinstance(result, list):
        return [dict(x) for x in result]
    if isinstance(result, dict):
        dss = result.get("DocumentSummarySet", result)
        if isinstance(dss, dict) and "DocumentSummary" in dss:
            return [dict(x) for x in dss["DocumentSummary"]]
    return []


def parse_bioproject(doc):
    text = "\n".join(all_text(doc))
    return {"accessions": unique(BIOPROJECT_RE.findall(text))}


def parse_sra(doc):
    text = "\n".join(all_text(doc))
    exp = str(doc.get("ExpXml", ""))
    runs_blob = str(doc.get("Runs", ""))

    runs = unique(SRA_RUN_RE.findall(runs_blob + "\n" + text))
    experiments = unique(SRA_EXP_RE.findall(exp + "\n" + text))
    studies = unique(SRA_STUDY_RE.findall(exp + "\n" + text))
    projects = unique(BIOPROJECT_RE.findall(exp + "\n" + text))

    platforms, instruments = [], []

    for m in re.finditer(
        r"<Platform[^>]*instrument_model=[\"']([^\"']+)[\"'][^>]*>([^<]+)</Platform>",
        exp, re.I
    ):
        instruments.append(m.group(1))
        platforms.append(m.group(2))

    instruments += re.findall(
        r"instrument_model=[\"']([^\"']+)[\"']", exp, re.I
    )
    platforms += re.findall(
        r"<Platform[^>]*>([^<]+)</Platform>", exp, re.I
    )

    return {
        "runs": unique(runs),
        "experiments": unique(experiments),
        "studies": unique(studies),
        "projects": unique(projects),
        "platforms": unique(platforms),
        "instruments": unique(instruments),
    }


def parse_nuccore_summary(doc):
    acc = []
    for k in ("AccessionVersion", "Caption"):
        if k in doc and clean(doc[k]):
            acc.append(str(doc[k]))
    return {"accessions": unique(acc)}


def extract_assembly_fields(text):
    fields = {"coverage": [], "assembly_method": [], "sequencing_technology": []}
    patterns = {
        "assembly_method": [
            r"Assembly\s+Method\s*(?:::|:|=)\s*([^\n\r<>]+)"
        ],
        "coverage": [
            r"(?:Genome\s+)?Coverage\s*(?:::|:|=)\s*([^\n\r<>]+)"
        ],
        "sequencing_technology": [
            r"Sequencing\s+Technology\s*(?:::|:|=)\s*([^\n\r<>]+)"
        ],
    }
    for key, pats in patterns.items():
        for pat in pats:
            for v in re.findall(pat, text, re.I):
                v = re.split(r"(?:##|;|\t)", v)[0].strip(" :")
                if v:
                    fields[key].append(v)
        fields[key] = unique(fields[key])
    return fields


def parse_assembly(doc):
    text = "\n".join(all_text(doc))
    out = extract_assembly_fields(text)
    out["accessions"] = unique(ASSEMBLY_RE.findall(text))
    return out


def fetch_summary_metadata(db, ids, cache, parser):
    data, todo = {}, []
    for uid in sorted(ids):
        old = cache.get_meta(db, uid)
        if old is not None:
            data[uid] = old
        else:
            todo.append(uid)

    print(f"    {db} summaries: {len(todo):,} uncached linked records")

    for batch_no, batch in enumerate(batches(todo, SUMMARY_BATCH), 1):
        docs = esummary(db, batch)
        mapped = {}
        for doc in docs:
            uid = document_uid(doc)
            if uid:
                mapped[uid] = parser(doc)

        # Safe positional fallback only when record counts match.
        if not mapped and len(docs) == len(batch):
            mapped = {uid: parser(doc) for uid, doc in zip(batch, docs)}

        for uid in batch:
            payload = mapped.get(uid, {})
            data[uid] = payload
            cache.put_meta(db, uid, payload)

        if batch_no % 10 == 0 or batch_no * SUMMARY_BATCH >= len(todo):
            done = min(batch_no * SUMMARY_BATCH, len(todo))
            print(f"      {db} summaries: processed {done:,}/{len(todo):,}")

    return data


def fetch_nuccore_assembly_data(ids, nuccore_summary, cache, full=False):
    """
    Recover GenBank Assembly-Data structured-comment fields.

    By default seq_start=1/seq_stop=1 are used to minimise sequence transfer
    while retaining record-level metadata. --full-genbank requests full records.
    """
    key = "nuccore_full" if full else "nuccore_metadata"
    data, todo = {}, []

    for uid in sorted(ids):
        old = cache.get_meta(key, uid)
        if old is not None:
            data[uid] = old
        else:
            todo.append(uid)

    print(f"    GenBank Assembly-Data: {len(todo):,} uncached linked records")

    # Accession -> UID map
    acc_to_uid = {}
    for uid, meta in nuccore_summary.items():
        for acc in meta.get("accessions", []):
            acc_to_uid[acc] = uid
            acc_to_uid[acc.split(".")[0]] = uid

    for i, batch in enumerate(batches(todo, NUCCORE_FETCH_BATCH), 1):
        def _fetch():
            kwargs = dict(
                db="nuccore", id=batch, rettype="gb", retmode="xml"
            )
            if not full:
                kwargs["seq_start"] = 1
                kwargs["seq_stop"] = 1
            h = Entrez.efetch(**kwargs)
            try:
                raw = h.read()
            finally:
                h.close()
            if isinstance(raw, bytes):
                raw = raw.decode("utf-8", errors="replace")
            return raw

        raw = retry(_fetch, "Nuccore EFetch")
        try:
            root = ET.fromstring(raw)
        except ET.ParseError:
            root = ET.fromstring("<root>" + raw + "</root>")

        records = root.findall(".//GBSeq") + root.findall(".//INSDSeq")
        if root.tag in ("GBSeq", "INSDSeq"):
            records = [root]

        mapped = {}
        for rec in records:
            accession = ""
            for tag in (
                "GBSeq_accession-version", "GBSeq_primary-accession",
                "INSDSeq_accession-version", "INSDSeq_primary-accession"
            ):
                el = rec.find(".//" + tag)
                if el is not None and clean(el.text):
                    accession = clean(el.text)
                    break

            text = "\n".join(t.strip() for t in rec.itertext() if t and t.strip())
            uid = acc_to_uid.get(accession) or acc_to_uid.get(accession.split(".")[0])
            if uid:
                mapped[uid] = extract_assembly_fields(text)

        for uid in batch:
            payload = mapped.get(uid, {})
            data[uid] = payload
            cache.put_meta(key, uid, payload)

        if i % 50 == 0:
            print(f"      processed ~{min(i*NUCCORE_FETCH_BATCH, len(todo)):,}/{len(todo):,}")

    return data


def gather(meta, ids, field):
    values = []
    for uid in ids:
        values.extend(meta.get(uid, {}).get(field, []) or [])
    return unique(values)


def local_or_link_status(local, linked):
    if clean(local):
        return "present_in_biosample"
    if linked:
        return "recovered_from_linked_resource"
    return "not_found_in_queried_resources"


def link_status(local, linked):
    if clean(local):
        return "present_in_biosample"
    if linked:
        return "recovered_from_entrez_link"
    return "not_found_in_queried_resources"


def make_summary(df):
    definitions = [
        ("P2", "BioProject linkage", "P2_crossrepo_status"),
        ("P3", "Nucleotide/GenBank linkage", "P3_crossrepo_status"),
        ("P4", "SRA linkage", "P4_crossrepo_status"),
        ("A1", "Sequencing platform/technology", "A1_crossrepo_status"),
        ("A3", "Reported coverage", "A3_crossrepo_status"),
        ("A4", "Assembly method", "A4_crossrepo_status"),
    ]
    rows = []
    n = len(df)
    for code, desc, col in definitions:
        c = df[col].value_counts().to_dict()
        local = int(c.get("present_in_biosample", 0))
        recovered = int(c.get("recovered_from_linked_resource", 0)
                        + c.get("recovered_from_entrez_link", 0))
        notfound = int(c.get("not_found_in_queried_resources", 0))
        rows.append({
            "check": code,
            "description": desc,
            "n": n,
            "present_in_biosample_n": local,
            "present_in_biosample_pct": round(100 * local / n, 3) if n else 0,
            "recovered_crossrepo_n": recovered,
            "recovered_crossrepo_pct": round(100 * recovered / n, 3) if n else 0,
            "not_found_after_crossrepo_query_n": notfound,
            "not_found_after_crossrepo_query_pct": round(100 * notfound / n, 3) if n else 0,
        })
    return pd.DataFrame(rows)


def parse_args():
    """Parse command-line arguments for the public reproducibility script."""
    parser = argparse.ArgumentParser(
        description=(
            "Enrich the APBT BioSample cohort with linked NCBI "
            "BioProject, SRA, and Nucleotide/GenBank metadata."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="CSV produced by Program 01.",
    )
    parser.add_argument(
        "--output",
        default="",
        help=(
            "Output enriched CSV. If omitted, a file is created under "
            "results/cross_repository/."
        ),
    )
    parser.add_argument(
        "--email",
        default=DEFAULT_ENTREZ_EMAIL,
        help=(
            "Contact email required by NCBI Entrez. Can also be supplied "
            "through the NCBI_EMAIL environment variable."
        ),
    )
    parser.add_argument(
        "--api-key",
        default=os.environ.get("NCBI_API_KEY", ""),
        help=(
            "Optional NCBI API key. Can also be supplied through the "
            "NCBI_API_KEY environment variable."
        ),
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=DEFAULT_RECORD_LIMIT,
        help="Number of rows to analyse; 0 means all rows.",
    )
    parser.add_argument(
        "--cache",
        default="",
        help="Optional SQLite cache path.",
    )
    parser.add_argument(
        "--skip-genbank-details",
        action="store_true",
        help=(
            "Skip detailed GenBank metadata retrieval. This prevents "
            "cross-repository recovery of A3/A4 from GenBank structured "
            "comments and is not recommended for the full study run."
        ),
    )
    parser.add_argument(
        "--full-genbank",
        action="store_true",
        help=(
            "Retrieve full linked GenBank records instead of the minimal "
            "record slice used for metadata extraction."
        ),
    )

    return parser.parse_args()


def prepare_args(args):
    """Validate paths and derive portable output locations."""
    if not args.email:
        args.email = input(
            "Entrez contact email (required by NCBI): "
        ).strip()

    if not args.email or "@" not in args.email:
        raise ValueError(
            "A valid NCBI Entrez contact email is required. "
            "Provide --email or set the NCBI_EMAIL environment variable."
        )

    input_path = Path(args.input).expanduser().resolve()

    if not input_path.exists():
        raise FileNotFoundError(
            f"Input dataset not found: {input_path}"
        )

    if args.limit < 0:
        raise ValueError("--limit must be 0 or a positive integer.")

    if args.output:
        output_path = Path(args.output).expanduser().resolve()
    else:
        results_dir = DEFAULT_RESULTS_DIR.resolve()
        results_dir.mkdir(parents=True, exist_ok=True)
        output_path = (
            results_dir
            / f"{input_path.stem}_cross_repository_audit.csv"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)

    args.input = str(input_path)
    args.output = str(output_path)

    print("=" * 72)
    print("APBT CROSS-REPOSITORY LINKAGE AND ENRICHMENT")
    print("=" * 72)
    print("Read-only retrieval from public NCBI Entrez resources.")
    print(f"Input dataset    : {args.input}")
    print(f"Output dataset   : {args.output}")
    print(
        f"Records to audit : "
        f"{'ALL' if args.limit == 0 else f'{args.limit:,}'}"
    )
    print(
        f"NCBI API key     : "
        f"{'provided' if args.api_key else 'not provided'}"
    )
    print("=" * 72)
    print()

    return args

def main():
    args = prepare_args(parse_args())

    Entrez.email = args.email
    Entrez.tool = TOOL
    if args.api_key:
        Entrez.api_key = args.api_key
    Entrez.max_tries = 5
    Entrez.sleep_between_tries = 10

    inp = Path(args.input)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    cache_file = (
        Path(args.cache).expanduser().resolve()
        if args.cache
        else out.with_name(f"{out.stem}_cache.sqlite")
    )
    cache = Cache(cache_file)

    try:
        df = pd.read_csv(inp, dtype=str, keep_default_na=False)
        if "biosample_accession" not in df.columns:
            raise ValueError("CSV must contain biosample_accession")

        if args.limit > 0:
            df = df.head(args.limit).copy()

        # Keep the same cohort as the original study file.
        accessions = unique(df["biosample_accession"].tolist())

        print("=" * 72)
        print("APBT CROSS-REPOSITORY AUDIT")
        print(f"Records: {len(df):,}")
        print(f"Unique BioSamples: {len(accessions):,}")
        print(f"Cache: {cache_file}")
        print("=" * 72)

        uid_map = resolve_biosample_uids(accessions, cache)
        print(f"Resolved {len(uid_map):,}/{len(accessions):,} BioSamples")

        print("\n[2/5] Following Entrez links")
        bp_links = get_links(accessions, uid_map, "bioproject", cache)
        sra_links = get_links(accessions, uid_map, "sra", cache)
        nuc_links = get_links(accessions, uid_map, "nuccore", cache)

        # The cross-repository method evaluates whether technical/provenance
        # information is recoverable through linked SRA/GenBank resources.
        # SARS-CoV-2 BioSample records in the pilot did not yield linked
        # NCBI Assembly records, while coverage and assembly-method metadata
        # were recoverable from GenBank structured comments. To avoid an
        # unnecessary 100,000-record Assembly ELink pass, Assembly linking is
        # intentionally skipped in this cross-repository run.
        asm_links = {acc: [] for acc in accessions}
        print("    BioSample -> assembly: SKIPPED (A3/A4 use BioSample and linked GenBank metadata)")

        print("\n[3/5] Retrieving linked metadata")
        bp_ids = {x for v in bp_links.values() for x in v}
        sra_ids = {x for v in sra_links.values() for x in v}
        nuc_ids = {x for v in nuc_links.values() for x in v}

        bp_meta = fetch_summary_metadata(
            "bioproject", bp_ids, cache, parse_bioproject
        )
        sra_meta = fetch_summary_metadata(
            "sra", sra_ids, cache, parse_sra
        )
        nuc_meta = fetch_summary_metadata(
            "nuccore", nuc_ids, cache, parse_nuccore_summary
        )

        # No separate Assembly summary metadata are used in this implementation.
        asm_meta = {}

        if args.skip_genbank_details:
            nuc_deep = {}
        else:
            nuc_deep = fetch_nuccore_assembly_data(
                nuc_ids, nuc_meta, cache, full=args.full_genbank
            )

        for col in [
            "bioproject", "sra_accession", "genbank_accession",
            "sequencing_platform", "sequencing_technology",
            "coverage_depth", "assembly_method"
        ]:
            if col not in df.columns:
                df[col] = ""

        print("\n[4/5] Building cross-repository status fields")
        rows = []

        for _, series in df.iterrows():
            r = series.to_dict()
            acc = clean(r["biosample_accession"])

            bpi = bp_links.get(acc, [])
            sri = sra_links.get(acc, [])
            nui = nuc_links.get(acc, [])
            asi = asm_links.get(acc, [])

            bp_acc = unique(
                gather(bp_meta, bpi, "accessions")
                + gather(sra_meta, sri, "projects")
            )

            sra_runs = gather(sra_meta, sri, "runs")
            sra_exps = gather(sra_meta, sri, "experiments")
            sra_studies = gather(sra_meta, sri, "studies")
            sra_platforms = gather(sra_meta, sri, "platforms")
            sra_instruments = gather(sra_meta, sri, "instruments")

            nuc_acc = gather(nuc_meta, nui, "accessions")

            asm_acc = []
            coverage = unique(
                gather(nuc_deep, nui, "coverage")
            )
            assembly_method = unique(
                gather(nuc_deep, nui, "assembly_method")
            )
            seq_tech = unique(
                gather(nuc_deep, nui, "sequencing_technology")
            )
            platform_values = unique(
                sra_platforms + sra_instruments + seq_tech
            )

            local_platform = clean(r.get("sequencing_platform")) or clean(
                r.get("sequencing_technology")
            )

            r["biosample_entrez_uid"] = uid_map.get(acc, "")

            r["linked_bioproject_accessions"] = joined(bp_acc)
            r["linked_sra_run_accessions"] = joined(sra_runs)
            r["linked_sra_experiment_accessions"] = joined(sra_exps)
            r["linked_sra_study_accessions"] = joined(sra_studies)
            r["linked_genbank_nucleotide_accessions"] = joined(nuc_acc)
            r["linked_assembly_accessions"] = joined(asm_acc)

            r["linked_sequencing_platforms"] = joined(sra_platforms)
            r["linked_instrument_models"] = joined(sra_instruments)
            r["linked_sequencing_technology"] = joined(seq_tech)
            r["linked_reported_coverage"] = joined(coverage)
            r["linked_assembly_method"] = joined(assembly_method)

            r["P2_crossrepo_status"] = link_status(
                r.get("bioproject"), bpi or bp_acc
            )
            r["P3_crossrepo_status"] = link_status(
                r.get("genbank_accession"), nui or nuc_acc
            )
            r["P4_crossrepo_status"] = link_status(
                r.get("sra_accession"), sri or sra_runs or sra_exps
            )
            r["A1_crossrepo_status"] = local_or_link_status(
                local_platform, platform_values
            )
            r["A3_crossrepo_status"] = local_or_link_status(
                r.get("coverage_depth"), coverage
            )
            r["A4_crossrepo_status"] = local_or_link_status(
                r.get("assembly_method"), assembly_method
            )

            for check in ("P2", "P3", "P4", "A1", "A3", "A4"):
                r[f"{check}_crossrepo_not_found"] = int(
                    r[f"{check}_crossrepo_status"]
                    == "not_found_in_queried_resources"
                )

            rows.append(r)

        out_df = pd.DataFrame(rows)
        out_df.to_csv(out, index=False, encoding="utf-8-sig")

        print("\n[5/5] Writing summary")
        summary = make_summary(out_df)
        summary_file = out.with_name(out.stem + "_summary.csv")
        summary.to_csv(summary_file, index=False, encoding="utf-8-sig")

        unresolved_file = out.with_name(out.stem + "_unresolved.csv")
        status_cols = [
            "P2_crossrepo_status", "P3_crossrepo_status", "P4_crossrepo_status",
            "A1_crossrepo_status", "A3_crossrepo_status", "A4_crossrepo_status"
        ]
        unresolved = out_df[
            (out_df["biosample_entrez_uid"] == "")
            | (out_df[status_cols] == "not_found_in_queried_resources").any(axis=1)
        ]
        unresolved[
            ["biosample_accession", "biosample_entrez_uid"] + status_cols
        ].to_csv(unresolved_file, index=False, encoding="utf-8-sig")

        print("\nCOMPLETE")
        print(f"Enriched file : {out}")
        print(f"Summary file  : {summary_file}")
        print(f"Unresolved    : {unresolved_file}")
        print(f"Resume cache  : {cache_file}")
        print()
        print(summary.to_string(index=False))
        print()
        print("Use the cross-repository status fields for the final")
        print("audit scoring. Do not describe not_found_in_queried_resources")
        print("as proof of global absence from all INSDC resources.")

    finally:
        cache.close()


if __name__ == "__main__":
    main()
