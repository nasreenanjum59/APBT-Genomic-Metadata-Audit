#!/usr/bin/env python3
"""
01_APBT_BioSample_Metadata_Retrieval.py
=======================================

APBT Study 1: Metadata retrieval from NCBI BioSample

Retrieves SARS-CoV-2 BioSample metadata from NCBI using the Entrez API
and saves the parsed metadata as a structured CSV file for the APBT
metadata exposure audit.

Public reproducibility version
------------------------------
This version removes machine-specific paths and personal contact details
from the source code. The NCBI contact email is supplied at runtime through
the --email option or the NCBI_EMAIL environment variable.

Example
-------
python 01_APBT_BioSample_Metadata_Retrieval.py \
    --email your_email@example.com \
    --n-records 100000 \
    --output data/SARS_CoV2_BioSample_metadata_100000.csv

Requirements
------------
pip install biopython
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import time
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

try:
    from Bio import Entrez
except ImportError:
    print("\nERROR: Biopython is not installed.")
    print("Install it with: pip install biopython")
    sys.exit(1)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

ENTREZ_TOOL_NAME = "APBT_Genomic_Metadata_Audit"
DEFAULT_BATCH_SIZE = 200
DEFAULT_OUTPUT_DIR = Path(
    os.getenv("APBT_OUTPUT_DIR", str(Path.cwd() / "data"))
)

# NCBI taxonomy ID 2697049:
# Severe acute respiratory syndrome coronavirus 2
BIOSAMPLE_QUERY_STRATEGIES = [
    'txid2697049[Organism]',
    'txid2697049[Organism] AND biosample[filter]',
    '"Severe acute respiratory syndrome coronavirus 2"[Organism]',
    'SARS-CoV-2[All Fields] AND (biosample[filter] OR "BioSample"[All Fields])',
    'txid2697049[Organism] AND "Homo sapiens"[Host]',
    'COVID-19[All Fields] AND coronavirus[All Fields]',
]

SRA_FALLBACK_QUERIES = [
    'txid2697049[Organism]',
    '"SARS-CoV-2"[Organism]',
    'COVID-19[All Fields]',
]

CSV_COLUMNS = [
    # Provenance / administrative
    "biosample_accession",
    "submission_date",
    "last_update",
    "submitter_handle",
    "bioproject",
    "sra_accession",
    "genbank_accession",

    # Biological
    "organism",
    "taxonomy_id",

    # Temporal
    "collection_date",

    # Geographic
    "geographic_location",
    "country",
    "lat_lon",

    # Host / specimen
    "host",
    "isolation_source",
    "tissue",
    "sample_type",

    # Technical
    "sequencing_platform",
    "sequencing_technology",
    "assembly_method",
    "coverage_depth",
    "quality_check",

    # Analytical
    "pango_lineage",
    "clade",
    "variant",

    # Laboratory / collection context
    "lab_name",
    "collected_by",
    "sequenced_by",
    "submitting_lab",

    # Reserved audit columns
    "audit_missing_date",
    "audit_date_anomaly",
    "audit_missing_geo",
    "audit_geo_inconsistent",
    "audit_host_mismatch",
    "audit_notes",
]


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def safe_text(element: Optional[ET.Element], default: str = "") -> str:
    """Return stripped element text, or a default value when absent."""
    if element is None:
        return default
    return (element.text or "").strip()


def configure_entrez(email: str, api_key: Optional[str] = None) -> None:
    """Configure Biopython Entrez without hard-coding personal credentials."""
    Entrez.email = email
    Entrez.tool = ENTREZ_TOOL_NAME

    if api_key:
        Entrez.api_key = api_key


# ---------------------------------------------------------------------------
# BioSample XML parser
# ---------------------------------------------------------------------------

def parse_biosample(xml_text: str) -> Dict[str, str]:
    """Parse one BioSample XML record into the audit CSV schema."""
    record = {column: "" for column in CSV_COLUMNS}

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return record

    # Top-level BioSample attributes
    record["biosample_accession"] = root.get("accession", "")
    record["submission_date"] = root.get("submission_date", "")
    record["last_update"] = root.get("last_update", "")

    # Organism block
    organism = root.find(".//Organism")
    if organism is not None:
        record["organism"] = organism.get("taxonomy_name", "")
        record["taxonomy_id"] = organism.get("taxonomy_id", "")

    # Owner / submitter
    owner = root.find(".//Owner")
    if owner is not None:
        record["submitter_handle"] = safe_text(owner.find("Name"))
        full_name = owner.find(".//Contact/Name/Full")
        if full_name is not None:
            record["lab_name"] = safe_text(full_name)

    # BioSample attributes
    for attribute in root.findall(".//Attribute"):
        name = (attribute.get("attribute_name") or "").lower().strip()
        harmonized_name = (
            attribute.get("harmonized_name") or ""
        ).lower().strip()
        value = safe_text(attribute)

        for key in (name, harmonized_name):
            if not key or not value:
                continue

            if key in (
                "collection_date",
                "collection date",
                "collection_date_precision",
            ):
                if not record["collection_date"]:
                    record["collection_date"] = value

            elif key in (
                "geo_loc_name",
                "geographic location",
                "geographic_location",
                "geo_loc_name_country_only",
                "country",
                "geographic location (country and/or sea)",
            ):
                if not record["geographic_location"]:
                    record["geographic_location"] = value

                    if ":" in value:
                        record["country"] = value.split(":", 1)[0].strip()
                    elif "/" in value:
                        record["country"] = value.split("/", 1)[0].strip()
                    else:
                        record["country"] = value

            elif key in (
                "lat_lon",
                "latitude and longitude",
                "latitude_and_longitude",
            ):
                record["lat_lon"] = value

            elif key in (
                "host",
                "host scientific name",
                "host organism",
                "specific host",
            ):
                if not record["host"]:
                    record["host"] = value

            elif key in ("isolation_source", "isolation source"):
                record["isolation_source"] = value

            elif key in ("tissue", "tissue type", "host tissue sampled"):
                record["tissue"] = value

            elif key in (
                "sample_type",
                "sample type",
                "biomaterial_provider",
            ):
                record["sample_type"] = value

            elif key in (
                "instrument_model",
                "sequencing instrument",
                "instrument model",
            ):
                record["sequencing_platform"] = value

            elif key in (
                "platform",
                "sequencing technology",
                "sequencing_technology",
                "library_strategy",
            ):
                record["sequencing_technology"] = value

            elif key in (
                "assembly_method",
                "assembly method",
                "de_novo_single_contig",
                "genome_assembly_method",
            ):
                record["assembly_method"] = value

            elif key in (
                "coverage",
                "depth of coverage",
                "sequencing depth",
                "genome_coverage",
            ):
                record["coverage_depth"] = value

            elif key in (
                "seq_quality_check",
                "quality check",
                "sequence_quality",
            ):
                record["quality_check"] = value

            elif key in (
                "lineage",
                "pango_lineage",
                "pango lineage",
                "pangolin_lineage",
            ):
                record["pango_lineage"] = value

            elif key in (
                "clade",
                "nextstrain_clade",
                "nextstrain clade",
            ):
                record["clade"] = value

            elif key in (
                "variant",
                "who_label",
                "who label",
                "variant_designation",
            ):
                record["variant"] = value

            elif key in ("collected_by", "collected by"):
                record["collected_by"] = value

            elif key in ("sequenced_by", "sequenced by"):
                record["sequenced_by"] = value

            elif key in (
                "submitting_lab",
                "submitting lab",
                "originating_lab",
                "originating lab",
            ):
                if not record["submitting_lab"]:
                    record["submitting_lab"] = value

            elif key in ("bioproject_accession", "bioproject"):
                record["bioproject"] = value

    # Linked identifiers
    for identifier in root.findall(".//Id"):
        database = (identifier.get("db") or "").upper()
        value = safe_text(identifier)

        if database == "SRA" and not record["sra_accession"]:
            record["sra_accession"] = value
        elif database in ("GENBANK", "INSDC") and not record["genbank_accession"]:
            record["genbank_accession"] = value

    # Links block
    for link in root.findall(".//Link"):
        label = (link.get("label") or "").upper()
        target = (link.get("target") or "").upper()

        if ("SRA" in label or "SRA" in target) and not record["sra_accession"]:
            record["sra_accession"] = link.get("value", "")

    return record


# ---------------------------------------------------------------------------
# NCBI search
# ---------------------------------------------------------------------------

def search_with_fallback(
    n_records: int,
    database: str = "biosample",
) -> Tuple[List[str], Optional[str], Optional[str], int, Optional[str]]:
    """
    Try BioSample query strategies in order until one returns results.

    Returns:
        ids, WebEnv, QueryKey, total_available, winning_query
    """
    for index, query in enumerate(BIOSAMPLE_QUERY_STRATEGIES, start=1):
        print(f"  Strategy {index}: {query}")

        try:
            with Entrez.esearch(
                db=database,
                term=query,
                retmax=n_records,
                usehistory="y",
            ) as handle:
                result = Entrez.read(handle)

            total = int(result["Count"])
            ids = list(result["IdList"])

            print(f"           -> {total:,} records found", end="")

            if total > 0:
                print("  [using this query]\n")
                return (
                    ids,
                    result["WebEnv"],
                    result["QueryKey"],
                    total,
                    query,
                )

            print("  (0 results, trying next...)")

        except Exception as exc:
            print(f"           -> Error: {exc}  (trying next...)")

        time.sleep(0.5)

    # Preserve the original script's SRA fallback behavior.
    print("\n  BioSample queries returned 0 results.")
    print("  Trying SRA as a fallback source...")

    for query in SRA_FALLBACK_QUERIES:
        try:
            with Entrez.esearch(
                db="sra",
                term=query,
                retmax=n_records,
                usehistory="y",
            ) as handle:
                result = Entrez.read(handle)

            total = int(result["Count"])

            if total > 0:
                print(f"  SRA: {total:,} records found with: {query}\n")
                return (
                    list(result["IdList"]),
                    result["WebEnv"],
                    result["QueryKey"],
                    total,
                    f"SRA:{query}",
                )

        except Exception as exc:
            print(f"  SRA error: {exc}")

        time.sleep(0.5)

    return [], None, None, 0, None


# ---------------------------------------------------------------------------
# Record retrieval
# ---------------------------------------------------------------------------

def fetch_records(
    ids: List[str],
    web_env: str,
    query_key: str,
    database: str = "biosample",
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> Iterable[Dict[str, str]]:
    """Fetch XML records in batches and yield parsed dictionaries."""
    total = len(ids)
    fetched = 0

    for start in range(0, total, batch_size):
        end = min(start + batch_size, total)

        print(
            f"  Fetching {start + 1:,}-{end:,} of {total:,}...",
            end=" ",
            flush=True,
        )

        raw = None

        for attempt in range(2):
            try:
                with Entrez.efetch(
                    db=database,
                    rettype="xml",
                    retmode="xml",
                    retstart=start,
                    retmax=batch_size,
                    webenv=web_env,
                    query_key=query_key,
                ) as handle:
                    raw = handle.read()
                break

            except Exception as exc:
                if attempt == 0:
                    print(f"fetch failed ({exc}); retrying in 5 seconds...")
                    time.sleep(5)
                else:
                    print(f"skipping batch after retry failure: {exc}")

        if raw is None:
            continue

        try:
            raw_text = (
                raw
                if isinstance(raw, str)
                else raw.decode("utf-8", errors="replace")
            )

            try:
                root = ET.fromstring(raw_text)
            except ET.ParseError:
                root = ET.fromstring(f"<root>{raw_text}</root>")

            samples = root.findall(".//BioSample")

            if not samples and root.tag == "BioSample":
                samples = [root]

            # Fallback parsing for SRA SAMPLE records.
            if not samples:
                samples = root.findall(".//SAMPLE")

            parsed = 0

            for sample in samples:
                xml_text = ET.tostring(sample, encoding="unicode")
                record = parse_biosample(xml_text)

                # Basic SRA fallback extraction if no BioSample accession
                # was parsed from the element.
                if not record["biosample_accession"]:
                    accession = sample.get("accession", "")

                    if accession:
                        record["biosample_accession"] = accession
                        record["organism"] = safe_text(
                            sample.find(".//TAXON/SCIENTIFIC_NAME")
                        )
                        record["taxonomy_id"] = safe_text(
                            sample.find(".//TAXON/TAXON_ID")
                        )

                        for attribute in sample.findall(".//SAMPLE_ATTRIBUTE"):
                            tag = safe_text(
                                attribute.find("TAG")
                            ).lower()
                            value = safe_text(attribute.find("VALUE"))

                            if "collection" in tag:
                                record["collection_date"] = value
                            elif (
                                "geo" in tag
                                or "country" in tag
                                or "location" in tag
                            ):
                                record["geographic_location"] = value
                            elif "host" in tag:
                                record["host"] = value
                            elif "isolation" in tag:
                                record["isolation_source"] = value

                if record["biosample_accession"]:
                    yield record
                    parsed += 1

            fetched += parsed
            print(f"parsed {parsed:,}  (total: {fetched:,})")

        except Exception as exc:
            print(f"\n  Parse error in batch: {exc}")

        # Keep requests within the standard NCBI rate limit when no API key
        # is supplied.
        time.sleep(0.4)


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def save_csv(records: Iterable[Dict[str, str]], filepath: Path) -> None:
    """Write parsed records to CSV."""
    filepath.parent.mkdir(parents=True, exist_ok=True)

    with filepath.open(
        "w",
        newline="",
        encoding="utf-8-sig",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=CSV_COLUMNS,
            extrasaction="ignore",
        )
        writer.writeheader()

        for record in records:
            writer.writerow(record)


def default_output_path(n_records: int) -> Path:
    """Return a portable default output path."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    filename = (
        f"SARS_CoV2_BioSample_metadata_"
        f"{n_records}records_{timestamp}.csv"
    )

    return DEFAULT_OUTPUT_DIR / filename


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Retrieve SARS-CoV-2 BioSample metadata from NCBI "
            "for the APBT metadata exposure audit."
        )
    )

    parser.add_argument(
        "--n-records",
        type=int,
        default=None,
        help="Number of records to retrieve.",
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output CSV path.",
    )

    parser.add_argument(
        "--email",
        default=os.getenv("NCBI_EMAIL"),
        help=(
            "Contact email required by NCBI. "
            "Can also be supplied through the NCBI_EMAIL environment variable."
        ),
    )

    parser.add_argument(
        "--api-key",
        default=os.getenv("NCBI_API_KEY"),
        help=(
            "Optional NCBI API key. "
            "Can also be supplied through the NCBI_API_KEY environment variable."
        ),
    )

    parser.add_argument(
        "--batch-size",
        type=int,
        default=DEFAULT_BATCH_SIZE,
        help=f"Entrez fetch batch size (default: {DEFAULT_BATCH_SIZE}).",
    )

    return parser.parse_args()


def ask_for_record_count() -> int:
    """Prompt for the number of records when not provided on the command line."""
    print("=" * 68)
    print("  APBT Study 1 - NCBI BioSample Metadata Retrieval")
    print("=" * 68)
    print()

    while True:
        try:
            value = int(
                input("  Number of records to retrieve [100000]: ").strip()
                or "100000"
            )

            if value < 1:
                print("  Please enter a number greater than zero.")
                continue

            return value

        except ValueError:
            print("  Please enter a whole number.")


def resolve_email(command_line_email: Optional[str]) -> str:
    """Obtain the NCBI contact email without storing it in source code."""
    email = (command_line_email or "").strip()

    if not email:
        print()
        email = input(
            "  Enter the contact email to use for NCBI Entrez: "
        ).strip()

    if not email or "@" not in email:
        print("\nERROR: A valid contact email is required by NCBI Entrez.")
        print(
            "Provide it with --email or the NCBI_EMAIL environment variable."
        )
        sys.exit(1)

    return email


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_arguments()

    n_records = (
        args.n_records
        if args.n_records is not None
        else ask_for_record_count()
    )

    if n_records < 1:
        print("ERROR: --n-records must be greater than zero.")
        sys.exit(1)

    if args.batch_size < 1:
        print("ERROR: --batch-size must be greater than zero.")
        sys.exit(1)

    email = resolve_email(args.email)
    configure_entrez(email=email, api_key=args.api_key)

    output_path = (
        args.output
        if args.output is not None
        else default_output_path(n_records)
    )

    output_path = output_path.expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print()
    print("=" * 68)
    print("  APBT NCBI METADATA RETRIEVAL")
    print("=" * 68)
    print(f"  Requested records : {n_records:,}")
    print(f"  Output file       : {output_path}")
    print(f"  Batch size        : {args.batch_size}")
    print(
        f"  NCBI API key      : "
        f"{'provided' if args.api_key else 'not provided'}"
    )
    print("=" * 68)
    print()

    start_time = time.time()

    try:
        (
            ids,
            web_env,
            query_key,
            total_available,
            winning_query,
        ) = search_with_fallback(n_records)

    except Exception as exc:
        print(f"\nERROR connecting to NCBI: {exc}")
        print("Check the internet connection and try again.")
        sys.exit(1)

    if not ids or web_env is None or query_key is None or winning_query is None:
        print("\nNo records were returned by the configured query strategies.")
        sys.exit(1)

    database = "sra" if winning_query.startswith("SRA:") else "biosample"
    actual_query = winning_query.replace("SRA:", "", 1)

    print(f"  Query used        : {actual_query}")
    print(f"  Database          : {database.upper()}")
    print(f"  Records available : {total_available:,}")
    print(f"  Records requested : {min(n_records, total_available):,}")
    print()

    records = list(
        fetch_records(
            ids=ids,
            web_env=web_env,
            query_key=query_key,
            database=database,
            batch_size=args.batch_size,
        )
    )

    if not records:
        print("\nNo records could be parsed from the XML response.")
        sys.exit(1)

    print(f"\nSaving {len(records):,} parsed records...")
    save_csv(records, output_path)

    elapsed = time.time() - start_time
    size_kb = output_path.stat().st_size // 1024

    print()
    print("=" * 68)
    print("  RETRIEVAL COMPLETE")
    print("=" * 68)
    print(f"  Records saved : {len(records):,}")
    print(f"  File          : {output_path}")
    print(f"  File size     : {size_kb:,} KB")
    print(f"  Time taken    : {elapsed:.0f} seconds")
    print()

    # Field-completeness preview
    with output_path.open(
        newline="",
        encoding="utf-8-sig",
    ) as handle:
        rows = list(csv.DictReader(handle))

    preview_fields = [
        "collection_date",
        "geographic_location",
        "host",
        "pango_lineage",
        "submission_date",
        "submitter_handle",
        "isolation_source",
        "sequencing_platform",
    ]

    print("  Field completeness preview:")
    print(
        f"  {'Field':<28} "
        f"{'Present':>9} "
        f"{'Missing':>9} "
        f"{'Complete':>10}"
    )
    print("  " + "-" * 60)

    for field in preview_fields:
        present = sum(
            1
            for row in rows
            if row.get(field, "").strip()
        )

        missing = len(rows) - present
        percentage = (
            present / len(rows) * 100
            if rows
            else 0
        )

        print(
            f"  {field:<28} "
            f"{present:>9,} "
            f"{missing:>9,} "
            f"{percentage:>9.1f}%"
        )

    print()
    print("  Continue with the next APBT audit program.")
    print("=" * 68)


if __name__ == "__main__":
    main()
