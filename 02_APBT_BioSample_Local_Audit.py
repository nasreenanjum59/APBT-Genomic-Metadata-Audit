#!/usr/bin/env python3
"""
02_APBT_BioSample_Local_Audit.py
================================

APBT Study 1: BioSample-local metadata audit

This program reproduces the BioSample-local audit stage used in the APBT
study before cross-repository enrichment was applied. It evaluates 23
metadata checks across five layers using fields available directly in the
BioSample-derived CSV.

IMPORTANT
---------
This is a methodological baseline only. It does NOT represent the final
cross-repository results reported in the revised manuscript. The final
analysis is produced by Program 04 after linked NCBI resources have been
examined.

Article
-------
"From data poisoning to advanced persistent biological threats (APBTs):
a taxonomy and large-scale metadata audit of genomic surveillance
infrastructure"

Usage
-----
python 02_APBT_BioSample_Local_Audit.py \
    --input data/SARS_CoV2_BioSample_metadata_100000.csv \
    --output-dir results/biosample_local

Requirements
------------
pip install pandas numpy matplotlib
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# ---------------------------------------------------------------------------
# Fixed study configuration
# ---------------------------------------------------------------------------

COVID_START = pd.Timestamp("2019-11-01")

# Fixed study cutoff for reproducibility.
# The audited retrieval window ended on 15 April 2026.
AUDIT_CUTOFF = pd.Timestamp("2026-04-15")

# Exact placeholder list used by the audit: 16 strings.
PLACEHOLDER_VALS = {
    "not applicable",
    "not provided",
    "na",
    "n/a",
    "none",
    "unknown",
    "missing",
    "not available",
    "undetermined",
    "not collected",
    "not determined",
    "restricted access",
    "missing: not provided",
    "missing: not collected",
    "missing: third party data",
    "not given",
}

# Exact H2 BioSample-local screen.
NON_HUMAN_HOST = {
    "vero",
    "vero-e6",
    "vero-e6-tmprss2",
    "bat",
    "mink",
    "cat",
    "dog",
    "not applicable",
    "animal",
    "mus musculus",
    "hamster",
}

# H5 isolation-source terms.
IMPLAUSIBLE_ISO = {
    "library",
    "plasmid",
    "pathogen bank",
    "cultured virus",
    "cell culture",
    "vero",
    "vero-e6",
    "laboratory",
}


# ---------------------------------------------------------------------------
# Figure palette and style
# ---------------------------------------------------------------------------

COLORS = {
    "temporal": "#E07B54",
    "geographic": "#E8A838",
    "host": "#5A9E6F",
    "provenance": "#4A86C8",
    "technical": "#8B6BB1",
    "none": "#B0B0B0",
    "low": "#A8CBAB",
    "medium": "#F5C97A",
    "high": "#EE9E82",
    "critical": "#C2596A",
    "bar_base": "#6BAED6",
    "line": "#2171B5",
    "grid": "#E8E8E8",
}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "axes.labelsize": 10,
    "axes.labelweight": "bold",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.axisbelow": True,
    "grid.color": COLORS["grid"],
    "grid.linewidth": 0.6,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.facecolor": "white",
    "figure.facecolor": "white",
})


LAYER_ORDER = [
    "Temporal",
    "Geographic",
    "Host",
    "Provenance",
    "Technical",
]

CATEGORY_ORDER = [
    "None",
    "Low",
    "Medium",
    "High",
    "Critical",
]


# ---------------------------------------------------------------------------
# Audit definitions
# ---------------------------------------------------------------------------

CHECK_META = {
    "T1": (
        "Temporal",
        "Missing date",
        "Missing or unparseable collection date",
    ),
    "T2": (
        "Temporal",
        "Year-only precision",
        "Year-only collection-date precision",
    ),
    "T3": (
        "Temporal",
        "Month-year precision",
        "Month-year collection-date precision",
    ),
    "T4": (
        "Temporal",
        "Pre-emergence date",
        "Collection date before 1 November 2019",
    ),
    "T5": (
        "Temporal",
        "Future date",
        "Collection date after the fixed study cutoff",
    ),
    "T6": (
        "Temporal",
        "After submission",
        "Collection date after submission date",
    ),
    "T7": (
        "Temporal",
        "Same-day submission",
        "Collection and submission on the same calendar day",
    ),
    "G1": (
        "Geographic",
        "Missing location",
        "Missing or functionally missing geographic location",
    ),
    "G2": (
        "Geographic",
        "Country mismatch",
        "Country field inconsistent with geographic_location",
    ),
    "G3": (
        "Geographic",
        "Placeholder geography",
        "Geographic field contains placeholder text",
    ),
    "H1": (
        "Host",
        "Missing host",
        "Missing or functionally missing host organism field",
    ),
    "H2": (
        "Host",
        "Non-human host",
        "Non-human or laboratory-associated host organism by exact-match screen",
    ),
    "H3": (
        "Host",
        "Missing isolation source",
        "Missing or functionally missing isolation source",
    ),
    "H4": (
        "Host",
        "Placeholder isolation source",
        "Isolation-source placeholder text not already captured under H3",
    ),
    "H5": (
        "Host",
        "Implausible isolation source",
        "Implausible isolation source for a human respiratory virus",
    ),
    "P1": (
        "Provenance",
        "Missing collector information",
        "Missing information on who collected the sample",
    ),
    "P2": (
        "Provenance",
        "No BioProject",
        "BioProject linkage unavailable in the BioSample-local record",
    ),
    "P3": (
        "Provenance",
        "No GenBank linkage",
        "Nucleotide/GenBank linkage unavailable in the BioSample-local record",
    ),
    "P4": (
        "Provenance",
        "No SRA linkage",
        "SRA linkage unavailable in the BioSample-local record",
    ),
    "A1": (
        "Technical",
        "No platform",
        "Sequencing platform unavailable in the BioSample-local record",
    ),
    "A2": (
        "Technical",
        "No lineage",
        "Missing Pango lineage assignment",
    ),
    "A3": (
        "Technical",
        "No coverage",
        "Coverage or sequencing depth unavailable in the BioSample-local record",
    ),
    "A4": (
        "Technical",
        "No assembly method",
        "Assembly method unavailable in the BioSample-local record",
    ),
}

LAYER_CHECKS = {
    layer: [
        check_id
        for check_id, (check_layer, _, _) in CHECK_META.items()
        if check_layer == layer
    ]
    for layer in LAYER_ORDER
}

ALL_CHECK_IDS = list(CHECK_META.keys())


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def is_blank(value) -> bool:
    """Return True for empty values or one of the predefined placeholders."""
    if pd.isna(value):
        return True

    text = str(value).strip().lower()
    return text == "" or text in PLACEHOLDER_VALS


def parse_date(value) -> Tuple[pd.Timestamp, str]:
    """
    Parse collection-date values while preserving their precision category.
    """
    if is_blank(value):
        return pd.NaT, "missing"

    text = str(value).strip()

    if re.match(r"^\d{4}-\d{2}-\d{2}", text):
        try:
            return pd.Timestamp(text[:10]), "full"
        except Exception:
            return pd.NaT, "unparseable"

    if re.match(r"^\d{4}-\d{2}$", text):
        try:
            return pd.Timestamp(text + "-01"), "month-year"
        except Exception:
            return pd.NaT, "unparseable"

    if re.match(r"^\d{4}$", text):
        try:
            return pd.Timestamp(text + "-01-01"), "year-only"
        except Exception:
            return pd.NaT, "unparseable"

    return pd.NaT, "unparseable"


def extract_country(geographic_location) -> str:
    """Extract a country-like prefix from a geographic-location field."""
    if is_blank(geographic_location):
        return ""

    text = str(geographic_location).strip()

    for separator in (":", "/", ","):
        if separator in text:
            return text.split(separator, 1)[0].strip()

    return text


def pct(number: int, total: int) -> float:
    """Percentage rounded to one decimal place."""
    return round(number / total * 100, 1) if total else 0.0


def ensure_required_columns(df: pd.DataFrame) -> None:
    """Fail clearly if the expected retrieval fields are absent."""
    required = {
        "biosample_accession",
        "submission_date",
        "collection_date",
        "geographic_location",
        "country",
        "host",
        "isolation_source",
        "collected_by",
        "bioproject",
        "genbank_accession",
        "sra_accession",
        "sequencing_platform",
        "pango_lineage",
        "coverage_depth",
        "assembly_method",
        "submitter_handle",
    }

    missing = sorted(required.difference(df.columns))

    if missing:
        raise ValueError(
            "Input CSV is missing required columns: "
            + ", ".join(missing)
        )


# ---------------------------------------------------------------------------
# Audit engine
# ---------------------------------------------------------------------------

def run_audit(df: pd.DataFrame) -> pd.DataFrame:
    """Run the 23-check BioSample-local audit."""
    df = df.copy()
    n_records = len(df)

    print(f"  Parsing dates for {n_records:,} records...")

    parsed = df["collection_date"].apply(parse_date)
    df["_collection_timestamp"] = parsed.apply(lambda item: item[0])
    df["_collection_precision"] = parsed.apply(lambda item: item[1])

    df["_submission_timestamp"] = pd.to_datetime(
        df["submission_date"],
        errors="coerce",
    )

    print(f"  Running {len(ALL_CHECK_IDS)} checks...")

    # Temporal
    df["T1"] = df["_collection_precision"].isin(
        ["missing", "unparseable"]
    )
    df["T2"] = df["_collection_precision"] == "year-only"
    df["T3"] = df["_collection_precision"] == "month-year"

    df["T4"] = (
        df["_collection_timestamp"].notna()
        & (df["_collection_timestamp"] < COVID_START)
    )

    df["T5"] = (
        df["_collection_timestamp"].notna()
        & (df["_collection_timestamp"] > AUDIT_CUTOFF)
    )

    df["T6"] = (
        df["_collection_timestamp"].notna()
        & df["_submission_timestamp"].notna()
        & (df["_collection_timestamp"] > df["_submission_timestamp"])
    )

    df["T7"] = (
        df["_collection_timestamp"].notna()
        & df["_submission_timestamp"].notna()
        & (
            df["_collection_timestamp"].dt.date
            == df["_submission_timestamp"].dt.date
        )
    )

    # Geographic
    df["G1"] = df["geographic_location"].apply(is_blank)

    df["G3"] = df["geographic_location"].apply(
        lambda value: (
            not pd.isna(value)
            and str(value).strip().lower() in PLACEHOLDER_VALS
        )
    )

    def country_mismatch(row) -> bool:
        geographic_location = str(
            row.get("geographic_location", "")
        ).strip()
        country = str(row.get("country", "")).strip()

        if is_blank(geographic_location) or is_blank(country):
            return False

        expected_country = extract_country(
            geographic_location
        ).lower()

        return (
            expected_country != country.lower()
            and not geographic_location.lower().startswith(
                country.lower()
            )
        )

    df["G2"] = df.apply(country_mismatch, axis=1)

    # Host and specimen
    df["H1"] = df["host"].apply(is_blank)

    df["H2"] = df["host"].apply(
        lambda value: (
            not is_blank(value)
            and str(value).strip().lower() in NON_HUMAN_HOST
        )
    )

    df["H3"] = df["isolation_source"].apply(is_blank)

    # Because is_blank already classifies placeholder values as missing,
    # H4 captures only placeholder text not already flagged by H3.
    df["H4"] = df["isolation_source"].apply(
        lambda value: (
            not is_blank(value)
            and str(value).strip().lower() in PLACEHOLDER_VALS
        )
    )

    df["H5"] = df["isolation_source"].apply(
        lambda value: (
            not is_blank(value)
            and any(
                term in str(value).strip().lower()
                for term in IMPLAUSIBLE_ISO
            )
        )
    )

    # Provenance: BioSample-local baseline only
    df["P1"] = df["collected_by"].apply(is_blank)
    df["P2"] = df["bioproject"].apply(is_blank)
    df["P3"] = df["genbank_accession"].apply(is_blank)
    df["P4"] = df["sra_accession"].apply(is_blank)

    # Technical: BioSample-local baseline only
    df["A1"] = df["sequencing_platform"].apply(is_blank)
    df["A2"] = df["pango_lineage"].apply(is_blank)
    df["A3"] = df["coverage_depth"].apply(is_blank)
    df["A4"] = df["assembly_method"].apply(is_blank)

    # Composite BioSample-local exposure score
    df["exposure_score"] = (
        df[ALL_CHECK_IDS].astype(int).sum(axis=1)
    )

    def category(score: int) -> str:
        if score == 0:
            return "None"
        if score <= 2:
            return "Low"
        if score <= 5:
            return "Medium"
        if score <= 9:
            return "High"
        return "Critical"

    df["exposure_category"] = df["exposure_score"].apply(category)

    return df


def build_summary(
    df: pd.DataFrame,
) -> Tuple[Dict[str, Dict[str, object]], int]:
    """Build per-check and per-layer summary statistics."""
    n_records = len(df)
    summary: Dict[str, Dict[str, object]] = {}

    for check_id in ALL_CHECK_IDS:
        layer, label, description = CHECK_META[check_id]
        count = int(df[check_id].sum())

        summary[check_id] = {
            "check_id": check_id,
            "layer": layer,
            "label": label,
            "description": description,
            "n": count,
            "pct": pct(count, n_records),
        }

    for layer in LAYER_ORDER:
        columns = LAYER_CHECKS[layer]
        count = int(df[columns].any(axis=1).sum())

        summary[f"L_{layer}"] = {
            "check_id": f"L_{layer}",
            "layer": layer,
            "label": f"Any {layer}",
            "description": f"At least one {layer} exposure indicator",
            "n": count,
            "pct": pct(count, n_records),
        }

    count = int((df["exposure_score"] > 0).sum())

    summary["OVERALL"] = {
        "check_id": "OVERALL",
        "layer": "All",
        "label": "Any exposure indicator",
        "description": "At least one audit indicator triggered",
        "n": count,
        "pct": pct(count, n_records),
    }

    return summary, n_records


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def save_figure(fig, path: Path) -> None:
    fig.savefig(
        path,
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)
    print(f"    {path.name}")


def figure_overview(
    df: pd.DataFrame,
    summary: Dict[str, Dict[str, object]],
    n_records: int,
    output_dir: Path,
) -> None:
    """Layer exposure rates and BioSample-local score categories."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    percentages = [
        summary[f"L_{layer}"]["pct"]
        for layer in LAYER_ORDER
    ]

    colours = [
        COLORS[layer.lower()]
        for layer in LAYER_ORDER
    ]

    bars = ax1.barh(
        LAYER_ORDER,
        percentages,
        color=colours,
        alpha=0.82,
        height=0.55,
        edgecolor="white",
        linewidth=0.8,
    )

    for bar, percentage in zip(bars, percentages):
        ax1.text(
            percentage + 0.5,
            bar.get_y() + bar.get_height() / 2,
            f"{percentage:.1f}%",
            va="center",
            fontsize=9,
            fontweight="bold",
        )

    ax1.set_xlabel("Records with >=1 exposure indicator (%)")
    ax1.set_title("BioSample-local Exposure by Metadata Layer")
    ax1.set_xlim(0, max(percentages) * 1.22)
    ax1.invert_yaxis()

    category_counts = {
        category: int(
            (df["exposure_category"] == category).sum()
        )
        for category in CATEGORY_ORDER
    }

    category_colours = [
        COLORS["none"],
        COLORS["low"],
        COLORS["medium"],
        COLORS["high"],
        COLORS["critical"],
    ]

    left = 0.0

    for category, colour in zip(
        CATEGORY_ORDER,
        category_colours,
    ):
        count = category_counts[category]
        percentage = count / n_records * 100

        ax2.barh(
            ["All records"],
            percentage,
            left=left,
            color=colour,
            alpha=0.88,
            height=0.45,
            edgecolor="white",
            linewidth=0.8,
            label=f"{category} ({percentage:.1f}%)",
        )

        if percentage > 4:
            ax2.text(
                left + percentage / 2,
                0,
                f"{percentage:.1f}%",
                ha="center",
                va="center",
                fontsize=8.5,
                fontweight="bold",
            )

        left += percentage

    ax2.set_xlabel("Percentage of records (%)")
    ax2.set_title("BioSample-local Composite Exposure Categories")
    ax2.set_xlim(0, 105)
    ax2.set_yticks([])

    ax2.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, -0.28),
        ncol=3,
        fontsize=8.5,
        frameon=False,
    )

    plt.tight_layout(pad=2.5)

    save_figure(
        fig,
        output_dir / "BioSample_Local_Exposure_Overview.png",
    )


def figure_per_check(
    summary: Dict[str, Dict[str, object]],
    output_dir: Path,
) -> None:
    """Plot all non-zero BioSample-local audit checks."""
    rows = [
        (
            check_id,
            CHECK_META[check_id][0],
            CHECK_META[check_id][1],
            summary[check_id]["pct"],
        )
        for check_id in ALL_CHECK_IDS
        if summary[check_id]["n"] > 0
    ]

    rows.sort(
        key=lambda item: (
            LAYER_ORDER.index(item[1]),
            -item[3],
        )
    )

    labels = [
        f"{item[0]}: {item[2]}"
        for item in rows
    ]

    percentages = [
        item[3]
        for item in rows
    ]

    colours = [
        COLORS[item[1].lower()]
        for item in rows
    ]

    fig, ax = plt.subplots(
        figsize=(10, max(6, len(rows) * 0.42))
    )

    bars = ax.barh(
        labels,
        percentages,
        color=colours,
        alpha=0.80,
        height=0.65,
        edgecolor="white",
        linewidth=0.6,
    )

    for bar, percentage in zip(bars, percentages):
        ax.text(
            bar.get_width() + 0.3,
            bar.get_y() + bar.get_height() / 2,
            f"{percentage:.1f}%",
            va="center",
            fontsize=8,
        )

    ax.set_xlabel("Records with indicator (%)")
    ax.set_title("BioSample-local Exposure Rate per Audit Check")
    ax.set_xlim(
        0,
        max(percentages) * 1.2 if percentages else 100,
    )
    ax.invert_yaxis()

    patches = [
        mpatches.Patch(
            color=COLORS[layer.lower()],
            alpha=0.80,
            label=layer,
        )
        for layer in LAYER_ORDER
    ]

    ax.legend(
        handles=patches,
        loc="lower right",
        fontsize=8.5,
        frameon=False,
    )

    plt.tight_layout()

    save_figure(
        fig,
        output_dir / "BioSample_Local_Per_Check_Exposure.png",
    )


# ---------------------------------------------------------------------------
# Tables and exports
# ---------------------------------------------------------------------------

def write_csv(
    rows: Iterable[Iterable[object]],
    path: Path,
) -> None:
    with path.open(
        "w",
        newline="",
        encoding="utf-8-sig",
    ) as handle:
        csv.writer(handle).writerows(rows)

    print(f"    {path.name}")


def build_tables(
    df: pd.DataFrame,
    summary: Dict[str, Dict[str, object]],
    n_records: int,
    output_dir: Path,
) -> None:
    """Write concise BioSample-local baseline tables."""
    category_counts = df["exposure_category"].value_counts()

    rows = [
        ["Metric", "Value", "Percentage (%)"],
        ["Total records audited", n_records, "100.0"],
        [
            "Records with >=1 exposure indicator",
            summary["OVERALL"]["n"],
            f"{summary['OVERALL']['pct']:.1f}",
        ],
        [
            "None",
            int(category_counts.get("None", 0)),
            f"{pct(int(category_counts.get('None', 0)), n_records):.1f}",
        ],
        [
            "Low",
            int(category_counts.get("Low", 0)),
            f"{pct(int(category_counts.get('Low', 0)), n_records):.1f}",
        ],
        [
            "Medium",
            int(category_counts.get("Medium", 0)),
            f"{pct(int(category_counts.get('Medium', 0)), n_records):.1f}",
        ],
        [
            "High",
            int(category_counts.get("High", 0)),
            f"{pct(int(category_counts.get('High', 0)), n_records):.1f}",
        ],
        [
            "Critical",
            int(category_counts.get("Critical", 0)),
            f"{pct(int(category_counts.get('Critical', 0)), n_records):.1f}",
        ],
        [
            "Mean BioSample-local exposure score",
            f"{df['exposure_score'].mean():.5f}",
            "",
        ],
        [
            "Median BioSample-local exposure score",
            f"{df['exposure_score'].median():.1f}",
            "",
        ],
        [
            "Maximum BioSample-local exposure score",
            int(df["exposure_score"].max()),
            "",
        ],
        [
            "Audit checks",
            len(ALL_CHECK_IDS),
            "",
        ],
    ]

    write_csv(
        rows,
        output_dir / "Table1_BioSample_Local_Summary.csv",
    )

    rows = [[
        "Check ID",
        "Layer",
        "Exposure condition",
        "N with indicator",
        "% with indicator",
        "N without indicator",
        "% without indicator",
    ]]

    for check_id in ALL_CHECK_IDS:
        item = summary[check_id]

        rows.append([
            check_id,
            item["layer"],
            item["description"],
            item["n"],
            f"{item['pct']:.3f}",
            n_records - item["n"],
            f"{(n_records - item['n']) / n_records * 100:.3f}",
        ])

    write_csv(
        rows,
        output_dir / "Table2_BioSample_Local_Checks.csv",
    )


def export_audit_csv(
    df: pd.DataFrame,
    output_dir: Path,
) -> None:
    """Export per-record BioSample-local indicator results."""
    keep = [
        "biosample_accession",
        "submission_date",
        "collection_date",
        "geographic_location",
        "country",
        "host",
        "isolation_source",
        "pango_lineage",
        "submitter_handle",
        "sra_accession",
        "exposure_score",
        "exposure_category",
    ] + ALL_CHECK_IDS

    keep = [
        column
        for column in keep
        if column in df.columns
    ]

    path = output_dir / "biosample_local_audit_results.csv"

    df[keep].to_csv(
        path,
        index=False,
        encoding="utf-8-sig",
    )

    print(f"    {path.name}")


def write_summary_report(
    df: pd.DataFrame,
    summary: Dict[str, Dict[str, object]],
    n_records: int,
    output_dir: Path,
) -> None:
    """
    Write a factual baseline summary.

    This report deliberately avoids presenting the BioSample-local values as
    the final study results. Program 04 produces the final cross-repository
    analysis.
    """
    category_counts = df["exposure_category"].value_counts()

    high_critical = (
        int(category_counts.get("High", 0))
        + int(category_counts.get("Critical", 0))
    )

    lines = [
        "APBT BIOSAMPLE-LOCAL BASELINE AUDIT",
        "=" * 72,
        "",
        "IMPORTANT:",
        "These results represent the BioSample-local baseline only.",
        "They are retained for methodological comparison and do not represent",
        "the final cross-repository results reported in the revised manuscript.",
        "",
        f"Records analysed: {n_records:,}",
        f"Audit checks: {len(ALL_CHECK_IDS)}",
        f"Mean exposure score: {df['exposure_score'].mean():.5f}",
        f"SD: {df['exposure_score'].std():.5f}",
        f"Median: {df['exposure_score'].median():.1f}",
        f"Maximum: {int(df['exposure_score'].max())}",
        "",
        f"Low: {int(category_counts.get('Low', 0)):,}",
        f"Medium: {int(category_counts.get('Medium', 0)):,}",
        f"High: {int(category_counts.get('High', 0)):,}",
        f"Critical: {int(category_counts.get('Critical', 0)):,}",
        (
            f"High/Critical: {high_critical:,} "
            f"({high_critical / n_records * 100:.3f}%)"
        ),
        "",
        "Final cross-repository analysis: see Program 04.",
    ]

    path = output_dir / "biosample_local_summary.txt"

    path.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )

    print(f"    {path.name}")


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the BioSample-local 23-check APBT metadata audit. "
            "This is the methodological baseline, not the final "
            "cross-repository analysis."
        )
    )

    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input CSV produced by Program 01.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results") / "biosample_local",
        help=(
            "Output directory "
            "(default: results/biosample_local)."
        ),
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_arguments()

    input_path = args.input.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()

    if not input_path.exists():
        print(f"ERROR: Input CSV not found: {input_path}")
        sys.exit(1)

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    figure_dir = output_dir / "figures"
    table_dir = output_dir / "tables"

    figure_dir.mkdir(exist_ok=True)
    table_dir.mkdir(exist_ok=True)

    print()
    print("=" * 72)
    print("  APBT BIOSAMPLE-LOCAL METADATA AUDIT")
    print("=" * 72)
    print(f"  Input       : {input_path}")
    print(f"  Output      : {output_dir}")
    print(f"  Study cutoff: {AUDIT_CUTOFF.date()}")
    print()
    print(
        "  NOTE: This program reproduces the BioSample-local "
        "baseline only."
    )
    print(
        "        Final cross-repository results are produced "
        "by Program 04."
    )
    print("=" * 72)
    print()

    start_time = time.time()

    try:
        df = pd.read_csv(
            input_path,
            encoding="utf-8-sig",
            low_memory=False,
        )

        ensure_required_columns(df)

    except Exception as exc:
        print(f"ERROR loading input CSV: {exc}")
        sys.exit(1)

    print(f"  {len(df):,} records loaded.")

    df = run_audit(df)
    summary, n_records = build_summary(df)

    print("\n  Generating baseline figures...")
    figure_overview(
        df,
        summary,
        n_records,
        figure_dir,
    )
    figure_per_check(
        summary,
        figure_dir,
    )

    print("\n  Writing baseline tables...")
    build_tables(
        df,
        summary,
        n_records,
        table_dir,
    )

    print("\n  Exporting per-record baseline results...")
    export_audit_csv(
        df,
        output_dir,
    )

    print("\n  Writing baseline summary...")
    write_summary_report(
        df,
        summary,
        n_records,
        output_dir,
    )

    elapsed = time.time() - start_time

    print()
    print("=" * 72)
    print("  BIOSAMPLE-LOCAL BASELINE COMPLETE")
    print("=" * 72)
    print(f"  Records analysed : {n_records:,}")
    print(f"  Total time       : {elapsed:.1f} seconds")
    print(f"  Output folder    : {output_dir}")
    print()
    print("  Next step: Program 03 performs cross-repository linkage.")
    print("=" * 72)


if __name__ == "__main__":
    main()
