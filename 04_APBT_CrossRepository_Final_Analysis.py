#!/usr/bin/env python3
"""
04_APBT_CrossRepository_Final_Analysis.py
=========================================

APBT Study 1: Final cross-repository metadata exposure analysis

This program applies the 23-check APBT metadata audit to the enriched
100,000-record SARS-CoV-2 BioSample cohort produced by Program 03.

For P2, P3, P4, A1, A3, and A4, an exposure indicator is assigned only when
the relevant information remains ``not_found_in_queried_resources`` after
the linked NCBI resources examined by Program 03 have been queried.
Information present directly in BioSample or recovered from linked records
is treated as available.

The phrase ``not_found_in_queried_resources`` refers only to the NCBI
resources examined in this study. It does not establish absence from all
NCBI or INSDC resources and does not demonstrate deliberate manipulation.

The script produces the final composite exposure scores, layer-level and
per-check summaries, publication figures, cross-repository status tables,
and the exploratory sensitivity analyses used in the manuscript.

Study cutoff
------------
The frozen retrieval window ended on 15 April 2026. This fixed cutoff is used
for temporal checks so that the analysis is reproducible when rerun later.

Usage
-----
python 04_APBT_CrossRepository_Final_Analysis.py \
    --input data/SARS_CoV2_BioSample_metadata_100000_cross_repository_audit.csv \
    --output-dir results/cross_repository

Requirements
------------
pip install pandas numpy matplotlib
"""

import argparse
import csv
import os
import re
import sys
import textwrap
import time
import warnings
from datetime import datetime
from pathlib import Path

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Dependency check
# ---------------------------------------------------------------------------
MISSING_PACKAGES = []
for package in ["pandas", "numpy", "matplotlib"]:
    try:
        __import__(package)
    except ImportError:
        MISSING_PACKAGES.append(package)

if MISSING_PACKAGES:
    print(f"\nERROR: Missing packages: {', '.join(MISSING_PACKAGES)}")
    print(f"Run: pip install {' '.join(MISSING_PACKAGES)}\n")
    sys.exit(1)

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator

# ===========================================================================
# CONFIGURATION
# ===========================================================================

DEFAULT_OUTPUT_DIR = Path("results") / "cross_repository"

COVID_START = pd.Timestamp("2019-11-01")
AUDIT_CUTOFF = pd.Timestamp("2026-04-15")

PLACEHOLDER_VALS = {
    "not applicable", "not provided", "na", "n/a", "none", "unknown",
    "missing", "not available", "undetermined", "not collected",
    "not determined", "restricted access", "missing: not provided",
    "missing: not collected", "missing: third party data", "not given",
}

NON_HUMAN_HOST = {
    "vero", "vero-e6", "vero-e6-tmprss2", "bat", "mink", "cat", "dog",
    "not applicable", "animal", "mus musculus", "hamster",
}

IMPLAUSIBLE_ISO = {
    "library", "plasmid", "pathogen bank", "cultured virus",
    "cell culture", "vero", "vero-e6", "laboratory",
}

HUMAN_HOST_TERMS_SENSITIVITY = {
    "homo sapiens",
    "homo sapiens [ncbitaxon:9606]",
}

CROSS_REPOSITORY_CHECKS = ["P2", "P3", "P4", "A1", "A3", "A4"]

CROSS_REPOSITORY_FLAG_COLUMNS = {
    "P2": "P2_crossrepo_not_found",
    "P3": "P3_crossrepo_not_found",
    "P4": "P4_crossrepo_not_found",
    "A1": "A1_crossrepo_not_found",
    "A3": "A3_crossrepo_not_found",
    "A4": "A4_crossrepo_not_found",
}

CROSS_REPOSITORY_STATUS_COLUMNS = {
    "P2": "P2_crossrepo_status",
    "P3": "P3_crossrepo_status",
    "P4": "P4_crossrepo_status",
    "A1": "A1_crossrepo_status",
    "A3": "A3_crossrepo_status",
    "A4": "A4_crossrepo_status",
}

# ===========================================================================
# PUBLICATION FIGURE SPECIFICATION
# ===========================================================================

MM_TO_IN = 1.0 / 25.4
W_SINGLE = 85 * MM_TO_IN
W_ONEHALF = 120 * MM_TO_IN
W_DOUBLE = 180 * MM_TO_IN
H_MAX = 220 * MM_TO_IN

FIG_DPI = 600

# ---------------------------------------------------------------------------
# Colour scheme – STRICTLY CONSISTENT WITH FIGURE 1
# This is the single source of truth for all colours in all figures.
# ---------------------------------------------------------------------------
C = {
    # Layer colours (exactly as used in Figure 1A)
    "temporal": "#0B5394",  # dark blue
    "geographic": "#3D85C6",  # mid blue
    "host": "#E69F00",  # orange
    "provenance": "#56B4E9",  # light sky blue
    "technical": "#4A8EC2",  # medium-light blue

    # Severity ramp colours (exactly as used in Figure 1B)
    "none": "#E0E0E0",  # light gray
    "low": "#C6DBEF",  # light blue
    "medium": "#6BAED6",  # medium blue
    "high": "#2171B5",  # dark blue
    "critical": "#08306B",  # very dark blue

    # Flag colour (same hex as Host, used for any flagged/missing/anomalous value)
    "flag": "#E69F00",  # orange

    # Structural colours
    "grid": "#E8E8E8",
    "text_muted": "#4D4D4D",
    "rule": "#999999",
}

SEVERITY_COLOURS = [
    C["none"], C["low"], C["medium"], C["high"], C["critical"]
]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans",
                        "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 8.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 8,
    "axes.labelweight": "normal",
    "axes.linewidth": 0.7,
    "axes.edgecolor": "#333333",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.axisbelow": True,
    "grid.color": C["grid"],
    "grid.linewidth": 0.5,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "legend.fontsize": 7,
    "lines.linewidth": 1.2,
    "figure.dpi": 150,
    "savefig.dpi": FIG_DPI,
    "savefig.facecolor": "white",
    "figure.facecolor": "white",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "figure.constrained_layout.use": True,
    "figure.constrained_layout.w_pad": 0.06,
    "figure.constrained_layout.h_pad": 0.06,
})

LAYER_ORDER = ["Temporal", "Geographic", "Host", "Provenance", "Technical"]
CATEGORY_ORDER = ["None", "Low", "Medium", "High", "Critical"]

# ===========================================================================
# CHECK DEFINITIONS
# ===========================================================================

CHECK_META = {
    "T1": ("Temporal", "Missing date",
           "Missing or unparseable collection date"),
    "T2": ("Temporal", "Year-only precision",
           "Year-only collection-date precision"),
    "T3": ("Temporal", "Month-year precision",
           "Month-year collection-date precision"),
    "T4": ("Temporal", "Pre-COVID date",
           "Collection date pre-dates COVID-19 emergence (<1 Nov 2019)"),
    "T5": ("Temporal", "Future date",
           "Collection date is in the future"),
    "T6": ("Temporal", "After submission",
           "Collection date after submission date"),
    "T7": ("Temporal", "Same-day submission",
           "Collection and submission on same calendar day"),

    "G1": ("Geographic", "Missing location",
           "Missing or functionally missing geographic location"),
    "G2": ("Geographic", "Country mismatch",
           "Country field inconsistent with geographic_location"),
    "G3": ("Geographic", "Placeholder geo",
           "Geographic field contains placeholder text"),

    "H1": ("Host", "Missing host",
           "Missing or functionally missing host organism field"),
    "H2": ("Host", "Non-human host",
           "Non-human or laboratory-associated host by exact-match screen"),
    "H3": ("Host", "Missing isolation src",
           "Missing or functionally missing isolation source"),
    "H4": ("Host", "Placeholder isolation",
           "Isolation source placeholder text not already captured under H3"),
    "H5": ("Host", "Implausible isolation",
           "Implausible isolation source for a human respiratory virus"),

    "P1": ("Provenance", "No collector",
           "Missing information on who collected the sample"),
    "P2": ("Provenance", "BioProject not recoverable",
           "BioProject association not recoverable from BioSample or linked NCBI records"),
    "P3": ("Provenance", "Sequence link not recoverable",
           "Nucleotide/GenBank linkage not recoverable from queried NCBI resources"),
    "P4": ("Provenance", "SRA link not recoverable",
           "SRA linkage not recoverable from BioSample or linked NCBI records"),

    "A1": ("Technical", "Platform not recoverable",
           "Sequencing platform or technology not recoverable from queried NCBI resources"),
    "A2": ("Technical", "No lineage",
           "Missing Pango lineage assignment in BioSample"),
    "A3": ("Technical", "Coverage not recoverable",
           "Coverage or sequencing depth not recoverable from queried NCBI resources"),
    "A4": ("Technical", "Assembly method not recoverable",
           "Assembly method not recoverable from queried NCBI resources"),
}

LAYER_CHECKS = {
    layer: [cid for cid, (check_layer, _, _) in CHECK_META.items()
            if check_layer == layer]
    for layer in LAYER_ORDER
}

ALL_CHECK_IDS = list(CHECK_META.keys())


# ===========================================================================
# HELPERS
# ===========================================================================

def is_blank(value):
    if pd.isna(value):
        return True
    text = str(value).strip().lower()
    return text == "" or text in PLACEHOLDER_VALS


def parse_date(value):
    if is_blank(value):
        return pd.NaT, "missing"
    text = str(value).strip()
    if re.match(r"^\d{4}-\d{2}-\d{2}", text):
        try:
            return pd.Timestamp(text[:10]), "full"
        except Exception:
            return pd.NaT, "unparseable"
    if re.match(r"^\d{2}/\d{2}/\d{4}$", text):
        try:
            return pd.to_datetime(text, format="%d/%m/%Y"), "full"
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


def extract_country(geo):
    if is_blank(geo):
        return ""
    text = str(geo).strip()
    for sep in [":", "/", ","]:
        if sep in text:
            return text.split(sep)[0].strip()
    return text


def pct(count, total):
    return round(count / total * 100, 1) if total else 0.0


def exposure_category(score):
    if score == 0:
        return "None"
    if score <= 2:
        return "Low"
    if score <= 5:
        return "Medium"
    if score <= 9:
        return "High"
    return "Critical"


def numeric_flag(series, column_name):
    numeric = pd.to_numeric(series, errors="coerce")
    if numeric.isna().any():
        bad = int(numeric.isna().sum())
        raise ValueError(
            f"{column_name} contains {bad:,} missing/non-numeric values. "
            "The enriched cross-repository dataset appears incomplete."
        )
    invalid = ~numeric.isin([0, 1])
    if invalid.any():
        values = sorted(numeric[invalid].unique().tolist())
        raise ValueError(
            f"{column_name} contains values other than 0/1: {values}"
        )
    return numeric.astype(bool)


def shorten(text, limit):
    text = str(text)
    if len(text) <= limit:
        return text
    return text[:limit - 1].rstrip() + "\u2026"


def panel_label(ax, letter, x=-0.10, y=1.04):
    ax.text(
        x, y, letter,
        transform=ax.transAxes,
        fontsize=9,
        fontweight="bold",
        va="bottom",
        ha="left",
        color="#111111",
    )


# ===========================================================================
# INPUT VALIDATION
# ===========================================================================

def validate_input_columns(df):
    base_required = [
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
    ]
    crossrepo_required = (
            list(CROSS_REPOSITORY_FLAG_COLUMNS.values())
            + list(CROSS_REPOSITORY_STATUS_COLUMNS.values())
    )
    missing = [c for c in base_required + crossrepo_required
               if c not in df.columns]
    if missing:
        raise ValueError(
            "\nThe selected file is not the required enriched "
            "cross-repository dataset.\nMissing columns:\n  "
            + "\n  ".join(missing)
        )
    if df["biosample_accession"].duplicated().any():
        duplicates = int(df["biosample_accession"].duplicated().sum())
        raise ValueError(
            f"Input contains {duplicates:,} duplicated BioSample accessions. "
            "Expected one row per BioSample."
        )


def validate_crossrepo_statuses(df):
    problems = []
    for code in CROSS_REPOSITORY_CHECKS:
        status_col = CROSS_REPOSITORY_STATUS_COLUMNS[code]
        flag_col = CROSS_REPOSITORY_FLAG_COLUMNS[code]
        flag = numeric_flag(df[flag_col], flag_col)
        expected = (
                df[status_col].fillna("").astype(str).str.strip()
                == "not_found_in_queried_resources"
        )
        mismatch = int((flag != expected).sum())
        if mismatch:
            problems.append(
                f"{code}: {mismatch:,} rows disagree between "
                f"{status_col} and {flag_col}"
            )
    if problems:
        raise ValueError(
            "\nCross-repository status validation failed:\n  "
            + "\n  ".join(problems)
        )


# ===========================================================================
# AUDIT ENGINE
# ===========================================================================

def run_audit(df):
    n = len(df)
    print(f"    Parsing dates for {n:,} records...")
    parsed = df["collection_date"].apply(parse_date)
    df["_col_ts"] = parsed.apply(lambda x: x[0])
    df["_col_prec"] = parsed.apply(lambda x: x[1])
    df["_sub_ts"] = pd.to_datetime(df["submission_date"], errors="coerce")
    print(f"    Running {len(ALL_CHECK_IDS)} checks...")
    # Temporal
    df["T1"] = df["_col_prec"].isin(["missing", "unparseable"])
    df["T2"] = df["_col_prec"] == "year-only"
    df["T3"] = df["_col_prec"] == "month-year"
    df["T4"] = df["_col_ts"].notna() & (df["_col_ts"] < COVID_START)
    df["T5"] = df["_col_ts"].notna() & (df["_col_ts"] > AUDIT_CUTOFF)
    df["T6"] = (
            df["_col_ts"].notna()
            & df["_sub_ts"].notna()
            & (df["_col_ts"] > df["_sub_ts"])
    )
    df["T7"] = (
            df["_col_ts"].notna()
            & df["_sub_ts"].notna()
            & (df["_col_ts"].dt.date == df["_sub_ts"].dt.date)
    )
    # Geographic
    df["G1"] = df["geographic_location"].apply(is_blank)
    df["G3"] = df["geographic_location"].apply(
        lambda v: (
                not pd.isna(v)
                and str(v).strip().lower() in PLACEHOLDER_VALS
        )
    )

    def _mismatch(row):
        geo = str(row.get("geographic_location", "")).strip()
        country = str(row.get("country", "")).strip()
        if is_blank(geo) or is_blank(country):
            return False
        expected = extract_country(geo).lower()
        return (
                expected != country.lower()
                and not geo.lower().startswith(country.lower())
        )

    df["G2"] = df.apply(_mismatch, axis=1)
    # Host
    df["H1"] = df["host"].apply(is_blank)
    df["H2"] = df["host"].apply(
        lambda v: (
                not is_blank(v)
                and str(v).strip().lower() in NON_HUMAN_HOST
        )
    )
    df["H3"] = df["isolation_source"].apply(is_blank)
    df["H4"] = df["isolation_source"].apply(
        lambda v: (
                not is_blank(v)
                and str(v).strip().lower() in PLACEHOLDER_VALS
        )
    )
    df["H5"] = df["isolation_source"].apply(
        lambda v: (
                not is_blank(v)
                and any(
            term in str(v).strip().lower()
            for term in IMPLAUSIBLE_ISO
        )
        )
    )
    # Provenance
    df["P1"] = df["collected_by"].apply(is_blank)
    df["P2_biosample_local"] = df["bioproject"].apply(is_blank)
    df["P3_biosample_local"] = df["genbank_accession"].apply(is_blank)
    df["P4_biosample_local"] = df["sra_accession"].apply(is_blank)
    df["P2"] = numeric_flag(df["P2_crossrepo_not_found"], "P2_crossrepo_not_found")
    df["P3"] = numeric_flag(df["P3_crossrepo_not_found"], "P3_crossrepo_not_found")
    df["P4"] = numeric_flag(df["P4_crossrepo_not_found"], "P4_crossrepo_not_found")
    # Technical
    df["A1_biosample_local"] = df["sequencing_platform"].apply(is_blank)
    df["A2"] = df["pango_lineage"].apply(is_blank)
    df["A3_biosample_local"] = df["coverage_depth"].apply(is_blank)
    df["A4_biosample_local"] = df["assembly_method"].apply(is_blank)
    df["A1"] = numeric_flag(df["A1_crossrepo_not_found"], "A1_crossrepo_not_found")
    df["A3"] = numeric_flag(df["A3_crossrepo_not_found"], "A3_crossrepo_not_found")
    df["A4"] = numeric_flag(df["A4_crossrepo_not_found"], "A4_crossrepo_not_found")
    # Composite
    df["exposure_score"] = df[ALL_CHECK_IDS].astype(int).sum(axis=1)
    df["exposure_category"] = df["exposure_score"].apply(exposure_category)
    return df


# ===========================================================================
# SUMMARY
# ===========================================================================

def build_summary(df):
    n = len(df)
    summary = {}
    for cid in ALL_CHECK_IDS:
        layer, label, description = CHECK_META[cid]
        count = int(df[cid].sum())
        summary[cid] = {
            "check_id": cid,
            "layer": layer,
            "label": label,
            "description": description,
            "n": count,
            "pct": pct(count, n),
        }
    for layer in LAYER_ORDER:
        cols = LAYER_CHECKS[layer]
        count = int(df[cols].any(axis=1).sum())
        summary[f"L_{layer}"] = {
            "check_id": f"L_{layer}",
            "layer": layer,
            "label": f"Any {layer}",
            "description": f"At least one {layer} exposure indicator",
            "n": count,
            "pct": pct(count, n),
        }
    count = int((df["exposure_score"] > 0).sum())
    summary["OVERALL"] = {
        "check_id": "OVERALL",
        "layer": "All",
        "label": "Any exposure indicator",
        "description": "At least one audit indicator triggered",
        "n": count,
        "pct": pct(count, n),
    }
    return summary, n


# ===========================================================================
# VALIDATION
# ===========================================================================

EXPECTED_LOCAL_COUNTS = {
    "T1": 1116, "T2": 1279, "T3": 1409, "T4": 0, "T5": 0, "T6": 0, "T7": 0,
    "G1": 479, "G2": 0, "G3": 479,
    "H1": 25581, "H2": 12, "H3": 50908, "H4": 0, "H5": 326,
    "P1": 51348, "P2": 99280, "P3": 100000, "P4": 3744,
    "A1": 99909, "A2": 98139, "A3": 100000, "A4": 100000,
}

EXPECTED_CROSS_REPOSITORY_COUNTS = {
    "P2": 28844,
    "P3": 54557,
    "P4": 3646,
    "A1": 30556,
    "A3": 99723,
    "A4": 73972,
}


def reconstructed_local_count(df, code):
    if code in CROSS_REPOSITORY_CHECKS:
        return int(df[f"{code}_biosample_local"].sum())
    return int(df[code].sum())


def write_validation_report(df, output_path):
    lines = []
    lines.append("=" * 78)
    lines.append("APBT FINAL CROSS-REPOSITORY ANALYSIS - VALIDATION REPORT")
    lines.append("=" * 78)
    lines.append(f"Records loaded: {len(df):,}")
    lines.append("")
    if len(df) != 100000:
        lines.append(
            "NOTE: Input is not the frozen 100,000-record manuscript cohort. "
            "Exact manuscript-count validation was skipped."
        )
        output_path.write_text("\n".join(lines), encoding="utf-8")
        return
    lines.append("A. Reconstructed BioSample-local check counts")
    lines.append("-" * 78)
    all_local_match = True
    for code in ALL_CHECK_IDS:
        actual = reconstructed_local_count(df, code)
        expected = EXPECTED_LOCAL_COUNTS[code]
        status = "PASS" if actual == expected else "FAIL"
        if status == "FAIL":
            all_local_match = False
        lines.append(
            f"{code:>2}  actual={actual:>7,}  expected={expected:>7,}  {status}"
        )
    lines.append("")
    lines.append(
        "BioSample-local reconstruction: "
        + ("PASS" if all_local_match else "FAIL")
    )
    lines.append("")
    lines.append("B. Cross-repository check counts")
    lines.append("-" * 78)
    all_crossrepo_match = True
    for code in CROSS_REPOSITORY_CHECKS:
        actual = int(df[code].sum())
        expected = EXPECTED_CROSS_REPOSITORY_COUNTS[code]
        status = "PASS" if actual == expected else "FAIL"
        if status == "FAIL":
            all_crossrepo_match = False
        lines.append(
            f"{code:>2}  actual={actual:>7,}  expected={expected:>7,}  {status}"
        )
    lines.append("")
    lines.append(
        "Cross-repository checks: "
        + ("PASS" if all_crossrepo_match else "FAIL")
    )
    lines.append("")
    mean_score = float(df["exposure_score"].mean())
    sd_score = float(df["exposure_score"].std())
    median_score = float(df["exposure_score"].median())
    high_critical = int(
        df["exposure_category"].isin(["High", "Critical"]).sum()
    )
    lines.append("C. Primary composite exposure score")
    lines.append("-" * 78)
    lines.append(f"Mean score       : {mean_score:.5f}")
    lines.append(f"SD               : {sd_score:.5f}")
    lines.append(f"Median           : {median_score:.1f}")
    lines.append(f"High/Critical n  : {high_critical:,}")
    lines.append(
        f"High/Critical %  : {high_critical / len(df) * 100:.3f}%"
    )
    score_pass = (
            abs(mean_score - 5.22374) < 0.00001
            and high_critical == 51492
    )
    lines.append("")
    lines.append(
        "Primary manuscript-cohort score: "
        + ("PASS" if score_pass else "FAIL")
    )
    lines.append("")
    lines.append(
        "If any validation line above "
        "reads FAIL, the input file or audit logic differs from the frozen manuscript cohort. "
        ""
    )
    output_path.write_text("\n".join(lines), encoding="utf-8")


# ===========================================================================
# FIGURE SAVING (STRICTLY PNG ONLY)
# ===========================================================================

def _save(fig, outdir, stem):
    """Saves the figure strictly as a high-DPI PNG. No PDF or other formats."""
    path = os.path.join(outdir, f"{stem}.png")
    fig.savefig(path, dpi=FIG_DPI, facecolor="white", format="png")
    plt.close(fig)
    kb = os.path.getsize(path) // 1024
    print(f"      {os.path.basename(path)}  ({kb} KB)")


# ---------------------------------------------------------------------------
# Figure 1
# ---------------------------------------------------------------------------
def figure_1(df, summary, n, outdir):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(W_DOUBLE, 2.95))
    layers = LAYER_ORDER
    percentages = [summary[f"L_{layer}"]["pct"] for layer in layers]
    colours = [C[layer.lower()] for layer in layers]
    bars = ax1.barh(layers, percentages, color=colours, height=0.58, edgecolor="white", linewidth=0.6)
    for bar, percentage in zip(bars, percentages):
        ax1.text(percentage + 1.0, bar.get_y() + bar.get_height() / 2, f"{percentage:.1f}",
                 va="center", fontsize=7, fontweight="bold", color="#222222")
    ax1.set_xlabel("Records with at least one indicator (%)")
    ax1.set_title("Exposure rate by metadata layer")
    ax1.set_xlim(0, min(118, max(percentages) * 1.24))
    ax1.invert_yaxis()
    ax1.set_ylim(len(layers) + 0.15, -0.55)
    overall_pct = summary["OVERALL"]["pct"]
    ax1.axvline(overall_pct, color=C["rule"], linestyle="--", linewidth=0.8)
    ax1.text(overall_pct - 1.5, len(layers) - 0.28, f"Overall {overall_pct:.1f}%",
             fontsize=6.5, color=C["text_muted"], va="center", ha="right")
    panel_label(ax1, "A", x=-0.24)

    category_counts = {cat: int((df["exposure_category"] == cat).sum()) for cat in CATEGORY_ORDER}
    left = 0
    legend_handles = []
    for category, colour in zip(CATEGORY_ORDER, SEVERITY_COLOURS):
        percentage = category_counts[category] / n * 100
        ax2.barh([0], percentage, left=left, color=colour, height=0.55, edgecolor="white", linewidth=0.6)
        legend_handles.append(
            mpatches.Patch(facecolor=colour, edgecolor="#8C8C8C", linewidth=0.4,
                           label=f"{category} ({percentage:.1f}%)")
        )
        if percentage > 6:
            text_colour = "white" if category in ("High", "Critical") else "#222222"
            ax2.text(left + percentage / 2, 0, f"{percentage:.1f}%",
                     ha="center", va="center", fontsize=7, fontweight="bold", color=text_colour)
        left += percentage
    ax2.set_xlabel("Records (%)")
    ax2.set_title("Composite exposure category")
    ax2.set_xlim(0, 100)
    ax2.set_ylim(-0.88, 0.48)
    ax2.set_yticks([])
    ax2.grid(False)
    ax2.legend(handles=legend_handles, loc="lower center", bbox_to_anchor=(0.5, 0.0),
               ncol=3, fontsize=6.5, frameon=False, handlelength=1.1, handleheight=0.9,
               columnspacing=1.0, labelspacing=0.5)
    high_critical = category_counts["High"] + category_counts["Critical"]
    ax2.text(0.99, 0.38,
             f"n = {n:,}\nMean score {df['exposure_score'].mean():.2f}\nHigh + Critical {high_critical / n * 100:.1f}%",
             transform=ax2.transAxes, ha="right", va="top", fontsize=6.5, color=C["text_muted"], linespacing=1.4)
    panel_label(ax2, "B", x=-0.06)
    _save(fig, outdir, "Exposure_Overview")


# ---------------------------------------------------------------------------
# Figure 2
# ---------------------------------------------------------------------------
def figure_2(df, summary, n, outdir):
    """
    Figure 2: exposure indicator rate for each audit check.
    Uses a discrete blue ramp (dark -> light) based on percentage bins,
    consistent with Figure 1B and Figure 3.
    """
    from matplotlib.colors import LinearSegmentedColormap

    # Build rows and sort descending by percentage
    rows = [(cid, CHECK_META[cid][0], CHECK_META[cid][1], summary[cid]["pct"])
            for cid in ALL_CHECK_IDS if summary[cid]["n"] > 0]
    rows.sort(key=lambda x: -x[3])

    labels = [f"{row[0]}  {shorten(row[2], 32)}" for row in rows]
    percentages = [row[3] for row in rows]

    # Define bins and corresponding colours (high percentage -> dark)
    # Bins: 0-5%, 5-20%, 20-40%, 40-70%, >70%
    bins = [0, 5, 20, 40, 70, 100]
    # Colours in order of increasing percentage: light to dark
    # We'll map each bin to a severity colour
    # For bin 0-5%: C["none"] (lightest), 5-20%: C["low"], 20-40%: C["medium"],
    # 40-70%: C["high"], >70%: C["critical"] (darkest)
    bin_colours = [C["none"], C["low"], C["medium"], C["high"], C["critical"]]

    # Assign colour to each percentage
    colours = []
    for pct in percentages:
        for i, (low, high) in enumerate(zip(bins[:-1], bins[1:])):
            if low <= pct < high:
                colours.append(bin_colours[i])
                break
        else:
            colours.append(C["critical"])  # fallback

    height = min(H_MAX, max(3.2, len(rows) * 0.235 + 0.95))
    fig, ax = plt.subplots(figsize=(W_DOUBLE, height))

    bars = ax.barh(labels, percentages, color=colours, height=0.68,
                   edgecolor="white", linewidth=0.5)

    # Add percentage labels
    for bar, pct in zip(bars, percentages):
        ax.text(bar.get_width() + 0.8, bar.get_y() + bar.get_height() / 2,
                f"{pct:.1f}", va="center", fontsize=6.5, color="#222222")

    ax.set_xlabel("Records with exposure indicator (%)")
    ax.set_title("Exposure indicator rate for each audit check")
    ax.set_xlim(0, min(112, max(percentages) * 1.16) if percentages else 100)
    ax.invert_yaxis()
    ax.tick_params(axis="y", length=0)

    # Optional: add a small legend showing the ramp (can be omitted)
    # We'll add a subtle colorbar-like patch or just leave it as is.
    # For clarity, we can add a note or keep the figure self-explanatory.

    _save(fig, outdir, "Per_Check_Exposure_Rates")
# ---------------------------------------------------------------------------
# Figure 3 – CONSISTENT WITH FIGURE 1 (HIGH VISIBILITY FIX)
# ---------------------------------------------------------------------------
def figure_3(df, summary, n, outdir):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(W_DOUBLE, 2.95))

    # Panel A: Date precision
    precision_order = ["full", "month-year", "year-only", "missing", "unparseable"]
    precision_labels = ["Full\nYYYY-MM-DD", "Month-year\nYYYY-MM", "Year only\nYYYY", "Missing", "Unparseable"]

    # FIX: Use Temporal layer color (Dark Blue) for 'Full', and Flag color (Orange) for issues.
    # This makes the chart highly visible and matches Figure 1's "Blue = Layer, Orange = Flag" rule.
    precision_colours = [
        C["temporal"],  # Full (Dark Blue #0B5394 - highly visible)
        C["flag"],  # Month-year (Orange #E69F00 - flagged as reduced)
        C["flag"],  # Year only (Orange)
        C["flag"],  # Missing (Orange)
        C["flag"]  # Unparseable (Orange)
    ]

    counts = df["_col_prec"].value_counts()
    values = [counts.get(p, 0) for p in precision_order]
    percentages = [value / n * 100 for value in values]

    bars = ax1.bar(range(len(precision_order)), percentages, color=precision_colours, width=0.62,
                   edgecolor="white", linewidth=0.6)
    for bar, count, percentage in zip(bars, values, percentages):
        ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1.2,
                 f"{count:,}\n{percentage:.1f}%", ha="center", va="bottom", fontsize=6.5, linespacing=1.3)

    ax1.set_xticks(range(len(precision_order)))
    ax1.set_xticklabels(precision_labels, fontsize=6.5, linespacing=1.3)
    ax1.set_ylabel("Records (%)")
    ax1.set_title("Collection date precision")
    ax1.set_ylim(0, 118)

    insufficient = sum(counts.get(p, 0) for p in ["month-year", "year-only", "missing", "unparseable"])
    ax1.text(0.99, 0.68, f"Reduced precision\n{insufficient:,} ({insufficient / n * 100:.1f}%)",
             transform=ax1.transAxes, ha="right", va="top", fontsize=6.5, color=C["text_muted"], linespacing=1.4)
    panel_label(ax1, "A", x=-0.14)

    # Panel B: Annual distribution (Already correct: Dark Blue for normal, Orange for anomalies)
    collection_years = df["_col_ts"].dt.year.dropna().astype(int)
    year_counts = collection_years.value_counts().sort_index()

    year_colours = []
    for year in year_counts.index:
        if year < 2019 or year > AUDIT_CUTOFF.year:
            year_colours.append(C["flag"])  # Orange for pre-emergence or future
        else:
            year_colours.append(C["temporal"])  # Dark Blue for normal years

    ax2.bar(year_counts.index, year_counts.values, color=year_colours, width=0.72,
            edgecolor="white", linewidth=0.5)

    ax2.set_xlabel("Collection year")
    ax2.set_ylabel("Records (n)")
    ax2.set_title("Annual distribution of collection dates")

    if len(year_counts) <= 14:
        ax2.set_xticks(list(year_counts.index))
        ax2.set_xticklabels([str(year) for year in year_counts.index], fontsize=6.5,
                            rotation=90 if len(year_counts) > 9 else 0)
    else:
        ax2.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=9))

    ax2.ticklabel_format(axis="y", style="plain")

    pre = int((collection_years < 2019).sum())
    if pre > 0:
        ax2.text(0.02, 0.97, f"Pre-emergence: {pre:,} records", transform=ax2.transAxes,
                 va="top", fontsize=6.5, color=C["flag"], fontweight="bold")

    panel_label(ax2, "B", x=-0.14)
    _save(fig, outdir, "Temporal_Analysis")


# ---------------------------------------------------------------------------
# Figure 4 – CONSISTENT WITH FIGURE 1 (3 DISTINCT COLORS & LEGENDS)
# ---------------------------------------------------------------------------
def figure_4(df, summary, n, outdir):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(W_DOUBLE, 3.35), gridspec_kw={"width_ratios": [1.25, 1.0]})

    # Panel A: Geographic labels
    country_counts = df["country"].value_counts().head(12)
    country_pcts = [value / n * 100 for value in country_counts.values]

    # Use Geographic layer color for valid locations, Flag color for missing/placeholder
    bar_colours = [
        C["geographic"] if str(country).strip().lower() not in PLACEHOLDER_VALS else C["flag"]
        for country in country_counts.index
    ]

    ax1.barh(range(len(country_counts)), country_pcts, color=bar_colours, height=0.62, edgecolor="white", linewidth=0.5)
    ax1.set_yticks(range(len(country_counts)))
    ax1.set_yticklabels([shorten(country, 22) for country in country_counts.index], fontsize=6.5)
    ax1.invert_yaxis()
    ax1.set_xlabel("Records (%)")
    ax1.set_title(f"{len(country_counts)} most frequent geographic labels")
    ax1.set_xlim(0, max(country_pcts) * 1.40)
    ax1.set_ylim(len(country_counts) + 0.45, -0.6)
    ax1.tick_params(axis="y", length=0)

    for i, (count, percentage) in enumerate(zip(country_counts.values, country_pcts)):
        ax1.text(percentage + max(country_pcts) * 0.02, i, f"{count:,} ({percentage:.1f}%)",
                 va="center", fontsize=6)

    # FIX: Corrected syntax error. 'in' already returns a boolean.
    has_missing = C["flag"] in bar_colours
    if has_missing:
        ax1.legend(
            handles=[
                mpatches.Patch(color=C["geographic"], label="Recorded location"),
                mpatches.Patch(color=C["flag"], label="Missing or placeholder")
            ],
            fontsize=6.5, frameon=False, loc="lower right", ncol=2,
            handlelength=1.1, handleheight=0.9, columnspacing=1.0, borderaxespad=0.2
        )
    panel_label(ax1, "A", x=-0.28)

    # Panel B: Host classification
    host_map = df["host"].fillna("Missing").apply(
        lambda value: ("Human" if str(value).lower() in {"homo sapiens", "human"} else
                       ("Missing" if is_blank(value) else "Non-human or\nambiguous"))
    )
    host_counts = host_map.value_counts()
    host_pcts = [value / n * 100 for value in host_counts.values]

    # Assign 3 distinct colors from the Figure 1 palette to the 3 distinct bars
    host_colours = []
    for host in host_counts.index:
        if str(host).lower() in ["human", "homo sapiens"]:
            host_colours.append(C["temporal"])  # Dark Blue (#0B5394) for Human (Valid/Expected)
        elif str(host).lower() == "missing":
            host_colours.append(C["flag"])  # Orange (#E69F00) for Missing (Flagged)
        else:
            host_colours.append(C["technical"])  # Medium-Light Blue (#4A8EC2) for Non-human (Distinct 3rd category)

    ax2.bar(range(len(host_counts)), host_pcts, color=host_colours, width=0.58, edgecolor="white", linewidth=0.6)
    for i, (count, percentage) in enumerate(zip(host_counts.values, host_pcts)):
        ax2.text(i, percentage + 1.2, f"{count:,}\n{percentage:.1f}%", ha="center", va="bottom", fontsize=6.5,
                 linespacing=1.3)

    ax2.set_xticks(range(len(host_counts)))
    ax2.set_xticklabels(host_counts.index, fontsize=6.5, linespacing=1.3)
    ax2.set_ylabel("Records (%)")
    ax2.set_title("Host organism classification")
    ax2.set_ylim(0, 118)

    # Legend now explicitly lists all 3 categories with their distinct colors
    ax2.legend(
        handles=[
            mpatches.Patch(color=C["temporal"], label="Human"),
            mpatches.Patch(color=C["flag"], label="Missing"),
            mpatches.Patch(color=C["technical"], label="Non-human or ambiguous")
        ],
        fontsize=6.5, frameon=False, loc="upper right"
    )
    panel_label(ax2, "B", x=-0.18)

    _save(fig, outdir, "Geographic_and_Host_Analysis")


# ---------------------------------------------------------------------------
# Figure 5 – CONSISTENT WITH FIGURE 1 & 3 (UNIFIED BLUE)
# ---------------------------------------------------------------------------
def figure_5(df, summary, n, outdir):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(W_DOUBLE, 3.0))

    def _draw_group(ax, check_ids, title, letter, label_x, no_indicator_colour):
        labels = [CHECK_META[cid][1] for cid in check_ids]
        exposure_pcts = [summary[cid]["pct"] for cid in check_ids]
        no_indicator_pcts = [100 - value for value in exposure_pcts]
        x = np.arange(len(check_ids))
        width = 0.36

        # 'No exposure indicator' bars
        ax.bar(x - width / 2, no_indicator_pcts, width,
               color=no_indicator_colour, edgecolor="white", linewidth=0.6)
        # 'Exposure indicator' bars (always Orange/Flag)
        ax.bar(x + width / 2, exposure_pcts, width,
               color=C["flag"], edgecolor="white", linewidth=0.6)

        # Add labels
        for i, (no_pct, exp_pct) in enumerate(zip(no_indicator_pcts, exposure_pcts)):
            ax.text(x[i] - width / 2, no_pct + 1.2, f"{no_pct:.1f}",
                    ha="center", va="bottom", fontsize=6, color=C["text_muted"])
            if exp_pct > 0:
                ax.text(x[i] + width / 2, exp_pct + 1.2, f"{exp_pct:.1f}",
                        ha="center", va="bottom", fontsize=6,
                        fontweight="bold", color=C["flag"])

        ax.set_xticks(x)
        ax.set_xticklabels([f"{cid}\n{textwrap.fill(label, 14)}" for cid, label in zip(check_ids, labels)],
                           fontsize=6.5, linespacing=1.25)
        ax.set_ylabel("Records (%)")
        ax.set_ylim(0, 124)
        ax.set_title(title)

        # Legend with unified colors
        handles = [mpatches.Patch(color=no_indicator_colour, label="No exposure indicator"),
                   mpatches.Patch(color=C["flag"], label="Exposure indicator")]
        ax.legend(handles=handles, fontsize=6.5, frameon=False, loc="upper center",
                  ncol=2, handlelength=1.1, handleheight=0.9, columnspacing=1.0)
        panel_label(ax, letter, x=label_x)

    # FIX: Use C["temporal"] (Dark Blue) for 'No exposure indicator' in BOTH panels.
    # This matches the primary blue from Figure 1A (Temporal) and Figure 3,
    # ensuring the "No exposure indicator" legend item is the same color everywhere.
    _draw_group(ax1, ["P1", "P2", "P3", "P4"], "Provenance layer", "A", -0.14, C["temporal"])
    _draw_group(ax2, ["A1", "A2", "A3", "A4"], "Technical and analytical layer", "B", -0.14, C["temporal"])

    _save(fig, outdir, "Provenance_and_Technical_Exposure")


# ---------------------------------------------------------------------------
# Figure 6
# ---------------------------------------------------------------------------
def figure_6(df, summary, n, outdir):
    submitter_counts = df["submitter_handle"].value_counts()
    top_n = min(12, len(submitter_counts))
    top = submitter_counts.head(top_n)
    top_pcts = [value / n * 100 for value in top.values]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(W_DOUBLE, 3.35), gridspec_kw={"width_ratios": [1.15, 1.0]})
    cumulative = 0
    bar_colours = []
    for value in top.values:
        cumulative += value / n * 100
        if cumulative <= 50:
            bar_colours.append(C["critical"])
        elif cumulative <= 80:
            bar_colours.append(C["high"])
        else:
            bar_colours.append(C["medium"])
    ax1.barh(range(top_n), top_pcts, color=bar_colours, height=0.64, edgecolor="white", linewidth=0.5)
    ax1.set_yticks(range(top_n))
    ax1.set_yticklabels([shorten(value, 28) for value in top.index], fontsize=6)
    ax1.invert_yaxis()
    ax1.set_xlabel("Share of dataset (%)")
    ax1.set_title(f"{top_n} largest submitting institutions")
    ax1.set_xlim(0, max(top_pcts) * 1.34)
    ax1.set_ylim(top_n + 0.9, -0.6)
    ax1.tick_params(axis="y", length=0)
    for i, (count, percentage) in enumerate(zip(top.values, top_pcts)):
        ax1.text(percentage + max(top_pcts) * 0.02, i, f"{count:,} ({percentage:.1f}%)", va="center", fontsize=6)
    ax1.legend(handles=[mpatches.Patch(color=C["critical"], label="Within first 50%"),
                        mpatches.Patch(color=C["high"], label="50-80%"),
                        mpatches.Patch(color=C["medium"], label="Above 80%")],
               fontsize=6, frameon=False, loc="lower right", ncol=3,
               title="Cumulative share of dataset", title_fontsize=6,
               handlelength=1.1, handleheight=0.9, columnspacing=0.9, borderaxespad=0.2)
    panel_label(ax1, "A", x=-0.34)
    sorted_values = np.sort(submitter_counts.values)[::-1]
    cumulative_pct = np.cumsum(sorted_values) / n * 100
    ranks = np.arange(1, len(cumulative_pct) + 1)
    ax2.plot(ranks, cumulative_pct, color=C["high"], linewidth=1.4, zorder=3)
    ax2.fill_between(ranks, cumulative_pct, color=C["low"], alpha=0.5)
    for threshold in (50, 80):
        idx = int(np.searchsorted(cumulative_pct, threshold))
        ax2.axhline(threshold, color=C["rule"], linestyle="--", linewidth=0.7)
        ax2.plot([idx + 1], [cumulative_pct[idx]], marker="o", markersize=3, color=C["high"], zorder=4)
        ax2.annotate(f"{idx + 1} submitters = {threshold}%", xy=(idx + 1, threshold),
                     xytext=(idx + 1 + len(ranks) * 0.04, threshold - 9),
                     fontsize=6.5, color=C["text_muted"])
    ax2.set_xlabel("Submitters ranked by volume (n)")
    ax2.set_ylabel("Cumulative share of dataset (%)")
    ax2.set_title("Submitter concentration")
    ax2.set_ylim(0, 104)
    ax2.set_xlim(0, len(ranks))
    panel_label(ax2, "B", x=-0.18)
    _save(fig, outdir, "Submitter_Concentration")


# ===========================================================================
# TABLES
# ===========================================================================

def write_csv(rows, path):
    with open(path, "w", newline="", encoding="utf-8-sig") as handle:
        csv.writer(handle).writerows(rows)


def table_path(outdir, stem):
    return os.path.join(outdir, f"{stem}.csv")


def build_tables(df, summary, n, outdir):
    categories = df["exposure_category"].value_counts()
    rows = [["Metric", "Value", "Percentage (%)"],
            ["Total records audited", f"{n:,}", "100.0"],
            ["Records with >=1 exposure indicator", f"{summary['OVERALL']['n']:,}", f"{summary['OVERALL']['pct']:.1f}"],
            ["Records with no detected exposure indicators", f"{n - summary['OVERALL']['n']:,}",
             f"{pct(n - summary['OVERALL']['n'], n):.1f}"],
            ["", "", ""],
            ["None (score = 0)", f"{categories.get('None', 0):,}", f"{pct(categories.get('None', 0), n):.1f}"],
            ["Low (score 1-2)", f"{categories.get('Low', 0):,}", f"{pct(categories.get('Low', 0), n):.1f}"],
            ["Medium (score 3-5)", f"{categories.get('Medium', 0):,}", f"{pct(categories.get('Medium', 0), n):.1f}"],
            ["High (score 6-9)", f"{categories.get('High', 0):,}", f"{pct(categories.get('High', 0), n):.1f}"],
            ["Critical (score >=10)", f"{categories.get('Critical', 0):,}",
             f"{pct(categories.get('Critical', 0), n):.1f}"],
            ["", "", ""],
            ["Mean exposure score", f"{df['exposure_score'].mean():.5f}", ""],
            ["SD exposure score", f"{df['exposure_score'].std():.5f}", ""],
            ["Median exposure score", f"{df['exposure_score'].median():.1f}", ""],
            ["Maximum exposure score", f"{int(df['exposure_score'].max())}", ""],
            ["Total audit checks performed", f"{len(ALL_CHECK_IDS)}", ""]]
    write_csv(rows, table_path(outdir, "Table1_Summary_Statistics"))

    apbt_relevance = {"Temporal": "Temporal metadata exposure relevant to molecular-clock and timeline interpretation",
                      "Geographic": "Geographic metadata exposure relevant to outbreak reconstruction",
                      "Host": "Host/specimen metadata exposure relevant to transmission-context interpretation",
                      "Provenance": "Provenance exposure relevant to attribution, linkage, and traceability",
                      "Technical": "Technical/analytical exposure relevant to interpretation and verification"}
    rows = [["Layer", "Records with >=1 exposure indicator (n)", "Percentage (%)", "Checks in layer",
             "APBT taxonomy relevance"]]
    for layer in LAYER_ORDER:
        layer_summary = summary[f"L_{layer}"]
        rows.append([layer, f"{layer_summary['n']:,}", f"{layer_summary['pct']:.1f}", len(LAYER_CHECKS[layer]),
                     apbt_relevance[layer]])
    write_csv(rows, table_path(outdir, "Table2_Layer_Exposure_Rates"))

    rows = [["Check ID", "Layer", "Label", "Description", "N flagged", "% flagged", "N not flagged", "% not flagged"]]
    for cid in ALL_CHECK_IDS:
        result = summary[cid]
        rows.append([cid, result["layer"], result["label"], result["description"],
                     result["n"], f"{result['pct']:.3f}", n - result["n"], f"{(n - result['n']) / n * 100:.3f}"])
    write_csv(rows, table_path(outdir, "Table3_Individual_Check_Results"))

    rows = [["Country", "N records", "% total", "Missing geo (G1)", "Country mismatch (G2)", "Mean exposure score"]]
    for country, group in sorted(df.groupby("country"), key=lambda x: -len(x[1]))[:20]:
        rows.append([country, len(group), f"{pct(len(group), n):.1f}",
                     int(group["G1"].sum()), int(group["G2"].sum()), f"{group['exposure_score'].mean():.2f}"])
    write_csv(rows, table_path(outdir, "Table4_Country_Breakdown"))

    rows = [
        ["Submitting institution", "N records", "% total", "Mean exposure score", "High/Critical n", "High/Critical %"]]
    for submitter, group in sorted(df.groupby("submitter_handle"), key=lambda x: -len(x[1]))[:20]:
        high_critical = int(group["exposure_category"].isin(["High", "Critical"]).sum())
        rows.append([str(submitter)[:65], len(group), f"{pct(len(group), n):.1f}",
                     f"{group['exposure_score'].mean():.2f}", high_critical, f"{pct(high_critical, len(group)):.1f}"])
    write_csv(rows, table_path(outdir, "Table5_Submitter_Analysis"))

    rows = [["Check", "BioSample-local flagged n", "BioSample-local %", "Cross-repository flagged n",
             "Cross-repository %", "Change in flagged n"]]
    for code in CROSS_REPOSITORY_CHECKS:
        old_n = int(df[f"{code}_biosample_local"].sum())
        new_n = int(df[code].sum())
        rows.append([code, old_n, f"{old_n / n * 100:.3f}", new_n, f"{new_n / n * 100:.3f}", new_n - old_n])
    write_csv(rows, table_path(outdir, "Table6_CrossRepository_Recovery_Comparison"))

    rows = [["Check", "Status", "N records", "Percentage (%)"]]
    for code in CROSS_REPOSITORY_CHECKS:
        status_col = CROSS_REPOSITORY_STATUS_COLUMNS[code]
        counts = df[status_col].fillna("").astype(str).str.strip().value_counts()
        for status, count in counts.items():
            rows.append([code, status, int(count), f"{int(count) / n * 100:.3f}"])
    write_csv(rows, table_path(outdir, "Table7_CrossRepository_Status_Breakdown"))
    print("      Main tables written.")


# ===========================================================================
# SENSITIVITY
# ===========================================================================

def run_sensitivity_analysis(df, outdir):
    n = len(df)
    primary_score = df["exposure_score"].copy()
    primary_category = df["exposure_category"].copy()
    host_lower = df["host"].fillna("").astype(str).str.strip().str.lower()
    expanded_h2 = (~df["host"].apply(is_blank) & ~host_lower.isin(HUMAN_HOST_TERMS_SENSITIVITY))
    expanded_checks = df[ALL_CHECK_IDS].copy()
    expanded_checks["H2"] = expanded_h2
    expanded_score = expanded_checks.astype(int).sum(axis=1)
    expanded_category = expanded_score.apply(exposure_category)
    no_a2_checks = expanded_checks.copy()
    no_a2_checks["A2"] = False
    no_a2_score = no_a2_checks.astype(int).sum(axis=1)
    no_a2_category = no_a2_score.apply(exposure_category)
    analyses = [
        ("Primary cross-repository 23-check score", primary_score, primary_category, int(df["H2"].sum())),
        ("Sensitivity: expanded H2 classification", expanded_score, expanded_category, int(expanded_h2.sum())),
        ("Exploratory sensitivity: expanded H2 and A2 excluded", no_a2_score, no_a2_category, int(expanded_h2.sum()))]
    rows = [
        ["Analysis", "N records", "Mean score", "SD", "Median", "None n", "Low n", "Medium n", "High n", "Critical n",
         "High/Critical n", "High/Critical %", "H2 flagged n"]]
    for name, score, category, h2_n in analyses:
        counts = category.value_counts()
        high_critical = int(category.isin(["High", "Critical"]).sum())
        rows.append([name, n, f"{score.mean():.5f}", f"{score.std():.5f}", f"{score.median():.1f}",
                     int(counts.get("None", 0)), int(counts.get("Low", 0)), int(counts.get("Medium", 0)),
                     int(counts.get("High", 0)), int(counts.get("Critical", 0)), high_critical,
                     f"{high_critical / n * 100:.3f}", h2_n])
    write_csv(rows, table_path(outdir, "Table8_Sensitivity_Analysis"))
    return {"expanded_h2_n": int(expanded_h2.sum()),
            "expanded_high_critical_n": int(expanded_category.isin(["High", "Critical"]).sum()),
            "expanded_high_critical_pct": float(expanded_category.isin(["High", "Critical"]).mean() * 100),
            "no_a2_high_critical_n": int(no_a2_category.isin(["High", "Critical"]).sum()),
            "no_a2_high_critical_pct": float(no_a2_category.isin(["High", "Critical"]).mean() * 100)}


# ===========================================================================
# REPORTS
# ===========================================================================

def write_results_report(df, summary, n, sensitivity, outdir):
    categories = df["exposure_category"].value_counts()
    high = int(categories.get("High", 0))
    critical = int(categories.get("Critical", 0))
    high_critical = high + critical
    lines = ["=" * 78,
             "APBT FINAL CROSS-REPOSITORY RESULTS",
             f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
             f"Dataset: {n:,} SARS-CoV-2 BioSample records",
             "=" * 78, "",
             "PRIMARY COMPOSITE EXPOSURE SCORE",
             "-" * 78,
             f"Mean exposure score : {df['exposure_score'].mean():.5f}",
             f"SD                  : {df['exposure_score'].std():.5f}",
             f"Median              : {df['exposure_score'].median():.1f}",
             f"Maximum             : {int(df['exposure_score'].max())}",
             "",
             f"None                : {int(categories.get('None', 0)):,}",
             f"Low                 : {int(categories.get('Low', 0)):,}",
             f"Medium              : {int(categories.get('Medium', 0)):,}",
             f"High                : {high:,}",
             f"Critical            : {critical:,}",
             f"High + Critical     : {high_critical:,} ({high_critical / n * 100:.3f}%)",
             "",
             "CROSS-REPOSITORY CHECKS",
             "-" * 78]
    for code in CROSS_REPOSITORY_CHECKS:
        old_n = int(df[f"{code}_biosample_local"].sum())
        new_n = int(df[code].sum())
        lines.append(
            f"{code}: BioSample-local {old_n:,} ({old_n / n * 100:.3f}%) -> cross-repository {new_n:,} ({new_n / n * 100:.3f}%)")
    lines.extend(["", "LAYER-LEVEL EXPOSURE", "-" * 78])
    for layer in LAYER_ORDER:
        result = summary[f"L_{layer}"]
        lines.append(f"{layer:<12}: {result['n']:,} ({result['pct']:.1f}%)")
    lines.extend(["", "SENSITIVITY ANALYSES", "-" * 78,
                  f"Expanded H2 flagged records: {sensitivity['expanded_h2_n']:,}",
                  f"Expanded-H2 High/Critical: {sensitivity['expanded_high_critical_n']:,} ({sensitivity['expanded_high_critical_pct']:.3f}%)",
                  f"Expanded-H2 + A2-excluded High/Critical: {sensitivity['no_a2_high_critical_n']:,} ({sensitivity['no_a2_high_critical_pct']:.3f}%)",
                  "", "INTERPRETATION NOTE", "-" * 78,
                  "For P2/P3/P4/A1/A3/A4, 'not found in the queried resources' means that the relevant information was not recovered from the NCBI resources examined by this analysis. It does not establish global absence from all INSDC resources and does not establish deliberate metadata manipulation.",
                  ""])
    path = os.path.join(outdir, "final_results_report.txt")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))


def write_figure_specification(outdir):
    lines = ["APBT PUBLICATION FIGURE SPECIFICATION",
             "=" * 60, "",
             "Size", "-" * 60,
             f"Two-column width : 180 mm ({W_DOUBLE:.2f} in)",
             f"Single column    : 85 mm ({W_SINGLE:.2f} in)",
             f"Maximum height   : 220 mm ({H_MAX:.2f} in)",
             "Layout engine    : constrained layout, exact saved width",
             "", "Resolution and format", "-" * 60,
             f"Raster resolution: {FIG_DPI} dpi",
             "Format written   : PNG ONLY (strictly enforced)",
             "Colour mode      : RGB", "", "Typography", "-" * 60,
             "Font stack       : Arial, Helvetica, Liberation Sans, DejaVu Sans",
             "Smallest text    : 6 pt data labels, 6.5 pt axis text, 8.5 pt titles",
             "", "Colour scheme - ONE palette, used in every figure", "-" * 60,
             "Layer colours (for categorical identity):",
             f"  Temporal    {C['temporal']}  (dark blue)",
             f"  Geographic  {C['geographic']}  (mid blue)",
             f"  Provenance  {C['provenance']}  (light sky blue)",
             f"  Technical   {C['technical']}  (medium-light blue)",
             f"  Host        {C['host']}  (orange - the one non-blue layer)",
             "", "Severity ramp (for graded scales, used in Figure 1B and Figure 3):",
             f"  None       {C['none']}", f"  Low        {C['low']}", f"  Medium     {C['medium']}",
             f"  High       {C['high']}", f"  Critical   {C['critical']}",
             "",
             "Figure 3 uses the temporal layer colour for full-precision dates and the flag accent for reduced, missing, or unparseable dates.",
             "Figure 5 uses a consistent blue for 'No exposure indicator' and the flag accent for exposure indicators.",
             f"Flag accent (missing / placeholder / exposure indicator): {C['flag']}",
             "Same hex as the Host layer colour - orange has exactly one job across the whole figure set.",
             "",
             "These swatches are the whole palette. Blue always reads as 'no problem' or a neutral category; orange always reads as 'flagged'.",
             "", "No green, purple, or pink hues are used anywhere in this figure set.",
             "", "Figure list", "-" * 60,
             "Figure 1  Exposure overview (layer rates, composite categories)",
             "Figure 2  Exposure rate for each audit check",
             "Figure 3  Temporal analysis (date precision, annual distribution)",
             "Figure 4  Geographic labels and host classification",
             "Figure 5  Provenance and technical exposure after cross-repository recovery",
             "Figure 6  Submitter concentration", ""]
    path = os.path.join(outdir, "figure_specification.txt")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))


# ===========================================================================
# FULL AUDIT EXPORT
# ===========================================================================

def export_audit_csv(df, outdir):
    keep = (["biosample_accession", "submission_date", "collection_date",
             "geographic_location", "country", "host", "isolation_source",
             "pango_lineage", "submitter_handle", "sra_accession",
             "exposure_score", "exposure_category"] + ALL_CHECK_IDS +
            ["P2_crossrepo_status", "P3_crossrepo_status", "P4_crossrepo_status",
             "A1_crossrepo_status", "A3_crossrepo_status", "A4_crossrepo_status"])
    keep = [c for c in keep if c in df.columns]
    path = os.path.join(outdir, "final_full_audit_results.csv")
    df[keep].to_csv(path, index=False, encoding="utf-8-sig")
    size_mb = os.path.getsize(path) / 1024 / 1024
    print(f"      final_full_audit_results.csv  ({size_mb:.1f} MB)")


# ===========================================================================
# COMMAND-LINE INTERFACE
# ===========================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run the final 23-check APBT cross-repository metadata exposure "
            "analysis on the enriched dataset produced by Program 03."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Enriched cross-repository CSV produced by Program 03.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Output directory (default: results/cross_repository).",
    )
    return parser.parse_args()


def resolve_paths(args):
    csv_path = args.input.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()

    if not csv_path.exists():
        raise FileNotFoundError(f"Input CSV not found: {csv_path}")

    return csv_path, output_dir


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    args = parse_args()
    csv_path, output_dir = resolve_paths(args)
    figure_dir = output_dir / "figures"
    table_dir = output_dir / "tables"
    for folder in [output_dir, figure_dir, table_dir]:
        folder.mkdir(parents=True, exist_ok=True)
    print()
    print(f"  Loading {csv_path.name} ...")
    start_time = time.time()
    df = pd.read_csv(csv_path, encoding="utf-8-sig", low_memory=False)
    n = len(df)
    print(f"  {n:,} records loaded ({time.time() - start_time:.1f}s)")
    print()
    print("  Validating enriched cross-repository fields ...")
    validate_input_columns(df)
    validate_crossrepo_statuses(df)
    print("  Input validation passed.")
    print()
    print("  Running 23-check cross-repository exposure audit ...")
    df = run_audit(df)
    summary, n = build_summary(df)
    print(f"  Audit complete ({time.time() - start_time:.1f}s)")
    print()
    print("  Writing manuscript-cohort validation report ...")
    write_validation_report(df, output_dir / "validation_report.txt")
    print("      validation_report.txt")
    print()
    print(f"  Generating 6 publication figures at {FIG_DPI} dpi (PNG ONLY)...")
    figure_1(df, summary, n, figure_dir)
    figure_2(df, summary, n, figure_dir)
    figure_3(df, summary, n, figure_dir)
    figure_4(df, summary, n, figure_dir)
    figure_5(df, summary, n, figure_dir)
    figure_6(df, summary, n, figure_dir)
    print()
    print("  Building tables ...")
    build_tables(df, summary, n, table_dir)
    print()
    print("  Running supplementary sensitivity analyses ...")
    sensitivity = run_sensitivity_analysis(df, table_dir)
    print("      Table8_Sensitivity_Analysis.csv")
    print()
    print("  Exporting per-record audit CSV ...")
    export_audit_csv(df, output_dir)
    print()
    print("  Writing reports ...")
    write_results_report(df, summary, n, sensitivity, output_dir)
    print("      final_results_report.txt")
    write_figure_specification(output_dir)
    print("      figure_specification.txt")
    print()
    categories = df["exposure_category"].value_counts()
    high = int(categories.get("High", 0))
    critical = int(categories.get("Critical", 0))
    high_critical = high + critical
    elapsed = time.time() - start_time
    print("=" * 78)
    print("  COMPLETE")
    print("=" * 78)
    print(f"  Records analysed    : {n:,}")
    print(f"  Mean exposure score : {df['exposure_score'].mean():.5f}")
    print(f"  SD                  : {df['exposure_score'].std():.5f}")
    print(f"  Median              : {df['exposure_score'].median():.1f}")
    print(f"  High                : {high:,}")
    print(f"  Critical            : {critical:,}")
    print(f"  High + Critical     : {high_critical:,} ({high_critical / n * 100:.3f}%)")
    print(f"  Total time          : {elapsed:.0f}s")
    print(f"  Output folder       : {output_dir.resolve()}")
    print("=" * 78)
    print()
    print("  Expected primary manuscript-cohort result:\n  Mean = 5.22374; High/Critical = 51,492 (51.492%).")
    print()
    print("  Open validation_report.txt and confirm all manuscript-cohort checks show PASS.")
    print()


if __name__ == "__main__":
    main()
