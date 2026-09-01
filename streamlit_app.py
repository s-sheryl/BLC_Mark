from __future__ import annotations

import base64
import html
from pathlib import Path
from typing import Optional
from io import BytesIO

import pandas as pd
import streamlit as st
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer


# ============================================================
# BLC MARK — STREAMLIT APPLICATION
# ============================================================
#
# IMPORTANT SCIENTIFIC BOUNDARY
# ------------------------------------------------------------
# This application is a READ-ONLY presentation layer.
#
# It consumes frozen Phase 5 and Phase 6 outputs.
# It does NOT:
#   - recalculate differential expression
#   - recalculate evidence scores
#   - recalculate final scores
#   - assign or modify ranks
#   - retrieve new external evidence during normal use
#   - replace missing scores with zero
#
# The scientific pipeline remains authoritative.
# ============================================================


# ============================================================
# PAGE CONFIGURATION
# ============================================================

st.set_page_config(
    page_title="BLC Mark",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)


# ============================================================
# PATHS
# ============================================================

ROOT = Path(__file__).resolve().parent

ASSETS_DIR = ROOT / "assets"
PHASE5_DIR = ROOT / "results" / "phase5"
PHASE6_DIR = ROOT / "results" / "phase6"

PHASE6_TABLE_DIR = PHASE6_DIR / "tables"
PHASE6_FIGURE_DIR = PHASE6_DIR / "figures"
PHASE6_REPORT_DIR = PHASE6_DIR / "report"
PHASE6_REPRO_DIR = PHASE6_DIR / "reproducibility"


# ============================================================
# COHORT CONFIGURATION
# ============================================================

COHORTS = {
    "TCGA-BRCA": {
        "short": "BRCA",
        "name": "Breast Cancer",
        "image": "breast_cancer.png",
        "top25": "brca_top_25_biomarkers.csv",
        "expected_top_gene": "CXCL2",
    },
    "TCGA-LUAD": {
        "short": "LUAD",
        "name": "Lung Adenocarcinoma",
        "image": "lung_cancer.png",
        "top25": "luad_top_25_biomarkers.csv",
        "expected_top_gene": "BUB1B",
    },
    "TCGA-COAD": {
        "short": "COAD",
        "name": "Colorectal Cancer",
        "image": "colorectal_cancer.png",
        "top25": "coad_top_25_biomarkers.csv",
        "expected_top_gene": "GOLGA7B",
    },
}


# ============================================================
# ASSET HELPERS
# ============================================================

def asset_path(filename: str) -> Path:
    return ASSETS_DIR / filename


@st.cache_data(show_spinner=False)
def image_data_uri(path_string: str) -> str:
    path = Path(path_string)

    if not path.exists():
        return ""

    suffix = path.suffix.lower()

    mime = {
        ".png": "image/png",
        ".jpg": "image/jpeg",
        ".jpeg": "image/jpeg",
        ".webp": "image/webp",
    }.get(suffix, "image/png")

    encoded = base64.b64encode(path.read_bytes()).decode("utf-8")
    return f"data:{mime};base64,{encoded}"


def image_uri(filename: str) -> str:
    return image_data_uri(str(asset_path(filename)))


def safe(value) -> str:
    if value is None:
        return ""

    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass

    return html.escape(str(value))


# ============================================================
# COLUMN HELPERS
# ============================================================

COLUMN_ALIASES = {
    "gene": [
        "gene_id",
        "gene",
        "gene_symbol",
        "symbol",
    ],
    "rank": [
        "rank",
        "final_rank",
        "biomarker_rank",
    ],
    "final_score": [
        "final_score",
        "score",
        "biomarker_score",
        "prioritization_score",
    ],
    "effect_size": [
        "effect_size",
        "log2_fold_change",
        "log2fc",
        "fold_change",
    ],
    "adjusted_p_value": [
        "adjusted_p_value",
        "adj_p_value",
        "padj",
        "fdr",
        "q_value",
    ],
    "clinical": [
        "clinical_category",
        "clinical_evidence_category",
        "clinical_evidence",
        "clinical_status",
    ],
    "cross_cancer": [
    "cross_cancer_cohort_count",
    "cross_cancer_count",
    "cross_cancer_cohorts",
    "cohort_count",
    "cancer_count",
    "cross_cancer",
],
    "pathway_count": [
        "pathway_count",
        "pathways",
        "pathway_evidence_count",
    ],
    "de_score": [
        "de_score",
        "differential_expression_score",
    ],
    "cancer_score": [
        "cancer_association_score",
        "cancer_score",
    ],
    "clinical_score": [
        "clinical_score",
        "clinical_evidence_score",
    ],
    "cross_score": [
        "cross_cancer_score",
        "cross_cohort_score",
    ],
}


def resolve_column(df: pd.DataFrame, logical_name: str) -> Optional[str]:
    candidates = COLUMN_ALIASES.get(logical_name, [])

    for candidate in candidates:
        if candidate in df.columns:
            return candidate

    lower_lookup = {str(col).lower(): col for col in df.columns}

    for candidate in candidates:
        if candidate.lower() in lower_lookup:
            return lower_lookup[candidate.lower()]

    return None


def require_gene_column(df: pd.DataFrame) -> str:
    gene_col = resolve_column(df, "gene")

    if gene_col is None:
        raise ValueError(
            "No recognized gene identifier column was found in this file."
        )

    return gene_col


# ============================================================
# DATA LOADING — PHASE 6 TOP 25
# ============================================================

@st.cache_data(show_spinner=False)
def load_top25(cohort: str) -> pd.DataFrame:
    filename = COHORTS[cohort]["top25"]
    path = PHASE6_TABLE_DIR / filename

    if not path.exists():
        raise FileNotFoundError(
            f"Required Phase 6 table was not found: {path}"
        )

    df = pd.read_csv(path)

    gene_col = require_gene_column(df)

    # Normalized helper field only for application lookup.
    # Scientific source values are not changed.
    df = df.copy()
    df["_gene_lookup"] = (
        df[gene_col]
        .astype(str)
        .str.strip()
        .str.upper()
    )

    return df


# ============================================================
# DATA LOADING — PHASE 5 FULL CANDIDATE SET
# ============================================================

def _candidate_phase5_csvs(cohort: str) -> list[Path]:
    cohort_dir = PHASE5_DIR / cohort

    if not cohort_dir.exists():
        return []

    return sorted(cohort_dir.rglob("*.csv"))


@st.cache_data(show_spinner=False)
def discover_phase5_candidate_file(cohort: str) -> str:
    """
    Locate the frozen Phase 5 ranked candidate table without assuming an
    unverified filename.

    A valid candidate file must contain:
        - a recognizable gene column
        - final score
        - rank

    Files are NEVER modified.
    """

    candidates = _candidate_phase5_csvs(cohort)

    valid_files: list[tuple[int, Path]] = []

    for path in candidates:
        try:
            sample = pd.read_csv(path, nrows=5)
        except Exception:
            continue

        gene_col = resolve_column(sample, "gene")
        score_col = resolve_column(sample, "final_score")
        rank_col = resolve_column(sample, "rank")

        if gene_col and score_col and rank_col:
            # Prefer files whose names explicitly indicate final/ranked/
            # prioritized candidates, but do not scientifically infer data.
            name = path.name.lower()

            priority = 0

            if "rank" in name:
                priority += 4
            if "prior" in name:
                priority += 4
            if "candidate" in name:
                priority += 3
            if "biomarker" in name:
                priority += 2
            if "result" in name:
                priority += 1

            valid_files.append((priority, path))

    if not valid_files:
        raise FileNotFoundError(
            f"No Phase 5 ranked candidate CSV containing gene, rank, and "
            f"final score columns was found under: {PHASE5_DIR / cohort}"
        )

    valid_files.sort(
        key=lambda item: (
            -item[0],
            str(item[1]).lower(),
        )
    )

    return str(valid_files[0][1])


@st.cache_data(show_spinner=False)
def load_phase5_candidates(cohort: str) -> pd.DataFrame:
    path = Path(discover_phase5_candidate_file(cohort))

    df = pd.read_csv(path)

    gene_col = require_gene_column(df)

    df = df.copy()
    df["_gene_lookup"] = (
        df[gene_col]
        .astype(str)
        .str.strip()
        .str.upper()
    )

    return df


# ============================================================
# FORMATTING
# ============================================================

def format_decimal(value, digits: int = 4) -> str:
    if value is None:
        return "Unavailable"

    try:
        if pd.isna(value):
            return "Unavailable"
    except (TypeError, ValueError):
        pass

    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return str(value)


def format_score(value) -> str:
    return format_decimal(value, digits=4)


def format_effect(value) -> str:
    return format_decimal(value, digits=4)


def format_pvalue(value) -> str:
    if value is None:
        return "Unavailable"

    try:
        if pd.isna(value):
            return "Unavailable"
    except (TypeError, ValueError):
        pass

    try:
        numeric = float(value)

        if numeric == 0:
            return "0.000e+00"

        return f"{numeric:.3e}"

    except (TypeError, ValueError):
        return str(value)


def format_rank(value) -> str:
    if value is None:
        return "Unavailable"

    try:
        if pd.isna(value):
            return "Unavailable"
    except (TypeError, ValueError):
        pass

    try:
        numeric = float(value)

        if numeric.is_integer():
            return str(int(numeric))

        return str(numeric)

    except (TypeError, ValueError):
        return str(value)


def format_count(value) -> str:
    if value is None:
        return "Unavailable"

    try:
        if pd.isna(value):
            return "Unavailable"
    except (TypeError, ValueError):
        pass

    try:
        numeric = float(value)

        if numeric.is_integer():
            return str(int(numeric))

        return str(numeric)

    except (TypeError, ValueError):
        return str(value)


# ============================================================
# APP CSS
# ============================================================

st.markdown(
    """
    <style>

    /* --------------------------------------------------------
       GLOBAL
    -------------------------------------------------------- */

    :root {
        --blc-bg: #050813;
        --blc-panel: rgba(12, 18, 36, 0.92);
        --blc-panel-2: rgba(17, 25, 49, 0.92);
        --blc-border: rgba(116, 145, 255, 0.18);
        --blc-text: #f3f7ff;
        --blc-muted: #96a3bd;
        --blc-blue: #46b9ff;
        --blc-purple: #9c6cff;
        --blc-cyan: #51f0ef;
        --blc-gold: #ffc95b;
    }

    html,
    body,
    [class*="css"] {
        font-family:
            Inter,
            ui-sans-serif,
            -apple-system,
            BlinkMacSystemFont,
            "Segoe UI",
            sans-serif;
    }

    .stApp {
        background:
            radial-gradient(
                circle at 78% 8%,
                rgba(56, 93, 255, 0.12),
                transparent 26%
            ),
            radial-gradient(
                circle at 15% 38%,
                rgba(132, 71, 255, 0.08),
                transparent 27%
            ),
            #050813;
        color: #f3f7ff;
    }

    .block-container {
        max-width: 1540px;
        padding-top: 1.15rem;
        padding-bottom: 3rem;
    }

    h1, h2, h3, h4 {
        color: #f4f7ff !important;
    }

    p {
        color: #dce4f5;
    }

    /* --------------------------------------------------------
       SIDEBAR
    -------------------------------------------------------- */

    section[data-testid="stSidebar"] {
        background:
            linear-gradient(
                180deg,
                rgba(8, 13, 29, 0.99),
                rgba(5, 8, 19, 0.99)
            );
        border-right: 1px solid rgba(115, 142, 255, 0.13);
    }

    section[data-testid="stSidebar"] > div {
        padding-top: 0.8rem;
    }

    .sidebar-brand {
        text-align: center;
        padding: 0.25rem 0.2rem 0.8rem;
    }

    .sidebar-brand img {
        display: block;
        width: 170px;
        height: 115px;
        margin: 0 auto 0.2rem;
        object-fit: contain;
        object-position: center;
        border-radius: 20px;
        filter: drop-shadow(0 0 22px rgba(69, 160, 255, 0.25));
    }

    .sidebar-brand .subtext {
        margin-top: 0.15rem;
        font-size: 0.76rem;
        letter-spacing: 0.12em;
        text-transform: uppercase;
        color: #8797b9;
    }

    .sidebar-version {
        margin-top: 1rem;
        padding: 0.9rem 0.95rem;
        border: 1px solid rgba(101, 131, 235, 0.2);
        border-radius: 16px;
        background:
            linear-gradient(
                135deg,
                rgba(56, 91, 181, 0.15),
                rgba(103, 64, 176, 0.08)
            );
    }

    .sidebar-version strong {
        color: #eef4ff;
        font-size: 0.87rem;
    }

    .sidebar-version span {
        display: block;
        color: #8f9db9;
        font-size: 0.75rem;
        line-height: 1.55;
        margin-top: 0.25rem;
    }

    .sidebar-network {
        margin-top: 1rem;
        border-radius: 16px;
        overflow: hidden;
        opacity: 0.88;
    }

    .sidebar-network img {
        width: 100%;
        max-height: 170px;
        display: block;
        object-fit: cover;
    }

    .sidebar-footer {
        margin-top: 0.75rem;
        color: #72809c;
        text-align: center;
        font-size: 0.7rem;
        line-height: 1.5;
    }

    /* Radio navigation */
    section[data-testid="stSidebar"] div[role="radiogroup"] label {
        background: transparent;
        border-radius: 10px;
        padding: 0.43rem 0.55rem;
        margin-bottom: 0.16rem;
        color: #cfd7e8;
    }

    section[data-testid="stSidebar"] div[role="radiogroup"] label:hover {
        background: rgba(91, 118, 224, 0.11);
    }

    /* --------------------------------------------------------
       HERO
    -------------------------------------------------------- */

    .hero-card {
        position: relative;
        display: grid;
        grid-template-columns: minmax(0, 1.04fr) minmax(380px, 0.96fr);
        min-height: 390px;
        overflow: hidden;
        border-radius: 28px;
        border: 1px solid rgba(104, 129, 239, 0.22);
        background:
            radial-gradient(
                circle at 78% 35%,
                rgba(50, 101, 255, 0.2),
                transparent 36%
            ),
            linear-gradient(
                135deg,
                rgba(12, 18, 39, 0.98),
                rgba(7, 10, 25, 0.98)
            );
        box-shadow:
            0 24px 65px rgba(0, 0, 0, 0.25),
            inset 0 0 50px rgba(54, 78, 167, 0.045);
        margin-bottom: 1.2rem;
    }

    .hero-copy {
        position: relative;
        z-index: 2;
        padding: 3rem 1rem 2.6rem 3rem;
        display: flex;
        flex-direction: column;
        justify-content: center;
    }

    .hero-eyebrow {
        color: #64e8ff;
        font-size: 0.76rem;
        font-weight: 700;
        letter-spacing: 0.16em;
        text-transform: uppercase;
        margin-bottom: 0.65rem;
    }

    .hero-title {
        font-size: clamp(3.3rem, 5.2vw, 5.3rem);
        line-height: 0.98;
        font-weight: 800;
        letter-spacing: -0.055em;
        margin: 0;
        background:
            linear-gradient(
                90deg,
                #f4f7ff 4%,
                #9a78ff 43%,
                #4ac8ff 86%
            );
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
    }

    .hero-subtitle {
        color: #aab5cc;
        font-size: 1rem;
        margin-top: 0.85rem;
        margin-bottom: 0.85rem;
        font-weight: 600;
    }

    .hero-description {
        max-width: 670px;
        color: #c4cee1;
        line-height: 1.7;
        font-size: 0.95rem;
        margin-bottom: 1.2rem;
    }

    .hero-warning {
        max-width: 680px;
        padding: 0.8rem 1rem;
        border-radius: 13px;
        border: 1px solid rgba(66, 226, 239, 0.22);
        background: rgba(34, 180, 197, 0.07);
        color: #bdeef3;
        font-size: 0.78rem;
        line-height: 1.5;
    }

    .hero-art {
        position: relative;
        min-height: 390px;
    }

    .hero-art img {
        position: absolute;
        width: 100%;
        height: 100%;
        inset: 0;
        object-fit: cover;
        object-position: center;
        filter: saturate(1.08);
    }

    .hero-art::after {
        content: "";
        position: absolute;
        inset: 0;
        background:
            linear-gradient(
                90deg,
                rgba(7, 10, 25, 0.92) 0%,
                rgba(7, 10, 25, 0.3) 32%,
                rgba(7, 10, 25, 0.02) 68%
            );
    }

    /* --------------------------------------------------------
       SECTION HEADERS
    -------------------------------------------------------- */

    .section-heading {
        margin-top: 1.65rem;
        margin-bottom: 0.85rem;
    }

    .section-heading h2 {
        margin: 0;
        font-size: 1.48rem;
        letter-spacing: -0.025em;
    }

    .section-heading p {
        margin: 0.3rem 0 0;
        color: #8593af;
        font-size: 0.84rem;
    }

    /* --------------------------------------------------------
       COHORT CARDS
    -------------------------------------------------------- */

    .cohort-card {
        position: relative;
        min-height: 195px;
        border: 1px solid rgba(112, 139, 238, 0.18);
        border-radius: 20px;
        overflow: hidden;
        background:
            linear-gradient(
                145deg,
                rgba(16, 23, 47, 0.96),
                rgba(8, 12, 26, 0.98)
            );
        box-shadow: 0 12px 28px rgba(0, 0, 0, 0.16);
    }

    .cohort-card.active {
        border-color: rgba(78, 193, 255, 0.6);
        box-shadow:
            0 0 0 1px rgba(78, 193, 255, 0.12),
            0 15px 34px rgba(17, 82, 136, 0.18);
    }

    .cohort-image {
        position: absolute;
        right: -3%;
        bottom: -8%;
        width: 51%;
        height: 116%;
        object-fit: contain;
        object-position: center;
        filter: drop-shadow(0 0 18px rgba(95, 125, 255, 0.16));
    }

    .cohort-content {
        position: relative;
        z-index: 2;
        width: 61%;
        padding: 1.45rem 0.5rem 1rem 1.4rem;
    }

    .cohort-kicker {
        color: #6dcfff;
        font-size: 0.69rem;
        text-transform: uppercase;
        letter-spacing: 0.13em;
        font-weight: 700;
        margin-bottom: 0.45rem;
    }

    .cohort-content h3 {
        margin: 0;
        font-size: 1.04rem;
        line-height: 1.3;
    }

    .cohort-content .tcga {
        margin-top: 0.45rem;
        color: #8997b2;
        font-size: 0.8rem;
    }

    .cohort-active-pill {
        display: inline-block;
        margin-top: 0.75rem;
        padding: 0.28rem 0.58rem;
        border-radius: 999px;
        background: rgba(59, 193, 243, 0.1);
        border: 1px solid rgba(77, 207, 255, 0.17);
        color: #77dfff;
        font-size: 0.67rem;
        font-weight: 700;
        letter-spacing: 0.04em;
    }

    /* --------------------------------------------------------
       SEARCH PANEL
    -------------------------------------------------------- */

    .search-intro {
        margin-top: 1.45rem;
        padding: 1.05rem 1.15rem 0.2rem;
        border-radius: 19px 19px 0 0;
        border: 1px solid rgba(108, 137, 240, 0.2);
        border-bottom: none;
        background:
            linear-gradient(
                135deg,
                rgba(17, 26, 51, 0.9),
                rgba(9, 14, 31, 0.92)
            );
    }

    .search-intro h3 {
        margin: 0;
        font-size: 1.03rem;
    }

    .search-intro p {
        color: #8795b2;
        font-size: 0.8rem;
        margin: 0.3rem 0 0.65rem;
    }

    /* --------------------------------------------------------
       METRIC CARDS — LARGE IMAGE VERSION
    -------------------------------------------------------- */

    .metric-card {
        min-height: 178px;
        border: 1px solid rgba(104, 130, 225, 0.18);
        border-radius: 20px;
        padding: 1rem 1rem 0.95rem;
        background:
            linear-gradient(
                145deg,
                rgba(16, 23, 46, 0.96),
                rgba(8, 12, 27, 0.97)
            );
        overflow: hidden;
        position: relative;
    }

    .metric-card-inner {
        display: grid;
        grid-template-columns: 92px minmax(0, 1fr);
        gap: 0.9rem;
        align-items: center;
        height: 100%;
    }

    .metric-art {
        width: 92px;
        height: 112px;
        border-radius: 16px;
        object-fit: cover;
        object-position: center;
        background: rgba(19, 28, 54, 0.8);
        filter: saturate(1.1);
    }

    .metric-label {
        color: #8f9db8;
        text-transform: uppercase;
        font-size: 0.67rem;
        font-weight: 700;
        letter-spacing: 0.1em;
        margin-bottom: 0.32rem;
    }

    .metric-value {
        color: #f5f8ff;
        font-size: 1.8rem;
        line-height: 1.05;
        font-weight: 800;
        letter-spacing: -0.035em;
        word-break: break-word;
    }

    .metric-note {
        color: #7785a1;
        font-size: 0.7rem;
        margin-top: 0.42rem;
        line-height: 1.4;
    }

    /* --------------------------------------------------------
       GENERIC PANEL
    -------------------------------------------------------- */

    .panel {
        border: 1px solid rgba(105, 132, 230, 0.18);
        border-radius: 20px;
        background:
            linear-gradient(
                145deg,
                rgba(15, 22, 43, 0.96),
                rgba(8, 12, 27, 0.98)
            );
        padding: 1.15rem;
    }

    /* --------------------------------------------------------
       SEARCH RESULT
    -------------------------------------------------------- */

    .gene-result-card {
        margin-top: 0.25rem;
        padding: 1.25rem 1.35rem;
        border-radius: 18px;
        border: 1px solid rgba(80, 196, 255, 0.23);
        background:
            radial-gradient(
                circle at 90% 12%,
                rgba(92, 77, 255, 0.13),
                transparent 27%
            ),
            linear-gradient(
                135deg,
                rgba(16, 25, 51, 0.96),
                rgba(8, 14, 32, 0.98)
            );
    }

    .gene-result-header {
        display: flex;
        align-items: baseline;
        gap: 0.7rem;
        margin-bottom: 1rem;
    }

    .gene-name {
        color: #f7f9ff;
        font-size: 1.55rem;
        font-weight: 800;
        letter-spacing: -0.025em;
    }

    .gene-cohort {
        color: #62d4ff;
        font-size: 0.73rem;
        font-weight: 700;
        letter-spacing: 0.08em;
        text-transform: uppercase;
    }

    .result-grid {
        display: grid;
        grid-template-columns: repeat(6, minmax(0, 1fr));
        gap: 0.75rem;
    }

    .result-item {
        background: rgba(20, 30, 57, 0.6);
        border: 1px solid rgba(105, 133, 231, 0.12);
        border-radius: 12px;
        padding: 0.72rem;
        min-width: 0;
    }

    .result-item span {
        display: block;
        color: #7684a1;
        font-size: 0.64rem;
        text-transform: uppercase;
        letter-spacing: 0.07em;
        margin-bottom: 0.25rem;
    }

    .result-item strong {
        color: #eef4ff;
        font-size: 0.9rem;
        word-break: break-word;
    }

    .evidence-disclaimer {
        margin-top: 0.85rem;
        color: #8594b1;
        font-size: 0.7rem;
        line-height: 1.55;
    }

    /* --------------------------------------------------------
       FEATURE CARDS — LARGE ARTWORK
    -------------------------------------------------------- */

    .feature-card {
        min-height: 470px;
        border: 1px solid rgba(104, 132, 230, 0.17);
        border-radius: 22px;
        overflow: hidden;
        background:
            linear-gradient(
                145deg,
                rgba(15, 22, 43, 0.96),
                rgba(8, 12, 27, 0.98)
            );
        box-shadow: 0 14px 32px rgba(0, 0, 0, 0.14);
    }

    .feature-image-wrap {
        width: 100%;
        height: 300px;
        display: flex;
        align-items: center;
        justify-content: center;
        padding: 1rem;
        box-sizing: border-box;
        background:
            radial-gradient(
                circle at center,
                rgba(69, 109, 255, 0.14),
                transparent 68%
            );
        overflow: hidden;
    }

    .feature-image-wrap img {
        width: auto;
        height: auto;
        max-width: 100%;
        max-height: 100%;
        object-fit: contain;
        object-position: center;
        display: block;
    }

    .feature-card:hover .feature-image-wrap img {
        transform: none;
    }

    .feature-copy {
        padding: 1rem 1.05rem 1.2rem;
    }

    .feature-copy h3 {
        font-size: 1.03rem;
        margin: 0 0 0.45rem;
    }

    .feature-copy p {
        font-size: 0.77rem;
        color: #8e9bb5;
        line-height: 1.6;
        margin: 0;
    }

    /* --------------------------------------------------------
       INFORMATIONAL PAGE
    -------------------------------------------------------- */

    .info-hero {
        padding: 1.6rem 1.7rem;
        border-radius: 23px;
        border: 1px solid rgba(106, 133, 230, 0.18);
        background:
            linear-gradient(
                135deg,
                rgba(19, 28, 55, 0.94),
                rgba(8, 12, 26, 0.97)
            );
        margin-bottom: 1rem;
    }

    .info-hero h1 {
        margin-bottom: 0.35rem;
    }

    .info-hero p {
        color: #94a2bd;
        line-height: 1.65;
        margin: 0;
    }

    .notice {
        border: 1px solid rgba(81, 217, 234, 0.18);
        background: rgba(34, 174, 195, 0.06);
        border-radius: 14px;
        padding: 0.8rem 0.95rem;
        color: #b8e5eb;
        font-size: 0.78rem;
        line-height: 1.55;
        margin-bottom: 1rem;
    }

    /* --------------------------------------------------------
       DATAFRAME
    -------------------------------------------------------- */

    div[data-testid="stDataFrame"] {
        border: 1px solid rgba(108, 137, 235, 0.17);
        border-radius: 16px;
        overflow: hidden;
    }

    div[data-testid="stDataFrame"] * {
        font-size: 0.82rem;
    }

    /* --------------------------------------------------------
       INPUTS / BUTTONS
    -------------------------------------------------------- */

    div[data-baseweb="select"] > div,
    div[data-baseweb="input"] > div {
        background-color: rgba(13, 20, 39, 0.94) !important;
        border-color: rgba(103, 131, 232, 0.22) !important;
    }

    .stButton > button {
        width: 100%;
        border-radius: 11px;
        border: 1px solid rgba(77, 185, 255, 0.25);
        background:
            linear-gradient(
                135deg,
                rgba(67, 124, 240, 0.95),
                rgba(105, 76, 220, 0.95)
            );
        color: white;
        font-weight: 700;
    }

    .stButton > button:hover {
        border-color: rgba(99, 218, 255, 0.55);
        color: white;
    }

    /* --------------------------------------------------------
       RESPONSIVE
    -------------------------------------------------------- */

    @media (max-width: 1100px) {

        .hero-card {
            grid-template-columns: 1fr;
        }

        .hero-art {
            min-height: 330px;
        }

        .hero-copy {
            padding: 2rem;
        }

        .result-grid {
            grid-template-columns: repeat(3, minmax(0, 1fr));
        }

    }

    @media (max-width: 760px) {

        .result-grid {
            grid-template-columns: repeat(2, minmax(0, 1fr));
        }

        .metric-card-inner {
            grid-template-columns: 78px minmax(0, 1fr);
        }

        .metric-art {
            width: 78px;
            height: 100px;
        }

    }


    /* --------------------------------------------------------
       SECONDARY PAGE CARDS
    -------------------------------------------------------- */

    .info-section {
        margin-top: 1rem;
        padding: 1.15rem 1.25rem;
        border: 1px solid rgba(105, 132, 230, 0.16);
        border-radius: 18px;
        background:
            linear-gradient(
                145deg,
                rgba(15, 22, 43, 0.92),
                rgba(8, 12, 27, 0.95)
            );
    }

    .info-section h3 {
        margin: 0 0 0.55rem;
        font-size: 1.08rem;
        letter-spacing: -0.015em;
    }

    .info-section p,
    .info-section li {
        color: #a7b3c9;
        line-height: 1.65;
        font-size: 0.88rem;
    }

    .info-section p {
        margin: 0;
    }

    .info-section ul {
        margin-top: 0.55rem;
        margin-bottom: 0.55rem;
    }

    .cohort-summary-card {
        min-height: 118px;
        padding: 1rem 1.1rem;
        border: 1px solid rgba(104, 132, 230, 0.18);
        border-radius: 18px;
        background:
            linear-gradient(
                145deg,
                rgba(16, 23, 46, 0.95),
                rgba(8, 12, 27, 0.97)
            );
    }

    .cohort-summary-card .kicker {
        color: #67d8ff;
        font-size: 0.68rem;
        text-transform: uppercase;
        letter-spacing: 0.1em;
        font-weight: 700;
    }

    .cohort-summary-card .value {
        color: #f5f7ff;
        font-size: 1.38rem;
        line-height: 1.15;
        font-weight: 800;
        margin-top: 0.28rem;
    }

    .cohort-summary-card .note {
        color: #8290aa;
        font-size: 0.73rem;
        margin-top: 0.4rem;
        line-height: 1.45;
    }

    .download-card {
        padding: 1rem 1.1rem;
        border: 1px solid rgba(105, 132, 230, 0.16);
        border-radius: 16px;
        background:
            linear-gradient(
                145deg,
                rgba(15, 22, 43, 0.92),
                rgba(8, 12, 27, 0.95)
            );
    }

    .download-card .file-title {
        color: #eef4ff;
        font-size: 0.9rem;
        font-weight: 700;
    }

    .download-card .file-name {
        color: #7f8da8;
        font-size: 0.74rem;
        margin-top: 0.25rem;
        font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
        word-break: break-all;
    }

    </style>
    """,
    unsafe_allow_html=True,
)


# ============================================================
# HTML RENDERING
# ============================================================

def render_html(content: str) -> None:
    """
    Streamlit 1.62-compatible HTML rendering.

    st.html() is intentionally used instead of multiline
    unsafe_allow_html markdown, which previously rendered literal HTML.
    """
    st.html(content.strip())


# ============================================================
# SIDEBAR
# ============================================================

logo_uri = image_uri("blc_mark_logo.png")
network_uri = image_uri("sidebar_network.png")

with st.sidebar:

    if logo_uri:
        render_html(
            f"""
            <div class="sidebar-brand">
                <img src="{logo_uri}" alt="BLC Mark logo">
                <div class="subtext">
                    Cancer Biomarker Discovery
                </div>
            </div>
            """
        )
    else:
        st.title("BLC Mark")

    page = st.radio(
        "Navigation",
        [
            "Overview",
            "Biomarker Explorer",
            "Compare Cohorts",
            "Methodology",
            "About BLC Mark",
            "Download Data",
        ],
        label_visibility="collapsed",
    )

    render_html(
        """
        <div class="sidebar-version">
            <strong>BLC Mark — Version 1</strong>
            <span>
                RNA-seq gene expression<br>
                TCGA-BRCA • TCGA-LUAD • TCGA-COAD
            </span>
        </div>
        """
    )

    if network_uri:
        render_html(
            f"""
            <div class="sidebar-network">
                <img
                    src="{network_uri}"
                    alt="Decorative molecular network artwork"
                >
            </div>
            """
        )

    render_html(
        """
        <div class="sidebar-footer">
            Read-only scientific results interface<br>
            BLC Mark scientific outputs
        </div>
        """
    )


# ============================================================
# SHARED UI COMPONENTS
# ============================================================

def section_heading(title: str, subtitle: str = "") -> None:
    subtitle_html = (
        f"<p>{safe(subtitle)}</p>"
        if subtitle
        else ""
    )

    render_html(
        f"""
        <div class="section-heading">
            <h2>{safe(title)}</h2>
            {subtitle_html}
        </div>
        """
    )



def info_section(title: str, body_html: str) -> None:
    render_html(
        f"""
        <div class="info-section">
            <h3>{safe(title)}</h3>
            <div>{body_html}</div>
        </div>
        """
    )


def render_cohort_card(
    cohort: str,
    active_cohort: str,
) -> None:

    config = COHORTS[cohort]
    uri = image_uri(config["image"])

    active_class = " active" if cohort == active_cohort else ""
    active_pill = (
        '<div class="cohort-active-pill">ACTIVE COHORT</div>'
        if cohort == active_cohort
        else ""
    )

    render_html(
        f"""
        <div class="cohort-card{active_class}">
            <div class="cohort-content">
                <div class="cohort-kicker">
                    {safe(config["short"])}
                </div>

                <h3>
                    {safe(config["name"])}
                </h3>

                <div class="tcga">
                    {safe(cohort)}
                </div>

                {active_pill}
            </div>

            {
                f'<img class="cohort-image" '
                f'src="{uri}" '
                f'alt="{safe(config["name"])} artwork">'
                if uri
                else ""
            }
        </div>
        """
    )


def metric_card(
    title: str,
    value: str,
    note: str,
    image_filename: str,
) -> None:

    uri = image_uri(image_filename)

    render_html(
        f"""
        <div class="metric-card">
            <div class="metric-card-inner">

                {
                    f'<img class="metric-art" '
                    f'src="{uri}" '
                    f'alt="{safe(title)}">'
                    if uri
                    else '<div class="metric-art"></div>'
                }

                <div>
                    <div class="metric-label">
                        {safe(title)}
                    </div>

                    <div class="metric-value">
                        {safe(value)}
                    </div>

                    <div class="metric-note">
                        {safe(note)}
                    </div>
                </div>

            </div>
        </div>
        """
    )


def feature_card(
    title: str,
    description: str,
    image_filename: str,
) -> None:

    uri = image_uri(image_filename)

    render_html(
        f"""
        <div class="feature-card">

            <div class="feature-image-wrap">
                {
                    f'<img src="{uri}" '
                    f'alt="{safe(title)} decorative artwork">'
                    if uri
                    else ""
                }
            </div>

            <div class="feature-copy">
                <h3>{safe(title)}</h3>
                <p>{safe(description)}</p>
            </div>

        </div>
        """
    )


def get_value(
    row: pd.Series,
    df: pd.DataFrame,
    logical_name: str,
):
    column = resolve_column(df, logical_name)

    if column is None:
        return None

    return row.get(column)


def render_gene_result(
    row: pd.Series,
    df: pd.DataFrame,
    cohort: str,
) -> None:

    gene_col = require_gene_column(df)

    gene = row[gene_col]

    rank = format_rank(
        get_value(row, df, "rank")
    )

    score = format_score(
        get_value(row, df, "final_score")
    )

    effect = format_effect(
        get_value(row, df, "effect_size")
    )

    pvalue = format_pvalue(
        get_value(row, df, "adjusted_p_value")
    )

    cross = format_count(
        get_value(row, df, "cross_cancer")
    )

    pathways = format_count(
        get_value(row, df, "pathway_count")
    )

    clinical_value = get_value(
        row,
        df,
        "clinical",
    )

    clinical = (
        safe(clinical_value)
        if clinical_value is not None
        and not pd.isna(clinical_value)
        else "Unavailable"
    )

    render_html(
        f"""
        <div class="gene-result-card">

            <div class="gene-result-header">
                <div class="gene-name">
                    {safe(gene)}
                </div>

                <div class="gene-cohort">
                    {safe(cohort)}
                </div>
            </div>

            <div class="result-grid">

                <div class="result-item">
                    <span>Rank</span>
                    <strong>{safe(rank)}</strong>
                </div>

                <div class="result-item">
                    <span>Final score</span>
                    <strong>{safe(score)}</strong>
                </div>

                <div class="result-item">
                    <span>Effect size</span>
                    <strong>{safe(effect)}</strong>
                </div>

                <div class="result-item">
                    <span>Adjusted p-value</span>
                    <strong>{safe(pvalue)}</strong>
                </div>

                <div class="result-item">
                    <span>Cross-cancer cohorts</span>
                    <strong>{safe(cross)}</strong>
                </div>

                <div class="result-item">
                    <span>Pathway count</span>
                    <strong>{safe(pathways)}</strong>
                </div>

            </div>

            <div style="
                margin-top:0.9rem;
                padding:0.75rem;
                border-radius:11px;
                background:rgba(20,30,57,0.56);
            ">
                <span style="
                    display:block;
                    color:#7684a1;
                    font-size:0.64rem;
                    text-transform:uppercase;
                    letter-spacing:0.07em;
                    margin-bottom:0.25rem;
                ">
                    Clinical evidence category
                </span>

                <strong style="
                    color:#eef4ff;
                    font-size:0.87rem;
                ">
                    {clinical}
                </strong>
            </div>

            <div class="evidence-disclaimer">
                Clinical evidence categories reproduce source-integrated
                evidence from the frozen pipeline. They do not constitute
                independent clinical validation by BLC Mark.
            </div>

        </div>
        """
    )


# ============================================================
# CURATED TABLE
# ============================================================

def build_curated_table(df: pd.DataFrame) -> pd.DataFrame:

    gene_col = resolve_column(df, "gene")
    rank_col = resolve_column(df, "rank")
    score_col = resolve_column(df, "final_score")
    effect_col = resolve_column(df, "effect_size")
    padj_col = resolve_column(df, "adjusted_p_value")
    clinical_col = resolve_column(df, "clinical")
    cross_col = resolve_column(df, "cross_cancer")
    pathway_col = resolve_column(df, "pathway_count")

    output = pd.DataFrame(index=df.index)

    if rank_col:
        output["Rank"] = df[rank_col].map(format_rank)

    if gene_col:
        output["Gene"] = df[gene_col].astype(str)

    if score_col:
        output["Final Score"] = df[score_col].map(format_score)

    if effect_col:
        output["Effect Size"] = df[effect_col].map(format_effect)

    if padj_col:
        output["Adjusted p-value"] = df[padj_col].map(format_pvalue)

    if clinical_col:
        output["Clinical Evidence Category"] = (
            df[clinical_col]
            .fillna("Unavailable")
            .astype(str)
        )

    if cross_col:
        output["Cross-Cancer Cohorts"] = df[cross_col].map(format_count)

    if pathway_col:
        output["Pathway Count"] = df[pathway_col].map(format_count)

    return output.reset_index(drop=True)


# ============================================================
# CHART
# ============================================================

def render_top_candidate_chart(df: pd.DataFrame) -> None:

    gene_col = resolve_column(df, "gene")
    score_col = resolve_column(df, "final_score")
    rank_col = resolve_column(df, "rank")

    if gene_col is None or score_col is None:
        st.info(
            "The frozen table does not expose the columns required "
            "for the candidate score chart."
        )
        return

    chart_df = df[[gene_col, score_col]].copy()

    if rank_col:
        chart_df["_rank"] = pd.to_numeric(
            df[rank_col],
            errors="coerce",
        )

        chart_df = chart_df.sort_values(
            "_rank",
            ascending=True,
            na_position="last",
        )

    else:
        chart_df[score_col] = pd.to_numeric(
            chart_df[score_col],
            errors="coerce",
        )

        chart_df = chart_df.sort_values(
            score_col,
            ascending=False,
            na_position="last",
        )

    chart_df = chart_df.head(15)

    chart_df[score_col] = pd.to_numeric(
        chart_df[score_col],
        errors="coerce",
    )

    chart_df = chart_df.dropna(subset=[score_col])

    if chart_df.empty:
        st.info("No scored candidates are available for this chart.")
        return

    chart_df = chart_df.set_index(gene_col)

    st.bar_chart(
        chart_df[[score_col]],
        height=430,
        use_container_width=True,
    )


# ============================================================
# OVERVIEW PAGE
# ============================================================

def page_overview() -> None:

    # --------------------------------------------------------
    # HERO
    # --------------------------------------------------------

    hero_uri = image_uri("blc_mark_hero.png")

    render_html(
        f"""
        <div class="hero-card">

            <div class="hero-copy">

                <div class="hero-eyebrow">
                    Cancer Biomarker Discovery & Evidence Integration
                </div>

                <div class="hero-title">
                    BLC Mark
                </div>

                <div class="hero-subtitle">
                    Transparent prioritization of computational
                    cancer biomarker candidates
                </div>

                <div class="hero-description">
                    BLC Mark integrates RNA-seq differential-expression
                    results with structured biological and clinical
                    evidence to transparently prioritize candidate
                    biomarkers across three TCGA cancer cohorts.
                </div>

                <div class="hero-warning">
                    Computationally prioritized candidates are research
                    outputs and should not be interpreted as clinically
                    validated biomarkers or medical recommendations.
                </div>

            </div>

            <div class="hero-art">
                {
                    f'<img src="{hero_uri}" '
                    f'alt="Decorative cancer cell and DNA artwork">'
                    if hero_uri
                    else ""
                }
            </div>

        </div>
        """
    )

    # --------------------------------------------------------
    # ACTIVE COHORT
    # --------------------------------------------------------

    section_heading(
        "Cancer Cohorts",
        "Choose the active cohort used throughout the Overview dashboard.",
    )

    selected_cohort = st.selectbox(
        "Active cohort",
        list(COHORTS.keys()),
        format_func=lambda c: (
            f"{COHORTS[c]['name']} — {c}"
        ),
        key="overview_cohort",
    )

    cohort_columns = st.columns(3)

    for column, cohort in zip(
        cohort_columns,
        COHORTS.keys(),
    ):
        with column:
            render_cohort_card(
                cohort,
                selected_cohort,
            )

    # --------------------------------------------------------
    # LOAD DATA
    # --------------------------------------------------------

    try:
        top25 = load_top25(selected_cohort)

    except Exception as exc:
        st.error(
            f"Unable to load the Top-25 candidate table for "
            f"{selected_cohort}.\n\n{exc}"
        )
        return

    # --------------------------------------------------------
    # SEARCH — MOVED HIGHER
    # --------------------------------------------------------

    render_html(
        """
        <div class="search-intro">
            <h3>Search for a gene</h3>
            <p>
                Search within the selected cohort's Top-25
                candidate table.
            </p>
        </div>
        """
    )

    search_col, button_col = st.columns(
        [5.8, 1.2],
        vertical_alignment="bottom",
    )

    with search_col:
        search_gene = st.text_input(
            "Gene symbol",
            placeholder="Search Top 25 — e.g. BRCA2, CXCL2, BUB1B",
            label_visibility="collapsed",
            key="overview_gene_search",
        )

    with button_col:
        search_clicked = st.button(
            "Search",
            key="overview_search_button",
            use_container_width=True,
        )

    if search_gene.strip():

        gene_query = search_gene.strip().upper()

        matches = top25[
            top25["_gene_lookup"] == gene_query
        ]

        if not matches.empty:
            render_gene_result(
                matches.iloc[0],
                top25,
                selected_cohort,
            )

        elif search_clicked or search_gene:

            st.info(
                f"**{search_gene.strip()}** is not present in the "
                f"Top 25 candidates for {selected_cohort}. "
                f"Use **Biomarker Explorer** to search the complete "
                f"candidate set."
            )

    # --------------------------------------------------------
    # METRICS
    # --------------------------------------------------------

    section_heading(
        "Cohort Snapshot",
        "Values are read directly from the pipeline-generated results.",
    )

    gene_col = resolve_column(top25, "gene")
    score_col = resolve_column(top25, "final_score")
    rank_col = resolve_column(top25, "rank")
    cross_col = resolve_column(top25, "cross_cancer")

    if rank_col:
        ranking_numeric = pd.to_numeric(
            top25[rank_col],
            errors="coerce",
        )

        ranked = top25.assign(
            _rank_numeric=ranking_numeric
        ).sort_values(
            "_rank_numeric",
            ascending=True,
            na_position="last",
        )
    else:
        ranked = top25.copy()

    top_row = ranked.iloc[0]

    top_gene = (
        str(top_row[gene_col])
        if gene_col
        else "Unavailable"
    )

    top_score = (
        format_score(top_row[score_col])
        if score_col
        else "Unavailable"
    )

    # Overview represents the Phase 6 Top-25 table.
    candidate_count = len(top25)

    top_cross = (
        format_count(top_row[cross_col])
        if cross_col
        else "Unavailable"
    )

    m1, m2, m3, m4 = st.columns(4)

    with m1:
        metric_card(
            "Candidates",
            str(candidate_count),
            "Top-25 prioritized candidates",
            "candidates_icon.png",
        )

    with m2:
        metric_card(
            "Top Gene",
            top_gene,
            "Rank 1 in the selected cohort",
            "top_gene_icon.png",
        )

    with m3:
        metric_card(
            "Top Score",
            top_score,
            "Frozen prioritization score",
            "top_score_icon.png",
        )

    with m4:
        metric_card(
            "Cross-Cancer",
            (
                f"{top_cross}/3"
                if top_cross not in {
                    "Unavailable",
                    "",
                }
                and "/" not in top_cross
                else top_cross
            ),
            "Cohorts represented for the top candidate",
            "cross_cancer_icon.png",
        )

    # --------------------------------------------------------
    # CHART
    # --------------------------------------------------------

    section_heading(
        "Top Candidate Profile",
        "Top 15 candidates displayed from the frozen prioritization output.",
    )

    render_top_candidate_chart(top25)

    # --------------------------------------------------------
    # LARGE FULL-WIDTH TABLE
    # --------------------------------------------------------

    section_heading(
        "Top 25 Biomarker Candidates",
        (
            "Full-width scientific results table. "
            "Values are formatted only for display; "
            "underlying results are unchanged."
        ),
    )

    curated = build_curated_table(top25)

    st.dataframe(
        curated,
        use_container_width=True,
        hide_index=True,
        height=860,
    )

    st.caption(
        "Clinical Evidence Category reproduces source-integrated evidence "
        "classification and does not represent independent clinical "
        "validation by BLC Mark."
    )

    # --------------------------------------------------------
    # FEATURE CARDS
    # --------------------------------------------------------

    section_heading(
        "Platform Foundations",
        (
            "Core principles underlying BLC Mark's scientific "
            "and engineering design."
        ),
    )

    f1, f2 = st.columns(2)

    with f1:
        feature_card(
            "Evidence Integration",
            (
                "Combines differential-expression evidence with "
                "structured cancer association, clinical evidence, "
                "and pathway evidence while preserving source "
                "traceability."
            ),
            "evidence_icon.png",
        )

    with f2:
        feature_card(
            "Transparent Scoring",
            (
                "Candidate prioritization is performed by the frozen "
                "pipeline using explicit evidence dimensions and "
                "deterministic ranking rules. The web application "
                "does not recalculate scores."
            ),
            "scoring_icon.png",
        )

    f3, f4 = st.columns(2)

    with f3:
        feature_card(
            "Cross-Cancer Analysis",
            (
                "Presents evidence across TCGA-BRCA, TCGA-LUAD, "
                "and TCGA-COAD without introducing additional "
                "cohorts or recalculating the underlying evidence."
            ),
            "cross_analysis_icon.png",
        )

    with f4:
        feature_card(
            "Reproducible & Traceable",
            (
                "Frozen outputs, explicit configuration, deterministic "
                "prioritization, provenance information, and release "
                "artifacts support transparent scientific review."
            ),
            "reproducibility_icon.png",
        )


# ============================================================
# BIOMARKER EXPLORER
# ============================================================

def page_biomarker_explorer() -> None:

    render_html(
        """
        <div class="info-hero">
            <h1>Biomarker Explorer</h1>
            <p>
                Search the complete prioritized candidate set for
                an individual cancer cohort. Scores and ranks shown here
                are loaded directly from the prioritization output and
                are never recalculated by the application.
            </p>
        </div>
        """
    )

    section_heading(
        "Select a Cancer Cohort",
        "The explorer loads the complete prioritized candidate set for the selected cohort.",
    )

    cohort = st.selectbox(
        "Cancer cohort",
        list(COHORTS.keys()),
        format_func=lambda c: f"{COHORTS[c]['name']} — {c}",
        key="explorer_cohort",
    )

    try:
        candidates = load_phase5_candidates(cohort)

    except Exception as exc:
        st.error(
            f"Unable to locate/load the ranked "
            f"candidate table for {cohort}.\n\n{exc}"
        )
        return

    score_col = resolve_column(candidates, "final_score")

    if score_col:
        scored_count = int(
            pd.to_numeric(candidates[score_col], errors="coerce")
            .notna()
            .sum()
        )
    else:
        scored_count = 0

    c1, c2, c3 = st.columns(3)

    with c1:
        render_html(
            f"""
            <div class="cohort-summary-card">
                <div class="kicker">Selected cohort</div>
                <div class="value">{safe(cohort)}</div>
                <div class="note">{safe(COHORTS[cohort]["name"])}</div>
            </div>
            """
        )

    with c2:
        render_html(
            f"""
            <div class="cohort-summary-card">
                <div class="kicker">Prioritized candidates</div>
                <div class="value">{len(candidates):,}</div>
                <div class="note">Frozen candidate rows loaded</div>
            </div>
            """
        )

    with c3:
        render_html(
            f"""
            <div class="cohort-summary-card">
                <div class="kicker">Scored candidates</div>
                <div class="value">{scored_count:,}</div>
                <div class="note">Missing scores remain unavailable</div>
            </div>
            """
        )

    section_heading(
        "Search Complete Candidate Set",
        "Enter a gene symbol to inspect its prioritization record.",
    )

    search_col, button_col = st.columns(
        [5.8, 1.2],
        vertical_alignment="bottom",
    )

    with search_col:
        gene = st.text_input(
            "Search complete candidate set",
            placeholder="Enter a gene symbol — e.g. BRCA2",
            label_visibility="collapsed",
            key="explorer_gene",
        )

    with button_col:
        st.button(
            "Search",
            key="explorer_search_button",
            use_container_width=True,
        )

    if gene.strip():

        query = gene.strip().upper()

        matches = candidates[
            candidates["_gene_lookup"] == query
        ]

        if matches.empty:
            st.warning(
                f"**{gene.strip()}** was not found in the frozen "
                f"candidate set for {cohort}."
            )
        else:
            render_gene_result(
                matches.iloc[0],
                candidates,
                cohort,
            )

    st.caption(
        f"Loaded {len(candidates):,} prioritized candidates for "
        f"{cohort}."
    )

    section_heading(
        "Browse Candidate Table",
        "Open the complete frozen candidate table for wider inspection.",
    )

    with st.expander(
        f"Browse all {len(candidates):,} candidates",
        expanded=False,
    ):
        preview = build_curated_table(candidates)

        if preview.empty:
            st.dataframe(
                candidates.drop(
                    columns=["_gene_lookup"],
                    errors="ignore",
                ),
                use_container_width=True,
                hide_index=True,
                height=700,
            )
        else:
            st.dataframe(
                preview,
                use_container_width=True,
                hide_index=True,
                height=760,
            )

# ============================================================
# COMPARE COHORTS
# ============================================================

def page_compare_cohorts() -> None:

    render_html(
        """
        <div class="info-hero">
            <h1>Compare Cohorts</h1>
            <p>
                Compare the frozen Top-25 prioritization outputs across
                TCGA-BRCA, TCGA-LUAD, and TCGA-COAD. This view presents
                existing results only; it does not create new
                cross-cohort scores.
            </p>
        </div>
        """
    )

    combined_path = (
        PHASE6_TABLE_DIR
        / "combined_top_25_biomarkers.csv"
    )

    section_heading(
        "Combined Top-25 Results",
        "Integrated Top-25 results across all three cancer cohorts.",
    )

    if combined_path.exists():
        combined = pd.read_csv(combined_path)

        st.dataframe(
            combined,
            use_container_width=True,
            hide_index=True,
            height=720,
        )
    else:
        st.info(
            "The combined Top-25 table was not found. "
            "Individual cohort outputs are shown below."
        )

    section_heading(
        "Cohort Top Candidates",
        "Top-ranked candidates from each cancer cohort.",
    )

    columns = st.columns(3)

    for column, cohort in zip(
        columns,
        COHORTS.keys(),
    ):

        with column:
            config = COHORTS[cohort]

            render_html(
                f"""
                <div class="cohort-summary-card">
                    <div class="kicker">{safe(cohort)}</div>
                    <div class="value">{safe(config["name"])}</div>
                    <div class="note">Top-10 ranked candidates</div>
                </div>
                """
            )

            try:
                df = load_top25(cohort)
                compact = build_curated_table(df)

                preferred = [
                    c
                    for c in [
                        "Rank",
                        "Gene",
                        "Final Score",
                    ]
                    if c in compact.columns
                ]

                if preferred:
                    compact = compact[preferred]

                st.dataframe(
                    compact.head(10),
                    use_container_width=True,
                    hide_index=True,
                    height=405,
                )

            except Exception as exc:
                st.error(str(exc))

# ============================================================
# METHODOLOGY
# ============================================================

def page_methodology() -> None:

    render_html(
        """
        <div class="info-hero">
            <h1>Methodology</h1>
            <p>
                BLC Mark Version 1 is a reproducible RNA-seq cancer
                biomarker discovery and evidence-integration workflow
                spanning three TCGA cohorts.
            </p>
        </div>
        """
    )

    render_html(
        """
        <div class="notice">
            This page summarizes the BLC Mark Version 1 methodology.
            The Streamlit interface itself performs no differential
            expression, evidence retrieval, score calculation, or ranking.
        </div>
        """
    )

    info_section(
        "1. Data scope",
        """
        <p>
            Version 1 operates on <strong>RNA-seq gene-expression data</strong>
            from <strong>TCGA-BRCA</strong> (Breast Cancer),
            <strong>TCGA-LUAD</strong> (Lung Adenocarcinoma), and
            <strong>TCGA-COAD</strong> (Colorectal Cancer).
            The application does not introduce additional cohorts.
        </p>
        """,
    )

    info_section(
        "2. Differential expression",
        """
        <p>
            For the normalized expression representation used in the executed
            Version 1 analysis, differential expression was performed using a
            <strong>two-sided Welch two-sample t-test</strong>, allowing
            unequal variance between groups.
        </p>
        <p style="margin-top:0.7rem;">
            Multiple testing correction uses the
            <strong>Benjamini-Hochberg false discovery rate procedure</strong>
            with the configured significance threshold. Expression
            representation and statistical method are explicit pipeline
            configuration choices rather than inferred silently.
        </p>
        """,
    )

    info_section(
        "3. Evidence integration",
        """
        <p>
            Significant differential-expression candidates are connected to
            structured evidence sources. Missing or unavailable evidence is
            preserved as unavailable rather than silently interpreted as
            negative evidence.
        </p>
        """,
    )

    info_section(
        "4. Biomarker prioritization",
        """
        <p>
            Version 1 uses four evidence dimensions with equal configured
            weights: differential-expression evidence (<strong>0.25</strong>),
            cancer-association evidence (<strong>0.25</strong>),
            clinical evidence (<strong>0.25</strong>), and
            cross-cancer evidence (<strong>0.25</strong>).
        </p>
        <p style="margin-top:0.7rem;">
            Final scores and ranks are generated by the prioritization
            pipeline. Exact score ties are resolved deterministically using
            gene identifier ordering.
        </p>
        """,
    )

    info_section(
        "5. Application boundary",
        """
        <p>
            This Streamlit application is a
            <strong>read-only scientific interface</strong>. It consumes
            pipeline-generated result files and does not recompute, overwrite,
            normalize, fill, or reinterpret prioritization scores.
        </p>
        <p style="margin-top:0.7rem;">
            A missing score remains unavailable. It is never replaced with
            zero and the application never assigns a rank manually.
        </p>
        """,
    )

    st.warning(
        "BLC Mark prioritizes computational research candidates. "
        "It is not a diagnostic system and does not establish clinical "
        "validity or clinical utility."
    )

# ============================================================
# ABOUT
# ============================================================

def page_about() -> None:

    render_html(
        """
        <div class="info-hero">
            <h1>About BLC Mark</h1>
            <p>
                A research-grade cancer biomarker discovery and evidence
                integration platform designed around scientific correctness,
                reproducibility, transparency, explainability, and
                traceable outputs.
            </p>
        </div>
        """
    )

    info_section(
        "Research question",
        """
        <p>
            Can publicly available transcriptomic datasets from major
            cancer types be integrated into a transparent, reproducible
            evidence framework to systematically prioritize clinically
            relevant cancer biomarker candidates?
        </p>
        """,
    )

    info_section(
        "Version 1 scope",
        """
        <p>
            The executed Version 1 release contains three cancer cohorts:
            <strong>TCGA-BRCA, TCGA-LUAD, and TCGA-COAD</strong>, using
            RNA-seq gene-expression data.
        </p>
        <p style="margin-top:0.7rem;">
            BLC Mark is intended to support research exploration and
            reproducible candidate prioritization rather than function as
            a biomarker database or clinical decision-support system.
        </p>
        """,
    )

    info_section(
        "Scientific interpretation",
        """
        <p>
            A high BLC Mark ranking indicates prioritization within the
            Version 1 computational evidence framework. It does
            <strong>not</strong> demonstrate diagnostic performance,
            prognostic validity, therapeutic utility, or prospective
            clinical validation.
        </p>
        <p style="margin-top:0.7rem;">
            Candidate rankings should therefore be interpreted as a basis
            for further research and independent validation.
        </p>
        """,
    )

# ============================================================
# PDF EXPORT HELPERS
# ============================================================

def dataframe_to_pdf_bytes(
    df: pd.DataFrame,
    title: str,
    subtitle: str = "",
) -> bytes:
    """Create a readable PDF summary from an existing result table.

    This helper is presentation-only. It does not recalculate, normalize,
    score, rank, or otherwise modify scientific results.
    """

    buffer = BytesIO()

    document = SimpleDocTemplate(
        buffer,
        pagesize=landscape(A4),
        rightMargin=10 * mm,
        leftMargin=10 * mm,
        topMargin=10 * mm,
        bottomMargin=10 * mm,
        title=title,
        author="BLC Mark",
    )

    styles = getSampleStyleSheet()

    title_style = ParagraphStyle(
        "BLCMarkTitle",
        parent=styles["Title"],
        fontName="Helvetica-Bold",
        fontSize=17,
        leading=20,
        spaceAfter=6,
    )

    subtitle_style = ParagraphStyle(
        "BLCMarkSubtitle",
        parent=styles["BodyText"],
        fontName="Helvetica",
        fontSize=8.5,
        leading=11,
        textColor=colors.HexColor("#4B5563"),
        spaceAfter=8,
    )

    cell_style = ParagraphStyle(
        "BLCMarkCell",
        parent=styles["BodyText"],
        fontName="Helvetica",
        fontSize=6.5,
        leading=8,
        wordWrap="CJK",
    )

    header_style = ParagraphStyle(
        "BLCMarkHeader",
        parent=cell_style,
        fontName="Helvetica-Bold",
        textColor=colors.white,
    )

    story = [Paragraph(title, title_style)]

    if subtitle:
        story.append(Paragraph(subtitle, subtitle_style))

    story.append(
        Paragraph(
            "BLC Mark presents computationally prioritized research candidates. "
            "These outputs should not be interpreted as independently clinically "
            "validated biomarkers or medical recommendations.",
            subtitle_style,
        )
    )
    story.append(Spacer(1, 4))

    pdf_df = df.copy()

    preferred_columns = [
        "gene_id",
        "cancer_cohort",
        "rank",
        "final_score",
        "effect_size",
        "adjusted_p_value",
        "clinical_category",
        "cross_cancer_cohort_count",
        "pathway_count",
    ]

    selected_columns = [
        column
        for column in preferred_columns
        if column in pdf_df.columns
    ]

    # If a table uses slightly different column names, keep the PDF useful
    # without changing the underlying CSV or scientific data.
    if not selected_columns:
        selected_columns = list(pdf_df.columns[:9])

    pdf_df = pdf_df[selected_columns].copy()
    pdf_df = pdf_df.where(pd.notna(pdf_df), "Unavailable")

    def format_pdf_value(value) -> str:
        if isinstance(value, float):
            if abs(value) != 0 and (abs(value) < 0.0001 or abs(value) >= 10000):
                return f"{value:.4e}"
            return f"{value:.6g}"
        return str(value)

    table_data = [
        [Paragraph(str(column), header_style) for column in pdf_df.columns]
    ]

    for _, row in pdf_df.iterrows():
        table_data.append(
            [
                Paragraph(
                    html.escape(format_pdf_value(value)),
                    cell_style,
                )
                for value in row
            ]
        )

    available_width = landscape(A4)[0] - (20 * mm)
    column_width = available_width / max(len(pdf_df.columns), 1)

    table = Table(
        table_data,
        colWidths=[column_width] * len(pdf_df.columns),
        repeatRows=1,
        hAlign="LEFT",
    )

    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#16213E")),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#CBD5E1")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 3),
                ("RIGHTPADDING", (0, 0), (-1, -1), 3),
                ("TOPPADDING", (0, 0), (-1, -1), 3),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
                (
                    "ROWBACKGROUNDS",
                    (0, 1),
                    (-1, -1),
                    [colors.white, colors.HexColor("#F8FAFC")],
                ),
            ]
        )
    )

    story.append(table)
    document.build(story)

    buffer.seek(0)
    return buffer.getvalue()


# ============================================================
# DOWNLOAD DATA
# ============================================================

def page_download_data() -> None:

    render_html(
        """
        <div class="info-hero">
            <h1>Download Data</h1>
            <p>
                Download BLC Mark scientific outputs directly from
                the local BLC Mark release. Download operations do not
                modify the underlying files.
            </p>
        </div>
        """
    )

    section_heading(
        "Top-25 Candidate Tables",
        "Download each cohort table as CSV or as a human-readable PDF summary.",
    )

    for cohort, config in COHORTS.items():

        path = PHASE6_TABLE_DIR / config["top25"]

        left, csv_col, pdf_col = st.columns(
            [4.2, 1.25, 1.25],
            vertical_alignment="center",
        )

        with left:
            render_html(
                f"""
                <div class="download-card">
                    <div class="file-title">{safe(config["name"])}</div>
                    <div class="file-name">{safe(path.name)}</div>
                </div>
                """
            )

        if path.exists():
            with csv_col:
                st.download_button(
                    "Download CSV",
                    data=path.read_bytes(),
                    file_name=path.name,
                    mime="text/csv",
                    key=f"download_csv_{cohort}",
                    use_container_width=True,
                )

            with pdf_col:
                cohort_df = pd.read_csv(path)
                pdf_bytes = dataframe_to_pdf_bytes(
                    cohort_df,
                    f"BLC Mark - {config['name']} Top-25 Prioritized Candidates",
                    f"Cohort: {cohort}",
                )

                st.download_button(
                    "Download PDF",
                    data=pdf_bytes,
                    file_name=path.with_suffix(".pdf").name,
                    mime="application/pdf",
                    key=f"download_pdf_{cohort}",
                    use_container_width=True,
                )
        else:
            with csv_col:
                st.caption("File unavailable")
            with pdf_col:
                st.caption("PDF unavailable")

    combined_path = (
        PHASE6_TABLE_DIR
        / "combined_top_25_biomarkers.csv"
    )

    section_heading(
        "Combined Top-25 Table",
        "Combined Top-25 results across all three cancer cohorts.",
    )

    left, csv_col, pdf_col = st.columns(
        [4.2, 1.25, 1.25],
        vertical_alignment="center",
    )

    with left:
        render_html(
            f"""
            <div class="download-card">
                <div class="file-title">Combined Top-25 Table</div>
                <div class="file-name">{safe(combined_path.name)}</div>
            </div>
            """
        )

    if combined_path.exists():
        with csv_col:
            st.download_button(
                "Download CSV",
                data=combined_path.read_bytes(),
                file_name=combined_path.name,
                mime="text/csv",
                key="download_combined_csv",
                use_container_width=True,
            )

        with pdf_col:
            combined_df = pd.read_csv(combined_path)
            combined_pdf = dataframe_to_pdf_bytes(
                combined_df,
                "BLC Mark - Combined Top-25 Prioritized Candidates",
                "TCGA-BRCA, TCGA-LUAD, and TCGA-COAD",
            )

            st.download_button(
                "Download PDF",
                data=combined_pdf,
                file_name="combined_top_25_biomarkers.pdf",
                mime="application/pdf",
                key="download_combined_pdf",
                use_container_width=True,
            )
    else:
        with csv_col:
            st.caption("File unavailable")
        with pdf_col:
            st.caption("PDF unavailable")

    report_path = (
        PHASE6_REPORT_DIR
        / "BLC_Mark_V1_scientific_report.md"
    )

    section_heading(
        "Scientific Report",
        "Generated Version 1 scientific report.",
    )

    if report_path.exists():
        st.download_button(
            "Download Scientific Report",
            data=report_path.read_bytes(),
            file_name=report_path.name,
            mime="text/markdown",
            use_container_width=False,
        )
    else:
        st.caption("Scientific report unavailable")

    manifest_path = (
        PHASE6_REPRO_DIR
        / "phase6_reproducibility_manifest.json"
    )

    section_heading(
        "Reproducibility Manifest",
        "Release and provenance metadata.",
    )

    if manifest_path.exists():
        st.download_button(
            "Download Reproducibility Manifest",
            data=manifest_path.read_bytes(),
            file_name=manifest_path.name,
            mime="application/json",
            use_container_width=False,
        )
    else:
        st.caption("Reproducibility manifest unavailable")

# ============================================================
# PAGE ROUTER
# ============================================================

if page == "Overview":
    page_overview()

elif page == "Biomarker Explorer":
    page_biomarker_explorer()

elif page == "Compare Cohorts":
    page_compare_cohorts()

elif page == "Methodology":
    page_methodology()

elif page == "About BLC Mark":
    page_about()

elif page == "Download Data":
    page_download_data()