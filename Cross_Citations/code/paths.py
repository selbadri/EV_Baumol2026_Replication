"""
==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    paths.py
 Purpose: Every input and output path of the Cross_Citations step, relative to the
          Cross_Citations folder. All scripts import this file, so the package runs on
          another machine without editing any path.
 ------------------------------------------------------------------------------------------
 Inputs:  none (definitions only)
 Outputs: none
==========================================================================================
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]          # Cross_Citations/
PACKAGE = ROOT.parent                                # the replication package
EMPIRICAL = PACKAGE / "Empirical"

INPUT = ROOT / "input"
INTERMEDIATE = ROOT / "intermediate"                 # extracted archives, not tracked
OUTPUT = ROOT / "output"

# ------------------------------------------------------------------------- inputs
# Concordances (tracked)
CONCORDANCES = INPUT / "concordances"
ALP_IPC = CONCORDANCES / "ipc4_to_isic_rev31_4.txt"            # IPC4 -> ISIC Rev 3.1
SECTOR_MAP = CONCORDANCES / "isic31_4digit_to_sector.csv"       # ISIC Rev 3.1 -> sector
ALP_USPC_ISIC = CONCORDANCES / "uspc_Class_to_isic_rev31_4.txt"  # USPC -> ISIC (robustness)

# USPTO Historical Patent Data Files: the annual counts are tracked here, the patent-level
# masterfile and the USPC -> NAICS crosswalk are shared with ../Empirical/input/
HPDF_ANNUAL = INPUT / "USPTO_HPDF" / "annual.csv"
HPDF_MASTER = EMPIRICAL / "input" / "historical_masterfile.dta"
ALP_USPC_NAICS = EMPIRICAL / "input" / "ALP_USPC" / "NAICS_2007" / "uspc_Class_to_naics07_2.txt"

# Large downloads (not tracked; see README)
PATENTSVIEW = INPUT / "PatentsView"
OECD_CITATIONS = INPUT / "OECD_Citations_Feb2025"
OECD_TPF = INPUT / "OECD_TPF_Feb2025"
OECD_QUALITY = INPUT / "OECD_Quality_Feb2025"
KPSS = INPUT / "KPSS"

# ------------------------------------------------------------------------- caches
CITATION_CELLS = OUTPUT / "citations" / "cell_counts_ipc4_exbm.parquet"   # step 06 cache

# ------------------------------------------------------------------------ outputs
OUT_PATENTS = OUTPUT / "patents"           # steps 01-05
OUT_CITATIONS = OUTPUT / "citations"       # step 06
OUT_ROBUSTNESS = OUTPUT / "robustness"     # code/robustness/
LOGS = OUTPUT / "logs"


def pv_file(stem):
    """A PatentsView table, zipped or not. Exits with a pointer to the README if absent."""
    for ext in (".tsv.zip", ".tsv"):
        p = PATENTSVIEW / (stem + ext)
        if p.exists():
            return p
    raise SystemExit(f"Missing {PATENTSVIEW / (stem + '.tsv.zip')}. "
                     "See 'Data availability' in the README for the download.")


def need(path):
    """Return path if it exists, otherwise exit with a pointer to the README."""
    if not Path(path).exists():
        raise SystemExit(f"Missing {path}. See the README for where it comes from.")
    return Path(path)
