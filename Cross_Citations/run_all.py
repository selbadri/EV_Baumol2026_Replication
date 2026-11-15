"""
==========================================================================================
 An Anatomy of the Baumol Cost Disease -- Elbadri & Vizcaino (2026)
 ------------------------------------------------------------------------------------------
 File:    run_all.py
 Purpose: Run the whole Cross_Citations step in order:
            01_business_methods.py        business-method flag (current CPC G06Q)
            02_class_counts.py            patents by IPC subclass and year
            03_partition_compare.py       ISIC vs NAICS sector partition check
            04_hpdf_sector_flows.py       USPC-route sector flows, 1790-2014
            05_patent_stocks.py           patent shares and knowledge stocks
            06_citation_probabilities.py  cross-sector citation probabilities
          and, with --robustness, the alternative patent stocks in code/robustness/.
          Each step runs in its own Python process; output also goes to output/logs/.
 ------------------------------------------------------------------------------------------
 Usage:   python run_all.py                  everything (needs the large downloads in input/)
          python run_all.py --quick          skip every pass over the large downloads:
                                             steps 03 and 05 from the tracked class counts and
                                             USPC-route flows, step 06 from the tracked citation
                                             cache. Needs nothing outside the repository
          python run_all.py --robustness     also run code/robustness/r1-r6
 ------------------------------------------------------------------------------------------
 Inputs:  see the README and code/paths.py
 Outputs: output/patents/, output/citations/, output/robustness/, output/logs/
==========================================================================================
"""

import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CODE = ROOT / "code"
LOGS = ROOT / "output" / "logs"

QUICK = "--quick" in sys.argv
ROBUSTNESS = "--robustness" in sys.argv

# (script, arguments, needs the large downloads)
STEPS = [
    ("01_business_methods.py", [], True),
    ("02_class_counts.py", [], True),
    ("03_partition_compare.py", [], False),
    ("04_hpdf_sector_flows.py", [], True),     # reads the 800 MB USPTO masterfile
    ("05_patent_stocks.py", [], False),
    ("06_citation_probabilities.py", ["--from-cache"] if QUICK else [], False),
]
ROBUSTNESS_STEPS = [
    ("robustness/r1_concordance_overlap.py", [], True),
    ("robustness/r2_ipc_sector_flows.py", [], True),
    ("robustness/r3_triadic_flows.py", [], True),
    ("robustness/r4_quality_flows.py", [], True),
    ("robustness/r5_kpss_flows.py", [], True),
    ("robustness/r5_kpss_flows.py",
     ["--drop", str(ROOT / "output" / "patents" / "bm_ids_cpc.csv"), "--tag", "_exbm"], True),
    ("robustness/r6_stock_variants.py", [], False),
]


def run(script, args):
    name = Path(script).stem + ("_exbm" if "--tag" in args else "")
    log = LOGS / f"{name}.log"
    print(f"\n=== {script} {' '.join(args)}", flush=True)
    t0 = time.time()
    with open(log, "w") as fh:
        proc = subprocess.Popen([sys.executable, str(CODE / script), *args], cwd=ROOT,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in proc.stdout:
            sys.stdout.write(line)
            fh.write(line)
        proc.wait()
    if proc.returncode != 0:
        sys.exit(f"\n{script} failed (exit {proc.returncode}); see {log}")
    print(f"=== done in {time.time() - t0:.0f}s, log {log.relative_to(ROOT)}")


def main():
    LOGS.mkdir(parents=True, exist_ok=True)
    steps = STEPS + (ROBUSTNESS_STEPS if ROBUSTNESS else [])
    for script, args, heavy in steps:
        if QUICK and heavy:
            print(f"\n[skip] {script} (--quick: reads the large downloads)")
            continue
        run(script, args)
    print("\nCross_Citations finished.")


if __name__ == "__main__":
    main()
