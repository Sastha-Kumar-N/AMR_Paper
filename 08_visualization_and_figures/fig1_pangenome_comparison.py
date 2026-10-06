#!/usr/bin/env python3
# Figure: Figure 1; Fig. 4 callout; panels a
# Purpose: Frequency classes across cohorts
# Inputs: config.yaml, checksum-verified files listed there.
# Tools/versions: ../requirements.txt; date: 2026-10-04.
from pathlib import Path
import sys
sys.dont_write_bytecode = True
HERE=Path(__file__).resolve().parent
# Works both in canonical script folder and self-contained figure copy.
COMMON=HERE.parent if (HERE.parent/"common").is_dir() else HERE
sys.path.insert(0,str(COMMON))
from common.plotters import run
if __name__ == "__main__": run(HERE/"config.yaml", Path(__file__))
