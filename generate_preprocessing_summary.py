#!/usr/bin/env python3
"""
generate_preprocessing_summary.py
==================================
Thin CLI wrapper for spatia.analysis.preprocessing_summary.generate_summary_report().

Every `run_preprocessing()` call (i.e. every `--steps preprocessing`
pipeline run) now generates this same summary_report/ automatically at the
end of the step -- see spatia/analysis/preprocessing.py. Run this script
directly only when you want to cheaply re-render the report WITHOUT
re-running image processing -- e.g. right after a report-code fix, or to
regenerate a report you deleted -- since it just reads
run_preprocessing()'s already-written outputs and never touches raw images.

Usage
-----
    python generate_preprocessing_summary.py --config experiments/crc_tma_full_pipeline.yaml
"""

from __future__ import annotations

import argparse
import os
import sys

import yaml

from spatia.analysis.preprocessing_summary import generate_summary_report


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", required=True, help="Same YAML config used for the pipeline run")
    args = parser.parse_args()

    if not os.path.exists(args.config):
        print(f"Config not found: {args.config}")
        sys.exit(2)

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    generate_summary_report(cfg)


if __name__ == "__main__":
    main()
