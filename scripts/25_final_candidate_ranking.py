#!/usr/bin/env python3
"""Downstream final ranking entry point (consumes pre-audit and robustness outputs)."""
import os, runpy
os.environ.pop('PHASE1_BUILD_PRE_ONLY', None)
runpy.run_path(os.path.join(os.path.dirname(__file__), '20_final_candidate_ranking.py'), run_name='__main__')
