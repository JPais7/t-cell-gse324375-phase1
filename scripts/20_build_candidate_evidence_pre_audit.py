#!/usr/bin/env python3
"""Single upstream construction entry point for candidate evidence."""
import os, runpy
os.environ['PHASE1_BUILD_PRE_ONLY']='1'
runpy.run_path(os.path.join(os.path.dirname(__file__), '20_final_candidate_ranking.py'), run_name='__main__')
