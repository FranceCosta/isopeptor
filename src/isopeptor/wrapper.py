#! /usr/env/python
# -*- coding: utf-8 -*-

import subprocess
from pathlib import Path

def run_jess(input_file):
    jess_exec = Path(__file__).parent / "resources" / "jess" / "jess"
    if not jess_exec.exists():
        raise FileNotFoundError("Jess executable not found.")
    
    # Run Jess with subprocess
    result = subprocess.run([str(jess_exec), input_file], capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"Jess error: {result.stderr}")
    
    return result.stdout