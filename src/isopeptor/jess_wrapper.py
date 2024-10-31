#! /usr/env/python
# -*- coding: utf-8 -*-

import subprocess
from pathlib import Path

def run_jess(pdb_dir : str, distance : float) -> str:
    """
    
        Runs Jess using stored isopeptide bond templates. Only intrachain matches are considered. In case of
        multiple models (delimited by ENDMDL) only the first is considered.

        Args:
        - pdb_dir: directory containing .pdb files
        - distance: used to set rmsd_threshold, distance_cutoff, max_allowed_distance jess parameters
                               higher values will cause a more permissive search. Note that this does not influence
                               probability prediction

        Returns:
        - jess stdout

        Rises:
        - FileNotFoundError if jess executable, template or pdb files are not found
        - RuntimeError if jess fails

    """
    jess_exec = Path(__file__).parent / "resources" / "jess" / "jess"
    isopep_templates = Path(__file__).parent / "resources" / "data" / "templates"
    jess_exec_dir = Path(__file__).parent / "resources" / "data"

    if not jess_exec.exists():
        raise FileNotFoundError("Jess executable not found.")

    if not isopep_templates.exists():
        raise FileNotFoundError("Jess templates not found.")

    pdb_files = [str(p.resolve()) for p in Path(pdb_dir).glob("*.pdb")]
    if not pdb_files:
        raise FileNotFoundError("No PDB files found in the specified directory.")

    rmsd_threshold, distance_cutoff, max_allowed_distance = [str(distance)]*3
    cmd = [str(jess_exec), str(isopep_templates), "-", rmsd_threshold, distance_cutoff, max_allowed_distance]
    
    process = subprocess.Popen(
        cmd, cwd=jess_exec_dir, text=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = process.communicate(input="\n".join(pdb_files))
    
    if process.returncode != 0:
        raise RuntimeError(f"Jess error: {result.stderr}")
    
    return stdout