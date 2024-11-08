#! /usr/env/python
# -*- coding: utf-8 -*-

from pyjess import Template
from pyjess import Jess
from pyjess import Molecule
from pathlib import Path
import os

def _run_jess(pdb_dir: str, distance: float, templates = Path(__file__).parent / "resources" / "data" / "template_structures") -> list:
    """
    
        Runs Jess using stored isopeptide bond templates (default templates). 
        Only intrachain matches are considered. In case of multiple models 
        (delimited by ENDMDL) only the first is considered.

        Args:
        - pdb_dir: directory containing .pdb files
        - distance: used to set rmsd_threshold, distance_cutoff, max_allowed_distance jess parameters
                               higher values will cause a more permissive search. Note that this does not influence
                               probability prediction
        - templates: directory containing templates to be used with .pdb extension

        Returns:
        - list of pyjess._jess.Hit 

        Rises:
        - FileNotFoundError if template or pdb files are not found

    """
    template_files = [str(p.resolve()) for p in Path(templates).glob("*.pdb")]
    if not template_files:
        raise FileNotFoundError("Jess templates not found.")

    pdb_files = [str(p) for p in Path(pdb_dir).glob("*.pdb")]
    if not pdb_files:
        raise FileNotFoundError("No PDB files found in the specified directory.")

    rmsd_threshold, distance_cutoff, max_dynamic_distance = [distance]*3
    
    # Load templates
    templates = []
    for path in list(template_files):
        templates.append(Template.load(path, id=os.path.basename(path).split(".pdb")[0]))
    jess = Jess(templates)
    
    # Run on structures
    hits = []
    for path in pdb_files:
        mol = Molecule.load(path, id=path)
        query = jess.query(mol, rmsd_threshold=rmsd_threshold, distance_cutoff=distance_cutoff, max_dynamic_distance=max_dynamic_distance)
        hits.extend(query)

    return hits