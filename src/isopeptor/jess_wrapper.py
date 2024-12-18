#! /usr/env/python
# -*- coding: utf-8 -*-

from pyjess import Template
from pyjess import Jess
from pyjess import Molecule
from pathlib import Path
import os
import biotite.structure.io.pdbx as pdbx
import biotite.structure.io.pdb as pdb
import warnings

def run_jess(structure_files: str, distance: float, 
            templates = Path(__file__).parent / "resources" / "data" / "template_structures") -> list:
    """
    
        Runs Jess using stored isopeptide bond templates (default templates). 
        Only intrachain matches are considered. In case of multiple models 
        (delimited by ENDMDL) only the first is considered.
        If .cif files are found they are converted into .pdb models (and stored in "./converted_models/").

        Args:
            structure_files (list): containing paths to .cif and .pdb files
            distance (float): used to set rmsd_threshold, distance_cutoff, max_allowed_distance jess parameters
                               higher values will cause a more permissive search. Note that this does not influence
                               probability prediction
            templates (str): directory containing templates to be used with .pdb extension

        Returns:
            list of pyjess._jess.Hit 

        Rises:
            FileNotFoundError if template or pdb/cif files are not found

    """
    template_files = [str(p.resolve()) for p in Path(templates).glob("*.pdb")]
    if not template_files:
        raise FileNotFoundError("Jess templates not found.")

    pdb_files = [p for p in structure_files if p[-4:]==".pdb"]
    cif_files = [p for p in structure_files if p[-4:]==".cif"]

    if not pdb_files+cif_files:
        raise FileNotFoundError("No PDB/CIF files found in the specified directory.")
    
    # Convert
    if cif_files:
        warnings.warn("CIF files detected. Converting them into PDB format.", UserWarning)
        converted_models_dir = os.path.join(os.getcwd(), "converted_models")
        os.makedirs(converted_models_dir, exist_ok=True)
        for cif_file_path in cif_files:
            cif_file = pdbx.CIFFile.read(cif_file_path)
            arr = pdbx.get_structure(cif_file)
            pdb_file = pdb.PDBFile()
            pdb_file.set_structure(arr)
            pdb_file_path = os.path.join(converted_models_dir, os.path.basename(cif_file_path).replace(".cif", ".pdb"))
            pdb_file.write(pdb_file_path)
            pdb_files.append(pdb_file_path)

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