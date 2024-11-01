#! /usr/env/python
# -*- coding: utf-8 -*-

import biotite.structure.io.pdb as pdb
import biotite.structure as struc
import os

def _get_structure_asa(pdb_file_path: str) -> tuple:
    """

        Calculate structure asa

        Args:
        - pdb_file_path

        Rises:
        - FileNotFoundError if pdb_file_path not found

    """
    if not os.path.isfile(pdb_file_path):
        raise FileNotFoundError("PDB file not found.")
    pdb_file = pdb.PDBFile.read(pdb_file_path)
    # Load structure
    structure = struc.array([atom for atom in pdb_file.get_structure()[0]])
    # Exlude hetero atoms, hydrogens, and select proper chain 
    structure = structure[(structure.hetero==False) & (structure.element != "H")]
    # Calc sasa
    structure_sasa = struc.sasa(structure, point_number=500)
    
    return (structure_sasa, structure)