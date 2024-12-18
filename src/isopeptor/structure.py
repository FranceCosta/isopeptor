#! /usr/env/python
# -*- coding: utf-8 -*-

import biotite.structure.io.pdb as pdb
import biotite.structure as struc
import os

def get_structure(pdb_file_path: str) -> tuple:
    """

        Load structure and remove hydrogens and hetero atoms

        Args:
            pdb_file_path (str): path to input pdb file

        Rises:
             FileNotFoundError if pdb_file_path not found

    """
    if not os.path.isfile(pdb_file_path):
        raise FileNotFoundError("PDB file not found.")
    pdb_file = pdb.PDBFile.read(pdb_file_path)
    # Load structure
    structure = struc.array([atom for atom in pdb_file.get_structure()[0]])
    # Exlude hetero atoms and hydrogens
    structure = structure[(structure.hetero==False) & (structure.element != "H")]

    return structure