#! /usr/env/python
# -*- coding: utf-8 -*-

import biotite.structure.io.pdb as pdb
import biotite.structure as struc
import os
from isopeptor.structure import get_structure

def get_structure_asa(pdb_file_path: str) -> tuple:
    """

        Return structure and ASA

        Args:
            pdb_file_path (str): path to PDB structure

    """
    # Load structure
    structure = get_structure(pdb_file_path)
    # Calc sasa
    structure_sasa = struc.sasa(structure, point_number=500)
    
    return (structure_sasa, structure)