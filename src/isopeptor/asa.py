#! /usr/env/python
# -*- coding: utf-8 -*-
import biotite.structure.io.pdb as pdb
import biotite.structure as struc
import os
import numpy as np

MAX_ASA = { "rost_sander": { "LYS": 205, "ASP": 163, "GLU": 194, "ASN": 157}}

def calculate_asa(pdb_file_path : str, residues : list, chain : str) -> float:
    """

        Calculate rASA

        Args:
        - pdb_file_path
        - residues: residues involved in isopeptide bond
        - chain

        Rises:
        - FileNotFoundError if pdb_file_path not found
        - ValueError if the amino acid is among expected ones

    """
    if not os.path.isfile(pdb_file_path):
        raise FileNotFoundError("PDB file not found.")
    pdb_file = pdb.PDBFile.read(pdb_file_path)
    # Load structure
    structure = struc.array([atom for atom in pdb_file.get_structure()[0]])
    # Exlude hetero atoms, hydrogens, and select proper chain 
    structure = structure[structure.hetero==False and structure.element != "H" and strcture.chain == chain]
    # Calc sasa
    structure_sasa = struc.sasa(structure, point_number=500)
    
    tot_r_asa = 0
    for res_id in residues:
        res_indeces = [i for i, atom in enumerate(structure) if atom.res_id == res_id]
        res_name = structure[res_indeces[0]].res_name
        if res_name not in MAX_ASA["rost_sander"].keys():
            raise ValueError(f"{res_name} not in {MAX_ASA['rost_sander'].keys()}")

        r_asa = sum([structure_sasa[i] for i in res_indeces]) / MAX_ASA["rost_sander"][res_name]
        tot_r_asa += r_asa

    return tot_r_asa / 3