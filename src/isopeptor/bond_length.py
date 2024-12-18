#! /usr/env/python
# -*- coding: utf-8 -*-

import biotite.structure as struc
from isopeptor.constants import MEAN_BOND_LENGTH
from isopeptor.constants import STD_DEV_BOND_LENGTH
import warnings
import numpy as np

def get_bond_length(structure:struc.AtomArray, chain:str, r1_bond:int, 
                    r2_bond:int, r2_bond_name:str):
    """

        Get bond length

        Args:
            structure (biotite.structure.AtomArray)
            chain (str)
            r1_bond (int)
            r2_bond (int)
            r2_bond_name (str)

    """
    bond_length = None
    n_atom = None
    c_atom = None

    # Get N atom (from lys or asn)
    n_res  = [atom for atom in structure if atom.res_id == r1_bond and atom.chain_id == chain and atom.atom_name == "NZ"]
    if len(n_res) == 0:
        n_res = [atom for atom in structure if atom.res_id == r2_bond and atom.chain_id == chain and atom.atom_name == "ND2"]
    if len(n_res) == 0:
        warnings.warn("Missing nitrogen from bond. Bond length not calculated", UserWarning)
    else:
        n_atom = n_res[0]

    if r2_bond_name == "ASN" or r2_bond_name == "ASP":
        # CG is the last
        c_res = [atom for atom in structure if atom.res_id == r2_bond and atom.chain_id == chain and atom.atom_name == "CG"]

    if r2_bond_name == "GLU":
        # CD is the last
        c_res = [atom for atom in structure if atom.res_id == r2_bond and atom.chain_id == chain and atom.atom_name == "CD"]
    
    if len(c_res) == 0:
        warnings.warn("Missing carbon from bond. Bond length not calculated", UserWarning)
    else:
        c_atom = c_res[0]
    
    if n_atom != None and c_atom != None:
        bond_length = round(struc.distance(n_atom, c_atom), 3)

    return bond_length

def get_bond_stats(structure:struc.AtomArray, chain:str, r1_bond:int, 
                    r2_bond:int, r2_bond_name:str) -> tuple:
    """

        Get bond zscore

        Args:
            structure (biotite.structure.AtomArray)
            chain (str)
            r1_bond (int)
            r2_bond (int)
            r2_bond_name (str)

        Returns:
            (bond_length:float, bond_length_zscore:float, bond_length_allowed:bool)

    """
    # Get bond length
    bond_length = get_bond_length(structure, chain, r1_bond, r2_bond, r2_bond_name)
    bond_length_zscore = None
    if bond_length != None:
        bond_length_zscore = round((bond_length - MEAN_BOND_LENGTH)/STD_DEV_BOND_LENGTH, 3)
        bond_length_allowed = False
        if np.abs(bond_length_zscore) < 4:
            bond_length_allowed = True

    return (bond_length, bond_length_zscore, bond_length_allowed)

