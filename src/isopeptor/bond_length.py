#! /usr/env/python
# -*- coding: utf-8 -*-

import biotite.structure as struc
from isopeptor.constants import MEAN_BOND_LENGTH
from isopeptor.constants import STD_DEV_BOND_LENGTH

def get_bond_length(structure:struc.AtomArray, chain:str, r1_bond:int, 
                    r2_bond:int, r2_bond_name:str):
    """

        Get bond length

        Args:
        - structure: biotite.structure.AtomArray
        - chain: str
        - r1_bond: int
        - r2_bond: int
        - r2_bond_name: str

    """
    lys_nz = [atom for atom in structure if atom.res_id == r1_bond and atom.chain_id == chain and atom.atom_name == "NZ"][0]
    if r2_bond_name == "N" or r2_bond_name == "D":
        # CG is the last
        c1 = [atom for atom in structure if atom.res_id == r2_bond and atom.chain_id == chain and atom.atom_name == "CG"][0]
    if r2_bond_name == "E":
        # CD is the last
        c1 = [atom for atom in structure if atom.res_id == r2_bond and atom.chain_id == chain and atom.atom_name == "CD"][0]

    bond_length = struc.distance(lys_nz, c1)

    return bond_length

def get_bond_zscore(structure:struc.AtomArray, chain:str, r1_bond:int, 
                    r2_bond:int, r2_bond_name:str) -> tuple:
    """

        Get bond zscore

        Args:
        - structure: biotite.structure.AtomArray
        - chain: str
        - r1_bond: int
        - r2_bond: int
        - r2_bond_name: str

        Returns:
        - (bond_length:float, zscore:float, bond_length_allowed:bool)

    """
    # Get bond length
    bond_length = get_bond_length(structure, chain, r1_bond, r2_bond, r2_bond_name)
    zscore = (bond_length - MEAN_BOND_LENGTH)/STD_DEV_BOND_LENGTH
    bond_length_allowed = False
    if zscore < 4:
        bond_length_allowed = True

    return (bond_length, zscore, bond_length_allowed)

