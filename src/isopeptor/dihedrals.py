#! /usr/env/python
# -*- coding: utf-8 -*-

import biotite.structure as struc
import numpy as np
from isopeptor.constants import DIHEDRAL_ANGLES_TRHESHOLDS
from pathlib import Path
import warnings
import os
import joblib

KDE = None
def get_dihedral_distrib_models():
    global KDE
    bond_types = ["CnaA-like", "CnaB-like"]
    bond_pairs = ["phi_psi", "omega_phi", "omega_psi"]
    if KDE is None:
        KDE = {bond_type:
                {bond_pair: None for bond_pair in bond_pairs
                } for bond_type in bond_types}

        for bond_type in bond_types:
            for bond_pair in bond_pairs:
                model_path = Path(__file__).parent / "resources" / "models" / f"{bond_type}_{bond_pair}.pkl"
                if not os.path.isfile(model_path):
                    raise ValueError(f"Model file not found in {model_path}.")
                KDE[bond_type][bond_pair] = joblib.load(model_path)

    return KDE

def get_dihedral_angles(structure:struc.AtomArray, chain:str, r1_bond:int, 
                    r2_bond:int, r2_bond_name:str) -> tuple:
    """

        Get dihedral torsion angles. Returns pseudo_omega, pseudo_psi and pseudo_phi in a tuple

        Args:
        - structure: biotite.structure.AtomArray
        - chain: str
        - r1_bond: int
        - r2_bond: int
        - r2_bond_name: str

        Returns:
        - (pseudo_omega, pseudo_psi, pseudo_phi)

    """
    
    lys_cd = [atom for atom in structure if atom.res_id == r1_bond and atom.chain_id == chain and atom.atom_name == "CD"]
    lys_ce = [atom for atom in structure if atom.res_id == r1_bond and atom.chain_id == chain and atom.atom_name == "CE"]
    lys_nz = [atom for atom in structure if atom.res_id == r1_bond and atom.chain_id == chain and atom.atom_name == "NZ"]
    asn_nz = [atom for atom in structure if atom.res_id == r2_bond and atom.chain_id == chain and atom.atom_name == "ND2"]

    if len(lys_cd) > 0:
        lys_cd = lys_cd[0]
    else:
        lys_cd = None

    if len(lys_ce) > 0:
        lys_ce = lys_ce[0]
    else:
        lys_ce = None
    
    if len(lys_nz) > 0:
        n_atom = lys_nz[0]
    elif len(asn_nz) > 0:
        n_atom = asn_nz[0]
    else:
        n_atom = None

    if r2_bond_name == "ASN" or r2_bond_name == "ASP":
        # CG is the last
        c1 = [atom for atom in structure if atom.res_id == r2_bond and atom.chain_id == chain and atom.atom_name == "CG"]
        if len(c1) > 0:
            c1 = c1[0]
        else: 
            c1 = None
        # CB
        c2 = [atom for atom in structure if atom.res_id == r2_bond and atom.chain_id == chain and atom.atom_name == "CB"]
        if len(c2) > 0:
            c2 = c2[0]
        else:
            c2 = None
        # CA
        c3 = [atom for atom in structure if atom.res_id == r2_bond and atom.chain_id == chain and atom.atom_name == "CA"]
        if len(c3) > 0:
            c3 = c3[0]
        else:
            c3 = None

    if r2_bond_name == "GLU":
        # CD is the last
        c1 = [atom for atom in structure if atom.res_id == r2_bond and atom.chain_id == chain and atom.atom_name == "CD"]
        if len(c1) > 0:
            c1 = c1[0]
        else:
            c1 = None
        # CG
        c2 = [atom for atom in structure if atom.res_id == r2_bond and atom.chain_id == chain and atom.atom_name == "CG"]
        if len(c2) > 0:
            c2 = c2[0]
        else:
            c2 = None
        # CB
        c3 = [atom for atom in structure if atom.res_id == r2_bond and atom.chain_id == chain and atom.atom_name == "CB"]
        if len(c3) > 0:
            c3 = c3[0]
        else:
            c3 = None

    if any([True if atom == None else False for atom in [lys_cd, lys_ce, n_atom, c1, c2, c3]]):
        pseudo_omega, pseudo_psi, pseudo_phi = [None]*3
        warnings.warn("Missing atom from bond. Dihedrals not calculated", UserWarning)
        #print([lys_cd, lys_ce, n_atom, c1, c2, c3])
    else:
        # This corresponds to the peptide omega angle (torsion angle on bond CN)
        pseudo_omega = struc.dihedral(c2, c1, n_atom, lys_ce)
        pseudo_omega = round(pseudo_omega*180/np.pi, 3)
        # This corresponds to psi (one carbon is used instead of a second N)
        pseudo_psi = struc.dihedral(c3, c2, c1, n_atom)
        pseudo_psi = round(pseudo_psi*180/np.pi, 3)
        # This corresponds to phi (lys_cd should be bound to oxygen and then to Nitrogrn to form the next peptide bond)
        pseudo_phi = struc.dihedral(c1, n_atom, lys_ce, lys_cd)
        pseudo_phi = round(pseudo_phi*180/np.pi, 3)
    
    return (pseudo_omega, pseudo_psi, pseudo_phi)

def get_dihedral_angles_stats(structure:struc.AtomArray, chain:str,
                                r1_bond:int, r2_bond:int, r2_bond_name:str, bond_type:str) -> tuple:
    """

    Get dihedral torsion angles and stats about their distribution

    Args:
        structure (biotite.structure.AtomArray)
        chain (str)
        r1_bond (int)
        r2_bond (int)
        r2_bond_name (str)
        bond_type (str)

    Returns:
            (
                pseudo_omega:float, pseudo_phi:float, pseudo_psi:float, phi_psi_likelihood:float, omega_psi_likelihood:float, 
                omega_phi_likelihood:float, phi_psi_allowed:bool, omega_psi_allowed:bool, omega_phi_allowed:bool
            )

    """
    pseudo_omega, pseudo_psi, pseudo_phi = get_dihedral_angles(structure, chain, r1_bond, 
                    r2_bond, r2_bond_name)
    if pseudo_omega == None or pseudo_psi == None or pseudo_phi == None:
        phi_psi_likelihood, omega_psi_likelihood, omega_phi_likelihood, phi_psi_allowed, omega_psi_allowed, omega_phi_allowed = [None]*6
    else:
        # Load models
        kde = get_dihedral_distrib_models()

        # Get likelihoods
        phi_psi_likelihood = round(kde[bond_type]["phi_psi"].score([[pseudo_phi, pseudo_psi]]), 3)
        omega_psi_likelihood = round(kde[bond_type]["omega_psi"].score([[pseudo_omega, pseudo_psi]]), 3)
        omega_phi_likelihood = round(kde[bond_type]["omega_phi"].score([[pseudo_omega, pseudo_phi]]), 3)

        phi_psi_allowed, omega_psi_allowed, omega_phi_allowed = [False]*3
        
        if phi_psi_likelihood > DIHEDRAL_ANGLES_TRHESHOLDS[bond_type]["phi_psi"]:
            phi_psi_allowed = True

        if omega_psi_likelihood > DIHEDRAL_ANGLES_TRHESHOLDS[bond_type]["omega_psi"]:
            omega_psi_allowed = True

        if omega_phi_likelihood > DIHEDRAL_ANGLES_TRHESHOLDS[bond_type]["omega_phi"]:
            omega_phi_allowed = True

    return (
                pseudo_omega, pseudo_phi, pseudo_psi, phi_psi_likelihood, omega_psi_likelihood, omega_phi_likelihood, 
                phi_psi_allowed, omega_psi_allowed, omega_phi_allowed
            )

