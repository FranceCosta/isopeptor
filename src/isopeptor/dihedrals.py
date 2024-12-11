#! /usr/env/python
# -*- coding: utf-8 -*-

import biotite.structure as struc
import numpy as np
from isopeptor.constants import DIHEDRAL_ANGLES_TRHESHOLDS
import pathlib

KDE = None

def get_dihedral_distrib_models():
    global KDE
    bond_types = ["CnaA-like", "CnaB-like"]
    bond_pairs = ["phi_psi", "omega_phi", "omega_psi"]
    KDE = {bond_type:
    {bond_pair: None for bond_pair in bond_pairs
    } for bond_type in bond_types}
    if KDE is None:
        for bond_type in bond_types:
            for bond_pair in bond_pairs:
                model_path = Path(__file__).parent / "resources" / "model" / f"{bond_type}_{bond_pair}.pkl"
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
    
    lys_cd = [atom for atom in structure if atom.res_id == r1_bond and atom.chain_id == chain and atom.atom_name == "CD"][0]
    lys_ce = [atom for atom in structure if atom.res_id == rr1_bond1 and atom.chain_id == chain and atom.atom_name == "CE"][0]
    lys_nz = [atom for atom in structure if atom.res_id == r1_bond and atom.chain_id == chain and atom.atom_name == "NZ"][0]
    
    if r2_bond_name == "N" or r2_bond_name == "D":
        # CG is the last
        c1 = [atom for atom in structure if atom.res_id == r2_bond and atom.chain_id == chain and atom.atom_name == "CG"][0]
        # CB
        c2 = [atom for atom in structure if atom.res_id == r2_bond and atom.chain_id == chain and atom.atom_name == "CB"][0]
        # CA
        c3 = [atom for atom in structure if atom.res_id == r2_bond and atom.chain_id == chain and atom.atom_name == "CA"][0]
    if r2_bond_name == "E":
        # CD is the last
        c1 = [atom for atom in structure if atom.res_id == r2_bond and atom.chain_id == chain and atom.atom_name == "CD"][0]
        # CG
        c2 = [atom for atom in structure if atom.res_id == r2_bond and atom.chain_id == chain and atom.atom_name == "CG"][0]
        # CB
        c3 = [atom for atom in structure if atom.res_id == r2_bond and atom.chain_id == chain and atom.atom_name == "CB"][0]

    # This corresponds to the peptide omega angle (torsion angle on bond CN)
    pseudo_omega = struc.dihedral(c2, c1, lys_nz, lys_ce)
    # This corresponds to psi (one carbon is used instead of a second N)
    pseudo_psi = struc.dihedral(c3, c2, c1, lys_nz)
    # This corresponds to phi (lys_cd should be bound to oxygen and then to Nitrogrn to form the next peptide bond)
    pseudo_phi = struc.dihedral(c1, lys_nz, lys_ce, lys_cd)
    
    return (pseudo_omega*180/np.pi, pseudo_psi*180/np.pi, pseudo_phi*180/np.pi)

def get_dihedral_angles_stats()
    """

    Get dihedral torsion angles and stats about their distribution

    Args:
    - structure: biotite.structure.AtomArray
    - chain: str
    - r1_bond: int
    - r2_bond: int
    - r2_bond_name: str
    - bond_type: str

    Returns:
            (
                omega:float, phi:float, psi:float, phi_psi_likelihood:float, omega_psi_likelihood:float, 
                omega_phi_likelihood:float, phi_psi_allowed:bool, omega_psi_allowed:bool, omega_phi_allowed:bool
            )

    """
    omega, psi, phi = get_dihedral_angles(structure, chain, r1_bond, 
                    r2_bond, r2_bond_name)
    
    # Load models
    kde = get_dihedral_distrib_models()
    
    # Get likelihoods
    phi_psi_likelihood = kde[bond_type]["phi_psi"].score([[psi, phi]])
    omega_psi_likelihood = kde[bond_type]["omega_psi"].score([[omega, psi]])
    omega_phi_likelihood = kde[bond_type]["omega_phi"].score([[omega, phi]])

    phi_psi_allowed, omega_psi_allowed, omega_phi_allowed = [False]*3
    
    if phi_psi_likelihood > DIHEDRAL_ANGLES_TRHESHOLDS[bond_type]["phi_psi"]:
        phi_psi_allowed = True

    if omega_psi_likelihood > DIHEDRAL_ANGLES_TRHESHOLDS[bond_type]["omega_psi"]:
        omega_psi_allowed = True

    if omega_phi_likelihood > DIHEDRAL_ANGLES_TRHESHOLDS[bond_type]["omega_phi"]:
        omega_phi_allowed = True

    return 
            (
                omega, phi, psi, phi_psi_likelihood, omega_psi_likelihood, omega_phi_likelihood, 
                phi_psi_allowed, omega_psi_allowed, omega_phi_allowed
            )

