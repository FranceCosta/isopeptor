#! /usr/env/python
# -*- coding: utf-8 -*-

class BondElement:
    """

        Stores isopeptide bond predictions

        Attributes:
            struct_file (str): path to structure file
            protein_name (str): protein name from file name
            rmsd (float): rmsd with template (in A)
            template (str): template name
            chain (str): protein chain
            r1_bond (int): residue number of residue 1 involved in bond
            r_cat (int): residue number of catalytic residue
            r2_bond (int): residue number of residue 2 involved in bond

    """
    def __init__(
        self, struct_file: str, protein_name: str, rmsd: float, template: str,
        chain: str, r1_bond: int, r_cat: int, r2_bond: int,
        r1_bond_name: str, r_cat_name: str, r2_bond_name: str
        ):
        self.struct_file = struct_file
        self.protein_name = protein_name
        self.rmsd = rmsd
        self.template = template
        self.chain = chain
        self.r1_bond = r1_bond
        self.r_cat = r_cat
        self.r2_bond = r2_bond
        self.r1_bond_name = r1_bond_name
        self.r_cat_name = r_cat_name
        self.r2_bond_name = r2_bond_name
        self.r_asa: float | None = None
        self.probability: float | None = None
        self.bond_length: float | None = None
        self.bond_length_zscore: float | None = None
        self.bond_length_allowed: bool | None = None
        self.pseudo_phi: float | None = None
        self.pseudo_psi: float | None = None
        self.pseudo_omega: float | None = None
        self.phi_psi_likelihood: float | None = None
        self.phi_psi_allowed: bool | None = None
        self.omega_psi_likelihood: float | None = None
        self.omega_psi_allowed: bool | None = None
        self.omega_phi_likelihood: float | None = None
        self.omega_phi_allowed: bool | None = None

    def __repr__(self):
        s = f"BondElement(struct_file={self.struct_file}, protein_name={self.protein_name}, rmsd={self.rmsd}, template={self.template}, "+\
            f"chain={self.chain}, r1_bond={self.r1_bond}, r_cat={self.r_cat}, r2_bond={self.r2_bond}, "+\
            f"r1_bond_name={self.r1_bond_name}, r_cat_name={self.r_cat_name}, r2_bond_name={self.r2_bond_name}, "+\
            f"bond_type={self.bond_type}, r_asa={self.r_asa}, probability={self.probability})"
        return s
