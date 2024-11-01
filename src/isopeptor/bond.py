#! /usr/env/python
# -*- coding: utf-8 -*-

class BondElement:
    """

        Stores isopeptide bond predictions

        Attributes:
            pdb_file: path to pdb file
            protein_name: protein name from file name
            rmsd: rmsd with template (in A)
            template: template name
            chain: protein chain
            r1_bond: residue number of residue 1 involved in bond
            r_cat: residue number of catalytic residue
            r2_bond: residue number of residue 2 involved in bond

    """
    def __init__(self, pdb_file: str, protein_name: str, rmsd: float, template: str, chain: str, r1_bond: int, r_cat: int, r2_bond: int):
        self.pdb_file = pdb_file
        self.protein_name = protein_name
        self.rmsd = rmsd
        self.template = template
        self.chain = chain
        self.r1_bond = r1_bond
        self.r_cat = r_cat
        self.r2_bond = r2_bond
        self.r_asa: float | None = None
        self.probability: float | None = None

    def __repr__(self):
        s = f"BondElement(pdb_file={self.pdb_file}, protein_name={self.protein_name}, rmsd={self.rmsd}, template={self.template}, "+\
            f"chain={self.chain}, r1_bond={self.r1_bond}, r_cat={self.r_cat}, r2_bond={self.r2_bond}, r_asa={self.r_asa}, "+\
            f"probability={self.probability})"
        return s
