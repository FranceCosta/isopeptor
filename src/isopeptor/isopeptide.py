#! /usr/env/python
# -*- coding: utf-8 -*-

import re
from typing import List
import warnings
import numpy as np
from isopeptor.jess_wrapper import _run_jess
from isopeptor.asa import _get_structure_asa
from isopeptor.bond import BondElement 
from isopeptor.constants import MAX_ASA
from isopeptor.constants import BOND_TYPE
from isopeptor.model import _predict

class Isopeptide:
    """
    
        Handles isopeptide bond prediction running jess via the jess_wrapper module
        and solvent accessible area via the asa module. Stores isopeptide bond predictions
        as a list of BondElement. Prediction is run for all structures from pdb_dir. 
        If multiple matches with an isopeptide bond signature are detected, only the one
        with the lowest RMSD is retained.

        Attributes:
            pdb_dir: str where pdb files are located
            distance: float that specifies permissivity of jess search
            jess_output: None | str that stores jess raw output for debug purposes
            isopeptide_bonds: list that stores isopeptide bonds as BondElement elements
            fixed_r_asa: float | None which ranges between 0 and 1 that fixes r_asa allowing to skip its calculation

        Usage:
            
            Use a fixed solvent accessible area for a quick prediction:
            >>> from isopeptor.isopeptide import Isopeptide
            >>> i = Isopeptide("tests/data/test_structures", distance=1.5, fixed_r_asa=0.1)
            >>> i.predict()
            >>> i.isopeptide_bonds[0]
            BondElement(pdb_file=/nfs/research/agb/research/francesco/projects/20241024_isopeptor_v1/tests/data/test_structures/8beg.pdb, protein_name=8beg, rmsd=0.0, template=8beg_A_590_636_729, chain=A, r1_bond=590, r_cat=636, r2_bond=729, r1_bond_name=LYS, r_cat_name=ASP, r2_bond_name=ASN, bond_type=CnaA-like, r_asa=0.1, probability=0.984)
            >>> i.print_tabular()
            protein_name        chain   r1_bond r_cat   r2_bond r1_bond_name    r_cat_name      r2_bond_name    bond_type       rmsd    r_asa   probability     template
            8beg                A       590     636     729     LYS             ASP             ASN             CnaA-like       0.0     0.1     0.984           8beg_A_590_636_729
            7woi                A       57      158     195     LYS             GLU             ASN             CnaB-like       0.001   0.1     0.984           7woi_A_57_158_195 
            5dz9                A       556     606     703     LYS             ASP             ASN             CnaA-like       0.0     0.1     0.984           4z1p_A_3_53_150   
            4z1p                A       3       53      150     LYS             ASP             ASN             CnaA-like       0.0     0.1     0.984           4z1p_A_3_53_150   
            6to1_af             A       13      334     420     LYS             ASP             ASN             CnaA-like       0.346   0.1     0.795           2x9z_A_193_241_318
                    
            Calculate solvent accessible area for a more accurate (and slow) prediction
            >>> i = Isopeptide("tests/data/test_structures", distance=1.5)
            >>> i.predict()
    """
    def __init__(self, pdb_dir: str, distance: float = 1.5, fixed_r_asa: float | None = None):
        """
        
            Raises
               ValueError if fixed_r_asa not between 0 and 1
        
        """
        self.pdb_dir: str = pdb_dir
        self.distance: float = distance
        self.jess_output: str | None = None
        self.isopeptide_bonds: List[BondElement] = []
        self.fixed_r_asa: float | None = fixed_r_asa
        self.bond_type: str | None = None
        if self.fixed_r_asa != None:
            if self.fixed_r_asa < 0 or self.fixed_r_asa > 1:
                raise ValueError(f"fixed_r_asa is not in 0-1 range. Found: {self.fixed_r_asa}")

    def predict(self):
        """

            Predict isopeptide bonds by 1. running jess template-based search,
            2. calculating asa, 3. predicting isopeptide bond probability with
            logistic regression.

        """

        self.jess_output = _run_jess(self.pdb_dir, self.distance)
        self._parse_jess_output()
        if len(self.isopeptide_bonds) == 0:
            warnings.warn("No isopeptide bond predictions detected. Try increasing the distance parameter.", UserWarning)
            return
        self._reduce_redundant()
        if self.fixed_r_asa != None:
            for bond in self.isopeptide_bonds:
                bond.r_asa = self.fixed_r_asa
        else:
            self._calc_rasa()
        # Make prediction with linear regression
        self._infer()
        # Sort based on probability
        self.isopeptide_bonds.sort(key=lambda x: (x.probability, x.protein_name), reverse=True)
        # Infer type of bond
        self._infer_type()

    def print_tabular(self):
        """
        
            Print isopeptide bonds in a tabular format
        
        """
        headers = [
        "protein_name", "chain", "r1_bond", "r_cat", "r2_bond",
        "r1_bond_name", "r_cat_name", "r2_bond_name", "bond_type",
        "rmsd", "r_asa", "probability", "template"
        ]
        column_widths = [max(len(header), max(len(str(getattr(bond, header))) for bond in self.isopeptide_bonds)) for header in headers]
        print("\t".join(headers))
        for bond in self.isopeptide_bonds:
            row = [
                bond.protein_name, bond.chain, bond.r1_bond, bond.r_cat, bond.r2_bond, 
                bond.r1_bond_name, bond.r_cat_name, bond.r2_bond_name, bond.bond_type,
                bond.rmsd, bond.r_asa, bond.probability, bond.template
            ]
            formatted_row = "\t".join(f"{str(item):<{column_widths[i]}}" for i, item in enumerate(row))
            print(formatted_row)

    def _parse_jess_output(self):
        """

            Parse jess output and stores as list of bondElement in self.isopeptide_bonds

            Raises:
                ValueError if jess output format or number of bond residues are different from expected

        """
        line_to_info = re.compile(r"REMARK (.*?\.pdb) ([0-9].[0-9]{3}) .+/(\w+_\w+_\d+_\d+_\d+).pdb Det= \d\.\d log\(E\)~ (-?\d+\.\d+)")
        for line in self.jess_output.split("\n"):
            if "REMARK" in line:
                match = line_to_info.findall(line)[0]
                if len(match) != 4:
                    raise ValueError(f"Jess output format different from expected: {line}")
                pdb_file = match[0]
                protein_name = pdb_file.split("/")[-1].replace(".pdb", "")
                rmsd = float(match[1])
                template = match[2].split("/")[-1].replace(".pdb", "")
                residues = []
                residue_names = []
                chain = None
                
            if line.startswith("ATOM"):
                chain = line[21]
                resi = int(line[22:26].strip())
                resi_name = line[17:20]
                if resi not in residues:
                    residues.append(resi)
                    residue_names.append(resi_name)
            
            if line.startswith("ENDMDL"):
                
                if len(residues) != 3:
                    raise ValueError(f"Number of residues found different from expected: 3!={len(residues)}")

                self.isopeptide_bonds.append(
                    BondElement(
                                pdb_file, protein_name, rmsd, template, chain, 
                                residues[0], residues[1], residues[2],
                                residue_names[0], residue_names[1], residue_names[2]
                    )
                )

    def _reduce_redundant(self):
        """

            Reduces redundant isopeptide bond predictions by keeping prediction with lower RMSD.

        """
        bonds = self.isopeptide_bonds
        filtered_bonds = []
        for protein_name in set([b.protein_name for b in bonds]):
            filtered_bonds.append(sorted([b for b in bonds if b.protein_name == protein_name], key = lambda b: (b.protein_name, b.r1_bond, b.r_cat, b.r2_bond, b.rmsd))[0])
        self.isopeptide_bonds = filtered_bonds

    def _calc_rasa(self):
        """

            Calculate r_asa for every isopeptide bond in every structure

        """
        bonds = self.isopeptide_bonds
        for pdb_file in set([b.pdb_file for b in bonds]):
            # Get asa and structure atom array pdb file
            structure_sasa, structure = _get_structure_asa(pdb_file)
            # Get asa of isopeptide residues
            for bond in [b for b in bonds if b.pdb_file == pdb_file]:
                isopep_residues = [bond.r1_bond, bond.r_cat, bond.r2_bond]
                isopep_residue_names = [bond.r1_bond_name, bond.r_cat_name, bond.r2_bond_name]
                tmp_r_asa = 0
                for res_id, res_name in zip(isopep_residues, isopep_residue_names):
                    res_indeces = [i for i, atom in enumerate(structure) if atom.res_id == res_id and atom.chain_id == bond.chain]
                    #res_name = structure[res_indeces[0]].res_name
                    if res_name not in MAX_ASA["rost_sander"].keys():
                        raise ValueError(f"{res_name} not in {MAX_ASA['rost_sander'].keys()}")
                    # Normalise ASA by residue surface area
                    r_asa = sum([structure_sasa[i] for i in res_indeces]) / MAX_ASA["rost_sander"][res_name]
                    tmp_r_asa += r_asa
                bond.r_asa = round(tmp_r_asa / 3, 3)
    
    def _infer(self):
        """

            Infer presence of isoeptide bond using logistic regression model

        """
        for bond in self.isopeptide_bonds:
            bond.probability = _predict(bond.rmsd, bond.r_asa)

    def _infer_type(self):
        """

            Infer isopeptide bond type (CnaA/B-like)

        """
        for bond in self.isopeptide_bonds:
            bond.bond_type = BOND_TYPE.get(bond.template, None)

if __name__ == "__main__":
    import doctest
    doctest.testmod(report=True, optionflags=doctest.NORMALIZE_WHITESPACE)