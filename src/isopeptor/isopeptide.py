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
        self._infer()

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
                chain = None
                
            if line.startswith("ATOM"):
                chain = line[21]
                resi = int(line[22:26].strip())
                residues.append(resi)
            
            if line.startswith("ENDMDL"):
                residues = list(set(residues))
                
                if len(residues) != 3:
                    raise ValueError(f"Number of residues found different from expected: 3!={len(residues)}")

                self.isopeptide_bonds.append(
                    BondElement(
                                pdb_file, protein_name, rmsd, template, chain, residues[0], residues[1], residues[2], 
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

            Calcolate r_asa for every isopeptide bond in every structure

        """
        bonds = self.isopeptide_bonds
        for pdb_file in set([b.pdb_file for b in bonds]):
            # Get asa and structure atom array pdb file
            structure_sasa, structure = _get_structure_asa(pdb_file)
            # Get asa of isopeptide residues
            for bond in [b for b in bonds if b.pdb_file == pdb_file]:
                isopep_residues = [bond.r1_bond, bond.r_cat, bond.r2_bond]
                tmp_r_asa = 0
                for res_id in isopep_residues:
                    res_indeces = [i for i, atom in enumerate(structure) if atom.res_id == res_id and atom.chain_id == bond.chain]
                    res_name = structure[res_indeces[0]].res_name
                    if res_name not in MAX_ASA["rost_sander"].keys():
                        raise ValueError(f"{res_name} not in {MAX_ASA['rost_sander'].keys()}")
                    # Normalise ASA by residue surface area
                    r_asa = sum([structure_sasa[i] for i in res_indeces]) / MAX_ASA["rost_sander"][res_name]
                    tmp_r_asa += r_asa
                bond.r_asa = "%.3f" % round(tmp_r_asa / 3, 3)
    
    def _infer(self):
        """

            Infer presence of isoeptide bond using logistic regression model

        """
        for bond in self.isopeptide_bonds:
            bond.probability = _predict(bond.rmsd, bond.r_asa)