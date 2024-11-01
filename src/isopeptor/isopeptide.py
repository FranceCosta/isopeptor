#! /usr/env/python
# -*- coding: utf-8 -*-

from isopeptor.jess_wrapper import _run_jess
from isopeptor.asa import calculate_asa
import re
from typing import List
from isopeptor.bond import BondElement 

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

    """
    def __init__(self, pdb_dir: str, distance: float = 1.5):
        self.pdb_dir: str = pdb_dir
        self.distance: float = distance
        self.jess_output: str | None = None
        self.isopeptide_bonds: List[BondElement] = []

    def predict(self):
        """

            Predict isopeptide bonds by 1. running jess template-based search,
            2. calculating asa, 3. predicting isopeptide bond probability with
            logistic regression.

        """

        self.jess_output = _run_jess(self.pdb_dir, self.distance)
        self._parse_jess_output()
        
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

            Raises:
                ValueError if no isopeptide bond predictions are detected in self.isopeptide_bonds

        """
        if len(self.isopeptide_bonds) == 0:
            raise ValueError("No isopeptide bond predictions detected.")
        