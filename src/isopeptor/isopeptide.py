#! /usr/env/python
# -*- coding: utf-8 -*-

from asa import calculate_asa
from jess_wrapper import run_jess

class Isopeptide:
    def __init__(self, pdb_dir : str, distance : float = 1.5):
        self.pdb_dir = pdb_dir
        self.distance = distance
        self.jess_output = None
        self.bonds = None

        def run_jess_and_parse(self):
            self.jess_output = run_jess(self.pdb_dir, self.distance)
            #self.parse_jess_output()

        def parse_jess_output(self):
            # Assuming jess_output lines have a format like "pdb_file bond_type bond_strength"
            for line in self.jess_output.splitlines():
                parts = line.split()
                if len(parts) == 3:
                    pdb_file, bond_type, bond_strength = parts
                    bond_strength = float(bond_strength)  # Convert strength to float
                    self.bonds.append(BondElement(pdb_file, bond_type, bond_strength))
                else:
                    raise ValueError(f"Unexpected line format in jess_output: '{line}'")

class BondElement:
    def __init__(self, pdb_file: str):
        self.pdb_file = pdb_file

    def __repr__(self):
        return f"BondElement(pdb_file={self.pdb_file}, bond_type={self.bond_type}, bond_strength={self.bond_strength})"
