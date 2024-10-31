#! /usr/env/python
# -*- coding: utf-8 -*-

from isopeptor.jess_wrapper import run_jess
from isopeptor.asa import calculate_asa
import re

class Isopeptide:
    def __init__(self, pdb_dir : str, distance : float = 1.5):
        self.pdb_dir = pdb_dir
        self.distance = distance
        self.jess_output = None
        self.isopeptide_bonds = []

    def run_jess_and_parse(self):
        self.jess_output = run_jess(self.pdb_dir, self.distance)
        self._parse_jess_output()

    def _parse_jess_output(self):
        line_to_info = re.compile(r"([A-Z0-9]{,15}).pdb ([0-9].[0-9]{3}) .+/(\w+_\w+_\d+_\d+_\d+).pdb Det= \d\.\d log\(E\)~ (-?\d+\.\d+)")
        for line in self.jess_output.split("\n"):
            if "REMARK" in line:
                line = line.split("REMARK")[1]
                match_ = line_to_info.findall(line)[0]
                pdb_file = match_[0]
                protein_name = pdb_file.split("/")[-1].replace(".pdb", "")
                rmsd = float(match_[1])
                template = match_[2].split("/")[-1].replace(".pdb", "")
                #loge = float(match_[3].strip())
                residues = []
                chain = None
                
            if line.startswith("ATOM"):
                chain = line[21]
                resi = int(line[22:26].strip())
                residues.append(resi)
            
            if line.startswith("ENDMDL"):
                self.isopeptide_bonds.append(BondElement(
                    pdb_file, 
                    protein_name,
                    rmsd, 
                    template,
                    chain,
                    residues[0],
                    residues[1],
                    residues[2], 
                    ))

class BondElement:
    def __init__(self, pdb_file : str, protein_name : str, rmsd : float, template : str, chain : str, r1_bond : int, r_cat : int, r2_bond : int):
        self.pdb_file = pdb_file
        self.protein_name = protein_name
        self.rmsd = rmsd
        self.template = template
        self.chain = chain
        self.r1_bond = r1_bond
        self.r_cat = r_cat
        self.r2_bond = r2_bond

    def __repr__(self):
        s = f"BondElement(pdb_file={self.pdb_file}, protein_name={self.protein_name}, rmsd={self.rmsd}, template={self.template},"+\
            f"template={self.template},)"
        return s
