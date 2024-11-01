#! /usr/env/python
# -*- coding: utf-8 -*-

"""

    Run isopeptide bond prediction from command line

"""

from isopeptor import Isopeptide 
import argparse

parser = argparse.ArgumentParser(
    prog = "isopeptor",
    description = "Run isopeptide bond prediction from command line. Usage: isopeptor path/to/*.pdb > isopeptide_bonds.csv"
)

parser.add_argument(
    "--path_to_pdb_files", 
    required=True, 
    help="Path to directory containing .pdb files.", 
    type=str
)

parser.add_argument(
    "--fixed_r_asa", 
    required=False, 
    help="Fixes the relative solvent accessible area using a value between 0 and 1 to speed up the prediction.", 
    type=float,
    default=None
)

def main():
    args = parser.parse_args()
    

if __name__ == "__main__":
    main()