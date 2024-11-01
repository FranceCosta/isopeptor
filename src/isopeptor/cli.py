#! /usr/env/python
# -*- coding: utf-8 -*-

"""

    Run isopeptide bond prediction from command line

"""

from isopeptor.isopeptide import Isopeptide 
import argparse

parser = argparse.ArgumentParser(
    prog = "isopeptor",
    description = "Run isopeptide bond prediction from command line. Usage: isopeptor path/to pdb files/ > isopeptide_bonds.csv"
)

parser.add_argument(
    "path_to_pdb_files", 
    help="Path to directory containing .pdb files.", 
    type=str
)

parser.add_argument(
    "--distance", 
    required=False, 
    help="Specifies permissivity of jess search. The higher, the more permissive.", 
    type=float,
    default=1.5
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
    i = Isopeptide(args.path_to_pdb_files, args.distance, args.fixed_r_asa)
    i.predict()
    i.print_tabular()

if __name__ == "__main__":
    main()