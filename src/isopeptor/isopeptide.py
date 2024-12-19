#! /usr/env/python
# -*- coding: utf-8 -*-

import re
from typing import List
import warnings
import numpy as np
import os
from pathlib import Path
from isopeptor.jess_wrapper import run_jess
from isopeptor.asa import get_structure_asa
from isopeptor.bond import BondElement 
from isopeptor.constants import MAX_ASA
from isopeptor.constants import BOND_TYPE
from isopeptor.logistic_regression import predict
from isopeptor.structure import get_structure
from isopeptor.bond_length import get_bond_stats
from isopeptor.dihedrals import get_dihedral_angles_stats

class Isopeptide:
    """
    
        Handles isopeptide bond prediction running pyjess via the jess_wrapper module
        and solvent accessible area via the asa module. Stores isopeptide bond predictions
        as a list of BondElement. Prediction is run for all structures from struct_dir. 
        Structures in .pdb format are directly analysed. If structures are in .cif format, they are first converted into .pdb.
        If multiple matches with an isopeptide bond signature are detected, only the one
        with the lowest RMSD is retained.

        Attributes:
            struct_dir (str): where pdb/cif files are located
            distance (float): that specifies permissivity of jess search
            jess_output (None | str): that stores jess raw output for debug purposes
            isopeptide_bonds (list): stores isopeptide bonds as BondElement elements
            fixed_r_asa (float | None): which ranges between 0 and 1 that fixes r_asa allowing to skip its calculation

        Example:
            
            >>> # Use a fixed solvent accessible area for a quick prediction:
            >>> from isopeptor.isopeptide import Isopeptide
            >>> i = Isopeptide("tests/data/test_structures", distance=1.5, fixed_r_asa=0.1)
            >>> i.predict()
            >>> i.isopeptide_bonds[0]
            BondElement(struct_file=tests/data/test_structures/8beg.pdb, protein_name=8beg, rmsd=0.0, template=8beg_A_590_636_729, chain=A, r1_bond=590, r_cat=636, r2_bond=729, r1_bond_name=LYS, r_cat_name=ASP, r2_bond_name=ASN, bond_type=CnaA-like, r_asa=0.1, probability=0.986)                     
            
            >>> # Calculate solvent accessible area for a more accurate (and slow) prediction:
            >>> i = Isopeptide("tests/data/test_structures", distance=1.5)
            >>> i.predict()
    """
    def __init__(self, struct_dir: str, distance: float = 1.5, fixed_r_asa: float | None = None):
        """
        
            Raises
               ValueError if fixed_r_asa not between 0 and 1
        
        """
        self.struct_dir: str = struct_dir
        self.distance: float = distance
        self.jess_hits: list | None = None
        self.isopeptide_bonds: List[BondElement] = []
        self.fixed_r_asa: float | None = fixed_r_asa
        if self.fixed_r_asa != None:
            if self.fixed_r_asa < 0 or self.fixed_r_asa > 1:
                raise ValueError(f"fixed_r_asa is not in 0-1 range. Found: {self.fixed_r_asa}")
        pdb_files = [str(p) for p in Path(self.struct_dir).glob("*.pdb")]
        cif_files = [str(p) for p in Path(self.struct_dir).glob("*.cif")]
        self.structure_files = pdb_files + cif_files
        self.geometry_calculated = False

    def predict(self):
        """

            Predict isopeptide bonds by 1. running jess template-based search,
            2. calculating asa, 3. predicting isopeptide bond probability with
            logistic regression.

        """

        self.jess_hits = run_jess(self.structure_files, self.distance)
        if len(self.jess_hits) == 0:
            warnings.warn("No isopeptide bond predictions detected. Try increasing the distance parameter.", UserWarning)
            return
        self._load_hits()
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

            Example:
                >>> from isopeptor.isopeptide import Isopeptide
                >>> i = Isopeptide("tests/data/test_structures", distance=1.5, fixed_r_asa=0.1)
                >>> i.predict()
                >>> i.print_tabular()
                protein_name        probability     chain   r1_bond r_cat   r2_bond r1_bond_name    r_cat_name      r2_bond_name    bond_type       rmsd    r_asa   template
                8beg                0.986           A       590     636     729     LYS             ASP             ASN             CnaA-like       0.0     0.1     8beg_A_590_636_729   
                8beg                0.986           A       756     806     894     LYS             ASP             ASN             CnaA-like       0.0     0.1     8beg_A_756_806_894   
                8beg                0.986           A       922     973     1049    LYS             ASP             ASN             CnaA-like       0.0     0.1     8beg_A_922_973_1049  
                8beg                0.986           A       1076    1123    1211    LYS             ASP             ASN             CnaA-like       0.0     0.1     8beg_A_1076_1123_1211
                5dz9                0.986           A       556     606     703     LYS             ASP             ASN             CnaA-like       0.0     0.1     4z1p_A_3_53_150      
                5dz9                0.986           A       730     776     861     LYS             ASP             ASN             CnaA-like       0.0     0.1     4z1p_A_177_223_308   
                4z1p                0.986           A       3       53      150     LYS             ASP             ASN             CnaA-like       0.0     0.1     4z1p_A_3_53_150      
                4z1p                0.986           A       177     223     308     LYS             ASP             ASN             CnaA-like       0.0     0.1     4z1p_A_177_223_308   
                7woi                0.935           B       57      158     195     LYS             GLU             ASN             CnaB-like       0.226   0.1     5j4m_A_47_139_172    
                7woi                0.89            A       203     246     318     LYS             ASP             ASN             CnaA-like       0.306   0.1     2woy_A_1259_1307_1393
                7woi                0.807           A       355     435     466     LYS             GLU             ASN             CnaB-like       0.398   0.1     6n0a_A_15_102_153    
                1amx                0.793           A       176     209     293     LYS             ASP             ASN             CnaA-like       0.411   0.1     2f68_X_176_209_293   
                6to1_af             0.495           A       13      334     420     LYS             ASP             ASN             CnaA-like       0.601   0.1     6to1_A_111_432_518   
                8beg                0.016           B       662     591     589     LYS             GLU             ASN             CnaB-like       1.176   0.1     3pg2_A_387_471_515
            
        """
        if not self.geometry_calculated:
            headers = [
            "protein_name", "probability", "chain", "r1_bond", "r_cat", "r2_bond",
            "r1_bond_name", "r_cat_name", "r2_bond_name", "bond_type",
            "rmsd", "r_asa", "template"
            ]
            # Print the header row using formatted widths
            column_widths = []
            for header in headers:
                # Find the maximum length among the bond attributes for the current header
                max_bond_length = 0
                for bond in self.isopeptide_bonds:
                    bond_value_length = len(str(getattr(bond, header)))
                    if bond_value_length > max_bond_length:
                        max_bond_length = bond_value_length
                # Compare it with the length of the header and append the maximum to column_widths
                column_widths.append(max(len(header), max_bond_length))
            formatted_header = "\t".join(f"{header:<{column_widths[i]}}" for i, header in enumerate(headers))
            print(formatted_header)
            
            if len(self.isopeptide_bonds) > 0:
                column_widths = [max(len(header), max(len(str(getattr(bond, header))) for bond in self.isopeptide_bonds)) for header in headers]
                for bond in self.isopeptide_bonds:
                    row = [
                        bond.protein_name, bond.probability, bond.chain, bond.r1_bond, bond.r_cat, bond.r2_bond, 
                        bond.r1_bond_name, bond.r_cat_name, bond.r2_bond_name, bond.bond_type,
                        bond.rmsd, bond.r_asa, bond.template
                    ]
                    formatted_row = "\t".join(f"{str(item):<{column_widths[i]}}" for i, item in enumerate(row))
                    print(formatted_row)
        else:
            headers = [
                "protein_name", "probability", "chain", "r1_bond", "r_cat", "r2_bond",
                "r1_bond_name", "r_cat_name", "r2_bond_name", "bond_type",
                "rmsd", "r_asa", "template", "bond_length", "bond_length_zscore", "bond_length_allowed",
                "pseudo_phi", "pseudo_psi", "pseudo_omega", "phi_psi_likelihood", "phi_psi_allowed", "omega_psi_likelihood",
                "omega_psi_allowed",  "omega_phi_likelihood", "omega_phi_allowed"
            ]
            # Print the header row using formatted widths
            column_widths = []
            for header in headers:
                # Find the maximum length among the bond attributes for the current header
                max_bond_length = 0
                for bond in self.isopeptide_bonds:
                    bond_value_length = len(str(getattr(bond, header)))
                    if bond_value_length > max_bond_length:
                        max_bond_length = bond_value_length
                # Compare it with the length of the header and append the maximum to column_widths
                column_widths.append(max(len(header), max_bond_length))
                
            formatted_header = "\t".join(f"{header:<{column_widths[i]}}" for i, header in enumerate(headers))
            print(formatted_header)
            if len(self.isopeptide_bonds) > 0:                
                # Print each data row using the same formatting
                for bond in self.isopeptide_bonds:
                    row = [
                        bond.protein_name, bond.probability, bond.chain, bond.r1_bond, bond.r_cat, bond.r2_bond, 
                        bond.r1_bond_name, bond.r_cat_name, bond.r2_bond_name, bond.bond_type,
                        bond.rmsd, bond.r_asa, bond.template, bond.bond_length, bond.bond_length_zscore, bond.bond_length_allowed,
                        bond.pseudo_phi, bond.pseudo_psi, bond.pseudo_omega, bond.phi_psi_likelihood, bond.phi_psi_allowed, bond.omega_psi_likelihood,
                        bond.omega_psi_allowed, bond.omega_phi_likelihood, bond.omega_phi_allowed
                    ]
                    formatted_row = "\t".join(f"{str(item):<{column_widths[i]}}" for i, item in enumerate(row))
                    print(formatted_row)
    
    def save_csv(self, output_table: str = "results.csv") -> None:
        """
        
            Save isopeptide bond results in .csv format.
            
            Args:
                output_table (str): path to output table
            
            Example:
                >>> from isopeptor.isopeptide import Isopeptide
                >>> i = Isopeptide("tests/data/test_structures", distance=1.5, fixed_r_asa=0.1)
                >>> i.predict()
                >>> i.save_csv()
        
        """
        with open(output_table, "wt") as fh:
            if not self.geometry_calculated:
                headers = [
                "protein_name", "probability", "chain", "r1_bond", "r_cat", "r2_bond",
                "r1_bond_name", "r_cat_name", "r2_bond_name", "bond_type",
                "rmsd", "r_asa", "template"
                ]
                fh.write(",".join(headers)+"\n")
                if len(self.isopeptide_bonds) > 0:
                    for bond in self.isopeptide_bonds:
                        row = [
                            bond.protein_name, str(bond.probability), str(bond.chain), str(bond.r1_bond), 
                            str(bond.r_cat), str(bond.r2_bond), 
                            bond.r1_bond_name, bond.r_cat_name, bond.r2_bond_name, bond.bond_type,
                            str(bond.rmsd), str(bond.r_asa), bond.template
                        ]
                        formatted_row = ",".join(row)
                        fh.write(formatted_row+"\n")
            else:
                headers = [
                    "protein_name", "probability", "chain", "r1_bond", "r_cat", "r2_bond",
                    "r1_bond_name", "r_cat_name", "r2_bond_name", "bond_type",
                    "rmsd", "r_asa", "template", "bond_length", "bond_length_zscore", "bond_length_allowed",
                    "pseudo_phi", "pseudo_psi", "pseudo_omega", "phi_psi_likelihood", "phi_psi_allowed", "omega_psi_likelihood",
                    "omega_psi_allowed",  "omega_phi_likelihood", "omega_phi_allowed"
                ]

                if len(self.isopeptide_bonds) > 0:
                    
                    formatted_header = ",".join(headers)
                    fh.write(formatted_header+"\n")
                    for bond in self.isopeptide_bonds:
                        row = [
                            bond.protein_name, str(bond.probability), str(bond.chain), 
                            str(bond.r1_bond), str(bond.r_cat), str(bond.r2_bond), 
                            bond.r1_bond_name, bond.r_cat_name, bond.r2_bond_name, bond.bond_type,
                            str(bond.rmsd), str(bond.r_asa), bond.template, str(bond.bond_length), 
                            str(bond.bond_length_zscore), bond.bond_length_allowed,
                            str(bond.pseudo_phi), str(bond.pseudo_psi), str(bond.pseudo_omega), 
                            str(bond.phi_psi_likelihood), str(bond.phi_psi_allowed), str(bond.omega_psi_likelihood),
                            bond.omega_psi_allowed, bond.omega_phi_likelihood, bond.omega_phi_allowed
                        ]
                        formatted_row = ",".join(row)
                        fh.write(formatted_row+"\n")

    def get_geometry(self):
        """
        
            Get geometry measures (bond length and dihedral angles and map to existing distribution of measures from 
            the database of PDB derived structures).

            Example:
                >>> from isopeptor.isopeptide import Isopeptide
                >>> i = Isopeptide("tests/data/test_structures", distance=1.5, fixed_r_asa=0.1)
                >>> i.predict()
                >>> i.get_geometry()
                >>> i.print_tabular()
                protein_name        probability     chain   r1_bond r_cat   r2_bond r1_bond_name    r_cat_name      r2_bond_name    bond_type       rmsd    r_asa      template                bond_length     bond_length_zscore      bond_length_allowed     pseudo_phi      pseudo_psi      pseudo_omega    phi_psi_likelihood phi_psi_allowed omega_psi_likelihood    omega_psi_allowed       omega_phi_likelihood    omega_phi_allowed
                8beg                0.986           A       590     636     729     LYS             ASP             ASN             CnaA-like       0.0     0.1        8beg_A_590_636_729      1.33            0.024                   True                    92.305          -123.286        -0.032          -8.706             True            -7.801                  True                    -8.797                  True             
                8beg                0.986           A       756     806     894     LYS             ASP             ASN             CnaA-like       0.0     0.1        8beg_A_756_806_894      1.328           -0.0                    True                    -114.423        -130.136        44.142          -10.503            False           -10.472                 True                    -10.503                 False            
                8beg                0.986           A       922     973     1049    LYS             ASP             ASN             CnaA-like       0.0     0.1        8beg_A_922_973_1049     1.333           0.061                   True                    70.944          -131.258        17.005          -9.737             True            -8.997                  True                    -9.882                  True             
                8beg                0.986           A       1076    1123    1211    LYS             ASP             ASN             CnaA-like       0.0     0.1        8beg_A_1076_1123_1211   1.341           0.159                   True                    102.317         -124.836        2.258           -8.286             True            -7.882                  True                    -8.233                  True             
                5dz9                0.986           A       556     606     703     LYS             ASP             ASN             CnaA-like       0.0     0.1        4z1p_A_3_53_150         1.291           -0.451                  True                    109.407         -112.845        -8.356          -7.72              True            -7.722                  True                    -7.884                  True             
                5dz9                0.986           A       730     776     861     LYS             ASP             ASN             CnaA-like       0.0     0.1        4z1p_A_177_223_308      1.332           0.049                   True                    97.334          -122.501        6.356           -8.391             True            -7.892                  True                    -8.633                  True             
                4z1p                0.986           A       3       53      150     LYS             ASP             ASN             CnaA-like       0.0     0.1        4z1p_A_3_53_150         1.291           -0.451                  True                    109.407         -112.845        -8.356          -7.72              True            -7.722                  True                    -7.884                  True             
                4z1p                0.986           A       177     223     308     LYS             ASP             ASN             CnaA-like       0.0     0.1        4z1p_A_177_223_308      1.332           0.049                   True                    97.334          -122.501        6.356           -8.391             True            -7.892                  True                    -8.633                  True             
                7woi                0.935           B       57      158     195     LYS             GLU             ASN             CnaB-like       0.226   0.1        5j4m_A_47_139_172       1.312           -0.195                  True                    -148.85         71.47           179.899         -9.376             True            -9.126                  True                    -9.675                  True             
                7woi                0.89            A       203     246     318     LYS             ASP             ASN             CnaA-like       0.306   0.1        2woy_A_1259_1307_1393   1.336           0.098                   True                    136.992         -158.356        4.495           -14.872            False           -12.529                 False                   -9.337                  True             
                7woi                0.807           A       355     435     466     LYS             GLU             ASN             CnaB-like       0.398   0.1        6n0a_A_15_102_153       1.328           -0.0                    True                    106.44          118.566         168.851         -11.023            False           -10.77                  False                   -10.955                 False            
                1amx                0.793           A       176     209     293     LYS             ASP             ASN             CnaA-like       0.411   0.1        2f68_X_176_209_293      3.345           24.598                  False                   96.887          -147.455        66.383          -10.918            False           -14.477                 False                   -25.835                 False            
                6to1_af             0.495           A       13      334     420     LYS             ASP             ASN             CnaA-like       0.601   0.1        6to1_A_111_432_518      2.928           19.512                  False                   83.132          -155.408        131.186         -13.496            False           -51.581                 False                   -15.775                 False            
                8beg                0.016           B       662     591     589     LYS             GLU             ASN             CnaB-like       1.176   0.1        3pg2_A_387_471_515      4.806           42.415                  False                   -175.107        131.767         146.384         -11.214            False           -13.479                 False                   -11.449                 False              

        """
        self.geometry_calculated = True
        bonds = self.isopeptide_bonds
        for struct_file in set([b.struct_file for b in bonds]):
            structure = get_structure(struct_file)
            for bond in [b for b in bonds if b.struct_file == struct_file]:
                # Get bond lenght and stats
                bond_length, bond_length_zscore, bond_length_allowed = get_bond_stats(structure, bond.chain, bond.r1_bond, 
                                bond.r2_bond, bond.r2_bond_name)
                
                # Get dihedrals and stats
                angle_stats = get_dihedral_angles_stats(structure, bond.chain, bond.r1_bond, 
                                bond.r2_bond, bond.r2_bond_name, bond.bond_type)
                pseudo_omega, pseudo_phi, pseudo_psi, phi_psi_likelihood, omega_psi_likelihood, omega_phi_likelihood, phi_psi_allowed, omega_psi_allowed, omega_phi_allowed = angle_stats

                # Load values into the bond instance
                bond.bond_length = bond_length
                bond.bond_length_zscore = bond_length_zscore
                bond.bond_length_allowed = bond_length_allowed
                bond.pseudo_omega = pseudo_omega
                bond.pseudo_phi = pseudo_phi
                bond.pseudo_psi = pseudo_psi
                bond.phi_psi_likelihood = phi_psi_likelihood
                bond.omega_psi_likelihood = omega_psi_likelihood
                bond.omega_phi_likelihood = omega_phi_likelihood
                bond.phi_psi_allowed = phi_psi_allowed
                bond.omega_psi_allowed = omega_psi_allowed
                bond.omega_phi_allowed = omega_phi_allowed

    def _load_hits(self):
        """

            Load pyjess._jess.Hit as list of bondElement in self.isopeptide_bonds.

            Raises:
                ValueError if number of residues found is not expected.

        """
        for hit in self.jess_hits:
            template, rmsd, struct_file, atoms = hit.template.id, hit.rmsd, hit.molecule.id, hit.atoms()
            protein_name = os.path.basename(struct_file).replace(".pdb", "").replace(".cif", "")
            residues, residue_names = [], []
            for atom in atoms:
                if atom.residue_number not in residues:
                    residues.append(atom.residue_number)
                    residue_names.append(atom.residue_name)
            chain = atom.chain_id

            if len(residues) == 3:
                self.isopeptide_bonds.append(
                    BondElement(
                                struct_file, protein_name, round(rmsd, 3), template, chain, 
                                residues[0], residues[1], residues[2],
                                residue_names[0], residue_names[1], residue_names[2]
                    )
                )
            else:
                raise ValueError(f"Found {len(residues)} residues.")

    def _reduce_redundant(self):
        """

            Reduces redundant isopeptide bond predictions by keeping prediction with lower RMSD.

        """
        bonds = self.isopeptide_bonds
        grouped_bonds = {}
        
        for bond in bonds:
            key = (bond.protein_name, bond.r1_bond, bond.r_cat, bond.r2_bond)
            if key not in grouped_bonds or bond.rmsd < grouped_bonds[key].rmsd:
                grouped_bonds[key] = bond

        # Keep only the best bond for each group
        self.isopeptide_bonds = list(grouped_bonds.values())

    def _calc_rasa(self):
        """

            Calculate r_asa for every isopeptide bond in every structure.

        """
        bonds = self.isopeptide_bonds
        for struct_file in set([b.struct_file for b in bonds]):
            # Get asa and structure atom array pdb file
            structure_sasa, structure = get_structure_asa(struct_file)
            # Get asa of isopeptide residues
            for bond in [b for b in bonds if b.struct_file == struct_file]:
                isopep_residues = [bond.r1_bond, bond.r_cat, bond.r2_bond]
                isopep_residue_names = [bond.r1_bond_name, bond.r_cat_name, bond.r2_bond_name]
                tmp_r_asa = 0
                for res_id, res_name in zip(isopep_residues, isopep_residue_names):
                    # Handle the case of two res bb
                    if res_name == None:
                        continue
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

            Infer presence of isoeptide bond using logistic regression model.

        """
        for bond in self.isopeptide_bonds:
            bond.probability = predict(bond.rmsd, bond.r_asa)

    def _infer_type(self):
        """

            Infer isopeptide bond type (CnaA/B-like).

        """
        for bond in self.isopeptide_bonds:
            bond.bond_type = BOND_TYPE.get(bond.template, None)

if __name__ == "__main__":
    import doctest
    doctest.testmod(report=True, optionflags=doctest.NORMALIZE_WHITESPACE)