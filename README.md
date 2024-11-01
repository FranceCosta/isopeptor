# ISOPEPtide bond detecTOR

Package for isopeptide bond prediction.

## Installation

`pip install isopeptor`

## Usage

On the command line:
`isopeptor tests/data/test_structures/`

Output:
```
protein_name    chain   r1_bond r_cat   r2_bond r1_bond_name    r_cat_name      r2_bond_name    bond_type       rmsd    r_asa   probability     template
6to1_af         A       13      334     420     LYS             ASP             ASN             CnaA-like       0.346   0.002   0.835           2x9z_A_193_241_318
7woi            A       57      158     195     LYS             GLU             ASN             CnaB-like       0.001   0.019   0.987           7woi_A_57_158_195 
5dz9            A       556     606     703     LYS             ASP             ASN             CnaA-like       0.0     0.009   0.987           4z1p_A_3_53_150   
8beg            A       590     636     729     LYS             ASP             ASN             CnaA-like       0.0     0.004   0.987           8beg_A_590_636_729
4z1p            A       3       53      150     LYS             ASP             ASN             CnaA-like       0.0     0.009   0.987           4z1p_A_3_53_150
```