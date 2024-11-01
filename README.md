# ISOPEPtide bond detecTOR

Package for isopeptide bond prediction.

## Installation

`pip install isopeptor`

## Usage

From the command line:
```
isopeptor tests/data/test_structures/
```

Output:
```
protein_name    chain   r1_bond r_cat   r2_bond r1_bond_name    r_cat_name      r2_bond_name    bond_type       rmsd    r_asa   probability     template
6to1_af         A       13      334     420     LYS             ASP             ASN             CnaA-like       0.346   0.002   0.835           2x9z_A_193_241_318
7woi            A       57      158     195     LYS             GLU             ASN             CnaB-like       0.001   0.019   0.987           7woi_A_57_158_195 
5dz9            A       556     606     703     LYS             ASP             ASN             CnaA-like       0.0     0.009   0.987           4z1p_A_3_53_150   
8beg            A       590     636     729     LYS             ASP             ASN             CnaA-like       0.0     0.004   0.987           8beg_A_590_636_729
4z1p            A       3       53      150     LYS             ASP             ASN             CnaA-like       0.0     0.009   0.987           4z1p_A_3_53_150
```

To redirect the output to a `.tsv` file use:

```
isopeptor tests/data/test_structures/ > output.tsv
```

## Advanced usage

### Full command line options

```
usage: isopeptor [-h] [--distance DISTANCE] [--fixed_r_asa FIXED_R_ASA] path_to_pdb_files

Run isopeptide bond prediction from command line. Usage: isopeptor path/to pdb files/ > isopeptide_bonds.csv

positional arguments:
  path_to_pdb_files     Path to directory containing .pdb files.

options:
  -h, --help            show this help message and exit
  --distance DISTANCE   Specifies permissivity of jess search. The higher, the more permissive.
  --fixed_r_asa FIXED_R_ASA
                        Fixes the relative solvent accessible area using a value between 0 and 1 to speed up the prediction.
```

### Python API

```


```

## Test

```
python -m unittest discover -s tests -p "test_isopeptide.py"
```
