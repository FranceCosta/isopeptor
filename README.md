# ISOPEPtide bond detecTOR

Isopeptide bond prediction from .pdb structure. Method described in "A global survey of intramolecular isopeptide bonds".

The method can be accessed via [this google colab]() or installed and run locally.

## Installation

```
git clone https://github.com/FranceCosta/isopeptor
cd isopeptor 
pip install .
```

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

### Full command line options:

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

Specifying --fixed_r_asa with a relative solvent accessibile area between 0 and 1 allows skipping its time-consuming calculation. The downside is that the prediction is less reliable. It is useful for very high throughput screenings where high precision is not required.

### Python API

Quick prediction:
```
>>> from isopeptor.isopeptide import Isopeptide
>>> i = Isopeptide("tests/data/test_structures", distance=1.5, fixed_r_asa=0.1)
>>> i.predict()
>>> i.isopeptide_bonds[0]
BondElement(pdb_file=/nfs/research/agb/research/francesco/projects/20241024_isopeptor_v1/tests/data/test_structures/8beg.pdb, protein_name=8beg, rmsd=0.0, template=8beg_A_590_636_729, chain=A, r1_bond=590, r_cat=636, r2_bond=729, r1_bond_name=LYS, r_cat_name=ASP, r2_bond_name=ASN, bond_type=CnaA-like, r_asa=0.1, probability=0.984)
```
```
>>> i.print_tabular()
protein_name        chain   r1_bond r_cat   r2_bond r1_bond_name    r_cat_name      r2_bond_name    bond_type       rmsd    r_asa   probability     template
8beg                A       590     636     729     LYS             ASP             ASN             CnaA-like       0.0     0.1     0.984           8beg_A_590_636_729
7woi                A       57      158     195     LYS             GLU             ASN             CnaB-like       0.001   0.1     0.984           7woi_A_57_158_195 
5dz9                A       556     606     703     LYS             ASP             ASN             CnaA-like       0.0     0.1     0.984           4z1p_A_3_53_150   
4z1p                A       3       53      150     LYS             ASP             ASN             CnaA-like       0.0     0.1     0.984           4z1p_A_3_53_150   
6to1_af             A       13      334     420     LYS             ASP             ASN             CnaA-like       0.346   0.1     0.795           2x9z_A_193_241_318
```

Calculate solvent accessible area for a more accurate (and slow) prediction:
```
>>> i = Isopeptide("tests/data/test_structures", distance=1.5)
>>> i.predict()
protein_name	probability	chain	r1_bond	r_cat	r2_bond	r1_bond_name	r_cat_name	r2_bond_name	bond_type	rmsd	r_asa	template
8beg        	0.984      	A    	590    	636  	729    	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.1  	8beg_A_590_636_729
7woi        	0.984      	A    	57     	158  	195    	LYS         	GLU       	ASN         	CnaB-like	0.001	0.1  	7woi_A_57_158_195 
5dz9        	0.984      	A    	556    	606  	703    	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.1  	4z1p_A_3_53_150   
4z1p        	0.984      	A    	3      	53   	150    	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.1  	4z1p_A_3_53_150   
6to1_af     	0.795      	A    	13     	334  	420    	LYS         	ASP       	ASN         	CnaA-like	0.346	0.1  	2x9z_A_193_241_318
```

## Test

```
python -m unittest discover -s tests -p "test_isopeptide.py"
```


## TODO
- modify _reduce_redundant() to select unique by protein match and not by protein!
- get asa by chain