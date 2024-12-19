# ISOPEPtide bond detecTOR

Python package for the detection of intamolecular isopeptide bonds in protein structures. 
The method is described in "Isopeptor: a tool for detecting intramolecular isopeptide bonds in protein structures".

Isopeptor can be accessed via [this google colab](https://colab.research.google.com/github/FranceCosta/Isopeptor_development/blob/main/notebooks/Isopeptide_finder.ipynb) or installed and run locally.

Read the [documentation](https://isopeptor.readthedocs.io/en/latest/index.html) for complete information on its usage and examples.

## Installation

```
pip install isopeptor
```

## Usage

From the command line:
```
isopeptor tests/data/test_structures/
```

Output:
```
protein_name	probability	chain	r1_bond	r_cat	r2_bond	r1_bond_name	r_cat_name	r2_bond_name	bond_type	rmsd	r_asa	template
8beg        	0.987      	A    	590    	636  	729    	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.004	8beg_A_590_636_729   
8beg        	0.987      	A    	756    	806  	894    	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.028	8beg_A_756_806_894   
8beg        	0.987      	A    	922    	973  	1049   	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.015	8beg_A_922_973_1049  
8beg        	0.987      	A    	1076   	1123 	1211   	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.015	8beg_A_1076_1123_1211
7woi        	0.987      	A    	57     	158  	195    	LYS         	GLU       	ASN         	CnaB-like	0.001	0.019	7woi_A_57_158_195    
7woi        	0.987      	A    	203    	246  	318    	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.012	7woi_A_203_246_318   
5dz9        	0.987      	A    	556    	606  	703    	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.009	4z1p_A_3_53_150      
5dz9        	0.987      	A    	730    	776  	861    	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.019	4z1p_A_177_223_308   
4z1p        	0.987      	A    	3      	53   	150    	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.009	4z1p_A_3_53_150      
4z1p        	0.987      	A    	177    	223  	308    	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.019	4z1p_A_177_223_308   
6to1_af     	0.505      	A    	13     	334  	420    	LYS         	ASP       	ASN         	CnaA-like	0.548	0.002	4uzg_A_187_225_330
```

To redirect the output to a `.tsv` file use:

```
isopeptor tests/data/test_structures/ > output.tsv
```

### Full command line options:

```
usage: isopeptor [-h] [--distance DISTANCE] [--fixed_r_asa FIXED_R_ASA] [--eval_geometry] path_to_structure_files

Run isopeptide bond prediction from command line. Usage: isopeptor path/to pdb files/ > isopeptide_bonds.csv

positional arguments:
  path_to_structure_files
                        Path to directory containing .pdb/.cif files.

options:
  -h, --help            show this help message and exit
  --distance DISTANCE   Specifies permissivity of jess search. The higher, the more permissive.
  --fixed_r_asa FIXED_R_ASA
                        Fixes the relative solvent accessible area using a value between 0 and 1 to speed up the prediction.
  --eval_geometry       Run geometric evaluation of isopeptide bonds.
```

## Test

```
python -m unittest discover -s tests -p "test_isopeptide.py"
```
