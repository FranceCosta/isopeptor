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
BondElement(struct_file=/nfs/research/agb/research/francesco/projects/20241024_isopeptor_v1/tests/data/test_structures/8beg.pdb, protein_name=8beg, rmsd=0.0, template=8beg_A_590_636_729, chain=A, r1_bond=590, r_cat=636, r2_bond=729, r1_bond_name=LYS, r_cat_name=ASP, r2_bond_name=ASN, bond_type=CnaA-like, r_asa=0.1, probability=0.984)
```
```
>>> i.print_tabular()
protein_name        probability     chain   r1_bond r_cat   r2_bond r1_bond_name    r_cat_name      r2_bond_name    bond_type       rmsd    r_asa   template
8beg                0.99            A       590     636     729     LYS             ASP             ASN             CnaA-like       0.0     0.1     8beg_A_590_636_729   
8beg                0.99            A       756     806     894     LYS             ASP             ASN             CnaA-like       0.0     0.1     8beg_A_756_806_894   
8beg                0.99            A       922     973     1049    LYS             ASP             ASN             CnaA-like       0.0     0.1     8beg_A_922_973_1049  
8beg                0.99            A       1076    1123    1211    LYS             ASP             ASN             CnaA-like       0.0     0.1     8beg_A_1076_1123_1211
5dz9                0.99            A       556     606     703     LYS             ASP             ASN             CnaA-like       0.0     0.1     4z1p_A_3_53_150      
5dz9                0.99            A       730     776     861     LYS             ASP             ASN             CnaA-like       0.0     0.1     4z1p_A_177_223_308   
4z1p                0.99            A       3       53      150     LYS             ASP             ASN             CnaA-like       0.0     0.1     4z1p_A_3_53_150      
4z1p                0.99            A       177     223     308     LYS             ASP             ASN             CnaA-like       0.0     0.1     4z1p_A_177_223_308   
7woi                0.907           B       57      158     195     LYS             GLU             ASN             CnaB-like       0.314   0.1     5j4m_A_47_139_172    
1amx                0.879           A       176     209     293     LYS             ASP             ASN             CnaA-like       0.353   0.1     2f68_X_176_209_293   
7woi                0.871           A       203     246     318     LYS             ASP             ASN             CnaA-like       0.363   0.1     4hss_A_187_224_299   
7woi                0.834           B       355     435     466     LYS             GLU             ASN             CnaB-like       0.403   0.1     8f70_A_299_386_437   
6to1_af             0.602           A       13      334     420     LYS             ASP             ASN             CnaA-like       0.565   0.1     5mkc_A_191_600_695                    

```

Calculate solvent accessible area for a more accurate (and slow) prediction:
```
>>> i = Isopeptide("tests/data/test_structures", distance=1.5)
>>> i.predict()
>>> i.print_tabular()
protein_name	probability	chain	r1_bond	r_cat	r2_bond	r1_bond_name	r_cat_name	r2_bond_name	bond_type	rmsd	r_asa	template
8beg        	0.992      	A    	590    	636  	729    	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.004	8beg_A_590_636_729   
8beg        	0.991      	A    	756    	806  	894    	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.028	8beg_A_756_806_894   
8beg        	0.991      	A    	922    	973  	1049   	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.015	8beg_A_922_973_1049  
8beg        	0.991      	A    	1076   	1123 	1211   	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.015	8beg_A_1076_1123_1211
5dz9        	0.991      	A    	556    	606  	703    	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.009	4z1p_A_3_53_150      
5dz9        	0.991      	A    	730    	776  	861    	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.019	4z1p_A_177_223_308   
4z1p        	0.991      	A    	3      	53   	150    	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.009	4z1p_A_3_53_150      
4z1p        	0.991      	A    	177    	223  	308    	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.019	4z1p_A_177_223_308   
7woi        	0.917      	B    	57     	158  	195    	LYS         	GLU       	ASN         	CnaB-like	0.314	0.027	5j4m_A_47_139_172    
1amx        	0.894      	A    	176    	209  	293    	LYS         	ASP       	ASN         	CnaA-like	0.353	0.016	2f68_X_176_209_293   
7woi        	0.887      	A    	203    	246  	318    	LYS         	ASP       	ASN         	CnaA-like	0.363	0.012	4hss_A_187_224_299   
7woi        	0.855      	B    	355    	435  	466    	LYS         	GLU       	ASN         	CnaB-like	0.403	0.009	8f70_A_299_386_437   
6to1_af     	0.642      	A    	13     	334  	420    	LYS         	ASP       	ASN         	CnaA-like	0.565	0.002	5mkc_A_191_600_695   
```

## Test

```
python -m unittest discover -s tests -p "test_isopeptide.py"
```

## TODO