# ISOPEPtide bond detecTOR

![Workflow](./workflow.png)

Isopeptide bond prediction from .pdb/.cif structure. Method described in "Isopeptor: a tool for detecting intramolecular isopeptide bonds in protein structures".

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

Specifying --fixed_r_asa with a relative solvent accessibile area between 0 and 1 allows skipping its time-consuming calculation. The downside is that the prediction is less reliable. It is useful for very high throughput screenings where high precision is not required.


To run a geometric evaluation use the `--eval_geometry` flag. This will report bond length, dihedral angles (pseudo-psi,phi and psi) and statistical measures. Two metrics are used to compare the parameters with the dataset of high quality structures: bond length Z-score and Kernel Density Estimate likelihood of dihedral angles. Bond length outliers are considered for values exceeding 4 standard deviations from the mean while angle outliers are calculated for each pair of dihedral angles for values of KDE likelihood exceeding the 95th percentile of the high wuality dataset distributions. Outlier parameters are marked under \<PARAMETER\>_allowed column.

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

To include geometric measures use:

```
>>> i = Isopeptide("tests/data/test_structures", distance=1.5, fixed_r_asa=0.1)
>>> i.predict()
>>> i.get_geometry()
>>> i.print_tabular()
protein_name	probability	chain	r1_bond	r_cat	r2_bond	r1_bond_name	r_cat_name	r2_bond_name	bond_type	rmsd 	r_asa	template             	bond_length	bond_length_zscore	bond_length_allowed	phi     	psi     	omega  	phi_psi_likelihood	phi_psi_allowed	omega_psi_likelihood	omega_psi_allowed	omega_phi_likelihood	omega_phi_allowed
8beg        	0.986      	A    	590    	636  	729    	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.1  	8beg_A_590_636_729   	1.33       	0.024             	True               	92.305  	-123.286	-0.032 	-8.706            	True           	-7.801              	True             	-8.797              	True             
8beg        	0.986      	A    	756    	806  	894    	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.1  	8beg_A_756_806_894   	1.328      	0.0               	True               	-114.423	-130.136	44.142 	-10.503           	False          	-10.472             	True             	-10.503             	False            
8beg        	0.986      	A    	922    	973  	1049   	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.1  	8beg_A_922_973_1049  	1.333      	0.061             	True               	70.944  	-131.258	17.005 	-9.737            	True           	-8.997              	True             	-9.882              	True             
8beg        	0.986      	A    	1076   	1123 	1211   	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.1  	8beg_A_1076_1123_1211	1.341      	0.159             	True               	102.317 	-124.836	2.258  	-8.286            	True           	-7.882              	True             	-8.233              	True             
5dz9        	0.986      	A    	556    	606  	703    	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.1  	4z1p_A_3_53_150      	1.291      	-0.451            	True               	109.407 	-112.845	-8.356 	-7.72             	True           	-7.722              	True             	-7.884              	True             
5dz9        	0.986      	A    	730    	776  	861    	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.1  	4z1p_A_177_223_308   	1.332      	0.049             	True               	97.334  	-122.501	6.356  	-8.391            	True           	-7.892              	True             	-8.633              	True             
4z1p        	0.986      	A    	3      	53   	150    	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.1  	4z1p_A_3_53_150      	1.291      	-0.451            	True               	109.407 	-112.845	-8.356 	-7.72             	True           	-7.722              	True             	-7.884              	True             
4z1p        	0.986      	A    	177    	223  	308    	LYS         	ASP       	ASN         	CnaA-like	0.0  	0.1  	4z1p_A_177_223_308   	1.332      	0.049             	True               	97.334  	-122.501	6.356  	-8.391            	True           	-7.892              	True             	-8.633              	True             
7woi        	0.935      	B    	57     	158  	195    	LYS         	GLU       	ASN         	CnaB-like	0.226	0.1  	5j4m_A_47_139_172    	1.312      	-0.195            	True               	-148.85 	71.47   	179.899	-9.376            	True           	-9.126              	True             	-9.675              	True             
7woi        	0.89       	A    	203    	246  	318    	LYS         	ASP       	ASN         	CnaA-like	0.306	0.1  	2woy_A_1259_1307_1393	1.336      	0.098             	True               	136.992 	-158.356	4.495  	-14.872           	False          	-12.529             	False            	-9.337              	True             
7woi        	0.807      	A    	355    	435  	466    	LYS         	GLU       	ASN         	CnaB-like	0.398	0.1  	6n0a_A_15_102_153    	1.328      	0.0               	True               	106.44  	118.566 	168.851	-11.023           	False          	-10.77              	False            	-10.955             	False            
1amx        	0.793      	A    	176    	209  	293    	LYS         	ASP       	ASN         	CnaA-like	0.411	0.1  	2f68_X_176_209_293   	3.345      	24.598            	False              	96.887  	-147.455	66.383 	-10.918           	False          	-14.477             	False            	-25.835             	False            
6to1_af     	0.495      	A    	13     	334  	420    	LYS         	ASP       	ASN         	CnaA-like	0.601	0.1  	6to1_A_111_432_518   	2.928      	19.512            	False              	83.132  	-155.408	131.186	-13.496           	False          	-51.581             	False            	-15.775             	False            
8beg        	0.016      	B    	662    	591  	589    	LYS         	GLU       	ASN         	CnaB-like	1.176	0.1  	3pg2_A_387_471_515   	4.806      	42.415            	False              	-175.107	131.767 	146.384	-11.214           	False          	-13.479             	False            	-11.449             	False            

```

## Test

```
python -m unittest discover -s tests -p "test_isopeptide.py"
```

## Notes

- .CIF file parsing is handled via .CIF to .PDB conversion at the moment at may not be suitable for very large protein structures.