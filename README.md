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
protein_name        probability     chain   r1_bond r_cat   r2_bond r1_bond_name    r_cat_name      r2_bond_name    bond_type       rmsd    r_asa   template             
8beg                0.991           A       590     636     729     LYS             ASP             ASN             CnaA-like       0.0     0.1     8beg_A_590_636_729   
8beg                0.991           A       756     806     894     LYS             ASP             ASN             CnaA-like       0.0     0.1     8beg_A_756_806_894   
8beg                0.991           A       922     973     1049    LYS             ASP             ASN             CnaA-like       0.0     0.1     8beg_A_922_973_1049  
8beg                0.991           A       1076    1123    1211    LYS             ASP             ASN             CnaA-like       0.0     0.1     8beg_A_1076_1123_1211
5dz9                0.991           A       556     606     703     LYS             ASP             ASN             CnaA-like       0.0     0.1     4z1p_A_3_53_150      
5dz9                0.991           A       730     776     861     LYS             ASP             ASN             CnaA-like       0.0     0.1     4z1p_A_177_223_308   
4z1p                0.991           A       3       53      150     LYS             ASP             ASN             CnaA-like       0.0     0.1     4z1p_A_3_53_150      
4z1p                0.991           A       177     223     308     LYS             ASP             ASN             CnaA-like       0.0     0.1     4z1p_A_177_223_308   
7woi                0.909           B       57      158     195     LYS             GLU             ASN             CnaB-like       0.314   0.1     5j4m_A_47_139_172    
1amx                0.882           A       176     209     293     LYS             ASP             ASN             CnaA-like       0.353   0.1     2f68_X_176_209_293   
7woi                0.875           A       203     246     318     LYS             ASP             ASN             CnaA-like       0.363   0.1     4hss_A_187_224_299   
7woi                0.838           B       355     435     466     LYS             GLU             ASN             CnaB-like       0.403   0.1     8f70_A_299_386_437   
6to1_af             0.607           A       13      334     420     LYS             ASP             ASN             CnaA-like       0.565   0.1     5mkc_A_191_600_695   

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

## Reference
If you use isopeptor please cite:

```bibtex
@article{10.1093/bioadv/vbaf049,
    author = {Costa, Francesco and Barringer, Rob and Riziotis, Ioannis and Andreeva, Antonina and Bateman, Alex},
    title = {Isopeptor: a tool for detecting intramolecular isopeptide bonds in protein structures},
    journal = {Bioinformatics Advances},
    volume = {5},
    number = {1},
    pages = {vbaf049},
    year = {2025},
    month = {03},
    issn = {2635-0041},
    doi = {10.1093/bioadv/vbaf049},
    url = {https://doi.org/10.1093/bioadv/vbaf049},
    eprint = {https://academic.oup.com/bioinformaticsadvances/article-pdf/5/1/vbaf049/62375500/vbaf049.pdf},
}
```
