Command line usage
==================

.. automodule:: isopeptor.cli
    :members:

Usage
-----
.. Example::
    
    >>>isopeptor -h
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
    