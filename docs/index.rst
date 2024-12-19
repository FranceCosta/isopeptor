
Isopeptor Documentation
=======================
Python package for the detection of intamolecular isopeptide bonds in protein structures.  
The method described in "Isopeptor: a tool for detecting intramolecular isopeptide bonds in protein structures".

Workflow
--------

1. Jess-based template scan is used to detect potential isopeptide bond signatures in target protein structures. A total of 140 templates coming from high quality structures are used for this purpose. Input structures can be in either PDB or CIF format. 
2. If a site on the target protein structure matches with multiple templates, only the template-site match with the lowest RMSD is retained. 
3. Relative Accessible Solvent Area (rASA) is calculated. 
4. Both RMSD with the closest template and rASA are taken as input by the logistic regression model for isopeptide bond classification which outputs a probability. 
5. An optional final step is used to evaluate the geometry of isopeptide bonds. Two metrics are used to compare the parameters with the dataset of high quality structures: bond length Z-score and Kernel Density Estimate likelihood of dihedral angles.

Installation
------------

.. code-block:: console

    $ pip install isopeptor

Usage
-----

.. toctree::

   commandline
   api




