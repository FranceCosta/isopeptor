Isopeptor Documentation
=======================
Python package for the detection of intamolecular isopeptide bonds in protein structures.  
The method is described in **Isopeptor: a tool for detecting intramolecular isopeptide bonds in protein structures**.

Intramolecular Isopeptide bonds
-------------------------------
Intramolecular isopeptide bonds typically form between lysine and asparagine/aspartate residues, 
typically catalyzed by a nearby aspartate or glutamate. 
When found in β-sandwich folds, they occur in CnaA domains (linking β-strands of opposing β-sheets) and 
CnaB domains (linking adjacent β-strands of the same β-sheet). 
They stabilize protein structures against various stresses.

.. note::
 
   *Inter*-molecular isopeptide bonds lack the asp/glu catalytic residue and therefore 
   can not be detected by isopeptor which has been designed to predict the presence of intramolecular isopeptide bonds
   by detecting all three residues involved in the formation of the bond.

Documentation
-------------

.. toctree::

   workflow
   installation
   commandline
   api
   reference




