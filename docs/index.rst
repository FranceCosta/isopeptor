
Isopeptor Documentation
=======================
Python package for the detection of intamolecular isopeptide bonds in protein structures.  
The method described in "Isopeptor: a tool for detecting intramolecular isopeptide bonds in protein structures".

Intramolecular Isopeptide bonds
-------------------------------
Intramolecular isopeptide bonds typically form between lysine and asparagine/aspartate residues, 
typically catalyzed by a nearby aspartate or glutamate. 
When found in β-sandwich folds, they occur in CnaA domains (linking β-strands of opposing β-sheets) and 
CnaB domains (linking adjacent β-strands of the same β-sheet). 
They stabilize protein structures against various stresses.

Input
-----
Protein structures in PDB or CIF format.

.. warning::
    CIF file parsing is handled by CIF to PDB conversion at the moment. 
    This is not suitable for very large protein structures.

Prediction
----------
Isopeptor employs `pyjess <https://pypi.org/project/pyjess/>`_ to perform a template-based search of input protein structures
from a database of 140 high resolution intramolecular isopeptide bond structures. For each triad of residues detected by pyjess,
only the template-match with the lowest `Root Mean Square Deviation (RMSD) <https://en.wikipedia.org/wiki/Root_mean_square_deviation>`
is reported. 
Isopeptor also calculates the `relative Accessible Solvent Area (rASA) <https://en.wikipedia.org/wiki/Relative_accessible_surface_area>`_
We decided to include rASA calculation to take into consideration the fact that
 intramolecular isopeptide bonds form only in a hydrophobic environment.
RMSD and rASA are used by a `logistic regression model <https://en.wikipedia.org/wiki/Logistic_regression>`_ to 
predict the presence of an **intramolecular isopeptide bond**.


Installation
------------

.. code-block:: console

    $ pip install isopeptor

Usage
-----

.. toctree::

   commandline
   api




