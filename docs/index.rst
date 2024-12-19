
Isopeptor Documentation
=======================
Python package for the detection of intamolecular isopeptide bonds in protein structures.  
The method described in **Isopeptor: a tool for detecting intramolecular isopeptide bonds in protein structures**.

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
only the template-match with the lowest `Root Mean Square Deviation (RMSD) <https://en.wikipedia.org/wiki/Root_mean_square_deviation>`_
is reported. 
Isopeptor also calculates the `relative Accessible Solvent Area (rASA) <https://en.wikipedia.org/wiki/Relative_accessible_surface_area>`_
We decided to include rASA calculation to take into consideration the fact that
intramolecular isopeptide bonds form only in a buiried hydrophobic environment.
RMSD and rASA are used by a `logistic regression model <https://en.wikipedia.org/wiki/Logistic_regression>`_ to 
predict the presence of **intramolecular isopeptide bonds**, which are assigned for probability values above 0.5.

.. image:: figures/f2.png
   :name: Figure 1. 
   :alt: Isopeptor workflow.
   :width: 600px

   Isopeptor workflow.

Geometric evaluation
--------------------

.. image:: figures/bond_angles.png
   :name: Figure 2. 
   :alt: Isopeptor workflow.
   :width: 300px

Isopeptor can optionally evaluate the quality of intramolecular isopeptide bonds. This is done using two metrics:

#. Bond lenght `Z-score <https://en.wikipedia.org/wiki/Standard_score>`_. Outliers are assigned for z-score values above 4 or below -4.
#. `Kernel Density Estimate (KDE) <https://en.wikipedia.org/wiki/Kernel_density_estimation>`_ likelihood of pseudo dihedral angles (pseudo φ, pseudo ψ, pseudo ω). 6 different KDE models are employed: one for each combination of of angle pairs  and isopeptide bond type. Outliers are assigned for values of KDE likelihood outside the 95th percentile distribution of values from our database.
    
.. note::
    Isopeptide bond angles have been named after the peptide-bond dihedral angles nomenclature: 
    pseudo φ is the angle around the bond between Asp/Asn\ :sub:`Cβ`\ and Asp/Asn\ :sub:`Cγ`\, ω between Asp/Asn\ :sub:`Cγ`\ 
    and Lys\ :sub:`Nζ`\ bond, and ψ between  Lys\ :sub:\`Nζ`\ and Lys\ :sub:`Cε`\ bond.

Output
------
Isopeptor output consists of the follwoing fields:

* *protein_name*
* *probability*: probability calculated with the logistic regression model. Intramolecular isopeptide bonds are assigned for values above 0.5. 
* *chain*: protein chain
* *r1_bond*: position of the first residue involved in the intramolecular isopeptide bond (lysine)
* *r_cat*: position of the intramolecular isopeptide bond catalytic residue (aspartate/glutamate)
* *r2_bond*: position of the second residue involved in the intramolecular isopeptide bond (asparagine/aspartate)
* *r1_bond_name*: name of the first residue involved in the intramolecular isopeptide bond (lysine)
* *r_cat_name*: name of the intramolecular isopeptide bond catalytic residue (aspartate/glutamate)
* *r2_bond_name*: name of the second residue involved in the intramolecular isopeptide bond (asparagine/aspartate)
* *bond_type*: CnaA-like or CnaB-like. This is assigned based on the type of the closest template.
* *rmsd*: RMSD (Å) with the closest template.
* *r_asa*: rASA (ranges between 0 and 1).
* *template*: name of the closest template.

And optional fields:

* *bond_length*: bond length calculated between Lys\ :sub:`Nζ`\ and Asp/Asn\ :sub:`Cγ`\. If  Lys\ :sub:`Nζ`\ atom is missing, isopeptor will attempt using Asn\ :sub:`Nδ`\ instead.
* *bond_length_zscore*: value of bond length Z-score.
* *bond_length_allowed*: True/False. False if Z-score is below -4 or above 4.
* *pseudo_phi*: value of pseudo ...
* *pseudo_psi*: 
* *pseudo_omega*: 
* *phi_psi_likelihood*:
* *phi_psi_allowed*:
* *omega_psi_likelihood*:
* *omega_psi_allowed*:
* *omega_phi_likelihood*:
* *omega_phi_allowed*:

Installation
------------

.. code-block:: console

    $ pip install isopeptor

Usage
-----

.. toctree::

   commandline
   api




