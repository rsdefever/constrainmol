
constrainmol
============
|License|
|CodeCov|
|Azure|

.. |Codecov| image:: https://codecov.io/gh/rsdefever/constrainmol/branch/main/graph/badge.svg?token=H7BBUYGNTU
             :target: https://codecov.io/gh/rsdefever/constrainmol 
.. |Azure| image:: https://dev.azure.com/rdefever/constrainmol/_apis/build/status/rsdefever.constrainmol?branchName=main
.. |License| image:: https://img.shields.io/github/license/rsdefever/constrainmol

Overview
~~~~~~~~

A package to update the coordinates of molecular systems to match bond-length constraints. The package
is designed to work with the `MoSDeF tools <https://mosdef.org>`_
and `ParmEd <https://parmed.github.io/ParmEd/html/index.html#>`_.


Warning
~~~~~~~

**constrainmol** is still in early development (0.x releases). The API may
change unexpectedly.

Usage
~~~~~

**constrainmol** usage is encapsulated in a single ``ConstrainedMolecule``
class. Here we demonstrate a simple example: load a single molecule
from file, create the constrained molecule, solve for the new coordinates,
update the coordinates of molecule to their constrained values, and save
the result to file.

.. code-block:: python

  from constrainmol import ConstrainedMolecule
  molecule = parmed.load_file("molecule.top", xyz="molecule.gro")  # load a molecule w/ bond length info
  constrained_mol = ConstrainedMolecule(molecule)  # create the constrained molecule
  constrained_mol.solve()  # solve for the constrained coordinates
  molecule.coordinates = constrained_mol.xyz  # update the coordinates to their constrained values
  molecule.save("constrained_molecule.gro")  # save the new coordinates to disk


If there is a system with many molecules (with the same desired geometry),
the coordinates can be updated and repeatedly solved. This is faster
than creating a new ``ConstrainedMolecule`` every time.

.. code-block:: python

    constrained_mol = ConstrainedMolecule(molecule)
    constrained_mol.solve()
    constrained_mol.update_xyz(new_coordinates)
    constrained_mol.solve()


More examples are provided in the ``examples/`` directory of this repository.

Installation
~~~~~~~~~~~~

Installation from source is the only supported option at this time:

.. code-block:: bash

    git clone git@github.com/rsdefever/constrainmol.git
    cd constrainmol/
    conda create --name constrain --file requirements.txt -c conda-forge -c mosdef -c omnia
    pip install .

If you want to test your installation you can replace the final two steps above with the following:

.. code-block:: bash

    conda create --name constrain --file requirements-dev.txt -c conda-forge -c mosdef -c omnia
    pip install .
    pytest constrainmol/tests


Credits
~~~~~~~

Development of constrainmol was supported by the National Science Foundation
under grant NSF Grant Number 1835874. Any opinions, findings, and conclusions or
recommendations expressed in this material are those of the author(s) and do
not necessarily reflect the views of the National Science Foundation.
