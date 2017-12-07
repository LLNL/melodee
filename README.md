# M E L O D E E
Modular Expression Language for Ordinary Differential Equation Editing
----------------------------------------------------------------------

Author: Robert Clayton Blake III, Lawrence Livermore National Laboratory <blake14@llnl.gov>

Melodee is a new language for ODE modeling and development. It has been designed to make it easy for domain experts to write modular equation systems.  Using a familiar Matlab-like syntax, modelers can describe their systems of interests and combine them in new and novel ways at run time.  Existing models can easily be converted to MELODEE because of its familiar semantics and syntax.

Included in this package are the following files:
* Parser.py
  * A parser for the language
  * Code-generation tools that can be used to write MELODEE source-to-source translation tools (i.e. convert a MELODEE model to code that solves the equations in a given software suite like Matlab or CVODE)
* Units.py
  * Python module handling unit conversions
* Cellml_converter.py
  * A CellML conversion tool that converts XML-based [CellML models](http://cellml.org) to the MELODEE syntax
* TT06.mel
  * A MELODEE file describing TT06, converted from the [following CellML model](https://models.physiomeproject.org/exposure/a7179d94365ff0c9c0e6eb7c6a787d3d/ten_tusscher_model_2006_IK1Ko_endo_units.cellml/view). The generated MELODEE code was then cleaned up and modified to show off more MELODEE features.

Requirements:
* Python - language everything is written in
* PLY - for LR parsing
* sympy - for symobolic expression manipulation

To get started, create a virtual environment and install the package.

```bash
virtualenv venv
source venv/bin/activate
pip install -e .
```

... then run the simulator
```bash
melodee -g matlab -t modifiedModel -t bclDefault HH.mel domains/electrophysiology/singleCell.mel
```

-----------------

This work was performed under the auspices of the United States Department of Energy by the Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344. Release Number LLNL-CODE-720003.

