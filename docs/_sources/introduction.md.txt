![PyPI - License](https://img.shields.io/pypi/l/vise?color=blue)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/vise)
[![CircleCI](https://circleci.com/gh/kumagai-group/vise/tree/master.svg?style=shield)](https://circleci.com/gh/kumagai-group/vise/tree/master)

vise
=========
Vasp Integrated Simulation Environment (vise) is a collection of tools that 
supports VASP users to prepare input files for the VASP calculations, handle its calculation errors, and analyze its results.

**Note: Units used in pydefect are eV for energy and angstrom for length following the vasp convention.**

Installation instructions
---------------------------------------------------------
1. Requirements
  - Python 3.6 or higher
  - pymatgen
  - matplotlib
  - spglib
  - seekpath
  - custodian
  - PyYAML
  - monty
  - numpy
  - scipy
  - tabulate
  - palettable


2. Latest stable version is released at PyPI repository, so one can download 
it using `pip install vise`.


Executing this software
---------------------------------------------------------

1. Command line method
  - execute ***vise -h*** for detailed description of available options

3. Usage as a module
  - vise can be imported as a python module
<!--
Input files for several materials can be found in the same example/inputs directory.
More information in the online manual at: http://xxxx
-->
Files and directories included in vise distribution
--------------------------------------------------------
~~~
  README                    : introduction
  LICENSE                   : the MIT license 
  setup.py                  : installation script
  requirements.txt          : list of required packages

  /vise/analyzer            : tools for VASP analysis especially for band figure and density of states
  /vise/chempotdiag         : tools for drawing the chemical potential diagram
  /vise/cli                 : command line interfaces
  /vise/custodian_extension : original custodian exention
  /vise/input_set           : tools for generating VASP input files
  /vise/test_files          : test files used mainly for unitests
  /vise/util                : useful tools 
~~~~

License
-----------------------
Python code is licensed under the MIT License.

Development notes
-----------------
### Bugs, requests and questions
Please use the [Issue Tracker](https://github.com/kumagai-group/vise/issues) to report bugs, request features.

### Code contributions
We are always welcome people who want to make vise better.
Please use the ["Fork and Pull"](https://guides.github.com/activities/forking/) workflow to make contributions and stick as closely as possible to the following:

- Code style follows [PEP8](http://www.python.org/dev/peps/pep-0008) and [Google's writing style](https://google.github.io/styleguide/pyguide.html).
- Add unittests wherever possible including scripts for command line interfaces.

### Tests
Run the tests using `pytest vise`.
We also use integrated testing on Github via [circleCI]().

Contributors
--------------------------------------------------------
Yu Kumagai<br>
Akira Takahashi

Contact info
---------------------------------------------------------
Yu Kumagai<br>
yuuukuma@gmail.co.jp

Tokyo Institute of Technology (Japan)

