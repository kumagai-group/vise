![PyPI - License](https://img.shields.io/pypi/l/vise?color=blue)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/vise)
[![CircleCI](https://circleci.com/gh/kumagai-group/vise/tree/master.svg?style=shield)](https://circleci.com/gh/kumagai-group/vise/tree/master)

vise
=========
VASP Integrated Supporting Environment (vise) is a collection of tools that 
supports VASP users to prepare input files for the VASP calculations, handle its calculation errors, and analyze its results.

**Note1: Only ensure VASP ver5.4.4 so far.**

**Note2: Units used in pydefect are eV for energy and angstrom for length following the vasp convention.**

**Note3: When transforming the structure to the primitive one, antiferromagnetic magnetization is not supported.**

Installation instructions
---------------------------------------------------------
1. Requirements
  - Python 3.7 or higher

Vise depends largely on the following packages, which should be acknowledged sincerely,
  - pymatgen
  - spglib
  - seekpath
  - see requirements.txt for others

2. Latest stable version is released at PyPI repository, so one can download 
it using `pip install vise`.


Executing this software
---------------------------------------------------------

1. Command line method
  - execute ***vise -h*** for detailed description of available options

2. Usage as a module
  - vise can be imported as a python module

More information in the online manual at: https://kumagai-group.github.io/vise/

Files and directories included in vise distribution
--------------------------------------------------------
~~~
  README            : introduction
  LICENSE           : the MIT license 
  setup.py          : installation script
  requirements.txt  : list of required packages

  /vise/analyzer    : tools for VASP analysis especially for band figure and density of states
  /vise/cli         : command line interfaces
  /vise/input_set   : tools for generating VASP input files
  /vise/tests       : test files used mainly for unitests
  /vise/util        : useful tools 
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
We also use integrated testing on Github via circleCI.

Citing vise
---------------
If vise has been used in your research, please consider citing our Github repo!

Contact info
------------
Yu Kumagai<br>
yuuukuma@gmail.com

Tokyo Institute of Technology (Japan)

