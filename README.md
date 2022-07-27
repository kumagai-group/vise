![PyPI - License](https://img.shields.io/pypi/l/vise?color=blue)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/vise)
[![CircleCI](https://circleci.com/gh/kumagai-group/vise/tree/master.svg?style=shield)](https://circleci.com/gh/kumagai-group/vise/tree/master)

vise
=========
VASP Integrated Supporting Environment (vise) is a collection of tools that 
supports VASP users to prepare input files for the VASP calculations, handle its calculation errors, and analyze its results.

**Note1: Only ensure VASP ver5.4.4 or later so far.**

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

Detailed information is provided in the online manual at: https://kumagai-group.github.io/vise/

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
If vise has been used in your research, please consider citing the following paper.

["Insights into oxygen vacancies from high-throughput first-principles calculations"](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.5.123803)<br>
Yu Kumagai, Naoki Tsunoda, Akira Takahashi, and Fumiyasu Oba<br>
Phys. Rev. Materials 5, 123803 (2021)

Contact info
------------
Yu Kumagai<br>
yuuukuma@gmail.com

Tohoku University (Japan)

