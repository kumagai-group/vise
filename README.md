![PyPI - License](https://img.shields.io/pypi/l/vise?color=blue)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/vise)
[![CircleCI](https://circleci.com/gh/kumagai-group/vise/tree/master.svg?style=shield)](https://circleci.com/gh/kumagai-group/vise/tree/master)

vise
=========
Vasp Integrated Simulation Environment (vise) is a collection of tools that 
supports VASP users to prepare input files for the VASP calculations, handle its calculation errors, and analyze its results.

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
    ```
    vise 
    ```

3. Usage as a module
  - vise can be imported as a python module
#  - The comments in the script makes it (hopefully) self explained.

#Input files for several materials can be found in the same example/inputs directory.
#More information in the online manual at: http://xxxx

Files and directories included in vise distribution
--------------------------------------------------------
~~~
  README             this file 
  LICENSE            the MIT license 
  setup.py           installation script
  requirements.txt   list of required packages

  /vise/analyzer     
  /vise/chempotdiag
  /vise/cli
  /vise/custodian_extension
  /vise/input_set
  /vise/test_files
  /vise/util
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

- Code style should comply with [PEP8](http://www.python.org/dev/peps/pep-0008) where possible. 
[Google's house style](https://google.github.io/styleguide/pyguide.html) is also helpful, including a good model for docstrings.
- Add tests wherever possible, and use the test suite to check if you broke anything.

### Tests
Run the tests using `pytest vise`.
We also use integrated testing on Github via [circleCI]().

Contributors
--------------------------------------------------------
Yu Kumagai
Akira Takahashi

Contact info
---------------------------------------------------------
Yu Kumagai
<br>yuuukuma@gmail.co.jp

Tokyo Institute of Technology
<br>Tokyo (Japan)

