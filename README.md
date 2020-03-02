
[![License](https://img.shields.io/apm/l/vise)](https://img.shields.io/apm/l/vise)
[![PyPI version]()]()
[![Build Status]()]()
[![Coverage Status]()]()

vise
=========
Vasp Integrated Simulation Environment (vise) is a collection of tools that


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

3. Scripting method (as a module)
  - vise can be imported as a python module
  - In examples/api_scripts directory an example script is available (script_silicon.py)
  - The comments in the script makes it (hopefully) self explained.

Input files for several materials can be found in the same example/inputs directory.
More information in the online manual at: http://xxxx

Files and directories included in DynaPhoPy distribution
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
### Bugs, features and questions
Please use the [Issue Tracker](https://github.com/kumagai-group/vise/issues) to report bugs or request features in the first instance. 
While we hope that most questions can be answered by searching [the docs](https://smact.readthedocs.io/en/latest/), 
we welcome new questions on the issue tracker, especially if they helps us improve the docs! 
For other queries about any aspect of the code, please contact Yu Kumagai by e-mail: yuuukuma@gmail.com.

### Code contributions
We are always looking for ways to make vise better and more useful to the wider community; contributions are very welcome. 
Please use the ["Fork and Pull"](https://guides.github.com/activities/forking/) workflow to make contributions and stick as closely as possible to the following:

- Code style should comply with [PEP8](http://www.python.org/dev/peps/pep-0008) where possible. 
[Google's house style](https://google.github.io/styleguide/pyguide.html) is also helpful, including a good model for docstrings.
- Please use comments liberally when adding nontrivial features, and take the chance to clean up other people's code while looking at it.
- Add tests wherever possible, and use the test suite to check if you broke anything.

### Tests
Run the tests using `pytest vise`.
We also use integrated testing on Github via [circleCI]().

Contributors
--------------------------------------------------------
Akira Takahashi


Contact info
---------------------------------------------------------
Yu Kumagai
<br>yuuukuma@gmail.co.jp

Tokyo Institute of Technology
<br>Tokyo (Japan)

