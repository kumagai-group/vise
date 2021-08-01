Tutorial -- preparation of POSCAR
---------------------------------

In this tutorial, we show how to prepare POSCAR file via the Materials Project (MP) database.

==============================
Preparation of the POSCAR file
==============================
Firstly, we obtain the :code:`POSCAR` file through the MP REST API.
(Of course, it's also fine to prepare POSCAR by another way.)
For this, we need to set the PMG_MAPI_KEY in the .pmgrc.yaml file, e.g.,
See `pymatgen web page 1 <https://pymatgen.org/usage.html>`_, for more details.

To confirm the pymatgen setting works properly, run the following python script snippet.
Note that it creates vasp files, so it would be better to move to a temporary directory.
If the :code:`VASP` files are not created, there should be a problem related to the pymatgen.

::

    from pymatgen.io.vasp.sets import VaspInputSet
    from pymatgen.core import Structure, Lattice
    from pymatgen.io.vasp.sets import MPRelaxSet
    s = Structure(Lattice.cubic(1), ["H", "He"], [[0.0]*3, [0.5]*3])
    vasp_set = MPRelaxSet(s)
    vasp_set.write_input(".")

Once we find the MP id (e.g., mp-2857 for ScN) via the MP web page,
:code:`Vise` allows one to automatically retrieve the POSCAR files
using the :code:`get_poscar` (= :code:`gp`) sub-command.
For example, we can get the crystal structure of ScN by typing as follows:

::

    vise gp -m mp-2857

