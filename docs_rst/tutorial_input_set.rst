Tutorial -- input set
---------------------

In this tutorial, we show how to use :code:`vise` to generate the :code:`VASP` input files.

============================
Preparation of the unit cell
============================
Firstly, we obtain the POSCAR file via Materials Project REST API.
(Of course, it's also fine to prepare POSCAR by another way instead.)
When we use the Materials Project REST API,
we need to set the PMG_MAPI_KEY in the .pmgrc.yaml file at the home directory, e.g.,
See `pymatgen web page 1 <https://pymatgen.org/usage.html>`_, for more details.

By checking the Materials Project web page, we know the id for ScN is mp-2857.
:code:`Vise` automatically retrieve the POSCAR files
using the :code:`get_poscar` (= :code:`gp`) sub-command.
For example, we can get the crystal structure of ScN. by typing as follows:

::

    vise gp -m mp-2857

===================================
Input files for the cell relaxation
===================================
Secondly, we prepare INCAR, POTCAR, and KPOINTS files.
In :code:`vise`, :code:`vasp_set` (= :code:`vs`) sub-command generates these files.
It includes various arguments, in which important ones are
:code:`--task` (or :code:`-t`) and :code:`--xc` (or :code:`-x`),
which determine the task and exchange-correlation (XC) functional.
The defaults are structure optimization with the PBE functional.
So, we can generate INCAR, POTCAR, and KPOINTS files,
by typing as follows at the directory where POSCAR exists,

::

    vise vs

Note that the structure optimization must be generally iterated with 1.3 times larger cutoff energy
until the forces and stresses converge at the first ionic step so as to reduce the errors caused
by the Pulay Stress to an acceptable accuracy.
See `vasp manual <https://cms.mpi.univie.ac.at/vasp/vasp/Volume_vs_energy_volume_relaxations_Pulay_Stress.html>`_
or `wikipedia <https://en.wikipedia.org/wiki/Pulay_stress>`_ for details.
Such increase of the cutoff energy is also done by :code:`vise`.

The :code:`vs` sub-command allows the :code:`POTCAR` file to be modified
from the default set via :code:`--potcar`.
The :code:`INCAR` setting is also controlled via :code:`--user_incar_setting` (= :code:`-uis`) argument.

It is also possible to control the :code:`POTCAR` and :code:`INCAR` setting
via :code:`vise.yaml` file, but the command line arguments are prioritized.
See :doc:`tutorial_vise_yaml` for details.

Next, let's move to :doc:`tutorial_band_dos`.

=============
KPOINTS files
=============
The kpoint mesh is determined to be proportional to the reciprocal lattice constants.
Let's consider the cubic lattice with a=10A, where the reciprocal lattice length in the "physics" definition is 2*pi/10.
When the density of the k-point mesh is set to 2.5 A,
the number of k points along this direction is ceil(2*pi/10*2.5) = ceil(pi/2) = 2.
Users can also control this density via `vise.yaml`.

There is an exception for body centered orthorhombic and tetragonal systems.
In these, although distances of reciprocal lattice constants are not the same,
but the number of k-points needs to be the same to keep the lattice symmetry.

