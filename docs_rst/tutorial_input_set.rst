Tutorial -- input set
---------------------

We here show how to generate the :code:`VASP` input files.

===================================
Input files for the cell relaxation
===================================
Let's assume we have a POSCAR file and are about to prepare INCAR, POTCAR, and KPOINTS files.
The sub-command for this is the :code:`vasp_set` (= :code:`vs`).
It includes various arguments, but important are
:code:`--task` (or :code:`-t`) and :code:`--xc` (or :code:`-x`)
that determine the target task and exchange-correlation (XC) functional adopted,
respectively.
The defaults are structure optimization and the PBE functional.
So, we can generate INCAR, POTCAR, and KPOINTS files by typing as follows at the directory where POSCAR exists,

::

    vise vs


The :code:`vs` sub-command allows the :code:`POTCAR` file to be modified
from the default set via :code:`--potcar`.
The :code:`INCAR` setting is also controlled via :code:`--user_incar_setting` (= :code:`-uis`) argument.
An example is

::

    vise vs --potcar Mg_pv -uis ALGO All

It is also possible to control the :code:`POTCAR` and :code:`INCAR` setting
via :code:`vise.yaml` file, but the command line arguments are prioritized.
See :doc:`tutorial_vise_yaml` for details.


=====
Tasks
=====
Here, we show a list of INCAR settings and some conditions for input files for each task. The VASP defaults at August, 1, 2021 are also shown.
Note that the structure optimization must be generally iterated
with 1.3 times larger cutoff energy until the forces and stresses converge
at the first ionic step so as to reduce the errors caused
by the Pulay Stress to an acceptable accuracy.
See `wikipedia <https://en.wikipedia.org/wiki/Pulay_stress>`_ for details.
Such increase of the cutoff energy is also considered in :code:`vise`.

.. csv-table:: Tasks
   :file: task.csv
   :header-rows: 1

**: Not set, thus same as default.

**1: Number of bands is determined from the sum of valence electron / 2 and unoccupied band number per element written in `this yaml <https://github.com/kumagai-group/vise/blob/master/vise/input_set/datasets/unoccupied_bands.yaml>`_.

**2: Set ISMEAR=-4 or -5 when there is a band gap and number of irreducible k-points is 4 or more, otherwise ISMEAR=0.

**3: Combination of IBRION=8 and NPAR>=2 is prohibited.

**4: The multiplication factor for the number of k-points. For example, the numbers of k-points for DOS are doubled along all the directions.

**5: For the band structure calculations, the structures must be primitive cells that are uniquely defined in the seekpath.

**6: See comments in `this script <https://github.com/kumagai-group/vise/blob/master/vise/input_set/kpoints_mode.py>`_.

Caveats:

- SIGMA=0.1, NELM=100, and LASPH=True are set in common.

- EMIN, EMAX, and NEDOS are set for dos and dielectric_function tasks.
  When vbm and cbm are given, EMIN=vbm - 15 and EMAX=cbm + 15, otherwise EMIN=20 and EMAX=20.
  NEDOS = round((EMAX - EMIN) / 0.01) + 1

- Combination of ISMEAR=-5 and non-Gamma-centered kpoint mesh is prohibited, so the centering
  is changed to Gamma when ISMEAR=-5.

=============
XC functional
=============
We next show a list of INCAR settings and some conditions for input files for each XC functional.

.. csv-table:: XC functional
   :file: xc.csv
   :header-rows: 1

**1: The LDAUU and LDAUL parameters are written in the `yaml file <https://github.com/kumagai-group/vise/blob/master/vise/input_set/datasets/u_parameter_set.yaml>`_.

=============
POTCAR files
=============
One can find the `POTCAR list <https://github.com/kumagai-group/vise/blob/master/vise/input_set/datasets/potcar_set.yaml>`_, where the first column shows the POTCAR used in vise as default and the second column those adopted in the materials project database.

=============
KPOINTS files
=============
The kpoint mesh is determined to be proportional to the reciprocal lattice constants.
Let's consider the cubic lattice with a=10Å,
where the reciprocal lattice length in the "physics" definition is 2π/10.
When the density of the k-point mesh is set to 2.5Å,
the number of k points along this direction is ceil(2*π/10*2.5) = ceil(π/2) = 2.

Body centered orthorhombic and tetragonal systems are exceptions;
although distances of reciprocal lattice constants are not the same,
the number of k-points needs to be the same to keep the lattice symmetry.
Therefore, we first calculate the geometric mean of the reciprocal lattice constants,
and estimate the average number of k-points based on it.

===================
--prev_dir argument
===================
The --prev_dir argument allows for parsing the VASP calculations performed in the designated directory
to extract information on the structure, charge, band-edge positions, and magnetization.
If one wants to copy, move, and/or link files from the directory, use the --file_transfer argument together.

==================
--options argument
==================
The options in `IncarSettingsGenerator <https://github.com/kumagai-group/vise/blob/master/vise/input_set/incar_settings_generator.py>`_
and `StructureKpointsGenerator <https://github.com/kumagai-group/vise/blob/master/vise/input_set/structure_kpoints_generator.py>`_
classes are set with the --options argument.
For example, those for IncarSettingsGenerator are,

* charge: float = 0.0,
* band_gap: Optional[float] = None,
* vbm_cbm: Optional[List[float]] = None,
* exchange_ratio: float = 0.25,
* set_hubbard_u: Optional[bool] = None,
* auto_npar_kpar: bool = True,
* cutoff_energy: Optional[float] = None,
* is_magnetization: bool = False,

and those for StructureKpointsGenerator are,

* kpt_density: Optional[float] = None,  # in Å
* gamma_centered: Optional[bool] = None,
* only_even_num_kpts: bool = False,  # Set when ceiling kpt numbers to be even.
* num_kpt_factor: Optional[int] = None,  # NKRED is set to this as well.

::

    vise vs --options cutoff_energy 1000 only_even_num_kpts True

Next, let's move on to :doc:`tutorial_band_dos`.


