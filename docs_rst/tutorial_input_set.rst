Tutorial -- input set
---------------------

In this tutorial, we show how to use :code:`vise` to generate the :code:`VASP` input files.

===================================
Input files for the cell relaxation
===================================
Let's assume we have POSCAR file are are about to prepare INCAR, POTCAR, and KPOINTS files
The used sub-command is the :code:`vasp_set` (= :code:`vs`).
It includes various arguments, among which important are
:code:`--task` (or :code:`-t`) and :code:`--xc` (or :code:`-x`)
that determine the target task and exchange-correlation (XC) functional,
respectively.
The defaults are structure optimization and the PBE functional.
So, we can generate INCAR, POTCAR, and KPOINTS files,
by typing as follows at the directory where POSCAR exists,

::

    vise vs

Note that the structure optimization must be generally iterated
with 1.3 times larger cutoff energy until the forces and stresses converge
at the first ionic step so as to reduce the errors caused
by the Pulay Stress to an acceptable accuracy.
See `vasp manual <https://cms.mpi.univie.ac.at/vasp/vasp/Volume_vs_energy_volume_relaxations_Pulay_Stress.html>`_
or `wikipedia <https://en.wikipedia.org/wiki/Pulay_stress>`_ for details.
Such increase of the cutoff energy is also considered in :code:`vise`.

The :code:`vs` sub-command allows the :code:`POTCAR` file to be modified
from the default set via :code:`--potcar`.
The :code:`INCAR` setting is also controlled via :code:`--user_incar_setting` (= :code:`-uis`) argument.
An example is

::

    vise vs vise vs --potcar Mg_pv -uis ALGO All

It is also possible to control the :code:`POTCAR` and :code:`INCAR` setting
via :code:`vise.yaml` file, but the command line arguments are prioritized.
See :doc:`tutorial_vise_yaml` for details.


=====
Tasks
=====

.. csv-table:: Tasks
   :file: task.csv
   :header-rows: 1

**: Not set, thus same as default.

**1: Number of bands is determined from the sum of valence electron / 2 and unoccupied band number per element written in `this yaml <https://github.com/kumagai-group/vise/blob/master/vise/input_set/datasets/unoccupied_bands.yaml>`_.

**2: Set -5 when there is a band gap and number of irreducible k-points >= 4, otherwise 0.

**3: Combination of IBRION=8 and NPAR>=2 is prohibited.

**4: The multiplication factor for the number of k-points. For example, the numbers of k-points for DOS are doubled along all the directions.

**5: For the band structure calculations, the structures must be primitive cells that are uniquely defined in the seekpath.

**6: See comments in `this script <https://github.com/kumagai-group/vise/blob/master/vise/input_set/kpoints_mode.py>`_.

Caveats:

- SIGMA=0.1, NELM=100, and LASPH=True are set in common.

- EMIN, EMAX, and NEDOS are set for dos and dielectric_function.
  When vbm and cbm are set, EMIN=vbm - 15 and EMAX=cbm + 15, otherwise EMIN=20 and EMAX=20.
  NEDOS = round((EMAX - EMIN) / 0.01) + 1

- Combination of ISMEAR=-5 and non-Gamma-centered kpoint mesh is prohibited, so the centering
  is changed to Gamma.

=============
XC functional
=============

.. csv-table:: XC functional
   :file: xc.csv
   :header-rows: 1

**1: The LDAUU and LDAUL parameters are written in the `yaml file <https://github.com/kumagai-group/vise/blob/master/vise/input_set/datasets/u_parameter_set.yaml>`_.

=============
POTCAR files
=============
Here, one can see the `POTCAR list <https://github.com/kumagai-group/vise/blob/master/vise/input_set/datasets/potcar_set.yaml>`_, where the first column shows the POTCAR used in vise as default and the second column those adopted in the materials project database.

=============
KPOINTS files
=============
The kpoint mesh is determined to be proportional to the reciprocal lattice constants.
Let's consider the cubic lattice with a=10Å,
where the reciprocal lattice length in the "physics" definition is 2π/10.
When the density of the k-point mesh is set to 2.5Å,
the number of k points along this direction is ceil(2*π/10*2.5) = ceil(π/2) = 2.
Users can also control this density via :code:`vise.yaml`.

There is an exception for body centered orthorhombic and tetragonal systems.
In these, although distances of reciprocal lattice constants are not the same,
the number of k-points needs to be the same to keep the lattice symmetry.
Therefore, we first calculate the geometric mean of the reciprocal lattice constants,
and estimate the averaged number of k-points based on it.

=================
prev_dir argument
=================
When one sets --prev_dir argument, vise parses the VASP calculations performed in that directory
to determine the structure, charge, band-edge positions, and magnetization.
If one wants to transfer files, use the --file_transfer argument together.

================
options argument
================
The options in `IncarSettingsGenerator <https://github.com/kumagai-group/vise/blob/master/vise/input_set/incar_settings_generator.py>`_ and `StructureKpointsGenerator <https://github.com/kumagai-group/vise/blob/master/vise/input_set/structure_kpoints_generator.py>`_
e.g.,

for IncarSettingsGenerator,

* charge: float = 0.0,
* band_gap: Optional[float] = None,
* vbm_cbm: Optional[List[float]] = None,
* exchange_ratio: float = 0.25,
* set_hubbard_u: Optional[bool] = None,
* auto_npar_kpar: bool = True,
* cutoff_energy: Optional[float] = None,
* is_magnetization: bool = False,

and for StructureKpointsGenerator,

* kpt_density: Optional[float] = None,  # in Å
* gamma_centered: Optional[bool] = None,  # Vasp definition
* only_even_num_kpts: bool = False,  # If ceil kpt numbers to be even.
* num_kpt_factor: Optional[int] = None,  # Set NKRED to this as well.

Next, let's move to :doc:`tutorial_band_dos`.


