Tutorial -- input set
---------------------

In this tutorial, the subparsers of :code:`get_poscars=(gp)` and :code:`vasp_set=(vs)` are introduced.

============================
Preparation of the unit cell
============================
Here, we show how to use vise using an example of ScN.
Firstly, we obtain the POSCAR file via Materials Project REST API.
(Of course, it's also fine to prepare POSCAR by another way by yourself instead.)
When we use the Materials Project REST API,
we need to set the PMG_MAPI_KEY in the .pmgrc.yaml file at the home directory, e.g.,
See `pymatgen web page 1 <https://pymatgen.org/usage.html>`_, for more details.

When we check the Materials Project web page, we know the id for ScN is mp-2857.
vise automatically retrieve the POSCAR files using the :code:`gp` (= :code:`get_poscar`) sub-command.
For example, we can get the crystal structure of ScN. by typing as follows,

::

    vise gp -m mp-2857

===================================
Input files for the cell relaxation
===================================
Let's begin with the relaxation of the unit cell using VASP.
For this purpose, we need to prepare INCAR, POTCAR, and KPOINTS files.
In vise, :code:`vs` (= :code:`vasp_set`) sub-option automatically generates these files.
:code:`vs` includes various arguments, and the most important ones are
:code:`--task` (or :code:`-t`) and :code:`--xc` (or :code:`-x`),
determining the task and exchange-correlation (XC) functional.
The defaults are structure optimization (structure_opt) for task and PBE functional (pbe) for XC functional.
So we can generate INCAR, POTCAR, and KPOINTS files, by typing as follows at the directory where POSCAR exists,

::

    vise vs

The :code:`vs` sub-option has a lot of arguments, such as :code:`--potcar` and :code:`--kpt_density` (or :code:`-k`).
Some users may want to use their favorite potcar set rather than the vise defaults set.
In this case, it is cumbersome to add :code:`--potcar` argument all the time.
To avoid such circumstance, users can use the :code:`vise.yaml` file.
See vise.yaml section for details.


Note that the structure optimization must be generally iterated with 1.3 times larger cutoff energy
until the forces and stresses converge at the first ionic step so as to reduce the errors caused by the Pulay Stress to an acceptable accuracy.
See [vasp manual](https://cms.mpi.univie.ac.at/vasp/vasp/Volume_vs_energy_volume_relaxations_Pulay_Stress.html) or [wikipedia](https://cms.mpi.univie.ac.at/vasp/vasp/Volume_vs_energy_volume_relaxations_Pulay_Stress.html) for details.
Such increase of the cutoff energy is done by :code:`vise`.

==============================
Advanced usage of the vasp_set
==============================
In the arguments of :code:`vs`, there are -uis(=user_incar_setting) and -auis(=additional_user_incar_setting) arguments.
The former INCAR setting is set by vise.yaml by default, and if one uses this option, the vise.yaml default setting is overwritten.
Conversely, the latter is set only via this option and it does not overwrite user_incar_setting written in vise.yaml
Please use both the options as the situation demands.







