Tutorial -- vise.yaml
---------------------

================================
Setting for the vise.yaml
================================
In :code:`vise`, users can control various types of parameters using :code:`vise.yaml` file.
For example, we provide the default :code:`POTCAR` set regularly we use, but users can adopt their favorite :code:`POTCAR` set,
by setting :code:`potcar_set` key as follows

.. code-block::

    potcar_set: Mg_pv O_h

Then, the Mg_pv and O_h :code:`POTCAR` files are used instead of the normal Mg and O :code:`POTCAR` files in vise.

Note that vise try to find the :code:`vise.yaml` file from the current working directory to the parent directly up to the home or root directory.
And once it is found only the file is read and others locating in the parent directories are ignored.
Therefore, when :code:`vise.yaml` is located only at the top directory of the project, the parameters are always used.

The keys for :code:`vise.yaml` are

=================== ======== =========================================================================================================================
name                type      explanation
=================== ======== =========================================================================================================================
vasp_cmd            str/list  Vasp commands
symprec             float     Distance precision in angstrom used for symmetry analysis.
angle_tolerance     float     Angle tolerance in angstrom used for symmetry analysis.
xc                  str       Exchange-correlation functional. Supported functions are seen in the help of :code:`vasp_set`.
kpt_density         float     K-point density in inverse of angstrom.
initial_kpt_density float     Initial k-point density in inverse of angstrom used for k-point convergence calculations.
vise_opts           dict      Options used for the ViseInputSet in vise. See docstrings of ViseInputSet for detailed options.
user_incar_setting  dict      User-specific INCAR settings.
ldauu               dict      LDAUU values in INCAR for each elements.
ldaul               dict      LDAUL values in INCAR for each elements.
outcar              str       OUTCAR file name, which should be OUTCAR.finish when using :code:`vasp_run` or :code:`kpt_conv`.
contcar             str       CONTCAR file name, which should be CONTCAR.finish when using :code:`vasp_run` or :code:`kpt_conv`.
vasprun             str       vasprun.xml file name, which should be vasprun.xml.finish when using :code:`vasp_run` or :code:`kpt_conv`.
procar              str       PROCAR file name, which should be PROCAR.finish when using :code:`vasp_run` or :code:`kpt_conv`.
potcar_set          dict      POTCAR names such as O_h and Mg_pv for each elements.
potcar_set_name     str       Default POTCAR set name. Supported sets are found in the help of :code:`vasp_set`.
max_relax_num       int       Max relaxation number for :code:`vasp_run` or :code:`kpt_conv` (see help of :code:`vasp_run`).
removed_files       list      Removed files for :code:`vasp_run` or :code:`kpt_conv` (see help of :code:`vasp_run`)..
left_files          list      Left files at the calculation directory for :code:`vasp_run` or :code:`kpt_conv` (see help of :code:`vasp_run`)..
timeout             int       Timeout for :code:`vasp_run` or :code:`kpt_conv` (see help of :code:`vasp_run`)..
=================== ======== =========================================================================================================================

For example, one can write in :code:`vise.yaml` as

.. code-block::

    xc: hse
    user_incar_setting:
        ENCUT: 550

