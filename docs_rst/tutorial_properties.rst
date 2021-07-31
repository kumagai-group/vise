Tutorial -- Calculations of various properties
---------------------------------------------

====================================
Band structure and density of states
====================================
We can create the input files for the calculations of the band structure (BS) and
the density of states (DOS) by simply typing

::

    vise vs -t band

and

::

    vise vs -t dos

In each directory, we run the vasp calculations.

Here, the band path is determined based upon the
`seekpath code <https://www.materialscloud.org/work/tools/seekpath>`_,
so if one uses the plot for publication or presentation, please cite the following paper.

- `Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, Band structure diagram paths based on crystallography, Comp. Mat. Sci. 128, 140 (2017). <https://www.sciencedirect.com/science/article/pii/S0927025616305110?via%3Dihub>`_ DOI: 10.1016/j.commatsci.2016.10.015.

:code:`Vise` also provides the plotters of BS and DOS based on with
:code:`plot_band` (= :code:`pb`) and :code:`plot_dos` (= :code:`pd`) sub-commands.
Type the following commands in the BS and DOS calculation directories,

::

    vise pb -f band.pdf

and

::

    vise pd -f dos.pdf

The band-edge positions including the band gap is also evaluated using
the :code:`be`(= :code:`band_edge`) sub-command.

::

    vise be

The effective masses are calculated using
`BoltzTrap2 <https://www.imc.tuwien.ac.at/forschungsbereich_theoretische_chemie/forschungsgruppen/prof_dr_gkh_madsen_theoretical_materials_chemistry/boltztrap2/>`_.
The subcommand is :code:`em`(= :code:`effective_mass`) sub-command.
The essential key parameter is the carrier concentration, and user need to
input exponential parts of base-10 in cm-3 (A typical value is between 17 -- 20).

::

    vise em -c 17 19

====================================
Static and ionic dielectric constant
====================================

We can also create the input files for the calculation of the static and ionic dielectric constant.

::

    vise vs -t dielectric_dfpt

Instead, when one wants to use the finite field method,


::

    vise vs -t dielectric_finite_field


======================
Absorption coefficient
======================

We can also create the input files for the calculation of the absorption coefficient.

::

    vise vs -t dielectric_function

::

    vise pa
