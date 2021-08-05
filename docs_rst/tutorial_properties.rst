Tutorial -- Calculations of various properties
----------------------------------------------
==================
List of properties
==================

.. csv-table:: properties
   :file: properties.csv
   :header-rows: 1

**1: The band path is determined based upon the
`seekpath code <https://www.materialscloud.org/work/tools/seekpath>`_,
so if one uses the plot for publication or presentation, please cite the related paper.

**2: For the carrier concentration, calculations for the density of states (DOS) or absorption coefficient would be fine.

====================
Plotter and analyzer
====================
:code:`Vise` provides plotters for band structures, DOS, and absorption coefficients with
:code:`plot_band` (= :code:`pb`), :code:`plot_dos` (= :code:`pd`), and :code:`plot_absorption` (= :code:`pa`) sub-commands.

In :code:`plot_absorption` sub-command, 
The -ckk option allows us to calculate the real part of the dielectric function explicitly from the imaginary part using the Kramers-Kronig transformation.
he complex shift η is then set via the --ita option.

The band-edge positions including the band gap and effective masses are also evaluated using the :code:`band_edge` (= :code:`be`) and :code:`effective_mass` (= :code:`em`) sub-commands, respectively.
The calculations of the effective masses use
`BoltzTrap2 <https://www.imc.tuwien.ac.at/forschungsbereich_theoretische_chemie/forschungsgruppen/prof_dr_gkh_madsen_theoretical_materials_chemistry/boltztrap2/>`_.
So, if one uses the effective masses, please cite the related paper.

The essential key parameter is the carrier concentration, and user need to
input exponential parts of base-10 in cm-3 (A typical value would be between 17 -- 20).
For example, carrier concentrations of 10**17 and 10**20 cm-3 are set as follows:

::

    vise em -c 17 19
