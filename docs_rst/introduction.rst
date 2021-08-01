=====================
Introduction of vise.
=====================

First-principles calculations are becoming increasingly prevalent in materials science.
However, preparing input files for the calculations and analyzing calculation
results are still not an easy task for researchers other than experts,
and human errors could come up.

To overcome the situation, we have developed :code:`vise` code,
which means *VASP Integrated Supporting Environment*.
:code:`Vise` supports :code:`VASP` users to generate its input files
for several tasks with suitable defaults, and allows the users to
tune some parameters depending on their purposes.

It also supports the analyses of the calculated results for electronic band structures, density of states (DOS), effective masses, and optical absorption coefficients.
Since simple APIs are provided, :code:`vise` can be easily implemented to python codes for automatic calculations.

**Note: Units used in vise are eV for energy and Å for length
following the vasp convention.**

