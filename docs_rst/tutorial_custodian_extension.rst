Tutorial -- custodian extension
-------------------------------

In this tutorial, the subparsers of :code:`vasp_run=(vr)` and :code:`kpt_conv=(kc)` are introduced.

========================================
VASP calculations with custodian wrapper
========================================
We assume that there are four vasp input files, namely INCAR, POSCAR, POTCAR, KPOINTS.
To generate these input files, please refer to tutorial of vise input set.

Here, We explain how to run VASP calculations with custodian.
For this purpose, vise provides an extension utility of :code:`custodian`.
One can run the typical structure optimization using the :code:`vr`(=:code:`vasp_run`) sub-option.
An example is shown as
.. code-block::

    vise vr -v mpirun -np 16 vasp_std

We provide some error handler groups, e.g., minimum, default, dielectric, and no_handler.
Note that the original error handlers are slightly modified and some handlers are added.
If you do not want to use custodian error handler, please set no_handler to -handler_name option.
To know the details, please look into the custodian_extension directory in vise.

In most cases, we run VASP with cluster nodes via submitting runshell script.
An example of runshell script is shown below.
.. code-block::

    #!/bin/zsh
    #$ -S /bin/zsh
    #$ -cwd
    #$ -V
    #$ -j y
    #$ -N pydefect
    #$ -o std.log
    #$ -pe all_pe* 36
    #============ Shell Script ============

    vise vr -v mpirun -np 16 vasp_std

The :code:`vr` sub-option is used also for static calculations such as band structure and dielectric constant.

In order to check the k-point convergence, one can use :code:`kc(=kpt_conv)` options.
In such cases, we can set various options such as -criteria.
.. code-block::

    vise vr --print


