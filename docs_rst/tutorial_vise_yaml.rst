Tutorial -- vise.yaml
---------------------

=========================
Setting for the vise.yaml
=========================
In :code:`vise`, users can control various types of parameters using :code:`vise.yaml` file.
For example, we provide the default :code:`POTCAR` set,
but users can override their favorite :code:`POTCAR` set as follows:

::

    overridden_potcar: Mg_pv O_h

The Mg_pv and O_h :code:`POTCAR` files are then used instead of the default normal Mg and O :code:`POTCAR` files.

Note that :code:`vise` tries to find the :code:`vise.yaml` files
from the current working directory to the parent directly up to the home or root directory.
If the same options exist in multiple `vise.yaml` files,
the value in the yaml file locating at the deeper directory is prioritized.
For example, when two different :code:`vise.yaml` are located
at the home directory of the project and :code:`$home/something/` directory,
the parameters written in :code:`$home/something/vise.yaml` are prioritized.

The properties of Defaults class in `defaults.py <https://github.com/kumagai-group/vise/blob/master/vise/defaults.py>`_ can be set as the keys
of :code:`vise.yaml`.
For example, one can write in :code:`vise.yaml` as

::

    xc: hse
    user_incar_settings:
        LASPH: False


for the XC functional and INCAR parameters.
