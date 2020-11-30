Tutorial -- vise.yaml
---------------------

===============================
Setting for the vise.yaml
===============================
In :code:`vise`, users can control various types of parameters using :code:`vise.yaml` file.
For example, we provide the default :code:`POTCAR` set regularly we use,
but of course, users can adopt their favorite :code:`POTCAR` set by setting :code:`potcar_set` key as follows

::

    overridden_potcar: Mg_pv O_h

The Mg_pv and O_h :code:`POTCAR` files are then used instead of the default normal Mg and O :code:`POTCAR` files.
(The default :code:`POTCAR` is described in :code:`potcar_set.yaml` file.
Note that :code:`vise` tries to find the :code:`vise.yaml` files
from the current working directory to the parent directly up to the home or root directory.
If the same options exist in multiple `vise.yaml` files,
the value written in the yaml file at deeper directory is prioritized.
Therefore, when two different :code:`vise.yaml` are located
at the home directory of the project and :code:`$home/something/` directory,
the parameters written in :code:`$home/something/vise.yaml` are prioritized.

The properties of Defaults class in :code:`defaults.py` can be set as the keys
of :code:`vise.yaml`.
For example, one can write in :code:`vise.yaml` as

::

    xc: hse
    potcar_set: gw

for the XC functional and :code:`POTCAR` default set.