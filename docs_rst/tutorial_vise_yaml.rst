Tutorial -- vise.yaml
---------------------

===============================
Setting for the vise.yaml
===============================
In :code:`vise`, users can control various types of parameters using :code:`vise.yaml` file.
For example, we provide the default :code:`POTCAR` set regularly we use,
but of course, users can adopt their favorite :code:`POTCAR` set by setting :code:`potcar_set` key as follows

::

    potcar_set: Mg_pv O_h

The Mg_pv and O_h :code:`POTCAR` files are then used instead of the default normal Mg and O :code:`POTCAR` files.
Note that :code:`vise` try to find the :code:`vise.yaml` files
from the current working directory to the parent directly up to the home or root directory.
If the same options exist in multiple `vise.yaml` files,
the value written in the yaml file at deeper directory is prioritized.
Therefore, when :code:`vise.yaml` is located at the home directory of the project
and :code:`$home/something/` directory,
the parameters written in :code:`$home/something/vise.yaml` are
always used for the :code:`vise` commands.

The keys for :code:`vise.yaml` are shown in defaults.py in :code:`vise` package.
For example, one can write in :code:`vise.yaml` as

::

    xc: hse

for the XC functional.