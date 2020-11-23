Tutorial -- vise.yaml
---------------------

================================
Setting for the vise.yaml
================================
In :code:`vise`, users can control various types of parameters using :code:`vise.yaml` file.
For example, we provide the default :code:`POTCAR` set regularly we use, but users can adopt their favorite :code:`POTCAR` set,
by setting :code:`potcar_set` key as follows

::

    overridden_potcar:
        Mg_pv O_h

Then, the Mg_pv and O_h :code:`POTCAR` files are used instead of the normal Mg and O :code:`POTCAR` files in vise.

Note that vise try to find the :code:`vise.yaml` file from the current working directory to the parent directly up to the home or root directory.
And, multiple files are read with the priority.

The keys for :code:`vise.yaml` are shown in defaults.py

For example, one can write in :code:`vise.yaml` as

::

    xc: hse
    user_incar_settings:
        ENCUT: 550

