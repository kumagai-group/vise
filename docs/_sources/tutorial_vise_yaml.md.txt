Tutorial -- Setting for the vise.yaml
-----------------------
In `VISE`, users can control various types of parameters using `vise.yaml` files.
For example, we provide the default `POTCAR` set regularly we use, but users can adopt their favorite `POTCAR` set,
by setting `potcar_set` key as follows
```
potcar_set: Mg_pv O_h
```
Then, the Mg_pv and O_h `POTCAR` files are used instead of the normal Mg and O `POTCAR` files.

Note that `VISE` try to find the `vise.yaml` file from the current working directory to the parent directly up to the home or root directory.
Therefore, when `vise.yaml` is located at the top directory of the project as shown above directory tree,
the parameters are always used for the `pydefect` commands.
The keys for `vise.yaml` are
- "vasp_cmd" (str)
- "symprec" (float)
- "angle_tolerance" (float)
- "xc" (str)
- "kpt_density" (float)
- "initial_kpt_density" (float)
- "vise_opts" (dict): 
- "user_incar_setting" (dict)
- "ldauu" (dict)
- "ldaul" (dict)
- "potcar_set"
- "potcar_set_name" (str)
- "relax_iter_num" (int)
- "removed_files" (list)

For vise_opts corresponds to the Options written in ViseInputSet class in vise.input_set/input_set.py

Same info is written at the top of the `main.py` file.
For example, one can write in `vise.yaml` 
```
xc: hse
```
for XC functional.


