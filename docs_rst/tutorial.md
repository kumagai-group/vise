# Tutorial  pydefect
-----------------------

### 1. Preparation of the unit cell
Here, we show how to use vise using an example of ScN.

Please create ScN directory and unitcell subdirectory in it, and enter there.
```
ScN
 │
 ├ vise.yaml
 │
 ├ unitcell/ ── structure_opt/
 │            ├ band/
 │            └ dos/
 │            
 └ optical_absorption
             ...
```

Firstly, we obtain the POSCAR file via Materials Project REST API.
(Of course, it's also fine to prepare POSCAR by another way by yourself instead.)
When we use the Materials Project REST API,
we need to set the PMG_MAPI_KEY in the .pmgrc.yaml file at the home directory, e.g.,
```
PMG_MAPI_KEY: xxxxxxxxxxxxxxxx
```
See [pymatgen web page 1](https://pymatgen.org/usage.html) or [2](https://pymatgen.org/_modules/pymatgen/io/vasp/inputs.html), for more details.

When we ckeck the Materials Project web page, we can know the id for ScN is mp-2857.
Therefore, by typing as follows,
```
python $PATH_TO_VISE/vise/vise/main.py gp -n 2857 
```
we can get the crystal structure of ScN.


### 2. Preparation of the input files for the unit cell relaxation.
We next relax the unit cell using VASP.
Please create `structure_opt`

For this purpose, we need to prepare INCAR, POTCAR, KPOINTS files.
In vise, `vs`(=`vasp_set`) automatically generate these files.
The default task and exchange-correlation (XC) functional of `vs` are structure optimization (structure_opt) and PBEsol functional (pbesol).
So, by typing as follows, at the directory 
```
python $PATH_TO_VISE/vise/vise/main.py vs
```

In options, there are -uis(=user_incar_setting) and -auis(=additional_user_incar_setting).
The former is set by vise.yaml by default, and if one uses this option, the vise.yaml default setting is overwritten.
Conversely, the latter is set only via this option.
Please use both the options as the situation demands.

### 3. Calculation of the unit cell relaxation.
We then run the vasp.
For this purpose, vise provides an extension utility of `custodian`.

`vr`(=`vasp_run`)
```
python $PATH_TO_VISE/vise/vise/main.py vr -v "mpirun -np xx vasp_std"
```

We provide four error handler groups, i.e., "minimum", "default", "dielectric", and "no_handler".
If you do not want to use custodian error handler, please set "no_handler" to -handler_name option.
See also vise.custodian_extension.handler_groups.py for details.

In most cases, we run VASP with cluster nodes via submitting runshell script.
An example of runshell script is shown below.
```
#!/bin/zsh
#$ -S /bin/zsh
#$ -cwd
#$ -V
#$ -j y
#$ -N pydefect
#$ -o std.log
#$ -pe all_pe* 36
#============ Shell Script ============

python ~/my_bin/vise/vise/main.py vr
```

In order to perform the k-point convergence, please use `-kc` options.
In such cases, we can set various options such as -criteria.


```
python $PATH_TO_VISE/vise/vise/main.py vr --print
```


### 4. Calculation of the band structure and the density of states
We then run the Calculation of the band structure and density of states.
Please create band and dos directories at the unitcell, and type 
```
python $PATH_TO_VISE/vise/vise/main.py vs -t band
```
```
python $PATH_TO_VISE/vise/vise/main.py vs -t dos
```
in each directory and run the vasp calculations as above.

Here, the band path is determined based upon the [seekpath code](https://www.materialscloud.org/work/tools/seekpath), 
so if one uses the plot for publication or presentation, please cite the following paper.
- [Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, Band structure diagram paths based on crystallography, Comp. Mat. Sci. 128, 140 (2017).](https://www.sciencedirect.com/science/article/pii/S0927025616305110?via%3Dihub) 
  DOI: 10.1016/j.commatsci.2016.10.015 (the "HPKOT" paper; arXiv version: arXiv:1602.06402).

`VISE` also provides the plotters of BS and DOS based on with `plot_band` (=`pb`) and `plot_dos` (=`pd`) sub-commands.
Type the following commands in `band/` and  `dos/`, respectively. 
```
python $PATH_TO_VISE/vise/vise/main.py pb -f band.pdf
```
```
python $PATH_TO_VISE/vise/vise/main.py pb -f dos.pdf
```

The band gap is also evaluated using `bg`(=`band_gap`) sub command.
```
python $PATH_TO_VISE/vise/vise/main.py bg -v vasprun.xml.finish -o OUTCAR.finish
```


### 5. Calculation of the dielectric constant

### Setting for the vise.yaml
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
- "vasp_cmd"
- "symprec"
- "angle_tolerance"
- "xc"
- "kpt_density"
- "initial_kpt_density"
- "vise_opts"
- "user_incar_setting"
- "ldauu"
- "ldaul"
- "potcar_set"
- "potcar_set_name"
- "relax_iter_num"
- "removed_files"

Same info is written at the top of the `main.py` file.
For example, one can write in `vise.yaml` 
```
xc: hse
```
for XC functional.


