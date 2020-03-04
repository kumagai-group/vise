Tutorial -- analyzer
-----------------------

`VISE` also provides the plotters of BS and DOS based on with `plot_band` (=`pb`) and `plot_dos` (=`pd`) sub-commands.
Type the following commands in `band/` and  `dos/`, respectively. 
```
vise pb -f band.pdf
```



```
vise pd -f dos.pdf
```
The band gap is also evaluated using `bg`(=`band_gap`) sub command.
```
vise bg -v vasprun.xml.finish -o OUTCAR.finish
```

For band structures, one can show two band structures in the same figure.


