Tutorial -- vise API
---------------------------------

In this tutorial, we show how to use vise in the python program.

==============================
Generation of VASP input files
============================== 

Here, an example to generate the VAPS input files is shown.

::

    from pathlib import Path
    
    from pymatgen.core import Structure, Lattice
    
    from vise.input_set.input_options import CategorizedInputOptions
    from vise.input_set.vasp_input_files import VaspInputFiles
    from vise.input_set.task import Task
    from vise.input_set.xc import Xc
    
    structure = Structure(Lattice.cubic(1), ["Mg", "O"], [[0.0]*3, [0.5]*3])
    
    categorized_input_options = CategorizedInputOptions(
                structure=structure,
                task=Task.band,
                xc=Xc.pbe, 
                overridden_potcar={"Mg": "Mg_pv"})
    
    input_files = VaspInputFiles(categorized_input_options, overridden_incar_settings={"NSW": 20})
    input_files.create_input_files(dirname=Path("."))

The CategorizedInputOptions class constructor takes the keyword arguments that 
are arguments of 
`generate_potcar <https://github.com/kumagai-group/vise/blob/master/vise/input_set/potcar_generator.py>`_,
`IncarSettingsGenerator <https://github.com/kumagai-group/vise/blob/master/vise/input_set/incar_settings_generator.py>`_,
and `StructureKpointsGenerator <https://github.com/kumagai-group/vise/blob/master/vise/input_set/structure_kpoints_generator.py>`_,
An example is overridden_potcar shown above.

The VaspInputFiles class constructor also takes the overridden_incar_settings, which can control the INCAR tags.

Note also that the vise.yaml files are also parsed.


Here, we show an example of FireTask in `FireWorks <https://materialsproject.github.io/fireworks/index.html>`_.


::

    from pathlib import Path
    
    from fireworks import FiretaskBase, explicit_serialize
    from vise.input_set.input_options import CategorizedInputOptions
    from vise.input_set.vasp_input_files import VaspInputFiles
    
    @explicit_serialize
    class WriteVaspInputsTask(FiretaskBase):
    
        required_params = ["task", "xc"]
        optional_params = ["input_options", "overridden_incar_settings"]
    
        def run_task(self, fw_spec):
            categorized_input_options = CategorizedInputOptions(
                structure=fw_spec["structure"],
                task=self["task"],
                xc=self["xc"],
                **self.get("input_options", {}))
    
            input_files = VaspInputFiles(categorized_input_options,
                                         self.get("overridden_incar_settings", {}))
            input_files.create_input_files(dirname=Path("."))
