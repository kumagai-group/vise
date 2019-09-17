from pathlib import Path

from pymatgen.io.vasp import VaspInput
from vise.input_set.incar import ViseIncar


class ViseVaspInput(VaspInput):

    @classmethod
    def from_directory(cls, input_dir, optional_files=None):
        """Construct ViseVspInput object from a directory.

        Note: Only Incar is overridden by ViseIncar.

        Args:
            input_dir (str):
            optional_files (dict):
                See docstrings of from_directory method of VaspInput
        """
        vasp_input = super().from_directory(input_dir, optional_files)
        path = Path(input_dir)
        vasp_input["INCAR"] = ViseIncar.from_file(path / "INCAR")

        return cls(vasp_input)