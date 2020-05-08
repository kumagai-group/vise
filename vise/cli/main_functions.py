# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

from argparse import Namespace
from pathlib import Path
from copy import deepcopy

from pymatgen.ext.matproj import MPRester
from vise.input_set.input_options import CategorizedInputOptions
from vise.input_set.vasp_input_files import VaspInputFiles
from vise.defaults import defaults
from vise.input_set.task import Task
from vise.input_set.xc import Xc
from vise.input_set.prior_info import PriorInfoFromCalcDir
from vise.input_set.kpoints_mode import KpointsMode


def get_poscar_from_mp(args: Namespace) -> None:
    s = MPRester().get_structure_by_material_id(args.mpid)
    s.to(fmt="poscar", filename=args.poscar)


def vasp_set(args: Namespace) -> None:

    option_kwargs = {}

    if args.prev_dir:
        pi = PriorInfoFromCalcDir(args.prev_dir, args.file_transfer_type)
        pi.file_transfers.transfer()
        option_kwargs.update(pi.input_options_kwargs)

    if args.options:
        option_kwargs.update(args.options)

    if args.uniform_kpt_mode:
        option_kwargs["kpt_mode"] = KpointsMode.uniform

    options = CategorizedInputOptions(
        structure=args.structure,
        task=args.task,
        xc=args.xc,
        kpt_density=args.kpt_density,
        override_potcar_set=args.override_potcar_set,
        charge=args.charge,
        **option_kwargs)

    overridden_incar_settings = deepcopy(defaults.user_incar_settings)
    if args.user_incar_settings:
        overridden_incar_settings.update(args.user_incar_settings)
    vif = VaspInputFiles(options, overridden_incar_settings)
    vif.create_input_files(Path.cwd())

