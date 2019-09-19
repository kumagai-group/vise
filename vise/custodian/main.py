# -*- coding: utf-8 -*-
from pydash import max_
from vise.custodian.jobs import ViseVaspJob
from custodian.custodian import Custodian
#from vise.custodian.vise_custodian_vasp_handler import ViseMemorySwapHandler

import argparse


__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-c", "--command", dest="vasp_cmd", nargs="+", type=str,
        help="VASP command. If you are using mpirun, set this to something like \"mpirun pvasp\".", )

    parser.add_argument(
        "-kpt", "--converge_kpt", dest="converge_kpt", action="store_true",
        help="Relax and converge KPOINTS")

    parser.add_argument(
        "-rw", "--remove_wavecar", dest="rm_wavecar", action="store_true",
        help="Remove finish/WAVECAR")

    parser.add_argument(
        "--max_relax_num", dest="max_relax_num", default=10, type=int,
        help="Maximum number of relaxations.")

    parser.add_argument(
        "-kcr", "--kpoints_criteria", dest="kpoints_criteria",
        default=0.03, type=float,
        help="convergence criteria of kpoints"
    )

    args = parser.parse_args()

    if len(args.vasp_cmd) == 1:
        vasp_cmd = args.vasp_cmd[0].split()
    # handlers
    handlers = []
    # handlers = [ObaVaspErrorHandler(), MeshSymmetryErrorHandler(),
    #             ObaUnconvergedErrorHandler(), NonConvergingErrorHandler(),
    #             PotimErrorHandler()]
#    kpoints_criteria = args.kpoints_criteria

    # if args.converge_kpt:
    #     c = Custodian(handlers,
    #                   ViseVaspJob.kpt_converge(
    #                       cmd, max_relax_number, kpoints_criteria,
    #                   removes_wavecar=removes_wavecar),
    #                   polling_time_step=5, monitor_freq=1,
    #                   max_errors=10, gzipped_output=False)
    # else:
        # st = Poscar.from_file("POSCAR").structure
        # oba_vis = ObaSet.make_input(st)
        # # oba_vis = \
        # #     ObaSet.from_prev_calc(".",
        # #                           parse_calc_results=False,
        # #                           standardize_structure=True,
        # #                           parse_potcar=True,
        # #                           parse_incar=True,
        # #                           parse_kpoints=True)
        # oba_vis.write_input(".")
    # job = ViseVaspJob.structure_optimization_run(vasp_cmd=args.vasp_cmd,
    #                                              max_relax_num=args.max_relax_num,
    #                                              removes_wavecar=args.rm_wavecar)
#    job =
    job = ViseVaspJob.kpt_converge(vasp_cmd=args.vasp_cmd)
    c = Custodian(handlers=handlers,
                  jobs=job,
                  polling_time_step=5,
                  monitor_freq=1,
                  max_errors=10,
                  gzipped_output=False)
    c.run()
