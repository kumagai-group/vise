# -*- coding: utf-8 -*-

from pymatgen.io.vasp import Poscar

from obadb.vasp.input_set import ObaSet
from obadb.custodian.oba_vaspjob import ObaVaspJob

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


if __name__ == "__main__":
    from custodian.custodian import Custodian
    from obadb.custodian.custodian_vasp_handler \
        import MeshSymmetryErrorHandler, NonConvergingErrorHandler, \
        PotimErrorHandler
    from obadb.custodian.oba_custodian_vasp_handler \
        import ObaUnconvergedErrorHandler, ObaVaspErrorHandler,\
        ObaMemorySwapHandler

    import argparse

    parser = argparse.ArgumentParser()

    mpi_exec = "/home/common/bin/openmpi-1.8.8_intel-16.0.2/bin/mpirun"
    vasp_exec = "/home/kuma/bin/vasp.5.4.4/bin/vasp_std"
    default_cmd = "{} {}".format(mpi_exec, vasp_exec)

    parser.add_argument(
        "-c", "--command", dest="command", nargs="?",
        default=default_cmd, type=str,
        help="VASP command. Defaults to pvasp. If you are using mpirun, "
             "set this to something like \"mpirun pvasp\".", )

    parser.add_argument(
        "-kpt", "--converge_kpt", dest="converge_kpt",
        action="store_true",
        help="Relax and converge KPOINTS")

    parser.add_argument(
        "-rw", "--remove_wavecar", dest="remove_wavecar",
        action="store_true",
        help="Remove finish/WAVECAR")

    parser.add_argument(
        "-mr", "--max_relax", dest="max_relax",
        default=10, type=int,
        help="Maximum number of relaxations to allow")

    parser.add_argument(
        "-kcr", "--kpoints_criteria", dest="kpoints_criteria",
        default=0.03, type=float,
        help="convergence criteria of kpoints"
    )

    args = parser.parse_args()

    # handlers
    handlers = [ObaMemorySwapHandler()]
    # handlers = [ObaVaspErrorHandler(), MeshSymmetryErrorHandler(),
    #             ObaUnconvergedErrorHandler(), NonConvergingErrorHandler(),
    #             PotimErrorHandler()]
    kpoints_criteria = args.kpoints_criteria
    cmd = args.command
    max_relax_number = args.max_relax
    removes_wavecar = args.remove_wavecar

    if args.converge_kpt:
        c = Custodian(handlers,
                      ObaVaspJob.kpt_converge(
                          cmd, max_relax_number, kpoints_criteria,
                      removes_wavecar=removes_wavecar),
                      polling_time_step=5, monitor_freq=1,
                      max_errors=10, gzipped_output=False)
    else:
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
        c = Custodian(handlers,
                      ObaVaspJob.get_runs_geom(
                          cmd, max_relax_number,
                          removes_wavecar=removes_wavecar,
                          removes_outputs_in_current_dir=True),
                      polling_time_step=5, monitor_freq=1,
                      max_errors=10, gzipped_output=False)
    c.run()
