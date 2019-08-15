# -*- coding: utf-8 -*-
import argparse
import os
import warnings

from obadb.analyzer.vasp_process_analyzer import check_vasp_output, vasp_convergence_ionic, vasp_convergence_electronic
from obadb.analyzer.band_plotter import PrettyBSPlotter
from obadb.analyzer.dos_plotter import get_dos_plot
from obadb.util.structure_handler import SYMPREC

__version__ = "0.0.1"
__date__ = "8.Aug.2018"

# Following defaults determine the condition of automatic defect calculations.
# electronegativity difference for antisites and substitutional impurities
_EN_DIFF = 1.0
# Maximum displacement distance
_DISTANCE = 0.2
# Cutoff radius in which atoms are perturbed.
_CUTOFF = 3.0
_SYMPREC = SYMPREC
_ANGLE_TOLERANCE = -1.0


def main():
    parser = argparse.ArgumentParser(
        description="""                            
    pydefect is a package for first-principles point defect calculations. It 
    allows us to construct input files, parse first-principles calculation 
    results, and analyze data.""",
        epilog="""                                 
    Author: Yu Kumagai, Akira Takahashi
    Version: {}                                                                 
    Last updated: {}""".format(__version__, __date__),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #        allow_abbrev=False)

    subparsers = parser.add_subparsers()

    # -- check_convergence -----------------------------------------------------
    parser_check_convergence = subparsers.add_parser(
        name="check_convergence",
        description="Tools for checking the vasp calculation results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['sr'])

    parser_check_convergence.add_argument(
        "-c", "--convergence", dest="convergence", action="store_true",
        help="Check convergence of vasp calc at electronic and ionic steps.")
    parser_check_convergence.add_argument(
        "--dirs", dest="dirs", nargs="+", type=str,
        help="Directory names.")
    parser_check_convergence.set_defaults(func=check_convergence)

    # -- plot_band -----------------------------------------------------------
    parser_plot_band = subparsers.add_parser(
        name="plot_band",
        description="Tools for plotting band structures",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['pb'])
    parser_plot_band.add_argument(
        "-v", dest="vasprun", nargs="+", type=str)
    parser_plot_band.add_argument(
        "-v2", dest="vasprun2", nargs="+", type=str)
    parser_plot_band.add_argument(
        "-k", dest="kpoints", nargs="+", type=str)
    parser_plot_band.add_argument(
        "-y", dest="y_range", nargs="+", type=float)
    parser_plot_band.add_argument(
        "-f", dest="filename", type=str, help="pdf file name.")
    parser_plot_band.add_argument(
        "-a", dest="absolute", action="store_false",
        help="Show in the absolute energy scale.")
    parser_plot_band.add_argument(
        "-l", dest="legend", action="store_false",
        help="Not show the legend.")

    parser_plot_band.set_defaults(func=plot_band)

    # -- plot_dos -----------------------------------------------------------
    parser_plot_dos = subparsers.add_parser(
        name="plot_dos",
        description="Tools for plotting density of states",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['pd'])
    parser_plot_dos.add_argument(
        "-v", dest="vasprun", type=str)
    parser_plot_dos.add_argument(
        "-cv", dest="cbm_vbm", type=float, nargs="+",
        help="Set CBM and VBM.")
    parser_plot_dos.add_argument(
        "-t", dest="pdos_type", type=str, default="element",
        help=".")
    parser_plot_dos.add_argument(
        "-s", dest="specific", type=str, nargs="+", default=None, help=".")
    parser_plot_dos.add_argument(
        "-o", dest="orbital", action="store_false",
        help="Switch off the orbital decomposition.")
    parser_plot_dos.add_argument(
        "-x", dest="x_range", nargs="+", type=float, default=None,
        help="Set energy minimum and maximum.")
    parser_plot_dos.add_argument(
        "-y", dest="ymaxs", nargs="+", type=float, default=None,
        help="Set max values of y ranges. Support two ways." +
             "1st: total_max, each_atom" +
             "2nd: total_max, 1st_atom, 2nd_atom, ...")
    parser_plot_dos.add_argument(
        "-f", dest="filename", type=str, help="pdf file name.")
    parser_plot_dos.add_argument(
        "-a", dest="absolute", action="store_false",
        help="Show in the absolute energy scale.")
    parser_plot_dos.add_argument(
        "-l", dest="legend", action="store_false",
        help="Not show the legend.")
    parser_plot_dos.add_argument(
        "-cfv", dest="cfv", action="store_false",
        help="Not crop the first value.")
    parser_plot_dos.add_argument(
        "--symprec", dest="symprec", type=float, default=_SYMPREC,
        help="Set precision used for symmetry analyzer [A].")

    parser_plot_dos.set_defaults(func=plot_dos)

    args = parser.parse_args()
    args.func(args)


def check_convergence(args):
    for d in args.dirs:
        if os.path.isdir(d) is False:
            warnings.warn(message=d + " does not exist, so nothing is done.")
        else:
            check = check_vasp_output(d,
                                      contcar_name=args.poscar,
                                      outcar_name=args.outcar,
                                      vasprun_name=args.vasprun)
            if check["all"]:
                if vasp_convergence_ionic(d, args.vasprun):
                    ionic = "Y"
                else:
                    ionic = "N"
                if vasp_convergence_electronic(d, args.vasprun):
                    electronic = "Y"
                else:
                    electronic = "N"
                print("{:>20}  ionic:{:>3}  electronic:{:>3}".
                      format(d, ionic, electronic))
            else:
                contcar = "Y"
                outcar = "Y"
                vasprun = "Y"

                if check["contcar"] is None:
                    contcar = "NA"
                elif check["contcar"] is False:
                    contcar = "N"

                if check["outcar"] is None:
                    outcar = "NA"
                elif check["outcar"] is False:
                    outcar = "N"

                if check["vasprun"] is None:
                    vasprun = "NA"
                elif check["vasprun"] is False:
                    vasprun = "N"

            print("{:>20}  CONTCAR:{:>3}  OUTCAR:{:>3}, vasprun:{:>3}".
                  format(d, contcar, outcar, vasprun))


def plot_band(args):

    p = PrettyBSPlotter(args.kpoints, args.vasprun, args.vasprun2,
                        args.absolute, args.y_range, args.legend)

    if args.filename:
        p.save_fig(args.filename, format_type="pdf")
    else:
        p.show_fig()


def plot_dos(args):
    dos = get_dos_plot(vasprun_file=args.vasprun,
                       cbm_vbm=args.cbm_vbm,
                       pdos_type=args.pdos_type,
                       specific=args.specific,
                       orbital=args.orbital,
                       xlim=args.x_range,
                       ymaxs=args.ymaxs,
                       zero_at_efermi=args.absolute,
                       legend=args.legend,
                       crop_first_value=args.cfv,
                       symprec=args.symprec)

    if args.filename:
        dos.savefig(args.filename, format="pdf")
    else:
        dos.show()


if __name__ == "__main__":
    main()
