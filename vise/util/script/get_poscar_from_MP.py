# -*- coding: utf-8 -*-

import argparse

from pymatgen.ext.matproj import MPRester

__author__ = "Yu Kumagai"
__copyright__ = "Copyright 2018, Oba group"
__version__ = "0.1"
__maintainer__ = "Yu Kumagai"
__email__ = "yuuukuma@gmail.com"
__status__ = "Development"
__date__ = "Aug 3, 2018"


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--number", dest="number",
                        type=int, help="MP entry number w/o mp-")
    parser.add_argument("-p", "--poscar", dest="poscar", default="POSCAR",
                        type=str,
                        help="POSCAR-type file name for the unitcell.",
                        metavar="FILE")

    opts = parser.parse_args()

    s = MPRester().get_structure_by_material_id("mp-" + str(opts.number))
    s.to(filename=opts.poscar)


if __name__ == "__main__":
    main()
