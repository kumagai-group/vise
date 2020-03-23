#!/usr/bin/env python
#  Copyright (c) 2020. Distributed under the terms of the MIT License.

# import sys
# import argparse
# from pylab import *

# parser= argparse.ArgumentParser()
# parser.add_argument("-o", "--outcar", dest="outcar", default="OUTCAR",
#                     type=str, help="OUTCAR name.")
# parser.add_argument("-p", "--plot", dest="plot",action="store_true",
#                     help="Show plot.")
# opts = parser.parse_args()

from pymatgen.io.vasp import Outcar

#outcar = open(opts.outcar).readlines()

# nr_ion = None

# for line in outcar:
#     if line.find('NION') != -1:
#         nr_ion = int(line.split()[11])
#     elif line.find('EDIFFG') != -1:
#         try:
#             force_convergence = float(line.split()[2])
#         except:
#             force_convergence = "NA"
#         break

# if nr_ion is None:
#     raise ValueError("NION could not be found.")

# nr_line_force = []

# for nr_line, line in enumerate(outcar):
#     if line.find('TOTAL-FORCE') != -1:
#         nr_line_force.append(nr_line)

# if not nr_line_force:
#     raise ValueError("TOTAL-FORCE line not in {}. "
#                      "Stop here. Bye.".format(opts.outcar))

# force_raws = [[outcar[i + j + 2].split()[3:6]
#                for j in range(nr_ion)] for i in nr_line_force]
# forces = [[sum([float(z) ** 2 for z in y]) ** 0.5 for y in x]
#           for x in force_raws]
# max_average_forces = [[max(x), sum(x) / nr_ion] for x in forces]

# print("EDIFF = {}".format(force_convergence))
# print("maximum force | average force")

# for i in max_average_forces:
#     print("{0:12.7f} {1:12.7f}".format(i[0], i[1]))
# #    print("%12.7f  %12.7f".format(i[0], i[1]))

# if opts.plot:
#     title('Force convergence')
#     xlabel('Nr of ionic step')
#     ylabel('Residual force')
#     ediff = [-force_convergence for i in range(len(max_average_forces))]
#     plot(ediff)
#     plot(max_average_forces)
#     plot(max_average_forces,'r*')
#     show()

# sys.exit(0)