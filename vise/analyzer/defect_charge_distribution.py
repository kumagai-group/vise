# -*- coding: utf-8 -*-

import argparse
import itertools

import matplotlib.pyplot as plt
import numpy as np
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp.outputs import Chgcar, VolumetricData

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"


class ViseChgcar(Chgcar):

    @classmethod
    def from_file(cls, filename):
        (poscar, data, data_aug) = VolumetricData.parse_file(filename)
        return cls(poscar, data, data_aug=data_aug)

    def max_value_positions(self, data_type):

        supported_type = ["total", "diff", "up", "down"]
        if data_type in ["total", "diff"]:
            vals = self.data[data_type]
        elif data_type is "up":
            vals = self.spin_data[Spin.up]
        elif data_type is "down":
            vals = self.spin_data[Spin.down]
        else:
            raise ValueError(
                "Supported data_type is {}.".format(str(supported_type)))

        indices = tuple(zip(np.where(vals == vals.max())[0],
                            np.where(vals == vals.max())[1],
                            np.where(vals == vals.max())[2]))
        coords = []
        for i in indices:
            coords.append([i[0] / self.dim[0],
                           i[1] / self.dim[1],
                           i[2] / self.dim[2]])
        return indices, coords

    def radial_distribution(self, data_type, ind, radius, nbins=1):
        """
        Get integrated difference of atom index up to radius. This can be
        an extremely computationally intensive process, depending on how many
        grid points are in the VolumetricData.

        Args:
            data_type (str):  Type name
            ind (int): Index of atom.
            radius (float): Radius of integration.
            nbins (int): Number of bins. Defaults to 1. This allows one to
                obtain the charge integration up to a list of the cumulative
                charge integration values for radii for [radius/nbins,
                2 * radius/nbins, ....].

        Returns:
            Differential integrated charge as a np array of [[radius, value],
            ...]. Format is for ease of plotting. E.g., plt.plot(data[:,0],
            data[:,1])
        """

        # For non-spin-polarized runs, this is zero by definition.
        if not self.is_spin_polarized:
            radii = [radius / nbins * (i + 1) for i in range(nbins)]
            data = np.zeros((nbins, 2))
            data[:, 0] = radii
            return data

        struct = self.structure
        a = self.dim
        if ind not in self._distance_matrix or \
                self._distance_matrix[ind]["max_radius"] < radius:
            coords = []
            for (x, y, z) in itertools.product(*[list(range(i)) for i in a]):
                coords.append([x / a[0], y / a[1], z / a[2]])
            sites_dist = struct.lattice.get_points_in_sphere(
                coords, struct[ind].coords, radius)
            self._distance_matrix[ind] = {"max_radius": radius,
                                          "data": np.array(sites_dist)}

        data = self._distance_matrix[ind]["data"]

        # Use boolean indexing to find all charges within the desired distance.
        # data[:, 0]: shifted_coords
        # data[:, 1]: distances
        # data[:, 2]: indices
        # data[:, 3]: images

        inds = data[:, 1] <= radius
        dists = data[inds, 1]
        data_inds = np.rint(np.mod(list(data[inds, 0]), 1) *
                            np.tile(a, (len(dists), 1))).astype(int)

        supported_type = ["total", "diff", "up", "down"]
        if data_type in ["total", "diff"]:
            vals = [self.data[data_type][x, y, z] for x, y, z in data_inds]
        elif data_type is "up":
            vals = [self.spin_data[Spin.up][x, y, z] for x, y, z in data_inds]
        elif data_type is "down":
            vals = [self.spin_data[Spin.down][x, y, z] for x, y, z in data_inds]
        else:
            raise ValueError(
                "Supported data_type is {}.".format(str(supported_type)))

        hist, edges = np.histogram(dists, bins=nbins, range=[0, radius],
                                   weights=vals)
        hist_numbers, _ = np.histogram(dists, bins=nbins, range=[0, radius])

        self.hist_data = np.zeros((nbins, 2))
        self.hist_data[:, 0] = [sum(edges[i:i + 2]) / 2 for i in range(nbins)]
        # 4pi * r^2 * rho
        self.hist_data[:, 1] = \
            4 * np.pi * self.hist_data[:, 0] ** 2 * hist / hist_numbers
        volume = self.structure.lattice.volume
        mesh_distance = edges[1] - edges[0]
        self.hist_data[:, 1] = self.hist_data[:, 1] / volume

        for i in range(self.hist_data[:, 1].size):
            if sum(self.hist_data[:i + 1, 1]) * mesh_distance > 0.5:
                # Obtain from the calculation of the area of a trapezoid
                x0 = self.hist_data[i - 1, 0]
                y0 = self.hist_data[i, 1]
                y1 = self.hist_data[i + 1, 1]
                x = (0.5 - sum(self.hist_data[:i, 1]) * mesh_distance) \
                    * 2 / (y0 + y1) + x0
                self.half_point = x
                break

        print("Central position")
        print(struct[ind].frac_coords)
        print("half_point")
        print(self.half_point)
        print("sum")
        print(sum(self.hist_data[:, 1]) * mesh_distance)

    def plot_radial_distribution(self, filename=None):
        plt.plot(self.hist_data[:, 0], self.hist_data[:, 1], 'o-')
        plt.axvline(self.half_point)

        if filename:
            plt.savefig(filename)
        else:
            plt.show()


def main():
    parser = argparse.ArgumentParser(
        description=""".""",
        epilog="""                                 
    Author: Yu Kumagai
    Version: {}                                                                 
    Last updated: {}""".format(__version__, __date__),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "-p", "--parchg", dest="parchg", default="PARCHG", type=str,
        help="PARCHG name to be parsed.")
    parser.add_argument(
        "-c", "--center", dest="center", type=float, nargs="+", default=None,
        help="Atom index beginning from 1 or fractional coordinates.")
    parser.add_argument(
        "-t", "--type", dest="type", type=str, default="up",
        help="Type from total, diff, up, down.")
    parser.add_argument(
        "-d", "--distance", dest="distance", type=float, default=4,
        help="Distance.")
    parser.add_argument(
        "-b", "--nbins", dest="nbins", type=int, default=40,
        help="Number of bins.")
    parser.add_argument(
        "-f", dest="filename", type=str, default=None, help="File name.")

    args = parser.parse_args()

    c = ViseChgcar.from_file(args.parchg)
    if args.center is not None:
        if len(args.center) == 3:
            c.structure.append("X", args.center)
            index = len(c.structure) - 1
        elif len(args.center) == 1:
            index = int(args.center) - 1
        else:
            assert ValueError("The input centering is not proper.")
    else:
        centers = c.max_value_positions(data_type=args.type)[1]
        if len(centers) > 1:
            assert ValueError("The number of positions with the largest charge"
                              "is more than 1.")
        center = centers[0]
        c.structure.append("X", center)
        index = len(c.structure) - 1

    c.radial_distribution(args.type, index, args.distance, nbins=args.nbins)
    c.plot_radial_distribution(filename=args.filename)


if __name__ == "__main__":
    main()
