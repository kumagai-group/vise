#!/usr/bin/env python

import argparse

from pymatgen.io.vasp.outputs import BSVasprun


def band_gap_properties(v, round_digit_number=1):
    """
    Args:
        v (Vasprun):
        round_digit_number (int):
    """

    for s in v.eigenvalues:
        occupation = [round(sum([i[1] for i in k]), round_digit_number)
                      for k in v.eigenvalues[s]]

        if len(set(occupation)) != 1:
            return None

    band_structure = v.get_band_structure()
    vbm = band_structure.get_vbm()
    cbm = band_structure.get_cbm()

    band_gap = band_structure.get_band_gap()

    if band_gap["energy"] == 0.0:
        return None

    vbm_info = {'energy': vbm['energy'],
                'band_index': dict(vbm['band_index']),
                'kpoints': [v.actual_kpoints[x] for x in vbm["kpoint_index"]]}
    cbm_info = {'energy': cbm['energy'],
                'band_index': dict(cbm['band_index']),
                'kpoints': [v.actual_kpoints[x] for x in cbm["kpoint_index"]]}

    return band_gap, vbm_info, cbm_info


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", dest="vasprun", type=str,
                         default="vasprun.xml", metavar="FILE")
    args = parser.parse_args()

    v = BSVasprun(args.vasprun)
    bgp = band_gap_properties(v)
    if bgp is None:
        print("Metallic system")
    else:
        band_gap, vbm_info, cbm_info = bgp
        print("CBM info {}".format(cbm_info))
        print("VBM info {}".format(vbm_info))
        print("band gap info {}".format(band_gap))


if __name__ == "__main__":
    main()
