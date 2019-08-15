# -*- coding: utf-8 -*-

import sys

from atomate.utils.utils import get_database

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

dbfile = sys.argv[1]

db = get_database(dbfile)

materials = db["materials"]

mp_data = db["mp_e_above_hull"]

for m in materials.find({"thermo": {"$exists": True}},
                        projection=["thermo.energy_per_atom", "task_info"]):

    if m["task_info"][0]["swapped_element"] == {}:
        mp_id = m["task_info"][0]["mp_ids"]
        energy = m["thermo"]["energy_per_atom"]
        single_mp_data = mp_data.find_one({"material_id": mp_id},
                                          projection=["energy_per_atom", "pretty_formula"])
        energy_mp = single_mp_data["energy_per_atom"]
        energy_difference = energy - energy_mp
        pretty_formula = single_mp_data["pretty_formula"]

        print(mp_id, pretty_formula, energy, energy_mp, energy_difference)
