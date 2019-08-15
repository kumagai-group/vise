# -*- coding: utf-8 -*-

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

import sys
from atomate.utils.utils import get_database

db_file = sys.argv[1]
db = get_database(config_file=db_file)

#for d in db:
for i in db["materials"].find({'dielectric.epsilon_ionic.0.0': {"$exists": 1}}):
    dielectric_constant_xx = i["dielectric"]["epsilon_ionic"][0][0]
    dielectric_constant_yy = i["dielectric"]["epsilon_ionic"][1][1]
    dielectric_constant_zz = i["dielectric"]["epsilon_ionic"][2][2]
    dielectric = [dielectric_constant_xx, dielectric_constant_yy, dielectric_constant_zz]
    dielectric_constant = sum(dielectric) / 3
    max_dielectric = max(dielectric)

    if max_dielectric > 50:
        m_id = i["material_id"]
        formula = i["formula_pretty"]
        sg = i["sg_number"]
        mpid = i["mp_ids"][0]
        print("{0:>10s} {1:>17s} {2:5d} {3:8.1f} {4:8.1f} {5:8.1f} {6:8.1f}".format(m_id, formula, sg, dielectric_constant, dielectric_constant_xx, dielectric_constant_yy, dielectric_constant_zz))
        print(mpid)

