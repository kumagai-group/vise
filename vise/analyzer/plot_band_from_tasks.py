# -*- coding: utf-8 -*-

import os
import sys

from atomate.utils.utils import get_database
from obadb.analyzer.band_plotter import PrettyBSPlotter

__author__ = "Yu Kumagai"
__maintainer__ = "Yu Kumagai"

db = get_database(config_file=sys.argv[1])
db_tasks = db["tasks"]

for t in db_tasks.find({"task_label": "band_hse"}):
    directory = t["dir_name"].split(":")[1]

    os.chdir(directory)

    kpoints = "KPOINTS"
    vasprun = "vasprun.xml"

    mp_id = t["db_ids"]["mp_id"]
#    directory2 = db_tasks.find_one({"db_ids.mp_id": mp_id})["dir_name"].split(":")[1]
#    vasprun2 = os.path.join(directory2, "vasprun.xml")

    print(t["formula_pretty"])
    print(directory)
#    p = PrettyBSPlotter(kpoints, vasprun, vasprun2)
    p = PrettyBSPlotter(kpoints, vasprun)
    p.save_fig("band.pdf", format_type="pdf")
