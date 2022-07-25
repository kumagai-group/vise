# -*- coding: utf-8 -*-
#  Copyright (c) 2020. Distributed under the terms of the MIT License.
import os

import pytest
from importlib import reload


@pytest.fixture
def modified_defaults(tmpdir):
    os.chdir(tmpdir)
    tmpdir.join("vise.yaml").write("""
kpoint_density: 0.5
overridden_potcar:
   O_h
dos_step_size: 0.01
band_mesh_distance: 0.25
str_opt_encut_factor: 1.5
ldauu:
   Mn: 10
ldaul:
   Mn:  4
potcar_set: mp   
""")

    # Need to import here when run multiple tests together.
    import vise.defaults
    reload(vise.defaults)
    return vise.defaults.defaults


@pytest.fixture
def simple_modified_defaults(tmpdir):
    os.chdir(tmpdir)
    tmpdir.join("vise.yaml").write("""
kpoint_density: 0.5
""")

    # Need to import here when run multiple tests together.
    import vise.defaults
    reload(vise.defaults)
    return vise.defaults.defaults


def test_user_settings_in_defaults(modified_defaults):
    assert modified_defaults.user_settings["ldauu"] == {"Mn": 10}
    assert modified_defaults.ldauu == {"Mn": 10}


@pytest.mark.skip(reason='default_override')
def test_integration(tmpdir):
    tmpdir.chdir()
    os.chdir(tmpdir)
    tmpdir.join("vise.yaml").write("""
kpoint_density: 0.5
overridden_potcar:
   O_h
dos_step_size: 0.01
band_mesh_distance: 0.25
str_opt_encut_factor: 1.5
ldauu:
   Mn: 10
ldaul:
   Mn:  4
potcar_set: mp
user_incar_settings:
    ISPIN: 2  
options:
    cutoff_energy: 1000
""")

    from vise.cli.main import parse_args
    tmpdir.join("POSCAR").write("""Mn1 O1
1.0
0.000000 2.123447 2.123447
2.123447 0.000000 2.123447
2.123447 2.123447 0.000000
Mn O
1 1
direct
0.000000 0.000000 0.000000 Mn
0.500000 0.500000 0.500000 O
""")
    args = parse_args(["vs", "-t", "band"])
    args.func(args)
    incar = tmpdir.join("INCAR")
    incar_expected = """# algorithm
ALGO  =  Normal

# accuracy
PREC   =  Normal
LREAL  =  False
EDIFF  =  1e-05
ENCUT  =  1000
LASPH  =  True
NELM   =  100

# ionic relaxation
ISIF    =  0
IBRION  =  2
NSW     =  0

# occupation
ISMEAR  =  0
SIGMA   =  0.1

# spin
ISPIN  =  2

# IO control
LWAVE   =  False
LCHARG  =  False

# analyzer
NBANDS  =  23
LORBIT  =  10

# hubbard u
LDAU       =  True
LDAUTYPE   =  2
LMAXMIX    =  4
LDAUPRINT  =  1
LDAUL      =  4 -1
LDAUU      =  10 0

# parallel
KPAR  =  2"""
    assert incar.read() == incar_expected

    potcar = tmpdir.join("POTCAR")
    potcar_expected = """  PAW_PBE Mn_pv 02Aug2007                
   13.0000000000000     
 parameters from PSCTR are:
   VRHFIN =Mn: 3p4s3d
   TITEL  = PAW_PBE Mn_pv 02Aug2007
   POMASS =   54.938; ZVAL   =   13.000    mass and valenz
   ENMAX  =  269.864; ENMIN  =  202.398 eV

   END of PSCTR-controll parameters
  PAW_PBE O_h 06Feb2004                  
   6.00000000000000     
 parameters from PSCTR are:
   VRHFIN =O: s2p4
   TITEL  = PAW_PBE O_h 06Feb2004
   POMASS =   16.000; ZVAL   =    6.000    mass and valenz
   ENMAX  =  700.000; ENMIN  =  500.000 eV

   END of PSCTR-controll parameters
"""

    assert potcar.read() == potcar_expected

    kpoints = tmpdir.join("KPOINTS")
    kpoints_expected = """k-path added by seekpath. Formula: MnO SG: 225 
23
Reciprocal
 0.00000000000000000  0.00000000000000000  0.00000000000000000   1 None
 0.50000000000000000  0.00000000000000000  0.00000000000000000   4 None
 0.50000000000000000  0.50000000000000000  0.00000000000000000   3 None
 0.00000000000000000  0.00000000000000000  0.00000000000000000   0 GAMMA
 0.12500000000000000  0.00000000000000000  0.12500000000000000   0 
 0.25000000000000000  0.00000000000000000  0.25000000000000000   0 
 0.37500000000000000  0.00000000000000000  0.37500000000000000   0 
 0.50000000000000000  0.00000000000000000  0.50000000000000000   0 X
 0.62500000000000000  0.25000000000000000  0.62500000000000000   0 U
 0.37500000000000000  0.37500000000000000  0.75000000000000000   0 K
 0.29999999999999999  0.29999999999999999  0.59999999999999998   0 
 0.22500000000000001  0.22500000000000001  0.45000000000000001   0 
 0.14999999999999999  0.14999999999999999  0.29999999999999999   0 
 0.07500000000000001  0.07500000000000001  0.15000000000000002   0 
 0.00000000000000000  0.00000000000000000  0.00000000000000000   0 GAMMA
 0.12500000000000000  0.12500000000000000  0.12500000000000000   0 
 0.25000000000000000  0.25000000000000000  0.25000000000000000   0 
 0.37500000000000000  0.37500000000000000  0.37500000000000000   0 
 0.50000000000000000  0.50000000000000000  0.50000000000000000   0 L
 0.50000000000000000  0.41666666666666669  0.58333333333333337   0 
 0.50000000000000000  0.33333333333333337  0.66666666666666663   0 
 0.50000000000000000  0.25000000000000000  0.75000000000000000   0 W
 0.50000000000000000  0.00000000000000000  0.50000000000000000   0 X
"""

    assert kpoints.read() == kpoints_expected

#
# def test_to_json_file(tmpdir, simple_modified_defaults):
#     os.chdir(tmpdir)
#     simple_modified_defaults.to_json_file()
#     kpoints = tmpdir.join("vise_defaults.json")
#
#     kpoints_expected