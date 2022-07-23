# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import numpy as np

# unit conversion and constants used in vasp.5.4.4 src
RYTOEV = 13.605826  # Rydberg to eV
AUTOA = 0.529177249  # atomic unit to angstrom
PI = 3.141592653589793238


def vasp_grid(encut: float, lattice_length: float, symprec: str) -> int:
    """ Calculate the numbers of FFT grids used in VASP

    << src from main.F of vasp.5.4.4>>
    XCUTOF =SQRT(INFO%ENMAX /RYTOEV)/(2*PI/(LATT_CUR%ANORM(1)/AUTOA))
    IF (INFO%SZPREC(1:1)=='h' .OR. INFO%SZPREC(1:1)=='a') THEN
      WFACT=4
    ELSE
      WFACT=3
    ENDIF
    GRID%NGPTAR(1)=XCUTOF*WFACT+0.5_q
    CALL FFTCHK(GRID%NGPTAR)
    <<end src>>
    FFTCHK update original grids based on specific algorithms,
     see setup_actual_grids() in this module file for details.
    """
    # In case of ACCURATE or HIGH, factor = 4
    factor = 4 if symprec[0].lower() in ["a", "h"] else 3
    cutoff_in_rydberg = encut / RYTOEV
    lattice_length_in_angstrom = lattice_length / AUTOA
    rec_lat_len = 2 * PI / lattice_length_in_angstrom
    return int(round(factor * np.sqrt(cutoff_in_rydberg) / rec_lat_len))


