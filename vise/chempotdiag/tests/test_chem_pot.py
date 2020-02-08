#  Copyright (c) Oba-group 
#  Distributed under the terms of the MIT License.

import unittest
import os

from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import Element

from chempotdiag.compound import Compound, CompoundsList
from chempotdiag.vertex import Vertex
from chempotdiag.chem_pot_diag import ChemPotDiag
from chempotdiag.config import (
    MOLECULE_SUFFIX, ROOM_TEMPERATURE, REFERENCE_PRESSURE)
from chempotdiag.gas import Gas

DISPLAY_DIAGRAM = False  # If True, plot diagrams. For check figure by eyes.

EXAMPLE_DIR = os.path.dirname(os.path.abspath(__file__)) + "/"

FILENAME_1D = EXAMPLE_DIR + "energy_1d.txt"
FILENAME_2D = EXAMPLE_DIR + "energy_2d.txt"
FILENAME_3D = EXAMPLE_DIR + "energy_MP-Ca-Al-O.txt"
FILENAME_4D = EXAMPLE_DIR + "energy_4d.txt"
# For read DFT test. We don't check these files are physically proper.
DFT_DIRECTORIES = [f"{EXAMPLE_DIR}/dft_data/O2{MOLECULE_SUFFIX}/",
                   f"{EXAMPLE_DIR}/dft_data/Mg/",
                   f"{EXAMPLE_DIR}/dft_data/MgO/"]
O2_OUTCAR_VAL = -10.26404723 / 2
MGO_OUTCAR_VAL = -12.52227962 / 2
MG_OUTCAR_VAL = -3.42786289 / 2
POSCAR_NAME = "POSCAR-finish"
OUTCAR_NAME = "OUTCAR-finish"
VASPRUN_NAME = "vasprun-finish.xml"


class TestChemPot(unittest.TestCase):

    #  Simple test to generate class object from text file without error.
    def test_read_from_file(self):
        _ = ChemPotDiag.from_file(FILENAME_1D)
        _ = ChemPotDiag.from_file(FILENAME_2D)
        _ = ChemPotDiag.from_file(FILENAME_3D)
        _ = ChemPotDiag.from_file(FILENAME_4D)

    def test_calculation(self):
        cl = CompoundsList([
            Compound(None, Composition("Al"), 0),
            Compound(None, Composition("Mg"), 0),
            Compound(None, Composition("MgAl"), -10),
            Compound(None, Composition("Mg2Al3"), 100)
        ])
        cpd = ChemPotDiag.from_compound_list(cl, temperature=None, pressure=None)

        # includes drawing boundary, -11 = -10 (minimum energy) * 1.1
        actual = cpd.vertices
        expected = [
            Vertex(None, {"Mg": -11, "Al": -11}),
            Vertex(None, {"Mg": 0, "Al": -11}),
            Vertex(None, {"Mg": -11, "Al": 0}),
            Vertex(None, {"Mg": -10, "Al": 0}),
            Vertex(None, {"Mg": 0, "Al": -10}),
        ]
        for v in actual:
            self.assertTrue(any(v.almost_equal(e) for e in expected))

        stable_comp: CompoundsList = cpd.stable_compounds
        expected = CompoundsList([
            Compound(None, Composition("Al"), 0),
            Compound(None, Composition("Mg"), 0),
            Compound(None, Composition("MgAl"), -10),
        ])
        self.assertTrue(expected.almost_equal(stable_comp))

        # unstable
        unstable: Compound = cpd.unstable_compounds[0]
        self.assertEqual(unstable.composition, Composition("Mg0.4Al0.6"))
        self.assertAlmostEqual(unstable.energy, 100/5, 5)

    def test_read_from_vasp_files(self):
        poscar_paths = [d+POSCAR_NAME for d in DFT_DIRECTORIES]
        outcar_paths = [d+OUTCAR_NAME for d in DFT_DIRECTORIES]
        vasprun_paths = [d+VASPRUN_NAME for d in DFT_DIRECTORIES]
        for temp, pres in [
            # only consider zero point energy
            (0, None),
            # default: REFERENCE_PRESSURE
            (ROOM_TEMPERATURE, None),
            # Specify both T and P
            (ROOM_TEMPERATURE, 1e+6)]:
            #  from outcar
            p = pres if pres else REFERENCE_PRESSURE
            cpd = \
                ChemPotDiag.from_vasp_calculations_files(
                    poscar_paths, outcar_paths,
                    temperature=temp, pressure={"O2": p})
            # element energy
            elem_en = cpd.element_energy
            if temp:
                self.assertAlmostEqual(elem_en[Element("O")],
                                       O2_OUTCAR_VAL +
                                       Gas.O2.energy_shift(temp, p))
            else:
                self.assertAlmostEqual(elem_en[Element("O")],
                                       O2_OUTCAR_VAL +
                                       Gas.O2.energy_shift(temp, p))
        _ = ChemPotDiag.from_vasp_calculations_files(poscar_paths,
                                                     vasprun_paths,
                                                     fmt="vasprun",
                                                     temperature=273,
                                                     pressure={"O2": 1e+6})
        # When calculation of any simple substance is not found,
        # program should be raise error.
        self.assertRaises(ValueError,
                          lambda: ChemPotDiag.from_vasp_calculations_files(
                              poscar_paths[1:], outcar_paths[1:]))

    #  Test to draw 1d diagram
    # @unittest.skipIf(not DISPLAY_DIAGRAM, "not display chempotdiag")
    def test_draw_diagram_1d(self):
        cp = ChemPotDiag.from_file(FILENAME_1D)
        cp.draw_diagram()

    #  Test to draw 2d diagram
    @unittest.skipIf(not DISPLAY_DIAGRAM, "not display chempotdiag")
    def test_draw_diagram_2d(self):
        cp = ChemPotDiag.from_file(FILENAME_2D)
        cp.draw_diagram()

    @unittest.skipIf(not DISPLAY_DIAGRAM, "not display chempotdiag")
    def test_draw_diagram_2d_with_vertex_name(self):
        cp = ChemPotDiag.from_file(FILENAME_2D)
        rc = "MgO"
        title_str = "With_vertex_name test, remarked_compound = {0}".format(rc)
        cp.draw_diagram(title=title_str, remarked_compound=rc)

        #  Maybe label of previous plot remains if you forget to delete.
        rc = "Mg2O"
        title_str = "With_vertex_name test, remarked_compound = {0}".format(rc)
        cp.draw_diagram(title=title_str, remarked_compound=rc)

    @unittest.skipIf(not DISPLAY_DIAGRAM, "not display chempotdiag")
    def test_draw_diagram_2d_without_label(self):
        cp = ChemPotDiag.from_file(FILENAME_2D)
        cp.draw_diagram(title="diagram_2d_without_label", with_label=False)

    @unittest.skipIf(not DISPLAY_DIAGRAM, "not display chempotdiag")
    def test_draw_range_diagram2d(self):
        cp = ChemPotDiag.from_file(FILENAME_2D)
        rc = "MgO"
        draw_range = -20
        title_str = "{0}, draw_range = {1}".format(rc, draw_range)
        cp.draw_diagram(title=title_str,
                        remarked_compound=rc,
                        draw_range=draw_range)

    @unittest.skipIf(not DISPLAY_DIAGRAM, "not display chempotdiag")
    def test_wrong_draw_range_diagram2d(self):
        cp = ChemPotDiag.from_file(FILENAME_2D)
        draw_range = 100
        rc = "MgO"
        title_str = "{0}, draw_range = {1}".format(rc, draw_range)
        self.assertRaises(ValueError,
                          lambda: cp.draw_diagram(title=title_str,
                                                  remarked_compound=rc,
                                                  draw_range=draw_range))

    #  Test to draw 3d diagram
    @unittest.skipIf(not DISPLAY_DIAGRAM, "not display chempotdiag")
    def test_draw_diagram_3d(self):
        cp = ChemPotDiag.from_file(FILENAME_3D)
        cp.draw_diagram()

    def test_draw_diagram_3d_with_vertex_name(self):
        cp = ChemPotDiag.from_file(FILENAME_3D)
        rc = "Ca"
        title_str = "With_vertex_name test, remarked_compound = {0}".format(rc)
        if DISPLAY_DIAGRAM:
            cp.draw_diagram(title=title_str, remarked_compound=rc)  # For debug
        #  Maybe label of previous plot remains if you forget to delete.
        rc = "Ca11Al14O32"
        title_str = "With_vertex_name test, remarked_compound = {0}".format(rc)
        if DISPLAY_DIAGRAM:
            cp.draw_diagram(title=title_str,
                            remarked_compound=rc,
                            elements=["Al", "Ca", "O"])
        got_vl = cp.get_neighbor_vertices(rc)
        self.assertTrue(
            got_vl[0].almost_equal(
                Vertex(None, {"Al": -8.72916, "Ca": -6.84027, "O": 0})))
        self.assertTrue(
            got_vl[1].almost_equal(
                Vertex(None, {"Al": -8.88137, "Ca": -6.64655, "O": 0})))
        self.assertTrue(
            got_vl[2].almost_equal(
                Vertex(None, {"Al": -0.70214, "Ca": -1.19373, "O": -5.45282})))
        self.assertTrue(
            got_vl[3].almost_equal(
                Vertex(None, {"Al": -0.25935, "Ca": -1.19373, "O": -5.64654})))

    @unittest.skipIf(not DISPLAY_DIAGRAM, "not display chempotdiag")
    def test_draw_diagram_3d_without_label(self):
        cp = ChemPotDiag.from_file(FILENAME_3D)
        cp.draw_diagram(title="diagram_3d_without_label",
                        with_label=False)

    @unittest.skipIf(not DISPLAY_DIAGRAM, "not display chempotdiag")
    def test_draw_range_diagram3d(self):
        cp = ChemPotDiag.from_file(FILENAME_3D)
        rc = "Ca11Al14O32"
        draw_range = -200
        title_str = f"{rc}, draw_range = {draw_range}"
        cp.draw_diagram(title=title_str,
                        remarked_compound=rc,
                        draw_range=draw_range)

    @unittest.skipIf(not DISPLAY_DIAGRAM, "not display chempotdiag")
    def test_wrong_draw_range_diagram3d(self):
        cp = ChemPotDiag.from_file(FILENAME_3D)
        draw_range = 100
        rc = "Ca11Al14O32"
        title_str = f"{rc}, draw_range = {draw_range}"
        self.assertRaises(ValueError,
                          lambda: cp.draw_diagram(title=title_str,
                                                  remarked_compound=rc,
                                                  draw_range=draw_range))

    def test_dim(self):
        cp1 = ChemPotDiag.from_file(FILENAME_1D)
        self.assertEqual(cp1.dim, 1)
        cp2 = ChemPotDiag.from_file(FILENAME_2D)
        self.assertEqual(cp2.dim, 2)
        cp3 = ChemPotDiag.from_file(FILENAME_3D)
        self.assertEqual(cp3.dim, 3)
        cp4 = ChemPotDiag.from_file(FILENAME_4D)
        self.assertEqual(cp4.dim, 4)

    @unittest.skipIf(not DISPLAY_DIAGRAM, "not display chempotdiag")
    def test_draw_diagram_3d_all_examples(self):
        for s1 in ("MP", "Oba"):
            for s2 in ("Ca-Al-O", "Mg-Ca-O", "Sr-Bi-N", "Sr-Fe-O", "Sr-Ti-O"):
                path = EXAMPLE_DIR + "energy_" + s1 + "-" + s2 + ".txt"
                cp = ChemPotDiag.from_file(path)
                rc = cp.stable_compounds[1].name
                cp.draw_diagram(remarked_compound=rc)

    @unittest.skipIf(not DISPLAY_DIAGRAM, "not display chempotdiag")
    def test_draw_diagram_2d_with_temperature_pressure(self):
        poscar_paths = [d+POSCAR_NAME for d in DFT_DIRECTORIES]
        outcar_paths = [d+OUTCAR_NAME for d in DFT_DIRECTORIES]
        for t in [298.15, 500, 1000]:
            for o2_p in [1e-5, 1e+100]:
                title_str = "T = {0} (K), P_O2 = {1} (Pa)".format(t, o2_p)
                p = {"O2": o2_p}
                cp = ChemPotDiag.from_vasp_calculations_files(poscar_paths,
                                                              outcar_paths,
                                                              temperature=t,
                                                              pressure=p)
                rc = cp.stable_compounds[1].name
                cp.draw_diagram(
                    title=title_str,
                    remarked_compound=rc)

        title_str = "temperature = 298.15, pressure = None"
        cp = ChemPotDiag.from_vasp_calculations_files(poscar_paths,
                                                      outcar_paths,
                                                      temperature=298.15,
                                                      pressure=None)
        rc = cp.stable_compounds[1].name
        cp.draw_diagram(
            title=title_str,
            remarked_compound=rc)

    def test_neighbor_vertices_as_dict(self):
        """
        Reference: result of commit ab86c47fc6957ffabcd942c5afe2dd8ee91ef82b
        """
        cp = ChemPotDiag.from_file(FILENAME_3D)  # Ca,Al,O

        def almost_equal_dict(o1, o2):
            self.assertTrue(type(o1), type(o2))
            if isinstance(o1, float):
                self.assertAlmostEqual(o1, o2, 5)
            if isinstance(o1, dict):
                self.assertSetEqual(set(o1.keys()),
                                    set(e for e in o2.keys()))
                for k in o1.keys():
                    almost_equal_dict(o1[k], o2[k])
            else:
                self.assertEqual(o1, o2)

        d = cp.get_neighbor_vertices_as_dict("Ca11Al14O32",
                                             elements=["Ca", "Al", "O"],
                                             comment="test_neighbor_1")
        cp.draw_diagram(remarked_compound="Ca11Al14O32", elements=["Ca", "Al", "O"])
        d_ref = {'compound': 'Ca11Al14O32',
                 'standard_energy': {Element('Ca'): -1.99870011,
                                     Element('Al'): -3.74799675,
                                     Element('O'): -4.93552791875},
                 'A': {Element('Ca'): -6.840267704583255,
                       Element('Al'): -8.729164283541706,
                       Element('O'): -1.4210854715202004e-14},
                 'B': {Element('Ca'): -6.646545131249994,
                       Element('Al'): -8.881374876874958,
                       Element('O'): 0.0},
                 'C': {Element('Ca'): -1.193725120000991,
                       Element('Al'): -0.702144860001539,
                       Element('O'): -5.45282001124896},
                 'D': {Element('Ca'): -1.1937251200010195,
                       Element('Al'): -0.25935040666836073,
                       Element('O'): -5.646542584582207},
                 'comment': 'test_neighbor_1'}
        almost_equal_dict(d, d_ref)
        d = cp.get_neighbor_vertices_as_dict("Ca",
                                             elements=["Ca", "Al", "O"],
                                             comment="test_neighbor_2")
        d_ref = \
            {'compound': 'Ca',
             'standard_energy': {Element('Ca'): -1.99870011,
                                 Element('Al'): -3.74799675,
                                 Element('O'): -4.93552791875},
             'A': {Element('Ca'): -1.4210854715202004e-14,
                   Element('Al'): -0.4947097849999835,
                   Element('O'): -6.646545131250008},
             'comment': 'test_neighbor_2'}
        almost_equal_dict(d, d_ref)

    def test_dump_yaml(self):
        dir_path = "./"

        filename = dir_path + "/vertices_MgO.yaml"
        cp = ChemPotDiag.from_file(FILENAME_2D)  # Mg, O
        comment = "This is from test_output_yaml in test_chem_pot.py"
        cp.dump_vertices_yaml(filename, "MgO", ["Mg", "O"], comment=comment)
        vertices, standard_energy = ChemPotDiag.load_vertices_yaml(filename)
        self.assertTrue(cp.get_neighbor_vertices("MgO", ["Mg", "O"]).
                        almost_equal(vertices))
        for elem in cp.elements:
            self.assertAlmostEqual(cp.element_energy[elem],
                                   standard_energy[elem])
        os.remove(filename)

        filename = dir_path + "/vertices_Ca11Al14O32.yaml"
        cp = ChemPotDiag.from_file(FILENAME_3D)  # Ca, Al, O
        comment = "This is from test_output_yaml in test_chem_pot.py"
        cp.dump_vertices_yaml(filename, "Ca11Al14O32", ["Ca", "Al", "O"],
                              comment=comment)
        vertices, standard_energy = ChemPotDiag.load_vertices_yaml(filename)
        self.assertTrue(cp.get_neighbor_vertices("Ca11Al14O32", ["Ca", "Al", "O"]).
                        almost_equal(vertices))
        for i, elem in enumerate(cp.elements):
            self.assertAlmostEqual(cp.element_energy[elem],
                                   standard_energy[elem])
        os.remove(filename)


if __name__ == "__main__":
    unittest.main()
