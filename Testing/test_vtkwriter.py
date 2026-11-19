
import unittest
from Pynite import FEModel3D
from Pynite.VTKWriter import VTKWriter
import sys, tempfile
from io import StringIO
import pytest
import vtk

class Test_VTKWriter(unittest.TestCase):
    '''
    Tests of writing VTK Data to disk.
    This Testsuit can only check for errors when writing to disk,
    not check if the data is correct or can even be read.
    '''

    def setUp(self):
        
        # Suppress printed output temporarily
        sys.stdout = StringIO()

    def tearDown(self):

        # Reset the print function to normal
        sys.stdout = sys.__stdout__
    
    def test_nodes(self):

        # Create a temporary file
        with tempfile.NamedTemporaryFile(suffix=".vtk") as tmp:
            temp_path = tmp.name
        
            model = FEModel3D()
            model.add_material("Steel", 210e9, 210e9, 0.3, 7850)
            model.add_section("test",1,1,1,1)
            model.add_node("A", 0, 0, 0)
            model.add_node("B", 1, 0, 0)

            model.add_member("AB", "A", "B", "Steel", "test")

            model.def_support("A",True,True,True,True,True,True)
            model.add_node_load("B","FY", 1)
            model.analyze()

            VTKWriter(model)._write_node_data(temp_path)

    def test_members(self):

        with tempfile.NamedTemporaryFile(suffix=".vtk") as tmp:
            temp_path = tmp.name
        
            model = FEModel3D()
            model.add_material("Steel", 210e9, 210e9, 0.3, 7850)
            model.add_section("test",1,1,1,1)
            model.add_node("A", 0, 0, 0)
            model.add_node("B", 1, 0, 0)

            model.add_member("AB", "A", "B", "Steel", "test")

            model.def_support("A",True,True,True,True,True,True)
            model.add_node_load("B","FY", 1)
            model.analyze()

            VTKWriter(model)._write_member_data(temp_path)

    def test_quads(self):

        with tempfile.NamedTemporaryFile(suffix=".vtk") as tmp:
            temp_path = tmp.name
        
            model = FEModel3D()
            model.add_material("Steel", 210e9, 210e9, 0.3, 7850)
            model.add_section("test",1,1,1,1)
            model.add_node("A", 0, 0, 0)
            model.add_node("B", 1, 0, 0)
            model.add_node("C", 1, 1, 0)
            model.add_node("D", 0, 1, 0)
            model.add_quad("ABCD","A", "B", "C", "D",1,"Steel")

            model.add_member("AB", "A", "B", "Steel", "test")

            model.def_support("A",True,True,True,True,True,True)
            model.add_node_load("B","FY", 1)
            model.analyze()

            VTKWriter(model)._write_quad_data(temp_path)


@pytest.mark.parametrize('factors, expected_force, expected_moment', [
    ({'D': 1.0}, (2, -12, 4), (6, 8, -2)),
    ({'D': 1.5, 'L': -0.5, 'W': 0.0}, (7, -16, 5), (7, 11, -6)),
    ({}, (0, 0, 0), (0, 0, 0)),
    ({'D': 1.0, 'L': 1.0, 'W': 1.0}, (94, 84, 106), (110, 110, 104)),
])
def test_nodal_load_combination_factors(tmp_path, factors, expected_force, expected_moment):
    """Read back exported nodal forces and moments for factored and omitted load cases."""
    model = FEModel3D()
    model.add_material('Steel', 210e9, 80e9, 0.3, 7850)
    model.add_section('Section', 1, 1, 1, 1)
    model.add_node('A', 0, 0, 0)
    model.add_node('B', 1, 0, 0)
    model.add_member('AB', 'A', 'B', 'Steel', 'Section')
    model.def_support('A', True, True, True, True, True, True)

    # Exercise all six directions, repeated loads, and a case excluded from some combinations.
    for case, magnitudes in [('D', (2, -10, 4, 6, 8, -2)),
                             ('L', (-8, -4, 2, 4, 2, 6)), ('W', (100,) * 6)]:
        for direction, magnitude in zip(('FX', 'FY', 'FZ', 'MX', 'MY', 'MZ'), magnitudes):
            model.add_node_load('B', direction, magnitude, case)
    model.add_node_load('B', 'FY', -2, 'D')
    model.add_load_combo('Test', factors)
    model.analyze_linear()

    # Use the public exporter and a real VTK reader, without opening a rendering window.
    path = tmp_path / 'model'
    VTKWriter(model).write_to_vtk(str(path))
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(str(path) + '_nodes.vtk')
    reader.ReadAllFieldsOn()
    reader.Update()
    data = reader.GetOutput().GetPointData()

    # Point 0 is unloaded; point 1 must contain only the factored loads for this combination.
    for label, expected in [('F', expected_force), ('M', expected_moment)]:
        loads = data.GetArray(f'Loads {label} - Test')
        assert loads.GetTuple3(0) == pytest.approx((0, 0, 0))
        assert loads.GetTuple3(1) == pytest.approx(expected)


@pytest.mark.parametrize('combos', [{'Service': 1.0}, {'Service': 1.0, 'Strength': 1.5}])
def test_nodal_displacements_per_combination(tmp_path, combos):
    """Keep separately named displacement arrays for each load combination."""
    model = FEModel3D()
    model.add_material('Material', 1000, 400, 0.25, 1)
    model.add_section('Section', 1, 1, 1, 1)
    model.add_node('A', 0, 0, 0)
    model.add_node('B', 3, 0, 0)
    model.add_member('AB', 'A', 'B', 'Material', 'Section')
    model.def_support('A', True, True, True, True, True, True)
    model.add_node_load('B', 'FY', -10, 'D')
    for combo, factor in combos.items():
        model.add_load_combo(combo, {'D': factor})
    model.analyze_linear()

    # Cantilever tip deflection is PL^3/(3EI) = -0.09 for the unit-factor combination.
    path = tmp_path / 'model'
    VTKWriter(model).write_to_vtk(str(path))
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(str(path) + '_nodes.vtk')
    reader.ReadAllFieldsOn()
    reader.Update()
    data = reader.GetOutput().GetPointData()
    for combo, factor in combos.items():
        assert model.nodes['B'].DY[combo] == pytest.approx(-0.09 * factor)
        displacements = data.GetArray(f'Displacements D - {combo}')
        assert displacements is not None
        assert displacements.GetTuple3(0) == pytest.approx((0, 0, 0))
        assert displacements.GetTuple3(1) == pytest.approx((0, -0.09 * factor, 0))
