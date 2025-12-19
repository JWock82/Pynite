
import unittest
from Pynite import FEModel3D
from Pynite.VTKWriter import VTKWriter
import sys, tempfile
from io import StringIO

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
