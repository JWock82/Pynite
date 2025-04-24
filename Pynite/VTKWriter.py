import vtk
import os
import tempfile
from pathlib import Path
import subprocess
from typing import Dict

from Pynite.FEModel3D import FEModel3D


class VTKWriter:
    def __init__(self, model: FEModel3D) -> None:
        self.model = model

    def write_to_vtk(self, path: str):
        """
        Writes model data into a VTK file using vtkUnstructuredGrid. The resulting file can
        then be postprocessed in ParaView. If multiple mesh types exist, they will be written as
        multiple .vtk files.

        Args:
            path (str): The path to the VTK file(s) to be created.
        """

        # remove Filetype if supplied explicitly
        path.removesuffix(".vtk")

        node_displacement_arrays: Dict[str, vtk.vtkFloatArray] = {}
        for combo_name in self.model._D.keys():
            node_displacement_arrays[combo_name] = vtk.vtkFloatArray()
            node_displacement_arrays[combo_name].SetName(f"Displacement_{combo_name}")
            node_displacement_arrays[combo_name].SetNumberOfComponents(3)

        node_ids: Dict[str, int] = {}
        node_name_array = vtk.vtkStringArray()
        node_name_array.SetName("Node_Names")

        points = vtk.vtkPoints()
        for node_name, node in self.model.nodes.items():
            point_id = points.InsertNextPoint(node.X, node.Y, node.Z)
            node_ids[node_name] = point_id
            node_name_array.InsertNextValue(node_name)

            for combo_name in self.model._D.keys():
                node_displacement_arrays[combo_name].InsertNextTuple3(node.DX[combo_name], node.DY[combo_name], node.DZ[combo_name])


        #### CREATE LINE CELLS ####
        lines = vtk.vtkCellArray()
        for member in self.model.members.values():
            for submember in member.sub_members.values():
                line = vtk.vtkLine()
                line.GetPointIds().SetId(0, node_ids[submember.i_node.name])
                line.GetPointIds().SetId(1, node_ids[submember.j_node.name])
                lines.InsertNextCell(line)

        #### CREATE QUAD CELLS ####
        quads = vtk.vtkCellArray()
        for quad in self.model.quads.values():
            vtkquad = vtk.vtkQuad()
            vtkquad.GetPointIds().SetId(0, node_ids[quad.i_node.name])
            vtkquad.GetPointIds().SetId(1, node_ids[quad.j_node.name])
            vtkquad.GetPointIds().SetId(2, node_ids[quad.m_node.name])
            vtkquad.GetPointIds().SetId(3, node_ids[quad.n_node.name])
            quads.InsertNextCell(vtkquad)

        ugrid_members = vtk.vtkUnstructuredGrid()
        ugrid_members.SetPoints(points)
        ugrid_members.GetPointData().AddArray(node_name_array)
        ugrid_members.SetCells(vtk.VTK_LINE, lines)

        ugrid_quads = vtk.vtkUnstructuredGrid()
        ugrid_quads.SetPoints(points)
        ugrid_quads.GetPointData().AddArray(node_name_array)
        ugrid_quads.SetCells(vtk.VTK_QUAD, quads)

        # writes displacement data into points
        for combo_name, array in node_displacement_arrays.items():
            ugrid_members.GetPointData().AddArray(array)

        #### EXTRACT QUAD DATA ####
        for combo in self.model.load_combos.keys():
            quad_data: Dict[str, vtk.vtkFloatArray] = {}
            for key in ("D", "F"): # Displacement, Force
                quad_data[key] = vtk.vtkFloatArray()
                quad_data[key].SetNumberOfComponents(3)
                quad_data[key].SetName(f'{key} "{combo}"')

            # fill zeros
            for node_id in node_ids.values():
                for array in quad_data.values():
                    array.InsertTuple3(node_id, 0, 0, 0)

            # get data and write into quad_data array
            for quad in self.model.quads.values():
                data = {"D": quad.D().reshape((4, 6)), "F": quad.F().reshape((4, 6))}

                for i,node in enumerate((quad.i_node, quad.j_node, quad.m_node, quad.n_node)):
                    for key in data.keys():
                        quad_data[key].InsertTuple3(node_ids[node.name], *data[key][i][:3])
            
            # add the data arrays to the points
            for array in quad_data.values():
                ugrid_quads.GetPointData().AddArray(array)

        #### WRITE DATA TO DISK ####
        if len(self.model.members)>0:
            member_writer = vtk.vtkUnstructuredGridWriter()
            member_writer.SetFileName(path + "_members.vtk")
            member_writer.SetInputData(ugrid_members)
            member_writer.Write()

        if len(self.model.quads)>0:
            quads_writer = vtk.vtkUnstructuredGridWriter()
            quads_writer.SetFileName(path + "_quads.vtk")
            quads_writer.SetInputData(ugrid_quads)
            quads_writer.Write()

    def open_in_paraview(self):
        """
        Open the Model inside Paraview, if installed. Paraview must be accessible with the `paraview` command for this method to work.
        """
        # Create a temporary file with .vtk extension
        with tempfile.NamedTemporaryFile(delete=False) as tmp:
            temp_path = tmp.name

        # Write the VTK file
        self.write_to_vtk(temp_path)
        files = []
        if len(self.model.members) > 0:
            files.append(temp_path + "_members.vtk")
        if len(self.model.quads) > 0:
            files.append(temp_path + "_quads.vtk")

        try:
            # Open paraview with the temporary file
            subprocess.run(
                ["paraview", *files],
                check=True,
            )
        finally:
            # Clean up the temporary file
            for f in files:
                os.remove(f)
