import vtk
import os
import tempfile
from pathlib import Path
import subprocess
from typing import Dict, Tuple

import numpy as np

from Pynite.FEModel3D import FEModel3D, Quad3D


class VTKWriter:
    """
    The `VTKWriter` Class allows for writing `FEModel3D` data into .vtk files. These can then be postprocessed
    using external Software i.e. Paraview.
    """

    def __init__(self, model: FEModel3D) -> None:
        self.model = model
        self.members_written = False
        self.quads_written = False

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

        self._write_member_data(path)
        self._write_quad_data(path)

    def _write_member_data(self, path: str):
        #### CREATE POINTS ####
        node_ids: Dict[str, int] = {}
        node_name_array = vtk.vtkStringArray()
        node_name_array.SetName("Node_Names")

        points = vtk.vtkPoints()
        for node_name, node in self.model.nodes.items():
            point_id = points.InsertNextPoint(node.X, node.Y, node.Z)
            node_ids[node_name] = point_id
            node_name_array.InsertNextValue(node_name)

        #### CREATE LINE CELLS ####
        lines = vtk.vtkCellArray()
        for member in self.model.members.values():
            for submember in member.sub_members.values():
                line = vtk.vtkLine()
                line.GetPointIds().SetId(0, node_ids[submember.i_node.name])
                line.GetPointIds().SetId(1, node_ids[submember.j_node.name])
                lines.InsertNextCell(line)

        ugrid_members = vtk.vtkUnstructuredGrid()
        ugrid_members.SetPoints(points)
        ugrid_members.GetPointData().AddArray(node_name_array)
        ugrid_members.SetCells(vtk.VTK_LINE, lines)

        if len(self.model.members) > 0:
            member_writer = vtk.vtkUnstructuredGridWriter()
            member_writer.SetFileName(path + "_members.vtk")
            member_writer.SetInputData(ugrid_members)
            member_writer.Write()
            self.members_written = True

    def _write_quad_data(self, path: str):
        #### CREATE QUAD CELLS ####
        quads = vtk.vtkCellArray()
        points = vtk.vtkPoints()
        for quad in self.model.quads.values():
            # xi and eta coordinates for the VTK_BIQUADRATIC_QUAD
            xis = (0,1,1,0,0.5,1,0.5,0,0.5)
            etas = (0,0,1,1,0,0.5,1,0.5,0.5)

            vtkquad = vtk.vtkBiQuadraticQuad()
            for i, (xi, eta) in enumerate(zip(xis, etas)):
                coords = self._interpolate_quad_coords(quad, xi, eta)
                vtkquad.GetPointIds().SetId(i, points.InsertNextPoint(coords))

            quads.InsertNextCell(vtkquad)

        ugrid_quads = vtk.vtkUnstructuredGrid()
        ugrid_quads.SetPoints(points)
        ugrid_quads.SetCells(vtk.VTK_BIQUADRATIC_QUAD, quads)

        #### WRITE DATA TO DISK ####

        if len(self.model.quads) > 0:
            quads_writer = vtk.vtkUnstructuredGridWriter()
            quads_writer.SetFileName(path + "_quads.vtk")
            quads_writer.SetInputData(ugrid_quads)
            quads_writer.Write()
            self.quads_written = True

    @staticmethod
    def _interpolate_quad_coords(quad: Quad3D, xi: float, eta: float) -> Tuple[float,float,float]:
        """
        Helper Method to return the global coordinates of a point on the quad, given the natural coordinates
        xi and eta. We should consider moving this over to Quad3D in the future.
        """
        pi = np.array([quad.i_node.X, quad.i_node.Y, quad.i_node.Z])
        pj = np.array([quad.j_node.X, quad.j_node.Y, quad.j_node.Z])
        pm = np.array([quad.m_node.X, quad.m_node.Y, quad.m_node.Z])
        pn = np.array([quad.n_node.X, quad.n_node.Y, quad.n_node.Z])

        return tuple(
            (1 - xi) * (1 - eta) * pi
            + xi * (1 - eta) * pj
            + xi * eta * pm
            + (1 - xi) * eta * pn
        )

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
        if self.members_written:
            files.append(temp_path + "_members.vtk")
        if self.quads_written:
            files.append(temp_path + "_quads.vtk")

        try:
            # Open paraview with the temporary file(s)
            subprocess.run(
                ["paraview", *files],
                check=True,
            )
        finally:
            # Clean up the temporary file
            for f in files:
                os.remove(f)
