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
        self._members_written = False
        self._quads_written = False

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

        #### Node Data ####
        for combo in self.model.load_combos.keys():
            node_D_array = vtk.vtkFloatArray()
            node_D_array.SetNumberOfComponents(3)
            node_D_array.SetName(f"D - {combo}")
            for node_name, node_id in node_ids.items():
                node = self.model.nodes[node_name]
                node_D_array.InsertTuple3(node_id, node.DX[combo], node.DY[combo], node.DZ[combo])
            ugrid_members.GetPointData().AddArray(node_D_array)

        if len(self.model.members) > 0:
            member_writer = vtk.vtkUnstructuredGridWriter()
            member_writer.SetFileName(path + "_members.vtk")
            member_writer.SetInputData(ugrid_members)
            member_writer.Write()
            self._members_written = True

    def _write_quad_data(self, path: str):

        # xi and eta natural coordinates [0-1] for the VTK_BIQUADRATIC_QUAD point positions
        xis = (0,1,1,0,0.5,1,0.5,0,0.5)
        etas = (0,0,1,1,0,0.5,1,0.5,0.5)

        #### CREATE QUAD CELLS ####
        quads = vtk.vtkCellArray()
        points = vtk.vtkPoints()
        quad_refs: Dict[str, vtk.vtkBiQuadraticQuad] = {}
        node_register: Dict[int,Tuple[float,float,float]] = {} # Contains all existing point ids and their positions
        for quad in self.model.quads.values():
            # Node corner coords
            pi = np.array([quad.i_node.X, quad.i_node.Y, quad.i_node.Z])
            pj = np.array([quad.j_node.X, quad.j_node.Y, quad.j_node.Z])
            pm = np.array([quad.m_node.X, quad.m_node.Y, quad.m_node.Z])
            pn = np.array([quad.n_node.X, quad.n_node.Y, quad.n_node.Z])

            vtkquad = vtk.vtkBiQuadraticQuad()
            for i, (xi, eta) in enumerate(zip(xis, etas)):
                coords = self._interpolate_quad_corner_data(pi, pj, pm, pn, xi, eta)
                # search existing nodes by position if node already exists (from other quad)
                for node_id, node_coords in node_register.items():
                    if sum(abs(coords[i] - node_coords[i]) for i in range(3)) <= 1e-10:
                        vtkquad.GetPointIds().SetId(i, node_id)
                        break
                else:
                    new_id = points.InsertNextPoint(*coords)
                    node_register[new_id] = coords
                    vtkquad.GetPointIds().SetId(i, new_id)

            quad_refs[quad.name] = vtkquad
            quads.InsertNextCell(vtkquad)

        ugrid_quads = vtk.vtkUnstructuredGrid()
        ugrid_quads.SetPoints(points)
        ugrid_quads.SetCells(vtk.VTK_BIQUADRATIC_QUAD, quads)

        #### READ QUAD DATA ####
        for quad_name, vtkquad in quad_refs.items():
            for combo in self.model.load_combos.keys():
                quad = self.model.quads[quad_name]

                # DISPLACEMENT
                D = vtk.vtkFloatArray()
                D.SetName(f"D - {combo}")
                D.SetNumberOfComponents(3)
                for i,(xi,eta) in enumerate(zip(xis,etas)):
                    di = np.array([quad.i_node.DX[combo], quad.i_node.DY[combo], quad.i_node.DZ[combo]])
                    dj = np.array([quad.j_node.DX[combo], quad.j_node.DY[combo], quad.j_node.DZ[combo]])
                    dm = np.array([quad.m_node.DX[combo], quad.m_node.DY[combo], quad.m_node.DZ[combo]])
                    dn = np.array([quad.n_node.DX[combo], quad.n_node.DY[combo], quad.n_node.DZ[combo]])

                    d = self._interpolate_quad_corner_data(di, dj, dm, dn, xi, eta)
                    D.InsertTuple3(vtkquad.GetPointId(i), *d)
                ugrid_quads.GetPointData().AddArray(D)

                # MEMBRANE STRESSES
                membrane = vtk.vtkFloatArray()
                membrane.SetName(f"Membrane - {combo}")
                membrane.SetNumberOfComponents(3)
                for i, (xi, eta) in enumerate(zip(xis, etas)):
                    membrane.InsertTuple3(vtkquad.GetPointId(i), *quad.membrane(xi, eta, False, combo)) # type: ignore
                ugrid_quads.GetPointData().AddArray(membrane)

        #### WRITE DATA TO DISK ####

        if len(self.model.quads) > 0:
            quads_writer = vtk.vtkUnstructuredGridWriter()
            quads_writer.SetFileName(path + "_quads.vtk")
            quads_writer.SetInputData(ugrid_quads)
            quads_writer.Write()
            self._quads_written = True

    @staticmethod
    def _interpolate_quad_corner_data(i:np.ndarray, j:np.ndarray, m:np.ndarray, n:np.ndarray, xi:float, eta:float) -> Tuple[float,float,float]:
        """
        Helper Method to return the linearly interpolated data of a point on the quad, given the natural coordinates
        xi and eta. We should consider moving this over to Quad3D in the future.
        """
        return tuple(
            (1 - xi) * (1 - eta) * i
            + xi * (1 - eta) * j
            + xi * eta * m
            + (1 - xi) * eta * n
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
        if self._members_written:
            files.append(temp_path + "_members.vtk")
        if self._quads_written:
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
