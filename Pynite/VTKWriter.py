import vtk
import os
import tempfile
from typing import Dict

from Pynite.FEModel3D import FEModel3D


class VTKWriter:
    def __init__(self, model: FEModel3D) -> None:
        self.model = model

    def write_to_vtk(self, path: str):
        """
        Writes model data into a VTK file using vtkUnstructuredGrid. The resulting file can
        then be postprocessed in ParaView.

        Args:
            path (str): The path to the VTK file to be created.
        """

        # Sanitize path if necessary
        if not path.endswith(".vtk"):
            path += ".vtk"

        # Initialize empty containers
        ugrid = vtk.vtkUnstructuredGrid()
        points = vtk.vtkPoints()
        verts = vtk.vtkCellArray()

        displacement_arrays: Dict[str, vtk.vtkFloatArray] = {}
        for combo_name in self.model._D.keys():
            displacement_arrays[combo_name] = vtk.vtkFloatArray()
            displacement_arrays[combo_name].SetName(f"Displacement_{combo_name}")
            displacement_arrays[combo_name].SetNumberOfComponents(3)

        node_ids = {}
        for i, (node_name, node) in enumerate(self.model.nodes.items()):
            point_id = points.InsertNextPoint(node.X, node.Y, node.Z)
            node_ids[node_name] = point_id

            verts.InsertNextCell(1)
            verts.InsertCellPoint(point_id)

            for combo_name in self.model._D.keys():
                try:
                    DX = self.model._D[combo_name][node.ID * 6 + 0, 0]
                    DY = self.model._D[combo_name][node.ID * 6 + 1, 0]
                    DZ = self.model._D[combo_name][node.ID * 6 + 2, 0]
                except:
                    DX = 0.0
                    DY = 0.0
                    DZ = 0.0
                displacement_arrays[combo_name].InsertNextTuple3(DX, DY, DZ)

        ugrid.SetPoints(points)
        ugrid.SetCells(vtk.VTK_VERTEX, verts)

        for combo_name, array in displacement_arrays.items():
            ugrid.GetPointData().AddArray(array)

        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(path)
        writer.SetInputData(ugrid)
        writer.Write()

    def open_in_paraview(self):
        raise NotImplementedError
