import vtk
import os
import tempfile
from pathlib import Path
import subprocess
from typing import Dict, Tuple, List

import numpy as np
from numpy.linalg import inv

from Pynite.FEModel3D import FEModel3D, Quad3D, Member3D


class VTKWriter:
    """
    The `VTKWriter` Class allows for writing `FEModel3D` data into .vtk files. These can then be postprocessed
    using external Software i.e. Paraview.
    """

    def __init__(self, model: FEModel3D, log=False) -> None:
        self.model = model
        self._members_written = False
        self._quads_written = False
        self.log = log

    def write_to_vtk(self, path: str):
        """
        Writes model data into a VTK file using vtkUnstructuredGrid. The resulting file can
        then be postprocessed in ParaView. If multiple mesh types exist, they will be written as
        multiple .vtk files.

        Args:
            path (str): The path to the VTK file(s) to be created.
        """

        # remove Filetype if supplied explicitly
        path = path.removesuffix(".vtk")
        if self.log:
            print(f"Writing Data to {path}...")

        self._write_member_data(path)
        self._write_quad_data(path)

    def _write_member_data(self, path: str):
        if self.log:
            print("- collecting member data...")

        points = vtk.vtkPoints()

        #### CREATE LINE CELLS ####
        # each (sub)member is subdivided into line segments
        lines = vtk.vtkCellArray()
        submembers: List[Tuple[Tuple[float, float], Member3D, vtk.vtkLine]] = []
        for member in self.model.members.values():
            for subm in member.sub_members.values():
                n = 11 # Number of submembers + 1
                for start,end in zip(np.linspace(0,1,n)[:-1], np.roll(np.linspace(0,1,n),-1)[:-1]):
                    # Calculate intermediary point positions
                    p1 = self._interpolate_member_data(
                        np.array([subm.i_node.X, subm.i_node.Y, subm.i_node.Z]),
                        np.array([subm.j_node.X, subm.j_node.Y, subm.j_node.Z]),
                        start,
                    )
                    p2 = self._interpolate_member_data(
                        np.array([subm.i_node.X, subm.i_node.Y, subm.i_node.Z]),
                        np.array([subm.j_node.X, subm.j_node.Y, subm.j_node.Z]),
                        end,
                    )

                    line = vtk.vtkLine()
                    line.SetObjectName(member.name)
                    line.GetPointIds().SetId(0, points.InsertNextPoint(*p1))
                    line.GetPointIds().SetId(1, points.InsertNextPoint(*p2))
                    submembers.append(((start, end), subm, line))
                    lines.InsertNextCell(line)

        ugrid_members = vtk.vtkUnstructuredGrid()
        ugrid_members.SetPoints(points)
        ugrid_members.SetCells(vtk.VTK_LINE, lines)

        #### MEMBER Data ####
        for combo in self.model.load_combos.keys():
            # Displacement
            D_array = vtk.vtkDoubleArray()
            D_array.SetNumberOfComponents(3)
            D_array.SetName(f"Displacement - {combo}")

            # Moments
            moment = vtk.vtkDoubleArray()
            moment.SetNumberOfComponents(3)
            moment.SetName(f"Moments - {combo}")

            # Forces
            force = vtk.vtkDoubleArray()
            force.SetNumberOfComponents(3)
            force.SetName(f"Forces - {combo}")

            # Bending Stress increase
            # Can be used with the Paraview Calculator Filter to get bending stresses at a given location by multiplying with the axial distance
            sigma_b = vtk.vtkDoubleArray()
            sigma_b.SetNumberOfComponents(3)
            sigma_b.SetName(f"Sigma/r - {combo}")

            for line_range,subm,line in submembers:
                for i,x in enumerate(line_range):
                    x = 1-x # go backwards
                    point_id = line.GetPointId(i)
                    T = inv(subm.T()[:3,:3]) # Transformation Matrix Local -> Global
                    # Displacement
                    deflection = T @ np.array([float(subm.deflection(direction,x*subm.L(),combo)) for direction in ("dx", "dy", "dz")]) # type: ignore
                    D_array.InsertTuple3(point_id, *deflection)

                    # moment
                    res = np.array([subm.torque(x, combo), subm.moment("My",x,combo), subm.moment("Mz",x,combo)])
                    if all(res):
                        m = T @ res
                    else:
                        m = (0,0,0)
                    moment.InsertTuple3(point_id, *m)

                    # forces
                    res = np.array([subm.axial(x, combo), subm.shear("Fy",x,combo), subm.shear("Fz",x,combo)])
                    if all(res):
                        s = T @ res
                    else:
                        s = (0,0,0)
                    force.InsertTuple3(point_id, *s)

                    # bending stress increase
                    sec = subm.section
                    my = subm.moment("My", x, combo)
                    mz = subm.moment("Mz", x, combo)
                    if all((my,mz)):
                        res = np.array([0, my/sec.Iy, mz/sec.Iz]) # type: ignore
                        sig_b = T @ res
                    else:
                        sig_b = (0,0,0)
                    sigma_b.InsertTuple3(point_id, *sig_b)

            ugrid_members.GetPointData().AddArray(D_array)
            ugrid_members.GetPointData().AddArray(moment)
            ugrid_members.GetPointData().AddArray(force)
            ugrid_members.GetPointData().AddArray(sigma_b)

        if len(self.model.members) > 0:
            member_writer = vtk.vtkUnstructuredGridWriter()
            member_writer.SetFileName(path + "_members.vtk")
            member_writer.SetInputData(ugrid_members)
            member_writer.Write()
            self._members_written = True
            if self.log:
                print("- members written!")
        elif self.log:
            print("- no members detected")

    def _write_quad_data(self, path: str):
        """
        Collects all quads in the model and writes their data into a VTK mesh.
        """
        if self.log:
            print("- collecting quad data...")
        # xi and eta natural coordinates [0-1] for the VTK_BIQUADRATIC_QUAD point positions
        xis = (0,1,1,0,0.5,1,0.5,0,0.5)
        etas = (0,0,1,1,0,0.5,1,0.5,0.5)

        #### CREATE QUAD CELLS ####
        quads = vtk.vtkCellArray()
        points = vtk.vtkPoints()
        quad_refs: Dict[str, vtk.vtkBiQuadraticQuad] = {}

        node_register = np.empty((0,4), dtype=np.float64) # (i,x,y,z)

        for quad in self.model.quads.values():
            # Node corner coords
            pi = np.array([quad.i_node.X, quad.i_node.Y, quad.i_node.Z])
            pj = np.array([quad.j_node.X, quad.j_node.Y, quad.j_node.Z])
            pm = np.array([quad.m_node.X, quad.m_node.Y, quad.m_node.Z])
            pn = np.array([quad.n_node.X, quad.n_node.Y, quad.n_node.Z])

            vtkquad = vtk.vtkBiQuadraticQuad()
            vtkquad.SetObjectName(quad.name)
            for i, (xi, eta) in enumerate(zip(xis, etas)):
                coords = np.array(self._interpolate_quad_corner_data(pi, pj, pm, pn, xi, eta))
                # search existing nodes by position if node already exists (from other quad)
                if node_register.shape[0] > 0:
                    # This was done in Numpy for necessary performance reasons, hence the unintuitive code.
                    # The (x,y,z) coords array gets subtracted from every entry in the register and summed up per entry,
                    # resulting in a vector of absolute distances. This is then checked for the smallest entry, that needs to
                    # be below a numeric threshold to be counted as a match. The index can then be used to get the node_id.
                    search = np.sum(np.abs(node_register[:,1:] - coords), axis=1)
                    query = node_register[:,0][np.where((search == np.min(search)) & (np.min(search)<1e-10))[0]]
                    if len(query)>0:
                        node_id = int(query[0])
                        # existing node found on position coords
                        vtkquad.GetPointIds().SetId(i, node_id)
                    else:
                        # no existing node found at positiion coords, create a new point
                        node_id = points.InsertNextPoint(*coords)
                        node_register = np.vstack((node_register,np.array([node_id, *coords])))
                else:
                    # in case of very first node
                    node_id = points.InsertNextPoint(*coords)
                    node_register = np.vstack((node_register,np.array([node_id, *coords])))

                vtkquad.GetPointIds().SetId(i, node_id)

            quad_refs[quad.name] = vtkquad
            quads.InsertNextCell(vtkquad)

        ugrid_quads = vtk.vtkUnstructuredGrid()
        ugrid_quads.SetPoints(points)
        ugrid_quads.SetCells(vtk.VTK_BIQUADRATIC_QUAD, quads)

        #### READ QUAD DATA ####
        for combo in self.model.load_combos.keys():
            # Displacement Data
            D = vtk.vtkDoubleArray()
            D.SetName(f"Displacement - {combo}")
            D.SetNumberOfComponents(3)

            # Membrane Data
            membrane = vtk.vtkDoubleArray()
            membrane.SetName(f"Membrane - {combo}")
            membrane.SetNumberOfComponents(3)

            # Moment Data
            moments = vtk.vtkDoubleArray()
            moments.SetName(f"Moments - {combo}")
            moments.SetNumberOfComponents(3)

            # SHEAR Data
            shear = vtk.vtkDoubleArray()
            shear.SetName(f"Shear - {combo}")
            shear.SetNumberOfComponents(3)

            for quad_name, vtkquad in quad_refs.items():
                quad = self.model.quads[quad_name]

                for i,(xi,eta) in enumerate(zip(xis,etas)):
                    # DISPLACEMENT
                    di = np.array([quad.i_node.DX[combo], quad.i_node.DY[combo], quad.i_node.DZ[combo]])
                    dj = np.array([quad.j_node.DX[combo], quad.j_node.DY[combo], quad.j_node.DZ[combo]])
                    dm = np.array([quad.m_node.DX[combo], quad.m_node.DY[combo], quad.m_node.DZ[combo]])
                    dn = np.array([quad.n_node.DX[combo], quad.n_node.DY[combo], quad.n_node.DZ[combo]])

                    d = self._interpolate_quad_corner_data(di, dj, dm, dn, xi, eta)
                    D.InsertTuple3(vtkquad.GetPointId(i), *d)

                    # MEMBRANE STRESSES
                    membrane.InsertTuple3(vtkquad.GetPointId(i), *quad.membrane(xi, eta, False, combo).flatten()) # type: ignore

                    # MOMENTS
                    m = quad.moment(xi, eta, False, combo) # type: ignore
                    moments.InsertTuple3(vtkquad.GetPointId(i), *m)

                    # SHEAR
                    s = quad.shear(xi, eta, False, combo) # type: ignore
                    shear.InsertTuple3(vtkquad.GetPointId(i), *s)

                ugrid_quads.GetPointData().AddArray(D)
                ugrid_quads.GetPointData().AddArray(membrane)
                ugrid_quads.GetPointData().AddArray(moments)
                ugrid_quads.GetPointData().AddArray(shear)

        #### WRITE DATA TO DISK ####

        if len(self.model.quads) > 0:
            quads_writer = vtk.vtkUnstructuredGridWriter()
            quads_writer.SetFileName(path + "_quads.vtk")
            quads_writer.SetInputData(ugrid_quads)
            quads_writer.Write()
            self._quads_written = True
            if self.log:
                print("- quads written!")
        elif self.log:
            print("- no quads detected")

    @staticmethod
    def _interpolate_member_data(p1:np.ndarray,p2:np.ndarray, x:float) -> np.ndarray:
        """
        Takes in vectors of size n from both ends of the member and interpolates linearly to relative position x on the member [0-1].
        """
        return x*p1 + (1-x)*p2

    @staticmethod
    def _interpolate_quad_corner_data(i:np.ndarray, j:np.ndarray, m:np.ndarray, n:np.ndarray, xi:float, eta:float) -> Tuple[float,float,float]:
        """
        Helper Method to return the linearly interpolated data of a point on the quad, given the natural coordinates
        xi and eta. We should consider moving this over to Quad3D in the future.
        """
        result = tuple(
            (1 - xi) * (1 - eta) * i
            + xi * (1 - eta) * j
            + xi * eta * m
            + (1 - xi) * eta * n
        )
        return tuple(float(res) for res in result) # type: ignore

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
