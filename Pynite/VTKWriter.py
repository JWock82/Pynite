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
        self._nodes_written = False
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

        self._write_node_data(path)
        self._write_member_data(path)
        self._write_quad_data(path)

    def _write_node_data(self, path:str):
        if self.log:
            print("- collecting node data...")

        ugrid = vtk.vtkUnstructuredGrid()
        points = vtk.vtkPoints()

        node_names = vtk.vtkStringArray()
        node_names.SetName("Name")
        
        node_cells = vtk.vtkCellArray()
        node_ids: Dict[str, int] = {}
        for node in self.model.nodes.values():
            point_id = points.InsertNextPoint(node.X, node.Y, node.Z)
            node_ids[node.name] = point_id
            node_names.InsertValue(point_id, node.name)
            vert = vtk.vtkVertex()
            vert.GetPointIds().SetId(0, point_id)
            node_cells.InsertNextCell(vert)

        ugrid.SetPoints(points)
        ugrid.SetCells(vtk.VTK_POINT_DATA, node_cells)
        ugrid.GetPointData().AddArray(node_names)
        
        #### LOAD SPECIFIC DATA ####
        for combo in self.model.load_combos.keys():
            reaction_constraints = vtk.vtkIntArray()
            reaction_constraints.SetName("Reaction Constraints")
            reaction_constraints.SetNumberOfComponents(6)

            displacements = vtk.vtkDoubleArray()
            displacements.SetName("Displacements D")
            displacements.SetNumberOfComponents(3)

            forces = vtk.vtkDoubleArray()
            forces.SetName(f"Force Reactions F - {combo}")
            forces.SetNumberOfComponents(3)

            moments = vtk.vtkDoubleArray()
            moments.SetName(f"Moment Reactions M - {combo}")
            moments.SetNumberOfComponents(3)

            for i, name in enumerate(["DX", "DY", "DZ", "RX", "RY", "RZ"]):
                reaction_constraints.SetComponentName(i,name)
            
            for node_name, node_id in node_ids.items():
                node = self.model.nodes[node_name]
                reaction_constraints.InsertTuple6(node_id, int(node.support_DX), int(node.support_DY), int(node.support_DZ), int(node.support_RX), int(node.support_RY), int(node.support_RZ))
                displacements.InsertTuple3(node_id, node.DX[combo], node.DY[combo], node.DZ[combo]) # type: ignore
                forces.InsertTuple3(node_id, node.RxnFX[combo], node.RxnFY[combo], node.RxnFZ[combo])
                moments.InsertTuple3(node_id, node.RxnMX[combo], node.RxnMY[combo], node.RxnMZ[combo])
                
            ugrid.GetPointData().AddArray(reaction_constraints)
            ugrid.GetPointData().AddArray(displacements)
            ugrid.GetPointData().AddArray(forces)
            ugrid.GetPointData().AddArray(moments)
        

        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(path + "_nodes.vtk")
        writer.SetInputData(ugrid)
        writer.Write()
        self._nodes_written = True
        

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
            D_array_G = vtk.vtkDoubleArray()
            D_array_G.SetNumberOfComponents(3)
            D_array_G.SetName(f"Displacement D - {combo}")

            D_array_lok = vtk.vtkDoubleArray()
            D_array_lok.SetNumberOfComponents(3)
            D_array_lok.SetName(f"Displacement d - {combo}")

            # Moments
            moment_G = vtk.vtkDoubleArray()
            moment_G.SetNumberOfComponents(3)
            moment_G.SetName(f"Moments M - {combo}")

            moment_lok = vtk.vtkDoubleArray()
            moment_lok.SetNumberOfComponents(3)
            moment_lok.SetName(f"Moments m - {combo}")

            # Forces
            force_G = vtk.vtkDoubleArray()
            force_G.SetNumberOfComponents(3)
            force_G.SetName(f"Forces F - {combo}")

            force_lok = vtk.vtkDoubleArray()
            force_lok.SetNumberOfComponents(3)
            force_lok.SetName(f"Forces f - {combo}")

            # Bending Stress increase
            # Can be used with the Paraview Calculator Filter to get bending stresses at a given location by multiplying with the axial distance
            sigma_b_G = vtk.vtkDoubleArray()
            sigma_b_G.SetNumberOfComponents(3)
            sigma_b_G.SetName(f"Sigma/r - {combo}")

            sigma_b_lok = vtk.vtkDoubleArray()
            sigma_b_lok.SetNumberOfComponents(3)
            sigma_b_lok.SetName(f"sigma/r - {combo}")

            for line_range,subm,line in submembers:
                for i,x in enumerate(line_range):
                    x = 1-x # go backwards
                    xl = x * subm.L()
                    point_id = line.GetPointId(i)
                    T = inv(subm.T()[:3,:3]) # Transformation Matrix Local -> Global

                    # Displacement
                    deflection = np.array([float(subm.deflection(direction,xl,combo)) for direction in ("dx", "dy", "dz")]) # type: ignore
                    D_array_lok.InsertTuple3(point_id, *deflection)
                    D_array_G.InsertTuple3(point_id, *(T @ deflection))

                    # moment
                    m = np.array([subm.torque(xl, combo), subm.moment("My",xl,combo), subm.moment("Mz",xl,combo)])
                    moment_lok.InsertTuple3(point_id, *m)
                    moment_G.InsertTuple3(point_id, *(T @ m))

                    # forces
                    s = np.array([subm.axial(xl, combo), subm.shear("Fy",xl,combo), subm.shear("Fz",xl,combo)])
                    force_lok.InsertTuple3(point_id, *s)
                    force_G.InsertTuple3(point_id, *(T @ s))

                    # bending stress increase
                    sec = subm.section
                    sig_b = np.array([0, subm.moment("My", xl, combo)/sec.Iy, subm.moment("Mz", xl, combo)/sec.Iz]) # type: ignore
                    sigma_b_lok.InsertTuple3(point_id, *sig_b)
                    sigma_b_G.InsertTuple3(point_id, *(T @ sig_b))

            ugrid_members.GetPointData().AddArray(D_array_G)
            ugrid_members.GetPointData().AddArray(D_array_lok)
            ugrid_members.GetPointData().AddArray(moment_G)
            ugrid_members.GetPointData().AddArray(moment_lok)
            ugrid_members.GetPointData().AddArray(force_G)
            ugrid_members.GetPointData().AddArray(force_lok)
            ugrid_members.GetPointData().AddArray(sigma_b_G)
            ugrid_members.GetPointData().AddArray(sigma_b_lok)

        # clean the data from duplicate points
        cleaner = vtk.vtkStaticCleanUnstructuredGrid()
        cleaner.SetInputData(ugrid_members)
        cleaner.SetToleranceIsAbsolute(True)
        cleaner.SetAbsoluteTolerance(0.01)
        cleaner.Update()
        ugrid_members = cleaner.GetOutput()

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
                node_id = points.InsertNextPoint(*coords)

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

        # clean the data from duplicate points
        cleaner = vtk.vtkStaticCleanUnstructuredGrid()
        cleaner.SetInputData(ugrid_quads)
        cleaner.SetToleranceIsAbsolute(True)
        cleaner.SetAbsoluteTolerance(0.01)
        cleaner.Update()
        ugrid_quads = cleaner.GetOutput()

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
        This can be used i.e. to calculate every point position on a member given by a relative coordinate.
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
        if self._nodes_written:
            files.append(temp_path + "_nodes.vtk")

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
