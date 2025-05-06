"""
This module contains the VTKWriter class, which fetches the FEModel3D data and writes it into one or multiple
`.vtk` files. The VTK fileformat is a commonly supported format designed for postprocessing of mesh results.
"""

from __future__ import annotations # Allows more recent type hints features

import vtk
import os
import tempfile
import subprocess

import numpy as np
from numpy.linalg import inv
from scipy.interpolate import RegularGridInterpolator

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from typing import Dict, Tuple, List
    from Pynite.FEModel3D import FEModel3D, Member3D

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
        
        if len(self.model.plates)>0:
            print("[WARNING] Plate export is not supported at this time. If needed, you can substitude them by Quad3D Elements.")
        
        self._write_node_data(path+"_nodes.vtk")

        if len(self.model.members)>0:
            self._write_member_data(path+"_members.vtk")
        elif self.log:
            print("- no members detected")

        if len(self.model.quads)>0:
            self._write_quad_data(path+"_quads.vtk")
        elif self.log:
            print("- no quads detected")

    def _write_node_data(self, path:str):
        """
        Collect all Nodes inside the model and write their data into a .vtk file.
        """
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

            force_loads = vtk.vtkDoubleArray()
            force_loads.SetName(f"Loads F - {combo}")
            force_loads.SetNumberOfComponents(3)

            moment_loads = vtk.vtkDoubleArray()
            moment_loads.SetName(f"Loads M - {combo}")
            moment_loads.SetNumberOfComponents(3)

            for i, name in enumerate(["DX", "DY", "DZ", "RX", "RY", "RZ"]):
                reaction_constraints.SetComponentName(i,name)

            for node_name, node_id in node_ids.items():
                node = self.model.nodes[node_name]
                reaction_constraints.InsertTuple6(node_id, int(node.support_DX), int(node.support_DY), int(node.support_DZ), int(node.support_RX), int(node.support_RY), int(node.support_RZ))
                displacements.InsertTuple3(node_id, node.DX[combo], node.DY[combo], node.DZ[combo]) # type: ignore
                forces.InsertTuple3(node_id, node.RxnFX[combo], node.RxnFY[combo], node.RxnFZ[combo])
                moments.InsertTuple3(node_id, node.RxnMX[combo], node.RxnMY[combo], node.RxnMZ[combo])

                # calculate the NodeLoad for each node
                fl: Dict[str,float] = {"X":0,"Y":0,"Z":0}
                ml: Dict[str,float] = {"X":0,"Y":0,"Z":0}
                for (f_or_m, direction), magnitude, case in node.NodeLoads:
                    if f_or_m == "F":
                        fl[direction] += magnitude
                    elif f_or_m == "M":
                        ml[direction] += magnitude

                force_loads.InsertTuple3(node_id, fl["X"], fl["Y"], fl["Z"])
                moment_loads.InsertTuple3(node_id, ml["X"], ml["Y"], ml["Z"])

            ugrid.GetPointData().AddArray(reaction_constraints)
            ugrid.GetPointData().AddArray(displacements)
            ugrid.GetPointData().AddArray(forces)
            ugrid.GetPointData().AddArray(moments)
            ugrid.GetPointData().AddArray(force_loads)
            ugrid.GetPointData().AddArray(moment_loads)

        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(path)
        writer.SetInputData(ugrid)
        writer.Write()
        self._nodes_written = True

        if self.log:
            print("- nodes written!")

    def _write_member_data(self, path: str):
        """
        Collects all members inside the model and writes their data into a .vtk file. Each member is represented by 10 lines
        that sample its data along its length.
        """
        if self.log:
            print(f"- collecting member data ({len(self.model.members)})...")

        points = vtk.vtkPoints()

        #### CREATE LINE CELLS ####
        # each (sub)member is further subdivided into line segments
        lines = vtk.vtkCellArray()
        submembers: List[Tuple[Tuple[float, float], Member3D, vtk.vtkLine]] = []
        for member in self.model.members.values():
            if len(member.sub_members) == 0:
                # The model has not been analyzed yet. Only add straight lines between the nodes
                line = vtk.vtkLine()
                line.SetObjectName(member.name)
                line.GetPointIds().SetId(0, points.InsertNextPoint(member.i_node.X, member.i_node.Y, member.i_node.Z))
                line.GetPointIds().SetId(1, points.InsertNextPoint(member.j_node.X, member.j_node.Y, member.j_node.Z))
                lines.InsertNextCell(line)
            else:
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

        member_writer = vtk.vtkUnstructuredGridWriter()
        member_writer.SetFileName(path)
        member_writer.SetInputData(ugrid_members)
        member_writer.Write()
        self._members_written = True
        if self.log:
            print("- members written!")

    def _write_quad_data(self, path: str):
        """
        Collects all quads in the model and writes their data into a VTK mesh. Each Pynite Quad gets subdivided in 4 vtkBiquadraticQuads,
        to avoid sampling datapoints only in the corners, where the quad results are least accurate.
        """
        if self.log:
            print(f"- collecting quad data ({len(self.model.quads)})...")

        # Define the natural coordinates xi,eta for every point-id inside a vtk_biquadratic_quad element as node_id:(xi,eta)
        # BiQuadraticQuad Element Vertex IDs
        #   3_____6_____2
        #    |   /|\   |
        #    |  / | \  |
        #    | /  |  \ |
        #    |/___|___\|
        #   7|\   |8  /|5
        #    | \  |  / |
        #    |  \ | /  |
        #    |   \|/   |
        #   0¯¯¯¯¯4¯¯¯¯¯1

        coordinates = {
            0:(-1,-1),
            1:(1,-1),
            2:(1,1),
            3:(-1,1),
            4:(0,-1),
            5:(1,0),
            6:(0,1),
            7:(-1,0),
            8:(0,0),
        }
        ugrid = vtk.vtkUnstructuredGrid() # this holds all data and gets written to .vtk at the end
        points = vtk.vtkPoints() # this holds all Point data
        quads = vtk.vtkCellArray() # this holds all cell data. A cell is in this case the vtkBiQuadraticQuads

        # keeps track of all vtk subquads for every Pynite Quad by name
        quad_references: Dict[str, List[vtk.vtkBiQuadraticQuad]] = {}

        # this loops through every Pynite quad and creates 4 vtkBiQuadraticQuad elements and their point data
        for quad in self.model.quads.values():
            quad_references[quad.name] = [vtk.vtkBiQuadraticQuad() for i in range(4)]
            # corner coordinates of the Pynite quad
            i_coords = np.array((quad.i_node.X, quad.i_node.Y, quad.i_node.Z))
            j_coords = np.array((quad.j_node.X, quad.j_node.Y, quad.j_node.Z))
            m_coords = np.array((quad.m_node.X, quad.m_node.Y, quad.m_node.Z))
            n_coords = np.array((quad.n_node.X, quad.n_node.Y, quad.n_node.Z))

            # goes through newly created subquads and calculates their vertex coordinates
            for i, subquad in enumerate(quad_references[quad.name]):
                for vert_id in range(9):
                    # this calculates the normalized coordinates (xi, eta) corresponding to the vertex by first looking up
                    # the natural coordinates in the whole grid and scaling/shifting them to the correct sector
                    xi  = coordinates[vert_id][0] * 0.5 + (i % 2) - 0.5  # Scale and shift xi
                    eta = coordinates[vert_id][1] * 0.5 + (i // 2) - 0.5 # Scale and shift eta
                    point_coords = self._interpolate_quad_corner_data(i_coords, j_coords, m_coords, n_coords, xi, eta)

                    # set the Vertex ID of the subquad to the global ID of the newly created point at the calculated coordinates
                    subquad.GetPointIds().SetId(vert_id, points.InsertNextPoint(point_coords)) # type: ignore
                subquad.SetObjectName(quad.name)
                # insert the quad into the cell array
                quads.InsertNextCell(subquad)

        #### READ QUAD DATA ####
        # After creating all the subquads, it is now necessary to fetch all relevant quad data from the Pynite model.
        # This is done for every vertex of every subquad, thus making this loop the bottleneck of the entire VTKWriter.
        # It needed to be heavily optimized to get the computation time to a reasonable level (~5s for 1000 Pynite Quads)
        for combo in self.model.load_combos.keys():
            # firstly, instantiate all data arrays for every load combo
            # Displacement Data
            D = vtk.vtkDoubleArray()
            D.SetName(f"Displacement - {combo}")
            D.SetNumberOfComponents(3)

            # Membrane Data
            membrane_loc = vtk.vtkDoubleArray()
            membrane_loc.SetName(f"Membrane Stresses sigma - {combo}")
            membrane_loc.SetNumberOfComponents(3)
            
            membrane_glob = vtk.vtkDoubleArray()
            membrane_glob.SetName(f"Membrane Stresses Sigma - {combo}")
            membrane_glob.SetNumberOfComponents(3)

            # Moment Data
            moments_loc = vtk.vtkDoubleArray()
            moments_loc.SetName(f"Moments m - {combo}")
            moments_loc.SetNumberOfComponents(3)

            moments_glob = vtk.vtkDoubleArray()
            moments_glob.SetName(f"Moments M - {combo}")
            moments_glob.SetNumberOfComponents(3)

            # SHEAR Data
            shear_loc = vtk.vtkDoubleArray()
            shear_loc.SetName(f"Forces f - {combo}")
            shear_loc.SetNumberOfComponents(3)

            shear_glob = vtk.vtkDoubleArray()
            shear_glob.SetName(f"Forces F - {combo}")
            shear_glob.SetNumberOfComponents(3)

            # go trough every pynite quad
            for quad_name, subquads in quad_references.items():
                quad = self.model.quads[quad_name]
                xi_range = np.linspace(-1,1,50)
                eta_range = np.linspace(-1,1,50)
                XI, ETA = np.meshgrid(xi_range,eta_range) # 2D Coordinate Grids
                T = quad.T()[:3,:3] # transformation matrix of the quad

                # Collect data inside the Pynite Quad at a bunch of points at once
                # to avoid calling quad.moments, .shear and .membrane inside the hotloop
                # for every points inside every subquad element. This data is then
                # interpolated to the requested vertex coordinate.
                MO = quad.moment(XI, ETA, True, combo) # type: ignore
                S = quad.shear(XI, ETA, True, combo) # type: ignore
                MB = quad.membrane(XI, ETA, True, combo) # type: ignore
                
                # Create Scipy interpolators to speed this up. Points are given reversed as (eta,xi) because
                # a transpose was needed for the correct dimensions 
                MO_interp = RegularGridInterpolator((xi_range, eta_range), MO.T)
                S_interp  = RegularGridInterpolator((xi_range, eta_range), S.T)
                MB_interp = RegularGridInterpolator((xi_range, eta_range), MB.T)

                # go through every subquad of the Pynite Quad
                for i, subquad in enumerate(subquads):
                    # calculate the (eta,xi) coordinates for every vertex in advance
                    xi  = [coordinates[vert_id][0] * 0.5 + (i % 2) - 0.5 for vert_id in range(9)]
                    eta = [coordinates[vert_id][1] * 0.5 + (i // 2) - 0.5 for vert_id in range(9)]

                    # Precalculate all vertex results at once to avoid calling the interpolator for every vertex.
                    # This gets the precomputed 2D Grid data from earlier and interpolates its data to the correct position.local -> global
                    # go trough every vertex of the subquad
                    p_membranes = MB_interp((eta,xi))
                    p_moments = MO_interp((eta,xi))
                    p_shears = S_interp((eta,xi))
                    for vert_id in range(9):
                        # DISPLACEMENT
                        di = np.array([quad.i_node.DX[combo], quad.i_node.DY[combo], quad.i_node.DZ[combo]])
                        dj = np.array([quad.j_node.DX[combo], quad.j_node.DY[combo], quad.j_node.DZ[combo]])
                        dm = np.array([quad.m_node.DX[combo], quad.m_node.DY[combo], quad.m_node.DZ[combo]])
                        dn = np.array([quad.n_node.DX[combo], quad.n_node.DY[combo], quad.n_node.DZ[combo]])

                        d = self._interpolate_quad_corner_data(di, dj, dm, dn, xi[vert_id], eta[vert_id])
                        # The fetched data now gets inserted into the corresponding data array by the id of the vertex in the
                        # points array
                        D.InsertTuple3(subquad.GetPointId(vert_id), *d) # type: ignore
                        
                        p_membrane = p_membranes[vert_id]
                        p_moment = p_moments[vert_id]
                        p_shear = np.insert(p_shears[vert_id].flatten(),0,(0))

                        # transform the local quad data into global coordinates
                        p_membrane_glob = T.T @ p_membrane
                        p_moment_glob  = T.T @ p_moment
                        p_shear_glob    = T.T @ p_shear

                        membrane_loc.InsertTuple3(subquad.GetPointId(vert_id), *p_membrane) # type: ignore
                        moments_loc.InsertTuple3(subquad.GetPointId(vert_id), *p_moment) # type: ignore
                        shear_loc.InsertTuple3(subquad.GetPointId(vert_id), *p_shear) # type: ignore

                        membrane_glob.InsertTuple3(subquad.GetPointId(vert_id), *p_membrane_glob) # type: ignore
                        moments_glob.InsertTuple3(subquad.GetPointId(vert_id), *p_moment_glob) # type: ignore
                        shear_glob.InsertTuple3(subquad.GetPointId(vert_id), *p_shear_glob) # type: ignore

                ugrid.GetPointData().AddArray(D)
                ugrid.GetPointData().AddArray(membrane_loc)
                ugrid.GetPointData().AddArray(membrane_glob)
                ugrid.GetPointData().AddArray(moments_loc)
                ugrid.GetPointData().AddArray(moments_glob)
                ugrid.GetPointData().AddArray(shear_loc)
                ugrid.GetPointData().AddArray(shear_glob)


        ugrid.SetPoints(points)
        ugrid.SetCells(vtk.VTK_BIQUADRATIC_QUAD, quads)

        # clean the data from duplicate points
        cleaner = vtk.vtkStaticCleanUnstructuredGrid()
        cleaner.SetInputData(ugrid)
        cleaner.SetToleranceIsAbsolute(True)
        cleaner.SetAbsoluteTolerance(0.01)
        cleaner.Update()
        ugrid = cleaner.GetOutput()

        #### WRITE DATA TO DISK ####
        quads_writer = vtk.vtkUnstructuredGridWriter()
        quads_writer.SetFileName(path)
        quads_writer.SetInputData(ugrid)
        quads_writer.Write()
        self._quads_written = True
        if self.log:
            print("- quads written!")

    @staticmethod
    def _interpolate_member_data(p1:np.ndarray,p2:np.ndarray, x:float) -> np.ndarray:
        """
        Takes in vectors of size n from both ends of the member and interpolates linearly to relative position x on the member [0-1].
        This can be used i.e. to calculate every point position on a member given by a relative coordinate.
        """
        return x*p1 + (1-x)*p2

    @staticmethod
    def _interpolate_quad_corner_data(i:float|np.ndarray, j:float|np.ndarray, m:float|np.ndarray, n:float|np.ndarray, xi:float, eta:float) -> float|np.ndarray:
        """
        Helper Method to return the linearly interpolated data of a point on the quad, given the natural coordinates
        xi and eta.
        """
        # Bilinear shape functions for natural coordinates [-1, 1]
        N1 = 0.25 * (1 - xi) * (1 - eta)  # Shape function for corner i @ (-1, -1)
        N2 = 0.25 * (1 + xi) * (1 - eta)  # Shape function for corner j @ ( 1, -1)
        N3 = 0.25 * (1 + xi) * (1 + eta)  # Shape function for corner m @ ( 1,  1)
        N4 = 0.25 * (1 - xi) * (1 + eta)  # Shape function for corner n @ (-1,  1)

        result_vector = N1 * i + N2 * j + N3 * m + N4 * n

        return result_vector

    def open_in_paraview(self):
        """
        Open the Model inside Paraview, if installed. Paraview must be accessible with the `paraview` command for this method to work.
        """
        # Create a temporary file
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
