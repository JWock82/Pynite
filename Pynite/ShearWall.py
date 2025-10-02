from __future__ import annotations # Allows more recent type hints features
from math import isclose
from numpy import average
import io
from typing import TYPE_CHECKING, Literal

from prettytable import PrettyTable

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

if TYPE_CHECKING:
    from typing import List, Dict, Tuple
    from Pynite.Quad3D import Quad3D
    import matplotlib.figure


class ShearWall():
    """Creates a new shear wall model that allows for modeling of complex shear walls. You can add openings and flanges (wall returns or intersections). Diaphragm levels can be defined in order to apply shear forces along the length of the wall. Diaphragms can be full or partial length. Supports can be applied at any level in the shear wall. Supports can also be full or partial length. A `ky_mod` factor is built in to account for cracking. Shear walls can automatically detect shear wall piers and coupling beams, and sum internal forces in those components.
    """

    def __init__(self, model, name, mesh_size, length, height, thickness, material_name, ky_mod=0.35, origin=[0, 0, 0], plane='XY') -> None:

        self.model = model
        self.name = name
        self.mesh_size = mesh_size
        self.L = length
        self.H = height
        self.ky_mod = ky_mod
        self.origin = origin
        self.plane = plane
        self._openings: List[List[str | float | None]] = []
        self._flanges: List[List[str | float]] = []
        self._supports: List[List[float]] = []
        self._stories: List[List[str | float]] = []
        self._shears: List[List[str | float]] = []
        self._axials: List[List[str | float]] = []
        self._materials: List[List[str | float]] = []
        self.piers: Dict[str, Pier] = {}
        self.coupling_beams: Dict[str, CouplingBeam] = {}
        self.asign_material(material_name, thickness)

    def asign_material(self, name: str, t: float, x_start: float | None = None, x_end: float | None = None, y_start: float | None = None, y_end: float | None = None) -> None:

        # Data validation: ensure the extents of the material fall within the wall's envelope
        if x_start is None: x_start = 0
        if x_end is None: x_end = self.L
        if y_start is None: y_start = 0
        if y_end is None: y_end = self.H

        # Store the material parameters
        self._materials.append([name, t, x_start, x_end, y_start, y_end])

    def add_opening(self, name: str, x_start: float, y_start: float, width: float, height: float, tie: float | None = None) -> None:
        self._openings.append([name, x_start, y_start, width, height, None])

    def add_flange(self, thickness: float, width: float, x: float, y_start: float, y_end: float, material: str, side: Literal['+z', '-z']) -> None:
        self._flanges.append([thickness, width, x, y_start, y_end, material, side])

    def add_support(self, elevation: float | None = None, x_start: float | None = None, x_end: float | None = None) -> None:

        # Default support settings
        if elevation is None: elevation = 0
        if x_start is None: x_start = 0
        if x_end is None: x_end = self.L

        # Handle invalid supports
        if elevation < 0 or elevation > self.H or x_start < 0 or x_start > self.L or x_end < 0 or x_end > self.L:
            raise Exception(f'Support for shear wall `{self.name}` falls outside the extents of the wall and is invalid.')

        # Add the supports to the wall
        self._supports.append([elevation, x_start, x_end])

    def add_story(self, story_name: str, elevation: float, x_start: float | None = None, x_end: float | None = None) -> None:

        # Validate input
        if elevation is None: elevation = self.H
        if x_start is None: x_start = 0
        if x_end is None: x_end = self.L

        # Add the story to the model
        self._stories.append([story_name, elevation, x_start, x_end])

        # Add a load combination to use when calculating the story's stiffness
        self.model.add_load_combo('Stiffness: ' + story_name, {story_name: 1.0}, 'stiffness')

        # Add a 100 kip story shear to the model to use when calculating the story's stiffness
        self.add_shear(story_name, 100, case=story_name)

    def add_shear(self, story_name: str, force: float, case: str = 'Case 1') -> None:
        self._shears.append([story_name, force, case])

    def add_axial(self, story_name: str, force: float, case: str = 'Case 1') -> None:
        self._axials.append([story_name, force, case])

    def generate(self) -> None:

        # Identify mesh control points
        x_control: List[float] = [0, self.L]
        y_control: List[float] = [0, self.H]

        for material in self._materials:
            x_control.append(material[2])
            x_control.append(material[3])
            y_control.append(material[4])
            y_control.append(material[5])

        z_control: List[float] = [0]
        for flg in self._flanges:
            if flg[6] == '+z': z_control.append(flg[1])
            else: z_control.append(-flg[1])
            x_control.append(flg[2])
            y_control.append(flg[3])
            y_control.append(flg[4])

        for support in self._supports:
            x_control.append(support[1])
            x_control.append(support[2])
            y_control.append(support[0])

        for story in self._stories:
            x_control.append(story[2])
            x_control.append(story[3])
            y_control.append(story[1])

        # While opening control points are auto-generated by the wall's mesh, we have no way of generating them for the flange meshes. We'll add some control points for the sake of the flanges. Duplicate control point values in the wall be be automatically resolved by Pynite.
        for opng in self._openings:
            y_control.append(opng[2])
            y_control.append(opng[2] + opng[4])

        # Add the wall mesh to the model
        # Note that the thickness can still be changed later
        material_name, thickness, _, _, _, _ = self._materials[0]
        self.model.add_rectangle_mesh(self.name, self.mesh_size, self.L, self.H, 1, material_name, thickness, self.ky_mod, self.origin, self.plane, x_control=x_control, y_control=y_control)

        # Add a tie material if it's not already in the model
        if 'Tie' not in self.model.materials.keys():   
            self.model.add_material('Tie', 1, 1, 0, 0)

        # Add the openings to the mesh
        for opng in self._openings:

            name, x_start, y_start, width, height, AE = opng
            self.model.meshes[self.name].add_rect_opening(name, x_start, y_start, width, height)

            # Add any ties over the opening
            if AE is not None:

                i_node_name = self.model.unique_name(self.model.nodes, 'N')
                self.model.add_node(i_node_name, *self._local2global(x_start, y_start + height, 0, self.plane))

                j_node_name = self.model.unique_name(self.model.nodes, 'N')
                self.model.add_node(j_node_name, *self._local2global(x_start + width, y_start + height, 0, self.plane))

                tie_name = self.model.unique_name(self.model.Members, 'Tie ')
                self.model.add_member(tie_name, i_node_name, j_node_name, 'Tie', 1, 1, 1, AE)
                self.model.def_releases(tie_name, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1)

        # Add the flanges to the mesh
        for i, flg  in enumerate(self._flanges):

            # Read in the flange's parameters
            t, b, x, y_start, y_end, material, side = flg

            # Determine which side of the wall to place the flange on and define control points for the flange mesh so nodes line up properly with other meshes
            if side == '+z':
                z = 0
                flg_x_control = [val for val in z_control if round(val, 10) >= 0 and round(val, 10) <= b]
            else:
                z = -b
                flg_x_control = [b - (-val) for val in z_control if round(val, 10) <= 0 and round(val, 10) >= -b]

            flg_y_control = [y - y_start for y in y_control if round(y, 10) >= round(y_start, 10) and round(y, 10) <= round(y_end, 10)]

            # Identify the global origin for the flange mesh
            if self.plane == 'XY':
                Xof = self.origin[0] + x
                Yof = self.origin[0] + y_start
                Zof = self.origin[2] + z
                flg_plane = 'YZ'
            elif self.plane == 'XZ':
                Xof = self.origin[0] + x
                Yof = self.origin[1] + z
                Zof = self.origin[2] + y_start
                flg_plane = 'YZ'
            elif self.plane == 'YZ':
                Xof = self.origin[0] + z
                Yof = self.origin[1] + y_start
                Zof = self.origin[2] + x
                flg_plane = 'XY'

            self.model.add_rectangle_mesh(self.name + ' Flg ' + str(i+1), self.mesh_size, b, y_end-y_start, t, material, 1, self.ky_mod, [Xof, Yof, Zof], flg_plane, flg_x_control, flg_y_control)

        # Generate the meshes
        self.model.meshes[self.name].generate()

        for i, flg in enumerate(self._flanges):
            self.model.meshes[self.name + ' Flg ' + str(i + 1)].generate()

        # Merge the flange nodes with the rest of the wall
        self.model.merge_duplicate_nodes()

        # Step through each plate in the model
        for plate in self.model.quads.values():

            # Step through each material in the wall
            for material in self._materials:

                # Get the material properties
                material_name, t, x_start, x_end, y_start, y_end = material

                xi, yi, zi = _global2local(plate.i_node.X, plate.i_node.Y, plate.i_node.Z, self.origin, self.plane)
                xj, yj, zj = _global2local(plate.j_node.X, plate.j_node.Y, plate.j_node.Z, self.origin, self.plane)
                xm, ym, zm = _global2local(plate.m_node.X, plate.m_node.Y, plate.m_node.Z, self.origin, self.plane)

                # Determine if the current plate is part of a flange
                if isclose(xi, xj):
                    # Flanges already have material properties and thicknesses assigned properly
                    pass
                else:
                    # Determine if the current plate is this material
                    if (round(xi, 10) >= round(x_start, 10)
                    and round(xm, 10) <= round(x_end, 10)
                    and round(yi, 10) >= round(y_start, 10)
                    and round(ym, 10) <= round(y_end, 10)):

                        # Assign material properties to the plate
                        plate.E = self.model.materials[material_name].E
                        plate.nu = self.model.materials[material_name].nu
                        plate.t = t

        # Add supports
        for support in self._supports:
            elevation, x_start, x_end = support
            for node in self.model.nodes.values():

                x, y, z = _global2local(node.X, node.Y, node.Z, self.origin, self.plane)

                if isclose(y, elevation) and round(x, 10) >= round(x_start, 10) and round(x, 10) <= round(x_end, 10):
                    self.model.def_support(node.name, True, True, True, True, True, True)

        # Add shear forces to the wall
        for story in self._stories:

            # Read in parameters for this story
            story_name, elevation, x_start, x_end = story

            # Initialize a list of story nodes
            node_list = []

            # Step through each node in the model
            for node in self.model.nodes.values():

                x, y, z = _global2local(node.X, node.Y, node.Z, self.origin, self.plane)

                # Check if this node belongs to this story
                if isclose(y, elevation) and x >= x_start and x <= x_end and isclose(z, 0):

                    # Add the node to the list of nodes in the current story
                    node_list.append(node)

            # Add shear and axial forces to all the nodes in the story
            for node in node_list:

                # Determine how many nodes are in the current story
                num_nodes = len(node_list)

                # Step through each shear force in the model
                for shear in self._shears:

                    # Read in parameters for this shear
                    story, force, case = shear

                    # Determine the global direction to apply the shear in
                    if self.plane == 'XY': direction = 'FX'
                    elif self.plane == 'YZ': direction = 'FZ'

                    # Determine if this shear acts on this story
                    if story == story_name:
                        self.model.add_node_load(node.name, direction, force/num_nodes, case)

                # Step through each axial force in the model
                for axial in self._axials:

                    # Read in parameters for this axial force
                    story, force, case = axial

                    # Determine if this axial force acts on this story
                    if story == story_name:
                        self.model.add_node_load(node.name, 'FY', -force/num_nodes, case)

        # Populate dictionaries of piers and coupling beams for the wall
        self._identify_piers()
        self._identify_coupling_beams()


    def _identify_piers(self) -> None:

        # Reset all piers in the wall
        self.piers = {}

        # Create a list of x and y coordinates that represent the edges of the wall
        x_vals: List[float] = [0, self.L]
        y_vals: List[float] = [0, self.H]

        # Add the edges of the openings to the lists
        for opng in self._openings:
            x_vals.append(opng[1])
            x_vals.append(opng[1] + opng[3])
            y_vals.append(opng[2])
            y_vals.append(opng[2] + opng[4])

        # Sort the lists (ascending)
        x_vals = sorted(x_vals)
        y_vals = sorted(y_vals)

        # Remove duplicate (or near duplicate) values
        unique_list: List[float] = []
        for i in range(len(x_vals) - 1):
            # Only keep the value at `i` if it's not a duplicate or near duplicate of the next value
            if not isclose(x_vals[i], x_vals[i+1]):
                unique_list.append(x_vals[i])
        unique_list.append(x_vals[-1])  # The last value will always be a keeper
        x_vals = unique_list

        unique_list = []
        for i in range(len(y_vals) - 1):
            # Only keep the value at `i` if it's not a duplicate or near duplicate of the next value
            if not isclose(y_vals[i], y_vals[i+1]):
                unique_list.append(y_vals[i])
        unique_list.append(y_vals[-1])  # The last value will always be a keeper
        y_vals = unique_list

        # Divide the wall into vertical strip piers using the left and right edges of each opening as strip boundaries
        self.piers = {}
        for i in range(len(x_vals) - 1):
            width = x_vals[i+1] - x_vals[i]
            height = self.H
            x = x_vals[i]
            y = 0
            self.piers['P' + str(i+1)] = Pier('P' + str(i+1), x, y, width, height, self)

        # Divide the strip piers further into rectanglular piers using the top and bottom of each opening as pier boundaries
        new_piers: Dict[str, Pier] = {}
        pier_count = 1
        for pier in self.piers.values():
            for i in range(len(y_vals) - 1):
                width = pier.width
                height = y_vals[i + 1] - y_vals[i]
                x = pier.x
                y = y_vals[i]
                new_piers['P' + str(pier_count)] = Pier('P' + str(pier_count), x, y, width, height, self)
                pier_count += 1
        self.piers = new_piers

        # Delete any piers that fall within an opening
        delete_list: List[str] = []
        for pier in self.piers.values():

            # Check if this pier is inside any of the openings
            for opng in self._openings:
                if (round(pier.x, 10) >= round(opng[1], 10)
                   and round(pier.x + pier.width, 10) <= round(opng[1] + opng[3], 10)
                   and round(pier.y, 10) >= round(opng[2], 10)
                   and round(pier.y + pier.height, 10) <= round(opng[2] + opng[4], 10)):
                    delete_list.append(pier.name)
                    break

        for pier in delete_list:
            del self.piers[pier]

        # Working horizontally (left to right), rejoin any rectangles that share a vertical edge to form a larger rectangle
        found_duplicate = True
        while found_duplicate == True:

            found_duplicate = False
            piers_copy = self.piers.copy()

            for key1, pier1 in piers_copy.items():

                for key2, pier2 in piers_copy.items():

                    # Check for piers that need to be merged
                    if (key1 != key2
                        and isclose(pier1.y, pier2.y)
                        and isclose(pier1.x + pier1.width, pier2.x)
                        and isclose(pier1.height, pier2.height)):

                        # Merge the piers in the `self.piers` dictionary
                        self.piers[key1].width = pier1.width + pier2.width

                        # Delete the 2nd pier from the `self.piers` dictionary
                        del self.piers[key2]

                        # Since the `self.piers` dictionary has changed we need `piers_copy` to get updated. Flag that we found a duplicate and break the loops.
                        found_duplicate = True
                        break

                # Break the `for` loop if a duplicate was found so we can get an updated copy of `self.piers`
                if found_duplicate == True:
                    break

        # Working vertically (bottom to top), rejoin any rectangles that share a horizontal edge to form a larger rectangle
        found_duplicate = True
        while found_duplicate == True:

            found_duplicate = False
            piers_copy = self.piers.copy()

            for key1, pier1 in piers_copy.items():

                for key2, pier2 in piers_copy.items():

                    if (key1 != key2
                        and isclose(pier1.x, pier2.x)
                        and isclose(pier1.y + pier1.height, pier2.y)
                        and isclose(pier1.width, pier2.width)):

                        # Merge the piers in the `self.piers` dictionary
                        self.piers[key1].height = pier1.height + pier2.height

                        # Delete the 2nd pier from the `self.piers` dictionary
                        del self.piers[key2]

                        # Since the `self.piers` dictionary has changed we need `piers_copy` to get updated. Flag that we found a duplicate and break the loops.
                        found_duplicate = True
                        break

                # Break the `for` loop if a duplicate was found so we can get an updated copy of `self.piers`
                if found_duplicate == True:
                    break

        # Generate a list of new keys in ascending order
        new_keys: List[str] = [f'P{i+1}' for i in range(len(self.piers))]

        # Replace the old dicionary with one that has updated keys
        self.piers = dict(zip(new_keys, self.piers.values()))
        for key, pier in self.piers.items():
            pier.name = key

        # Assign plates to each pier
        for plate in self.model.quads.values():

            xi, yi, zi = _global2local(plate.i_node.X, plate.i_node.Y, plate.i_node.Z, self.origin, self.plane)
            xj, yj, zj = _global2local(plate.j_node.X, plate.j_node.Y, plate.j_node.Z, self.origin, self.plane)
            xm, ym, zm = _global2local(plate.m_node.X, plate.m_node.Y, plate.m_node.Z, self.origin, self.plane)

            y_avg = (yi + ym)/2
            x_avg = (xi + xm)/2

            for pier in self.piers.values():

                if (round(x_avg, 10) >= round(pier.x, 10)
                    and round(x_avg, 10) <= round(pier.x + pier.width, 10)
                    and round(y_avg, 10) >= round(pier.y, 10)
                    and round(y_avg, 10) <= round (pier.y + pier.height, 10)):

                    pier.plates.append(plate)

    def _identify_coupling_beams(self) -> None:
        
        # Reset all coupling beams in the wall
        self.coupling_beams = {}

        # Create a list of x and y coordinates that represent the edges of the wall
        x_vals: List[float] = [0, self.L]
        y_vals: List[float] = [0, self.H]

        # Add the edges of the openings to the lists
        for opng in self._openings:
            x_vals.append(opng[1])
            x_vals.append(opng[1] + opng[3])
            y_vals.append(opng[2])
            y_vals.append(opng[2] + opng[4])

        # Sort the lists (ascending)
        x_vals = sorted(x_vals)
        y_vals = sorted(y_vals)

        # Remove duplicate (or near duplicate) values
        unique_list: List[float] = []
        for i in range(len(x_vals) - 1):
            # Only keep the value at `i` if it's not a duplicate or near duplicate of the next value
            if not isclose(x_vals[i], x_vals[i+1]):
                unique_list.append(x_vals[i])
        unique_list.append(x_vals[-1])  # The last value will always be a keeper
        x_vals = unique_list

        unique_list = []
        for i in range(len(y_vals) - 1):
            # Only keep the value at `i` if it's not a duplicate or near duplicate of the next value
            if not isclose(y_vals[i], y_vals[i+1]):
                unique_list.append(y_vals[i])
        unique_list.append(y_vals[-1])  # The last value will always be a keeper
        y_vals = unique_list

        # Divide the wall into horizontal strips using the bottom and top edges of each opening as strip boundaries
        self.coupling_beams = {}
        for i in range(len(y_vals) - 1):
            height = y_vals[i + 1] - y_vals[i]
            length = self.L
            y = y_vals[i]
            x = 0
            self.coupling_beams['B' + str(i+1)] = CouplingBeam('B' + str(i+1), x, y, length, height, self)

        # Divide the strips further into rectanglular beams using the left and right of each opening as beam boundaries
        new_beams: Dict[str, CouplingBeam] = {}
        beam_count = 1
        for beam in self.coupling_beams.values():
            for i in range(len(x_vals) - 1):
                height = beam.height
                length = x_vals[i+1] - x_vals[i]
                y = beam.y
                x = x_vals[i]
                new_beams['B' + str(beam_count)] = CouplingBeam('B' + str(beam_count), x, y, length, height, self)
                beam_count += 1
        self.coupling_beams = new_beams

        # Delete any beams that fall within an opening
        delete_list: List[str] = []
        for beam in self.coupling_beams.values():
           
           # Check if this beam is inside any of the openings
           for opng in self._openings:
               
               if (round(beam.x, 10) >= round(opng[1], 10)
                   and round(beam.x + beam.length, 10) <= round(opng[1] + opng[3], 10)
                   and round(beam.y, 10) >= round(opng[2], 10)
                   and round(beam.y + beam.height, 10) <= round(opng[2] + opng[4], 10)):
                    delete_list.append(beam.name)
                    break
        
        for beam in delete_list:
            del self.coupling_beams[beam]
        
        # Working vertically (bottom to top), rejoin any rectangles that share a horizontal edge to form a larger rectangle
        found_duplicate = True
        while found_duplicate == True:

            found_duplicate = False
            beams_copy = self.coupling_beams.copy()

            for key1, beam1 in beams_copy.items():

                for key2, beam2 in beams_copy.items():

                    # Check for beams that need to be merged
                    if (key1 != key2
                        and isclose(beam1.x, beam2.x)
                        and isclose(beam1.y + beam1.height, beam2.y)
                        and isclose(beam1.length, beam2.length)):

                        # Merge the beams in the `self.coupling_beams` dictionary
                        self.coupling_beams[key1].height = beam1.height + beam2.height

                        # Delete the 2nd beam from the `self.coupling_beams` dictionary
                        del self.coupling_beams[key2]

                        # Since the `self.coupling_beams` dictionary has changed we need `beams_copy` to get updated. Flag that we found a duplicate and break the loops.
                        found_duplicate = True
                        break

                # Break the `for` loop if a duplicate was found so we can get an updated copy of `self.coupling_beams`
                if found_duplicate == True:
                    break

        # Working horizontally (left to right), rejoin any rectangles that share a vertical edge to form a larger rectangle
        found_duplicate = True
        while found_duplicate == True:
            
            found_duplicate = False
            beams_copy = self.coupling_beams.copy()

            for key1, beam1 in beams_copy.items():

                for key2, beam2 in beams_copy.items():

                    if (key1 != key2
                        and isclose(beam1.y, beam2.y)
                        and isclose(beam1.x + beam1.length, beam2.x)
                        and isclose(beam1.height, beam2.height)):

                        # Merge the beams in the `self.coupling_beams` dictionary
                        self.coupling_beams[key1].length = beam1.length + beam2.length

                        # Delete the 2nd beam from the `self.coupling_beams` dictionary
                        del self.coupling_beams[key2]

                        # Since the `self.couping_beams` dictionary has changed we need `beams_copy` to get updated. Flag that we found a duplicate and break the loops.
                        found_duplicate = True
                        break

                # Break the `for` loop if a duplicate was found so we can get an updated copy of `self.coupling_beams`
                if found_duplicate == True:
                    break

        # Check for any coupling beams at the bottom of the wall. There should not be any. Delete them as they are found
        delete_list = []
        for beam in self.coupling_beams.values():
            if beam.y == 0:
                delete_list.append(beam.name)
        
        for beam in delete_list:
            del self.coupling_beams[beam]

        # Generate a list of new keys in ascending order
        new_keys: List[str] = [f'B{i + 1}' for i in range(len(self.coupling_beams))]

        # Replace the old dicionary with one that has updated keys
        self.coupling_beams = dict(zip(new_keys, self.coupling_beams.values()))
        for key, beam in self.coupling_beams.items():
            beam.name = key

        # Assign plates to each beam
        for plate in self.model.quads.values():

            xi, yi, zi = _global2local(plate.i_node.X, plate.i_node.Y, plate.i_node.Z, self.origin, self.plane)
            xm, ym, zm = _global2local(plate.m_node.X, plate.m_node.Y, plate.m_node.Z, self.origin, self.plane)

            y_avg = (yi + ym)/2
            x_avg = (xi + xm)/2

            for beam in self.coupling_beams.values():

                if (round(x_avg, 10) >= round(beam.x, 10)
                    and round(x_avg, 10) <= round(beam.x + beam.length, 10)
                    and round(y_avg, 10) >= round(beam.y, 10)
                    and round(y_avg, 10) <= round(beam.y + beam.height, 10)):

                    beam.plates.append(plate)

    def draw_piers(self, show: bool = False) -> None | matplotlib.figure.Figure:
        
        fig, ax = plt.subplots()

        ax.patch.set_facecolor((0.8, 0.8, 0.8))

        for pier in self.piers.values():
            self._add_rectangle(ax, pier.x, pier.y, pier.width, pier.height, pier.name)
        
        # Adjust the aspect ratio of the plot
        ax.set_aspect('equal')

        # Slim down the margins
        plt.tight_layout()

        # show plot or return it
        if show == True: plt.show()
        else: return plt

    def draw_coupling_beams(self, show: bool = False) -> None | matplotlib.figure.Figure:
        
        fig, ax = plt.subplots()

        ax.patch.set_facecolor((0.8, 0.8, 0.8))

        # Draw the overall Wall
        self._add_rectangle(ax, 0, 0, self.L, self.H, "", 'white')

        # Draw the openings
        for opng in self._openings:
            self._add_rectangle(ax, opng[1], opng[2], opng[3], opng[4], '', 'grey')

        for beam in self.coupling_beams.values():
            self._add_rectangle(ax, beam.x, beam.y, beam.length, beam.height, beam.name, 'white')
        
        # Adjust the aspect ratio of the plot
        ax.set_aspect('equal')

        # Slim down the margins
        plt.tight_layout()

        # show plot or return it
        if show == True: plt.show()
        else: return plt

    def _add_rectangle(self, ax: matplotlib.axes.Axes, x: float, y: float, w: float, h: float, name: str, color: str = 'white') -> None:
        """Adds a rectangle to the pyplot
        """

        # create rectangle
        rect = Rectangle((x, y), w, h, linewidth=1, edgecolor='r', facecolor=color)
        ax.add_patch(rect)

        # add name to center of rectangle
        ax.text(x + w/2, y + h/2, name, ha='center', va='center')

        # set plot limits
        ax.set_xlim(0, max(ax.get_xlim()[1], x + w))
        ax.set_ylim(0, max(ax.get_ylim()[1], y + h))

    def _sort_openings(self) -> None:

        # Sort the openings based on y-coordinates
        n = len(self._openings)
        for i in range(n):
            for j in range(0, n-i-1):
                if self._openings[j][2] > self._openings[j+1][2]:
                    self._openings[j], self._openings[j+1] = self._openings[j+1], self._openings[j]
        
        # Sort the openings based on x-coordinates
        n = len(self._openings)
        for i in range(n):
            for j in range(0, n-i-1):
                if self._openings[j][1] > self._openings[j+1][1]:
                    self._openings[j], self._openings[j+1] = self._openings[j+1], self._openings[j]

    def stiffness(self, story_name: str) -> float:

        # Validate that the specified story exists in the shear wall
        if not any(story[0] == story_name for story in self._stories):
            raise Exception(f'Unable to obtain stiffness. Story {story_name} does not exist.')

        # Step through each story in the model to find the one we're looking for
        for story in self._stories:

            # Determine if this story is the one we are interested in
            if story[0] == story_name:

                # Exit the loop
                break

        # 100 kips is being applied to the story for the purpose of determining stiffness
        V = 100

        # Initialize the maximum wall deflection to zero
        d_max = 0

        # Step through each node in the model
        for node in self.model.nodes.values():

            local_coords = _global2local(node.X, node.Y, node.Z, self.origin, self.plane)
            x, y, z = local_coords[0], local_coords[1], local_coords[2]

            # Determine if this node is in this story
            if (round(x, 10) >= round(story[2], 10)
            and round(x, 10) <= round(story[3], 10)
            and isclose(story[1], y)
            and isclose(z, 0)):

                # Check if this deflection is the largest in the story
                if node.DX['Stiffness: ' + story_name] > d_max:

                    d_max = node.DX['Stiffness: ' + story_name]

        # Return the story's stiffness:
        return V/(d_max)

    def screenshots(self, combo_name: str = 'Combo 1', dir_path: str = './') -> None:

        from Pynite.Rendering import Renderer

        renderer = Renderer(self.model)
        renderer.window_width = 750
        renderer.window_height = 750
        renderer.annotation_size = self.mesh_size/6
        renderer.deformed_shape = True
        renderer.deformed_scale = 400
        renderer.render_loads = True
        renderer.scalar_bar = True
        renderer.combo_name = combo_name
        renderer.labels = False
        
        # Save the shear plot screenshot to this file's directory
        renderer.color_map = 'Txy'
        renderer.screenshot(dir_path + '/shear_wall_screenshot1.png', interact=True)
        
        # Save the shear plot screenshot to this file's directory
        renderer.color_map = 'Sy'
        renderer.screenshot(dir_path + '/shear_wall_screenshot2.png', interact=False, reset_camera=False)

        # Save the pier screenshot to this file's directory
        pier_sketch = self.draw_piers(show=False)
        pier_sketch.savefig(dir_path + '/shear_wall_piers.png', format='png')

    def print_piers(self, combo_name: str = 'Combo 1') -> None:
        """Tabulates and prints pier results for the shear wall
        """

        # Create a PrettyTable object
        table = PrettyTable()

        # Define the headers
        table.field_names = ["ID", "Length", "Height", "M/(VL)", "V", "M", "P"]

        # Add rows to the table
        for pier_id, pier in self.piers.items():
            P, M, V, M_VL = pier.sum_forces(combo_name)
            table.add_row([pier.name, pier.width, pier.height, M_VL, V, M, P])

        # Print the table
        print('+-------------------+')
        print('| Wall Pier Results |')
        print('+-------------------+')
        print(table)

    def print_coupling_beams(self, combo_name: str = 'Combo 1') -> None:
        """Tabulates and prints coupling beam results for the shear wall
        """

        # Create a PrettyTable object
        table = PrettyTable()

        # Define the headers
        table.field_names = ["ID", "Length", "Height", "M/(VH)", "V", "M", "P"]

        # Add rows to the table
        for beam_id, beam in self.coupling_beams.items():
            P, M, V, M_VL = beam.sum_forces(combo_name)
            table.add_row([beam.name, beam.length, beam.height, M_VL, V, M, P])

        # Print the table
        print('+----------------------------+')
        print('| Wall Coupling Beam Results |')
        print('+----------------------------+')
        print(table)

    def _local2global(self, x: float, y: float, z: float, plane: Literal['XY', 'XZ', 'YZ'] = 'XY') -> List[float]:

        Xo, Yo, Zo = self.origin[0], self.origin[1], self.origin[2]

        if plane == 'XY':
            X = Xo + x
            Y = Yo + y
            Z = Zo
        elif plane == 'XZ':
            X = Xo + x
            Y = Yo
            Z = Zo + y
        elif plane == 'YZ':
            X = Xo
            Y = Yo + y
            Z = Zo + x

        return [X, Y, Z]


def _global2local(X: float, Y: float, Z: float, origin: List[float] = [0, 0, 0], plane: Literal['XY', 'YZ', 'XZ'] = 'XY') -> List[float]:

    Xo, Yo, Zo = origin[0], origin[1], origin[2]

    if plane == 'XY':
        x = X - Xo
        y = Y - Yo
        z = Z - Zo
    elif plane == 'YZ':
        x = Z - Zo
        y = Y - Yo
        z = -(X - Xo)
    elif plane == 'XZ':
        x = X - Xo
        y = Z - Zo
        z = -(Y - Yo)
    else:
        raise Exception('Invalid plane specified for coordinate transformation')

    return [x, y, z]


# %%
class Pier():

    def __init__(self, name: str, x: float, y: float, width: float, height: float, shear_wall: ShearWall) -> None:
        self.name: str = name
        self.x: float = x  # The location of the left side of the pier
        self.y: float = y  # The height of the bottom of the pier
        self.width: float = width
        self.height: float = height

        # Read in relevant data from the parent shear wall
        self.shear_wall = shear_wall
        self.plane = shear_wall.plane
        self.origin = shear_wall.origin

        # This list will be used by the parent shear wall to store a list of only the plates in this pier
        self.plates: List[Quad3D] = []

    def sum_forces(self, combo_name: str = 'Combo 1') -> Tuple[float, float, float, float]:

        # Initialize the forces in the plate
        P, M, V = 0, 0, 0

        # Step through each plate in the pier
        for plate in self.plates:

            xi, yi, zi = _global2local(plate.i_node.X, plate.i_node.Y, plate.i_node.Z, self.origin, self.plane)
            xj, yj, zj = _global2local(plate.j_node.X, plate.j_node.Y, plate.j_node.Z, self.origin, self.plane)

            # Determine if this plate is at the bottom of the pier
            if isclose(yi, self.y):

                # Find and sum the axial forces in this plate
                Pi = plate.f(combo_name)[1][0]
                Pj = plate.f(combo_name)[7][0]
                P += -Pi - Pj

                # Find and sum the moments about the pier's center in this plate
                xi_p = xi - (self.x + self.width/2)
                xj_p = xj - (self.x + self.width/2)
                Mi = plate.f(combo_name)[1][0]*xi_p
                Mj = plate.f(combo_name)[7][0]*xj_p
                M += -Mi - Mj

                # Find and sum the shear forces in this plate
                # Check if this is a flange plate or a web plate
                if isclose(xi, xj):
                    Vi = -plate.f(combo_name)[2][0]
                    Vj = -plate.f(combo_name)[8][0]
                else:
                    Vi = -plate.f(combo_name)[0][0]
                    Vj = -plate.f(combo_name)[6][0]
                V += -Vi - Vj

        # Calculate the shear span ratio
        M_VL = M/(V*self.width)

        # Return the summed forces and shear span ratio
        return P, M, V, M_VL


# %%
class CouplingBeam():

    def __init__(self, name: str, x: float, y: float, length: float, height: float, shear_wall: ShearWall) -> None:
        self.name: str = name
        self.x: float = x  # The location of the left side of the coupling beam
        self.y: float = y  # The height to the bottom of the coupling beam
        self.length: float = length
        self.height: float = height

        # Read in relevant data from the parent shear wall
        self.shear_wall = shear_wall
        self.plane = shear_wall.plane
        self.origin = shear_wall.origin

        self.plates: List[Quad3D] = []

    def sum_forces(self, combo_name: str = 'Combo 1') -> Tuple[float, float, float, float]:

        # Initialize plate forces to zero
        P, M, V = 0, 0, 0

        # Step through each plate in the coupling beam
        for plate in self.plates:

            xi, yi, zi = _global2local(plate.i_node.X, plate.i_node.Y, plate.i_node.Z, self.origin, self.plane)
            xj, yj, zj = _global2local(plate.j_node.X, plate.j_node.Y, plate.j_node.Z, self.origin, self.plane)
            xn, yn, zn = _global2local(plate.n_node.X, plate.n_node.Y, plate.n_node.Z, self.origin, self.plane)

            # Determine if this plate is at the left edge of the coupling beam
            if isclose(xi, self.x):

                # Check if this is a wall flange plate or a wall web plate
                if isclose(xi, xj):

                    # Plates that form wall flanges should not affect coupling beams, so forces will not be summed
                    pass

                else:

                    # Find and sum the axial forces in this plate
                    Pi = plate.f(combo_name)[0][0]
                    Pn = plate.f(combo_name)[18][0]
                    P += -Pi - Pn

                    # Find and sum the moments about the coupling beam's center in this plate
                    xi_cb = yi - (self.y + self.height/2)
                    xn_cb = yn - (self.y + self.height/2)
                    Mi = plate.f(combo_name)[0][0]*xi_cb
                    Mn = plate.f(combo_name)[18][0]*xn_cb
                    M += -Mi - Mn

                    # Find and sum the shear forces in this plate
                    Vi = -plate.f(combo_name)[1][0]
                    Vn = -plate.f(combo_name)[19][0]

                    V += -Vi - Vn

        # Calculate the shear span ratio
        M_VH = M/(V*self.height)

        # Return the summed forces and shear span ratio
        return P, M, V, M_VH
