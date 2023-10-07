from math import isclose
from PyNite.LoadCombo import LoadCombo
from numpy import array, atleast_2d, zeros

def _prepare_model(model):
    """Prepares a model for analysis by ensuring at least one load combination is defined, generating all meshes that have not already been generated, activating all non-linear members, and internally numbering all nodes and elements.

    :param model: The model being prepared for analysis.
    :type model: FEModel3D
    """

    # Ensure there is at least 1 load combination to solve if the user didn't define any
    if model.LoadCombos == {}:
        # Create and add a default load combination to the dictionary of load combinations
        model.LoadCombos['Combo 1'] = LoadCombo('Combo 1', factors={'Case 1':1.0})
    
    # Generate all meshes
    for mesh in model.Meshes.values():
        if mesh.is_generated == False:
            mesh.generate()

    # Activate all springs and members for all load combinations
    for spring in model.Springs.values():
        for combo_name in model.LoadCombos.keys():
            spring.active[combo_name] = True
    
    # Activate all physical members for all load combinations
    for phys_member in model.Members.values():
        for combo_name in model.LoadCombos.keys():
            phys_member.active[combo_name] = True
    
    # Assign an internal ID to all nodes and elements in the model. This number is different from the name used by the user to identify nodes and elements.
    _renumber(model)

def _identify_combos(model, combo_tags):
    """Returns a list of load combinations that are to be run based on tags given by the user.

    :param model: The model being analyzed.
    :type model: FEModel3D
    :param combo_tags: A list of tags used for the load combinations to be evaluated.
    :type combo_tags: list
    :return: A list containing the load combinations to be analyzed.
    :rtype: list
    """
    
    # Identify which load combinations to evaluate
    if combo_tags is None:
        combo_list = model.LoadCombos.values()
    else:
        combo_list = []
        for combo in model.LoadCombos.values():
            if any(tag in combo.combo_tags for tag in combo_tags):
                combo_list.append(combo)
    
    return combo_list

def _check_stability(model, K):
    """
    Identifies nodal instabilities in a model's stiffness matrix.
    """

    # Initialize the `unstable` flag to `False`
    unstable = False

    # Step through each diagonal term in the stiffness matrix
    for i in range(K.shape[0]):
        
        # Determine which node this term belongs to
        node = [node for node in model.Nodes.values() if node.ID == int(i/6)][0]

        # Determine which degree of freedom this term belongs to
        dof = i%6

        # Check to see if this degree of freedom is supported
        if dof == 0:
            supported = node.support_DX
        elif dof == 1:
            supported = node.support_DY
        elif dof == 2:
            supported = node.support_DZ
        elif dof == 3:
            supported = node.support_RX
        elif dof == 4:
            supported = node.support_RY
        elif dof == 5:
            supported = node.support_RZ

        # Check if the degree of freedom on this diagonal is unstable
        if isclose(K[i, i], 0) and not supported:

            # Flag the model as unstable
            unstable = True

            # Identify which direction this instability effects
            if i%6 == 0: direction = 'for translation in the global X direction.'
            if i%6 == 1: direction = 'for translation in the global Y direction.'
            if i%6 == 2: direction = 'for translation in the global Z direction.'
            if i%6 == 3: direction = 'for rotation about the global X axis.'
            if i%6 == 4: direction = 'for rotation about the global Y axis.'
            if i%6 == 5: direction = 'for rotation about the global Z axis.'

            # Print a message to the console
            print('* Nodal instability detected: node ' + node.name + ' is unstable ' + direction)

    if unstable:
        raise Exception('Unstable node(s). See console output for details.')

    return

def _store_displacements(model, D1, D2, D1_indices, D2_indices, combo):
    """Stores calculated displacements from the solver into the model's displacement vector `_D` and into each node object in the model.

    :param model: The finite element model being evaluated.
    :type model: FEModel3D
    :param D1: An array of calculated displacements
    :type D1: array
    :param D2: An array of enforced displacements
    :type D2: array
    :param D1_indices: A list of the degree of freedom indices for each displacement in D1
    :type D1_indices: list
    :param D2_indices: A list of the degree of freedom indices for each displacement in D2
    :type D2_indices: list
    :param combo: The load combination to store the displacements for
    :type combo: LoadCombo
    """
    
    D = zeros((len(model.Nodes)*6, 1))

    # Step through each node in the model
    for node in model.Nodes.values():
        
        if node.ID*6 + 0 in D2_indices:
            # Get the enforced displacement
            D[(node.ID*6 + 0, 0)] = D2[D2_indices.index(node.ID*6 + 0), 0]
        else:
            # Get the calculated displacement
            D[(node.ID*6 + 0, 0)] = D1[D1_indices.index(node.ID*6 + 0), 0]

        if node.ID*6 + 1 in D2_indices:
            # Get the enforced displacement
            D[(node.ID*6 + 1, 0)] = D2[D2_indices.index(node.ID*6 + 1), 0]
        else:
            # Get the calculated displacement
            D[(node.ID*6 + 1, 0)] = D1[D1_indices.index(node.ID*6 + 1), 0]

        if node.ID*6 + 2 in D2_indices:
            # Get the enforced displacement
            D[(node.ID*6 + 2, 0)] = D2[D2_indices.index(node.ID*6 + 2), 0]
        else:
            # Get the calculated displacement
            D[(node.ID*6 + 2, 0)] = D1[D1_indices.index(node.ID*6 + 2), 0]

        if node.ID*6 + 3 in D2_indices:
            # Get the enforced rotation
            D[(node.ID*6 + 3, 0)] = D2[D2_indices.index(node.ID*6 + 3), 0]
        else:
            # Get the calculated rotation
            D[(node.ID*6 + 3, 0)] = D1[D1_indices.index(node.ID*6 + 3), 0]

        if node.ID*6 + 4 in D2_indices:
            # Get the enforced rotation
            D[(node.ID*6 + 4, 0)] = D2[D2_indices.index(node.ID*6 + 4), 0]
        else:
            # Get the calculated rotation
            D[(node.ID*6 + 4, 0)] = D1[D1_indices.index(node.ID*6 + 4), 0]

        if node.ID*6 + 5 in D2_indices:
            # Get the enforced rotation
            D[(node.ID*6 + 5, 0)] = D2[D2_indices.index(node.ID*6 + 5), 0]
        else:
            # Get the calculated rotation
            D[(node.ID*6 + 5, 0)] = D1[D1_indices.index(node.ID*6 + 5), 0]

    # Save the global displacement vector to the model
    model._D[combo.name] = D

    # Store the calculated global nodal displacements into each node object
    for node in model.Nodes.values():

        node.DX[combo.name] = D[node.ID*6 + 0, 0]
        node.DY[combo.name] = D[node.ID*6 + 1, 0]
        node.DZ[combo.name] = D[node.ID*6 + 2, 0]
        node.RX[combo.name] = D[node.ID*6 + 3, 0]
        node.RY[combo.name] = D[node.ID*6 + 4, 0]
        node.RZ[combo.name] = D[node.ID*6 + 5, 0]

def _check_TC_convergence(model, combo_name='Combo 1', log=True):
    
    # Assume the model has converged until we find out otherwise
    convergence = True
    
    # Provide an update to the console if requested by the user
    if log: print('- Checking for tension/compression-only support spring convergence')
    for node in model.Nodes.values():

        # Check convergence of tension/compression-only spring supports and activate/deactivate them as necessary
        if node.spring_DX[1] is not None:
            if ((node.spring_DX[1] == '-' and node.DX[combo_name] > 0)
            or (node.spring_DX[1] == '+' and node.DX[combo_name] < 0)):
                # Check if the spring is switching from active to inactive
                if node.spring_DX[2] == True: convergence = False
                # Make sure the spring is innactive
                node.spring_DX[2] = False
            elif ((node.spring_DX[1] == '-' and node.DX[combo_name] < 0)
            or (node.spring_DX[1] == '+' and node.DX[combo_name] > 0)):
                # Check if the spring is switching from inactive to active
                if node.spring_DX[2] == False: convergence = False
                # Make sure the spring is active
                node.spring_DX[2] = True
        if node.spring_DY[1] is not None:
            if ((node.spring_DY[1] == '-' and node.DY[combo_name] > 0)
            or (node.spring_DY[1] == '+' and node.DY[combo_name] < 0)):
                # Check if the spring is switching from active to inactive
                if node.spring_DY[2] == True: convergence = False
                # Make sure the spring is innactive
                node.spring_DY[2] = False
            elif ((node.spring_DY[1] == '-' and node.DY[combo_name] < 0)
            or (node.spring_DY[1] == '+' and node.DY[combo_name] > 0)):
                # Check if the spring is switching from inactive to active
                if node.spring_DY[2] == False: convergence = False
                # Make sure the spring is active
                node.spring_DY[2] = True
        if node.spring_DZ[1] is not None:
            if ((node.spring_DZ[1] == '-' and node.DZ[combo_name] > 0)
            or (node.spring_DZ[1] == '+' and node.DZ[combo_name] < 0)):
                # Check if the spring is switching from active to inactive
                if node.spring_DZ[2] == True: convergence = False
                # Make sure the spring is innactive
                node.spring_DZ[2] = False
            elif ((node.spring_DZ[1] == '-' and node.DZ[combo_name] < 0)
            or (node.spring_DZ[1] == '+' and node.DZ[combo_name] > 0)):
                # Check if the spring is switching from inactive to active
                if node.spring_DZ[2] == False: convergence = False
                # Make sure the spring is active
                node.spring_DZ[2] = True
        if node.spring_RX[1] is not None:
            if ((node.spring_RX[1] == '-' and node.RX[combo_name] > 0)
            or (node.spring_RX[1] == '+' and node.RX[combo_name] < 0)):
                # Check if the spring is switching from active to inactive
                if node.spring_RX[2] == True: convergence = False
                # Make sure the spring is innactive
                node.spring_RX[2] = False
            elif ((node.spring_RX[1] == '-' and node.RX[combo_name] < 0)
            or (node.spring_RX[1] == '+' and node.RX[combo_name] > 0)):
                # Check if the spring is switching from inactive to active
                if node.spring_RX[2] == False: convergence = False
                # Make sure the spring is active
                node.spring_RX[2] = True
        if node.spring_RY[1] is not None:
            if ((node.spring_RY[1] == '-' and node.RY[combo_name] > 0)
            or (node.spring_RY[1] == '+' and node.RY[combo_name] < 0)):
                # Check if the spring is switching from active to inactive
                if node.spring_RY[2] == True: convergence = False
                # Make sure the spring is innactive
                node.spring_RY[2] = False
            elif ((node.spring_RY[1] == '-' and node.RY[combo_name] < 0)
            or (node.spring_RY[1] == '+' and node.RY[combo_name] > 0)):
                # Check if the spring is switching from inactive to active
                if node.spring_RY[2] == False: convergence = False
                # Make sure the spring is active
                node.spring_RY[2] = True
        if node.spring_RZ[1] is not None:
            if ((node.spring_RZ[1] == '-' and node.RZ[combo_name] > 0)
            or (node.spring_RZ[1] == '+' and node.RZ[combo_name] < 0)):
                # Check if the spring is switching from active to inactive
                if node.spring_RZ[2] == True: convergence = False
                # Make sure the spring is innactive
                node.spring_RZ[2] = False
            elif ((node.spring_RZ[1] == '-' and node.RZ[combo_name] < 0)
            or (node.spring_RZ[1] == '+' and node.RZ[combo_name] > 0)):
                # Check if the spring is switching from inactive to active
                if node.spring_RZ[2] == False: convergence = False
                # Make sure the spring is active
                node.spring_RZ[2] = True
    
    # TODO: Adjust the code below to allow elements to reactivate on subsequent iterations if deformations at element nodes indicate the member goes back into an active state. This will lead to a less conservative and more realistic analysis. Nodal springs (above) already do this.

    # Check tension/compression-only springs
    if log: print('- Checking for tension/compression-only spring convergence')
    for spring in model.Springs.values():

        if spring.active[combo_name] == True:

            # Check if tension-only conditions exist
            if spring.tension_only == True and spring.axial(combo_name) > 0:
                spring.active[combo_name] = False
                convergence = False
            
            # Check if compression-only conditions exist
            elif spring.comp_only == True and spring.axial(combo_name) < 0:
                spring.active[combo_name] = False
                convergence = False

    # Check tension/compression only members
    if log: print('- Checking for tension/compression-only member convergence')
    for phys_member in model.Members.values():

        # Only run the tension/compression only check if the member is still active
        if phys_member.active[combo_name] == True:

            # Check if tension-only conditions exist
            if phys_member.tension_only == True and phys_member.max_axial(combo_name) > 0:
                phys_member.active[combo_name] = False
                convergence = False

            # Check if compression-only conditions exist
            elif phys_member.comp_only == True and phys_member.min_axial(combo_name) < 0:
                phys_member.active[combo_name] = False
                convergence = False

    # Return whether the TC analysis has converged
    return convergence

def _calc_reactions(model, log=False, combo_tags=None):
    """
    Calculates reactions internally once the model is solved.

    Parameters
    ----------
    log : bool, optional
        Prints updates to the console if set to True. Default is False.
    """

    # Print a status update to the console
    if log: print('- Calculating reactions')

    # Identify which load combinations to evaluate
    if combo_tags is None:
        combo_list = model.LoadCombos.values()
    else:
        combo_list = []
        for combo in model.LoadCombos.values():
            if any(tag in combo.combo_tags for tag in combo_tags):
                combo_list.append(combo)

    # Calculate the reactions node by node
    for node in model.Nodes.values():
        
        # Step through each load combination
        for combo in combo_list:
            
            # Initialize reactions for this node and load combination
            node.RxnFX[combo.name] = 0.0
            node.RxnFY[combo.name] = 0.0
            node.RxnFZ[combo.name] = 0.0
            node.RxnMX[combo.name] = 0.0
            node.RxnMY[combo.name] = 0.0
            node.RxnMZ[combo.name] = 0.0

            # Determine if the node has any supports
            if (node.support_DX or node.support_DY or node.support_DZ 
            or  node.support_RX or node.support_RY or node.support_RZ):

                # Sum the spring end forces at the node
                for spring in model.Springs.values():

                    if spring.i_node == node and spring.active[combo.name] == True:
                        
                        # Get the spring's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        spring_F = spring.F(combo.name)

                        node.RxnFX[combo.name] += spring_F[0, 0]
                        node.RxnFY[combo.name] += spring_F[1, 0]
                        node.RxnFZ[combo.name] += spring_F[2, 0]
                        node.RxnMX[combo.name] += spring_F[3, 0]
                        node.RxnMY[combo.name] += spring_F[4, 0]
                        node.RxnMZ[combo.name] += spring_F[5, 0]

                    elif spring.j_node == node and spring.active[combo.name] == True:
                    
                        # Get the spring's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        spring_F = spring.F(combo.name)
                    
                        node.RxnFX[combo.name] += spring_F[6, 0]
                        node.RxnFY[combo.name] += spring_F[7, 0]
                        node.RxnFZ[combo.name] += spring_F[8, 0]
                        node.RxnMX[combo.name] += spring_F[9, 0]
                        node.RxnMY[combo.name] += spring_F[10, 0]
                        node.RxnMZ[combo.name] += spring_F[11, 0]

                # Step through each physical member in the model
                for phys_member in model.Members.values():

                    # Sum the sub-member end forces at the node
                    for member in phys_member.sub_members.values():
                        
                        if member.i_node == node and phys_member.active[combo.name] == True:
                        
                            # Get the member's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            member_F = member.F(combo.name)

                            node.RxnFX[combo.name] += member_F[0, 0]
                            node.RxnFY[combo.name] += member_F[1, 0]
                            node.RxnFZ[combo.name] += member_F[2, 0]
                            node.RxnMX[combo.name] += member_F[3, 0]
                            node.RxnMY[combo.name] += member_F[4, 0]
                            node.RxnMZ[combo.name] += member_F[5, 0]

                        elif member.j_node == node and phys_member.active[combo.name] == True:
                        
                            # Get the member's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            member_F = member.F(combo.name)
                        
                            node.RxnFX[combo.name] += member_F[6, 0]
                            node.RxnFY[combo.name] += member_F[7, 0]
                            node.RxnFZ[combo.name] += member_F[8, 0]
                            node.RxnMX[combo.name] += member_F[9, 0]
                            node.RxnMY[combo.name] += member_F[10, 0]
                            node.RxnMZ[combo.name] += member_F[11, 0]

                # Sum the plate forces at the node
                for plate in model.Plates.values():

                    if plate.i_node == node:

                        # Get the plate's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        plate_F = plate.F(combo.name)
                
                        node.RxnFX[combo.name] += plate_F[0, 0]
                        node.RxnFY[combo.name] += plate_F[1, 0]
                        node.RxnFZ[combo.name] += plate_F[2, 0]
                        node.RxnMX[combo.name] += plate_F[3, 0]
                        node.RxnMY[combo.name] += plate_F[4, 0]
                        node.RxnMZ[combo.name] += plate_F[5, 0]

                    elif plate.j_node == node:

                        # Get the plate's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        plate_F = plate.F(combo.name)
                
                        node.RxnFX[combo.name] += plate_F[6, 0]
                        node.RxnFY[combo.name] += plate_F[7, 0]
                        node.RxnFZ[combo.name] += plate_F[8, 0]
                        node.RxnMX[combo.name] += plate_F[9, 0]
                        node.RxnMY[combo.name] += plate_F[10, 0]
                        node.RxnMZ[combo.name] += plate_F[11, 0]

                    elif plate.m_node == node:

                        # Get the plate's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        plate_F = plate.F(combo.name)
                
                        node.RxnFX[combo.name] += plate_F[12, 0]
                        node.RxnFY[combo.name] += plate_F[13, 0]
                        node.RxnFZ[combo.name] += plate_F[14, 0]
                        node.RxnMX[combo.name] += plate_F[15, 0]
                        node.RxnMY[combo.name] += plate_F[16, 0]
                        node.RxnMZ[combo.name] += plate_F[17, 0]

                    elif plate.n_node == node:

                        # Get the plate's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        plate_F = plate.F(combo.name)
                
                        node.RxnFX[combo.name] += plate_F[18, 0]
                        node.RxnFY[combo.name] += plate_F[19, 0]
                        node.RxnFZ[combo.name] += plate_F[20, 0]
                        node.RxnMX[combo.name] += plate_F[21, 0]
                        node.RxnMY[combo.name] += plate_F[22, 0]
                        node.RxnMZ[combo.name] += plate_F[23, 0]

                # Sum the quad forces at the node
                for quad in model.Quads.values():

                    if quad.m_node == node:

                        # Get the quad's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        quad_F = quad.F(combo.name)

                        node.RxnFX[combo.name] += quad_F[0, 0]
                        node.RxnFY[combo.name] += quad_F[1, 0]
                        node.RxnFZ[combo.name] += quad_F[2, 0]
                        node.RxnMX[combo.name] += quad_F[3, 0]
                        node.RxnMY[combo.name] += quad_F[4, 0]
                        node.RxnMZ[combo.name] += quad_F[5, 0]

                    elif quad.n_node == node:

                        # Get the quad's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        quad_F = quad.F(combo.name)
                
                        node.RxnFX[combo.name] += quad_F[6, 0]
                        node.RxnFY[combo.name] += quad_F[7, 0]
                        node.RxnFZ[combo.name] += quad_F[8, 0]
                        node.RxnMX[combo.name] += quad_F[9, 0]
                        node.RxnMY[combo.name] += quad_F[10, 0]
                        node.RxnMZ[combo.name] += quad_F[11, 0]

                    elif quad.i_node == node:

                        # Get the quad's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        quad_F = quad.F(combo.name)
                
                        node.RxnFX[combo.name] += quad_F[12, 0]
                        node.RxnFY[combo.name] += quad_F[13, 0]
                        node.RxnFZ[combo.name] += quad_F[14, 0]
                        node.RxnMX[combo.name] += quad_F[15, 0]
                        node.RxnMY[combo.name] += quad_F[16, 0]
                        node.RxnMZ[combo.name] += quad_F[17, 0]

                    elif quad.j_node == node:

                        # Get the quad's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        quad_F = quad.F(combo.name)
                
                        node.RxnFX[combo.name] += quad_F[18, 0]
                        node.RxnFY[combo.name] += quad_F[19, 0]
                        node.RxnFZ[combo.name] += quad_F[20, 0]
                        node.RxnMX[combo.name] += quad_F[21, 0]
                        node.RxnMY[combo.name] += quad_F[22, 0]
                        node.RxnMZ[combo.name] += quad_F[23, 0]
                
                # Sum the joint loads applied to the node
                for load in node.NodeLoads:

                    for case, factor in combo.factors.items():
                        
                        if load[2] == case:

                            if load[0] == 'FX':
                                node.RxnFX[combo.name] -= load[1]*factor
                            elif load[0] == 'FY':
                                node.RxnFY[combo.name] -= load[1]*factor
                            elif load[0] == 'FZ':
                                node.RxnFZ[combo.name] -= load[1]*factor
                            elif load[0] == 'MX':
                                node.RxnMX[combo.name] -= load[1]*factor
                            elif load[0] == 'MY':
                                node.RxnMY[combo.name] -= load[1]*factor
                            elif load[0] == 'MZ':
                                node.RxnMZ[combo.name] -= load[1]*factor
            
            # Calculate reactions due to active spring supports at the node
            else:
                
                if node.spring_DX[0] != None and node.spring_DX[2] == True:
                    sign = node.spring_DX[1]
                    k = node.spring_DX[0]
                    if sign != None: k = float(sign + str(k))
                    DX = node.DX[combo.name]
                    node.RxnFX[combo.name] += k*DX
                if node.spring_DY[0] != None and node.spring_DY[2] == True:
                    sign = node.spring_DY[1]
                    k = node.spring_DY[0]
                    if sign != None: k = float(sign + str(k))
                    DY = node.DY[combo.name]
                    node.RxnFY[combo.name] += k*DY
                if node.spring_DZ[0] != None and node.spring_DZ[2] == True:
                    sign = node.spring_DZ[1]
                    k = node.spring_DZ[0]
                    if sign != None: k = float(sign + str(k))
                    DZ = node.DZ[combo.name]
                    node.RxnFZ[combo.name] += k*DZ
                if node.spring_RX[0] != None and node.spring_RX[2] == True:
                    sign = node.spring_RX[1]
                    k = node.spring_RX[0]
                    if sign != None: k = float(sign + str(k))
                    RX = node.RX[combo.name]
                    node.RxnMX[combo.name] += k*RX
                if node.spring_RY[0] != None and node.spring_RY[2] == True:
                    sign = node.spring_RY[1]
                    k = node.spring_RY[0]
                    if sign != None: k = float(sign + str(k))
                    RY = node.RY[combo.name]
                    node.RxnMY[combo.name] += k*RY
                if node.spring_RZ[0] != None and node.spring_RZ[2] == True:
                    sign = node.spring_RZ[1]
                    k = node.spring_RZ[0]
                    if sign != None: k = float(sign + str(k))
                    RZ = node.RZ[combo.name]
                    node.RxnMZ[combo.name] += k*RZ

def _check_statics(model, combo_tags=None):
    '''
    Checks static equilibrium and prints results to the console.

    Parameters
    ----------
    precision : number
        The number of decimal places to carry the results to.
    '''

    print('+----------------+')
    print('| Statics Check: |')
    print('+----------------+')
    print('')

    from prettytable import PrettyTable

    # Start a blank table and create a header row
    statics_table = PrettyTable()
    statics_table.field_names = ['Load Combination', 'Sum FX', 'Sum RX', 'Sum FY', 'Sum RY', 'Sum FZ', 'Sum RZ', 'Sum MX', 'Sum RMX', 'Sum MY', 'Sum RMY', 'Sum MZ', 'Sum RMZ']

    # Identify which load combinations to evaluate
    if combo_tags is None:
        combo_list = model.LoadCombos.values()
    else:
        combo_list = []
        for combo in model.LoadCombos.values():
            if any(tag in combo.combo_tags for tag in combo_tags):
                combo_list.append(combo)

    # Step through each load combination
    for combo in combo_list:

        # Initialize force and moment summations to zero
        SumFX, SumFY, SumFZ = 0.0, 0.0, 0.0
        SumMX, SumMY, SumMZ = 0.0, 0.0, 0.0
        SumRFX, SumRFY, SumRFZ = 0.0, 0.0, 0.0
        SumRMX, SumRMY, SumRMZ = 0.0, 0.0, 0.0

        # Get the global force vector and the global fixed end reaction vector
        P = model.P(combo.name)
        FER = model.FER(combo.name)

        # Step through each node and sum its forces
        for node in model.Nodes.values():

            # Get the node's coordinates
            X = node.X
            Y = node.Y
            Z = node.Z

            # Get the nodal forces
            FX = P[node.ID*6+0][0] - FER[node.ID*6+0][0]
            FY = P[node.ID*6+1][0] - FER[node.ID*6+1][0]
            FZ = P[node.ID*6+2][0] - FER[node.ID*6+2][0]
            MX = P[node.ID*6+3][0] - FER[node.ID*6+3][0]
            MY = P[node.ID*6+4][0] - FER[node.ID*6+4][0]
            MZ = P[node.ID*6+5][0] - FER[node.ID*6+5][0]

            # Get the nodal reactions
            RFX = node.RxnFX[combo.name]
            RFY = node.RxnFY[combo.name]
            RFZ = node.RxnFZ[combo.name]
            RMX = node.RxnMX[combo.name]
            RMY = node.RxnMY[combo.name]
            RMZ = node.RxnMZ[combo.name]

            # Sum the global forces
            SumFX += FX
            SumFY += FY
            SumFZ += FZ
            SumMX += MX - FY*Z + FZ*Y
            SumMY += MY + FX*Z - FZ*X
            SumMZ += MZ - FX*Y + FY*X

            # Sum the global reactions
            SumRFX += RFX
            SumRFY += RFY
            SumRFZ += RFZ
            SumRMX += RMX - RFY*Z + RFZ*Y
            SumRMY += RMY + RFX*Z - RFZ*X
            SumRMZ += RMZ - RFX*Y + RFY*X   

        # Add the results to the table
        statics_table.add_row([combo.name, '{:.3g}'.format(SumFX), '{:.3g}'.format(SumRFX),
                                            '{:.3g}'.format(SumFY), '{:.3g}'.format(SumRFY),
                                            '{:.3g}'.format(SumFZ), '{:.3g}'.format(SumRFZ),
                                            '{:.3g}'.format(SumMX), '{:.3g}'.format(SumRMX),
                                            '{:.3g}'.format(SumMY), '{:.3g}'.format(SumRMY),
                                            '{:.3g}'.format(SumMZ), '{:.3g}'.format(SumRMZ)])

    # Print the static check table
    print(statics_table)
    print('')
    
def _partition_D(model):
    """Builds a list with known nodal displacements and with the positions in global stiffness
        matrix of known and unknown nodal displacements

    :return: A list of the global matrix indices for the unknown nodal displacements (D1_indices). A
                list of the global matrix indices for the known nodal displacements (D2_indices). A list
                of the known nodal displacements (D2).
    :rtype: list, list, list
    """

    D1_indices = [] # A list of the indices for the unknown nodal displacements
    D2_indices = [] # A list of the indices for the known nodal displacements
    D2 = []         # A list of the values of the known nodal displacements (D != None)

    # Create the auxiliary table
    for node in model.Nodes.values():
        
        # Unknown displacement DX
        if node.support_DX==False and node.EnforcedDX == None:
            D1_indices.append(node.ID*6 + 0)
        # Known displacement DX
        elif node.EnforcedDX != None:
            D2_indices.append(node.ID*6 + 0)
            D2.append(node.EnforcedDX)
        # Support at DX
        else:
            D2_indices.append(node.ID*6 + 0)
            D2.append(0.0)

        # Unknown displacement DY
        if node.support_DY == False and node.EnforcedDY == None:
            D1_indices.append(node.ID*6 + 1)
        # Known displacement DY
        elif node.EnforcedDY != None:
            D2_indices.append(node.ID*6 + 1)
            D2.append(node.EnforcedDY)
        # Support at DY
        else:
            D2_indices.append(node.ID*6 + 1)
            D2.append(0.0)

        # Unknown displacement DZ
        if node.support_DZ == False and node.EnforcedDZ == None:
            D1_indices.append(node.ID*6 + 2)
        # Known displacement DZ
        elif node.EnforcedDZ != None:
            D2_indices.append(node.ID*6 + 2)
            D2.append(node.EnforcedDZ)
        # Support at DZ
        else:
            D2_indices.append(node.ID*6 + 2)
            D2.append(0.0)

        # Unknown displacement RX
        if node.support_RX == False and node.EnforcedRX == None:
            D1_indices.append(node.ID*6 + 3)
        # Known displacement RX
        elif node.EnforcedRX != None:
            D2_indices.append(node.ID*6 + 3)
            D2.append(node.EnforcedRX)
        # Support at RX
        else:
            D2_indices.append(node.ID*6 + 3)
            D2.append(0.0)

        # Unknown displacement RY
        if node.support_RY == False and node.EnforcedRY == None:
            D1_indices.append(node.ID*6 + 4)
        # Known displacement RY
        elif node.EnforcedRY != None:
            D2_indices.append(node.ID*6 + 4)
            D2.append(node.EnforcedRY)
        # Support at RY
        else:
            D2_indices.append(node.ID*6 + 4)
            D2.append(0.0)

        # Unknown displacement RZ
        if node.support_RZ == False and node.EnforcedRZ == None:
            D1_indices.append(node.ID*6 + 5)
        # Known displacement RZ
        elif node.EnforcedRZ != None:
            D2_indices.append(node.ID*6 + 5)
            D2.append(node.EnforcedRZ)
        # Support at RZ
        else:
            D2_indices.append(node.ID*6 + 5)
            D2.append(0.0)
    
    # Legacy code on the next line. I will leave it here until the line that follows has been proven over time.
    # D2 = atleast_2d(D2)
    
    # Convert D2 from a list to a matrix
    D2 = array(D2, ndmin=2).T

    # Return the indices and the known displacements
    return D1_indices, D2_indices, D2

def _renumber(model):
    """
    Assigns node and element ID numbers to be used internally by the program. Numbers are
    assigned according to the order in which they occur in each dictionary.
    """
    
    # Number each node in the model
    for id, node in enumerate(model.Nodes.values()):
        node.ID = id
    
    # Number each spring in the model
    for id, spring in enumerate(model.Springs.values()):
        spring.ID = id

    # Descritize all the physical members and number each member in the model
    id = 0
    for phys_member in model.Members.values():
        phys_member.descritize()
        for member in phys_member.sub_members.values():
            member.ID = id
            id += 1
    
    # Number each plate in the model
    for id, plate in enumerate(model.Plates.values()):
        plate.ID = id
    
    # Number each quadrilateral in the model
    for id, quad in enumerate(model.Quads.values()):
        quad.ID = id