
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