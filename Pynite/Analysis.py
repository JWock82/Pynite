from __future__ import annotations  # Allows more recent type hints features
from typing import TYPE_CHECKING
from math import isclose

from numpy import array, atleast_2d, zeros, subtract, matmul, divide, seterr, nanmax
from numpy.linalg import solve

from Pynite.LoadCombo import LoadCombo

if TYPE_CHECKING:
    from typing import List, Tuple
    from Pynite.FEModel3D import FEModel3D
    from numpy import float64
    from numpy.typing import NDArray
    from scipy.sparse import lil_matrix


def _prepare_model(model: FEModel3D, n_modes: int = 0) -> None:
    """Prepares a model for analysis by ensuring at least one load combination is defined, generating all meshes that have not already been generated, activating all non-linear members, and internally numbering all nodes and elements.

    :param model: The model being prepared for analysis.
    :type model: FEModel3D
    :param n_modes: The number of modes to be used for modal analysis. Default is 0.
    :type n_modes: int
    """

    # Reset any nodal displacements
    model._D = {}
    for node in model.nodes.values():
        node.DX = {}
        node.DY = {}
        node.DZ = {}
        node.RX = {}
        node.RY = {}
        node.RZ = {}

    # Ensure there is at least 1 load combination to solve if the user didn't define any
    if model.load_combos == {}:
        # Create and add a default load combination to the dictionary of load combinations
        model.load_combos['Combo 1'] = LoadCombo('Combo 1', factors={'Case 1':1.0})

    # Remove any modal load combinations from prior modal analyses
    to_remove = [name for name, combo in model.load_combos.items() if combo.combo_tags != None and 'modal' in combo.combo_tags]
    for name in to_remove:
        model.load_combos.pop(name)

    # Set up new load combinations for modal analysis if necessary
    for i in range(n_modes):
        model.add_load_combo(f'Mode {i + 1}', {}, ['modal'])

    # Generate all basic meshes
    for mesh in model.meshes.values():
        if mesh.is_generated == False:
            mesh.generate()

    # Generate all shear wall meshes
    for shear_wall in model.shear_walls.values():
        shear_wall.generate()

    # Generate all mat foundation meshes
    for mat in model.mats.values():
        mat.generate()

    # Activate all springs and members for all load combinations
    for spring in model.springs.values():
        for combo_name in model.load_combos.keys():
            spring.active[combo_name] = True

    # Activate all physical members for all load combinations
    for phys_member in model.members.values():
        for combo_name in model.load_combos.keys():
            phys_member.active[combo_name] = True

    # Assign an internal ID to all nodes and elements in the model. This number is different from the name used by the user to identify nodes and elements.
    _renumber(model)


def _identify_combos(model: FEModel3D, combo_tags: List[str] | None = None) -> List[LoadCombo]:
    """Returns a list of load combinations that are to be run based on tags given by the user.

    :param model: The model being analyzed.
    :type model: FEModel3D
    :param combo_tags: A list of tags used for the load combinations to be evaluated. Defaults to `None` in which case all load combinations will be added to the list of load combinations to be run.
    :type combo_tags: list, optional
    :return: A list containing the load combinations to be analyzed.
    :rtype: list
    """

    # Identify which load combinations to evaluate
    if combo_tags is None:
        # Evaluate all load combinations if not tags have been provided
        combo_list = list(model.load_combos.values())
    else:
        # Initialize the list of load combinations to be evaluated
        combo_list = []
        # Step through each load combination in the model
        for combo in model.load_combos.values():
            # Check if this load combination is tagged with any of the tags we're looking for
            if combo.combo_tags is not None and any(tag in combo.combo_tags for tag in combo_tags):
                # Add the load combination to the list of load combinations to be evaluated
                combo_list.append(combo)

    # Return the list of load combinations to be evaluated
    return combo_list


def _check_stability(model: FEModel3D, K: NDArray[float64]) -> None:
    """
    Identifies nodal instabilities in a model's stiffness matrix.
    """

    # Initialize the `unstable` flag to `False`
    unstable = False

    # Step through each diagonal term in the stiffness matrix
    for i in range(K.shape[0]):

        # Determine which node this term belongs to
        node = [node for node in model.nodes.values() if node.ID == int(i/6)][0]

        # Determine which degree of freedom this term belongs to
        dof = i % 6

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

            # Identify which direction this instability affects
            if i%6 == 0:
                direction = 'for translation in the global X direction.'
            if i%6 == 1:
                direction = 'for translation in the global Y direction.'
            if i%6 == 2:
                direction = 'for translation in the global Z direction.'
            if i%6 == 3:
                direction = 'for rotation about the global X axis.'
            if i%6 == 4:
                direction = 'for rotation about the global Y axis.'
            if i%6 == 5:
                direction = 'for rotation about the global Z axis.'

            # Print a message to the console
            print('* Nodal instability detected: node ' + node.name + ' is unstable ' + direction)

    if unstable:
        raise Exception('Unstable node(s). See console output for details.')

    return


def _PDelta(model: FEModel3D, combo_name: str, P1: NDArray[float64], FER1: NDArray[float64], D1_indices: List[int], D2_indices: List[int], D2: NDArray[float64], log: bool = True, sparse: bool = True, check_stability: bool = False, max_iter: int = 30) -> None:
    """Performs second order (P-Delta) analysis. This type of analysis is appropriate for most models using beams, columns and braces. Second order analysis is usually required by material-specific codes. Models with slender members and/or members with combined bending and axial loads will generally have more significant P-Delta effects. P-Delta effects in plates/quads are not considered by Pynite at this time.

    :param model: The finite element model to be solved.
    :type: FEModel3D
    :param combo_name: The name of the load combination to evaluate.
    :type combo_name: str
    :param P1: An array of the loads to apply.
    :type P1: numpy array
    :param FER1: An array of the fixed end reactions.
    :type FER1: numpy array
    :param log: Prints updates to the console if set to True. Default == False.
    :type log: bool, optional
    :param check_stability: When set to True, checks the stiffness matrix for any unstable degrees of freedom and reports them back to the console. This does add to the solution time. Defaults to True.
    :type check_stability: bool, optional
    :param max_iter: The maximum number of iterations permitted. If this value is exceeded the program will report divergence. Defaults to 30.
    :type max_iter: int, optional
    :param sparse: Indicates whether the sparse matrix solver should be used. A matrix can be considered sparse or dense depening on how many zero terms there are. Structural stiffness matrices often contain many zero terms. The sparse solver can offer faster solutions for such matrices. Using the sparse solver on dense matrices may lead to slower solution times. Be sure ``scipy`` is installed to use the sparse solver. Default is True.
    :type sparse: bool, optional
    :param check_stability: Indicates whether nodal stability should be checked. This slows down the analysis considerably, but can be useful for small models or for debugging. Default is `False`.
    :type check_stability: bool, optional
    :raises ValueError: Occurs when there is a singularity in the stiffness matrix, which indicates an unstable structure.
    :raises Exception: Occurs when a model fails to converge.
    """

    # Import `scipy` features if the sparse solver is being used
    if sparse == True:
        from scipy.sparse.linalg import spsolve

    convergence_TC = False  # Tracks tension/compression-only convergence
    divergence_TC = False   # Tracks tension/compression-only divergence
    iter_count_TC = 1

    # Iterate until either T/C convergence or divergence occurs. Perform at least 2 iterations for the P-Delta analysis.
    while convergence_TC == False and divergence_TC == False:

        # Inform the user which iteration we're on
        if log:
            print('- Beginning tension/compression-only iteration #' + str(iter_count_TC))

        # There will be 2 solution steps in the P-Delta analysis:
        # Step 1 - Analyze based on initial stiffness
        # Step 2 - Calculate geometric stiffness and add its effects
        for solution_step in [1, 2]:

            # Determine if the user has selected a sparse solution
            if sparse == True:

                if solution_step == 1:

                    # Calculate the partitioned initial stiffness matrices. These matrices must be recalculated on each T/C iteration due to tension/compression-only members deactivating or reactivating.
                    if log:
                        print('- Calculating initial stiffness matrix')
                    K11, K12, K21, K22 = _partition(model, model.K(combo_name, log, check_stability, sparse).tolil(), D1_indices, D2_indices)

                    # The initial stiffness matrices are currently `lil` format which is great for memory, but slow for mathematical operations. They will be converted to `csr` format.
                    K11 = K11.tocsr()
                    K12 = K12.tocsr()
                    K21 = K21.tocsr()
                    K22 = K22.tocsr()

                # Check if we are ready to calculate the geometric stiffness
                if solution_step == 2:

                    # After the first iteration, the geometric stiffness matrix will be added to the linear elastic stiffness matrix.
                    if log: print('- Calculating geometric stiffness matrix')
                    Kg11, Kg12, Kg21, Kg22 = _partition(model, model.Kg(combo_name, log, sparse, False), D1_indices, D2_indices)

                    # The Kg stiffness matrices are currently `lil` format which is great for memory, but slow for mathematical operations. They will be converted to `csr` format. Note that the `+` operator performs matrix addition on `csr` matrices.
                    if log: print('- Summing initial & geometric stiffness matrices')
                    K11 = K11 + Kg11.tocsr()
                    K12 = K12 + Kg12.tocsr()
                    K21 = K21 + Kg21.tocsr()
                    K22 = K22 + Kg22.tocsr()

            # Determine if the user has selected a dense solution
            else:

                if solution_step == 1:

                    # Calculate the partitioned initial stiffness matrices. These matrices must be recalculated on each T/C iteration due to tension/compression-only members deactivating or reactivating.
                    K11, K12, K21, K22 = _partition(model, model.K(combo_name, log, check_stability, sparse), D1_indices, D2_indices)

                # Check if we are ready to calculate the geometric stiffness
                if solution_step == 2:

                    # After the first iteration, the geometric stiffness matrix will be added to the linear elastic stiffness matrix.
                    Kg11, Kg12, Kg21, Kg22 = _partition(model.Kg(combo_name, log, sparse, False), D1_indices, D2_indices)

                    K11 = K11 + Kg11
                    K12 = K12 + Kg12
                    K21 = K21 + Kg21
                    K22 = K22 + Kg22

            # Calculate the global displacement vector
            if log:
                print('- Calculating the global displacement vector')
            if K11.shape == (0, 0):
                # All displacements are known, so D1 is an empty vector
                D1 = []
            else:
                try:
                    # Calculate the displacements, `D1`
                    if sparse == True:
                        # The partitioned stiffness matrix is already in `csr` format. The `@` operator performs matrix multiplication on sparse matrices.
                        D1 = spsolve(K11.tocsr(), subtract(subtract(P1, FER1), K12.tocsr() @ D2))
                        D1 = D1.reshape(len(D1), 1)
                    else:
                        # The partitioned stiffness matrix is in `csr` format. It will be converted to a 2D dense array for mathematical operations.
                        D1 = solve(K11, subtract(subtract(P1, FER1), matmul(K12, D2)))

                except:
                    # Return out of the method if 'K' is singular and provide an error message
                    raise ValueError('The stiffness matrix is singular, which indicates that the structure is unstable.')

            # Store the calculated displacements
            _store_displacements(model, D1, D2, D1_indices, D2_indices, model.load_combos[combo_name])

        # Check whether the tension/compression-only analysis has converged and deactivate any members that are showing forces they can't hold
        convergence_TC = _check_TC_convergence(model, combo_name, log)

        # Report on convergence of tension/compression only analysis
        if convergence_TC == False:

            if log:
                print('- Tension/compression-only analysis did not converge on this iteration')
                print('- Tension/compression-only members will be deactivated or reactivated as necessary')
                print('- P-Delta analysis will be rerun')

            # Increment the tension/compression-only iteration count
            iter_count_TC += 1

        else:
            if log:
                print('- Tension/compression-only analysis converged after ' + str(iter_count_TC) + ' iteration(s)')

        # Check for divergence in the tension/compression-only analysis
        if iter_count_TC > max_iter:
            divergence_TC = True
            raise Exception('- Model diverged during tension/compression-only analysis')

    # Flag the model as solved
    model.solution = 'P-Delta'


def _pushover_step(model: FEModel3D, combo_name: str, push_combo: str, step_num: int, P1: NDArray[float64], FER1: NDArray[float64], D1_indices: List[int], D2_indices: List[int], D2: NDArray[float64], log: bool = True, sparse: bool = True, check_stability: bool = False) -> None:

    # Run at least one iteration
    run_step = True

    # Run/rerun the load step until convergence occurs
    while run_step == True:

        # Calculate the partitioned global stiffness matrices
        # Sparse solver
        if sparse == True:

            from scipy.sparse.linalg import spsolve

            # Calculate the initial stiffness matrix
            if log: print('- Calculating elastic stiffness matrix [Ke]')
            K11, K12, K21, K22 = _partition(model, model.K(combo_name, log, check_stability, sparse).tolil(), D1_indices, D2_indices)

            # Calculate the geometric stiffness matrix
            # The `combo_name` variable in the code below is not the name of the pushover load combination. Rather it is the name of the primary combination that the pushover load will be added to. Axial loads used to develop Kg are calculated from the displacements stored in `combo_name`.
            if log: print('- Calculating geometric stiffness matrix [Kg]')
            Kg11, Kg12, Kg21, Kg22 = _partition(model, model.Kg(combo_name, log, sparse, False).tolil(), D1_indices, D2_indices)

            # Calculate the stiffness reduction matrix
            if log: print('- Calculating plastic reduction matrix [Km]')
            Km11, Km12, Km21, Km22 = _partition(model, model.Km(combo_name, push_combo, step_num, log, sparse).tolil(), D1_indices, D2_indices)

            # The stiffness matrices are currently `lil` format which is great for
            # memory, but slow for mathematical operations. They will be converted to
            # `csr` format. The `+` operator performs matrix addition on `csr`
            # matrices.
            K11 = K11.tocsr() + Kg11.tocsr() + Km11.tocsr()
            K12 = K12.tocsr() + Kg12.tocsr() + Km12.tocsr()
            K21 = K21.tocsr() + Kg21.tocsr() + Km21.tocsr()
            K22 = K22.tocsr() + Kg22.tocsr() + Km22.tocsr()

        # Dense solver
        else:

            # Initial stiffness matrix
            if log: print('- Calculating elastic stiffness matrix [Ke]')
            K11, K12, K21, K22 = _partition(model, model.K(combo_name, log, check_stability, sparse), D1_indices, D2_indices)

            # Geometric stiffness matrix
            # The `combo_name` variable in the code below is not the name of the pushover load combination. Rather it is the name of the primary combination that the pushover load will be added to. Axial loads used to develop Kg are calculated from the displacements stored in `combo_name`.
            if log: print('Calculating geometric stiffness matrix [Kg]')
            Kg11, Kg12, Kg21, Kg22 = _partition(model.Kg(combo_name, log, sparse, False), D1_indices, D2_indices)

            # Calculate the stiffness reduction matrix
            if log: print('Calculating plastic reduction matrix [Km]')
            Km11, Km12, Km21, Km22 = _partition(model, model.Km(combo_name, push_combo, step_num, log, sparse), D1_indices, D2_indices)

            K11 = K11 + Kg11 + Km11
            K12 = K12 + Kg12 + Km12
            K21 = K21 + Kg21 + Km21
            K22 = K22 + Kg22 + Km22

        # Calculate the changes to the global displacement vector
        if log: print('- Calculating changes to the global displacement vector')
        if K11.shape == (0, 0):
            # All displacements are known, so D1 is an empty vector
            Delta_D1 = []
        else:
            try:
                # Calculate the change in the displacements Delta_D1
                if sparse == True:
                    # The partitioned stiffness matrix is already in `csr` format. The `@`
                    # operator performs matrix multiplication on sparse matrices.
                    Delta_D1 = spsolve(K11.tocsr(), subtract(subtract(P1, FER1), K12.tocsr() @ D2))
                    Delta_D1 = Delta_D1.reshape(len(Delta_D1), 1)
                else:
                    # The partitioned stiffness matrix is in `csr` format. It will be
                    # converted to a 2D dense array for mathematical operations.
                    Delta_D1 = solve(K11, subtract(subtract(P1, FER1), matmul(K12, D2)))

            except:
                # Return out of the method if 'K' is singular and provide an error message
                raise ValueError('The structure is unstable. Unable to proceed any further with analysis.')

        # Unpartition the displacement results from the analysis step
        Delta_D = _unpartition_disp(model, Delta_D1, D2, D1_indices, D2_indices)

        # Assume no need to rerun this load step due to plastic load reversal until we prove it otherwise
        run_step = False

        # Step through each member in the model and check for plastic load reversal
        for phys_member in model.members.values():

            for sub_member in phys_member.sub_members.values():

                # print(f'Member {sub_member.name} lambda = {sub_member.lamb(Delta_D, combo_name, push_combo, step_num)}')

                # Check for plastic load reversal at the i-node in this load step
                if sub_member.lamb(Delta_D, combo_name, push_combo, step_num)[0, 0] < 0:

                    # Flag the load step for reanalysis
                    sub_member.i_reversal = True
                    run_step = True
                    print(f'-Plastic load reversal encountered at member {sub_member.name} i-node')

                else:

                    sub_member.i_reversal = False

                # Check for plastic load reversal at the j-node in this load step
                if sub_member.lamb(Delta_D, combo_name, push_combo, step_num)[1, 0] < 0:

                    # Flag the load step for reanalysis
                    sub_member.j_reversal = True
                    run_step = True
                    print(f'-Plastic load reversal encountered at member {sub_member.name} j-node')

                else:

                    sub_member.j_reversal = False

        # Undo the last loadstep if plastic load reversal was discovered. We'll rerun it with the corresponding gradients set to zero vectors.
        if run_step == True:
            _sum_displacements(model, -Delta_D1, D2, D1_indices, D2_indices, model.load_combos[combo_name])
            if log: print('- Restarting load step due to plastic load reversal')

    # Sum the calculated displacements
    _sum_displacements(model, Delta_D1, D2, D1_indices, D2_indices, model.load_combos[combo_name])


def _unpartition_disp(model: FEModel3D, D1: NDArray[float64], D2: NDArray[float64], D1_indices: List[int], D2_indices: List[int]) -> NDArray[float64]:
    """Unpartitions displacements from the solver and returns them as a global displacement vector

    :param model: The finite element model being evaluated
    :type model: FEModel3D
    :param D1: An array of calculated displacements
    :type D1: array
    :param D2: An array of enforced displacements
    :type D2: array
    :param D1_indices: A list of the degree of freedom indices for each displacement in D1
    :type D1_indices: list
    :param D2_indices: A list of the degree of freedom indices for each displacement in D2
    :type D2_indices: list
    :return: Global displacement matrix
    :rtype: array
    """

    D = zeros((len(model.nodes)*6, 1))

    # Step through each node in the model
    for node in model.nodes.values():

        # Step through each degree of freedom at the node
        for i in range(6):

            # Check if the dof is in the list of enforced displacements
            if node.ID*6 + i in D2_indices:
                # Get the enforced displacement
                D[(node.ID*6 + i, 0)] = D2[D2_indices.index(node.ID*6 + i), 0]
            else:
                # Get the calculated displacement
                D[(node.ID*6 + i, 0)] = D1[D1_indices.index(node.ID*6 + i), 0]

    # Return the displacement vector
    return D


def _store_displacements(model: FEModel3D, D1: NDArray[float64], D2: NDArray[float64], D1_indices: List[int], D2_indices: List[int], combo: LoadCombo) -> None:
    """Stores calculated displacements from the solver into the model's displacement vector `_D` and into each node object in the model

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
    :param combo: The load combination to store the displacements for.
    :type combo: LoadCombo
    """

    # The raw results from the solver are partitioned. Unpartition them.
    D = _unpartition_disp(model, D1, D2, D1_indices, D2_indices)

    if combo != None:

        # Store the displacements in the model's global displacement vector
        model._D[combo.name] = D

        # Store the calculated global nodal displacements into each node object
        for node in model.nodes.values():

            node.DX[combo.name] = D[node.ID*6 + 0, 0]
            node.DY[combo.name] = D[node.ID*6 + 1, 0]
            node.DZ[combo.name] = D[node.ID*6 + 2, 0]
            node.RX[combo.name] = D[node.ID*6 + 3, 0]
            node.RY[combo.name] = D[node.ID*6 + 4, 0]
            node.RZ[combo.name] = D[node.ID*6 + 5, 0]


def _sum_displacements(model: FEModel3D, Delta_D1: NDArray[float64], Delta_D2: NDArray[float64], D1_indices: List[int], D2_indices: List[int], combo: LoadCombo) -> None:
    """Sums calculated displacements for a load step from the solver into the model's displacement vector `_D` and into each node object in the model.

    :param model: The finite element model being evaluated.
    :type model: FEModel3D
    :param Delta_D1: An array of calculated displacements for a load step
    :type Delta_D1: array
    :param Delta_D2: An array of enforced displacements for a load step
    :type Delta_D2: array
    :param D1_indices: A list of the degree of freedom indices for each displacement in D1
    :type D1_indices: list
    :param D2_indices: A list of the degree of freedom indices for each displacement in D2
    :type D2_indices: list
    :param combo: The load combination to store the displacements for
    :type combo: LoadCombo
    """

    # The raw results from the solver are partitioned. Unpartition them.
    Delta_D = _unpartition_disp(model, Delta_D1, Delta_D2, D1_indices, D2_indices)

    # Sum the load step's global displacement vector with the model's global displacement vector
    model._D[combo.name] += Delta_D

    # Sum the load step's calculated global nodal displacements to each node object's global displacement
    for node in model.nodes.values():

        node.DX[combo.name] += Delta_D[node.ID*6 + 0, 0]
        node.DY[combo.name] += Delta_D[node.ID*6 + 1, 0]
        node.DZ[combo.name] += Delta_D[node.ID*6 + 2, 0]
        node.RX[combo.name] += Delta_D[node.ID*6 + 3, 0]
        node.RY[combo.name] += Delta_D[node.ID*6 + 4, 0]
        node.RZ[combo.name] += Delta_D[node.ID*6 + 5, 0]

def _check_TC_convergence(model: FEModel3D, combo_name: str = "Combo 1", log: bool = True, spring_tolerance: float = 0, member_tolerance: float = 0) -> bool:
    """Checks for convergence in tension-only and compression-only analysis.

    This function evaluates the status of tension-only and compression-only springs and members
    within the finite element model for a given load combination. Its primary purpose is to
    determine if the non-linear analysis has converged by checking if any adjustments to the
    active status of these elements are required.

    The function performs the following checks:
    *   **Nodal Spring Convergence**: It iterates through each nodal spring and assesses whether
        its active state (whether it is resisting force) aligns with the current displacement
        and its tension-only or compression-only designation. If a spring's active
        status needs to change based on the displacement and a specified tolerance, the analysis
        is flagged as not converged, and the spring's `active` status is updated.
    *   **Member Convergence**: It checks each physical member in the model. If a
        tension-only member is active but has a maximum axial force greater than the
        `member_tolerance` (indicating compression), or if a compression-only member is active
        but has a minimum axial force less than `-member_tolerance` (indicating tension),
        that physical member and all its sub-members are deactivated. This action
        flags the analysis as not converged. A future enhancement is noted to allow elements to
        reactivate if deformations indicate they should return to an active state.
    *   **Sub-member Reset**: After checking, the `_solved_combo` flag for all sub-members is
        reset to `None`. This ensures that they will be resegmented and re-evaluated in subsequent
        iterations of the analysis, allowing for further changes as needed for convergence.

    :param model: The finite element model currently being evaluated.
    :type model: FEModel3D
    :param combo_name: The name of the load combination for which the convergence check is
        being performed. Defaults to "Combo 1".
    :type combo_name: str, optional
    :param log: A boolean flag. If `True`, the function prints status updates and information
        about the convergence check to the console. Defaults to `True`.
    :type log: bool, optional
    :param spring_tolerance: A float value representing the displacement tolerance for
        tension/compression-only support springs. Springs will switch active state if their
        displacement is beyond this tolerance in the disallowed direction. Defaults to 0.
    :type spring_tolerance: float, optional
    :param member_tolerance: A float value representing the axial force tolerance for
        tension/compression-only members. Members will be deactivated if their axial force
        is beyond this tolerance in the disallowed direction. Defaults to 0.
    :type member_tolerance: float, optional
    :return: **True** if all tension-only and compression-only elements have converged (no
        changes in active status are needed), **False** otherwise (indicating that further
        iterations of the analysis are required).
    :rtype: bool
    """

    # Assume the model has converged until we find out otherwise
    convergence = True

    # Provide an update to the console if requested by the user
    if log:
        print("- Checking for tension/compression-only support spring convergence")

    # Loop through each node and each directional spring to check and update their active status
    for node in model.nodes.values():

        for direction in ["DX", "DY", "DZ", "RX", "RY", "RZ"]:
            spring = getattr(node, f"spring_{direction}")
            displacement = getattr(node, direction)[combo_name]

            if spring[1] is not None:
                # Determine if the spring should be active based on its designation ('+' or '-') and the displacement
                should_be_active = (
                    spring[1] == "-" and displacement <= -spring_tolerance
                ) or (spring[1] == "+" and displacement >= spring_tolerance)

                # Check if there's a need to switch the active state of the spring
                if spring[2] != should_be_active:
                    spring[2] = should_be_active
                    convergence = False

    # TODO: Adjust the code below to allow elements to reactivate on subsequent iterations if deformations at element nodes indicate the member goes back into an active state. This will lead to a less conservative and more realistic analysis. Nodal springs (above) already do this.

    # Check tension/compression-only springs
    if log:
        print('- Checking for tension/compression-only spring convergence')

    for spring in model.springs.values():

        if spring.active[combo_name] == True:

            # Check if tension-only conditions exist
            if (spring.tension_only == True) and (spring.axial(combo_name) > spring_tolerance):
                if log:
                    print(f'- Deactivating spring {spring.name}')
                spring.active[combo_name] = False
                convergence = False

            # Check if compression-only conditions exist
            elif (spring.comp_only == True) and (spring.axial(combo_name) < -spring_tolerance):
                if log:
                    print(f'- Deactivating spring {spring.name}')
                spring.active[combo_name] = False
                convergence = False

    # Check tension/compression only members
    if log:
        print('- Checking for tension/compression-only member convergence')
    for phys_member in model.members.values():

        # Only run the tension/compression only check if the member is still active
        if phys_member.active[combo_name] == True:

            # Check if a tension-only conditions exist
            if phys_member.tension_only == True and phys_member.max_axial(combo_name) > member_tolerance:

                # Deactivate the physical member
                if log:
                    print(f'- Deactivating member {phys_member.name}')
                phys_member.active[combo_name] = False

                # Deactivate all the sub-members
                for sub_member in phys_member.sub_members.values():
                    sub_member.active[combo_name] = False

                # Flag the analysis as not converged
                convergence = False

            # Check if a compression-only conditions exist
            elif phys_member.comp_only == True and phys_member.min_axial(combo_name) < - member_tolerance:

                # Deactivate the physical member
                if log:
                    print(f'- Deactivating member {phys_member.name}')
                phys_member.active[combo_name] = False

                # Deactivate all the sub-members
                for sub_member in phys_member.sub_members.values():
                    sub_member.active[combo_name] = False

                # Flag the analysis as not converged
                convergence = False

        # Reset the sub-member's flag to unsolved. This will allow it to resolve for the same load combination after subsequent iterations have made further changes.
        for sub_member in phys_member.sub_members.values():
            sub_member._solved_combo = None

    # Return whether the TC analysis has converged
    return convergence


def _calc_reactions(model: FEModel3D, log: bool = False, combo_tags: List[str] | None = None) -> None:
    """
    Calculates reactions internally once the model is solved.

    Parameters
    ----------
    model : FEModel3D
        The finite element model to calculate reactions for.
    log : bool, optional
        Prints updates to the console if set to True. Default == False.
    combo_tags : string, optional
        A list of tags that will be used to identify which load combinations need their reactions calculated. If set to `None` then all load combinations will have their reactions calculated. Default is `None`.
    """

    # Print a status update to the console
    if log:
        print('- Calculating reactions')

    # Identify which load combinations to evaluate
    combo_list = _identify_combos(model, combo_tags)

    # Calculate the reactions node by node
    for node in model.nodes.values():

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
                for spring in model.springs.values():

                    if spring.i_node == node and spring.active[combo.name] == True:

                        # Get the spring's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        spring_F = spring.F(combo.name)

                        if node.support_DX: node.RxnFX[combo.name] += spring_F[0, 0]
                        if node.support_DY: node.RxnFY[combo.name] += spring_F[1, 0]
                        if node.support_DZ: node.RxnFZ[combo.name] += spring_F[2, 0]
                        if node.support_RX: node.RxnMX[combo.name] += spring_F[3, 0]
                        if node.support_RY: node.RxnMY[combo.name] += spring_F[4, 0]
                        if node.support_RZ: node.RxnMZ[combo.name] += spring_F[5, 0]

                    elif spring.j_node == node and spring.active[combo.name] == True:

                        # Get the spring's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        spring_F = spring.F(combo.name)

                        if node.support_DX: node.RxnFX[combo.name] += spring_F[6, 0]
                        if node.support_DY: node.RxnFY[combo.name] += spring_F[7, 0]
                        if node.support_DZ: node.RxnFZ[combo.name] += spring_F[8, 0]
                        if node.support_RX: node.RxnMX[combo.name] += spring_F[9, 0]
                        if node.support_RY: node.RxnMY[combo.name] += spring_F[10, 0]
                        if node.support_RZ: node.RxnMZ[combo.name] += spring_F[11, 0]

                # Step through each physical member in the model
                for phys_member in model.members.values():

                    # Sum the sub-member end forces at the node
                    for member in phys_member.sub_members.values():

                        if member.i_node == node and phys_member.active[combo.name] == True:

                            # Get the member's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            member_F = member.F(combo.name)

                            if node.support_DX: node.RxnFX[combo.name] += member_F[0, 0]
                            if node.support_DY: node.RxnFY[combo.name] += member_F[1, 0]
                            if node.support_DZ: node.RxnFZ[combo.name] += member_F[2, 0]
                            if node.support_RX: node.RxnMX[combo.name] += member_F[3, 0]
                            if node.support_RY: node.RxnMY[combo.name] += member_F[4, 0]
                            if node.support_RZ: node.RxnMZ[combo.name] += member_F[5, 0]

                        elif member.j_node == node and phys_member.active[combo.name] == True:

                            # Get the member's global force matrix
                            # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                            member_F = member.F(combo.name)

                            if node.support_DX: node.RxnFX[combo.name] += member_F[6, 0]
                            if node.support_DY: node.RxnFY[combo.name] += member_F[7, 0]
                            if node.support_DZ: node.RxnFZ[combo.name] += member_F[8, 0]
                            if node.support_RX: node.RxnMX[combo.name] += member_F[9, 0]
                            if node.support_RY: node.RxnMY[combo.name] += member_F[10, 0]
                            if node.support_RZ: node.RxnMZ[combo.name] += member_F[11, 0]

                # Sum the plate forces at the node
                for plate in model.plates.values():

                    if plate.i_node == node:

                        # Get the plate's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        plate_F = plate.F(combo.name)

                        if node.support_DX: node.RxnFX[combo.name] += plate_F[0, 0]
                        if node.support_DY: node.RxnFY[combo.name] += plate_F[1, 0]
                        if node.support_DZ: node.RxnFZ[combo.name] += plate_F[2, 0]
                        if node.support_RX: node.RxnMX[combo.name] += plate_F[3, 0]
                        if node.support_RY: node.RxnMY[combo.name] += plate_F[4, 0]
                        if node.support_RZ: node.RxnMZ[combo.name] += plate_F[5, 0]

                    elif plate.j_node == node:

                        # Get the plate's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        plate_F = plate.F(combo.name)

                        if node.support_DX: node.RxnFX[combo.name] += plate_F[6, 0]
                        if node.support_DY: node.RxnFY[combo.name] += plate_F[7, 0]
                        if node.support_DZ: node.RxnFZ[combo.name] += plate_F[8, 0]
                        if node.support_RX: node.RxnMX[combo.name] += plate_F[9, 0]
                        if node.support_RY: node.RxnMY[combo.name] += plate_F[10, 0]
                        if node.support_RZ: node.RxnMZ[combo.name] += plate_F[11, 0]

                    elif plate.m_node == node:

                        # Get the plate's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        plate_F = plate.F(combo.name)

                        if node.support_DX: node.RxnFX[combo.name] += plate_F[12, 0]
                        if node.support_DY: node.RxnFY[combo.name] += plate_F[13, 0]
                        if node.support_DZ: node.RxnFZ[combo.name] += plate_F[14, 0]
                        if node.support_RX: node.RxnMX[combo.name] += plate_F[15, 0]
                        if node.support_RY: node.RxnMY[combo.name] += plate_F[16, 0]
                        if node.support_RZ: node.RxnMZ[combo.name] += plate_F[17, 0]

                    elif plate.n_node == node:

                        # Get the plate's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        plate_F = plate.F(combo.name)

                        if node.support_DX: node.RxnFX[combo.name] += plate_F[18, 0]
                        if node.support_DY: node.RxnFY[combo.name] += plate_F[19, 0]
                        if node.support_DZ: node.RxnFZ[combo.name] += plate_F[20, 0]
                        if node.support_RX: node.RxnMX[combo.name] += plate_F[21, 0]
                        if node.support_RY: node.RxnMY[combo.name] += plate_F[22, 0]
                        if node.support_RZ: node.RxnMZ[combo.name] += plate_F[23, 0]

                # Sum the quad forces at the node
                for quad in model.quads.values():

                    if quad.i_node == node:

                        # Get the quad's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        quad_F = quad.F(combo.name)

                        if node.support_DX: node.RxnFX[combo.name] += quad_F[0, 0]
                        if node.support_DY: node.RxnFY[combo.name] += quad_F[1, 0]
                        if node.support_DZ: node.RxnFZ[combo.name] += quad_F[2, 0]
                        if node.support_RX: node.RxnMX[combo.name] += quad_F[3, 0]
                        if node.support_RY: node.RxnMY[combo.name] += quad_F[4, 0]
                        if node.support_RZ: node.RxnMZ[combo.name] += quad_F[5, 0]

                    elif quad.j_node == node:

                        # Get the quad's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        quad_F = quad.F(combo.name)

                        if node.support_DX: node.RxnFX[combo.name] += quad_F[6, 0]
                        if node.support_DY: node.RxnFY[combo.name] += quad_F[7, 0]
                        if node.support_DZ: node.RxnFZ[combo.name] += quad_F[8, 0]
                        if node.support_RX: node.RxnMX[combo.name] += quad_F[9, 0]
                        if node.support_RY: node.RxnMY[combo.name] += quad_F[10, 0]
                        if node.support_RZ: node.RxnMZ[combo.name] += quad_F[11, 0]

                    elif quad.m_node == node:

                        # Get the quad's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        quad_F = quad.F(combo.name)

                        if node.support_DX: node.RxnFX[combo.name] += quad_F[12, 0]
                        if node.support_DY: node.RxnFY[combo.name] += quad_F[13, 0]
                        if node.support_DZ: node.RxnFZ[combo.name] += quad_F[14, 0]
                        if node.support_RX: node.RxnMX[combo.name] += quad_F[15, 0]
                        if node.support_RY: node.RxnMY[combo.name] += quad_F[16, 0]
                        if node.support_RZ: node.RxnMZ[combo.name] += quad_F[17, 0]

                    elif quad.n_node == node:

                        # Get the quad's global force matrix
                        # Storing it as a local variable eliminates the need to rebuild it every time a term is needed                    
                        quad_F = quad.F(combo.name)

                        if node.support_DX: node.RxnFX[combo.name] += quad_F[18, 0]
                        if node.support_DY: node.RxnFY[combo.name] += quad_F[19, 0]
                        if node.support_DZ: node.RxnFZ[combo.name] += quad_F[20, 0]
                        if node.support_RX: node.RxnMX[combo.name] += quad_F[21, 0]
                        if node.support_RY: node.RxnMY[combo.name] += quad_F[22, 0]
                        if node.support_RZ: node.RxnMZ[combo.name] += quad_F[23, 0]

                # Sum the joint loads applied to the node
                for load in node.NodeLoads:

                    for case, factor in combo.factors.items():

                        if load[2] == case:

                            if load[0] == 'FX' and node.support_DX:
                                node.RxnFX[combo.name] -= load[1]*factor
                            elif load[0] == 'FY' and node.support_DY:
                                node.RxnFY[combo.name] -= load[1]*factor
                            elif load[0] == 'FZ' and node.support_DZ:
                                node.RxnFZ[combo.name] -= load[1]*factor
                            elif load[0] == 'MX' and node.support_RX:
                                node.RxnMX[combo.name] -= load[1]*factor
                            elif load[0] == 'MY' and node.support_RY:
                                node.RxnMY[combo.name] -= load[1]*factor
                            elif load[0] == 'MZ' and node.support_RZ:
                                node.RxnMZ[combo.name] -= load[1]*factor

            # Calculate any reactions due to active spring supports at the node
            if node.spring_DX[0] is not None and node.spring_DX[2] == True:
                k = float(node.spring_DX[0])
                DX = node.DX[combo.name]
                node.RxnFX[combo.name] -= k*DX
            if node.spring_DY[0] is not None and node.spring_DY[2] == True:
                k = float(node.spring_DY[0])
                DY = node.DY[combo.name]
                node.RxnFY[combo.name] -= k*DY
            if node.spring_DZ[0] is not None and node.spring_DZ[2] == True:
                k = float(node.spring_DZ[0])
                DZ = node.DZ[combo.name]
                node.RxnFZ[combo.name] -= k*DZ
            if node.spring_RX[0] is not None and node.spring_RX[2] == True:
                k = float(node.spring_RX[0])
                RX = node.RX[combo.name]
                node.RxnMX[combo.name] -= k*RX
            if node.spring_RY[0] is not None and node.spring_RY[2] == True:
                k = float(node.spring_RY[0])
                RY = node.RY[combo.name]
                node.RxnMY[combo.name] -= k*RY
            if node.spring_RZ[0] is not None and node.spring_RZ[2] == True:
                k = float(node.spring_RZ[0])
                RZ = node.RZ[combo.name]
                node.RxnMZ[combo.name] -= k*RZ


def _check_statics(model: FEModel3D, combo_tags: List[str] | None = None) -> None:
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
        combo_list = model.load_combos.values()
    else:
        combo_list = []
        for combo in model.load_combos.values():
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
        for node in model.nodes.values():

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


def _partition_D(model: FEModel3D) -> Tuple[List[int], List[int], NDArray[float64]]:
    """Builds a list with known nodal displacements and with the positions in global stiffness matrix of known and unknown nodal displacements

    :return: A list of the global matrix indices for the unknown nodal displacements (D1_indices). A list of the global matrix indices for the known nodal displacements (D2_indices). A list of the known nodal displacements (D2).
    :rtype: list, list, list
    """

    D1_indices = []  # A list of the indices for the unknown nodal displacements
    D2_indices = []  # A list of the indices for the known nodal displacements
    D2 = []          # A list of the values of the known nodal displacements

    # Create the auxiliary table
    for node in model.nodes.values():

        # Unknown displacement DX
        if node.support_DX == False and node.EnforcedDX == None:
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


def _partition(model: FEModel3D, unp_matrix: NDArray[float64] | lil_matrix, D1_indices: List[int], D2_indices: List[int]) -> Tuple[NDArray[float64], NDArray[float64]] | Tuple[NDArray[float64], NDArray[float64], NDArray[float64], NDArray[float64]]:
    """Partitions a matrix (or vector) into submatrices (or subvectors) based on degree of freedom boundary conditions.

    :param unp_matrix: The unpartitioned matrix (or vector) to be partitioned.
    :type unp_matrix: ndarray or lil_matrix
    :param D1_indices: A list of the indices for degrees of freedom that have unknown displacements.
    :type D1_indices: list
    :param D2_indices: A list of the indices for degrees of freedom that have known displacements.
    :type D2_indices: list
    :return: Partitioned submatrices (or subvectors) based on degree of freedom boundary conditions.
    :rtype: array, array, array, array
    """

    # Determine if this is a 1D vector or a 2D matrix

    # 1D vectors
    if unp_matrix.shape[1] == 1:
        # Partition the vector into 2 subvectors
        m1 = unp_matrix[D1_indices, :]
        m2 = unp_matrix[D2_indices, :]
        return m1, m2
    # 2D matrices
    else:
        # Partition the matrix into 4 submatrices
        m11 = unp_matrix[D1_indices, :][:, D1_indices]
        m12 = unp_matrix[D1_indices, :][:, D2_indices]
        m21 = unp_matrix[D2_indices, :][:, D1_indices]
        m22 = unp_matrix[D2_indices, :][:, D2_indices]
        return m11, m12, m21, m22


def _renumber(model: FEModel3D) -> None:
    """
    Assigns node and element ID numbers to be used internally by the program. Numbers are
    assigned according to the order in which they occur in each dictionary.
    """
    
    # Number each node in the model
    for id, node in enumerate(model.nodes.values()):
        node.ID = id
    
    # Number each spring in the model
    for id, spring in enumerate(model.springs.values()):
        spring.ID = id

    # Descritize all the physical members and number each member in the model
    id = 0
    for phys_member in model.members.values():
        phys_member.descritize()
        for member in phys_member.sub_members.values():
            member.ID = id
            id += 1
    
    # Number each plate in the model
    for id, plate in enumerate(model.plates.values()):
        plate.ID = id
    
    # Number each quadrilateral in the model
    for id, quad in enumerate(model.quads.values()):
        quad.ID = id
