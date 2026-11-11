"""Tests for `Analysis._unpartition`, which rebuilds a global vector from its partitioned pieces."""

# Standard library imports
import random

# Third party imports
from numpy import array, zeros, allclose

# Pynite imports
from Pynite import FEModel3D
from Pynite import Analysis


def _reference_unpartition(num_nodes, V1, V2, V1_indices, V2_indices):
    """The original loop-based implementation, kept here as the reference to compare against."""

    # Start with a global vector of zeros
    V = zeros((num_nodes*6, 1))

    # Step through every global dof and pull its term from whichever partition holds it
    for dof in range(num_nodes*6):
        if dof in V2_indices:
            V[dof, 0] = V2[V2_indices.index(dof), 0]
        else:
            V[dof, 0] = V1[V1_indices.index(dof), 0]

    # Return the rebuilt vector
    return V


def _model_with_nodes(num_nodes):
    """Builds a model holding `num_nodes` free nodes, which is all `_unpartition` needs."""

    # Create the model and add the nodes along a line
    model = FEModel3D()
    for i in range(num_nodes):
        model.add_node(f'N{i + 1}', i, 0, 0)

    # Return the model
    return model


def test_unpartition_matches_reference():
    """The vectorized version must give exactly the same vector as the original loop."""

    # Use a fixed seed so the test is repeatable
    rng = random.Random(42)

    # Build a model with 50 nodes (300 dofs)
    num_nodes = 50
    model = _model_with_nodes(num_nodes)

    # Shuffle the dofs and split them into unknown (V1) and known (V2) partitions
    dofs = list(range(num_nodes*6))
    rng.shuffle(dofs)
    V2_indices = sorted(dofs[:40])
    V1_indices = sorted(dofs[40:])

    # Give every dof a distinct value so any misplaced term would be caught
    V1 = array([[float(i + 1)] for i in range(len(V1_indices))])
    V2 = array([[-float(i + 1)] for i in range(len(V2_indices))])

    # Compare the new implementation against the reference
    new = Analysis._unpartition(model, V1, V2, V1_indices, V2_indices)
    ref = _reference_unpartition(num_nodes, V1, V2, V1_indices, V2_indices)
    assert new.shape == (num_nodes*6, 1)
    assert allclose(new, ref, rtol=0, atol=0)


def test_unpartition_empty_partitions():
    """Either partition may be empty: all dofs restrained, or no dof restrained."""

    # Build a small model
    num_nodes = 3
    model = _model_with_nodes(num_nodes)
    all_dofs = list(range(num_nodes*6))
    values = array([[float(i)] for i in all_dofs])

    # Every dof known (V1 empty, as when the stiffness matrix partition K11 is empty)
    V = Analysis._unpartition(model, [], values, [], all_dofs)
    assert allclose(V, values)

    # No dof known (V2 empty)
    V = Analysis._unpartition(model, values, zeros((0, 1)), all_dofs, [])
    assert allclose(V, values)
