from Pynite import FEModel3D
import numpy as np
import pytest


def test_spring_elements():
    # A First Course in the Finite Element Method, 4th Edition
    # Daryl L. Logan
    # Example 2.1
    # Units for this model are pounds and inches

    system = FEModel3D()

    system.add_node('1', 0, 0, 0)
    system.add_node('2', 30, 0, 0)
    system.add_node('3', 10, 0, 0)
    system.add_node('4', 20, 0, 0)

    # Add spring members
    system.add_spring('S1', '1', '3', 1000)
    system.add_spring('S2', '3', '4', 2000)
    system.add_spring('S3', '4', '2', 3000)

    # Define supports
    system.def_support('1', True, True, True, True, True, True)
    system.def_support('2', True, True, True, True, True, True)
    system.def_support('3', False, True, True, True, True, True)
    system.def_support('4', False, True, True, True, True, True)

    # Add node loads
    system.add_node_load('4', 'FX', 5000)

    system.analyze(True)

    # Check the elastic stiffness matrix
    Ke = system.Ke(sparse=False)[[0, 6, 12, 18]][:, [0, 6, 12, 18]]

    assert np.array_equal(Ke, np.array([[1000, 0, -1000, 0],
                                       [0, 3000, 0, -3000],
                                       [-1000, 0, 3000, -2000],
                                       [0, -3000, -2000, 5000]]))

    # Check results against the problem's known solution
    n3_DX = system.nodes['3'].DX['Combo 1']
    assert round(n3_DX, 4) == round(10/11, 4), 'Failed spring displacement test.'

    n4_DX = system.nodes['4'].DX['Combo 1']
    assert round(n4_DX, 4) == round(15/11, 4), 'Failed spring displacement test.'

    n1_rxn = system.nodes['1'].RxnFX['Combo 1']
    assert round(n1_rxn, 0) == round(-10000/11, 0), 'Failed node reaction test.'

    n2_rxn = system.nodes['2'].RxnFX['Combo 1']
    assert round(n2_rxn, 0) == round(-45000/11, 0), 'Failed node reaction test'

def test_nodal_springs():

    model = FEModel3D()

    model.add_node('N1', 0, 0, 0)

    model.def_support_spring('N1', 'DX', 1000)
    model.def_support_spring('N1', 'DY', 2000)
    model.def_support_spring('N1', 'DZ', 3000)
    model.def_support_spring('N1', 'RX', 4000)
    model.def_support_spring('N1', 'RY', 5000)
    model.def_support_spring('N1', 'RZ', 6000)

    model.add_node_load('N1', 'FX', -1000)
    model.add_node_load('N1', 'FY', -2000)
    model.add_node_load('N1', 'FZ', -3000)
    model.add_node_load('N1', 'MX', -4000)
    model.add_node_load('N1', 'MY', -5000)
    model.add_node_load('N1', 'MZ', -6000)

    model.analyze(check_statics=True)

    assert model.nodes['N1'].DX['Combo 1'] == -1, 'Faied nodal spring displacement test in X direction.'
    assert model.nodes['N1'].DY['Combo 1'] == -1, 'Faied nodal spring displacement test in Y direction.'
    assert model.nodes['N1'].DZ['Combo 1'] == -1, 'Faied nodal spring displacement test in Z direction.'
    assert model.nodes['N1'].RX['Combo 1'] == -1, 'Faied nodal spring displacement test in RX direction.'
    assert model.nodes['N1'].RY['Combo 1'] == -1, 'Faied nodal spring displacement test in RY direction.'
    assert model.nodes['N1'].RZ['Combo 1'] == -1, 'Faied nodal spring displacement test in RZ direction.'

@pytest.mark.parametrize('analyze_first', [False, True])
def test_delete_node_removes_connected_springs(analyze_first):
    """Remove springs at either end of a deleted node and preserve an unrelated spring."""
    model = FEModel3D()
    for name, x in [('A', 0), ('C', 2), ('B', 1)]:
        model.add_node(name, x, 0, 0)
        model.def_support(name, name == 'A', True, True, True, True, True)
    model.add_spring('AB', 'A', 'B', 100)
    model.add_spring('BC', 'B', 'C', 100)
    model.add_spring('AC', 'A', 'C', 100)
    model.add_node_load('C', 'FX', 10)
    if analyze_first:
        model.analyze_linear()
        assert model.nodes['C'].DX['Combo 1'] == pytest.approx(10 / 150)

    # B is the j-node of AB and the i-node of BC; AC must remain usable.
    model.delete_node('B')
    assert set(model.springs) == {'AC'}
    model.analyze_linear()
    assert model.nodes['C'].DX['Combo 1'] == pytest.approx(10 / 100)
    assert model.nodes['A'].RxnFX['Combo 1'] == pytest.approx(-10)


if __name__ == '__main__':
    test_spring_elements()
    test_nodal_springs()
