# -*- coding: utf-8 -*-
from Pynite import FEModel3D
import numpy as np


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

    # Check the stiffness matrix
    K = system.K(sparse=False)[[0, 6, 12, 18]][:, [0, 6, 12, 18]]

    assert np.array_equal(K, np.array([[1000, 0, -1000, 0],
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


if __name__ == '__main__':
    test_spring_elements()
