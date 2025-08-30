# -*- coding: utf-8 -*-

from Pynite import FEModel3D


def test_spring_reactions():

    model = FEModel3D()

    model.add_node('N1', 0, 0, 0)
    model.add_node_load('N1', 'FY', -5)
    model.def_support('N1', True, False, True, True, True, True)
    model.def_support_spring('N1', 'DY', '5', '-')
    model.analyze()

    assert model.nodes['N1'].RxnFY['Combo 1'] == 5, 'Failed spring reaction test.'


if __name__ == '__main__':
    test_spring_reactions()
