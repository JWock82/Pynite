# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2020 D. Craig Brinck, SE; tamalone1
"""

import unittest
from PyNite import FEModel3D
import math
import sys
from io import StringIO

class Test_2D_Frame(unittest.TestCase):
    """Tests for reaction results

    :param unittest: _description_
    :type unittest: _type_
    """

    def setUp(self):
        # Suppress printed output temporarily
        sys.stdout = StringIO()

    def tearDown(self):
        # Reset the print function to normal
        sys.stdout = sys.__stdout__

    def test_spring_reactions(self):

        model = FEModel3D()

        model.add_node('N1', 0, 0, 0)
        model.add_node_load('N1', 'FY', -5)
        model.def_support('N1', True, False, True, True, True, True)
        model.def_support_spring('N1', 'DY', '5', '-')
        model.analyze()

        self.assertEqual(model.Nodes['N1'].RxnFY['Combo 1'], 5)
