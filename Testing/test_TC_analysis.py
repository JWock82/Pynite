"""
MIT License

Copyright (c) 2020 D. Craig Brinck, SE; tamalone1
"""

import unittest
from Pynite import FEModel3D
import sys
from io import StringIO

class Test_2D_Frame(unittest.TestCase):
    """Tests for tension/compression-only analysis"""

    def setUp(self):
        # Suppress printed output temporarily
        sys.stdout = StringIO()

    def tearDown(self):
        # Reset the print function to normal
        sys.stdout = sys.__stdout__

    def test_TC_members(self):

        # Create a new finite element model
        tc_model = FEModel3D()
        tc_model.add_node('N1', 0, 0, 0)
        tc_model.add_node('N2', 100, 0, 0)
        tc_model.add_node('N3', 0, 10, 0)
        tc_model.add_node('N4', 0, -10, 0)

        E = 29000 # ksi
        G = 11400 # ksi
        nu = 0.3  # Poisson's ratio
        rho = 0.490/12**2  # Density (kci)
        tc_model.add_material('Steel', E, G, nu, rho)

        Iy = 3 # in^4
        Iz = 3 # in^4
        J = 0.0438 # in^4
        A = 1.94 # in^2
        tc_model.add_section('Section', A, Iy, Iz, J)

        tc_model.add_member('both-ways', 'N1', 'N2', 'Steel', 'Section')
        tc_model.add_member('t-only top', 'N3', 'N2', 'Steel', 'Section', tension_only=True)
        tc_model.def_releases('t-only top', Ryi=True, Rzi=True, Ryj=True, Rzj=True)
        tc_model.def_releases('both-ways', Ryi=True, Rzi=True, Ryj=True, Rzj=True)

        tc_model.def_support('N2', *[False]*2, *[True]*4)
        tc_model.def_support('N1', *[True]*6)
        tc_model.def_support('N3', *[True]*6)
        tc_model.def_support('N4', *[True]*6)
        tc_model.add_node_load('N2', 'FY', -10)

        tc_model.add_member('t-only bott', 'N4', 'N2', 'Steel', 'Section', tension_only=True)
        tc_model.def_releases('t-only bott', Ryi=True, Rzi=True, Ryj=True, Rzj=True)

        tc_model.analyze()

        self.assertAlmostEqual(tc_model.members['t-only top'].max_axial(), -100.499, 3)
        self.assertAlmostEqual(tc_model.members['both-ways'].max_axial(), 100, 3)
        self.assertEqual(tc_model.members['t-only bott'].max_axial(), 0, 3)
        self.assertFalse(tc_model.members['t-only bott'].active['Combo 1'])

    def test_inactive_member_deflection(self):
        """An inactive (slack) tension-only member should still ride along with
        its nodes: its deflection is the linear interpolation of its end-node
        displacements rather than zero (issue #317)."""

        from numpy import allclose, diff

        tc_model = FEModel3D()
        tc_model.add_node('N1', 0, 0, 0)
        tc_model.add_node('N2', 100, 0, 0)
        tc_model.add_node('N3', 0, 10, 0)
        tc_model.add_node('N4', 0, -10, 0)

        tc_model.add_material('Steel', 29000, 11400, 0.3, 0.490/12**2)
        tc_model.add_section('Section', 1.94, 3, 3, 0.0438)

        tc_model.add_member('both-ways', 'N1', 'N2', 'Steel', 'Section')
        tc_model.add_member('t-only top', 'N3', 'N2', 'Steel', 'Section', tension_only=True)
        tc_model.def_releases('t-only top', Ryi=True, Rzi=True, Ryj=True, Rzj=True)
        tc_model.def_releases('both-ways', Ryi=True, Rzi=True, Ryj=True, Rzj=True)

        tc_model.def_support('N2', *[False]*2, *[True]*4)
        tc_model.def_support('N1', *[True]*6)
        tc_model.def_support('N3', *[True]*6)
        tc_model.def_support('N4', *[True]*6)
        tc_model.add_node_load('N2', 'FY', -10)

        tc_model.add_member('t-only bott', 'N4', 'N2', 'Steel', 'Section', tension_only=True)
        tc_model.def_releases('t-only bott', Ryi=True, Rzi=True, Ryj=True, Rzj=True)

        tc_model.analyze()

        member = tc_model.members['t-only bott']
        L = member.L()

        # The bottom member goes slack (compression) and is removed from the
        # stiffness matrix, so it carries no internal force.
        self.assertFalse(member.active['Combo 1'])
        self.assertEqual(member.max_axial(), 0)

        # The member's local end displacements: the i-end (N4) is fully fixed,
        # the j-end (N2) deflects under the applied load.
        d = member._inactive_local_disp('Combo 1')
        dyi, dyj = d[1, 0], d[7, 0]

        # The j-end actually moves, so this is a non-trivial check.
        self.assertEqual(dyi, 0.0)
        self.assertNotAlmostEqual(dyj, 0.0)

        # Deflection at the ends must equal the (local) end-node displacements,
        # and the interior must be the linear interpolation between them rather
        # than zero (the old, incorrect behaviour).
        self.assertAlmostEqual(member.deflection('dy', 0.0), dyi, 9)
        self.assertAlmostEqual(member.deflection('dy', L), dyj, 9)
        self.assertAlmostEqual(member.deflection('dy', L/2), 0.5*(dyi + dyj), 9)
        self.assertNotAlmostEqual(member.deflection('dy', L/2), 0.0)

        # deflection_array must agree with deflection() and be perfectly linear.
        arr = member.deflection_array('dy', 11)
        expected = dyi + (dyj - dyi)*arr[0]/L
        self.assertTrue(allclose(arr[1], expected))
        self.assertTrue(allclose(diff(arr[1], 2), 0.0, atol=1e-12))

    def test_TC_member_reactivation(self):
        """A tension-only member deactivated on one iteration is reactivated when the rest of
        the structure pulls it taut again.

        A node in the XY plane hangs from an ordinary bar straight up and from three
        tension-only ties, ten times stiffer, towards 0, 135 and 225 degrees, and is pulled
        right and down. With every member active the ties at 0 and 225 degrees are both
        compressed and both deactivated. Without them the node moves right, away from the
        anchor of the tie at 0 degrees, which is stretched and has to come back in: the answer
        is the bar and the ties at 0 and 135 degrees, the tie at 225 degrees slack.
        """
        from math import cos, sin, radians, isclose

        L = 2.0
        ties = {'T0': 0.0, 'T135': 135.0, 'T225': 225.0}

        def model(ties, tension_only=True):
            m = FEModel3D()
            m.add_material('Steel', 200e6, 77e6, 0.3, 0.0)
            m.add_section('Bar', 1e-4, 1e-8, 1e-8, 1e-8)
            m.add_section('Tie', 1e-3, 1e-8, 1e-8, 1e-8)
            m.add_node('N', 0, 0, 0)
            m.def_support('N', False, False, True, True, True, True)
            members = {'Bar': (90.0, 'Bar', False)}
            members.update({name: (angle, 'Tie', tension_only) for name, angle in ties.items()})
            for name, (angle, section, t_only) in members.items():
                m.add_node('A' + name, L*cos(radians(angle)), L*sin(radians(angle)), 0)
                m.def_support('A' + name, True, True, True, True, True, True)
                m.add_member(name, 'N', 'A' + name, 'Steel', section, tension_only=t_only)
                m.def_releases(name, Ryi=True, Rzi=True, Ryj=True, Rzj=True)
            m.add_node_load('N', 'FX', 1.0)
            m.add_node_load('N', 'FY', -2.0)
            return m

        tc_model = model(ties)
        tc_model.analyze(check_statics=False)
        active = {name: tc_model.members[name].active['Combo 1'] for name in ties}
        self.assertEqual(active, {'T0': True, 'T135': True, 'T225': False})

        # The same displacements as a linear model of the members that stay active
        linear = model({'T0': 0.0, 'T135': 135.0}, tension_only=False)
        linear.analyze_linear(check_statics=False)
        for d in ('DX', 'DY'):
            self.assertTrue(isclose(getattr(tc_model.nodes['N'], d)['Combo 1'],
                                    getattr(linear.nodes['N'], d)['Combo 1'], rel_tol=1e-9))

        # The active ties pull (compression is positive), and the slack one is not stretched
        for name in ('T0', 'T135'):
            self.assertLess(tc_model.members[name].max_axial('Combo 1'), 0)
        n = tc_model.nodes['N']
        stretch = -(n.DX['Combo 1']*cos(radians(225)) + n.DY['Combo 1']*sin(radians(225)))
        self.assertLessEqual(stretch, 0)
