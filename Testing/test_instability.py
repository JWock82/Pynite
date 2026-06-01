# Tests for detection of global instability (a singular stiffness matrix).
#
# These cover GitHub issues #255 ("No error is thrown on unstable model in some scenarios") and
# #275 ("Model produces results when unstable with a hinge"), where a globally unstable model
# (a mechanism or a structure with insufficient supports) was silently solved and returned
# erroneous reactions/displacements instead of raising. The nodal stability check only catches
# degrees of freedom with zero diagonal stiffness; these models are globally rank-deficient while
# every diagonal term is non-zero.

import unittest

from Pynite import FEModel3D


class TestInstabilityDetection(unittest.TestCase):
    """Globally unstable models must raise, while stable models must still solve."""

    def test_hinge_mechanism_raises_sparse(self):
        # Issue #275: a two-span beam joined by an internal hinge with no support at the hinge is a
        # mechanism. Both the dense and sparse solvers previously returned erroneous finite results.
        model = self._hinge_mechanism()
        with self.assertRaises(Exception):
            model.analyze(sparse=True, check_statics=False)

    def test_hinge_mechanism_raises_dense(self):
        model = self._hinge_mechanism()
        with self.assertRaises(Exception):
            model.analyze(sparse=False, check_statics=False)

    def test_unsupported_member_raises(self):
        # Issue #255, scenario 1: a single member with no supports (six rigid body modes).
        for sparse in (True, False):
            model = FEModel3D()
            model.add_node('N1', 0, 0, 0)
            model.add_node('N2', 0, 1, 0)
            model.add_material('Steel', 29000, 11400, 0.5, 490 / 1000 / 12**3)
            model.add_section('Section', 10, 100, 150, 250)
            model.add_member('M1', 'N1', 'N2', 'Steel', 'Section')
            with self.assertRaises(Exception):
                model.analyze(sparse=sparse, check_statics=False)

    def test_unsupported_frame_raises(self):
        # Issue #255, scenario 2: an unsupported frame loaded in-plane. The sparse solver previously
        # returned finite-but-meaningless displacements with no warning at all.
        for sparse in (True, False):
            model = FEModel3D()
            model.add_node('N1', 0, 0, 0)
            model.add_node('N2', 0, 144, 0)
            model.add_node('N3', 180, 144, 0)
            model.add_node('N4', 180, 0, 0)
            model.add_material('Steel', 29000, 11400, 0.5, 490 / 1000 / 12**3)
            model.add_section('Section', 10, 100, 150, 250)
            model.add_section('Section2', 15, 100, 250, 250)
            model.add_member('M1', 'N1', 'N2', 'Steel', 'Section')
            model.add_member('M2', 'N4', 'N3', 'Steel', 'Section')
            model.add_member('M3', 'N2', 'N3', 'Steel', 'Section2')
            model.add_node_load('N2', 'FX', 50)
            with self.assertRaises(Exception):
                model.analyze(sparse=sparse, check_statics=False)

    def test_stable_model_still_solves(self):
        # A stable propped cantilever must analyze normally and give the correct reaction (the fixed
        # end of a propped cantilever under a uniform load w over length L carries 5*w*L/8).
        for sparse in (True, False):
            model = self._propped_cantilever()
            model.analyze(sparse=sparse, check_statics=False)
            self.assertAlmostEqual(model.nodes['n1'].RxnFY['Combo 1'], 6.25, places=4)

    def test_check_stability_false_skips_check(self):
        # With ``check_stability=False`` the stability check is bypassed (the historical fast path),
        # so an unstable model does not raise from the stability check.
        model = self._hinge_mechanism()
        try:
            model.analyze(sparse=True, check_stability=False, check_statics=False)
        except Exception:
            # It is acceptable for the bare solver to fail on its own, but it must NOT be blocked by
            # our stability check. Re-running the stable model with the flag off must succeed.
            pass
        stable = self._propped_cantilever()
        stable.analyze(sparse=True, check_stability=False, check_statics=False)
        self.assertAlmostEqual(stable.nodes['n1'].RxnFY['Combo 1'], 6.25, places=4)

    @staticmethod
    def _hinge_mechanism():
        model = FEModel3D()
        model.add_section('W10x33', 0.006264504, 0.000015234, 0.0000712, 0.000000241)
        model.add_material('steel', 200000000, 77000000, 0.3, 77)
        model.add_node('n1', 0, 0, 0)
        model.add_node('n2', 1, 0, 0)
        model.add_node('n3', 2, 0, 0)
        model.add_member('m1', 'n1', 'n2', 'steel', 'W10x33')
        model.add_member('m2', 'n2', 'n3', 'steel', 'W10x33')
        model.def_support('n1', 1, 1, 1, 1, 1, 0)
        model.def_support('n3', 0, 1, 1, 1, 1, 0)
        model.add_member_dist_load('m1', 'Fy', 10, 10)
        model.add_member_dist_load('m2', 'Fy', 10, 10)
        model.def_releases('m2', Rzi=True)
        return model

    @staticmethod
    def _propped_cantilever():
        model = FEModel3D()
        model.add_node('n1', 0, 0, 0)
        model.add_node('n2', 10, 0, 0)
        model.add_material('steel', 29000, 11200, 0.3, 0.284)
        model.add_section('S', 10, 100, 150, 250)
        model.add_member('m', 'n1', 'n2', 'steel', 'S')
        model.def_support('n1', 1, 1, 1, 1, 1, 1)
        model.def_support('n2', 0, 1, 1, 0, 0, 0)
        model.add_member_dist_load('m', 'Fy', -1, -1)
        return model


if __name__ == '__main__':
    unittest.main()
