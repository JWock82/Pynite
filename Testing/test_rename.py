"""Tests for renaming model objects without overwriting existing names."""

import pytest

from Pynite import FEModel3D


@pytest.mark.parametrize('kind, prefix', [
    ('nodes', 'N'), ('springs', 'S'), ('members', 'M'),
    ('plates', 'P'), ('quads', 'Q'),
])
def test_rename_preserves_objects_with_overlapping_names(kind, prefix):
    model = FEModel3D()
    model.add_material('Material', 1000, 400, 0.25, 1)
    model.add_section('Section', 1, 1, 1, 1)
    if kind == 'nodes':
        model.add_node('N2', 0, 0, 0)
        model.add_node('N1', 1, 0, 0)
    else:
        for name, x, y in [('A', 0, 0), ('B', 1, 0), ('C', 1, 1), ('D', 0, 1)]:
            model.add_node(name, x, y, 0)
        for name in (prefix + '2', prefix + '1'):
            if kind == 'springs':
                model.add_spring(name, 'A', 'B', 100)
            elif kind == 'members':
                model.add_member(name, 'A', 'B', 'Material', 'Section')
            elif kind == 'plates':
                model.add_plate(name, 'A', 'B', 'C', 'D', 0.1, 'Material')
            else:
                model.add_quad(name, 'A', 'B', 'C', 'D', 0.1, 'Material')
    objects = getattr(model, kind)
    originals = list(objects.values())

    model.rename()
    assert list(objects) == [prefix + '1', prefix + '2']
    assert list(objects.values()) == originals
    assert [obj.name for obj in originals] == list(objects)
    # Already-numbered names must remain safe on a repeated rename.
    model.rename()
    assert list(objects.values()) == originals
    assert [obj.name for obj in originals] == list(objects)


def test_renamed_spring_model_still_solves():
    model = FEModel3D()
    model.add_node('N2', 0, 0, 0)
    model.add_node('N1', 1, 0, 0)
    model.add_spring('S2', 'N2', 'N1', 100)
    model.def_support('N2', True, True, True, True, True, True)
    model.def_support('N1', False, True, True, True, True, True)
    model.add_node_load('N1', 'FX', 10)

    model.rename()
    model.analyze_linear()
    # The tip displacement is F/k and the fixed reaction balances the load.
    assert model.nodes['N2'].DX['Combo 1'] == pytest.approx(0.1)
    assert model.nodes['N1'].RxnFX['Combo 1'] == pytest.approx(-10)
