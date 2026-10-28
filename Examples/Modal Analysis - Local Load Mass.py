"""Convert equivalent local/global member loads to mass on an inclined cantilever.

Run from the repository: uv run --locked python "Examples/Modal Analysis - Local Load Mass.py"
Units: N, m, kg, s. This example checks coordinate representation consistency.
"""

import numpy as np

from Pynite import FEModel3D


def cantilever(local_loads):
    """Build the same physical model using either local or global force components."""
    model = FEModel3D()
    model.add_node('Fixed', 0, 0, 0)
    model.add_node('Free', 8, 6, 0)
    model.def_support('Fixed', True, True, True, True, True, True)
    # Density is entered as weight per volume, consistent with member self-weight.
    model.add_material('Steel', 200e9, 80e9, 0.3, 7800*9.81)
    model.add_section('Beam', 0.01, 8.33e-6, 8.33e-6, 1.67e-5)
    model.add_member('M', 'Fixed', 'Free', 'Steel', 'Beam')

    # The local y unit axis is (-0.6, 0.8, 0), obtained from the 8-by-6 geometry.
    # A -1000 N local y load therefore has global components (600, -800, 0) N.
    loads = [('Fy', -1000)] if local_loads else [('FX', 600), ('FY', -800)]
    for direction, force in loads:
        model.add_member_pt_load('M', direction, force, 5, case='Equipment')

    # Include the member's consistent self-weight mass in both models.
    model.add_member_self_weight('FY', -1, case='Self weight')
    model.add_load_combo('Mass', {'Equipment': 1.5, 'Self weight': 1})
    model.analyze_modal(num_modes=2, mass_combo_name='Mass', mass_direction='Y', gravity=9.81)
    return model


if __name__ == '__main__':
    local_model = cantilever(local_loads=True)
    global_model = cantilever(local_loads=False)

    # The equipment load contributes |1.5*(-800)|/9.81 kg in global Y.
    equipment_mass, location = local_model.members['M']._calc_load_mass('Mass', 'Y', 9.81)
    np.testing.assert_allclose(equipment_mass, 1200/9.81)
    np.testing.assert_allclose(location, 5)
    np.testing.assert_allclose(local_model.frequencies, global_model.frequencies, rtol=1e-9)
    print(f'Equipment mass from global Y: {equipment_mass:.6f} kg at x = {location:.3f} m')
    print(f'Local load frequencies:  {local_model.frequencies} Hz')
    print(f'Global load frequencies: {global_model.frequencies} Hz')
