==========
Quickstart
==========

Simple Beam
-----------

Here's a simple example of how to analyze a simple beam. Many more examples are available here:
`Examples <https://github.com/JWock82/Pynite/tree/main/Examples>`_.

.. code-block:: python

    # Example of a simply supported beam with a uniform distributed load.
    # Units used in this example are inches and kips
    # This example does not use load combinations. The program will create a
    # default load combindation called 'Combo 1'

    # Import `FEModel3D` from `Pynite`
    from Pynite import FEModel3D

    # Create a new finite element model
    beam = FEModel3D()

    # Add nodes (14 ft = 168 inches apart)
    beam.add_node('N1', 0, 0, 0)
    beam.add_node('N2', 168, 0, 0)

    # Define a material
    E = 29000       # Modulus of elasticity (ksi)
    G = 11200       # Shear modulus of elasticity (ksi)
    nu = 0.3        # Poisson's ratio
    rho = 2.836e-4  # Density (kci)
    beam.add_material('Steel', E, G, nu, rho)

    # Add a section with the following properties:
    # Iy = 100 in^4, Iz = 150 in^4, J = 250 in^4, A = 20 in^2
    beam.add_section('MySection', 20, 100, 150, 250)

    #Add a member
    beam.add_member('M1', 'N1', 'N2', 'Steel', 'MySection')

    # Provide simple supports
    beam.def_support('N1', True, True, True, False, False, False)
    beam.def_support('N2', True, True, True, True, False, False)

    # Add a uniform load of 200 lbs/ft to the beam (from 0 in to 168 in)
    beam.add_member_dist_load('M1', 'Fy', -200/1000/12, -200/1000/12, 0, 168)

    # Alternatively the following line would do apply the load to the full
    # length of the member as well
    # beam.add_member_dist_load('M1', 'Fy', -200/1000/12, -200/1000/12)

    # Analyze the beam
    beam.analyze()

    # Print the shear, moment, and deflection diagrams
    beam.members['M1'].plot_shear('Fy')
    beam.members['M1'].plot_moment('Mz')
    beam.members['M1'].plot_deflection('dy')

    # Print reactions at each end of the beam
    print(f"Left Support Reaction: { {k: float(v) for k, v in beam.nodes['N1'].RxnFY.items()} }")
    print(f"Right Support Reaction: { {k: float(v) for k, v in beam.nodes['N2'].RxnFY.items()} }")

    # Render the deformed shape of the beam magnified 100 times, with a text
    # height of 5 inches
    from Pynite.Visualization import Renderer
    renderer = Renderer(beam)
    renderer.annotation_size = 6
    renderer.deformed_shape = True
    renderer.deformed_scale = 100
    renderer.render_loads = True
    renderer.render_model()

3D Frame
--------

This example builds a simple portal frame in 3D: two columns and a beam, with
a lateral point load at the top. It demonstrates load combinations, support
definitions, and extracting member forces — a typical starting point for frame
analysis.

.. code-block:: python

    from Pynite import FEModel3D

    # Create the model
    frame = FEModel3D()

    # Nodes: two column bases and two column tops (units: inches, kips)
    #   N1 ---M3--- N2      (beam at top, 20 ft span)
    #   |           |
    #   M1          M2      (columns, 12 ft tall)
    #   |           |
    #   N3          N4      (fixed bases)
    frame.add_node('N1', 0, 12*12, 0)
    frame.add_node('N2', 20*12, 12*12, 0)
    frame.add_node('N3', 0, 0, 0)
    frame.add_node('N4', 20*12, 0, 0)

    # Material: structural steel
    frame.add_material('Steel', 29000, 11200, 0.3, 2.836e-4)

    # Section: W10x33 (approximate properties)
    # A=9.71 in^2, Iy=36.6 in^4, Iz=171 in^4, J=0.583 in^4
    frame.add_section('W10x33', 9.71, 36.6, 171, 0.583)

    # Members: two columns and one beam
    frame.add_member('M1', 'N3', 'N1', 'Steel', 'W10x33')
    frame.add_member('M2', 'N4', 'N2', 'Steel', 'W10x33')
    frame.add_member('M3', 'N1', 'N2', 'Steel', 'W10x33')

    # Fixed supports at column bases
    frame.def_support('N3', True, True, True, True, True, True)
    frame.def_support('N4', True, True, True, True, True, True)

    # Loads
    # Dead load: uniform gravity on beam (500 plf = 0.0417 kip/in)
    frame.add_member_dist_load('M3', 'Fy', -0.5/12, -0.5/12, case='D')
    # Wind load: lateral point load at top-left node
    frame.add_node_load('N1', 'FX', 10, case='W')

    # Load combinations (ASCE 7 LRFD)
    frame.add_load_combo('1.4D', {'D': 1.4})
    frame.add_load_combo('1.2D+1.0W', {'D': 1.2, 'W': 1.0})

    # Run analysis
    frame.analyze(check_statics=True)

    # Print base reactions for the wind combination
    print('--- Base reactions (1.2D+1.0W) ---')
    for node_name in ['N3', 'N4']:
        node = frame.nodes[node_name]
        Rx = node.RxnFX['1.2D+1.0W']
        Ry = node.RxnFY['1.2D+1.0W']
        Mz = node.RxnMZ['1.2D+1.0W']
        print(f'{node_name}: Rx={Rx:.2f} k, Ry={Ry:.2f} k, Mz={Mz:.1f} k-in')

    # Print peak beam forces
    beam = frame.members['M3']
    print(f'\nBeam max moment: {beam.max_moment("Mz", "1.2D+1.0W")/12:.1f} k-ft')
    print(f'Beam max shear:  {beam.max_shear("Fy", "1.2D+1.0W"):.2f} k')
    print(f'Beam max defl:   {beam.min_deflection("dy", "1.2D+1.0W"):.4f} in')

    # Render the deformed shape
    from Pynite.Visualization import Renderer
    rndr = Renderer(frame)
    rndr.deformed_shape = True
    rndr.deformed_scale = 50
    rndr.render_loads = True
    rndr.combo_name = '1.2D+1.0W'
    rndr.render_model()
