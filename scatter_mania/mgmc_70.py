import numpy as np

import openmc


def create_random_ray_model():
    ###############################################################################
    # Create multigroup data

    ###############################################################################
    # Create materials for the problem

    # Instantiate some Macroscopic Data
    void_data   = openmc.Macroscopic('heavy water')
    shield_data = openmc.Macroscopic('borated water')

    void_mat = openmc.Material(name='heavy water')
    void_mat.set_density('macro', 1.0)
    void_mat.add_macroscopic(void_data)
    
    shield_mat = openmc.Material(name='borated water')
    shield_mat.set_density('macro', 1.0)
    shield_mat.add_macroscopic(shield_data)

    # Instantiate a Materials collection and export to XML
    materials_file = openmc.Materials([void_mat, shield_mat])
    materials_file.cross_sections = "mgxs.h5"

    ###############################################################################
    # Define problem geometry
    
    source_cell = openmc.Cell(fill=void_mat, name='infinite source region')
    void_cell = openmc.Cell(fill=void_mat, name='infinite void region')
    shield_cell = openmc.Cell(fill=shield_mat, name='infinite shield region')
    
    su = openmc.Universe()
    su.add_cells([source_cell])
    
    vu = openmc.Universe()
    vu.add_cells([void_cell])
    
    au = openmc.Universe()
    au.add_cells([shield_cell])

    z_base = [
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [vu, vu, vu, vu, au, au],
            [vu, au, au, au, au, au],
            [vu, au, au, au, au, au],
            [vu, au, au, au, au, au],
            [vu, au, au, au, au, au],
            [su, au, au, au, au, au]
                         ]
    
    z_col = [
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, vu, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au]
                         ]
    
    z_high = [
            [au, au, au, vu, au, au],
            [au, au, au, vu, au, au],
            [au, au, au, vu, au, au],
            [au, au, au, vu, au, au],
            [au, au, au, vu, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au]
                         ]
    
    z_cap = [
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au],
            [au, au, au, au, au, au]
                         ]

    dogleg_pattern = [
            z_base,
            z_col,
            z_col,
            z_high,
            z_cap,
            z_cap
            ]
    
    x = 60.0
    x_dim = 6

    y = 100.0
    y_dim = 10

    z = 60.0
    z_dim = 6
    
    lattice = openmc.RectLattice()
    lattice.lower_left = [0.0, 0.0, 0.0]
    lattice.pitch = [x/x_dim, y/y_dim, z/z_dim]
    lattice.universes = dogleg_pattern
    
    lattice_cell = openmc.Cell(fill=lattice)

    lattice_uni = openmc.Universe()
    lattice_uni.add_cells([lattice_cell])

    x_low  = openmc.XPlane(x0=0.0,boundary_type='reflective') 
    x_high = openmc.XPlane(x0=x,boundary_type='vacuum') 
    
    y_low  = openmc.YPlane(y0=0.0,boundary_type='reflective') 
    y_high = openmc.YPlane(y0=y,boundary_type='vacuum') 
    
    z_low  = openmc.ZPlane(z0=0.0,boundary_type='reflective') 
    z_high = openmc.ZPlane(z0=z,boundary_type='vacuum') 
    
    full_domain = openmc.Cell(fill=lattice_uni, region=+x_low & -x_high & +y_low & -y_high & +z_low & -z_high, name='full domain')

    root = openmc.Universe(name='root universe')
    root.add_cell(full_domain)

    # Create a geometry with the two cells and export to XML
    geometry = openmc.Geometry(root)

    ###############################################################################
    # Define problem settings

    # Instantiate a Settings object, set all runtime parameters, and export to XML
    settings = openmc.Settings()
    settings.energy_mode = "multi-group"
    settings.batches = 100
    #settings.inactive = 5
    settings.particles = 100000
    settings.run_mode = 'fixed source'

    #settings.random_ray_distance_active = 100.0
    #settings.random_ray_distance_inactive = 20.0
    #settings.solver_type = 'random ray'
    #settings.run_mode = 'fixed source'

    # Create an initial uniform spatial source for ray integration
    #lower_left = (-pitch, -pitch, -1)
    #upper_right = (pitch, pitch, 1)
    #uniform_dist = openmc.stats.Box(lower_left, upper_right, only_fissionable=False)
    #rr_source = openmc.IndependentSource(space=uniform_dist, particle="random_ray")
    
    # Create the neutron source in the bottom right of the moderator
    strengths = [1.0] # Good - fast group appears largest (besides most thermal)
    midpoints = [7.5e6]
    energy_distribution = openmc.stats.Discrete(x=midpoints,p=strengths)
    
    lower_left_src = [0.0, 0.0, 0.0]
    upper_right_src = [10.0, 10.0, 10.0]
    spatial_distribution = openmc.stats.Box(lower_left_src, upper_right_src, only_fissionable=False)

    #source = openmc.IndependentSource(energy=energy_distribution, domains=[source_mat], strength=2.0) # works
    source = openmc.IndependentSource(space=spatial_distribution, energy=energy_distribution,  strength=1.0) # works
    
    #settings.source = [source, rr_source]
    settings.source = [source]
    #settings.export_to_xml()

    ###############################################################################
    # Define tallies

    # Create a mesh that will be used for tallying
    #mesh = openmc.RegularMesh()
    #mesh.dimension = (x_dim, y_dim, z_dim)
    #mesh.lower_left = (0.0, 0.0, 0.0)
    #mesh.upper_right = (x, y, z)

    # Create a mesh filter that can be used in a tally
    #mesh_filter = openmc.MeshFilter(mesh)

    # Now use the mesh filter in a tally and indicate what scores are desired
    #tally = openmc.Tally(name="Mesh tally")
    #tally.filters = [mesh_filter]
    #tally.scores = ['flux']
    #tally.estimator = 'collision'
    #tally.estimator = 'analog'

    estimator = 'tracklength'

    # Case 3A
    mesh_3A = openmc.RegularMesh()
    mesh_3A.dimension = (1, y_dim, 1)
    mesh_3A.lower_left = (0.0, 0.0, 0.0)
    mesh_3A.upper_right = (10.0, y, 10.0)
    mesh_filter_3A = openmc.MeshFilter(mesh_3A)
    
    tally_3A = openmc.Tally(name="Case 3A")
    tally_3A.filters = [mesh_filter_3A]
    tally_3A.scores = ['flux']
    tally_3A.estimator = estimator
    
    # Case 3B
    mesh_3B = openmc.RegularMesh()
    mesh_3B.dimension = (x_dim, 1, 1)
    mesh_3B.lower_left = (0.0, 50.0, 0.0)
    mesh_3B.upper_right = (x, 60.0, 10.0)
    mesh_filter_3B = openmc.MeshFilter(mesh_3B)
    
    tally_3B = openmc.Tally(name="Case 3B")
    tally_3B.filters = [mesh_filter_3B]
    tally_3B.scores = ['flux']
    tally_3B.estimator = estimator
    
    # Case 3C
    mesh_3C = openmc.RegularMesh()
    mesh_3C.dimension = (x_dim, 1, 1)
    mesh_3C.lower_left = (0.0, 90.0, 30.0)
    mesh_3C.upper_right = (x, 100.0, 40.0)
    mesh_filter_3C = openmc.MeshFilter(mesh_3C)
    
    tally_3C = openmc.Tally(name="Case 3C")
    tally_3C.filters = [mesh_filter_3C]
    tally_3C.scores = ['flux']
    tally_3C.estimator = estimator

    # Instantiate a Tallies collection and export to XML
    tallies = openmc.Tallies([tally_3A, tally_3B, tally_3C])


    ###############################################################################
    #                   Exporting to OpenMC plots.xml file
    ###############################################################################

    plot = openmc.Plot()
    plot.origin = [x/2.0, y/2.0, z/2.0]
    plot.width = [x, y, z]
    plot.pixels = [60, 100, 60]
    plot.type = 'voxel'

    # Instantiate a Plots collection and export to XML
    plot_file = openmc.Plots([plot])
    #plot_file.export_to_xml()
    
    model = openmc.model.Model()
    model.geometry = geometry
    model.materials = materials_file
    model.settings = settings
    #model.xs_data = mg_cross_sections_file
    model.tallies = tallies
    model.plots = plot_file

    return model

model = create_random_ray_model()
model.export_to_model_xml()
