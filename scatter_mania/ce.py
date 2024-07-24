import numpy as np

import openmc

def fill_3d_list(n, val):
    """
    Generates a 3D list of dimensions nxnxn filled with copies of val.

    Parameters:
    n (int): The dimension of the 3D list.
    val (any): The value to fill the 3D list with.

    Returns:
    list: A 3D list of dimensions nxnxn filled with val.
    """
    return [[[val for _ in range(n)] for _ in range(n)] for _ in range(n)]

def create_random_ray_model():
    ###############################################################################
    # Create multigroup data

    heavy_water_mat = openmc.Material(name='heavy water')

    # Set the density of heavy water, approximately 1.11 g/cm3 at room temperature
    heavy_water_mat.set_density('g/cm3', 1.11)

    # Add deuterium (Hydrogen-2)
    heavy_water_mat.add_nuclide('H2', 2.0 / 20.027, 'ao')  # Deuterium atomic fraction

    # Add Oxygen-16
    heavy_water_mat.add_nuclide('O16', 1.0 / 20.027, 'ao')  # Oxygen-16 atomic fraction


    air_mat = openmc.Material(name='air')

    # Set the density of air_mat, approximately 0.001225 g/cm3 at sea level and 15 degrees Celsius
    air_mat.set_density('g/cm3', 0.001225)

    # Add nitrogen (N2)
    air_mat.add_element('N', 0.78, 'ao')  # 78% nitrogen by volume

    # Add oxygen (O2)
    air_mat.add_element('O', 0.21, 'ao')  # 21% oxygen by volume

    # Add argon (Ar)
    air_mat.add_element('Ar', 0.01, 'ao')  # 1% argon by volume


    borated_water_mat = openmc.Material(name='borated water')

    # Set the density of borated water, approximately 1.00 g/cm3 (similar to regular water)
    borated_water_mat.set_density('g/cm3', 1.00)

    # Add regular light water (H2O)
    # Add hydrogen
    borated_water_mat.add_nuclide('H1', 2.0, 'ao')  # Two hydrogen atoms per molecule of water

    # Add oxygen
    borated_water_mat.add_nuclide('O16', 1.0, 'ao')  # One oxygen atom per molecule of water

    # Add boric acid (H3BO3) - Here we focus on the boron content, typically 1000 ppm by mass
    # This requires calculating the atomic fraction of boron in the solution
    # Boron's atomic mass is about 10.81, and water's molar mass is about 18.015

    # Assuming 1000 ppm of boron by mass:
    # Mass fraction of boron in borated water is 0.001 (1000 ppm = 0.001 by weight)
    # Boron concentration in moles per gram of borated water:
    boron_concentration = 0.001 / 10.81  # moles of boron per gram of borated water

    # Total moles in 1 gram of borated water:
    total_moles = boron_concentration + (1 - 0.001) / 18.015

    # Atomic fraction of boron:
    boron_af = boron_concentration / total_moles

    # Add boron
    borated_water_mat.add_nuclide('B10', 0.2 * boron_af, 'ao')  # 20% B-10 (common isotopic composition)
    borated_water_mat.add_nuclide('B11', 0.8 * boron_af, 'ao')  # 80% B-11

    cadmium_mat = openmc.Material(name='cadmium')

    # Set the density of cadmium, approximately 8.65 g/cm3
    cadmium_mat.set_density('g/cm3', 8.65)

    # Add cadmium (Cd). Natural cadmium consists of several isotopes, but Cd-113 is the major absorber.
    # You can add natural cadmium directly, or specify isotopes if needed for your simulation accuracy.
    cadmium_mat.add_element('Cd', 1.0, 'ao')  # Adding natural cadmium

    # Collect the materials together and export to XML
    #materials_file = openmc.Materials([heavy_water_mat, air_mat, borated_water_mat, cadmium_mat])
    materials_file = openmc.Materials([heavy_water_mat, borated_water_mat])
    #materials.export_to_xml()

    ###############################################################################
    # Create materials for the problem

    # Instantiate some Macroscopic Data
    #source_data = openmc.Macroscopic('source')
    #void_data   = openmc.Macroscopic('void')
    #shield_data = openmc.Macroscopic('shield')

    # Instantiate some Materials and register the appropriate Macroscopic objects
    #source_mat = openmc.Material(name='source')
    #source_mat.set_density('macro', 1.0)
    #source_mat.add_macroscopic(source_data)
    
    #void_mat = openmc.Material(name='void')
    #void_mat.set_density('macro', 1.0)
    #void_mat.add_macroscopic(void_data)
    
    #shield_mat = openmc.Material(name='shield')
    #shield_mat.set_density('macro', 1.0)
    #shield_mat.add_macroscopic(shield_data)

    # Instantiate a Materials collection and export to XML
    #materials_file = openmc.Materials([source_mat, void_mat, shield_mat])
    #materials_file.cross_sections = "mgxs.h5"

    ###############################################################################
    # Define problem geometry
    
    source_cell = openmc.Cell(fill=heavy_water_mat, name='infinite source region')
    void_cell = openmc.Cell(fill=heavy_water_mat, name='infinite void region')
    shield_cell = openmc.Cell(fill=borated_water_mat, name='infinite shield region')
    
    sub = openmc.Universe()
    sub.add_cells([source_cell])
    
    vub = openmc.Universe()
    vub.add_cells([void_cell])
    
    aub = openmc.Universe()
    aub.add_cells([shield_cell])
    
    # n controls the dimension of subdivision within each outer lattice element
    # E.g., n = 10 results in 1cm cubic FSRs
    n = 1
    delta = 10.0 / n
    ll = [-5.0, -5.0, -5.0]
    pitch = [delta, delta, delta]

    source_lattice = openmc.RectLattice()
    source_lattice.lower_left = ll
    source_lattice.pitch = pitch
    source_lattice.universes = fill_3d_list(n, sub)
    
    void_lattice = openmc.RectLattice()
    void_lattice.lower_left = ll
    void_lattice.pitch = pitch
    void_lattice.universes = fill_3d_list(n, vub)
    
    shield_lattice = openmc.RectLattice()
    shield_lattice.lower_left = ll
    shield_lattice.pitch = pitch
    shield_lattice.universes = fill_3d_list(n, aub)
    
    source_lattice_cell = openmc.Cell(fill=source_lattice, name='source lattice cell')
    su = openmc.Universe()
    su.add_cells([source_lattice_cell])
    
    void_lattice_cell = openmc.Cell(fill=void_lattice, name='void lattice cell')
    vu = openmc.Universe()
    vu.add_cells([void_lattice_cell])
    
    shield_lattice_cell = openmc.Cell(fill=shield_lattice, name='shield lattice cell')
    au = openmc.Universe()
    au.add_cells([shield_lattice_cell])

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
    
    lattice_cell = openmc.Cell(fill=lattice, name='dogleg lattice cell')

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
    
    # Create an initial uniform spatial source for ray integration
    lower_left = (0.0, 0.0, 0.0)
    upper_right = (x, y, z)
    uniform_dist = openmc.stats.Box(lower_left, upper_right, only_fissionable=False)
    rr_source = openmc.IndependentSource(space=uniform_dist)

    # Instantiate a Settings object, set all runtime parameters, and export to XML
    settings = openmc.Settings()
    #settings.energy_mode = "multi-group"
    settings.batches = 100
    settings.inactive = 0
    settings.particles = 100000
    #settings.particles = 2810
    #settings.particles = 4200
    #settings.particles = 10000
    settings.run_mode = 'fixed source'
    #settings.random_ray['distance_active'] = 400.0
    #settings.random_ray['distance_inactive'] = 100.0
    #settings.random_ray['ray_source'] = rr_source
    #settings.random_ray['volume_normalized_flux_tallies'] = False
    #settings.random_ray['source_shape'] = 'linear'
    
    # Create the neutron source in the bottom right of the moderator
    strengths = [1.0] # Good - fast group appears largest (besides most thermal)
    midpoints = [7.5e6]
    energy_distribution = openmc.stats.Discrete(x=midpoints,p=strengths)
    
    lower_left_src = [0.0, 0.0, 0.0]
    upper_right_src = [10.0, 10.0, 10.0]
    spatial_distribution = openmc.stats.Box(lower_left_src, upper_right_src, only_fissionable=False)

    #source = openmc.IndependentSource(energy=energy_distribution, domains=[source_mat], strength=2.0) # works
    source = openmc.IndependentSource(space=spatial_distribution, energy=energy_distribution,  strength=1.0) # works
    
    #source = openmc.IndependentSource(energy=energy_distribution, constraints={'domains':[source_cell]}, strength=1.0)
    #source = openmc.IndependentSource(energy=energy_distribution, domains=[source_mat], strength=1.0) # base source material
    #source = openmc.IndependentSource(energy=energy_distribution, domains=[sub], strength=1.0) # universe containing source cell
    #source = openmc.IndependentSource(energy=energy_distribution, domains=[source_cell], strength=1.0) # Material-filled cell
    #source = openmc.IndependentSource(energy=energy_distribution, domains=[source_lattice_cell], strength=1.0) # Higher level cell
    
    #source = openmc.IndependentSource(energy=energy_distribution, domains=[source_lattice], strength=1.0) # Lattice (BAD)


    # These are full problem sources
    #source = openmc.IndependentSource(energy=energy_distribution, domains=[full_domain], strength=1.0) # Highest level cell (making entire domain a source)
    #source = openmc.IndependentSource(energy=energy_distribution, domains=[root], strength=1.0) # Root Universe
    
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
    #tally.estimator = 'analog'

    # Instantiate a Tallies collection and export to XML
    #tallies = openmc.Tallies([tally])
    
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
    
    model = openmc.model.Model()
    model.geometry = geometry
    model.materials = materials_file
    model.settings = settings
    #model.xs_data = mg_cross_sections_file
    #model.tallies = tallies
    model.plots = plot_file

    return model

model = create_random_ray_model()
model.export_to_xml()
