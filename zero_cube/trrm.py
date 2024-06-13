import numpy as np

import openmc

def fill_cube(N, n_1, n_2, fill_1, fill_2, fill_3):
    # Create an empty NxNxN cube as a list of lists of lists
    cube = [[[0 for _ in range(N)] for _ in range(N)] for _ in range(N)]

    # Loop through each element in the cube to assign values based on the indices
    for i in range(N):
        for j in range(N):
            for k in range(N):
                if i < n_1 or j < n_1 or k < n_1:
                    cube[i][j][k] = fill_1
                elif (i >= n_1 and i < n_2) or (j >= n_1 and j < n_2) or (k >= n_1 and k < n_2):
                    cube[i][j][k] = fill_2
                else:
                    cube[i][j][k] = fill_3

    return cube

###############################################################################
# Create multigroup data

# Instantiate the energy group data
ebins = [1e-5, 20.0e6]
groups = openmc.mgxs.EnergyGroups(group_edges=ebins)

void_sigma_a = 4.0e-6
void_sigma_s = 3.0e-4
void_mat_data = openmc.XSdata('void', groups)
void_mat_data.order = 0
void_mat_data.set_total([void_sigma_a + void_sigma_s])
void_mat_data.set_absorption([void_sigma_a])
void_mat_data.set_scatter_matrix(np.rollaxis(np.array([[[void_sigma_s]]]),0,3))

absorber_sigma_a = 0.75
absorber_sigma_s = 0.25
absorber_mat_data = openmc.XSdata('absorber', groups)
absorber_mat_data.order = 0
absorber_mat_data.set_total([absorber_sigma_a + absorber_sigma_s])
absorber_mat_data.set_absorption([absorber_sigma_a])
absorber_mat_data.set_scatter_matrix(np.rollaxis(np.array([[[absorber_sigma_s]]]),0,3))

multiplier = 1.0
source_sigma_a = void_sigma_a * multiplier
source_sigma_s = void_sigma_s * multiplier
source_mat_data = openmc.XSdata('source', groups)
source_mat_data.order = 0
source_mat_data.set_total([source_sigma_a + source_sigma_s])
source_mat_data.set_absorption([source_sigma_a])
source_mat_data.set_scatter_matrix(np.rollaxis(np.array([[[source_sigma_s]]]),0,3))

mg_cross_sections_file = openmc.MGXSLibrary(groups)
mg_cross_sections_file.add_xsdatas([source_mat_data, void_mat_data, absorber_mat_data])
mg_cross_sections_file.export_to_hdf5()

###############################################################################
# Create materials for the problem

# Instantiate some Macroscopic Data
source_data = openmc.Macroscopic('source')
void_data   = openmc.Macroscopic('void')
absorber_data = openmc.Macroscopic('absorber')

# Instantiate some Materials and register the appropriate Macroscopic objects
source_mat = openmc.Material(name='source')
source_mat.set_density('macro', 1.0)
source_mat.add_macroscopic(source_data)

void_mat = openmc.Material(name='void')
void_mat.set_density('macro', 1.0)
void_mat.add_macroscopic(void_data)

absorber_mat = openmc.Material(name='absorber')
absorber_mat.set_density('macro', 1.0)
absorber_mat.add_macroscopic(absorber_data)

# Instantiate a Materials collection and export to XML
materials_file = openmc.Materials([source_mat, void_mat, absorber_mat])
materials_file.cross_sections = "mgxs.h5"

###############################################################################
# Define problem geometry

# Basically, I want to define the three cell and universe types, and 
# then define ONE large lattice that has fills specified by the geometry.

source_cell = openmc.Cell(fill=source_mat, name='infinite source region')
void_cell = openmc.Cell(fill=void_mat, name='infinite void region')
absorber_cell = openmc.Cell(fill=absorber_mat, name='infinite absorber region')

source_universe = openmc.Universe()
source_universe.add_cells([source_cell])

void_universe = openmc.Universe()
void_universe.add_cells([void_cell])

absorber_universe = openmc.Universe()
absorber_universe.add_cells([absorber_cell])

# Here we define the 3D lattice pattern that will be used

source_width = 5.0
void_width = 25.0
absorber_width = 30.0

n_base = 6

# This variable can be increased above 1 to refine the FSR mesh resolution further
refinement_level = 1

n = n_base * refinement_level
pitch = absorber_width / n
print(f"pitch = {pitch}")

pattern = fill_cube(n, 1*refinement_level, 5*refinement_level, source_universe, void_universe, absorber_universe)

lattice = openmc.RectLattice()
lattice.lower_left = [0.0, 0.0, 0.0]
lattice.pitch = [pitch, pitch, pitch]
lattice.universes = pattern

lattice_cell = openmc.Cell(fill=lattice)

lattice_uni = openmc.Universe()
lattice_uni.add_cells([lattice_cell])

x_low  = openmc.XPlane(x0=0.0,boundary_type='reflective') 
x_high = openmc.XPlane(x0=absorber_width,boundary_type='vacuum') 

y_low  = openmc.YPlane(y0=0.0,boundary_type='reflective') 
y_high = openmc.YPlane(y0=absorber_width,boundary_type='vacuum') 

z_low  = openmc.ZPlane(z0=0.0,boundary_type='reflective') 
z_high = openmc.ZPlane(z0=absorber_width,boundary_type='vacuum') 

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
settings.batches = 200
settings.inactive = 100
settings.particles = 1000
settings.run_mode = 'fixed source'

# Create an initial uniform spatial source for ray integration
lower_left_ray = [0.0, 0.0, 0.0]
upper_right_ray = [absorber_width, absorber_width, absorber_width]
uniform_dist_ray = openmc.stats.Box(lower_left_ray, upper_right_ray, only_fissionable=False)
rr_source = openmc.IndependentSource(space=uniform_dist_ray)

settings.random_ray['distance_active'] = 500.0
settings.random_ray['distance_inactive'] = 100.0
settings.random_ray['ray_source'] = rr_source
settings.random_ray['volume_estimator'] = 'naive'

# Create the neutron source in the bottom right of the moderator
strengths = [1.0] # Good - fast group appears largest (besides most thermal)
midpoints = [100.0]
energy_distribution = openmc.stats.Discrete(x=midpoints,p=strengths)

lower_left_src = [0.0, 0.0, 0.0]
upper_right_src = [source_width, source_width, source_width]
spatial_distribution = openmc.stats.Box(lower_left_src, upper_right_src, only_fissionable=False)

#source = openmc.IndependentSource(energy=energy_distribution, domains=[source_mat], strength=2.0) # works
source = openmc.IndependentSource(energy=energy_distribution, constraints={'domains':[source_mat]}, strength=1.0) # works

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

absorber_filter = openmc.MaterialFilter(absorber_mat)
absorber_tally = openmc.Tally(name="Absorber Tally")
absorber_tally.filters = [absorber_filter]
absorber_tally.scores = ['flux']
absorber_tally.estimator = estimator

void_filter = openmc.MaterialFilter(void_mat)
void_tally = openmc.Tally(name="Void Tally")
void_tally.filters = [void_filter]
void_tally.scores = ['flux']
void_tally.estimator = estimator

source_filter = openmc.MaterialFilter(source_mat)
source_tally = openmc.Tally(name="Source Tally")
source_tally.filters = [source_filter]
source_tally.scores = ['flux']
source_tally.estimator = estimator


# Instantiate a Tallies collection and export to XML
tallies = openmc.Tallies([source_tally, void_tally, absorber_tally])

###############################################################################
#                   Exporting to OpenMC plots.xml file
###############################################################################

plot = openmc.Plot()
plot.origin = [absorber_width/2.0, absorber_width/2.0, absorber_width/2.0]
plot.width = [absorber_width, absorber_width, absorber_width]
plot.pixels = [n, n, n]
plot.type = 'voxel'

# Instantiate a Plots collection and export to XML
plot_file = openmc.Plots([plot])
#plot_file.export_to_xml()

model = openmc.model.Model()
model.geometry = geometry
model.materials = materials_file
model.settings = settings
model.xs_data = mg_cross_sections_file
model.tallies = tallies
model.plots = plot_file

model.export_to_model_xml()
