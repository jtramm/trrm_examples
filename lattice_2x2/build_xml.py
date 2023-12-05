from math import log10

import numpy as np

import openmc
import openmc.mgxs

###############################################################################
# Create multigroup data

# Instantiate the energy group data
groups = openmc.mgxs.EnergyGroups(group_edges=[
    1e-5, 0.0635, 10.0, 1.0e2, 1.0e3, 0.5e6, 1.0e6, 20.0e6])

# Instantiate the 7-group (C5G7) cross section data
uo2_xsdata = openmc.XSdata('UO2', groups)
uo2_xsdata.order = 0
uo2_xsdata.set_total(
    [0.1779492, 0.3298048, 0.4803882, 0.5543674, 0.3118013, 0.3951678,
     0.5644058])
uo2_xsdata.set_absorption([8.0248E-03, 3.7174E-03, 2.6769E-02, 9.6236E-02,
                           3.0020E-02, 1.1126E-01, 2.8278E-01])
scatter_matrix = np.array(
    [[[0.1275370, 0.0423780, 0.0000094, 0.0000000, 0.0000000, 0.0000000, 0.0000000],
      [0.0000000, 0.3244560, 0.0016314, 0.0000000, 0.0000000, 0.0000000, 0.0000000],
      [0.0000000, 0.0000000, 0.4509400, 0.0026792, 0.0000000, 0.0000000, 0.0000000],
      [0.0000000, 0.0000000, 0.0000000, 0.4525650, 0.0055664, 0.0000000, 0.0000000],
      [0.0000000, 0.0000000, 0.0000000, 0.0001253, 0.2714010, 0.0102550, 0.0000000],
      [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0012968, 0.2658020, 0.0168090],
      [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0085458, 0.2730800]]])
scatter_matrix = np.rollaxis(scatter_matrix, 0, 3)
uo2_xsdata.set_scatter_matrix(scatter_matrix)
uo2_xsdata.set_fission([7.21206E-03, 8.19301E-04, 6.45320E-03,
                        1.85648E-02, 1.78084E-02, 8.30348E-02,
                        2.16004E-01])
uo2_xsdata.set_nu_fission([2.005998E-02, 2.027303E-03, 1.570599E-02,
                           4.518301E-02, 4.334208E-02, 2.020901E-01,
                           5.257105E-01])
uo2_xsdata.set_chi([5.8791E-01, 4.1176E-01, 3.3906E-04, 1.1761E-07, 0.0000E+00,
                    0.0000E+00, 0.0000E+00])

h2o_xsdata = openmc.XSdata('LWTR', groups)
h2o_xsdata.order = 0
h2o_xsdata.set_total([0.15920605, 0.412969593, 0.59030986, 0.58435,
                      0.718, 1.2544497, 2.650379])
h2o_xsdata.set_absorption([6.0105E-04, 1.5793E-05, 3.3716E-04,
                           1.9406E-03, 5.7416E-03, 1.5001E-02,
                           3.7239E-02])
scatter_matrix = np.array(
    [[[0.0444777, 0.1134000, 0.0007235, 0.0000037, 0.0000001, 0.0000000, 0.0000000],
      [0.0000000, 0.2823340, 0.1299400, 0.0006234, 0.0000480, 0.0000074, 0.0000010],
      [0.0000000, 0.0000000, 0.3452560, 0.2245700, 0.0169990, 0.0026443, 0.0005034],
      [0.0000000, 0.0000000, 0.0000000, 0.0910284, 0.4155100, 0.0637320, 0.0121390],
      [0.0000000, 0.0000000, 0.0000000, 0.0000714, 0.1391380, 0.5118200, 0.0612290],
      [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0022157, 0.6999130, 0.5373200],
      [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.1324400, 2.4807000]]])
scatter_matrix = np.rollaxis(scatter_matrix, 0, 3)
h2o_xsdata.set_scatter_matrix(scatter_matrix)

mg_cross_sections_file = openmc.MGXSLibrary(groups)
mg_cross_sections_file.add_xsdatas([uo2_xsdata, h2o_xsdata])
mg_cross_sections_file.export_to_hdf5()

###############################################################################
# Create materials for the problem

# Instantiate some Macroscopic Data
uo2_data = openmc.Macroscopic('UO2')
h2o_data = openmc.Macroscopic('LWTR')

# Instantiate some Materials and register the appropriate Macroscopic objects
uo2 = openmc.Material(name='UO2 fuel')
uo2.set_density('macro', 1.0)
uo2.add_macroscopic(uo2_data)

water = openmc.Material(name='Water')
water.set_density('macro', 1.0)
water.add_macroscopic(h2o_data)

# Instantiate a Materials collection and export to XML
materials_file = openmc.Materials([uo2, water])
materials_file.cross_sections = "mgxs.h5"
materials_file.export_to_xml()

###############################################################################
# Define problem geometry

########################################
# Define an unbounded pincell universe

pitch = 1.26

# Create a surface for the fuel outer radius
fuel_or =      openmc.ZCylinder(r=0.54, name='Fuel OR')
inner_ring_a = openmc.ZCylinder(r=0.33, name='inner ring a')
inner_ring_b = openmc.ZCylinder(r=0.45, name='inner ring b')
outer_ring_a = openmc.ZCylinder(r=0.60, name='outer ring a')
outer_ring_b = openmc.ZCylinder(r=0.69, name='outer ring b')

# Instantiate Cells
fuel_a = openmc.Cell(fill=uo2, region=-inner_ring_a, name='fuel inner a')
fuel_b = openmc.Cell(fill=uo2, region=+inner_ring_a & -inner_ring_b, name='fuel inner b')
fuel_c = openmc.Cell(fill=uo2, region=+inner_ring_b & -fuel_or, name='fuel inner c')
moderator_a = openmc.Cell(fill=water, region=+fuel_or & -outer_ring_a, name='moderator inner a')
moderator_b = openmc.Cell(fill=water, region=+outer_ring_a & -outer_ring_b, name='moderator outer b')
moderator_c = openmc.Cell(fill=water, region=+outer_ring_b, name='moderator outer c')

# Create pincell universe
pincell_base = openmc.Universe()

# Register Cells with Universe
pincell_base.add_cells([fuel_a, fuel_b, fuel_c, moderator_a, moderator_b, moderator_c])

# Create planes for azimuthal sectors
azimuthal_planes = []
for i in range(8):
    angle = 2 * i * openmc.pi / 8
    normal_vector = (-openmc.sin(angle), openmc.cos(angle), 0)
    azimuthal_planes.append(openmc.Plane(a=normal_vector[0], b=normal_vector[1], c=normal_vector[2], d=0))

# Create a cell for each azimuthal sector
azimuthal_cells = []
for i in range(8):
    azimuthal_cell = openmc.Cell(name=f'azimuthal_cell_{i}')
    azimuthal_cell.fill = pincell_base
    azimuthal_cell.region = +azimuthal_planes[i] & -azimuthal_planes[(i+1) % 8]
    azimuthal_cells.append(azimuthal_cell)

# Create a geometry with the azimuthal universes
pincell = openmc.Universe(cells=azimuthal_cells)

########################################
# Define a moderator lattice universe

moderator_infinite = openmc.Cell(fill=water, name='moderator infinite')
mu = openmc.Universe()
mu.add_cells([moderator_infinite])

lattice = openmc.RectLattice()
lattice.lower_left = [-pitch/2.0, -pitch/2.0]
lattice.pitch = [pitch/10.0, pitch/10.0]
lattice.universes = [
        [mu, mu, mu, mu, mu, mu, mu, mu, mu, mu],
        [mu, mu, mu, mu, mu, mu, mu, mu, mu, mu],
        [mu, mu, mu, mu, mu, mu, mu, mu, mu, mu],
        [mu, mu, mu, mu, mu, mu, mu, mu, mu, mu],
        [mu, mu, mu, mu, mu, mu, mu, mu, mu, mu],
        [mu, mu, mu, mu, mu, mu, mu, mu, mu, mu],
        [mu, mu, mu, mu, mu, mu, mu, mu, mu, mu],
        [mu, mu, mu, mu, mu, mu, mu, mu, mu, mu],
        [mu, mu, mu, mu, mu, mu, mu, mu, mu, mu],
        [mu, mu, mu, mu, mu, mu, mu, mu, mu, mu]
                     ]

mod_lattice_cell = openmc.Cell(fill=lattice)

mod_lattice_uni = openmc.Universe()

mod_lattice_uni.add_cells([mod_lattice_cell])

########################################
# Define 2x2 outer lattice
lattice2x2 = openmc.RectLattice()
lattice2x2.lower_left = [-pitch, -pitch]
lattice2x2.pitch = [pitch, pitch]
lattice2x2.universes = [
        [pincell, pincell],
        [pincell, mod_lattice_uni]
        ]

########################################
# Define cell containing lattice and other stuff
box = openmc.rectangular_prism(pitch*2, pitch*2, boundary_type='reflective')

assembly = openmc.Cell(fill=lattice2x2, region=box, name='assembly')

root = openmc.Universe(name='root universe')
root.add_cell(assembly)

# Create a geometry with the two cells and export to XML
geometry = openmc.Geometry(root)
geometry.export_to_xml()

###############################################################################
# Define problem settings

# Instantiate a Settings object, set all runtime parameters, and export to XML
settings = openmc.Settings()
settings.energy_mode = "multi-group"
settings.batches = 1000
settings.inactive = 100
settings.particles = 1000
settings.random_ray_distance_active = 1000.0
settings.random_ray_distance_inactive = 100.0
settings.solver_type = 'random ray'

# Create an initial uniform spatial source distribution over fissionable zones
lower_left = (-pitch, -pitch, -1)
upper_right = (pitch, pitch, 1)
uniform_dist = openmc.stats.Box(lower_left, upper_right, only_fissionable=False)
settings.source = openmc.IndependentSource(space=uniform_dist)
settings.export_to_xml()

###############################################################################
# Define tallies

# Create a mesh that will be used for tallying
mesh = openmc.RegularMesh()
mesh.dimension = (100, 100)
mesh.lower_left = (-pitch/2, -pitch/2)
mesh.upper_right = (pitch/2, pitch/2)

# Create a mesh filter that can be used in a tally
mesh_filter = openmc.MeshFilter(mesh)

# Now use the mesh filter in a tally and indicate what scores are desired
mesh_tally = openmc.Tally(name="Mesh tally")
mesh_tally.filters = [mesh_filter]
mesh_tally.scores = ['flux', 'fission', 'nu-fission']
mesh_tally.estimator = 'analog'

# Let's also create a tally to get the flux energy spectrum. We start by
# creating an energy filter
e_min, e_max = 1e-5, 20.0e6
groups = 500
energies = np.logspace(log10(e_min), log10(e_max), groups + 1)
energy_filter = openmc.EnergyFilter(energies)

spectrum_tally = openmc.Tally(name="Flux spectrum")
spectrum_tally.filters = [energy_filter]
spectrum_tally.scores = ['flux']
spectrum_tally.estimator = 'analog'

# Instantiate a Tallies collection and export to XML
tallies = openmc.Tallies([mesh_tally, spectrum_tally])
tallies.export_to_xml()

###############################################################################
#                   Exporting to OpenMC plots.xml file
###############################################################################

plot = openmc.Plot()
plot.origin = [0, 0, 0]
plot.width = [2*pitch, 2*pitch, 2*pitch]
plot.pixels = [400, 400, 1]
plot.color_by = 'material'
plot.type = 'voxel'

# Instantiate a Plots collection and export to XML
plot_file = openmc.Plots([plot])
plot_file.export_to_xml()
