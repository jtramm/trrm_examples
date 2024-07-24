import openmc
import openmc.mgxs as mgxs

summary = openmc.Summary('summary.h5')
geom = summary.geometry
mats = summary.materials

statepoint_filename = 'statepoint.100.h5'
sp = openmc.StatePoint(statepoint_filename)


# MGXS
groups = mgxs.EnergyGroups(mgxs.GROUP_STRUCTURES['CASMO-70'])
mgxs_lib = openmc.mgxs.Library(geom)
mgxs_lib.energy_groups = groups
mgxs_lib.correction = None
mgxs_lib.mgxs_types = ['total', 'absorption', 'nu-fission', 'fission',
                       'nu-scatter matrix', 'multiplicity matrix', 'chi']

# Specify a "cell" domain type for the cross section tally filters
mgxs_lib.domain_type = "material"

# Specify the cell domains over which to compute multi-group cross sections
mgxs_lib.domains = geom.get_all_materials().values()

# Do not compute cross sections on a nuclide-by-nuclide basis
mgxs_lib.by_nuclide = False


# Check the library - if no errors are raised, then the library is satisfactory.
mgxs_lib.check_library_for_openmc_mgxs()


# Construct all tallies needed for the multi-group cross section library
mgxs_lib.build_library()


# Create a "tallies.xml" file for the MGXS Library
tallies = openmc.Tallies()
mgxs_lib.add_to_tallies_file(tallies, merge=True)

# Export
tallies.export_to_xml()
