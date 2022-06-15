import units

# Choose a unit system
unit_type = 'atomic'

umap = units.UnitSystem(unit_type).factors

# System parameters
qho_freq = 0.05*umap['eV']
mass = 39.9480*umap['amu']  # Argon mass in atomic units

# Simulation parameters
dt = 1  # In atomic units of time
steps = 2e6
threshold = 0.1
save_freq = 1000
gamma = 1/(100*dt)  # Friction term of the Langevin thermostat
estimator = 'default'  # 'default' (primitive), 'centroid_virial', 'virial'
