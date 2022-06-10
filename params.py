# Physical constants
hbar = 0.6582  # meV*ps
kB = 0.08617  # meV/K

# System parameters
qho_freq = 50 / hbar  # QHO frequency (THz)
mass = 0.04778 * hbar / qho_freq  # Argon mass expressed through other constants

gamma = 100
dt = 1e-3
steps = 1e5
save_freq = 1000
estimator = 'default'