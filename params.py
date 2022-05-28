# Physical constants
hbar = 0.6582  # meV*ps
kB = 0.08617  # meV/K

# System parameters
qho_freq = 50 / hbar  # QHO frequency (THz)
gamma = 0.3  # Friction term for the Langevin thermostat

mass = 0.04778 * hbar / qho_freq  # Argon mass expressed through other constants
