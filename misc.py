import numpy as np
from scipy.special import erf
from scipy.interpolate import interp1d


# Generate a 1D Maxwellian CDF
def mb_cdf(v, beta, mass):
    a = np.sqrt(1 / (beta * mass))
    return 0.5 * (1 + erf(v / (np.sqrt(2) * a)))


# Invert the CDF for sampling purposes
def inv_cdf(beta, mass):
    v = np.arange(0, 3000, 0.1)  # Velocities in the range 0-3000 m/s
    cdf = mb_cdf(v, beta, mass)

    return interp1d(cdf, v, fill_value="extrapolate")


# Analytic expression for average energy in a 1D QHO vs thermodynamic beta
def energy_vs_beta(beta):
    return 0.5 + 1 / (np.exp(beta) - 1)
