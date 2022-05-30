# PIMD Simulation for 1D QHO

Implementation of PIMD for a one-dimensional argon particle in a QHO potential with $\hbar \omega_0 = 50 \~ \mathrm{meV}$.

## Configuring the script

* `params.py` contains system parameters such as the frequency of the oscillator, its mass and the friction term of the Langevin thermostat

## Energy vs beta

For the default (primitive) energy estimator:

<p align="center"><img src="https://i.imgur.com/vJcDdrw.png" alt="Primitive energy estimator (8 beads)" width="450" /></p>

For the centroid-virial kinetic energy estimator (two runs):

<p align="center"><img src="https://i.imgur.com/aOTR4xK.png" alt="CVKE (20 beads)" width="450" /><img src="https://i.imgur.com/BVSy8k3.png" alt="CVKE (20 beads)" width="450" /></p>

## Energy vs number of beads

For the default (primitive) energy estimator ($\beta\hbar\omega_0 = 5$):

<p align="center"><img src="https://i.imgur.com/XNHC0VR.png" alt="Primitive energy estimator (8 beads)" width="450" /></p>

For the centroid-virial kinetic energy estimator ($\beta\hbar\omega_0 = 5$):

<p align="center"><img src="https://i.imgur.com/vwUNzgX.png" width="450" /></p>

For the centroid-virial kinetic energy estimator with $\beta\hbar\omega_0 = 6$ and averaging over $10$ simulations per beach bead:

<p align="center"><img src="https://i.imgur.com/ElCvpF4.png" width="450" /></p>

## Mean energy computation

The primitive (or *thermodynamic*) energy estimator is calculated using Eqns. (2.53) and (3.3) from [CMMR]:

$$E^{\mathrm{TD}}=\frac{P}{2\beta}\left[1-\sum_{i=1}^{P}\frac{m}{\beta\hbar^{2}}\left(q_{i}-q_{i+1}\right)^{2}\right]+\frac{1}{P}\sum_{i=1}^{P}V\left(q_{i}\right)$$

where

$$V\left(x\right)=\frac{1}{2}m\omega_{0}^{2}x^{2}$$

is the potential energy associated with the QHO. 

The centroid-virial energy estimator is calculated using Eqn. (3.15) from [CMMR]:

$$E^{\mathrm{CV}} = \frac{1}{2\beta} + \frac{1}{2 P}\sum_{i=1}^{P} \left( q_i - \bar{q} \right) \frac{\partial V}{\partial q_i} + \frac{1}{P}\sum_{i=1}^{P}V\left(q_{i}\right),\quad \bar{q} = \frac{1}{P} \sum_{i=1}^{P} q_i $$

In all cases the expectation value is evaluated using Eqn. (2.37) from [CMMR]

$$\left< \hat{E} \right> = \frac{1}{T} \int_{0}^{T} dt \~ E_P(t)$$

where $T$ is the total simulation time and $E_P$ is a PIMD estimator for the total energy. This expression is equivalent to the ensemble average if the system is ergodic. In practice, the continuous integral is replaced by the discrete sum over all the time steps.

## References

* [CMMR] "*Path Integral Methods in Atomistic Modelling - An Introduction*", Ceriotti M., Manolopoulos D. E., Markland T. E., Rossi M., 2021
