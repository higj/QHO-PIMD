# PIMD Simulation for 1D QHO

Implementation of PIMD for a one-dimensional argon atom in a QHO potential with $\hbar \omega_0 = 50 \~ \mathrm{meV}$.

## Configuring the script

* `params.py` contains system parameters such as the frequency of the oscillator, its mass and the friction term of the Langevin thermostat

* `units.py` contains a class which handles units. By default, Hartree atomic units are used

## Calibrating the timestep

The optimal value of $\Delta t$ for a given $\beta$ is defined to be the largest $\Delta t$ for which the
deviations of the classical energies of the polymer do not exceed $0.1%$. In other words, the criterion is

$$
\left | \frac{\Delta E}{E_0} \cdot 100% \right | \leq 0.1 
$$

Notice that for different temperatures (betas) the optimal $\Delta t$ will be different. Therefore,
if one insists on using the same value of $\Delta t$ for different temperatures, one must ensure 
convergence in the desired temperature range.

The calibration procedure described above is performed in `cl-e-vs-t.py`. 
There we take $P=10$ and $\beta \hbar \omega_0 = 6$ and perform a classical simulation for different 
time steps:

<p align="center"><img src="https://i.imgur.com/rJOpSmz.png" width="450" /></p>

In this particular case the times are measured in atomic units ($\sim 0.0242 \, \mathrm{fs}$).

## Energy vs β

After choosing the optimal time step (`dt=1` in this case) we set the Langevin friction term to 
$\gamma = \frac{1}{100 \Delta t}$, calculate the mean *quantum* energy of the harmonic oscillator
as a function of the thermodynamic beta, and then compare it to the theoretically predicted curve:

<p align="center"><img src="https://i.imgur.com/m6e4LNn.png" alt="Primitive energy estimator (15 beads)" width="450" /></p>

For the virial energy estimator:

<p align="center"><img src="https://i.imgur.com/Jg9sWnx.png" alt="Virial (15 beads)" width="450" />

## Energy vs number of beads

The optimal time step is higher in the low temperature regime (high $\beta$). 
Therefore, the dependence of the average energy for $\beta \hbar \omega_0 = 6$ can
be calculated with `dt=15`. 

For the default (primitive) energy estimator ($\beta\hbar\omega_0 = 6$):

<p align="center"><img src="https://i.imgur.com/oy2wcBJ.png" alt="Primitive energy estimator" width="450" /></p>

For the virial kinetic energy estimator ($\beta\hbar\omega_0 = 6$):

<p align="center"><img src="" width="450" /></p>


## Mean energy computation

The **primitive** (or *thermodynamic*) energy estimator is calculated using Eqns. (2.53) and (3.3) from [CMMR]:

$$E^{\mathrm{TD}}=\frac{P}{2\beta}\left[1-\sum_{i=1}^{P}\frac{m}{\beta\hbar^{2}}\left(q_{i}-q_{i+1}\right)^{2}\right]+\frac{1}{P}\sum_{i=1}^{P}V\left(q_{i}\right)$$

where

$$V\left(x\right)=\frac{1}{2}m\omega_{0}^{2}x^{2}$$

is the potential energy associated with the QHO. 

The **virial** energy estimator is calculated using Eqn. (12.6.34) in [T]:

$$E^{\mathrm{vir}} = \frac{1}{P} \sum_{k=1}^P \left [ \frac{1}{2} q_k \frac{\partial V}{\partial q_k} + V(q_k) \right ]$$

The **centroid-virial** energy estimator is calculated using Eqn. (3.15) from [CMMR]:

$$E^{\mathrm{CV}} = \frac{1}{2\beta} + \frac{1}{2 P}\sum_{i=1}^{P} \left( q_i - \bar{q} \right) \frac{\partial V}{\partial q_i} + \frac{1}{P}\sum_{i=1}^{P}V\left(q_{i}\right),\quad \bar{q} = \frac{1}{P} \sum_{i=1}^{P} q_i $$

In all cases the expectation value is evaluated using Eqn. (2.37) from [CMMR]

$$\left< \hat{E} \right> = \frac{1}{T} \int_{0}^{T} dt \~ E_P(t) \approx \frac{n_{\mathrm{freq}}}{(1-t) n_{\mathrm{steps}}} \sum_{\mathrm{iterations}} E_P (t)$$

where $T$ is the total simulation time, $E_P$ is a PIMD estimator for the total energy, $n_{\mathrm{steps}}$ is the number of steps, $n_{\mathrm{freq}}$ is the recording step and $t$ is the threshold (between 0 and 1). This expression is equivalent to the ensemble average if the system is ergodic. In practice, the continuous integral is replaced by the discrete sum over all the time steps.

## References

* [CMMR] "*Path Integral Methods in Atomistic Modelling - An Introduction*", Ceriotti M., Manolopoulos D. E., Markland T. E., Rossi M., 2021
* [T] "*Statistical Mechanics: Theory and Molecular Simulation*", Tuckerman M.