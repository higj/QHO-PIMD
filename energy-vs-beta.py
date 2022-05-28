import numpy as np
import matplotlib.pyplot as plt
import misc
import params
from md import MDSimulation

energy_arr = []
betas = np.linspace(0.1, 6, 20)  # Range of betas for the simulation
betas_cont = np.linspace(0.1, 6, 100)  # Range of betas for the analytic function

beads_num = 20
dt = 0.00001
steps = 20000

for beta in betas:
    md_sim = MDSimulation(beta / (params.hbar * params.qho_freq), dt=dt,
                          mass=params.mass, bead_num=beads_num, tsteps=steps,
                          estimator_type='centroid_virial')
    md_sim.run()

    # Add the (normalised) energy to the array
    energy_arr.append(md_sim.mean_energy / (params.hbar * params.qho_freq))

fig = plt.figure()
plt.plot(betas, np.array(energy_arr), '-o', color='black', label='PIMD')
plt.plot(betas_cont, misc.energy_vs_beta(betas_cont), color='r', label='Theory')
plt.xlabel(r"$\beta \hbar \omega_0$", fontsize=15)
plt.ylabel(r"$\left\langle E\right\rangle / \hbar \omega_0$", fontsize=15)
plt.title(r"PIMD for QHO ($P={0}$)".format(beads_num))
plt.legend()
plt.show()