import numpy as np
import matplotlib.pyplot as plt
import misc
import params
from md import MDSimulation

energy_arr = []
betas = np.linspace(0.1, 6, 30)
betas_cont = np.linspace(0.1, 6, 100)

beads_num = 20

for beta in betas:
    md_sim = MDSimulation(beta / (params.hbar * params.qho_freq), dt=params.dt,
                          mass=params.mass, bead_num=beads_num, tsteps=params.steps,
                          estimator_type=params.estimator, save_freq=params.save_freq)
    md_sim.run()
    energy_arr.append(md_sim.mean_energy / (params.hbar * params.qho_freq))

fig = plt.figure()
plt.plot(betas, np.array(energy_arr), '-o', color='black', label='PIMD')
plt.plot(betas_cont, misc.energy_vs_beta(betas_cont), color='r', label='Theory')
plt.xlabel(r"$\beta \hbar \omega_0$", fontsize=15)
plt.ylabel(r"$\left\langle E\right\rangle / \hbar \omega_0$", fontsize=15)
plt.title(misc.print_title(beads_num, params.dt, params.steps, params.gamma, params.estimator))
plt.legend()
plt.show()
