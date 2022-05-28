import numpy as np
import matplotlib.pyplot as plt
import params
from md import MDSimulation

energy_arr = []
beads_num_arr = np.linspace(1, 60, 10)

for beads_num in beads_num_arr:
    md_sim = MDSimulation(5 / (params.hbar * params.qho_freq), dt=0.00001,
                          mass=params.mass, bead_num=int(beads_num), tsteps=20000,
                          estimator_type='default')
    md_sim.run()
    energy_arr.append(md_sim.mean_energy / (params.hbar * params.qho_freq))

fig = plt.figure()
plt.plot(beads_num_arr, np.array(energy_arr), '-o', color='black')
plt.xlabel(r"$P$", fontsize=15)
plt.ylabel(r"$\left\langle E\right\rangle / \hbar \omega_0$", fontsize=15)
plt.show()