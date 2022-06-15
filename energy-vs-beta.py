import numpy as np
import matplotlib.pyplot as plt
import misc
import params
from md import MDSimulation
import time

energy_arr = []
bhw_arr = np.linspace(0.1, 6, 10)
bhw_cont_arr = np.linspace(0.1, 6, 100)

beads_num = 15
avg_per_beta = 1

t0 = time.time()

for bhw in bhw_arr:
    bead_energy = []
    for i in range(avg_per_beta):
        beta = bhw / (params.umap['hbar'] * params.qho_freq)
        md_sim = MDSimulation(beta, dt=params.dt,
                              mass=params.mass, bead_num=beads_num, tsteps=params.steps,
                              estimator_type=params.estimator, save_freq=params.save_freq,
                              threshold=params.threshold)
        md_sim.run()

        bead_energy.append(md_sim.mean_energy / (params.umap['hbar'] * params.qho_freq))

    energy_arr.append(np.mean(np.array(bead_energy)))

t1 = time.time()

fig = plt.figure()
plt.plot(bhw_arr, np.array(energy_arr), '-o', color='black', label='PIMD')
plt.plot(bhw_cont_arr, misc.energy_vs_beta(bhw_cont_arr), color='r', label='Theory')
plt.xlabel(r"$\beta \hbar \omega_0$", fontsize=15)
plt.ylabel(r"$\left\langle E\right\rangle / \hbar \omega_0$", fontsize=15)
plt.title(misc.print_title(beads_num, params.dt, params.steps, params.gamma, params.estimator))
plt.legend()
plt.show()

dir = 'out/e-vs-beta'
timestamp = int(time.time())
fig.savefig("{}/{}.png".format(dir, timestamp), bbox_inches='tight')

misc.write_report(dir, bhw_arr, energy_arr, timestamp,
                  {'Number of beads (P)' : beads_num,
                   'Threshold' : params.threshold,
                   'Averaging per beta' : avg_per_beta},
                  runtime=t1 - t0)
