import numpy as np
import matplotlib.pyplot as plt
import misc
import params
from md import MDSimulation
import time

energy_arr = []
beads_num_arr = [1, 5, 10, 20, 25, 30, 35, 40, 45, 50, 55, 60]

avg_per_bead = 1
bhw = 6  # beta*hbar*omega

t0 = time.time()

for beads_num in beads_num_arr:
    energy_p = []
    for i in range(avg_per_bead):
        md_sim = MDSimulation(bhw / (params.umap['hbar'] * params.qho_freq), dt=params.dt,
                              mass=params.mass, bead_num=int(beads_num), tsteps=params.steps,
                              estimator_type=params.estimator, threshold=params.threshold)
        md_sim.run()
        energy_p.append(md_sim.mean_energy / (params.umap['hbar'] * params.qho_freq))

    energy_arr.append(np.mean(np.array(energy_p)))

t1 = time.time()

fig = plt.figure()
plt.plot(beads_num_arr, np.array(energy_arr), '-o', color='black')
plt.xlabel(r"$P$", fontsize=15)
plt.ylabel(r"$\left\langle E\right\rangle / \hbar \omega_0$", fontsize=15)
plt.title(misc.print_title(bhw, params.dt, params.steps, params.gamma, params.estimator, type='e-vs-p'))
plt.show()

dir = 'e-vs-p'
timestamp = int(time.time())
fig.savefig("{}/{}.png".format(dir, timestamp), bbox_inches='tight')

misc.write_report(dir, beads_num_arr, energy_arr, timestamp,
                  {'Beta' : bhw,
                   'Threshold' : params.threshold,
                   'Averaging per bead' : avg_per_bead},
                  runtime=t1 - t0)
