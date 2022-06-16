import numpy as np
import matplotlib.pyplot as plt
import params
from md import MDSimulation

beads_num = 10
save_freq = 1

dts = np.array([0.5, 1, 5, 10, 20])
steps = np.array([1e5]*len(dts))
sim_time = dts * steps

fig = plt.figure()

pos = []
momenta = []

for i in range(len(dts)):
    md_sim = MDSimulation(6 / (params.umap['hbar'] * params.qho_freq), dt=dts[i],
                          mass=params.mass, bead_num=beads_num, tsteps=steps[i],
                          estimator_type=params.estimator, save_freq=save_freq,
                          enable_thermostat=False, threshold=params.threshold, classical=True)
    if i == 0:
        pos = md_sim.pos_arr
        momenta = md_sim.momenta
    else:
        md_sim.force_initial_conditions(pos, momenta)

    md_sim.run()

    time_interval = np.arange(params.threshold * sim_time[i], sim_time[i], dts[i])
    energies = md_sim.classical_energy

    skip_points = 500

    mean = np.mean(energies)
    e_out = np.array(energies[::skip_points])/mean - 1

    plt.plot(time_interval[::skip_points], e_out, marker='.',
             label=r"$\Delta t$={}".format(dts[i]), markersize=5)

    plt.axhline(y=0, color='r', linestyle='--')
    plt.ticklabel_format(useMathText=True, scilimits=(0,0))

    print("ΔE/E: {}%".format(md_sim.get_classical_de()))

plt.xlabel(r"$t$", fontsize=15)
plt.ylabel(r"$\Delta E / \bar{E}$", fontsize=15)
plt.legend()
plt.show()
fig.savefig("out/e-vs-t.png", bbox_inches='tight')