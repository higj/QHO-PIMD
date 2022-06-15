import numpy as np
import params
import time
import os
import shutil


# Analytic expression for average energy in a 1D QHO vs thermodynamic beta
def energy_vs_beta(beta):
    return 0.5 + 1 / (np.exp(beta) - 1)


def print_title(specific, dt, steps, gamma, est, type='e-vs-beta'):
    est_dic = {'default': 'Primitive', 'centroid_virial': 'CVKE', 'virial' : 'Virial', 'simple_virial_pe': '2V', 'simple_virial_ke': '2K'}

    specific_label = ''

    if type == 'e-vs-beta':
        specific_label = 'P={}'.format(specific)
    elif type == 'e-vs-p':
        specific_label = r'\beta \hbar \omega_0={}'.format(specific)

    return r"${0}, \Delta t={1:.1e}, n={2:.1e}, \gamma = {3:.1e}$, Estimator: {4}".format(specific_label, dt, steps, gamma, est_dic[est])


def write_report(dir, x_arr, y_arr, timestamp = -1, extra = {}, runtime = -1):
    report = [
        '--- System parameters ---',
        'mass = {}'.format(params.mass),
        'oscillator frequency = {}'.format(params.qho_freq),
        'dt = {}'.format(params.dt),
        'gamma = {}'.format(params.gamma),
        'steps = {}'.format(params.steps),
        'save_freq = {}'.format(params.save_freq),
        'estimator = {}'.format(params.estimator)
    ]

    report.extend([
        '--- Units ---',
        'Unit system name = {}'.format(params.unit_type),
        'hbar = {}'.format(params.umap['hbar']),
        #'kB = {}'.format(params.umap['kB']),
        'amu = {}'.format(params.umap['amu']),
        '-- Conversion factors',
        'eV = {}'.format(params.umap['eV']),
        #'Kelvin = {}'.format(params.umap['kelvin']),
        'femtoseconds = {}'.format(params.umap['femtosec']),
        #'Hertz = {}'.format(params.hertz),
        #'Hertz*Radian (rad/s) = {}'.format(params.umap['hertz*rad']),
        'picometer = {}'.format(params.umap['picometer']),
        #'m/s = {}'.format(params.umap['m/s'])
    ])

    report.append('--- Miscellaneous ---')

    if runtime > 0:
        report.append('Runtime = {}'.format(time.strftime('%H:%M:%S', time.gmtime(runtime))))

    for key, value in extra.items():
        report.append('{} = {}'.format(key, value))

    report.append('--- Results ---')

    for x, y in zip(x_arr, y_arr):
        report.append('f[{}]={}'.format(x, y))

    if timestamp == -1:
        timestamp = int(time.time())

    with open("{}/{}.txt".format(dir, timestamp), "w") as file:
        file.write('\n'.join(report))

    is_bad = input("Is the result bad (Y/N)? ")

    if is_bad.lower() == 'y':
        shutil.move("{}/{}.png".format(dir, timestamp), "{}/bad/{}.png".format(dir, timestamp))
        shutil.move("{}/{}.txt".format(dir, timestamp), "{}/bad/{}.txt".format(dir, timestamp))
    elif is_bad.lower() == 'd':
        os.remove("{}/{}.png".format(dir, timestamp))
        os.remove("{}/{}.txt".format(dir, timestamp))