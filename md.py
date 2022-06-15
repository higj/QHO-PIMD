import numpy as np
import params


class MDSimulation:
    """Single-particle Path-Integral Molecular Dynamics simulation for a 1D quantum harmonic oscillator.

    Attributes:
        beta: Thermodynamic beta
        dt: Time increment of the MD simulation
        mass: The mass of the quantum particle
        bead_num: The number of beads on the classical cyclic polymer chain
        size: The size of the system
        save_freq: Frequency of calculations (e.g., average energy)
        tsteps: Total number of time increments
        enable_thermostat: A boolean indicating if a Langevin thermostat should be enabled
        estimator_type:
            Specifies the type of estimator to be used for calculating the average energy.
            'default' corresponds to the primitive energy estimator (Eqn. 12.3.20 in Tuckerman)
            'virial' corresponds to the virial energy estimator (Eqn. 12.6.34 in Tuckerman)
            'centroid_virial' corresponds to the centroid-virial (kinetic) energy estimator (J. Chem. Phys., 76, 5150, 1982)
        threshold: Specifies the threshold value for thermalization (between 0 and 1)
        mean_energy: The mean energy calculated during the simulation using an estimator specified by est
    """

    def __init__(self, beta, dt, mass, bead_num=6, size=1, save_freq=1, tsteps=1000, enable_thermostat=True, estimator_type='default', threshold=0.1):
        self.beta = beta
        self.p = bead_num
        self.dt = dt
        # Create an array of masses (in case we want to modify individual masses later)
        self.mass_arr = np.full(bead_num, mass)
        self.size = size
        self.save_freq = save_freq
        self.tsteps = int(tsteps)
        self.enable_thermostat = enable_thermostat
        self.est_type = estimator_type
        self.threshold = threshold*tsteps

        self.sim_time = self.dt * self.tsteps  # Total simulation time

        # Generate initial conditions
        self.pos_arr = self.gen_pos()  # Current positions of the beads
        self.momenta = self.mass_arr * self.sample_mb()  # Beads momenta (velocities sampled from the MB distribution)

        # Initialize forces
        self.force_arr = np.zeros_like(self.pos_arr)

        self.mean_energy = 0

    def run(self):
        """Run an MD simulation by performing a series of velocity Verlet steps."""

        for step in range(self.tsteps):
            if self.enable_thermostat:
                self.langevin_step()  # Perform a Langevin step (thermostat)

            self.vv_step()  # Perform a velocity Verlet step

            if self.enable_thermostat:
                self.langevin_step()  # Perform a Langevin step (thermostat)

            if step < self.threshold:
                continue

            # Perform calculations every nth step (where n=save_freq)
            if step % self.save_freq == 0:
                # Accumulate contributions to the average energy
                # (Based on Eqn. 2.37 from "Path Integral Methods in Atomistic Modelling"
                # by Ceriotti, Manolopoulos, Markland and Rossi; "CMMR" for short)
                self.mean_energy += self.total_energy_estimator()

        # Integration begins after thermalization
        self.mean_energy *= self.save_freq / (self.tsteps - self.threshold)

    def gen_pos(self):
        """Generate random positions for the initial conditions by sampling a uniform distribution."""

        return np.random.uniform(0, self.size*params.umap['picometer'], self.p)

    def sample_mb(self):
        """Sample the Maxwell-Boltzmann distribution and return an array of velocities."""

        return np.random.normal(scale=1/np.sqrt(self.beta * self.mass_arr))

    def langevin_step(self):
        """Implement a thermostat according to the Langevin scheme. This is necessary in order to
            sample a canonical ensemble, as opposed to the NVE ensemble."""

        noise = np.random.normal(size=self.p)  # The noise term (time derivative of a Wiener process)
        a = np.exp(-0.5 * params.gamma * self.dt)
        b = np.sqrt((1 - a**2) * self.mass_arr / self.beta)
        self.momenta = a * self.momenta + b * noise

    def vv_step(self):
        """Propagate position and velocity according to the symmetric-split velocity Verlet algorithm."""

        # First step: momenta are propagated half a step
        self.momenta += 0.5 * self.dt * self.force_arr
        # Second step: positions are propagated using the new velocities
        self.pos_arr += self.dt * self.momenta / self.mass_arr
        # Third step: forces are updated using the new positions
        self.update_forces()
        # Fourth step: momenta are propagated once more
        self.momenta += 0.5 * self.dt * self.force_arr

    def update_forces(self):
        """Update forces for all the beads."""

        for i in range(self.p):
            self.force_arr[i] = self.get_total_force(i)

    def get_spring_force(self, bead_num):
        """Returns the force exerted on a specific bead by the two neighboring beads (the beads are enumerated
        from 0 to P-1). The force is given by Eqn. (12.6.4) in Tuckerman.
        """

        chain_frequency = np.sqrt(self.p)/(self.beta*params.umap['hbar'])  # Self-explanatory
        prefactor = self.mass_arr[bead_num] * chain_frequency**2  # Nearest-neighbor coupling constant
        # Only interactions with the two closest neighbors are included.
        # The modulo operation reflects the periodic boundary conditions.
        return prefactor * (self.pos_arr[(bead_num + 1) % self.p] + self.pos_arr[bead_num - 1]
                            - 2 * self.pos_arr[bead_num])

    def get_external_force(self, bead_num):
        """Returns the force exerted on a specific bead by the external potential.
        In this particular case the potential is that of a simple harmonic oscillator.
        """

        prefactor = self.mass_arr[bead_num] * params.qho_freq**2 / self.p
        return -prefactor*self.pos_arr[bead_num]

    def get_total_force(self, bead_num):
        """Calculate the total force acting on a given bead."""

        return self.get_spring_force(bead_num) + self.get_external_force(bead_num)

    def get_bead_pe(self, bead_num):
        """Get the (external) potential energy associated with a specific bead. (Currently not used)"""
        return 0.5 * self.mass_arr[bead_num] * (params.qho_freq * self.pos_arr[bead_num])**2

    def potential_energy_estimator(self):
        """Returns the potential energy estimator. Based on Eqn. 2.53 in CMMR.
        """
        pe_arr = 0.5 * params.qho_freq**2 * (self.mass_arr * self.pos_arr**2)
        return np.mean(pe_arr)

    def kinetic_energy_estimator(self):
        """Returns the primitive kinetic energy estimator. Based on Eqn. (3.3) in CMMR.
        """

        prefactor = 0.5 * self.p / self.beta

        c = self.mass_arr / (self.beta * params.umap['hbar']**2)

        corr_term = 0  # Corresponds to the term describing the cross-correlations between different replicas

        for bead_i in range(self.p):
            corr_term += c[bead_i] * (self.pos_arr[bead_i] - self.pos_arr[(bead_i + 1) % self.p])**2

        return prefactor * (1 - corr_term)

    def centroid_virial_ke(self):
        """Returns the centroid virial kinetic energy estimator. Based on Eqn. (3.15) in CMMR.
        """

        a = 1 / (2 * self.beta)  # The first term in the equation
        qbar = np.mean(self.pos_arr)
        b = 0  # Holds the second term

        for bead_i in range(self.p):
            dq = self.pos_arr[bead_i] - qbar
            force = - self.p * self.get_external_force(bead_i)
            b += dq * force

        return a + b / (2 * self.p)

    def virial_energy_estimator(self):
        """Returns the virial energy estimator. Based on Eqn. (12.6.34) in Tuckerman.
        """

        virial_e = 0
        for bead_i in range(self.p):
            grad = - self.p * self.get_external_force(bead_i)
            pe = 0.5 * self.mass_arr[bead_i] * (params.qho_freq**2) * (self.pos_arr[bead_i])**2
            virial_e += 0.5 * self.pos_arr[bead_i] * grad + pe

        return virial_e / self.p

    def total_energy_estimator(self):
        """Returns the total energy estimator of the desired type."""

        if self.est_type == 'centroid_virial':
            return self.potential_energy_estimator() + self.centroid_virial_ke()
        elif self.est_type == 'virial':
            return self.virial_energy_estimator()
        elif self.est_type == 'simple_virial_pe':
            # Only for QHO: the mean of E is twice the PE
            return 2 * self.potential_energy_estimator()
        elif self.est_type == 'simple_virial_ke':
            # Only for QHO: the mean of E is twice the KE
            return 2 * self.kinetic_energy_estimator()

        return self.potential_energy_estimator() + self.kinetic_energy_estimator()
