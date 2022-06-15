class UnitSystem:
    def __init__(self, unit_type='default'):
        self.type = unit_type
        self.factors = self.get_conversion_factors()

    def get_conversion_factors(self):
        if self.type == 'atomic':
            return {
                'hbar': 1.0,
                'kB': 1.0,
                'amu': 1822.8885,
                'eV': 0.036749326,
                'kelvin': 3.1668152e-06,
                'femtosec': 41.341373,
                'hertz*rad': 2.4188843e-17,
                'picometer': 1.8897261e-2,
                'm/s': 4.5710289e-7
            }
        elif self.type == 'old':
            return {
                'hbar': 0.6582,  # meV*ps
                # kB = 1.38e-23  # J/K
                'kB': 0.08617,  # meV/K
                'amu': 1.0,
            }