import math
from constants import PLANCK, LIGHT_SPEED

# turns a linear ratio into decibels
def convert_to_db(linear):
    return 10.0 * math.log10(linear)

# turns a decibel value back into a linear ratio
def convert_from_db(db):
    return 10**(db/10)

# turns a power in watts into dBm (dB relative to 1 mW)
def watts_to_dbm(power_watts):
    return 10 * math.log10(power_watts/ 1e-3)

# returns the energy of one photon in joules at the given wavelength in metres, E = hc/lambda
def photon_energy(wavelength):
    return PLANCK*LIGHT_SPEED /wavelength
