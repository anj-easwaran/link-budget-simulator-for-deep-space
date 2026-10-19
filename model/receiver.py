import math


# holds the receive telescope settings and works out its gain
class Receiver:
    # stores the aperture size in m, optics efficiency, pointing error in rad and obscuration
    # sensitivity_w is the published power in W needed for the target BER at the data rate,
    # which OOK links use, photon counting links use the photon flux instead
    # Reference: Liang et al. (2022)
    def __init__(self,
                 aperture_diameter_m,
                 optics_eff=1.0,
                 pointing_error_rad=0.0,
                 sensitivity_w=None,
                 obscuration_ratio=0.0):
        self.aperture_diameter_m = aperture_diameter_m
        self.optics_eff = optics_eff
        self.pointing_error_rad = pointing_error_rad
        self.sensitivity_w = sensitivity_w
        self.obscuration_ratio = obscuration_ratio

    # returns the linear receive antenna gain (pi D / lambda)^2 minus the blocked area
    def get_gain(self, wavelength_m):
        # the secondary mirror blocks this fraction of the collecting area (ratio squared)
        blocked = self.obscuration_ratio ** 2
        return ((math.pi * self.aperture_diameter_m / wavelength_m) ** 2
                * (1.0 - blocked))


# Academic Citations & Literature References:
#
# Links Closed Against a Published Receiver Sensitivity:
# Liang, J., Chaudhry, A. U., Erdogan, E., & Yanikomeroglu, H. (2022). Link Budget
# Analysis for Free-Space Optical Satellite Networks. Proc. IEEE 23rd Int. Symp. on a
# World of Wireless, Mobile and Multimedia Networks (WoWMoM). (Specifically the OOK links
# checked against a published receiver sensitivity).
