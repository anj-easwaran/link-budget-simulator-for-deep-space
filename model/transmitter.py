import math

# holds the laser and transmit telescope settings and works out gain and divergence
class Transmitter:
    # stores laser power in W, wavelength in m, optics efficiency and pointing error in rad
    # the beam is set either by the aperture diameter or by a full divergence angle
    def __init__(self,power_w,wavelength_m,optics_eff=1,pointing_error_rad=0,
                 aperture_diameter_m=None, obscuration_ratio=0,full_divergence_rad=None,
                 truncation_ratio=None, divergence_convention="liang"):
        # if the divergence convention is not one of the two known ones, stop with an error
        if divergence_convention not in ("liang", "1/e2"):
            raise ValueError("divergence_convention is 'liang' or '1/e2'")
        # checks that exactly one of aperture or divergence was given, not both or neither
        if (aperture_diameter_m is None) == (full_divergence_rad is None):
            raise ValueError("give exactly one of aperture_diameter_m or full_divergence_rad")
        self.power_w=power_w
        self.wavelength_m=wavelength_m
        self.optics_eff=optics_eff
        self.pointing_error_rad = pointing_error_rad
        self.aperture_diameter_m = aperture_diameter_m
        self.obscuration_ratio = obscuration_ratio
        self.full_divergence_rad = full_divergence_rad
        self.truncation_ratio = truncation_ratio
        self.divergence_convention = divergence_convention

    # returns the linear transmit antenna gain, from the divergence if one was given
    # or from the aperture and beam truncation otherwise
    def get_gain(self):
        # if a divergence angle was given, work the gain out straight from it
        if self.full_divergence_rad is not None:
            # if it is the full width between the 1/e^2 points, 8 over the half angle squared
            # which is 32 over the full angle squared (docs/equations.md 2.5)
            if self.divergence_convention == "1/e2":
                return 32.0 / self.full_divergence_rad ** 2
            # otherwise use Liang's 16/theta^2, only right for widths quoted their way (2.4)
            # Reference: Liang et al. (2022)
            return 16.0 / self.full_divergence_rad ** 2
        # radius and area of the transmit aperture
        a = self.aperture_diameter_m / 2.0
        area = math.pi * a ** 2
        # obscuration ratio, secondary mirror size over primary mirror size
        gamma = self.obscuration_ratio
        # if a measured truncation ratio is published, use that one
        if self.truncation_ratio is not None:
            alpha = self.truncation_ratio
        # otherwise use the optimum truncation ratio for this obscuration (docs/equations.md 2.2)
        else:
            alpha = 1.12 - 1.30 * gamma ** 2 + 2.12 * gamma ** 4
        # gain factor for a truncated Gaussian beam with a central obscuration
        g_t = (2.0 / alpha ** 2) * (math.exp(-alpha ** 2)
                                    - math.exp(-alpha ** 2 * gamma ** 2)) ** 2
        # aperture gain 4 pi A / lambda^2 times the truncation factor (docs/equations.md 2.1)
        return (4.0 * math.pi * area / self.wavelength_m ** 2) * g_t

    # returns the full 1/e^2 divergence angle of the beam in radians
    def get_divergence(self):
        # if a divergence was given, just return it
        if self.full_divergence_rad is not None:
            return self.full_divergence_rad
        # the half angle goes from 0.94 lambda/D with no obscuration to 0.79 lambda/D at 0.4,
        # and the paper says to interpolate in between
        # Reference: Klein & Degnan (1974), Figure 10
        gamma = self.obscuration_ratio
        factor = 0.94 + (0.79 - 0.94) * (gamma / 0.4)
        # doubles the half angle to get the full divergence
        return 2.0 * factor * self.wavelength_m / self.aperture_diameter_m


# Academic Citations & Literature References:
#
# Beam Divergence of the Transmit Aperture:
# Klein, B. J., & Degnan, J. J. (1974). Optical Antenna Gain. 1: Transmitting Antennas.
# Applied Optics, 13(9), 2134-2141. (Specifically Figure 10 for the far-field half-angle
# of an unobscured and an obscured aperture).
#
# Divergence-Based Gain Convention:
# Liang, J., Chaudhry, A. U., Erdogan, E., & Yanikomeroglu, H. (2022). Link Budget
# Analysis for Free-Space Optical Satellite Networks. Proc. IEEE 23rd Int. Symp. on a
# World of Wireless, Mobile and Multimedia Networks (WoWMoM). (Specifically the 16/theta^2
# gain for divergences quoted in their convention).
