import math

from scipy.special import j1


# linear pointing loss for a Gaussian beam, exp(-G theta^2), the worse case of the two models
# Reference: Valentini et al. (2021), Equation 3
def gaussian_pointing_loss(gain_linear, error_rad):
    return math.exp(-gain_linear * error_rad ** 2)


# linear pointing loss for a circular (Airy) aperture using the Bessel J1 function
# Reference: Valentini et al. (2021), Equation 4
def airy_pointing_loss(gain_linear, error_rad):
    # pointing error scaled by the square root of the gain
    x = math.sqrt(gain_linear) * error_rad
    # if there is no pointing error, return 1 so there is no loss (and no divide by zero)
    if x == 0.0:
        return 1.0
    return (2.0 * j1(x) / x) ** 2


# best linear antenna gain for a given worst pointing error in radians, 1/theta^2
# Reference: Valentini et al. (2021), Equation 17
def optimum_gain(max_error_rad):
    return 1.0 / max_error_rad ** 2


# Academic Citations & Literature References:
#
# Pointing Loss & Optimum Gain:
# Valentini, L., Faedi, A., Paolini, E., & Chiani, M. (2021). Analysis of Pointing Loss
# Effects in Deep Space Optical Links. Proc. IEEE Global Communications Conf. (GLOBECOM).
# (Specifically Equations 3 and 4 for the two pointing losses and Equation 17 for the
# optimum gain).
