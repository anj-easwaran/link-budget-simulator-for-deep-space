import math

import numpy as np
from scipy.special import erfc, i0e
from scipy.stats import ncx2

# picks numpy's trapezoid function, or the older trapz name on older numpy versions
_trapezoid = getattr(np, "trapezoid", None) or np.trapz

# offsets either side of the Rician peak used for the DPSK envelope integral
_DPSK_EPS_OFFSETS = np.linspace(-10.0, 10.0, 241)


# BER of homodyne BPSK limited only by shot noise, 1/2 erfc(sqrt(2 n_s)) for n_s photons per bit
# Reference: Kikuchi (2016), Equation 60
def bpsk_homodyne_ber(photons_per_bit):
    return 0.5 * erfc(math.sqrt(2.0 * photons_per_bit))


# BER of heterodyne BPSK, 1/2 erfc(sqrt(n_s)), which also covers intradyne digital coherent
# Reference: Kikuchi (2016), Equation 74
def bpsk_heterodyne_ber(photons_per_bit):
    return 0.5 * erfc(math.sqrt(photons_per_bit))


# phase error variance sigma_phi^2 in rad^2 for delay demodulated DPSK over one bit period
# Equation 62 as printed is off from Equation 60 by (2 pi)^2 and only Equation 60 gives
# the paper's 3 Gbit/s result, so Equation 60 is used (docs/equations.md 8.9)
# Reference: Kikuchi et al. (1984), Equation 60
def phase_error_variance_dpsk(tx_linewidth_hz, lo_linewidth_hz, bit_rate_hz):
    return 2.0 * math.pi * (tx_linewidth_hz + lo_linewidth_hz) / bit_rate_hz


# error floor at infinite SNR, the chance that the phase error is bigger than pi/2
# Reference: Kikuchi et al. (1984), Equations 19 and 20, Figure 4
def psk_ber_floor(sigma_phi):
    return erfc(math.pi / (2.0 * math.sqrt(2.0) * sigma_phi))


# DPSK BER at a fixed phase error phi before it is averaged over phi, as 4 times a double
# integral with alpha = 2 gamma cos^2(phi/2) and beta = 2 gamma sin^2(phi/2)
# Reference: Kikuchi et al. (1984), Equations 36 to 38
def _dpsk_conditional_ber(snr_gamma, phi):
    # phi as a column so every phase value is worked out across all the eps points at once
    phi = np.atleast_1d(np.asarray(phi, dtype=float))[:, None]
    # alpha and beta terms for each phase error
    alpha = 2.0 * snr_gamma * np.cos(phi / 2.0) ** 2
    beta = 2.0 * snr_gamma * np.sin(phi / 2.0) ** 2
    # position of the Rician peak
    peak = np.sqrt(alpha)
    # eps points spread out around the peak
    eps = peak + _DPSK_EPS_OFFSETS[None, :]
    # only eps above zero is part of the integral
    inside = eps > 0.0
    # sets the eps points below zero to zero so the maths below does not break
    eps = np.where(inside, eps, 0.0)
    # inner xi integral, a Marcum Q function written as a non-central chi-square survival
    # function (docs/equations.md 8.11). scipy overflows for tiny x but the answer is 1 there
    # whatever the non-centrality, so x is floored at 1e-6
    xi_integral = 0.5 * ncx2.sf(np.maximum(2.0 * eps ** 2, 1e-6), 2,
                                np.maximum(2.0 * beta, 1e-12))
    # eps I0(2 sqrt(alpha) eps) exp(-eps^2 - alpha), using the scaled Bessel function
    # so it cannot overflow
    envelope = eps * i0e(2.0 * peak * eps) * np.exp(-(eps - peak) ** 2)
    # integrand, set to zero outside the valid eps range
    integrand = np.where(inside, envelope * xi_integral, 0.0)
    # spacing between the eps points
    step = _DPSK_EPS_OFFSETS[1] - _DPSK_EPS_OFFSETS[0]
    # trapezoid rule over eps, times 4 from the formula
    return 4.0 * _trapezoid(integrand, dx=step, axis=1)


# heterodyne DPSK BER with a Gaussian phase error of standard deviation sigma_phi, averaged
# over phi >= 0 on a fine grid since the errors come from the tail past pi/2
# with sigma_phi at 0 it comes down to 1/2 exp(-gamma) (docs/equations.md 8.11)
# Reference: Kikuchi et al. (1984), Equation 36
def dpsk_ber_with_phase_noise(snr_gamma, sigma_phi, points=401):
    # if the SNR is very high, only the phase error causes errors so the floor is used
    # (within 0.1% of the full integral at gamma = 1e5)
    # Reference: Kikuchi et al. (1984), Equation 39
    if snr_gamma > 1e5:
        return psk_ber_floor(sigma_phi) if sigma_phi > 0.0 else 0.0
    # if there is no phase noise, there is nothing to average over
    if sigma_phi == 0.0:
        return float(_dpsk_conditional_ber(snr_gamma, 0.0)[0])
    # grid of phase errors out to 8 standard deviations
    phi = np.linspace(0.0, 8.0 * sigma_phi, points)
    # Gaussian probability density, doubled since only phi >= 0 is used
    pdf = 2.0 * np.exp(-phi ** 2 / (2.0 * sigma_phi ** 2)) / (
        sigma_phi * math.sqrt(2.0 * math.pi))
    # averages the fixed phase BER over the phase error
    return float(_trapezoid(pdf * _dpsk_conditional_ber(snr_gamma, phi), phi))


# highest FEC coded heterodyne DPSK information rate in bits/s for a detected flux and laser pair
# slower rates pick up more phase noise per bit and faster rates get fewer photons per bit,
# so only a window of rates works. it steps down from max_rate_hz to the first rate under
# the FEC threshold, then bisects the top edge of the window, and returns 0.0 if none work
# Reference: Kikuchi et al. (1984), Equations 36 and 62
def dpsk_max_data_rate(detected_flux_hz, tx_linewidth_hz, lo_linewidth_hz,
                       target_ber=3.8e-3, fec_overhead=0.07,
                       max_rate_hz=1e12, min_rate_hz=1e3):
    # BER at a given bit rate, gamma is the detected photons per bit
    def ber(rate_hz):
        # phase error standard deviation at this rate
        sigma = math.sqrt(phase_error_variance_dpsk(
            tx_linewidth_hz, lo_linewidth_hz, rate_hz))
        return dpsk_ber_with_phase_noise(detected_flux_hz / rate_hz, sigma)

    # 80 rates spaced evenly in log from the fastest down to the slowest
    rates = np.geomspace(max_rate_hz, min_rate_hz, 80)
    # if even the fastest rate is under the 3.8e-3 FEC threshold, return it minus the overhead
    # Reference: Li et al. (2018)
    if ber(rates[0]) < target_ber:
        return max_rate_hz * (1.0 - fec_overhead)
    # steps down through the rates until one meets the threshold
    for hi, lo in zip(rates[:-1], rates[1:]):
        if ber(lo) < target_ber:
            break
    # if the loop never breaks, no rate meets the threshold
    else:
        return 0.0
    # bisects between the last failing rate and the first passing one to find the top edge
    for _ in range(30):
        mid = math.sqrt(lo * hi)
        # if mid still meets the threshold, the edge is above it
        if ber(mid) < target_ber:
            lo = mid
        # otherwise the edge is below mid
        else:
            hi = mid
    # information rate after taking off the 7% FEC overhead
    return lo * (1.0 - fec_overhead)


# largest combined TX plus LO linewidth in Hz that keeps the DPSK phase error under the budget
# the default 0.2 rad keeps the sensitivity loss at BER 1e-9 under 3 dB, and the formula
# is Equation 60 turned around
# Reference: Kikuchi et al. (1984), Section V-A and Equation 60
def max_linewidth_for_dpsk(bit_rate_hz, sigma_phi_budget=0.2):
    return sigma_phi_budget ** 2 * bit_rate_hz / (2.0 * math.pi)


# Academic Citations & Literature References:
#
# Shot-Noise-Limited Coherent BER:
# Kikuchi, K. (2016). Fundamentals of Coherent Optical Fiber Communications. Journal of
# Lightwave Technology, 34(1), 157-179. (Specifically Equation 60 for homodyne BPSK and
# Equation 74 for heterodyne BPSK).
#
# Laser Phase Noise & DPSK BER:
# Kikuchi, K., Okoshi, T., Nagamatsu, M., & Henmi, N. (1984). Degradation of Bit-Error Rate
# in Coherent Optical Communications Due to Spectral Spread of the Transmitter and the Local
# Oscillator. Journal of Lightwave Technology, 2(6), 1024-1033. (Specifically Equation 60 for
# the phase error variance, Equations 19, 20 and 36 to 39 for the BER, and Section V-A for
# the 0.2 rad phase error budget).
#
# FEC Threshold:
# Li, B., et al. (2018). BER Analysis of a Deep Space Optical Communication System Based on
# SNSPD Over Double Generalized Gamma Channel. IEEE Photonics Journal, 10(5).
# (Specifically the 3.8e-3 FEC threshold).
