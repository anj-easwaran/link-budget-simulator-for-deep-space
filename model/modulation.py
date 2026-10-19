import math

import numpy as np
from scipy.stats import poisson


# bit error rate of OOK with photon counting and a maximum likelihood threshold,
# the inputs are mean photon counts per slot
# Reference: Hemmati et al. (2011), Equations A1 to A3
def ook_poisson_ber(signal_photons, background_photons):
    ks, kb = signal_photons, background_photons
    # if there is no background light, only a missed pulse causes an error
    if kb == 0.0:
        return 0.5 * math.exp(-ks)
    # count threshold that decides between a 0 and a 1
    threshold = math.floor(ks / math.log1p(ks / kb))
    # every count from 0 up to the threshold
    k = np.arange(threshold + 1)
    # chance of each count when a 0 is sent (background only)
    f0 = poisson.pmf(k, kb)
    # chance of each count when a 1 is sent (signal plus background)
    f1 = poisson.pmf(k, ks + kb)
    return 0.5 - 0.5 * float(np.sum(f0 - f1))


# exact symbol error of uncoded M-PPM, maximum likelihood with ties broken at random
# signal_photons is K_s in the pulsed slot and background_photons is K_b in every slot,
# both already include the detection efficiency and dark counts
# Reference: Hemmati et al. (2011), Equations A5 to A7
def ppm_symbol_error_poisson(signal_photons, background_photons, order):
    ks, kb, m = signal_photons, background_photons, order
    # if there is no background, an error only happens when the pulse gives zero counts
    if kb == 0.0:
        return (m - 1.0) / m * math.exp(-ks)
    # highest count worth summing to, well past the mean so the tail is tiny
    k_max = int((ks + kb) + 12.0 * math.sqrt(ks + kb) + 30)
    k = np.arange(k_max + 1)
    # chance of each count in the pulsed slot
    pmf1 = poisson.pmf(k, ks + kb)
    # chance a background only slot gets k counts or fewer
    cdf0 = poisson.cdf(k, kb)
    # same thing shifted by one, so k-1 counts or fewer
    cdf0_prev = np.concatenate(([0.0], cdf0[:-1]))
    # tie term written as a sum of products that are all at most 1, so it stays stable
    tie_sum = np.zeros_like(cdf0)
    # loops through j and adds each F^j times F_prev^(M-1-j) term
    for j in range(m):
        tie_sum += cdf0 ** j * cdf0_prev ** (m - 1 - j)
    # chance the pulsed slot is picked correctly
    p_correct = float(np.sum(pmf1 * tie_sum / m))
    return 1.0 - p_correct


# turns a PPM symbol error into a bit error rate, M/(2(M-1)) times P_s
# Reference: Li et al. (2018), Equation 3
def ppm_ber_from_symbol_error(symbol_error, order):
    return order / (2.0 * (order - 1.0)) * symbol_error


# soft decision capacity of M-PPM on the Poisson channel in bits per slot, done by monte carlo
# Reference: Faedi (2021), Equation 1.4
# Reference: Moision & Hamkins (2005)
def ppm_capacity_bits_per_slot(signal_photons, background_photons, order,
                               samples=200_000, seed=0):
    ks, kb, m = signal_photons, background_photons, order
    # if there is no signal at all, nothing can be sent
    if ks == 0.0:
        return 0.0
    # if there is no background, use the exact noiseless answer instead of simulating
    if kb == 0.0:
        return (1.0 - math.exp(-ks)) * math.log2(m) / m
    # random number generator with a fixed seed so results repeat
    rng = np.random.default_rng(seed)
    # log likelihood per count, the -Ks part cancels so only the y log(1+Ks/Kb) part is kept
    log_g = math.log1p(ks / kb)
    # random counts in the pulsed slot for every sample
    y1 = rng.poisson(ks + kb, samples)
    # random counts in the other M-1 background only slots
    y0 = rng.poisson(kb, (samples, m - 1))
    # log likelihoods for the pulsed slot and for the empty slots
    log_l1 = y1 * log_g
    log_l0 = y0 * log_g
    # puts all M slots side by side for each sample
    log_all = np.concatenate((log_l1[:, None], log_l0), axis=1)
    # log of the sum over slots, taking the biggest out first so exp does not overflow
    peak = log_all.max(axis=1)
    log_sum = peak + np.log(np.exp(log_all - peak[:, None]).sum(axis=1))
    # information in bits per symbol for each sample
    bits = (math.log2(m) + (log_l1 - log_sum) / math.log(2.0))
    # averages over the samples and divides by M to get bits per slot
    return float(bits.mean()) / m


# Academic Citations & Literature References:
#
# OOK & M-PPM Error Rates on the Poisson Channel:
# Hemmati, H., Biswas, A., & Djordjevic, I. B. (2011). Deep-Space Optical Communications:
# Future Perspectives and Applications. Proceedings of the IEEE, 99(11), 2020-2039.
# (Specifically Equations A1 to A3 for OOK and Equations A5 to A7 for the PPM symbol error).
#
# PPM Symbol to Bit Error Mapping:
# Li, B., et al. (2018). BER Analysis of a Deep Space Optical Communication System Based on
# SNSPD Over Double Generalized Gamma Channel. IEEE Photonics Journal, 10(5).
# (Specifically Equation 3).
#
# Soft-Decision PPM Capacity:
# Faedi, A. (2021). Design and Analysis of Optical Links for Space Communications (Master's
# thesis). Università di Bologna.
# Moision, B., & Hamkins, J. (2005). Coded Modulation for the Deep-Space Optical Channel:
# Serially Concatenated Pulse-Position Modulation (IPN Progress Report 42-161). JPL.
# (Specifically Equation 1.4 for the capacity worked out by Monte Carlo).
