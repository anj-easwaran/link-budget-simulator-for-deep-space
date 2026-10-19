import functools
import math

from modulation import ppm_capacity_bits_per_slot

# SCPPM settings, the code rates and PPM orders, codeword length and overhead in bits,
# and the guard time as a fraction of the signalling time
# Reference: CCSDS (2019)
SCPPM_CODE_RATES = (1.0 / 3.0, 1.0 / 2.0, 2.0 / 3.0)
SCPPM_PPM_ORDERS = (4, 8, 16, 32, 64, 128, 256)
SCPPM_CODEWORD_BITS = 15120
SCPPM_OVERHEAD_BITS = 34
SCPPM_GUARD_FRACTION = 0.25
# gap in dB between the code and capacity
# Reference: Hemmati et al. (2011), Section III
SCPPM_CAPACITY_GAP_DB = 1.0

# slot widths in s that get searched when the caller does not give any. DSOC shows 0.5 to 2 ns
# and this goes on to 8 ns. 0.25 ns and below is ASSUMED not available since no flight modem
# is cited at that width. searching only every other order (4, 16, 64, 256) can be about
# 3 dB too pessimistic, so the full order set is the default and the threshold is cached
SCPPM_SLOT_TIMES_S = (0.5e-9, 1.0e-9, 2.0e-9, 4.0e-9, 8.0e-9)


# information bits per second of one SCPPM frame, with the guard time added as (1 + guard)
# which gives the same rates as the paper's own Tables I to III
# Reference: Valentini et al. (2021), Equation 20
def scppm_data_rate(code_rate, ppm_order, slot_time_s,
                    guard_fraction=SCPPM_GUARD_FRACTION):
    # number of PPM symbols in one codeword
    symbols = SCPPM_CODEWORD_BITS / math.log2(ppm_order)
    # time taken by one frame including the guard slots
    frame_time = symbols * ppm_order * slot_time_s * (1.0 + guard_fraction)
    # information bits (codeword bits times code rate minus overhead) over the frame time
    return (SCPPM_CODEWORD_BITS * code_rate - SCPPM_OVERHEAD_BITS) / frame_time


# checks if the coded link closes, the code rate has to fit under the PPM capacity with the
# signal cut down by the code's gap to capacity (about 1 dB for SCPPM class codes)
# returns (supportable, data rate in bps, capacity in bits per slot)
# Reference: Hemmati et al. (2011), Section III
def scppm_supportable(code_rate, ppm_order, slot_time_s,
                      signal_photons_per_slot, background_photons_per_slot,
                      capacity_gap_db=SCPPM_CAPACITY_GAP_DB,
                      samples=200_000):
    # linear factor that cuts the signal down by the gap
    derate = 10.0 ** (-capacity_gap_db / 10.0)
    # PPM capacity with the reduced signal
    capacity = ppm_capacity_bits_per_slot(
        signal_photons_per_slot * derate, background_photons_per_slot,
        ppm_order, samples=samples)
    # data rate for this code rate, order and slot width
    rate_bps = scppm_data_rate(code_rate, ppm_order, slot_time_s)
    # code bits per signalling slot. guard slots cost data rate but not capacity because
    # they are dead time, so the check is done per signalling slot
    code_bits_per_slot = code_rate * math.log2(ppm_order) / ppm_order
    return code_bits_per_slot <= capacity, rate_bps, capacity


# smallest signal counts per PPM symbol K_s where the coded link still decodes, found by
# bisection in log K_s (capacity only goes up with K_s), cached so sweeps stay quick
# the counts are detected counts, so receiver efficiencies cancel out of any comparison
# this is the decode threshold flight papers quote, DSOC's design point is 0.93 at 20.6 Mbps
# Reference: Biswas & Srinivasan (2026)
@functools.lru_cache(maxsize=None)
def scppm_decode_threshold_ks(ppm_order, code_rate,
                              background_photons_per_slot,
                              capacity_gap_db=SCPPM_CAPACITY_GAP_DB,
                              samples=20_000, iterations=40):
    # bits per slot the code needs to carry
    target = code_rate * math.log2(ppm_order) / ppm_order
    # linear factor for the gap to capacity
    derate = 10.0 ** (-capacity_gap_db / 10.0)

    # checks if the link decodes at signal counts ks, the capacity has to reach the target
    def decodes(ks):
        return ppm_capacity_bits_per_slot(
            ks * derate, background_photons_per_slot, ppm_order,
            samples=samples) >= target

    # starting lower and upper limits for K_s
    lo, hi = 1e-3, 100.0
    # loops the bisection, closing in on the threshold each time
    for _ in range(iterations):
        # midpoint of the two limits on a log scale
        mid = math.sqrt(lo * hi)
        # if it decodes at mid, the threshold is at or below mid
        if decodes(mid):
            hi = mid
        # otherwise the threshold is above mid
        else:
            lo = mid
    # returns the upper limit so the answer always decodes
    return hi


# smallest detected photon flux in photons/s that decodes rate_bps, the lowest over every
# signalling setup that can reach that rate. it is the reverse of scppm_max_data_rate,
# which asks how fast for a given flux, and it returns None if nothing reaches the rate
def scppm_flux_needed_for_rate(rate_bps, dark_count_rate_hz=0.0,
                               capacity_gap_db=SCPPM_CAPACITY_GAP_DB,
                               ppm_orders=SCPPM_PPM_ORDERS,
                               slot_times_s=SCPPM_SLOT_TIMES_S,
                               samples=20_000):
    # no answer found yet
    best = None
    # loops through every code rate, PPM order and slot width
    for code_rate in SCPPM_CODE_RATES:
        for order in ppm_orders:
            for slot in slot_times_s:
                # if this setup cannot reach the rate at all, skip it
                if scppm_data_rate(code_rate, order, slot) < rate_bps:
                    continue
                # counts per symbol needed to decode, with dark counts as the background
                ks = scppm_decode_threshold_ks(
                    order, code_rate, dark_count_rate_hz * slot,
                    capacity_gap_db, samples=samples)
                # the pulsed slot carries the whole symbol's average power, so the flux is the
                # counts over the full symbol time with the guard time included
                flux = ks / (order * slot * (1.0 + SCPPM_GUARD_FRACTION))
                # keeps this flux if it is the lowest so far
                if best is None or flux < best:
                    best = flux
    return best


# slot widths going up in half octave steps from min_slot_time_s, longest first
# for a receiver free to pick any slot width, a real modem has a fixed set instead
# which should be passed in as slot_times_s
def scppm_slot_ladder(min_slot_time_s, steps=25):
    return tuple(min_slot_time_s * 2.0 ** (step / 2.0)
                 for step in reversed(range(steps)))


# searches every code rate, PPM order and slot width for the fastest SCPPM rate that closes
# ppm_orders and slot_times_s are what the receiver can really do. left as defaults the
# search is free (any order, any slot from min_slot_time_s up), and passing a mission's
# demonstrated set keeps it to real hardware, which matters at long range where the free
# best case needs slot widths no flight modem has
# returns (rate_bps, code_rate, ppm_order, slot_time_s) or (0.0, None, None, None) if none close
# Reference: Xie et al. (2016), Algorithm 1
def scppm_max_data_rate(signal_photon_flux_hz, background_flux_hz,
                        detection_eff, min_slot_time_s,
                        dark_count_rate_hz=0.0,
                        capacity_gap_db=SCPPM_CAPACITY_GAP_DB,
                        ppm_orders=SCPPM_PPM_ORDERS,
                        slot_times_s=None,
                        slot_time_steps=25,
                        samples=20_000):
    # if no slot widths are given, make the half octave ladder
    if slot_times_s is None:
        slot_times_s = scppm_slot_ladder(min_slot_time_s, slot_time_steps)
    # otherwise sorts the given widths longest first
    else:
        slot_times_s = tuple(sorted(slot_times_s, reverse=True))
    # best result so far, starts with nothing closing
    best = (0.0, None, None, None)
    for code_rate in SCPPM_CODE_RATES:
        for order in ppm_orders:
            # slots go longest first. the rate goes up as slots shrink but closing only gets
            # harder, so setups that cannot beat the best skip the monte carlo check and
            # the first failure ends the scan
            for slot in slot_times_s:
                # if this setup cannot beat the best rate so far, skip it
                if scppm_data_rate(code_rate, order, slot) <= best[0]:
                    continue
                # signal counts in the pulsed slot, which carries the average power of the
                # whole symbol including the guard time
                ks = (detection_eff * signal_photon_flux_hz
                      * slot * order * (1.0 + SCPPM_GUARD_FRACTION))
                # background and dark counts in each slot
                kb = (detection_eff * background_flux_hz
                      + dark_count_rate_hz) * slot
                # checks if this setup closes
                ok, rate, _ = scppm_supportable(
                    code_rate, order, slot, ks, kb, capacity_gap_db,
                    samples=samples)
                # if it does not close, shorter slots will not either so stop this order
                if not ok:
                    break
                # if it closes faster than the best so far, keep it
                if rate > best[0]:
                    best = (rate, code_rate, order, slot)
    return best


# Academic Citations & Literature References:
#
# SCPPM Coded Modulation:
# CCSDS. (2019). Optical Communications Coding and Synchronization, Recommended Standard
# CCSDS 142.0-B-1 (Blue Book). Consultative Committee for Space Data Systems.
# (Specifically the SCPPM coded modulation scheme).
#
# SCPPM Frame Data Rate:
# Valentini, L., Faedi, A., Paolini, E., & Chiani, M. (2021). Analysis of Pointing Loss
# Effects in Deep Space Optical Links. Proc. IEEE Global Communications Conf. (GLOBECOM).
# (Specifically Equation 20, with the guard time applied as (1 + guard fraction)).
#
# Code Gap to Capacity:
# Hemmati, H., Biswas, A., & Djordjevic, I. B. (2011). Deep-Space Optical Communications:
# Future Perspectives and Applications. Proceedings of the IEEE, 99(11), 2020-2039.
# (Specifically Section III, about 1 dB for SCPPM class codes on the Poisson PPM channel).
#
# DSOC Decode Threshold:
# Biswas, A., & Srinivasan, M. (2026). Deep Space Optical Communications (DSOC) System
# Description and Performance. IEEE Journal of Selected Topics in Quantum Electronics, 32(1).
# (Specifically the design point of 0.93 signal counts per symbol at 20.6 Mbps).
#
# Rate Search Procedure:
# Xie, H., Heckman, D., & Breidenthal, J. (2016). Link Characterization for Deep-Space
# Optical Communications (IPN Progress Report 42-205). JPL. (Specifically Algorithm 1, the
# SOLT procedure).
