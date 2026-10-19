# puts model/ on the import path, it is only imported for that side effect
import _bootstrap

import math

from scipy.optimize import brentq

from coherent import (bpsk_heterodyne_ber, bpsk_homodyne_ber,
                      max_linewidth_for_dpsk, phase_error_variance_dpsk,
                      psk_ber_floor)
from modulation import ppm_ber_from_symbol_error, ppm_symbol_error_poisson

# lasers to compare, as (name, linewidth in Hz) for a single laser
# transmitter and local oscillator are a matched pair, so the combined linewidth is double
# (equation register 8.9)
LASERS = [
    # ordinary DFB laser at the wide end of the published range
    # Reference: Kikuchi (2016)
    ("DFB, 10 MHz", 1e7),
    # DFB laser at the narrow end of that same range
    # Reference: Kikuchi (2016)
    ("DFB, 100 kHz", 1e5),
    # stabilised laser from the project brief, the source planned for Objective 4
    ("stabilised, 100 Hz", 1e2),
]


# lowest bit rate (bit/s) a matched laser pair can run DPSK at
def lowest_dpsk_rate(linewidth_hz):
    # the linewidth limit grows in step with rate, so it scales from the limit at 1 bit/s
    return 2.0 * linewidth_hz / max_linewidth_for_dpsk(1.0)


# finds how many detected photons per bit (or per symbol) give the target BER
def _photons_for_ber(ber_of_n, target=1e-5):
    # brentq searches between 0.1 and 100 photons for where the BER meets the target
    return brentq(lambda n: ber_of_n(n) - target, 0.1, 100.0)


# prints the sensitivity comparison, the DPSK feasibility table and the lowest DPSK rates
# deep space rates are low, so the laser phase drifts a lot within one symbol
def print_summary():
    # photons per bit for BER 1e-5 with homodyne and heterodyne BPSK, worked out not typed in
    # the paper gives 4.55 and 9.0 for these two, which is what they get compared with
    # Reference: Kikuchi (2016), Equations 60 and 74
    homodyne = _photons_for_ber(bpsk_homodyne_ber)
    heterodyne = _photons_for_ber(bpsk_heterodyne_ber)
    # photons per symbol for uncoded 16-PPM with no background light
    # Reference: Hemmati et al. (2011), Equation A7
    ks = _photons_for_ber(lambda n: ppm_ber_from_symbol_error(
        ppm_symbol_error_poisson(n, 0.0, 16), 16))
    # prints the three sensitivities, PPM is divided by 4 bits per symbol to get per bit
    print("Shot-noise-limited sensitivity at BER 1e-5 (uncoded, detected "
          "photons/bit) [Kikuchi 2016 Eqs 60, 74]:")
    print(f"  BPSK homodyne     {homodyne:.2f}")
    print(f"  BPSK heterodyne   {heterodyne:.2f} (also intradyne/digital "
          "coherent; the paper rounds to 9.0)")
    print(f"  16-PPM counting   {ks / 4:.2f} (noiseless Poisson, "
          f"Hemmati 2011 A7: {ks:.1f} photons per symbol)")
    print()
    # prints the title and header row of the DPSK table, one column per laser
    print("DPSK phase-noise feasibility across deep-space data rates")
    print(f"{'rate':>10}  {'max TX+LO linewidth':>20}  "
          + "  ".join(f"{name:>20}" for name, _ in LASERS))
    # loops through deep space data rates from 10 kb/s up to 100 Mb/s
    for rate in (1e4, 1e5, 1e6, 1e8):
        # widest combined linewidth DPSK can cope with at this rate
        lw_max = max_linewidth_for_dpsk(rate)
        # starts the row with the rate in Mb/s and that linewidth limit
        row = f"{rate/1e6:>8.2f} M  {lw_max:>18.3g} Hz"
        # checks each laser pair against the limit
        for _, lw in LASERS:
            # rms phase error in radians that the laser pair builds up over one bit
            sigma = math.sqrt(phase_error_variance_dpsk(lw, lw, rate))
            # BER floor that phase error causes, even with loads of signal
            floor = psk_ber_floor(sigma)
            # OK if inside the limit, marginal if the floor is still under 1e-4, otherwise FAILS
            ok = "OK" if 2 * lw <= lw_max else (
                "marginal" if floor < 1e-4 else "FAILS")
            # adds this laser's verdict to the row
            row += f"  {ok:>20}"
        # prints the finished row for this rate
        print(row)
    print()
    # prints the lowest DPSK rate each laser pair manages with 0.2 rad of phase error or less
    print("Lowest DPSK bit rate each matched pair supports "
          "(sigma_phi <= 0.2 rad):")
    for name, lw in LASERS:
        print(f"  {name:20s} {lowest_dpsk_rate(lw):.3g} bit/s")


# runs the summary when the file is started directly
if __name__ == "__main__":
    print_summary()


# Academic Citations & Literature References:
#
# Coherent Receiver Sensitivity & Laser Linewidths:
# Kikuchi, K. (2016). Fundamentals of Coherent Optical Fiber Communications. Journal of
# Lightwave Technology, 34(1), 157-179. (Specifically Equations 60 and 74 for the BPSK
# homodyne and heterodyne BER, and the carrier-phase estimation discussion for the
# 100 kHz to 10 MHz DFB laser linewidths).
#
# Photon-Counting PPM Sensitivity:
# Hemmati, H., Biswas, A., & Djordjevic, I. B. (2011). Deep-Space Optical Communications:
# Future Perspectives and Applications. Proceedings of the IEEE, 99(11), 2020-2039.
# (Specifically Equation A7 for uncoded 16-PPM on the noiseless Poisson channel).
