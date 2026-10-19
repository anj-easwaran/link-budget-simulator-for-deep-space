# puts model/ on the import path, it is only imported for that side effect
import _bootstrap

import math

from constants import AU
from channel import Atmosphere, Turbulence
from coding import scppm_flux_needed_for_rate
from link_budget import LinkBudget
from missions import SurreyOdk16, DsocDownlink
from receiver import Receiver
from transmitter import Transmitter

# settings to change before running
# "down" is spacecraft to ground, "up" is ground to spacecraft
DIRECTION = "up"
# where the spacecraft is: "moon", "mars_close", "mars_far" or a distance in km
WHERE = "mars_close"
# data rate wanted, in bits per second
HOW_FAST = 10e6
# spacecraft telescope diameters to try on the downlink, in metres
TX_TELESCOPES_M = [0.01, 0.05, 0.10]
# ground transmit aperture diameters to try on the uplink, in metres
UPLINK_TX_TELESCOPES_M = [0.05, 0.15, 0.40]
# receiving telescope diameter on the spacecraft for the uplink, in metres
SPACECRAFT_TELESCOPE_M = 0.10
# elevation angle of the link at the ground station, in degrees
ELEVATION_DEG = 40.0
# scintillation fading is budgeted for a 1% outage, same as station_scenarios.py
OUTAGE_PROBABILITY = 0.01

# detector and seeing values come from the assumed values block of SurreyOdk16 in missions.py
_station = SurreyOdk16
# seeing angle in microradians, how much the atmosphere spreads the beam
SEEING_URAD = _station.seeing_urad
# fraction of arriving photons the detector actually counts
DETECTOR_EFF = _station.detector_eff
# false counts per second the detector gives with no light on it
DARK_COUNTS_HZ = _station.detector_dark_counts_hz


# distance in metres for each named target
DISTANCES = {
    "moon": 400e6,
    "mars_close": 0.37 * AU,
    "mars_far": 2.68 * AU,
}


# how far the coded link sits from capacity, DSOC's measured 1 dB code gap plus 0.8 dB
# lost in flight, used as the typical figure for any SCPPM style receiver
# Reference: Biswas & Srinivasan (2026), Figures 14 and 15
GAP_DB = DsocDownlink.total_gap_db


# works out the laser power in W needed to reach rate_bps, or None if it can't be reached
def best_config_power(flux_per_watt, rate_bps):
    # detected photon flux needed, picking the cheapest PPM signalling setup
    flux = scppm_flux_needed_for_rate(rate_bps, DARK_COUNTS_HZ, GAP_DB)
    # divides by the flux each watt delivers to get watts, None if no setup reaches the rate
    return None if flux is None else flux / flux_per_watt


# downlink case, returns the spacecraft laser power and beam width, or None
# same link model as the closure figure in thesis_figures.py, so both Moon answers agree
def needed_power_down(distance_m, rate_bps, station, tx_aperture_m):
    # DSOC wavelength for the spacecraft laser
    wavelength = DsocDownlink.wavelength_m
    # 1 W spacecraft transmitter with DSOC's optics loss, pointing error and the chosen telescope
    # Reference: Biswas & Srinivasan (2026), Table III
    tx = Transmitter(power_w=1.0, wavelength_m=wavelength,
                     optics_eff=DsocDownlink.tx_optics_eff,
                     pointing_error_rad=DsocDownlink.tx_pointing_error_rad,
                     aperture_diameter_m=tx_aperture_m)
    # link down to the ground station with clear sky loss and the circular aperture
    # pointing loss (equation register 4.3)
    link = LinkBudget(tx, station.build(), distance_m,
                      atmosphere=Atmosphere(
                          station_alt_m=station.site_alt_m),
                      elevation_rad=math.radians(ELEVATION_DEG),
                      pointing_model="airy")
    # extra loss from scintillation fading at the chosen outage probability
    fade = Turbulence(station_alt_m=station.site_alt_m).get_fading_loss(
        math.radians(90.0 - ELEVATION_DEG), wavelength, station.aperture_m,
        OUTAGE_PROBABILITY)
    # detected photons per second for every watt of laser power
    flux_per_watt = link.get_photon_flux() * fade * DETECTOR_EFF
    # laser power needed for the target rate
    p = best_config_power(flux_per_watt, rate_bps)
    # gives back the power and beam width if the rate can be reached
    return (p, tx.get_divergence()) if p else None


# uplink case, returns the ground laser power and the beam width after seeing, or None
def needed_power_up(distance_m, rate_bps, station, ground_aperture_m):
    # beam width from diffraction alone for this ground aperture
    diffraction = Transmitter(power_w=1.0, wavelength_m=1550e-9,
                              aperture_diameter_m=ground_aperture_m
                              ).get_divergence()
    # seeing widens the beam, diffraction and seeing add in quadrature (GOPEX model,
    # equation register 5.9), a model checked against the 1992 Galileo measurements
    effective = math.hypot(diffraction, SEEING_URAD * 1e-6)
    # ground transmitter using that width as a full 1/e^2 width (register 2.3)
    # seeing is assumed to be quoted the same way, so the gain is 32/theta^2 (register 2.5)
    # and not Liang's 16/theta^2
    # Reference: Liang et al. (2022), text after Equation 1
    tx = Transmitter(power_w=1.0, wavelength_m=1550e-9,
                     optics_eff=station.optics_eff,
                     pointing_error_rad=1e-6,
                     full_divergence_rad=effective,
                     divergence_convention="1/e2")
    # spacecraft receiver with the chosen telescope and no pointing error
    # optics efficiency of 0.63 (-2 dB) is assumed, taken the same as the LLCD transmit optics
    rx = Receiver(aperture_diameter_m=SPACECRAFT_TELESCOPE_M,
                  optics_eff=0.63,
                  pointing_error_rad=0.0)
    # link up through the atmosphere to the spacecraft
    link = LinkBudget(tx, rx, distance_m,
                      atmosphere=Atmosphere(
                          station_alt_m=station.site_alt_m),
                      elevation_rad=math.radians(ELEVATION_DEG))
    # detected photons per second for every watt, there's no fading term on the uplink
    flux_per_watt = link.get_photon_flux() * DETECTOR_EFF
    # laser power needed for the target rate
    p = best_config_power(flux_per_watt, rate_bps)
    # gives back the power and the widened beam if the rate can be reached
    return (p, effective) if p else None


# only runs when the file is started directly
if __name__ == "__main__":
    # uses the named distance if WHERE is one of the keys in DISTANCES
    if WHERE in DISTANCES:
        distance = DISTANCES[WHERE]
    # a plain number means kilometres, so it gets turned into metres
    elif isinstance(WHERE, (int, float)):
        distance = WHERE * 1e3
    # otherwise WHERE is text holding a number of km
    else:
        distance = float(WHERE) * 1e3
    # the Surrey station is used for both directions
    station = SurreyOdk16

    # downlink, prints the spacecraft laser power for each telescope size
    if DIRECTION == "down":
        print(f"DOWNLINK: to send {HOW_FAST/1e6:g} Mb/s from {WHERE} "
              f"({distance/1e3:,.0f} km) to {station.name}, "
              f"the spacecraft needs:\n")
        # loops through each spacecraft telescope size
        for ap in TX_TELESCOPES_M:
            # works out the power and beam width for this telescope
            r = needed_power_down(distance, HOW_FAST, station, ap)
            # if nothing reaches the rate, says so and moves on to the next size
            if r is None:
                print(f"  {ap*100:>4.0f} cm telescope:  cannot reach this rate")
                continue
            p, div = r
            # shows the power in mW when it's under 1 W, otherwise in W
            power = f"{p*1e3:.2f} mW" if p < 1 else f"{p:.2f} W"
            print(f"  {ap*100:>4.0f} cm telescope needs  {power:>10} of laser"
                  f" power   (beam width {div*1e6:.0f} microradians)")
    # uplink, prints the ground laser power for each aperture size
    else:
        print(f"UPLINK: to send {HOW_FAST/1e6:g} Mb/s from {station.name} "
              f"to a {SPACECRAFT_TELESCOPE_M*100:.0f} cm spacecraft telescope"
              f" at {WHERE} ({distance/1e3:,.0f} km), the ground needs:\n")
        # loops through each ground aperture size
        for ap in UPLINK_TX_TELESCOPES_M:
            # works out the power and widened beam for this aperture
            r = needed_power_up(distance, HOW_FAST, station, ap)
            # if nothing reaches the rate, says so and moves on to the next size
            if r is None:
                print(f"  {ap*100:>4.0f} cm aperture:  cannot reach this rate")
                continue
            p, div = r
            # shows the power in mW when it's under 1 W, otherwise in W
            power = f"{p*1e3:.2f} mW" if p < 1 else f"{p:.2f} W"
            print(f"  {ap*100:>4.0f} cm aperture  ->  {power:>10} of laser"
                  f" power   (beam after seeing: {div*1e6:.0f} microradians)")
        # points out that the beam stops shrinking with bigger apertures, seeing sets the floor
        print(f"\n  (seeing assumed {SEEING_URAD:.0f} microradians -- note"
              f" how the beam width stops\n   shrinking as the aperture"
              f" grows: the atmosphere sets the floor)")


# Academic Citations & Literature References:
#
# DSOC Code Gap, Transmit Optics & Pointing:
# Biswas, A., & Srinivasan, M. (2026). Deep Space Optical Communications (DSOC) System
# Description and Performance. IEEE Journal of Selected Topics in Quantum Electronics,
# 32(1). (Specifically Figures 14-15 for the 1 dB code gap and 0.8 dB flight impairment,
# and Table III for the transmit optics and pointing bound).
#
# Divergence-Based Transmit Gain:
# Liang, J., Chaudhry, A. U., Erdogan, E., & Yanikomeroglu, H. (2022). Link Budget
# Analysis for Free-Space Optical Satellite Networks. Proc. IEEE 23rd Int. Symp. on a
# World of Wireless, Mobile and Multimedia Networks (WoWMoM). (Specifically the
# 16/theta^2 gain given after Equation 1, which the uplink here deliberately does not use).
