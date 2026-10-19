# puts model/ on the import path, it is only imported for that side effect
import _bootstrap

import math
from contextlib import contextmanager

from constants import AU, PLANCK, LIGHT_SPEED
from channel import (Atmosphere, Turbulence, background_light_power,
                     SKY_RADIANCE_1550)
from coding import scppm_max_data_rate
from geometry import slant_distance
from link_budget import LinkBudget
from missions import SurreyOdk16, DsocDownlink
from transmitter import Transmitter
from utility import convert_from_db, convert_to_db, watts_to_dbm


# settings for every scenario
# elevation angle of the link at the station, in degrees
ELEVATION_DEG = 40.0
# scintillation fading is budgeted for a 1% outage
OUTAGE_PROBABILITY = 0.01
# daytime sky for the second rate, blue sky at 30 deg elevation, the table entry closest
# to the 40 deg operating point (equation register 5.8)
# it goes through the station's assumed filter width and field of view, night is the normal case
# Reference: Giggenbach (2022), Tables 1 and 2
DAY_SKY = "blue sky 30 deg elevation"


# the Surrey ODK16 telescope, all of its values come from missions.py
station = SurreyOdk16
# elevation angle turned into radians
ELEVATION = math.radians(ELEVATION_DEG)
# light collecting area in m^2, minus the part blocked by the central obstruction
APERTURE_AREA = (math.pi * (station.aperture_m / 2) ** 2
                 * (1.0 - station.obscuration_ratio ** 2))


# makes a clear sky atmosphere model at the station's altitude
def atmosphere():
    return Atmosphere(station_alt_m=station.site_alt_m)


# works out the background photons per second reaching the detector for a named sky
def background_flux(sky):
    # background light power in W from the sky radiance, telescope area, view angle and filter
    p = background_light_power(SKY_RADIANCE_1550[sky], APERTURE_AREA,
                               station.detector_fov_urad * 1e-6,
                               station.filter_bandwidth_nm)
    # divides by the energy of one 1550 nm photon to get photons per second
    return p / (PLANCK * LIGHT_SPEED / 1550e-9)


# scintillation fading loss for the set outage probability at this elevation
def fading_loss(elevation_rad):
    # turbulence model at the station's altitude
    turb = Turbulence(station_alt_m=station.site_alt_m)
    # zenith angle is 90 deg minus the elevation, done in radians
    return turb.get_fading_loss(math.pi / 2.0 - elevation_rad, 1550e-9,
                                station.aperture_m, OUTAGE_PROBABILITY)


# finds the best coded PPM rate the link supports under the chosen sky
# returns (rate, code rate, PPM order, slot time)
def _best_rate(link, sky="night"):
    # gap to capacity is DSOC's measured code gap plus flight impairment, as in the closure study
    # Reference: Biswas & Srinivasan (2026), Figures 14 and 15
    return scppm_max_data_rate(
        link.get_photon_flux() * fading_loss(ELEVATION),
        background_flux(sky), station.detector_eff,
        0.5e-9, station.detector_dark_counts_hz,
        capacity_gap_db=DsocDownlink.total_gap_db)


# temporarily moves every assumed station value to its pessimistic or optimistic end
# case is "pessimistic", "optimistic" or None for normal values (equation register section 10)
@contextmanager
def assumed_case(case):
    # gets the (pessimistic, optimistic) pair for each assumed value, empty if there are none
    ranges = getattr(station, "assumed_ranges", {})
    # if no case is picked or there are no ranges, runs the block with the normal values
    if case is None or not ranges:
        yield
        return
    # index 0 picks the pessimistic end and 1 the optimistic end
    idx = 0 if case == "pessimistic" else 1
    # saves the current values so they can be put back afterwards
    saved = {k: getattr(station, k) for k in ranges}
    # sets each assumed value to the chosen end, then runs the block
    try:
        for k, pair in ranges.items():
            setattr(station, k, pair[idx])
        yield
    # puts the saved values back, even if something went wrong in the block
    finally:
        for k, v in saved.items():
            setattr(station, k, v)


# works out the (pessimistic, optimistic) rates across the assumed value ranges
# build_link is a function so the link gets rebuilt with the edge values in place
def _rate_bounds(build_link, sky="night"):
    out = []
    # loops through both ends and stores the best rate for each one
    for case in ("pessimistic", "optimistic"):
        with assumed_case(case):
            out.append(_best_rate(build_link(), sky)[0])
    return tuple(out)


# LEO CubeSat at 600 km sending 10 Gb/s OOK, prints its link margin and whether it closes
def scenario_leo_ook():
    # 1 W transmitter with Liang's optics efficiency, pointing error and 15 microradian beam
    # Reference: Liang et al. (2022), Table 2
    tx = Transmitter(power_w=1.0, wavelength_m=1550e-9,
                     optics_eff=0.8, pointing_error_rad=1e-6,
                     full_divergence_rad=15e-6)
    # station receiver with Liang's -35.5 dBm sensitivity turned into watts
    # Reference: Liang et al. (2022), Table 2
    rx = station.build(sensitivity_w=convert_from_db(-35.5) * 1e-3)
    # slant distance from the station up to a 600 km orbit at this elevation
    d = slant_distance(ELEVATION, 600e3, station.site_alt_m)
    # link budget through the clear sky atmosphere
    link = LinkBudget(tx, rx, d, atmosphere=atmosphere(),
                      elevation_rad=ELEVATION)
    # link margin in dB, above zero means the link closes
    margin = convert_to_db(link.get_margin())
    # prints the margin and the verdict
    print(f"LEO 600 km CubeSat, 10 Gb/s OOK   margin {margin:+6.1f} dB  "
          f"-> {'CLOSES' if margin > 0 else 'DOES NOT CLOSE'}")


# LLCD class terminal at the Moon, prints the rate, the assumed range bounds and the day rate
def scenario_moon_ppm():
    # 0.5 W laser with a 10 cm telescope, the LLCD class values
    # Reference: Boroson et al. (2009)
    tx = Transmitter(power_w=0.5, wavelength_m=1550e-9,
                     optics_eff=0.63, pointing_error_rad=1e-6,
                     aperture_diameter_m=0.10)

    # builds the 400,000 km link to the station, kept as a function so the bounds can rebuild it
    def build_link():
        return LinkBudget(tx, station.build(), 400e6, atmosphere=atmosphere(),
                          elevation_rad=ELEVATION)

    # normal link and its best night time rate, code rate, PPM order and slot time
    link = build_link()
    rate, rc, m, slot = _best_rate(link)
    # pessimistic and optimistic rates from the assumed value ranges
    lo, hi = _rate_bounds(build_link)
    # best rate under the daytime sky instead
    day = _best_rate(link, DAY_SKY)[0]
    # received power in dBm
    p = watts_to_dbm(link.get_received_power())
    # if the link closes, prints the rate and its settings
    if rate > 0:
        print(f"Moon, LLCD-class TX, coded PPM    P_rx {p:6.1f} dBm  "
              f"-> CLOSES at {rate/1e6:7.2f} Mb/s "
              f"(R={rc:.2f}, {m}-PPM, {slot*1e9:.1f} ns; "
              f"bounds {lo/1e6:.2f}-{hi/1e6:.2f}; day {day/1e6:.2f})")
    # otherwise prints that it doesn't close, along with the bounds
    else:
        print(f"Moon, LLCD-class TX               P_rx {p:6.1f} dBm  "
              f"-> DOES NOT CLOSE at any rate "
              f"(bounds {lo/1e6:.2f}-{hi/1e6:.2f} Mb/s)")


# DSOC class terminal at distance_au, prints the rate, the bounds and the day rate
# gain from the 22 cm aperture with measured truncation (register 2.1) as DSOC's own budget
def scenario_deep_space_ppm(distance_au):
    # 4 W DSOC class flight transmitter
    tx = _dsoc_class_transmitter()

    # builds the link out to distance_au with circular aperture pointing loss (register 4.3)
    def build_link():
        return LinkBudget(tx, station.build(),
                          distance_au * AU,
                          atmosphere=atmosphere(), elevation_rad=ELEVATION,
                          pointing_model="airy")

    # normal link and its best night time rate, code rate, PPM order and slot time
    link = build_link()
    rate, rc, m, slot = _best_rate(link)
    # pessimistic and optimistic rates from the assumed value ranges
    lo, hi = _rate_bounds(build_link)
    # best rate under the daytime sky instead
    day = _best_rate(link, DAY_SKY)[0]
    # received power in dBm
    p = watts_to_dbm(link.get_received_power())
    # if the link closes, prints the rate with the day rate, or says the day link doesn't close
    if rate > 0:
        print(f"{distance_au:4.2f} AU, DSOC-class TX, PPM      "
              f"P_rx {p:6.1f} dBm  -> CLOSES at {rate/1e3:7.1f} kb/s "
              f"(R={rc:.2f}, {m}-PPM, {slot*1e9:.1f} ns; "
              f"bounds {lo/1e3:.1f}-{hi/1e3:.1f}; "
              f"day {day/1e3:.1f} kb/s)" if day > 0 else
              f"{distance_au:4.2f} AU, DSOC-class TX, PPM      "
              f"P_rx {p:6.1f} dBm  -> CLOSES at {rate/1e3:7.1f} kb/s "
              f"(R={rc:.2f}, {m}-PPM, {slot*1e9:.1f} ns; "
              f"bounds {lo/1e3:.1f}-{hi/1e3:.1f}; day: does not close)")
    # otherwise prints that it doesn't close, along with the bounds
    else:
        print(f"{distance_au:4.2f} AU, DSOC-class TX           "
              f"P_rx {p:6.1f} dBm  -> DOES NOT CLOSE at any rate "
              f"(bounds {lo/1e3:.1f}-{hi/1e3:.1f} kb/s)")


# makes the 4 W DSOC class flight transmitter with DSOC's optics, pointing, aperture and truncation
# Reference: Biswas & Srinivasan (2026)
def _dsoc_class_transmitter():
    return Transmitter(power_w=4.0, wavelength_m=DsocDownlink.wavelength_m,
                       optics_eff=DsocDownlink.tx_optics_eff,
                       pointing_error_rad=DsocDownlink.tx_pointing_error_rad,
                       aperture_diameter_m=DsocDownlink.tx_aperture_m,
                       truncation_ratio=DsocDownlink.tx_truncation_ratio)


# link from the DSOC class transmitter to the station at any distance and elevation
# not used in this file, thesis_figures.py and check_thesis_numbers.py call it
def _dsoc_class_link(distance_m, elevation_rad):
    return LinkBudget(_dsoc_class_transmitter(), station.build(), distance_m,
                      atmosphere=atmosphere(), elevation_rad=elevation_rad,
                      pointing_model="airy")


# only runs when the file is started directly
if __name__ == "__main__":
    # heading line with the telescope size, elevation and which sky counts as day
    print(f"{station.name} ({station.aperture_m*100:.0f} cm) as ground"
          f" receiver, {ELEVATION_DEG:.0f} deg elevation; rates for a clear"
          f" night, 'day' = {DAY_SKY}\n")
    # runs the LEO and Moon scenarios
    scenario_leo_ook()
    scenario_moon_ppm()
    # runs the deep space scenario at 0.1, 0.5, 1 and 2.68 AU
    for au in (0.1, 0.5, 1.0, 2.68):
        scenario_deep_space_ppm(au)


# Academic Citations & Literature References:
#
# LEO Terminal & Receiver Sensitivity:
# Liang, J., Chaudhry, A. U., Erdogan, E., & Yanikomeroglu, H. (2022). Link Budget
# Analysis for Free-Space Optical Satellite Networks. Proc. IEEE 23rd Int. Symp. on a
# World of Wireless, Mobile and Multimedia Networks (WoWMoM). (Specifically Table 2 for
# the transmit optics, pointing error, 15 microradian divergence and -35.5 dBm sensitivity).
#
# Daytime Sky Radiance:
# Giggenbach, D. (2022). Free-Space Optical Data Receivers with Avalanche Detectors for
# Satellite Downlinks Regarding Background Light. Sensors, 22(18), 6773. (Specifically
# Tables 1 and 2 for the blue-sky radiance at 1550 nm and 30 degrees elevation).
#
# DSOC Flight Terminal & Code Gap:
# Biswas, A., & Srinivasan, M. (2026). Deep Space Optical Communications (DSOC) System
# Description and Performance. IEEE Journal of Selected Topics in Quantum Electronics,
# 32(1). (Specifically the 4 W flight terminal, and Figures 14-15 for the measured code gap
# plus flight impairment).
#
# LLCD Class Moon Terminal:
# Boroson, D. M., Scozzafava, J. J., Murphy, D. V., Robinson, B. S., & Shaw, H. (2009). The
# Lunar Laser Communications Demonstration (LLCD). Proc. Third IEEE Int. Conf. on Space
# Mission Challenges for Information Technology (SMC-IT), 23-28.
# https://doi.org/10.1109/SMC-IT.2009.57 (Specifically the 0.5 W laser and 10 cm telescope).
