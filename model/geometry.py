import math
from constants import EARTH_RADIUS


# works out the slant range in metres from the station to the satellite
# station altitude is 0 unless the user gives one, and the Earth radius can be
# swapped because different missions and papers use slightly different values
def slant_distance(elevation_rad, sat_alt_m, station_alt_m=0, earth_radius_m=EARTH_RADIUS):
    # distance from the Earth centre up to the station
    r = earth_radius_m + station_alt_m
    # height of the satellite above the station rather than above sea level
    h = sat_alt_m - station_alt_m
    # square root part of the slant range formula
    # Reference: Liang et al. (2022), Equation 4
    root = math.sqrt(((r + h)/ r)**2-math.cos(elevation_rad)**2)
    # slant range from the station to the satellite
    slant_range_m =r*(root-math.sin(elevation_rad))
    return slant_range_m

# gets the length in metres of the path through the atmosphere layer at this elevation
def atmos_path_length(elevation_rad, layer_height_m, station_alt_m=0.0):
    # layer thickness above the station stretched out by the slant angle
    # Reference: Liang et al. (2022), text under Equation 10
    atmos_path_length_m = (layer_height_m - station_alt_m)/math.sin(elevation_rad)
    return atmos_path_length_m


# Academic Citations & Literature References:
#
# Slant Range & Atmospheric Path Length:
# Liang, J., Chaudhry, A. U., Erdogan, E., & Yanikomeroglu, H. (2022). Link Budget
# Analysis for Free-Space Optical Satellite Networks. Proc. IEEE 23rd Int. Symp. on a
# World of Wireless, Mobile and Multimedia Networks (WoWMoM). (Specifically Equation 4
# for the slant range and the text under Equation 10 for the path length).
