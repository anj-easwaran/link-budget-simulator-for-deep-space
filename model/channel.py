import math

from scipy.integrate import quad
from scipy.special import erfinv

from geometry import atmos_path_length


# free space path loss as a linear value, (lambda / 4 pi d)^2 with both lengths in m
# Reference: Liang et al. (2022), Equations 2 and 5
def free_space_loss(wavelength_m, distance_m):
    return (wavelength_m / (4.0 * math.pi * distance_m)) ** 2


# clear air Mie scattering plus Kim model fog and cloud scattering above the ground station
# only valid for stations 0 to 5 km up (the range the Mie polynomial was fitted over),
# turbulence and background light are handled separately further down
class Atmosphere:

    # stores the station altitude and top of the scattering layer in m, cloud liquid water
    # in g/m^3, cloud droplet concentration per cm^3 and the particle size coefficient
    def __init__(self,
                 station_alt_m=0.0,
                 layer_height_m=20e3,
                 liquid_water_g_m3=3.128e-4,
                 cloud_concentration_cm3=0.5,
                 particle_size_coefficient=1.6):
        self.station_alt_m = station_alt_m
        self.layer_height_m = layer_height_m
        self.liquid_water_g_m3 = liquid_water_g_m3
        self.cloud_concentration_cm3 = cloud_concentration_cm3
        self.particle_size_coefficient = particle_size_coefficient

    # linear Mie scattering loss along the slant path at this elevation
    # Reference: Liang et al. (2022), Equations 6 and 7
    def get_mie_loss(self, elevation_rad, wavelength_m):
        # wavelength in micrometres and station height in km, the units the polynomial wants
        lam = wavelength_m * 1e6
        h = self.station_alt_m / 1e3
        # polynomial coefficients, each one depends on the wavelength
        # Reference: ITU-R (2003)
        a = -0.000545 * lam ** 2 + 0.002 * lam - 0.0038
        b = 0.00628 * lam ** 2 - 0.0232 * lam + 0.00439
        c = -0.028 * lam ** 2 + 0.101 * lam - 0.18
        d = -0.228 * lam ** 3 + 0.922 * lam ** 2 - 1.26 * lam + 0.719
        # extinction ratio at the station height
        rho = a * h ** 3 + b * h ** 2 + c * h + d
        # the fitted polynomial goes negative (a gain) above about 1.5 km at 1550 nm,
        # so it is clamped at zero extinction (docs/equations.md 5.1 validity note)
        rho = max(rho, 0.0)
        # the loss grows as the path gets longer at low elevation
        return math.exp(-rho / math.sin(elevation_rad))

    # linear Kim model scattering loss over the path through the lower atmosphere
    # Reference: Liang et al. (2022), Equations 8 to 10
    def get_geometric_loss(self, elevation_rad, wavelength_m):
        # visibility in km worked out from the cloud liquid water and droplet concentration
        visibility_km = 1.002 / (self.liquid_water_g_m3
                                 * self.cloud_concentration_cm3) ** 0.6473
        # wavelength in nm
        lam_nm = wavelength_m * 1e9
        # attenuation per km, scaled from the 550 nm visibility by the particle size coefficient
        att_per_km = ((3.91 / visibility_km)
                      * (lam_nm / 550.0) ** (-self.particle_size_coefficient))
        # length of the path through the scattering layer in km
        path_km = atmos_path_length(elevation_rad, self.layer_height_m,
                                          self.station_alt_m) / 1e3
        # loss falls off exponentially with attenuation times path length
        return math.exp(-att_per_km * path_km)

    # total linear atmospheric loss, the Mie loss times the scattering loss
    # Reference: Liang et al. (2022), Equation 11
    def get_loss(self, elevation_rad, wavelength_m):
        return (self.get_mie_loss(elevation_rad, wavelength_m)
                * self.get_geometric_loss(elevation_rad, wavelength_m))


# refractive index structure constant C_n^2 in m^(-2/3) at a height, Hufnagel-Valley profile
# the default wind speed and ground strength give the standard HV5/7 profile
# Reference: Casanovas Ventura (2021), Equation 47
def hufnagel_valley_cn2(height_m, wind_speed_ms=21.0, a0=1.7e-14):
    # height above the local ground, not sea level, which is why Turbulence passes in h - h0
    # Reference: ITU-R (2015), Equation 6
    h = height_m
    # adds the wind driven high altitude term, the middle term and the ground layer term
    # the ground term uses the usual 100 m scale height. Casanovas prints 1000 m but that
    # is a misprint, only 100 m gives the HV5/7 Fried parameter r_0 of 5 cm
    # (docs/equations.md 5.4 erratum, checked in validation.py through the Fried parameter)
    # Reference: Andrews & Phillips (2005)
    # Reference: ITU-R (2015)
    return (0.00594 * (wind_speed_ms / 27.0) ** 2
            * (1e-5 * h) ** 10 * math.exp(-h / 1000.0)
            + 2.7e-16 * math.exp(-h / 1500.0)
            + a0 * math.exp(-h / 100.0))


# downlink scintillation using Andrews plane wave theory over a Hufnagel-Valley C_n^2 profile,
# with aperture averaging and a lognormal fading loss for a target outage probability
# Reference: Casanovas Ventura (2021), Equations 43 to 58
class Turbulence:

    # stores the station altitude in m, wind speed in m/s, ground C_n^2 strength
    # and the height in m where the profile integrals stop
    def __init__(self,
                 station_alt_m=0.0,
                 wind_speed_ms=21.0,
                 cn2_ground=1.7e-14,
                 profile_top_m=30e3):
        self.station_alt_m = station_alt_m
        self.wind_speed_ms = wind_speed_ms
        self.cn2_ground = cn2_ground
        self.profile_top_m = profile_top_m

    # C_n^2 at a height h above sea level, turned into height above the station first
    # the integrals use sea level heights but the profile is defined above local ground,
    # so without this the ground layer would go missing at any raised site
    # Reference: ITU-R (2015), Equation 6
    def _cn2(self, h):
        return hufnagel_valley_cn2(h - self.station_alt_m, self.wind_speed_ms,
                                   self.cn2_ground)

    # Rytov variance sigma_1^2 for a plane wave on a slant path
    # Reference: Casanovas Ventura (2021), Equations 44 to 46
    def get_rytov_variance(self, zenith_rad, wavelength_m):
        # station altitude, where the integral starts
        h0 = self.station_alt_m
        # with xi = (h-h0)/(H-h0) the (H-h0)^(5/6) factor cancels, so only this integral of
        # C_n^2 times (h-h0)^(5/6) is needed and the satellite altitude drops out. the
        # integrals die off within about 30 km, which is why this works for deep space too
        integral = quad(lambda h: self._cn2(h) * (h - h0) ** (5.0 / 6.0),
                        h0, self.profile_top_m, limit=200)[0]
        # optical wavenumber in rad/m
        k = 2.0 * math.pi / wavelength_m
        # the path gets longer away from zenith, which is where the cos term comes in
        return (2.25 * k ** (7.0 / 6.0) * integral
                / math.cos(zenith_rad) ** (11.0 / 6.0))

    # Fried parameter r_0 in m for a plane wave, [0.423 k^2 sec(zenith) integral C_n^2 dh]^(-3/5)
    # (docs/equations.md 5.10)
    # Reference: Andrews & Phillips (2005), Chapter 12
    def get_fried_parameter(self, zenith_rad, wavelength_m):
        # integral of C_n^2 from the station to the top of the profile
        mu0 = quad(self._cn2, self.station_alt_m, self.profile_top_m,
                   limit=200)[0]
        # optical wavenumber in rad/m
        k = 2.0 * math.pi / wavelength_m
        return (0.423 * k ** 2 * mu0 / math.cos(zenith_rad)) ** (-3.0 / 5.0)

    # scintillation index sigma_I^2 from the weak to strong turbulence theory
    # Reference: Casanovas Ventura (2021), Equation 43
    def get_scintillation_index(self, zenith_rad, wavelength_m):
        # Rytov variance for this zenith angle and wavelength
        s1 = self.get_rytov_variance(zenith_rad, wavelength_m)
        return math.exp(0.49 * s1 / (1.0 + 1.11 * s1 ** 1.2) ** (7.0 / 6.0)
                        + 0.51 * s1 / (1.0 + 0.69 * s1 ** 1.2) ** (5.0 / 6.0)
                        ) - 1.0

    # turbulence scale height h_s in m, a C_n^2 weighted average height of the profile
    # Reference: Casanovas Ventura (2021), Equation 49
    def get_turbulence_scale_height(self):
        h0 = self.station_alt_m
        # integral of C_n^2 times the height above the station squared
        num = quad(lambda h: self._cn2(h) * (h - h0) ** 2,
                   h0, self.profile_top_m, limit=200)[0]
        # integral of C_n^2 times the height above the station to the power 5/6
        den = quad(lambda h: self._cn2(h) * (h - h0) ** (5.0 / 6.0),
                   h0, self.profile_top_m, limit=200)[0]
        return (num / den) ** (6.0 / 7.0)

    # aperture averaging factor, how much a big receiver smooths out the scintillation
    # uses sec of the zenith angle like the source does. Casanovas Ventura Equation 48 prints
    # cos there, which is an erratum (docs/equations.md 5.6)
    # Reference: Biswas & Piazzolla (2006), Equation 3.4-18
    # Reference: Yura & McKinley (1983)
    def get_aperture_averaging(self, zenith_rad, wavelength_m,
                               rx_diameter_m):
        # scale height of the turbulence profile
        hs = self.get_turbulence_scale_height()
        # aperture size compared with the Fresnel zone. a longer slant path widens the zone,
        # so the averaging gets weaker towards the horizon
        x = rx_diameter_m ** 2 * math.cos(zenith_rad) / (hs * wavelength_m)
        return 1.0 / (1.0 + 1.11 * x ** (7.0 / 6.0))

    # power scintillation index sigma_P^2, the aperture averaging times the scintillation index
    # Reference: Casanovas Ventura (2021), Equation 51
    def get_power_scintillation_index(self, zenith_rad, wavelength_m,
                                      rx_diameter_m):
        return (self.get_aperture_averaging(zenith_rad, wavelength_m,
                                            rx_diameter_m)
                * self.get_scintillation_index(zenith_rad, wavelength_m))

    # linear fading loss to allow for so the received power stays above the threshold except
    # for outage_probability of the time, using lognormal fading
    # Reference: Casanovas Ventura (2021), Equations 52 to 58
    def get_fading_loss(self, zenith_rad, wavelength_m, rx_diameter_m,
                        outage_probability):
        # power scintillation index at the receiver
        sp2 = self.get_power_scintillation_index(zenith_rad, wavelength_m,
                                                 rx_diameter_m)
        # lognormal variance of the received power
        s_ld2 = math.log(sp2 + 1.0)
        # loss from the inverse error function at the chosen outage probability
        return (math.exp(erfinv(2.0 * outage_probability - 1.0)
                         * math.sqrt(2.0 * s_ld2))
                / math.sqrt(sp2 + 1.0))


# background sky spectral radiance at 1550 nm in W/(m^2 nm sr) for different sky conditions
# Reference: Giggenbach (2022), Table 1
SKY_RADIANCE_1550 = {
    "night": 0.0,
    "blue sky zenith": 1.2e-3,
    "blue sky 30 deg elevation": 2.0e-3,
    "blue sky horizon": 25e-3,
    "sunlit cloud": 20e-3,
    "overcast cloud": 2e-3,
    "moon disk in FOV": 20e-3,
}


# background light power in W on the detector, radiance x area x solid angle x filter width
# with the solid angle as 4 pi sin^2(omega)
# Reference: Giggenbach (2022), Equation 22
def background_light_power(radiance_w_m2_nm_sr, aperture_area_m2,
                           fov_full_flat_rad, filter_bandwidth_nm):
    # solid angle of the field of view, written as printed in the paper since that gives
    # all four worked values in its Table 2
    # Reference: Giggenbach (2022), Table 2
    omega_sr = 4.0 * math.pi * math.sin(fov_full_flat_rad) ** 2
    return (radiance_w_m2_nm_sr * aperture_area_m2 * omega_sr
            * filter_bandwidth_nm)


# Academic Citations & Literature References:
#
# Free-Space Path Loss & Atmospheric Scattering:
# Liang, J., Chaudhry, A. U., Erdogan, E., & Yanikomeroglu, H. (2022). Link Budget
# Analysis for Free-Space Optical Satellite Networks. Proc. IEEE 23rd Int. Symp. on a
# World of Wireless, Mobile and Multimedia Networks (WoWMoM). (Specifically Equations 2
# and 5 for path loss, 6 and 7 for Mie loss, 8 to 10 for Kim model scattering and 11 for
# the total atmospheric loss).
#
# Mie Scattering Polynomial:
# ITU-R. (2003). Recommendation P.1622: Prediction methods required for the design of
# Earth-space systems operating between 20 THz and 375 THz. International
# Telecommunication Union. (Specifically the extinction polynomial used in the Mie loss).
#
# Turbulence Profile & Downlink Scintillation:
# Casanovas Ventura, M. (2021). An Assessment on the Requirements for Deep Space Optical
# Communications (Master's thesis). Universitat Politècnica de Catalunya. (Specifically
# Equation 47 for the Hufnagel-Valley profile and Equations 43 to 58 for scintillation,
# aperture averaging and fading loss).
# ITU-R. (2015). Propagation data required for the design of Earth-space systems operating
# between 20 THz and 375 THz (Recommendation ITU-R P.1621-2). International
# Telecommunication Union. (Specifically Equation 6, the profile height above local ground).
# Andrews, L. C., & Phillips, R. L. (2005). Laser Beam Propagation through Random Media
# (2nd ed.). SPIE Press. (Specifically Chapter 12 for the Fried parameter, and the 100 m
# ground layer scale height).
#
# Aperture Averaging:
# Biswas, A., & Piazzolla, S. (2006). The Atmospheric Channel. In H. Hemmati (Ed.), Deep
# Space Optical Communications (JPL DESCANSO Series, Vol. 7, Chapter 3). Wiley.
# (Specifically Equation 3.4-18).
# Yura, H. T., & McKinley, W. G. (1983). Aperture averaging of scintillation for
# space-to-ground optical communication applications. Applied Optics, 22(11), 1608-1609.
# (The origin of the aperture averaging formula).
#
# Background Light:
# Giggenbach, D. (2022). Free-Space Optical Data Receivers with Avalanche Detectors for
# Satellite Downlinks Regarding Background Light. Sensors, 22(18), 6773. (Specifically
# Table 1 for sky radiance, Equation 22 for background power and Table 2 to check it).
