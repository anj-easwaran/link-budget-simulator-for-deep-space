import math

from transmitter import Transmitter
from receiver import Receiver
from link_budget import LinkBudget
from utility import convert_from_db, convert_to_db


# Psyche DSOC downlink at 1550 nm down to the Palomar Hale telescope
# every value comes from Biswas and Srinivasan 2026 unless another source is given
class DsocDownlink:

    # laser wavelength in m
    # Reference: Biswas & Srinivasan (2026), Abstract
    wavelength_m = 1550e-9
    # laser power in W before solar conjunction
    # Reference: Biswas & Srinivasan (2026), Section III.E
    power_pre_conjunction_w = 2.0
    # laser power in W after conjunction
    # Reference: Biswas & Srinivasan (2026), Abstract and Section III.E
    power_post_conjunction_w = 4.0
    # flight laser transceiver aperture in m, 22 cm
    # Reference: Biswas & Srinivasan (2026), Abstract
    tx_aperture_m = 0.22
    # design full 1/e^2 divergence in rad, the measured values were 13.5 to 15.1 urad
    # Reference: Biswas & Srinivasan (2026), Section III.E
    tx_divergence_rad = 13.8e-6
    # measured design truncation ratio, then the favourable (1.18) and adverse (1.51) values,
    # which give divergences of 13.8, 13.5 and 15.1 urad
    # Reference: Biswas & Srinivasan (2026), Section III.E
    tx_truncation_ratio = 1.23
    tx_truncation_range = (1.18, 1.51)
    # transmit pointing error bound in rad, the sub-urad bias and 0.5 urad 1 sigma jitter added up
    # Reference: Biswas & Srinivasan (2026), Section III.E
    tx_pointing_error_rad = 1.0e-6
    # Hale telescope receive aperture in m, 5 m
    # Reference: Biswas & Srinivasan (2026), Abstract
    rx_aperture_m = 5.0
    # background counts per 1 ns slot, the range was 0.0003 to 0.006 (Figure 16 note)
    # Reference: Biswas & Srinivasan (2026), Figure 14 and Figure 16
    background_per_ns_slot = 0.003
    # shortest slot in s, 267 Mbps flew with M = 16, rate 2/3 and 0.5 ns slots
    # Reference: Biswas & Srinivasan (2026)
    min_slot_time_s = 0.5e-9
    # signalling set the rate search can pick from. the four setups in dsoc_config_for_rate()
    # flew (M = 16, 32, 64 and slots 0.5, 1, 2 ns), the rest are ASSUMED to fill out the modem.
    # the limit matters since a free search at far range picks slot widths no flight modem
    # has and overstates the rate
    # PPM orders, ASSUMED past the demonstrated 16, 32 and 64
    ppm_orders = (16, 32, 64, 128)
    # slot widths in s, ASSUMED past the demonstrated 0.5, 1 and 2 ns
    slot_times_s = (0.5e-9, 1.0e-9, 2.0e-9, 4.0e-9, 8.0e-9)
    # furthest range in AU
    # Reference: Biswas & Srinivasan (2026), Abstract
    max_range_au = 2.68
    # middle range in AU, where the 133 Mbps decoding test was done
    # Reference: Biswas & Srinivasan (2026), Section III.F
    mid_range_au = 0.56
    # closest range in AU
    # Reference: Biswas & Srinivasan (2026), Section I
    min_range_au = 0.37
    # measured rate at the closest range in bps at 2 W, the highest rate reached
    # Reference: Biswas & Srinivasan (2026), Abstract
    rate_at_min_range_bps = 267e6
    # measured rate at the middle range in bps at 2 W
    # Reference: Biswas & Srinivasan (2026), Section III.F and Figure 14
    rate_at_mid_range_bps = 133e6
    # measured rate at the furthest range in bps at 2 W before conjunction, since
    # Figures 13 and 17(a) show no 4 W pass past 2.4 AU
    # Reference: Biswas & Srinivasan (2026), Section I and Figures 13 and 17(a)
    rate_at_max_range_bps = 8.33e6
    # flight transceiver optical transmission, -1.72 dB design, then the adverse (-2.34 dB)
    # and favourable (-1.42 dB) range
    # Reference: Biswas & Srinivasan (2026), Table III
    tx_optics_eff = 0.673
    tx_optics_eff_range = (0.58, 0.72)
    # receive chain behind the Hale aperture from the ground losses. optical transmission
    # -2.60 dB, laser, detector and receiver jitter -0.8 dB, detection efficiency -1.43 dB,
    # and the Hale obscuration -0.45 dB (net gain 139.67 dB against 140.12 dB for the
    # unblocked 5 m aperture used here). total -5.28 dB design, -4.18 favourable, -5.96 adverse
    # Reference: Biswas & Srinivasan (2026), Table III
    rx_chain_eff = 0.30
    rx_chain_eff_range = (0.25, 0.38)
    # code gap to capacity in dB, about 1 dB from the 133 Mbps discussion
    # (the sentence itself points to Figure 15)
    # Reference: Biswas & Srinivasan (2026), Section III.F and Figure 14
    code_gap_db = 1.0
    # extra loss in dB of the flight decoder against an ideal one at 133 Mbps, same sentence
    # Reference: Biswas & Srinivasan (2026), Section III.F
    flight_impairment_db = 0.8
    # total gap in dB that the rate search takes off the signal
    total_gap_db = code_gap_db + flight_impairment_db

    # builds the DSOC link budget for a distance in m, laser power in W and receive chain efficiency
    # rx_chain_eff covers everything behind the Hale aperture (optics, obscuration, jitter and
    # SNSPD detection), DsocDownlink.rx_chain_eff is the design value and its range the tolerance
    @staticmethod
    def build(distance_m, power_w, rx_chain_eff):
        # flight transmitter with the measured truncation ratio
        tx = Transmitter(
            power_w=power_w,
            wavelength_m=DsocDownlink.wavelength_m,
            optics_eff=DsocDownlink.tx_optics_eff,
            pointing_error_rad=DsocDownlink.tx_pointing_error_rad,
            aperture_diameter_m=DsocDownlink.tx_aperture_m,
            truncation_ratio=DsocDownlink.tx_truncation_ratio,
        )
        # Hale receiver, with the receive chain efficiency put in as its optics efficiency
        rx = Receiver(
            aperture_diameter_m=DsocDownlink.rx_aperture_m,
            optics_eff=rx_chain_eff,
            pointing_error_rad=0.0,
        )
        # gain comes from the aperture (docs/equations.md 2.1) like DSOC's own budget does,
        # since 16/theta^2 on the 1/e^2 divergence gives too little (2.4), and the circular
        # aperture pointing loss (4.3) matches the budget table
        # Reference: Biswas & Srinivasan (2026), Table III
        # the -0.5 dB other loss is an ASSUMED night time clear air allowance, rounded from the
        # -0.44 dB design attenuation in Table III (mean -0.12 dB). DSOC budgets use MODTRAN
        # (Section III.A), which is not published as a model, and the Liang polynomial (5.1)
        # is not used because Palomar at 1712 m is outside its fitted range
        return LinkBudget(tx, rx, distance_m,
                          other_losses=convert_from_db(-0.5),
                          pointing_model="airy")


# LLCD downlink at 1550 nm from LADEE down to the LLGT array of four 40 cm telescopes
# every value comes from Boroson et al. 2009 unless another source is given
class LlcdDownlink:

    # laser wavelength in m, the band is 1550 to 1570 nm
    # Reference: Boroson et al. (2009), Section III
    wavelength_m = 1550e-9
    # laser power in W from the 0.5 W MOPA
    # Reference: Boroson et al. (2009), Section III
    power_w = 0.5
    # transmit reflector aperture in m, 10 cm
    # Reference: Boroson et al. (2009), Section III
    tx_aperture_m = 0.10
    # transmit pointing error bound in rad, ASSUMED since the MIRU stabilised pointing is
    # not published in this paper
    tx_pointing_error_rad = 1.0e-6
    # receive aperture of one telescope in m, and the number of telescopes in the array
    # Reference: Boroson et al. (2009), Section III
    rx_aperture_m = 0.40
    rx_telescopes = 4
    # SNSPD detection efficiency of 50%
    # Reference: Boroson et al. (2009), Section III
    snspd_eff = 0.50
    # PPM order, 16-ary PPM
    # Reference: Boroson et al. (2009), Section III
    ppm_order = 16
    # slot time in s from the 5 GHz slot clock
    # Reference: Boroson et al. (2009), Section III
    slot_time_s = 0.2e-9
    # code rate, ASSUMED 1/2 from MLCD heritage. Section I gives MLCD a rate 1/2 near capacity
    # code, but for LLCD Section III only says the code nearly reaches Shannon capacity, no rate
    code_rate = 0.5
    # highest data rate in bps
    # Reference: Boroson et al. (2009), Abstract
    max_rate_bps = 622e6
    # link range in m, about 400,000 km
    # Reference: Boroson et al. (2009), Abstract
    range_m = 400e3 * 1e3
    # photon efficiency, about one bit per photon
    # Reference: Boroson et al. (2009), Section III
    photon_efficiency_bpp = 1.0
    # receive optics efficiency, ASSUMED -3 dB for the telescope and multimode fibre coupling
    rx_optics_eff = 0.50
    # night sky and dark counts per 0.2 ns slot, ASSUMED because it is not published
    background_per_slot = 1e-4
    # code gap in dB, ASSUMED because the turbo code gap plus implementation loss is not
    # published, for comparison SCPPM class codes are about 1 dB
    # Reference: Hemmati et al. (2011), Section III
    code_gap_db = 2.5

    # builds the LLCD link budget at its set range
    @staticmethod
    def build():
        # flight transmitter, the 0.63 optics efficiency is ASSUMED (-2 dB) as it is not published
        tx = Transmitter(
            power_w=LlcdDownlink.power_w,
            wavelength_m=LlcdDownlink.wavelength_m,
            optics_eff=0.63,
            pointing_error_rad=LlcdDownlink.tx_pointing_error_rad,
            aperture_diameter_m=LlcdDownlink.tx_aperture_m,
        )
        # one telescope only, multiply the flux by rx_telescopes to get the whole array
        rx = Receiver(
            aperture_diameter_m=LlcdDownlink.rx_aperture_m,
            optics_eff=LlcdDownlink.rx_optics_eff,
            pointing_error_rad=0.0,
        )
        # the -1 dB other loss is the atmosphere, ASSUMED to be clear sky
        return LinkBudget(tx, rx, LlcdDownlink.range_m,
                          other_losses=convert_from_db(-1.0))


# GOPEX (December 1992), a pulsed 532 nm uplink from Table Mountain to the Galileo SSI camera
# at 0.7 to 6.0 million km, the only measured received signal dataset over many ranges in
# Literature/. the camera aperture and QE are not published, so everything is compared to
# day 1 where those constants cancel out
class GopexUplink:

    # laser wavelength in m
    # Reference: Wilson & Lesh (1993), Table 1
    wavelength_m = 532e-9
    # pulse energy in J at the Table Mountain Facility
    # Reference: Wilson & Lesh (1993), Table 1
    pulse_energy_j = 0.250
    # clear sky zenith transmittance at 532 nm, read off the plot to +/-0.01 (3 km site, 23 km
    # visual range) and ASSUMED to hold at TMF (2.2 km) on the GOPEX nights (equations.md 5.9)
    # Reference: Biswas & Piazzolla (2006), Figure 3-18
    zenith_transmittance = 0.825
    # one tuple per day: range in km, divergence in urad, seeing in urad, zenith angle in deg,
    # pulses received, lognormal mean m_l and std sigma_N from the data, and the paper's
    # expected photoelectrons per pulse
    # Reference: Levine et al. (1993), Table 2
    days = {
        1: (0.7e6, 110.0, 20.0, 50.0, 33, 10.86, 0.7263, 67116.0),
        2: (1.45e6, 110.0, 7.5, 51.0, 59, 9.27, 0.5356, 14860.0),
        3: (2.21e6, 110.0, 20.0, 53.0, 41, 8.45, 0.940, 7267.0),
        6: (4.48e6, 60.0, 40.0, 55.0, 24, 7.70, 1.047, 4352.0),
        7: (5.25e6, 60.0, 20.0, 56.0, 60, 8.36, 0.8060, 3295.0),
        8: (6.01e6, 60.0, 30.0, 56.0, 51, 7.41, 0.8619, 2629.0),
    }

    # clear air transmittance along the slant path on a day, T_zenith^sec(zenith) (5.9)
    @staticmethod
    def _transmittance(day):
        # that night's zenith angle in radians
        zenith = math.radians(GopexUplink.days[day][3])
        return GopexUplink.zenith_transmittance ** (1.0 / math.cos(zenith))

    # model predicted received signal on a day relative to day 1 in dB, 1/(Theta_eff^2 d^2)
    # with the clear air transmittance at that night's average zenith angle (equations.md 5.9)
    @staticmethod
    def model_delta_db(day):
        # range, divergence and seeing on day 1 and on the chosen day
        d1, div1, see1 = GopexUplink.days[1][:3]
        d, div, see = GopexUplink.days[day][:3]
        # beam spread ratio, with divergence and seeing added in quadrature
        spread = math.hypot(div1, see1) / math.hypot(div, see)
        # clear air transmittance ratio between the two nights
        air = GopexUplink._transmittance(day) / GopexUplink._transmittance(1)
        # range ratio squared times spread ratio squared times the air ratio, in dB
        return convert_to_db((d1 / d) ** 2 * spread ** 2 * air)

    # measured mean photoelectrons per pulse relative to day 1 in dB, using
    # <N> = exp(m_l + sigma_N^2 / 2) since that is what a link budget predicts (5.11)
    @staticmethod
    def measured_delta_db(day):
        # lognormal mean and std on day 1 and on the chosen day
        m1, s1 = GopexUplink.days[1][5:7]
        m, s = GopexUplink.days[day][5:7]
        return convert_to_db(math.exp((m - m1) + (s * s - s1 * s1) / 2.0))

    # plus or minus one standard error of measured_delta_db in dB. each night's mean of n pulses
    # has relative error sqrt(exp(sigma_N^2)-1)/sqrt(n) and the two nights add in quadrature
    @staticmethod
    def measured_se_db(day):
        # day 1 is the reference so it has no error
        if day == 1:
            return 0.0

        # relative error of one night's mean signal
        def relative_error(d):
            # number of pulses and lognormal std for that night
            n, _, s = GopexUplink.days[d][4:7]
            return math.sqrt(math.expm1(s * s) / n)

        # combines both nights in quadrature and turns it into dB
        return convert_to_db(math.exp(math.hypot(relative_error(day),
                                                 relative_error(1))))

    # the paper's own Kiasaleh and Yan prediction relative to day 1 in dB
    # Reference: Levine et al. (1993), Table 2
    @staticmethod
    def paper_expected_delta_db(day):
        return convert_to_db(GopexUplink.days[day][7] / GopexUplink.days[1][7])


# University of Surrey ground station with the Orion Optics ODK16 telescope
# telescope values are from the Orion Optics and Teleskop-Express ODK16 spec sheets
class SurreyOdk16:

    # station name for the printed tables
    name = "Surrey ODK16"
    # telescope aperture in m, 400 mm
    # Reference: Orion Optics UK (2026)
    aperture_m = 0.400
    # central obscuration ratio, the 152 mm secondary over the 400 mm primary
    # Reference: Orion Optics UK (2026)
    obscuration_ratio = 0.38
    # site height in m and latitude in deg for Guildford, which is close to sea level
    site_alt_m = 0.0
    site_latitude_deg = 51.24

    # every value below is ASSUMED, a placeholder until it is measured in the lab
    # changing one here updates every script (scenarios, inverse solver, closure maps)
    # telescope and relay optics throughput, ASSUMED since it has not been measured
    optics_eff = 0.7
    # detection efficiency, ASSUMED 0.5 as an SNSPD class placeholder since no comms detector
    # is confirmed at Surrey
    detector_eff = 0.5
    # detector false counts per second, ASSUMED. the counter is also ASSUMED to be an array
    # with no dead time blocking at the worked out count rates. one SNSPD or GM-APD pixel
    # saturates well below the GHz detected rates of the near Earth (Moon and LEO) scenarios,
    # which is why DSOC uses arrays, and deep space rates far below 1 photon per slot are fine
    detector_dark_counts_hz = 1000.0
    # detector full field of view in urad, ASSUMED
    detector_fov_urad = 100.0
    # optical background filter width in nm, ASSUMED
    filter_bandwidth_nm = 1.0
    # site seeing in urad, ASSUMED (the GOPEX days went from 7.5 to 40)
    seeing_urad = 15.0
    # the mount model and max tracking rate are unknown and not in the model yet, and the
    # stabilised laser arriving in 2026 has unknown wavelength, power and linewidth
    # (pessimistic, optimistic) range for each ASSUMED value above (equations.md section 10)
    # each covers what the unknown could really be and pushes rate or margin one way only,
    # so running both ends together brackets the answer. the ranges are ASSUMED bounds too
    assumed_ranges = {
        # uncoated relay optics up to clean optics
        "optics_eff": (0.5, 0.85),
        # InGaAs GM-APD up to a good SNSPD
        "detector_eff": (0.2, 0.85),
        "detector_dark_counts_hz": (1e4, 1e2),
        # a wide field of view lets in more background light
        "detector_fov_urad": (200.0, 50.0),
        "filter_bandwidth_nm": (2.0, 0.2),
        # the GOPEX day to day spread
        "seeing_urad": (40.0, 7.5),
    }

    # builds the ODK16 as a Receiver, with an optional sensitivity in W
    @staticmethod
    def build(sensitivity_w=None):
        return Receiver(
            aperture_diameter_m=SurreyOdk16.aperture_m,
            obscuration_ratio=SurreyOdk16.obscuration_ratio,
            optics_eff=SurreyOdk16.optics_eff,
            pointing_error_rad=0.0,
            sensitivity_w=sensitivity_w,
        )


# returns (PPM order, code rate, slot time in s) for a DSOC data rate that flew
# these setups are decoded from the SCPPM frame maths with the guard time as (1 + guard),
# and (32, 1/3, 2 ns) is the only combination whose frame rate lands on 20.7 Mbps
# Reference: Valentini et al. (2021), Equation 20
# Reference: Biswas & Srinivasan (2026)
def dsoc_config_for_rate(rate_bps):
    # data rate in bps matched to its (order, code rate, slot time) setup
    table = {
        267e6: (16, 2.0 / 3.0, 0.5e-9),
        133e6: (16, 2.0 / 3.0, 1.0e-9),
        25e6: (64, 1.0 / 3.0, 1.0e-9),
        20.6e6: (32, 1.0 / 3.0, 2.0e-9),
    }
    return table[rate_bps]


# Academic Citations & Literature References:
#
# DSOC Downlink:
# Biswas, A., & Srinivasan, M. (2026). Deep Space Optical Communications (DSOC) System
# Description and Performance. IEEE Journal of Selected Topics in Quantum Electronics, 32(1).
# (Specifically the abstract, Sections I, III.A, III.E and III.F, Table III and Figures 13,
# 14, 16 and 17(a)).
# Valentini, L., Faedi, A., Paolini, E., & Chiani, M. (2021). Analysis of Pointing Loss
# Effects in Deep Space Optical Links. Proc. IEEE Global Communications Conf. (GLOBECOM).
# (Specifically Equation 20, the SCPPM frame rate used to decode the DSOC configurations).
#
# LLCD Downlink:
# Boroson, D. M., Scozzafava, J. J., Murphy, D. V., Robinson, B. S., & Shaw, H. (2009). The
# Lunar Laser Communications Demonstration (LLCD). Proc. Third IEEE Int. Conf. on Space
# Mission Challenges for Information Technology (SMC-IT), 23-28.
# https://doi.org/10.1109/SMC-IT.2009.57 (Specifically the abstract and Sections I and III).
# Hemmati, H., Biswas, A., & Djordjevic, I. B. (2011). Deep-Space Optical Communications:
# Future Perspectives and Applications. Proceedings of the IEEE, 99(11), 2020-2039.
# (Specifically Section III for the gap to capacity of SCPPM class codes).
#
# GOPEX Uplink:
# Wilson, K. E., & Lesh, J. R. (1993). An Overview of the Galileo Optical Experiment (GOPEX).
# TDA Progress Report, 42-114. (Specifically Table 1).
# Levine, B. M., Shaik, K. S., & Yan, T.-Y. (1993). Data Analysis for GOPEX Image Frames
# (TDA Progress Report 42-114). JPL. (Specifically Table 2, the TMF columns).
# Biswas, A., & Piazzolla, S. (2006). The Atmospheric Channel. In H. Hemmati (Ed.), Deep
# Space Optical Communications (JPL DESCANSO Series, Vol. 7, Chapter 3). Wiley.
# (Specifically Figure 3-18 for the clear sky zenith transmittance at 532 nm).
#
# Surrey ODK16 Ground Station:
# Orion Optics UK. (2026). ODK16: 16-inch corrected Dall-Kirkham telescope [Manufacturer
# product page]. https://www.orionoptics.co.uk/product/odk16/ (Specifically the aperture
# and secondary mirror size).
