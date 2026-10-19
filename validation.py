import math

from constants import AU
from transmitter import Transmitter
from receiver import Receiver
from channel import (Atmosphere, background_light_power, hufnagel_valley_cn2,
                     Turbulence)
from coding import (SCPPM_GUARD_FRACTION, scppm_data_rate,
                    scppm_max_data_rate, scppm_decode_threshold_ks)
from geometry import slant_distance
from link_budget import LinkBudget
from missions import DsocDownlink, LlcdDownlink, dsoc_config_for_rate
from modulation import ppm_capacity_bits_per_slot
from pointing import optimum_gain
from utility import convert_from_db, convert_to_db, watts_to_dbm

# every check gets saved in this list so run() can print them and other scripts can read them
CHECKS = []


# saves one check with its kind, source, computed and published values and the tolerance
def check(name, kind, source, computed, expected, tol, unit=""):
    # adds it to CHECKS, and it passes if the computed value is within tol of the published one
    CHECKS.append({
        "name": name,
        "kind": kind,
        "source": source,
        "computed": computed,
        "expected": expected,
        "tol": tol,
        "unit": unit,
        "ok": abs(computed - expected) <= tol,
    })


# checks the transmit antenna efficiency g_T for a Gaussian beam with no central blockage
# Reference: Klein & Degnan (1974), Section II.B
def klein_degnan_efficiency():
    # alpha is the truncation ratio (aperture over beam size), gamma is the obscuration ratio
    alpha, gamma = 1.12, 0.0
    # works out g_T from alpha and gamma
    g_t = (2.0 / alpha ** 2) * (math.exp(-alpha ** 2)
                                - math.exp(-alpha ** 2 * gamma ** 2)) ** 2
    # compares with the 0.8145 given in the paper
    check("g_T at alpha=1.12, gamma=0", "reproduction",
          "Klein & Degnan 1974, Sec. II.B", g_t, 0.8145, 5e-4)


# checks the inter-satellite link path loss and needed transmit power at three distances
# Reference: Liang et al. (2022), Table 2
def liang_isl_table2():
    # sets up the paper's 1 W, 1550 nm transmitter with 1 urad pointing error and 15 urad beam
    tx = Transmitter(
        power_w=1.0,
        wavelength_m=1550e-9,
        optics_eff=0.8,
        pointing_error_rad=1e-6,
        full_divergence_rad=15e-6,
    )
    # sets up the 8 cm receiver, with its -35.5 dBm sensitivity turned into watts
    rx = Receiver(
        aperture_diameter_m=0.080,
        optics_eff=0.8,
        pointing_error_rad=1e-6,
        sensitivity_w=convert_from_db(-35.5) * 1e-3,
    )
    # published path loss [dB], needed power [dBm] and loss tolerance at 1000, 5000, 10000 km
    # the 5000 km row gets a looser tolerance because the table prints -272.18 dB there, but
    # the exact sum is -272.16 dB and the paper's own 29.30 dBm only fits that, so it's a slip
    expected = {1000e3: (-258.18, 15.32, 0.01),
                5000e3: (-272.18, 29.30, 0.03),
                10000e3: (-278.18, 35.32, 0.01)}
    # loops through each distance and builds a link budget with no atmosphere
    for d, (lps_db, pt_dbm, lps_tol) in expected.items():
        budget = LinkBudget(tx, rx, d)
        # compares the free-space path loss in dB with the table
        check(f"ISL free-space loss at {d/1e3:.0f} km", "reproduction",
              "Liang 2022, Table 2", convert_to_db(budget.get_path_loss()), lps_db,
              lps_tol, "dB")
        # works out the transmit power needed for a 3 dB margin and compares it in dBm
        pt = budget.get_required_transmit_power(margin_linear=convert_from_db(3.0))
        check(f"ISL required P_T at {d/1e3:.0f} km", "reproduction",
              "Liang 2022, Table 2", watts_to_dbm(pt), pt_dbm, 0.02, "dBm")


# checks slant distance, Mie loss, cloud loss and needed power at three elevation angles
# Reference: Liang et al. (2022), Table 5
def liang_updown_table5():
    # atmosphere with a 1 km high station, a 20 km layer and the paper's cloud values
    atmosphere = Atmosphere(
        station_alt_m=1e3,
        layer_height_m=20e3,
        liquid_water_g_m3=3.128e-4,
        cloud_concentration_cm3=0.5,
        particle_size_coefficient=1.6,
    )
    # same 1 W, 1550 nm transmitter as the Table 2 check
    tx = Transmitter(
        power_w=1.0,
        wavelength_m=1550e-9,
        optics_eff=0.8,
        pointing_error_rad=1e-6,
        full_divergence_rad=15e-6,
    )
    # same 8 cm receiver with -35.5 dBm sensitivity
    rx = Receiver(
        aperture_diameter_m=0.080,
        optics_eff=0.8,
        pointing_error_rad=1e-6,
        sensitivity_w=convert_from_db(-35.5) * 1e-3,
    )
    # the paper uses an Earth radius of 6378.1 km, so the same value is used here
    liang_earth_radius = 6378.1e3
    # published slant distance [km], Mie loss [dB], geometric loss [dB] and needed power [dBm]
    expected = {10.0: (1813.4, -0.57, -1.22, 22.28),
                40.0: (810.7, -0.15, -0.33, 13.98),
                90.0: (549.0, -0.10, -0.21, 10.42)}
    # loops through each elevation angle
    for el_deg, (dgs_km, im_db, ig_db, pt_dbm) in expected.items():
        el = math.radians(el_deg)
        # line of sight distance from the 1 km high station up to the 550 km satellite
        d = slant_distance(el, 550e3, 1e3, earth_radius_m=liang_earth_radius)
        # compares the slant distance in km
        check(f"Slant distance at {el_deg:.0f} deg", "reproduction",
              "Liang 2022, Table 5", d / 1e3, dgs_km, 0.1, "km")
        # compares the Mie scattering loss (small particles in the air) in dB
        check(f"Mie loss at {el_deg:.0f} deg", "reproduction",
              "Liang 2022, Table 5",
              convert_to_db(atmosphere.get_mie_loss(el, 1550e-9)), im_db, 0.005, "dB")
        # compares the geometric scattering loss (bigger fog and cloud droplets) in dB
        check(f"Geometric loss at {el_deg:.0f} deg", "reproduction",
              "Liang 2022, Table 5",
              convert_to_db(atmosphere.get_geometric_loss(el, 1550e-9)), ig_db,
              0.005, "dB")
        # builds the link budget through the atmosphere at this elevation
        budget = LinkBudget(tx, rx, d, atmosphere=atmosphere, elevation_rad=el)
        # works out the power needed for a 3 dB margin and compares it in dBm
        pt = budget.get_required_transmit_power(margin_linear=convert_from_db(3.0))
        check(f"Up/downlink required P_T at {el_deg:.0f} deg", "reproduction",
              "Liang 2022, Table 5", watts_to_dbm(pt), pt_dbm, 0.02, "dBm")


# checks the Mars telemetry link budget at 2.68 AU line by line
# Reference: Valentini et al. (2021), Table III
def valentini_mars_table3():
    # the paper's 129 dB antenna gain turned into a plain ratio
    gain = convert_from_db(129.0)
    # 5 W laser at 1064 nm with -5 dB optics, and a full divergence of 4/sqrt(G) to give 129 dB
    tx = Transmitter(
        power_w=5.0,
        wavelength_m=1064e-9,
        optics_eff=convert_from_db(-5.0),
        pointing_error_rad=0.35e-6,
        full_divergence_rad=4.0 / math.sqrt(gain),
    )
    # receiver aperture picked so its gain (pi D / lambda)^2 is also 129 dB
    rx = Receiver(
        aperture_diameter_m=math.sqrt(gain) * 1064e-9 / math.pi,
        optics_eff=convert_from_db(-5.0),
        pointing_error_rad=0.35e-6,
    )
    # link budget at 2.68 AU with 4 dB of other losses
    budget = LinkBudget(tx, rx, 2.68 * AU,
                        other_losses=convert_from_db(-4.0))
    # compares the transmit antenna gain in dB
    check("Mars TM transmit gain", "reproduction",
          "Valentini 2021, Table III", convert_to_db(tx.get_gain()), 129.0, 0.01, "dB")
    # total pointing loss is the transmit and receive pointing losses multiplied together
    total_pointing = (budget.get_transmit_pointing_loss()
                      * budget.get_receive_pointing_loss())
    # compares the total pointing loss in dB
    check("Mars TM total pointing loss", "reproduction",
          "Valentini 2021, Table III", convert_to_db(total_pointing), -8.45, 0.01, "dB")
    # compares the space loss (the beam spreading out over the distance) in dB
    check("Mars TM space loss", "reproduction",
          "Valentini 2021, Table III", convert_to_db(budget.get_path_loss()), -373.49,
          0.05, "dB")
    # compares the received power in dBW
    check("Mars TM received power", "reproduction",
          "Valentini 2021, Table III", convert_to_db(budget.get_received_power()),
          -130.96, 0.05, "dBW")
    # received photons per nanosecond, in dB
    flux_db_ns = convert_to_db(budget.get_photon_flux() / 1e9)
    # compares the received photon flux
    check("Mars TM received photon flux", "reproduction",
          "Valentini 2021, Table III", flux_db_ns, -33.67, 0.05, "dB(phe/ns)")
    # margin is the received flux above the -35.76 dB(phe/ns) the receiver needs
    margin_db = flux_db_ns - (-35.76)
    # compares the link margin in dB
    check("Mars TM link margin", "reproduction",
          "Valentini 2021, Table III", margin_db, 2.09, 0.05, "dB")
    # checks the best antenna gain for a 0.35 urad pointing error also comes out near 129 dB
    # Reference: Valentini et al. (2021), Equation 17
    check("Optimum gain at theta_max=0.35 urad", "reproduction",
          "Valentini 2021, Eq. (17)", convert_to_db(optimum_gain(0.35e-6)), 129.0,
          0.15, "dB")


# checks the background light power reaching the detector for four worked cases
# Reference: Giggenbach (2022), Table 2
def giggenbach_background_table2():
    # each row is sky radiance [W/m^2/nm/sr], aperture area [m^2], field of view [rad],
    # filter bandwidth [nm] and the published power [W]
    rows = [(25e-3, 0.1, 200e-6, 40.0, 50e-9),
            (25e-3, 0.7, 100e-6, 40.0, 88e-9),
            (1.2e-3, 0.1, 200e-6, 40.0, 2.4e-9),
            (1.2e-3, 0.7, 100e-6, 40.0, 4.2e-9)]
    # loops through each case and works out the background power in watts
    for radiance, area, fov, bw, expected in rows:
        p = background_light_power(radiance, area, fov, bw)
        # compares in nW, allowing 1% of the value plus 0.35 nW
        check(f"P_BGL {expected*1e9:.1f} nW case", "reproduction",
              "Giggenbach 2022, Table 2", p * 1e9, expected * 1e9,
              0.01 * expected * 1e9 + 0.35, "nW")


# checks the HV5/7 turbulence strength C_n^2 at ground level equals the profile's ground terms
# Reference: Casanovas Ventura (2021), Equation 47
def hufnagel_valley_ground_value():
    # both sides are in units of 1e-14 m^(-2/3) so the tolerance is easy to read
    check("HV5/7 C_n^2 at h=0", "reproduction", "Casanovas 2021, Eq. (47)",
          hufnagel_valley_cn2(0.0) * 1e14, (2.7e-16 + 1.7e-14) * 1e14,
          1e-4, "1e-14")


# checks the SCPPM data rate for 4, 64 and 256-PPM with a rate 1/3 code and 256 ns slots
# Reference: Valentini et al. (2021), Tables I and II
def valentini_scppm_rates():
    # each row is the PPM order and the published data rate in bits per second
    rows = [(4, 517.32e3), (64, 97.00e3), (256, 32.33e3)]
    # loops through each PPM order and works out the data rate
    for order, expected in rows:
        rate = scppm_data_rate(1.0 / 3.0, order, 256e-9)
        # compares in kbps with a 1% tolerance
        check(f"SCPPM data rate M={order}", "reproduction",
              "Valentini 2021, Tables I-II", rate / 1e3, expected / 1e3,
              0.01 * expected / 1e3, "kbps")


# checks the Fried parameter r_0 (roughly how wide the wavefront stays undisturbed by
# turbulence) for the HV5/7 profile, where r_0 = 5 cm at 0.5 um at zenith by definition
# Reference: Andrews & Phillips (2005), Chapter 12
def hv57_fried_parameter():
    # this check also caught a misprint in the profile's ground term, where the printed
    # 1000 m scale height only gives 1.44 cm (see section 5.4 of the equation register)
    # Reference: Casanovas Ventura (2021), Equation 47
    turb = Turbulence(station_alt_m=0.0)
    # Fried parameter looking straight up at 0.5 um, turned from metres into cm
    r0_cm = turb.get_fried_parameter(0.0, 0.5e-6) * 100.0
    # compares with 5 cm, allowing 0.1 cm
    check("HV5/7 Fried parameter at 0.5 um, zenith", "reproduction",
          "definitional (Fried 1966; Andrews & Phillips Ch. 12)",
          r0_cm, 5.0, 0.1, "cm")


# compares the turbulence model with measured data from three real ground stations
# each tolerance is the measurement's own spread, so it tests consistency and is not a fit
def atmosphere_measured():
    # turbulence at the DLR ground station in Oberpfaffenhofen, 600 m above sea level
    ogs = Turbulence(station_alt_m=600.0)
    # loops through 10 and 20 degrees elevation with the measured median fade and spread in dB
    for el_deg, median_db, spread_db in ((10.0, 4.2, 1.8), (20.0, 1.9, 0.9)):
        # fading loss for a 40 cm telescope at 1550 nm and 1% outage, made positive in dB
        # the zenith angle is 90 degrees minus the elevation
        fade_db = -convert_to_db(ogs.get_fading_loss(
            math.radians(90.0 - el_deg), 1550e-9, 0.40, 0.01))
        # compares with the median fade measured over 15 OSIRISv1 satellite passes
        # Reference: Moll et al. (2022)
        check(f"Fading loss at 1% outage, 40 cm, {el_deg:.0f} deg", "measured",
              "Moll et al. 2022, SPIE 12266 (15 passes, OSIRISv1)",
              fade_db, median_db, spread_db, "dB")

    # aperture averaging factor (how much a big telescope smooths out the twinkling)
    # at 5 and 30 degrees elevation, for 847 nm and a 40 cm telescope
    avg5 = ogs.get_aperture_averaging(math.radians(85.0), 847e-9, 0.40)
    avg30 = ogs.get_aperture_averaging(math.radians(60.0), 847e-9, 0.40)
    # KIODO measured the factor falling from about 0.2 at 5 deg to about 0.1 at 30 deg, so it
    # passes if the factor is bigger near the horizon (a cos misprint would flip this round)
    # Reference: Perlot et al. (2007)
    check("Aperture averaging weakens toward horizon", "measured",
          "Perlot et al. 2007, SPIE 6457 (KIODO, 847 nm, 40 cm)",
          1.0 if avg5 > avg30 else 0.0, 1.0, 0.1)

    # turbulence at the Table Mountain Facility, 2286 m above sea level
    tmf = Turbulence(station_alt_m=2286.0)
    # compares the Rytov variance (how strong the scintillation is) at zenith and 1550 nm
    # r_0 isn't checked here, as HV5/7 starts from local ground so r_0 is the same at any site
    # Reference: Birch et al. (2024), Table 1
    check("Rytov variance at 1550 nm, Table Mountain", "measured",
          "Birch et al. 2024, Opt. Express 32, Table 1 (zenith-corrected)",
          tmf.get_rytov_variance(0.0, 1550e-9), 0.063, 0.02)


# checks the model's decode threshold for DSOC's 20.6 Mbps downlink against the flight value
# Reference: Biswas & Srinivasan (2026), Table III and Figure 16
def dsoc_decode_threshold():
    # gets the PPM order, code rate and slot time DSOC used at 20.6 Mbps
    order, code_rate, slot = dsoc_config_for_rate(20.6e6)
    # background counts per slot, scaled up from the per nanosecond value to this slot length
    kb = DsocDownlink.background_per_ns_slot * (slot / 1e-9)
    # gap to capacity is 1 dB for the code plus 0.8 dB of flight impairments
    # Reference: Biswas & Srinivasan (2026), Figures 14 and 15
    gap = DsocDownlink.code_gap_db + DsocDownlink.flight_impairment_db
    # signal counts per symbol needed to decode, using only detected counts so every
    # unpublished efficiency cancels out
    ks_model = scppm_decode_threshold_ks(order, code_rate, kb, gap)
    # compares with the design 0.93 counts per symbol, which matched flight to tenths of a dB
    check("DSOC 20.6 Mbps decode threshold Ks", "measured",
          "Biswas & Srinivasan 2026, Table III / Fig. 16",
          convert_to_db(ks_model), convert_to_db(0.93), 2.0, "dB(counts)")


# checks DSOC's received power at the Table III design point and the counts seen in flight
def dsoc_link_budget_table3():
    # builds DSOC's link at 2.29 AU with 4 W after conjunction and the table's receive chain
    link = DsocDownlink.build(2.29 * AU, DsocDownlink.power_post_conjunction_w,
                              DsocDownlink.rx_chain_eff)
    # optics and receive chain come from the table, so this tests the rest: gain from the
    # 22 cm aperture, circular aperture pointing, space loss, 5 m receive gain, -0.5 dB air
    # Reference: Biswas & Srinivasan (2026), Table III
    check("DSOC link budget, received power at Table III point",
          "reproduction", "Biswas & Srinivasan 2026, Table III (design column)",
          convert_to_db(link.get_received_power()), -118.29, 0.5, "dBW")
    # gets the PPM order and slot time for 20.6 Mbps (the code rate isn't needed here)
    order, _, slot = dsoc_config_for_rate(20.6e6)
    # one symbol lasts 32 slots of 2 ns plus the 25% guard time
    symbol_s = order * slot * (1.0 + SCPPM_GUARD_FRACTION)
    # detected counts per symbol against the 0.93 the paper says matched flight, with the
    # printed 3-sigma of 0.29 as the tolerance
    # Reference: Biswas & Srinivasan (2026), Section III.F and Table III
    check("DSOC signal counts per symbol in flight, Table III point",
          "measured",
          "Biswas & Srinivasan 2026, Sec. III.F / Table III (Ks 0.93, 3-sigma 0.29)",
          link.get_photon_flux() * symbol_s, 0.93, 0.29, "counts")


# checks DSOC's beam width against the 13 to 15 urad measured in flight before conjunction
# Reference: Biswas & Srinivasan (2026), Section III.E and Figure 11
def dsoc_flight_divergence():
    # DSOC's transmitter from its aperture and the truncation ratio measured before launch
    tx = Transmitter(power_w=1.0, wavelength_m=DsocDownlink.wavelength_m,
                     aperture_diameter_m=DsocDownlink.tx_aperture_m,
                     truncation_ratio=DsocDownlink.tx_truncation_ratio)
    # full 1/e^2 beam width that matches this gain, sqrt(32/G), turned into urad
    width_urad = math.sqrt(32.0 / tx.get_gain()) * 1e6
    # compares with the middle of the measured range, and since truncation was measured on
    # the ground but the width in flight, the two numbers are independent
    check("DSOC beam divergence measured in flight", "measured",
          "Biswas & Srinivasan 2026, Sec. III.E, Fig. 11: 13-15 urad",
          width_urad, 14.0, 1.0, "urad")


# works out the best rate DSOC's own PPM orders and slot times allow at this distance
def _dsoc_best_rate(distance_m, power_w, rx_chain_eff):
    # builds DSOC's link budget for this distance, power and receive chain efficiency
    budget = DsocDownlink.build(distance_m, power_w, rx_chain_eff)
    # flux is already detected photons so detection efficiency is 1, and only the rate is kept
    return scppm_max_data_rate(
        budget.get_photon_flux(),
        DsocDownlink.background_per_ns_slot * 1e9,
        1.0, DsocDownlink.min_slot_time_s,
        capacity_gap_db=DsocDownlink.total_gap_db,
        ppm_orders=DsocDownlink.ppm_orders,
        slot_times_s=DsocDownlink.slot_times_s)[0]


# compares the model's rate ratio between 0.56 AU and 2.68 AU with the rates DSOC flew
# Reference: Biswas & Srinivasan (2026), Sections I and III.F, Figures 13 and 17a
def dsoc_rate_vs_distance():
    # 267 Mbps at 0.37 AU is left out as it is DSOC's top rate, so it can't test the rate model
    # measured ratio of 133 Mbps at 0.56 AU over 8.33 Mbps at 2.68 AU (both at 2 W), in dB
    measured_ratio_db = convert_to_db(DsocDownlink.rate_at_mid_range_bps
                                      / DsocDownlink.rate_at_max_range_bps)
    # low and high ends of the receive chain efficiency, from the link budget's tolerance
    lo, hi = DsocDownlink.rx_chain_eff_range
    # loops through low, nominal and high efficiency, since efficiency moves both rates
    # together and the ratio is the fair thing to compare
    for eff in (lo, DsocDownlink.rx_chain_eff, hi):
        # best rate at 0.56 AU with 2 W before conjunction
        near = _dsoc_best_rate(DsocDownlink.mid_range_au * AU,
                               DsocDownlink.power_pre_conjunction_w, eff)
        # best rate at 2.68 AU with the same power
        far = _dsoc_best_rate(DsocDownlink.max_range_au * AU,
                              DsocDownlink.power_pre_conjunction_w, eff)
        # compares the model's ratio with the flown one in dB
        check(f"DSOC rate ratio 0.56/2.68 AU (eff={eff})", "measured",
              "Biswas & Srinivasan 2026, Sec. III.F / Sec. I (133 vs 8.33 Mbps, 2 W)",
              convert_to_db(near / far), measured_ratio_db, 2.5, "dB")


# rough check that LLCD's 622 Mbps Moon downlink closes and works at about one bit per photon
def llcd_photon_budget():
    # the optics efficiencies are assumed, so both tolerances are kept wide on purpose
    # builds the LLCD downlink budget
    budget = LlcdDownlink.build()
    # detected photon flux over all four 40 cm telescopes with their 50% efficient SNSPDs
    flux_detected = (budget.get_photon_flux() * LlcdDownlink.rx_telescopes
                     * LlcdDownlink.snspd_eff)
    # 16-PPM with 5 GHz slots
    m = LlcdDownlink.ppm_order
    slot = LlcdDownlink.slot_time_s
    # signal counts per symbol, with no guard time at LLCD's top rate
    ks = flux_detected * m * slot
    # counts per symbol needed to decode the rate 1/2 code
    ks_th = scppm_decode_threshold_ks(m, LlcdDownlink.code_rate,
                                      LlcdDownlink.background_per_slot,
                                      LlcdDownlink.code_gap_db)
    # margin is how many times more counts arrive than are needed, in dB
    margin_db = convert_to_db(ks / ks_th)
    # passes if the link closes with a margin somewhere between 0 and 15 dB
    # Reference: Boroson et al. (2009), Section III
    check("LLCD 622 Mbps link closes (margin in 0..15 dB)", "measured",
          "Boroson 2009, Sec. III", 1.0 if 0.0 < margin_db < 15.0 else 0.0,
          1.0, 0.1)
    # bits per detected photon at the threshold, code rate times log2(M) over the counts needed
    bits_per_photon_at_threshold = (LlcdDownlink.code_rate * math.log2(m)
                                    / ks_th)
    # compares with the roughly one bit per photon stated for LLCD, within 4 dB
    check("LLCD photon efficiency at threshold", "measured",
          "Boroson 2009: 'approximately one bit per photon'",
          convert_to_db(bits_per_photon_at_threshold),
          convert_to_db(LlcdDownlink.photon_efficiency_bpp), 4.0, "dB(b/ph)")


# checks the coherent detection results: BPSK photons per bit, DPSK linewidth limit and BER
def kikuchi_coherent_sensitivity():
    from coherent import bpsk_homodyne_ber, bpsk_heterodyne_ber
    from scipy.optimize import brentq
    # finds the photons per bit where homodyne BPSK hits a BER of 1e-5, which should be 4.55
    # Reference: Kikuchi (2016), Equation 60
    ns_homo = brentq(lambda n: bpsk_homodyne_ber(n) - 1e-5, 0.1, 50.0)
    check("BPSK homodyne photons/bit at 1e-5", "reproduction",
          "Kikuchi 2016, after Eq. (60): 4.55", ns_homo, 4.55, 0.01)
    # same for heterodyne BPSK, which should need about 9 photons per bit
    # Reference: Kikuchi (2016), Equation 74
    ns_het = brentq(lambda n: bpsk_heterodyne_ber(n) - 1e-5, 0.1, 50.0)
    check("BPSK heterodyne photons/bit at 1e-5", "reproduction",
          "Kikuchi 2016, after Eq. (74): 9.0", ns_het, 9.0, 0.1)
    # with 10 MHz transmit and local oscillator lasers, heterodyne DPSK keeps the phase error
    # within 0.2 rad only above 3 Gbit/s
    # Reference: Kikuchi et al. (1984), Section V-B
    from coherent import max_linewidth_for_dpsk
    # allowed linewidth goes up in step with bit rate, so the rate for the 20 MHz pair is
    # 20 MHz over the linewidth allowed at 1 bit/s
    # the tolerance is 0.5 Gbit/s because the paper only gives one significant figure
    check("DPSK bit rate needed by a 10 MHz laser pair", "reproduction",
          "Kikuchi 1984, Sec. V-B: (1/T) > 3 Gbit/s",
          2.0 * 10e6 / max_linewidth_for_dpsk(1.0) / 1e9, 3.0, 0.5,
          "Gbit/s")
    from coherent import dpsk_ber_with_phase_noise
    # finds the S/N where heterodyne DPSK with a 0.20 rad phase error reaches a BER of 1e-9
    # Reference: Kikuchi et al. (1984), Equation 36
    gamma_1e9 = brentq(lambda g: math.log10(
        dpsk_ber_with_phase_noise(g, 0.2)) + 9.0, 10.0, 1e3)
    # compares in dB with the 16.1 dB read off the figure, which is only good to about 0.5 dB
    # Reference: Kikuchi et al. (1984), Figure 3
    check("DPSK S/N for BER 1e-9 at sigma_phi 0.20 rad", "reproduction",
          "Kikuchi 1984, Fig. 3 (digitised)",
          10.0 * math.log10(gamma_1e9), 16.1, 0.5, "dB")


# checks the M-PPM capacity per slot with Kb = 1 against points read off a published figure
# Reference: Hemmati et al. (2011), Figure 5 and Equation 10
def hemmati_capacity_fig5():
    # each point is PPM order, capacity level and the Ks/M in dB where that curve's rising edge
    # crosses the level, read off seven curves to about 0.1 dB (the flat tops test nothing)
    points = [(2, 0.1, -3.14), (2, 0.01, -8.99), (4, 0.01, -11.39),
              (8, 0.01, -13.19), (16, 0.003, -17.52), (32, 0.003, -19.06),
              (64, 0.003, -20.56), (256, 0.003, -23.26)]
    # biggest miss in dB over all eight points
    worst = 0.0
    # loops through each point
    for order, level, read_db in points:
        # 50000 Monte Carlo samples up to 64-PPM and 20000 for 256-PPM
        samples = 50_000 if order <= 64 else 20_000
        # search range for Ks/M in dB
        lo, hi = -35.0, 10.0
        # bisection, halving the range 25 times to find where the capacity reaches the level
        for _ in range(25):
            mid = 0.5 * (lo + hi)
            # capacity at this Ks/M, with signal counts M times Ks/M and 1 background count
            cap = ppm_capacity_bits_per_slot(order * convert_from_db(mid), 1.0,
                                             order, samples=samples)
            # if the capacity is already high enough the crossing is lower, otherwise higher
            lo, hi = (lo, mid) if cap >= level else (mid, hi)
        # keeps the largest gap between the crossing found and the value read off the figure
        worst = max(worst, abs(0.5 * (lo + hi) - read_db))
    # 0.3 dB tolerance covers the reading error plus about 0.15 dB of Monte Carlo scatter
    check("PPM capacity at Kb=1, worst Ks/M miss over 8 points",
          "reproduction", "Hemmati 2011, Fig. 5 (digitised)", worst, 0.0, 0.3,
          "dB")


# checks GOPEX, the laser uplink to Galileo, night by night compared with the first night
# Reference: Levine et al. (1993), Table 2
def gopex_uplink_scaling():
    from missions import GopexUplink
    # everything is relative to night 1, so the unpublished camera constants all cancel
    # the model goes as 1/(effective beam width squared x distance squared) with clear-air loss
    for day in (2, 3, 6, 7, 8):
        # model change in dB against the paper's own predicted change, with a tight 0.6 dB
        check(f"GOPEX day {day} vs paper theory", "reproduction",
              "Levine et al. TDA 42-114, Table 2 (expected photoelectrons)",
              GopexUplink.model_delta_db(day),
              GopexUplink.paper_expected_delta_db(day), 0.6, "dB")
        # model change against the measured mean photoelectrons per pulse from the fitted
        # lognormal, with a wide 2.5 dB because the measurements scatter from night to night
        # Reference: Levine et al. (1993), Table 2 and Equation 6
        check(f"GOPEX day {day} vs flight measurement", "measured",
              "Levine et al. TDA 42-114, Table 2 (mean from m_l, sigma_N; Eq. 6)",
              GopexUplink.model_delta_db(day),
              GopexUplink.measured_delta_db(day), 2.5, "dB")


# checks the NODE terminal parts list adds up to the published mass and power totals
# Reference: Kingsbury (2015), Tables 3.5 and 3.6
def node_terminal_totals():
    from terminal import total_average_power_w, total_mass_g
    # total mass should come to 550 g
    check("NODE terminal mass total", "reproduction",
          "Kingsbury 2015, Table 3.5", total_mass_g(), 550.0, 0.5, "g")
    # average power should come to 8.52 W, allowing 0.02 W because the printed total is rounded
    check("NODE terminal average power total", "reproduction",
          "Kingsbury 2015, Table 3.6 (8.52 printed, 0.01 rounding)",
          total_average_power_w(), 8.52, 0.02, "W")


# runs every check, prints a table of results and returns True if they all pass
def run():
    # runs each group of checks, which fills up CHECKS
    klein_degnan_efficiency()
    liang_isl_table2()
    liang_updown_table5()
    valentini_mars_table3()
    giggenbach_background_table2()
    hufnagel_valley_ground_value()
    hv57_fried_parameter()
    valentini_scppm_rates()
    atmosphere_measured()
    dsoc_decode_threshold()
    dsoc_link_budget_table3()
    dsoc_flight_divergence()
    dsoc_rate_vs_distance()
    llcd_photon_budget()
    kikuchi_coherent_sensitivity()
    gopex_uplink_scaling()
    node_terminal_totals()
    hemmati_capacity_fig5()

    # length of the longest check name, so the table columns line up
    width = max(len(c["name"]) for c in CHECKS)
    # counts how many checks pass
    n_pass = 0
    # prints the table heading
    print(f"{'CHECK':<{width}}  {'KIND':<12} {'COMPUTED':>12} {'EXPECTED':>12}"
          f"  {'UNIT':<10} RESULT")
    # loops through each check and prints its row
    for c in CHECKS:
        # PASS or FAIL depending on the tolerance test in check()
        status = "PASS" if c["ok"] else "FAIL"
        # True counts as 1, so this adds one for every pass
        n_pass += c["ok"]
        # prints the row, and for a failure also shows which source to look at
        print(f"{c['name']:<{width}}  {c['kind']:<12} {c['computed']:>12.4f}"
              f" {c['expected']:>12.4f}  {c['unit']:<10} {status}"
              + ("" if c["ok"] else f"   <-- {c['source']}"))
    # prints how many passed out of the total
    print(f"\n{n_pass}/{len(CHECKS)} checks passed")
    return n_pass == len(CHECKS)


# when run directly, exits with 0 if every check passed and 1 if any failed
if __name__ == "__main__":
    raise SystemExit(0 if run() else 1)


# Academic Citations & Literature References:
#
# Transmit Antenna Efficiency:
# Klein, B. J., & Degnan, J. J. (1974). Optical Antenna Gain. 1: Transmitting Antennas.
# Applied Optics, 13(9), 2134-2141.
# (Specifically Section II.B for the transmit efficiency at alpha = 1.12).
#
# Link Budget, Slant Distance & Atmospheric Loss:
# Liang, J., Chaudhry, A. U., Erdogan, E., & Yanikomeroglu, H. (2022). Link Budget
# Analysis for Free-Space Optical Satellite Networks. Proc. IEEE 23rd Int. Symp. on a
# World of Wireless, Mobile and Multimedia Networks (WoWMoM).
# (Specifically Table 2 for the inter-satellite link and Table 5 for the up/downlink).
#
# Deep Space Link Budget & SCPPM Data Rates:
# Valentini, L., Faedi, A., Paolini, E., & Chiani, M. (2021). Analysis of Pointing Loss
# Effects in Deep Space Optical Links. Proc. IEEE Global Communications Conf. (GLOBECOM).
# (Specifically Table III for the Mars link, Tables I-II for the SCPPM rates and
# Equation 17 for the optimum gain).
#
# Background Light:
# Giggenbach, D. (2022). Free-Space Optical Data Receivers with Avalanche Detectors for
# Satellite Downlinks Regarding Background Light. Sensors, 22(18), 6773.
# (Specifically Table 2 for the four worked background power values).
#
# Turbulence Profile & Fried Parameter:
# Casanovas Ventura, M. (2021). An Assessment on the Requirements for Deep Space Optical
# Communications. Master's thesis, Universitat Politecnica de Catalunya.
# (Specifically Equation 47 for the Hufnagel-Valley ground term).
# Andrews, L. C., & Phillips, R. L. (2005). Laser Beam Propagation through Random Media
# (2nd ed.). SPIE Press.
# (Specifically Chapter 12 for the Fried parameter).
#
# Measured Turbulence at Real Ground Stations:
# Moll, F., Giggenbach, D., Schmidt, C., & Fuchs, C. (2022). Analysis of power scintillation
# and fading margin in the LEO-ground downlink with the OSIRISv1 laser terminal on Flying
# Laptop and the DLR optical ground station Oberpfaffenhofen. Proc. SPIE 12266, 122660F.
# (Specifically the fading loss measured over 15 satellite passes).
# Perlot, N., et al. (2007). Results of the optical downlink experiment KIODO from OICETS
# satellite to optical ground station Oberpfaffenhofen (OGS-OP). Proc. SPIE 6457, 645704.
# (Specifically the aperture averaging measured at 847 nm with a 40 cm telescope).
# Birch, M., Piazzolla, S., Hooser, P., Bennet, F., Travouillon, T., & Buehlman, W. (2024).
# Optical turbulence profiling at the Table Mountain Facility with the Laser Communication
# Relay Demonstration GEO downlink. Optics Express, 32(12), 21962-21976.
# (Specifically Table 1 for the Rytov variance, corrected to zenith).
#
# DSOC Flight Mission:
# Biswas, A., & Srinivasan, M. (2026). Deep Space Optical Communications (DSOC) System
# Description and Performance. IEEE Journal of Selected Topics in Quantum Electronics, 32(1).
# (Specifically Table III, Sections I, III.E and III.F, and Figures 11 and 13 to 17).
#
# LLCD Moon Downlink:
# Boroson, D. M., Scozzafava, J. J., Murphy, D. V., Robinson, B. S., & Shaw, H. (2009). The
# Lunar Laser Communications Demonstration (LLCD). Proc. Third IEEE Int. Conf. on Space
# Mission Challenges for Information Technology (SMC-IT), 23-28.
# https://doi.org/10.1109/SMC-IT.2009.57 (Specifically Section III for the 622 Mbps
# downlink and its approximately one bit per photon efficiency).
#
# Coherent Detection:
# Kikuchi, K. (2016). Fundamentals of Coherent Optical Fiber Communications. Journal of
# Lightwave Technology, 34(1), 157-179.
# (Specifically the values after Equations 60 and 74 for BPSK photons per bit).
# Kikuchi, K., Okoshi, T., Nagamatsu, M., & Henmi, N. (1984). Degradation of Bit-Error Rate
# in Coherent Optical Communications Due to Spectral Spread of the Transmitter and the Local
# Oscillator. Journal of Lightwave Technology, 2(6), 1024-1033.
# (Specifically Section V-B, Equation 36 and Figure 3 for DPSK with laser phase noise).
#
# PPM Channel Capacity:
# Hemmati, H., Biswas, A., & Djordjevic, I. B. (2011). Deep-Space Optical Communications:
# Future Perspectives and Applications. Proceedings of the IEEE, 99(11), 2020-2039.
# (Specifically Figure 5 and Equation 10 for the M-PPM capacity per slot).
#
# GOPEX Uplink to Galileo:
# Levine, B. M., Shaik, K. S., & Yan, T.-Y. (1993). Data Analysis for GOPEX Image Frames.
# TDA Progress Report 42-114, JPL.
# (Specifically Table 2 and Equation 6 for the expected and measured photoelectrons).
#
# Terminal Mass & Power:
# Kingsbury, R. W. (2015). Optical Communications for Small Satellites. PhD thesis,
# Massachusetts Institute of Technology.
# (Specifically Tables 3.5 and 3.6 for the NODE terminal mass and power totals).
