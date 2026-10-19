from constants import PLANCK, LIGHT_SPEED
from channel import free_space_loss
from pointing import gaussian_pointing_loss, airy_pointing_loss
from utility import convert_to_db, watts_to_dbm


# puts a transmitter, receiver and the channel together into one link
class LinkBudget:
    # stores the link parts, distance in m, optional atmosphere and elevation in rad
    # other_losses is any extra linear loss, pointing_model is "gaussian" or "airy"
    def __init__(self,
                 transmitter,
                 receiver,
                 distance_m,
                 atmosphere=None,
                 elevation_rad=None,
                 other_losses=1.0,
                 pointing_model="gaussian"):
        # if an atmosphere is given without an elevation angle, stop since the loss needs it
        if atmosphere is not None and elevation_rad is None:
            raise ValueError("an atmosphere needs an elevation angle")
        self.transmitter = transmitter
        self.receiver = receiver
        self.distance_m = distance_m
        self.atmosphere = atmosphere
        self.elevation_rad = elevation_rad
        self.other_losses = other_losses
        self.pointing_model = pointing_model

    # returns the linear pointing loss using whichever pointing model was picked
    def _pointing_loss(self, gain, error_rad):
        # if the airy model was chosen, use the circular aperture version
        if self.pointing_model == "airy":
            return airy_pointing_loss(gain, error_rad)
        # otherwise fall back to the Gaussian beam version
        return gaussian_pointing_loss(gain, error_rad)

    # linear pointing loss at the transmitter from its gain and pointing error
    def get_transmit_pointing_loss(self):
        return self._pointing_loss(self.transmitter.get_gain(),
                                   self.transmitter.pointing_error_rad)

    # linear pointing loss at the receiver, its gain depends on the laser wavelength
    def get_receive_pointing_loss(self):
        gain = self.receiver.get_gain(self.transmitter.wavelength_m)
        return self._pointing_loss(gain, self.receiver.pointing_error_rad)

    # linear free space path loss over the link distance
    def get_path_loss(self):
        return free_space_loss(self.transmitter.wavelength_m, self.distance_m)

    # linear atmospheric loss at the set elevation
    def get_atmospheric_loss(self):
        # if no atmosphere is given, return 1 so it adds no loss
        if self.atmosphere is None:
            return 1.0
        return self.atmosphere.get_loss(self.elevation_rad,
                                        self.transmitter.wavelength_m)

    # received power in W, every term of the link equation multiplied together (register 3.2)
    def get_received_power(self):
        return (self.transmitter.power_w
                * self.transmitter.optics_eff
                * self.receiver.optics_eff
                * self.transmitter.get_gain()
                * self.receiver.get_gain(self.transmitter.wavelength_m)
                * self.get_transmit_pointing_loss()
                * self.get_receive_pointing_loss()
                * self.get_atmospheric_loss()
                * self.get_path_loss()
                * self.other_losses)

    # received signal photons per second, the received power over the energy of one photon
    def get_photon_flux(self):
        # energy of one photon at the laser wavelength
        # Reference: Valentini et al. (2021), Equation 18
        photon_energy = PLANCK * LIGHT_SPEED / self.transmitter.wavelength_m
        return self.get_received_power() / photon_energy

    # link margin, how many times more power arrives than the receiver needs
    # Reference: Liang et al. (2022), Equation 12
    def get_margin(self):
        # if the receiver has no sensitivity set there is nothing to compare against, so stop
        if self.receiver.sensitivity_w is None:
            raise ValueError("receiver has no sensitivity_w set")
        return self.get_received_power() / self.receiver.sensitivity_w

    # transmit power in W needed to get the given linear margin at the receiver sensitivity
    def get_required_transmit_power(self, margin_linear=1.0):
        # if the receiver has no sensitivity set, stop with an error
        if self.receiver.sensitivity_w is None:
            raise ValueError("receiver has no sensitivity_w set")
        # received power for every watt sent, since the link scales straight with transmit power
        per_watt = self.get_received_power() / self.transmitter.power_w
        # power needed is the margin times the sensitivity divided by the power per watt
        return margin_linear * self.receiver.sensitivity_w / per_watt

    # builds the link budget table as rows of (name, value in dB or dBm, unit)
    def report(self):
        # transmitter side and path loss rows, these are always shown
        rows = [
            ("Transmit power", watts_to_dbm(self.transmitter.power_w), "dBm"),
            ("Transmit antenna gain", convert_to_db(self.transmitter.get_gain()), "dB"),
            ("Transmit optics efficiency", convert_to_db(self.transmitter.optics_eff), "dB"),
            ("Transmit pointing loss", convert_to_db(self.get_transmit_pointing_loss()), "dB"),
            ("Free-space path loss", convert_to_db(self.get_path_loss()), "dB"),
        ]
        # if there is an atmosphere, add its loss row
        if self.atmosphere is not None:
            rows.append(("Atmospheric loss", convert_to_db(self.get_atmospheric_loss()), "dB"))
        # receiver side rows
        rows += [
            ("Receive antenna gain",
             convert_to_db(self.receiver.get_gain(self.transmitter.wavelength_m)), "dB"),
            ("Receive optics efficiency", convert_to_db(self.receiver.optics_eff), "dB"),
            ("Receive pointing loss", convert_to_db(self.get_receive_pointing_loss()), "dB"),
        ]
        # only add the other losses row if some extra loss was actually set
        if self.other_losses != 1.0:
            rows.append(("Other losses", convert_to_db(self.other_losses), "dB"))
        # received power row in dBm
        rows.append(("Received power", watts_to_dbm(self.get_received_power()), "dBm"))
        # if the receiver has a sensitivity, add it and the link margin at the bottom
        if self.receiver.sensitivity_w is not None:
            rows.append(("Receiver sensitivity",
                         watts_to_dbm(self.receiver.sensitivity_w), "dBm"))
            rows.append(("Link margin", convert_to_db(self.get_margin()), "dB"))
        return rows


# Academic Citations & Literature References:
#
# Received Power Chain & Link Margin:
# Liang, J., Chaudhry, A. U., Erdogan, E., & Yanikomeroglu, H. (2022). Link Budget
# Analysis for Free-Space Optical Satellite Networks. Proc. IEEE 23rd Int. Symp. on a
# World of Wireless, Mobile and Multimedia Networks (WoWMoM). (Specifically Equation 12
# for the link margin).
#
# Received Photon Flux:
# Valentini, L., Faedi, A., Paolini, E., & Chiani, M. (2021). Analysis of Pointing Loss
# Effects in Deep Space Optical Links. Proc. IEEE Global Communications Conf. (GLOBECOM).
# (Specifically Equation 18 for the photon rate at the receiver).
