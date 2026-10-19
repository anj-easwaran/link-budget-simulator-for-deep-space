# one part of the terminal with its mass, power and where the numbers came from
class Component:
    # stores the name, subsystem, mass in g, power in W, duty cycle (0 to 1),
    # provenance (measured, estimated or vendor) and cost in USD if known
    def __init__(self, name, subsystem, mass_g, power_w, duty_cycle,
                 provenance, cost_usd=None):
        self.name = name
        self.subsystem = subsystem
        self.mass_g = mass_g
        self.power_w = power_w
        self.duty_cycle = duty_cycle
        self.provenance = provenance
        self.cost_usd = cost_usd

    # average power in W, the power when on times the fraction of time it is on
    def get_average_power(self):
        return self.power_w * self.duty_cycle


# every part of the NODE terminal in transmit-active mode, mass in g and power in W
# a part that is only listed in one of the two tables gets 0 for the other value
# Reference: Kingsbury (2015), Table 3.5 and Table 3.6
NODE_COMPONENTS = [
    Component("High-speed electronics (PCB, FPGA, memory)", "electronics",
              60, 0.0, 0.0, "estimated"),
    Component("Modulator/framer (Spartan-6 FPGA)", "electronics",
              0, 1.00, 1.0, "estimated"),
    Component("FEC encoder (offline on microcontroller)", "electronics",
              0, 0.25, 0.0, "estimated"),
    Component("Non-volatile telemetry buffer (SSR)", "electronics",
              0, 0.25, 1.0, "estimated"),
    Component("Fine-steering mechanism + driver", "pointing",
              20, 0.25, 1.0, "estimated"),
    Component("PAT processor", "pointing", 0, 0.25, 1.0, "estimated"),
    Component("Structural mounts (optics bench)", "structure",
              100, 0.0, 0.0, "estimated"),
    Component("Seed laser", "transmitter", 20, 0.0, 0.0, "measured"),
    Component("Seed laser driver", "transmitter", 0, 0.25, 1.0, "estimated"),
    Component("Seed laser TEC", "transmitter", 0, 0.25, 1.0, "measured"),
    Component("Extinction filter", "transmitter", 20, 0.0, 0.0, "measured"),
    Component("EDFA", "transmitter", 150, 5.70, 1.0, "measured",),
    Component("TX collimator", "transmitter", 10, 0.0, 0.0, "measured"),
    Component("Fiber service loops, splices", "transmitter",
              20, 0.0, 0.0, "estimated"),
    Component("Beacon lens and filter assembly", "beacon receiver",
              100, 0.0, 0.0, "estimated"),
    Component("Beacon focal plane array", "beacon receiver",
              50, 2.00, 1.0, "estimated"),
    Component("Radio modem RX (Astrodev Lithium-1)", "radio",
              0, 0.10, 0.95, "vendor"),
    Component("Radio modem TX/RX (Astrodev Lithium-1)", "radio",
              0, 1.50, 0.05, "vendor"),
]

# the table prints the EDFA average as 4.10 W against 5.70 W at full duty (measured in the
# lab while transmitting), so that printed average is used to match the published total
# Reference: Kingsbury (2015), Table 3.6
NODE_PRINTED_AVERAGES = {"EDFA": 4.10}


# adds up the mass of every part in grams
def total_mass_g():
    return sum(c.mass_g for c in NODE_COMPONENTS)


# adds up the average power of every part in W
def total_average_power_w():
    total = 0.0
    # loops through the parts, using the printed average where the table has one (the EDFA)
    for c in NODE_COMPONENTS:
        total += NODE_PRINTED_AVERAGES.get(c.name.split(" (")[0],
                                           c.get_average_power())
    return total


# prints the parts table with the totals and the published numbers to check against
def report():
    # column headings
    print(f"{'Component':<44} {'Subsys':<16} {'Mass g':>7} {'P avg W':>8}"
          f" {'Provenance':<10} {'Cost'}")
    # loops through each part and prints one row
    for c in NODE_COMPONENTS:
        # average power, taking the printed table value where there is one
        p_avg = NODE_PRINTED_AVERAGES.get(c.name.split(" (")[0],
                                          c.get_average_power())
        # if no cost is known, say it is unsourced instead of printing a number
        cost = "n/a (unsourced)" if c.cost_usd is None else f"${c.cost_usd}"
        print(f"{c.name:<44} {c.subsystem:<16} {c.mass_g:>7.0f}"
              f" {p_avg:>8.2f} {c.provenance:<10} {cost}")
    # totals row
    print(f"{'TOTAL':<44} {'':<16} {total_mass_g():>7.0f}"
          f" {total_average_power_w():>8.2f}")
    # published mass and power totals to compare the table with
    # Reference: Kingsbury (2015), Table 3.5 and Table 3.6
    print("\nPublished check: 550 g of a 600 g budget; 8.52 W printed total"
          " (component averages sum to 8.53 - 0.01 W rounding in the"
          " duty-cycled radio rows) of a 10 W budget. Kingsbury 2015,"
          " Tables 3.5-3.6.")


# only prints the table when the file is run directly, not when it is imported
if __name__ == "__main__":
    report()


# Academic Citations & Literature References:
#
# NODE Terminal Mass & Power Budget:
# Kingsbury, R. W. (2015). Optical Communications for Small Satellites (Doctoral
# dissertation). Massachusetts Institute of Technology. (Specifically Table 3.5 for the
# mass and Table 3.6 for the power in transmit-active mode).
