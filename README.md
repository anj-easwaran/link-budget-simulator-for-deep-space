# Deep-Space Free-Space Optical Link Simulator

![Python](https://img.shields.io/badge/Python-3-3776AB?logo=python&logoColor=white)
![NumPy](https://img.shields.io/badge/NumPy-013243?logo=numpy&logoColor=white)
![SciPy](https://img.shields.io/badge/SciPy-8CAAE6?logo=scipy&logoColor=white)
![Validation](https://img.shields.io/badge/validation-65%20checks-brightgreen)
![MSc](https://img.shields.io/badge/MSc%20dissertation-University%20of%20Surrey%202026-blue)

An open-source Python simulator for laser (free-space optical) communication links between spacecraft and ground
stations, from low Earth orbit to deep space. Built for my MSc dissertation in Astronautics and Space Engineering at
the University of Surrey (2026).

**[What it models](#what-it-models) · [Validation](#validation) · [Results](#results) ·
[Running it](#running-it) · [Structure](#repository-structure)**

Laser links can carry far more data than radio, but the link budgets behind the missions that have flown are not
public. This simulator rebuilds them from published equations, checks every step against published calculations and
flight data, and then uses the model to answer design questions.

---

## What it models

The simulator starts from a first-principles link budget and adds effects in stages:

| Stage | What is modelled |
|---|---|
| **Geometry** | slant range and atmospheric path length between a ground station and a satellite |
| **Transmitter and receiver** | telescope gains and beam divergence (Klein and Degnan) |
| **Channel** | free-space path loss, atmospheric attenuation (Mie / Kim), Hufnagel-Valley turbulence and fading, background light |
| **Pointing** | Gaussian and Airy pointing loss, optimum transmit gain for a given pointing error |
| **Modulation and detection** | OOK on a photon-counting channel; PPM with photon-counting detectors such as SNSPDs, with PPM capacity by Monte Carlo; coherent BPSK and DPSK with laser phase noise |
| **Channel coding** | SCPPM data rates, decode thresholds, photon flux needed for a given rate |
| **Missions** | NASA's DSOC and LLCD, the GOPEX uplink experiment, and Surrey's 40 cm ODK16 telescope, one cited source per parameter |

Every formula is written in [`docs/equations.md`](docs/equations.md), with its statement, symbols, units and the
paper and equation number it comes from, before it appears in the code. Nothing is tuned to match measured mission
data.

---

## Validation

`validation.py` runs **65 automated pass/fail checks**:

| Type | Count | What it tests | Result |
|---|---|---|---|
| Reproduction | 48 | equations are implemented correctly | published link budgets reproduced within **0.1 dB** |
| Measurement | 17 | the physics against flight data | see below |

- The **GOPEX** laser uplink is matched to within **0.6 dB** of the experiment's own theory.
- The modelled **DSOC** decode threshold lies **0.25 to 1.6 dB** below the mission's published operating point
  across its background range.
- The modelled **turbulence fading** lies inside the spread measured on 40 cm, 1550 nm ground-station downlinks at
  10° and 20° elevation.

Each check has its own tolerance, set by the type of check. The dissertation also tests which measured checks a
deliberately simple model would pass too, to show which comparisons genuinely test the physics.

---

## Results

| Study | Finding |
|---|---|
| Sensitivity analysis | telescope apertures affect link performance more than any other parameter |
| Coherent modulation | the laser linewidth needed for coherent links at deep-space data rates |
| CubeSat terminal | component-level laser terminal sized at **550 g** and **8.5 W** |
| Surrey ground station | the 40 cm telescope with an assumed photon-counting receiver could receive about **129 kb/s** (32 to 275 kb/s) from a DSOC-class transmitter at 2.68 AU on a clear night |
| Solar activity | NRLMSISE-00 shows solar minimum to maximum changes clear-air loss by less than 10⁻⁴ dB, so it is left out |

### Uncertainty

Parameters that no source publishes are marked `ASSUMED` in `model/missions.py`, with a pessimistic and an
optimistic value recorded beside each one. Results are bracketed by running the model at both ends, and the
scenario printouts report these bounds alongside the nominal values.

---

## Running it

Requires Python 3 with NumPy and SciPy:

```bash
pip install numpy scipy
```

Each script prints its results; the maths stays in `model/`.

```bash
PYTHONPATH=model python validation.py   # 65-check pass/fail suite
python studies/station_scenarios.py     # Surrey ground-station scenarios (LEO, Moon, deep space)
python studies/inverse_solver.py        # hardware needed to close a link into Surrey
python studies/stabilised_laser.py      # coherent modulation and laser linewidth
python model/terminal.py                # CubeSat terminal mass and power breakdown
```

> [!TIP]
> On Windows PowerShell, set the path first with `$env:PYTHONPATH="model"`, then run `python validation.py`.

---

## Repository structure

<details>
<summary>Show files</summary>

| File | Contents |
|---|---|
| `model/constants.py`, `model/utility.py` | SI constants and the only dB conversions |
| `model/geometry.py` | slant range and atmospheric path |
| `model/transmitter.py`, `model/receiver.py` | telescope gains and divergence |
| `model/channel.py` | path loss, atmosphere, turbulence, background light |
| `model/pointing.py` | pointing loss and optimum gain |
| `model/modulation.py` | OOK and PPM error rates, PPM capacity by Monte Carlo |
| `model/coding.py` | SCPPM rates, decode thresholds, rate search |
| `model/coherent.py` | BPSK and DPSK with laser phase noise |
| `model/link_budget.py` | the forward model: received power, photon flux, margin |
| `model/missions.py` | DSOC, LLCD, GOPEX and Surrey parameter sets, one citation per line |
| `model/terminal.py` | CubeSat terminal mass and power budget |
| `validation.py` | the 65-check test suite |
| `studies/` | Surrey scenarios, inverse sizing and the coherent laser study |
| `docs/equations.md` | every equation with its source |

</details>

---

> [!NOTE]
> This is a static link budget: it gives the data rate a link supports under stated conditions, not how often those
> conditions occur. Cloud, seasons, daytime turbulence and solar conjunction geometry are outside its scope.

## Author

**Anchaly Easwaran**, MSc Astronautics and Space Engineering, University of Surrey, 2026.
The full write-up is in [`FSO_PROJ_SUBMISSION.pdf`](FSO_PROJ_SUBMISSION.pdf).
