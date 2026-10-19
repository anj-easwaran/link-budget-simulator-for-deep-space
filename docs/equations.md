# Equation Register

Every formula implemented in the simulator is recorded here first: statement,
symbol definitions, units, and the exact source (paper + equation number).
**No code is written that is not backed by an entry in this file.**

Rules:
- A value or formula is only used if it comes from a paper in `Literature/`
  (or a web-searched source, cited with enough detail to find it again).
- Nothing is tuned to match measured mission data without an independent source.
- Where two sources disagree, both are recorded and the choice is justified.

Sources are cited by short name; full references at the bottom.

Section numbers are stable identifiers (code comments cite them), so entries
removed in the 2026-09-11 cleanup leave gaps: 1.3 (point-ahead), 1.4-1.6
(pass geometry), 4.4 (pointing outage), 7.1-7.4 (APD noise, OOK Q-factor)
and 8.1 (OOK Gaussian BER) covered code that no thesis result used.

---

## 1. Geometry

### 1.1 Slant distance, ground station to satellite
d_GS = R · ( sqrt( ((R+H)/R)^2 − cos^2(θ_E) ) − sin(θ_E) )

- R = R_E + h_E : Earth radius + ground-station altitude [m]
- H = h_S − h_E : satellite altitude above the station [m]
- θ_E : elevation angle [rad]
- At θ_E = 90°, d_GS = h_S − h_E.

Source: Liang 2022, Eq. (4) (they cite Cakaj 2011).

### 1.2 Atmospheric path length through a layer
d_A = (h_A − h_E) · sec(θ_Z) = (h_A − h_E) · csc(θ_E)

- h_A : height of the atmospheric layer [m]; θ_Z = 90° − θ_E : zenith angle.
- This is the plane-parallel (flat-Earth secant) airmass; valid away from the
  horizon.

Source: Liang 2022, text below Eq. (10) (they cite Guo 2010).

---

## 2. Transmitter

### 2.1 Transmit antenna gain (Klein & Degnan)
G_T = (4π·A / λ²) · g_T(α, γ),  with A = π a²

Far-field, on-axis, collimated (β = 0):
g_T(α, γ) = (2/α²) · [ exp(−α²) − exp(−α²γ²) ]²

- a : primary aperture radius [m]; b : central obscuration radius [m]
- α = a/w : aperture radius / beam 1/e² intensity radius at the aperture
- γ = b/a : obscuration ratio
- The gain definition charges truncation and obscuration losses to G_T itself
  (power normalised over the whole incident plane, not the transmitted part).
- Used for any terminal that publishes its aperture. DSOC: D = 22 cm with the
  measured truncation ratio α = 1.23 (1.18 favourable, 1.51 adverse; Biswas &
  Srinivasan 2026, Sec. IV.E) gives g_T = −0.95 dB, the η_gain allocation of
  their Table III, and G_T = 112.03 dB. (Table III prints the ideal gain
  112.98 dB and η_gain −0.95 dB as separate lines; 112.03 dB is their sum,
  not a printed value. Likewise the 140.12 dB Hale gain quoted in
  `missions.py` is (πD/λ)² for 5 m at 1550.12 nm; the table prints the net
  139.67 dB with obscuration included.) Without a published α the optimum of
  2.2 is used.

NOTE on transcription: the PDF text layer garbles Eq. (9); the squared-bracket
form above is confirmed numerically against the paper's own stated values
(α = 1.12, γ = 0 gives g_T = 0.8145 — reproduced as 0.8147).

Source: Klein & Degnan 1974, Eqs. (6), (9), (13).

### 2.2 Optimum beam truncation
α_opt = 1.12 − 1.30 γ² + 2.12 γ⁴   (accurate to ~1% for γ < 0.4)

Source: Klein & Degnan 1974, Eq. (12).

### 2.3 Far-field beam half-angle of the optimised antenna
θ_1/e² = 0.94 λ/D (γ = 0)  …  0.79 λ/D (γ = 0.4)

- D = 2a : aperture diameter; half-angle to the 1/e² intensity point.

Source: Klein & Degnan 1974, Sec. III (Fig. 10 limits).

### 2.4 Divergence-based transmit gain (simple model)
G_T = 16 / Θ_T²

- Θ_T : FULL transmitting divergence angle [rad], in Liang's convention.
- CORRECTION (2026-09-10): this is NOT the ideal-Gaussian relation for a width
  quoted between the 1/e² intensity points. A Gaussian far field
  I ∝ exp(−2θ²/θ_h²), with θ_h the 1/e² half-angle, has peak gain 8/θ_h² =
  32/Θ² for the full 1/e² width; 16/Θ² corresponds to a full width between the
  1/e points. Applied to DSOC's 1/e² divergence (13.8 µrad) it put the
  transmit gain 2.8 dB below DSOC's own link budget (Table III). It is kept
  only to reproduce Liang's arithmetic (and Liang-class terminals quoted in
  Liang's convention); DSOC is built from its aperture with 2.1, and a width
  known to be a full 1/e² width uses 2.5.

Source: Liang 2022, text after Eq. (1) (they cite Polishuk & Arnon 2004).

### 2.5 Divergence-based transmit gain, full 1/e² width
G_T = 8 / θ_h² = 32 / Θ²

- θ_h : far-field half-angle to the 1/e² intensity point; Θ = 2θ_h.
- Derivation from 2.1: for an untruncated Gaussian (α = a/w → ∞) Klein &
  Degnan's g_T → 2/α² = 2w²/a², so G_T = (4π·πa²/λ²)(2w²/a²) = 8π²w²/λ².
  With the Gaussian far-field half-angle θ_h = λ/(πw), G_T = 8/θ_h².
- Check against a published budget: DSOC's 13.8 µrad 1/e² divergence gives
  32/Θ² = 112.25 dB, against the 112.03 dB of its truncated aperture
  (Biswas & Srinivasan 2026, Table III); truncation widens the beam slightly.
- Used for the Surrey uplink solver, where the width is the diffraction
  width of 2.3 combined with seeing (5.9), taken in the same convention.

Source: Klein & Degnan 1974, Eq. (9), in the untruncated limit.

---

## 3. Free-space channel

### 3.1 Free-space path loss (space loss)
L_PS = ( λ / (4π d) )²    (linear, < 1)

- d : link distance [m]. Same expression for inter-satellite (Eq. 2) and
  ground-satellite (Eq. 5) links.

Source: Liang 2022, Eqs. (2), (5); identical in Klein & Degnan 1974 range
equation (p. 2134) and Valentini 2021, Eq. (18).

### 3.2 Received power (link equation)
P_R = P_T · η_T · η_R · G_T · G_R · L_T · L_R · L_A · L_P

- η_T, η_R : transmit/receive optics efficiencies [-]
- L_T, L_R : transmit/receive pointing losses [-]
- L_A : atmospheric attenuation (unity for ISL / exo-atmospheric)
- L_P : free-space path loss [-]

Source: Liang 2022, Eqs. (1), (3); same structure as Klein & Degnan 1974
range equation P_S = P_T G_T G_R T_T T_R T_A (λ/4πR)².

### 3.3 Received signal photon flux
n_s = P_R · λ / (h c)    [photons/s]

- P_R from 3.2 with all losses (incl. detection/implementation losses) applied.

Source: Valentini 2021, Eq. (18).

### 3.4 Link margin
LM = P_R / P_req

- P_req : receiver sensitivity — power required for the target BER at the
  operating data rate. In photon-flux terms: LM[dB] = n_s[dB] − n_s_min[dB].

Source: Liang 2022, Eq. (12); Valentini 2021, Eq. (19).

---

## 4. Pointing

### 4.1 Angular error statistics
Rayleigh (zero-mean i.i.d. Gaussian azimuth/elevation errors):
f_θ(θ) = (θ/σ_θ²) · exp(−θ²/2σ_θ²),  θ ≥ 0

Rician (non-zero-mean, bias η):
f_θ(θ) = (θ/σ_θ²) · exp(−(θ²+η²)/2σ_θ²) · I₀(θη/σ_θ²)

Source: Valentini 2021, Eqs. (1), (2).

### 4.2 Pointing loss, Gaussian beam
L_p(θ) = exp(−G θ²)

- G : linear on-axis antenna gain; θ : angular miss-pointing [rad].
- Self-consistency: with G = 16/Θ_T², L_p = exp(−16 θ²/Θ_T²).
- Liang applies this at both ends: L_T = exp(−G_T θ_T²), L_R = exp(−G_R θ_R²).

Source: Valentini 2021, Eq. (3); Liang 2022 (after Eq. (1)); both cite
Chen & Gardner 1989.

### 4.3 Pointing loss, circular (Airy) aperture
L_p(θ) = [ 2 J₁(√G θ) / (√G θ) ]²

- Valentini's Eq. (5) approximates this for small angles as exp(−ζ G θ²)
  with ζ = 0.188 (their fitted value). The series expansion of the exact
  form is exp(−G θ²/4) as θ → 0 (2J₁(x)/x ≈ 1 − x²/8). The two differ by
  under 0.3 dB for √G θ ≤ 1; the exact form is what the code implements.

- The Gaussian-beam model (4.2) is the WORST CASE of the two (loss factor
  e^{−Gθ²} vs e^{−0.188Gθ²}); which applies depends on whether the far-field
  pattern is set by Gaussian illumination or a uniformly-filled aperture.
- Used for the DSOC-class transmitter and the closure study (spacecraft
  transmitters built from an aperture). The exact form above is implemented;
  for small θ it tends to exp(−Gθ²/4). With DSOC's aperture gain (2.1) and a
  1 µrad bound it gives −0.18 dB, against the −0.19 dB design pointing loss
  of Biswas & Srinivasan 2026, Table III. The optimum gain under this model
  is found numerically (the flux peaks inside the main lobe, where
  √G θ ≈ 1.84). LLCD-class and Liang-class terminals keep the Gaussian worst
  case of 4.2, which also underlies G_opt = 1/θ_max² (Valentini Eq. 17).

Source: Valentini 2021, Eqs. (4), (5), Fig. 1.

### 4.5 Optimum antenna gain under pointing error (deterministic)
G_opt = 1/θ_max²

- Beyond this gain, narrower beams lose more to jitter than they add on axis.

Source: Valentini 2021, Eq. (17).

---

## 5. Atmosphere

### 5.1 Mie scattering (elevation-dependent)
I_m = exp(−ρ / sin(θ_E)),
ρ = a h_E³ + b h_E² + c h_E + d  (extinction ratio; h_E in km, valid 0–5 km)

Wavelength-dependent empirical coefficients (λ in µm):
a = −0.000545λ² + 0.002λ − 0.0038
b =  0.00628λ² − 0.0232λ + 0.00439
c = −0.028λ² + 0.101λ − 0.18
d = −0.228λ³ + 0.922λ² − 1.26λ + 0.719

- VALIDITY: ρ is a fitted polynomial and turns negative (a gain) above
  about 1.5 km site altitude at 1550 nm (1.7 km at 532 nm). The code clamps
  ρ ≥ 0, so elevated sites (Palomar 1712 m, Table Mountain 2286 m) receive
  zero Mie extinction rather than a gain. Results at such sites carry only
  the Kim term of 5.2 and are a stated limitation; DSOC's Palomar budget
  uses a −0.5 dB allowance instead (`missions.py`).

Source: Liang 2022, Eqs. (6), (7) (they cite ITU-R P.1622).

### 5.2 Geometrical (fog/cloud) scattering — Kim model
V = 1.002 / (L_W · N)^0.6473          (visibility, km)
θ_A = (3.91/V) · (λ/550)^(−φ)          (attenuation coefficient; λ in nm)
I_g = exp(−θ_A · d_A)                  (Beer–Lambert over path d_A)

- L_W : liquid water content [g/m³]; N : cloud number concentration [cm⁻³]
- φ : Kim particle-size coefficient.

Source: Liang 2022, Eqs. (8)–(10).

### 5.3 Total atmospheric attenuation (baseline model)
L_A = I_m · I_g

Source: Liang 2022, Eq. (11).

### 5.4 Refractive-index structure profile (Hufnagel–Valley)
C_n²(h) = 0.00594 (w/27)² (10⁻⁵h)¹⁰ e^(−h/1000)
         + 2.7×10⁻¹⁶ e^(−h/1500) + A e^(−h/100)   [m^(−2/3)], h in m

- HV5/7: w = 21 m/s (rms pseudowind), A = 1.7×10⁻¹⁴ m^(−2/3)
  (gives Fried r₀ ≈ 5 cm, isoplanatic angle ≈ 7 µrad at 0.5 µm — the
  values that name the "5/7" parameterisation).
- ERRATUM (sources disagree): Casanovas 2021 Eq. (47) prints the
  ground-layer term as A·e^(−h/1000); the standard Hufnagel–Valley form
  (Andrews & Phillips, "Laser Beam Propagation through Random Media,"
  2nd ed., Ch. 12; ITU-R P.1621) has A·e^(−h/100). Only the 100 m scale
  height reproduces the definitional r₀ = 5 cm (the printed 1000 m gives
  1.44 cm — checked in validation via §5.10), so the standard form is
  used and Casanovas's is recorded as a misprint.

- HEIGHT REFERENCE: h is height above LOCAL GROUND, and the integrals of
  §5.5, §5.6 and §5.10 run from the station upward. ITU-R P.1621-2 (2015),
  §5.1.1, Eq. (6), p. 8 defines h as "height above ground level (m)", and
  its Eqs. (8a)/(8b) take h₀ as the station height above ground level.
  Casanovas 2021 is inconsistent: Fig. 3.16 plots the profile against height
  above ground, while the text at Eq. (46) calls h₀ the station height above
  sea level. Reading h as height above sea level and integrating from the
  site altitude removes the ground layer at any elevated site (Palomar,
  Table Mountain, Tenerife), giving r₀ ≈ 18-20 cm against 2.7-7.4 cm
  measured.
  Disclosed caveats: ITU-R P.1621-2 itself labels its Fig. 5 axis "height
  above sea level" and states its closed-form r₀ approximation (Eqs. 9-13)
  for station altitudes above sea level. Mata Calvo (DLR/ESA manuscript,
  2023, unpublished) argues the standard profile is sea-level referenced and
  proposes moving only the ground layer to the site. The convention follows
  the ITU definition; the measured r₀ comparison prompted the check but does
  not justify the choice.

Source: Casanovas 2021, Eq. (47), ground term corrected per Andrews &
Phillips / ITU-R P.1621; height reference per ITU-R P.1621-2 (2015), Eq. (6).

### 5.5 Downlink scintillation index (Andrews plane-wave theory)
σ_I² = exp[ 0.49σ₁²/(1+1.11σ₁^(12/5))^(7/6)
          + 0.51σ₁²/(1+0.69σ₁^(12/5))^(5/6) ] − 1

Rytov variance (slant path, plane wave):
σ₁² = 2.25 µ_pl k^(7/6) (H−h₀)^(5/6) sec^(11/6)(φ)
µ_pl = ∫_{h₀}^{H} C_n²(h) ξ^(5/6) dh,   ξ = (h−h₀)/(H−h₀)

- k = 2π/λ wavenumber; φ zenith angle; H satellite altitude; h₀ station.
- Valid for plane waves, any zenith angle; Casanovas notes the model is
  validated for LEO/GEO H — using it for deep space assumes the turbulence
  is all in the low atmosphere (H-independent µ_pl integral converges).

Source: Casanovas 2021, Eqs. (43)-(46), citing Andrews et al.

### 5.6 Aperture averaging (plane wave, slant path)
Λ_avg = [1 + 1.11 (D_Rx²/(h_s λ sec φ))^(7/6)]⁻¹
h_s = [ ∫C_n²(h)(h−h₀)² dh / ∫C_n²(h)(h−h₀)^(5/6) dh ]^(6/7)
σ_P² = Λ_avg · σ_I²   (power scintillation index)

- ERRATUM (sources disagree): Casanovas 2021 Eq. (48) prints cos φ in the
  denominator, and codes it that way in his appendix. Biswas & Piazzolla
  (JPL DESCANSO vol. 7, Ch. 3, Eq. (3.4-18)), citing the same original
  (Yura & McKinley, Appl. Opt. 22(11), 1608-1609, 1983), have sec φ in the
  denominator. The sec form is used: a longer slant path puts the
  turbulence further from the receiver, widening the Fresnel zone, so
  averaging weakens toward the horizon. Two independent measured campaigns
  show that trend (KIODO, Perlot 2007; SOTA/OPALS, Moll 2016); the cos form
  reverses it.
- UNRESOLVED: the coefficient is 1.11 in Casanovas but A₀⁻¹ with A₀ ≈ 1.1
  (≈ 0.91) in DESCANSO. 1.11 is kept; the difference moves the §5.7 fading
  loss by 0.1-0.4 dB for a 40 cm aperture. The 1983 original was not read.

Source: Biswas & Piazzolla (DESCANSO vol. 7, Ch. 3), Eqs. (3.4-18)-(3.4-19),
after Yura & McKinley 1983; scale height and σ_P² per Casanovas 2021,
Eqs. (49)-(51).

### 5.7 Scintillation (fading) loss for a target outage probability
Lognormal received-power fading with σ_LD² = ln(σ_P²+1), µ_LD = −σ_LD²/2:
p_thr = ½[1 + erf( ln(P_min/⟨P⟩ · √(σ_P²+1)) / √(2 ln(σ_P²+1)) )]

Inverted for the loss to budget at outage fraction p_thr:
L_turb(dB) = 4.343 [ erf⁻¹(2p_thr−1)·√(2 ln(σ_P²+1)) − ½ ln(σ_P²+1) ]

Source: Casanovas 2021, Eqs. (52)-(58).

### 5.8 Background light power on the detector
P_BGL = L_e,Ω,λ · A · Ω · Δλ,   Ω = 4π sin²(ω)

- L_e,Ω,λ : background spectral radiance [W/(m² nm sr)] — tabulated values
  (Giggenbach Table 1), e.g. at 1550 nm: blue sky zenith 1.2×10⁻³,
  horizon 25×10⁻³, sunlit cloud 20×10⁻³.
- A : aperture area [m²]; Δλ : filter bandwidth [nm]; ω : detector full-flat
  field-of-view angle [rad].
- NOTE: Ω = 4π sin²ω is transcribed as printed; it reproduces all four
  worked P_BGL values in the paper's own Table 2 (50/88/2.4/4.2 nW), which
  is the authority for how ω is meant to be used.

Source: Giggenbach 2022, Eq. (22), Tables 1-2.

### 5.9 Uplink effective divergence with turbulence (GOPEX model)
Θ_eff² = Θ² + θ_seeing²
S_day / S_1 = (d_1/d)² · (Θ_eff,1/Θ_eff)² · T_z^(sec ζ_day − sec ζ_1)

- The GOPEX analysis models atmosphere-induced angular jitter as two
  independent zero-mean Gaussian random variables; independent Gaussian
  angular spreads add in variance, broadening the effective beam. The
  received signal then scales as 1/(Θ_eff² d²).
- Clear-air transmittance along the slant path is T_z^(sec ζ): optical depth
  scales with path length (Beer–Lambert, DESCANSO Eqs. 3.2-1 to 3.2-3), over
  the same plane-parallel path as §5.1. ζ is each night's average zenith
  angle (Levine Table 2, 50–56°).
- T_z = 0.825 ± 0.01 at 532 nm, digitised from DESCANSO Fig. 3-18 (model
  transmittance, 3 km site, zenith, 23 km visual range). ASSUMED to hold at
  TMF (2.2 km) on the GOPEX nights; a lower-altitude site or the elevated
  post-1991 stratospheric aerosol would lower it. Sensitivity: sec ζ spans
  at most 0.23 across the nights, so T_z = 0.825 moves a ratio by at most
  0.19 dB, and T_z = 0.70 would move it by at most 0.36 dB.
- §5.1–5.3 are not used here: the §5.1 Mie polynomial gives negative
  extinction (a gain) above about 1.7 km site altitude.
- Convention caveat: the GOPEX tables do not state whether "seeing" and
  "divergence" share the same width convention; treated as same-convention
  widths. The resulting day-to-day ratios agree with the paper's own
  Kiasaleh & Yan "expected photoelectrons" to <0.6 dB, which supports the
  reading.

Source: Levine, Shaik & Yan (GOPEX data analysis), TDA PR 42-114, 1993,
Sec. IV and Table 2 (`19940009915.pdf`); transmitter parameters from
Wilson & Lesh, same volume (`19940009913.pdf`), Table 1; clear-air
transmittance from Biswas & Piazzolla 2006 (DESCANSO vol. 7, ch. 3),
Eqs. (3.2-1)–(3.2-3) and Fig. 3-18.

### 5.10 Fried parameter (atmospheric coherence length), plane wave
r₀ = [ 0.423 k² sec(ζ) ∫ C_n²(h) dh ]^(−3/5)

- k = 2π/λ; ζ zenith angle; integral over the station-to-space path.
- Validation anchor: HV5/7 profile (§5.4) at λ = 0.5 µm, zenith →
  r₀ = 5 cm by definition of the parameterisation.

Source: Fried 1966 (definition); standard form as in Andrews & Phillips,
Ch. 12 (chapter-level cite, textbook-universal); anchor value also quoted
by Casanovas 2021 with the HV5/7 definition.

### 5.11 GOPEX measured mean signal and its sampling error
⟨N⟩ = exp(m_l + σ_N²/2)
Δ_dB(day) = 10 log₁₀( ⟨N⟩_day / ⟨N⟩_1 )
CV = σ_I/⟨I⟩ = √(exp(σ_N²) − 1)
u(ln⟨N⟩) = CV / √n
u(Δ_dB) = (10/ln 10) · √( u_day² + u_1² )

- m_l, σ_N: lognormal parameters fitted to each night's received pulses
  ("data" rows of Levine Table 2); n: total received pulses (same table).
- The link budget, like the paper's own theory, predicts the MEAN
  photoelectrons per pulse, whereas m_l is the logarithm of the median. The
  first line is Levine Eq. (6) solved for ⟨N⟩. The appendix m-file
  (Fig. A-2) computes m_l and σ_N from each night's sample mean and standard
  deviation in exactly this way, so the relation recovers the measured mean.
- CV: Levine Eq. (12) solved for the ratio of standard deviation to mean.
- Sampling error: standard deviation of the mean (GUM 4.2.3, Eq. 5), carried
  to the logarithm to first order and combined for the two independent
  nights in quadrature (GUM 5.1.2, Eq. 10). Day 1's term is common to every
  ratio; on day 1 itself the ratio is 1 by definition and carries no bar.
- ASSUMED: pulses within a night are independent (15–30 Hz repetition).
  Slow changes in seeing through a night correlate them, so u is a lower
  bound on the uncertainty of a night's mean.

Source: Levine, Shaik & Yan 1993, Eqs. (6), (12), Table 2, Fig. A-2
(`19940009915.pdf`); JCGM 100:2008 (GUM), Sec. 4.2.3 Eq. (5), Sec. 5.1.2
Eq. (10).

---

## 6. Receiver

### 6.1 Receive antenna gain (ideal filled aperture)
G_R = (π D_R / λ)²

- D_R : receiver telescope diameter [same length units as λ].

Source: Liang 2022, text after Eq. (1) (they cite Polishuk & Arnon 2004);
consistent with the diffraction limit 4πA/λ² of Klein & Degnan 1974 for
A = π(D/2)²: (4π/λ²)·π D²/4 = (πD/λ)².

### 6.2 Central obscuration of the receive aperture
G_R multiplied by (1 − γ²), γ = b/a.

- Pure geometry: the secondary mirror blocks the area fraction γ² of the
  collecting aperture (energy collection, direct detection). Klein &
  Degnan's Part 2 (NASA TM X-524-73-185) treats the full receiver case.

---

## 7. Detection (photon counting)

### 7.5 Photon-counting statistics (Poisson channel)
f(k|0) = K_b^k e^(−K_b) / k!,   f(k|1) = (K_s+K_b)^k e^(−(K_s+K_b)) / k!

- K_s, K_b : mean detected signal/background photons per slot (detection
  efficiency and dark counts already folded into the means:
  K_b includes dark counts per slot, K_s = η_det × incident photons).

Source: Hemmati 2011, Eq. (3) and Appendix; same statistics as
Valentini 2021 Eq. (14).

---

## 8. Modulation and BER

### 8.2 Poisson PPM channel (photon counting)
P{k | slot pulsed}     = (n_s + n_b)^k e^{−(n_s+n_b)} / k!
P{k | slot not pulsed} = n_b^k e^{−n_b} / k!

- M-ary PPM: symbol = M = 2^m slots, one pulsed; m bits/symbol.
- n_s, n_b : mean signal/background photons per slot.

Source: Valentini 2021, Eq. (14).

### 8.3 SCPPM information data rate
B_r = (15120·R_c − 34) / ( M·T_s·(15120/log₂M)·α_gt )

- R_c : convolutional code rate; T_s : slot time; α_gt : guard time factor
  (CCSDS: M/4 idle slots between symbols, α_gt = 25%);
  15120 = SCPPM codeword length, 34 = CRC + termination bits.

- NOTE: superseded by §9.3. Eq. (20) prints α_gt alone, but only (1+α_gt)
  reproduces the paper's own data-rate rows; the code uses §9.3.

Source: Valentini 2021, Eq. (20); parameters from CCSDS 142.0-B-1.

### 8.4 OOK on the Poisson channel (photon counting)
ML threshold:  y_tsh = K_s / ln(1 + K_s/K_b)
P_b = ½ − ½ Σ_{k=0}^{⌊y_tsh⌋} [f(k|0) − f(k|1)]
No background:  P_b = ½ e^(−K_s)

Source: Hemmati 2011, Appendix Eqs. (A1)-(A3).

### 8.5 Uncoded M-PPM symbol error on the Poisson channel (exact, ML)
P_s = 1 − Σ_{k=0}^{∞} [ (1+K_s/K_b)^k e^(−K_s) / M ]
          × [ F(k|0)^M − F(k−1|0)^M ]
F(k|0) = Σ_{m=0}^{k} K_b^m e^(−K_b) / m!     (background CDF)
No background:  P_s = (M−1) e^(−K_s) / M

- Verified analytically: (A6) equals the ML-with-random-tie-breaking
  expression via the identity Σ_j C(M−1,j)/(j+1) p^j F^(M−1−j)
  = [(p+F)^M − F^M]/(M p).

Source: Hemmati 2011, Appendix Eqs. (A5)-(A7).

### 8.6 PPM symbol-to-bit error mapping
P_b ≈ M/(2(M−1)) · P_s

Source: Li 2018, Eq. (3), citing Hemmati (ed.), "Deep Space Optical
Communications," Wiley 2006, ch. 4.

---

### 8.7 Coherent detection, shot-noise-limited BER
Q = |d1 − d0|/(σ1 + σ0),  BER = ½ erfc(Q/√2)

- BPSK homodyne:   Q = 2√n̄s → BER = ½ erfc(√(2 n̄s));  n̄s = 4.55 ph/bit
  at BER 10⁻⁵ (paper's own anchor value).
- BPSK heterodyne (or phase/polarisation-diversity intradyne = digital
  coherent): Q = √(2 n̄s) → BER = ½ erfc(√n̄s); n̄s = 9.0 at 10⁻⁵.
  The 3 dB penalty is the image-band vacuum noise, and the phase-diversity
  homodyne (intradyne) receiver has identical S/N.

- n̄s : mean detected signal photons per bit.

Source: Kikuchi 2016, Eqs. (57)-(60), (74).

### 8.8 Laser phase noise (Lorentzian line)
White FM noise → Lorentzian field spectrum of FWHM δf, and the phase
drift accumulated over an interval T is Gaussian with variance
⟨Δφ²(T)⟩ = 2π · δf · T

Source: Kikuchi et al. 1984, Eqs. (40)-(46).

### 8.9 Phase-error variance by receiver scheme
The IF carrier linewidth is δf_T + δf_LO (transmitter + LO, Eq. 52).

- DPSK (delay demodulation, bit period T = 1/f_b):
  σ_φ² = 2π (δf_T + δf_LO) T               [Eq. 60]
  (Eq. (62) as printed, σ_φ² = (δf_T + δf_L)/f_s with Eq. (61) f_s = 2π/T,
  contradicts Eq. (60) by a factor (2π)²; only the Eq. (60) form gives the
  paper's own "(1/T) > 3 Gbit/s" below, so Eq. (60) is cited.)
- Phase-error budget: σ_φ ≤ 0.2 rad keeps the receiver-sensitivity loss at
  BER 10⁻⁹ below 3 dB (Sec. V-A, from their Figs. 2–3). For DPSK this gives
  the largest combined linewidth δf_T + δf_LO ≤ σ_φ² f_b / (2π)
  (`max_linewidth_for_dpsk`).
- Reproduction anchor: with δf_T = δf_LO = 10 MHz the paper requires
  "(1/T) > 3 Gbit/s" for heterodyne DPSK (Sec. V-B); the relation above
  gives 3.14 Gbit/s (validation check).
- Laser linewidths used for comparison, per laser (a matched transmitter/LO
  pair doubles them): DFB semiconductor lasers 100 kHz – 10 MHz (Kikuchi
  2016, "Carrier-Phase Estimation" subsection); the stabilised source
  100 Hz, the specification in the project brief for the hardware planned
  for Objective 4 (a specification, not a measured linewidth).
- Not modelled: a direct-detection DPSK receiver with an optical delay-line
  interferometer uses no local oscillator, so only δf_T enters.

Source: Kikuchi et al. 1984, Section III and Section V; Kikuchi 2016
(DFB linewidth range).

### 8.10 PSK error floor from residual phase error
Conditioned on Gaussian phase error Δφ ~ N(0, σ_φ²) the PSK BER is
E_Δφ[ ½ erfc( √γ · cos Δφ ) ] (Eqs. (19)-(20)); only its infinite-SNR
limit is implemented, the error floor where the decision flips for
|Δφ| > π/2:
BER_floor = erfc( π / (2√2 σ_φ) )
(the paper presents it as Fig. 4; used for the DPSK high-SNR limit of
§8.11 and the linewidth feasibility table in `stabilised_laser.py`)

Source: Kikuchi et al. 1984, Eqs. (18)-(20), Fig. 4.

### 8.11 Heterodyne DPSK BER with residual phase error
P_e = (4/√π) ∫₋∞^∞ e^(−t²) ∫₀^∞ ε I₀(2√α ε) e^(−ε²−α)
          · ∫_ε^∞ ξ I₀(2√β ξ) e^(−ξ²−β) dξ dε dt
α = 2γ cos²(σ_φ t/√2),   β = 2γ sin²(σ_φ t/√2)

- DPSK uses the previous bit as its phase reference, so its BER is not the
  PSK expression of §8.10; the two envelopes R₁, R₂ are Rician (Eqs.
  (30)–(31)) and averaging over the Gaussian phase error gives Eq. (36).
- γ is the same S/N ratio as §8.10; on the shot-noise-limited heterodyne
  channel γ = detected photons per bit (Kikuchi 2016 Eq. (74) for PSK).
- Evaluation: the ξ integral is a Marcum Q function,
  ∫_ε^∞ ξ I₀(2√β ξ) e^(−ξ²−β) dξ = ½ Q₁(√(2β), √2 ε), written with the
  non-central χ² survival function (2 degrees of freedom, non-centrality
  2β, at 2ε²). The ε integral runs over ±10 about the Rician peak √α. The
  average over the phase error φ = √2 σ_φ t uses a dense grid on
  0 ≤ φ ≤ 8σ_φ, because the errors come from the tail |φ| > π/2, which a
  coarse Gauss–Hermite rule under-samples. Above γ = 10⁵ the BER is set to
  its Eq. (39) limit, which the integral is within 0.1% of there and which
  keeps the special functions finite.
- Checks: with σ_φ = 0 the triple integral reduces analytically to the
  standard P_e = ½ e^(−γ) (verified analytically); the
  σ_φ = 0.20 rad curve reaches BER 10⁻⁹ at γ = 16.03 dB against 16.1 dB
  digitised from Fig. 3 (validation check, reproduction); large γ
  approaches the Eq. (39) floor of §8.10 from above. At σ_φ = 0.2 rad the
  phase-noise penalty at the 3.8 × 10⁻³ threshold is 0.34 dB.
- At the 3.8 × 10⁻³ FEC threshold (§9.3a) DPSK needs γ = 4.88 against 3.56
  for PSK, 1.37 dB more. `dpsk_max_data_rate` previously used the PSK
  expression and so overstated the DPSK rate by that factor.

Source: Kikuchi et al. 1984, Sec. II-E, Eqs. (23)–(38), Fig. 3.

## 9. Coding

### 9.1 Soft-decision capacity of M-ary PPM on the Poisson channel
C(M) = (1/M) · E_{Y₁..Y_M}[ log₂( M Λ(Y₁) / Σ_j Λ(Y_j) ) ]   [bits/slot]

- Λ(y) = f(y|1)/f(y|0) = (1+K_s/K_b)^y e^(−K_s) is the slot likelihood
  ratio; Y₁ ~ f(·|1) (pulsed slot), Y_j ~ f(·|0) for j>1.
- NOTE: Faedi prints L = p₀/p₁ in his Eq. (1.2) but the capacity formula
  only reproduces the known noiseless limit C → (1−e^(−K_s))·log₂M / M
  with the Λ = p₁/p₀ orientation (verified analytically); recorded so.
- Evaluated numerically (Monte Carlo over the Poisson statistics).
- Same expression as Hemmati 2011, Eq. (10), which is written per symbol:
  C_sym = log₂M − E[ log₂ Σ_j (1+K_s/K_b)^(Y_j−Y₁) ] = M·C, because
  Λ(Y_j)/Λ(Y₁) = (1+K_s/K_b)^(Y_j−Y₁) (the e^(−K_s) factors cancel).
  Hemmati's Fig. 5 plots C per slot at K_b = 1; eight rising-edge crossings
  read off it are a reproduction check (validation.py,
  hemmati_capacity_fig5; worst 0.17 dB against a 0.3 dB tolerance).

Source: Faedi (Bologna thesis), Eq. (1.4), citing Moision & Hamkins,
IPN Progress Report 42-161, 2005.

### 9.2 Code performance as a gap to capacity
A serially concatenated turbo/SCPPM-class code on the 64-PPM Poisson
channel operates ≈ 1 dB from channel capacity (at K_b = 0.2, BER 10⁻⁵);
LDPC ≈ 1.4 dB. The link model treats coded performance as
R_b ≤ C evaluated with the received flux reduced by this gap.

Source: Hemmati 2011, Sec. III (citing [53]); CCSDS 142.0-B-1 defines the
SCPPM code itself (rates 1/3, 1/2, 2/3; M = 4..256).

### 9.3a Generic hard-decision FEC threshold (OOK / coherent links)
A standard-FEC-coded link closes when the UNCODED (pre-FEC) BER is below
3.8 × 10⁻³; the information rate is the channel rate reduced by the code
overhead (7%, the standard FEC family this threshold belongs to —
overhead value is the industry convention, marked as such).

Source: Li 2018, Sec. 4, Conclusion ("the FEC limit of 3.8 × 10⁻³").

### 9.3 SCPPM frame data rate
R_b = (15120·R_c − 34) / [ (15120/log₂M) · M·T_s · (1+α_gt) ]

- NOTE: Valentini's Eq. (20) prints the guard factor as α_gt, but only
  (1+α_gt) reproduces the paper's own data-rate rows (97.00 kbps at M=64,
  32.33 kbps at M=256, both with R_c=1/3, T_s=256 ns, α_gt=25%): the
  guard slots extend each symbol by M/4 slots. Recorded as verified.

Source: Valentini 2021, Eq. (20) + Tables I-II; CCSDS 142.0-B-1.

---

## 10. Uncertainty treatment

Model outputs are point values computed from published parameters plus the
values marked ASSUMED. The uncertainty carried by the ASSUMED values is
propagated by INTERVAL EVALUATION: each ASSUMED parameter is assigned a
(pessimistic, optimistic) range, recorded next to its nominal value in
`missions.py`, and an output bound is obtained by evaluating the full model
once with every parameter at its pessimistic edge and once at its
optimistic edge. Because every ASSUMED parameter acts monotonically on the
outputs used here (rate and margin increase with efficiency, decrease with
dark counts, background field of view and filter width), the two joint
evaluations bracket the output exactly.

Interval bounds are used instead of probabilistic (Monte Carlo)
propagation deliberately: no distributions are published for these
parameters, and assigning them would contradict the register's no-tuning
rule. This is the "worst-case" treatment of JCGM 100 (GUM) Sec. 4.3.7 for
quantities with known limits and unknown distribution (there, a rectangular
distribution's limits; here reported as raw bounds without dividing by
sqrt(3), i.e. bounds, not a standard uncertainty).

Separate, smaller: the PPM capacity estimator (register 9.1) is itself a
Monte Carlo mean over N samples, so quoted rates carry a sampling standard
error of order sigma/sqrt(N) (JCGM 101 Sec. 7.9); with the N = 20,000 used
for figures this is small against the ASSUMED-parameter intervals and is
not propagated further.

Source: JCGM 100:2008 (GUM), Sec. 4.3.7; JCGM 101:2008 (GUM Supplement 1,
Monte Carlo method), Sec. 7.9.

---

## Validation anchors (published worked examples)

1. **Liang 2022, Table 2** — ISL: λ=1550 nm, η_T=η_R=0.8, D_R=80 mm,
   θ_T=θ_R=1 µrad, Θ_T=15 µrad, P_req=−35.5 dBm, LM=3 dB.
   d_SS=1000 km → L_PS=−258.18 dB, P_T=15.32 dBm; …; d_SS=10000 km →
   L_PS=−278.18 dB, P_T=35.32 dBm. (Genuine end-to-end arithmetic check.)
2. **Liang 2022, Tables 4–5** — up/downlink with atmosphere: h_E=1 km,
   L_W=3.128e-4 g/m³ (paper's Table 3 swaps the L_W/N labels vs Eq. 8 usage),
   N=0.5 cm⁻³, φ=1.6, h_A=20 km. θ_E=40° → I_m=−0.15 dB, I_g=−0.33 dB,
   L_A=−0.48 dB; full elevation sweep 10°–90° in Table 5.
3. **Valentini 2021, Table III** — deep-space Mars TM budget: **λ = 1064 nm**
   (the −373.49 dB space loss only reproduces at 1064 nm; at 1550 nm it
   would be −370.24 dB), M=64, R_c=1/3, T_s=256 ns, θ_max=0.35 µrad,
   G_T=G_R=129 dB, η=−5 dB each end, other losses −4 dB, pointing loss
   −8.45 dB, range 2.68 AU, P_av=5 W → n_s=−33.67 dB phe/ns against
   n_s_min=−35.76 dB phe/ns, LM=2.09 dB, 100 kbps.
4. **Klein & Degnan 1974** — α=1.12 optimum for γ=0 with g_T=0.8145.

---

## References (short name → file)

- **Liang 2022** — J. Liang, A.U. Chaudhry, E. Erdogan, H. Yanikomeroglu,
  "Link Budget Analysis for Free-Space Optical Satellite Networks,"
  IEEE WoWMoM 2022. `Link_Budget_Analysis_for_Free-Space_Optical_Satellite_Networks.pdf`
- **Valentini 2021** — L. Valentini, A. Faedi, E. Paolini, M. Chiani,
  "Analysis of Pointing Loss Effects in Deep Space Optical Links,"
  IEEE GLOBECOM 2021. `Analysis_of_Pointing_Loss_Effects_in_Deep_Space_Optical_Links.pdf`
- **Klein & Degnan 1974** — B.J. Klein, J.J. Degnan, "Optical Antenna Gain 1:
  Transmitting Antennas," Applied Optics 13(9), 1974. `OpticalAntennaGain1.pdf`
  (Extended NASA report: `19730019291.pdf`, TM X-524-73-185.)
- **Giggenbach 2022** — D. Giggenbach, "Free-Space Optical Data Receivers with
  Avalanche Detectors for Satellite Downlinks Regarding Background Light,"
  Sensors 22, 6773, 2022. `sensors-22-06773-v5 (2).pdf`
- **Casanovas 2021** — M. Casanovas Ventura, "An Assessment on the Requirements
  for Deep Space Optical Communications," MSc thesis, UPC, 2021.
  `Casanovas_Ventura_Marc_TFM_2021.pdf`
- **Hemmati 2011** — H. Hemmati, A. Biswas, I.B. Djordjevic, "Deep-Space
  Optical Communications: Future Perspectives and Applications,"
  Proc. IEEE 99(11), 2011. `ProcIEEE2011.pdf`
- **CCSDS 142.0-B-1** — "Optical Communications Coding and Synchronization,"
  Blue Book, Aug 2019. `142x0b1.pdf`
- **Li 2018** — B. Li et al., "BER Analysis of a Deep Space Optical
  Communication System Based on SNSPD Over Double Generalized Gamma
  Channel," IEEE Photonics J. 10(5), 2018.
  `ber-analysis-of-a-deep-space-optical-communication-system-1ctizazzc0.pdf`
- **Faedi** — A. Faedi, "Design and Analysis of Optical Links for Space
  Communications," MSc thesis, Università di Bologna. `masters_FAEIDO.pdf`
- **Xie 2016** — H. Xie, D. Heckman, J. Breidenthal, "Link Characterization
  for Deep-Space Optical Communications," IPN Progress Report 42-205, 2016.
  `205D.pdf`
- **Cakaj 2011** — S. Cakaj, W. Keim, K. Malaric, "The Range and Horizon
  Plane Simulation for Ground Stations of Low Earth Orbiting (LEO)
  Satellites," Int. J. Communications, Network and System Sciences, 2011.
  (Cited by Liang 2022 for §1.1.)
