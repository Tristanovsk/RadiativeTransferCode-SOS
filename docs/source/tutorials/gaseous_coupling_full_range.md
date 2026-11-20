# Gaseous coupling over 364-2500 nm: results

This page summarises the follow-up of the tutorial {doc}`gaseous_coupling`, which tested the
**decoupled** gaseous transmittance of [RadCalNet-OC](https://github.com/Tristanovsk/radcalnet_oc) on
five absorption windows of 684-990 nm, with one water vapour column and two solar angles. The study was
extended with SOS-ABS to the whole spectral range of the CKD coefficients (364-2500 nm), to the water
vapour column, the solar angle and the aerosol load, and a photon path-length model of the coupling was
fitted. The simulations and the analysis are in RadCalNet-OC
(`scripts/sos_abs_coupling`, notebook `notebook/research/radcalnet_oc_gaseous_coupling_350_2500nm.ipynb`);
the figures below are taken from that notebook.

## Simulation sets

All sets are made with {func}`sosabs.build_lut` (script `make_luts.py` of RadCalNet-OC): TOA field with
and without absorption (CKD 5 cm$^{-1}$), over a rough sea (Cox-Munk, 5 m s$^{-1}$) plus a Lambertian
reflector of albedo $\rho$ standing for the water-leaving signal, and the direct-beam transmittance
versus air mass.

| set | spectral range | grid | status |
|---|---|---|---|
| v1 | 5 windows of 684-990 nm | H$_2$O 2 g cm$^{-2}$, SZA 30° and 60°, AOT(550) 0.05 and 0.2 | tutorial {doc}`gaseous_coupling` |
| v2 | 684-990 nm | H$_2$O 0.5-4 g cm$^{-2}$, SZA 20-70°, AOT(550) 0.05-0.4 | done |
| v3 | 364-2500 nm outside 684-990 nm, absorbing CKD intervals ($0.02 < T_g < 0.99$ for $m = 2$) | as v2 | done |
| v4 | wavelengths of v2 + v3 | OPAC aerosols (MOPSMAP phase matrices) and one condition changed at a time: aerosol model and scale height, gas profile, pressure, polarisation, wind, direct glint | nominal set done, variants to run |
| v5 | saturated cores of the H$_2$O 1380 and 1880 nm bands, skipped by v3 ($T_g < 0.02$) | as v3 | running |

Atmosphere of v1-v3: WMO maritime aerosols (scale height 2 km), US62 profiles, 1013 hPa, scalar radiative
transfer. Together, v2 and v3 give 678 absorbing channels in 19 windows.

### Notes on running SOS-ABS for these sets

- The CKD coefficients (`fic/COEFF_CKD/ckd_5cmm1.nc`) start at 364 nm.
- The simulated wavelengths are the centres of the 5 cm$^{-1}$ CKD intervals and avoid the edges of the
  250 cm$^{-1}$ CKD files, where SOS-ABS stops (`ERROR_3500`, e.g. 1.6 µm = 6250 cm$^{-1}$).
- The direct-transmittance runs need `-SURF.Alb 0` even without a surface (`ERROR_2401` otherwise).
- With external aerosol data (`-AER.Model 4`, see {doc}`opac_aerosol_luts`), `-AER.Waref` must equal
  `-SOS_Main.Wa`: the runs are made wavelength by wavelength, with the AOT at the wavelength computed from
  the OPAC extinction.
- The Mie and surface-matrix cache is shared by all the workers; the surface matrices are cached per solar
  zenith angle (LUTs made with older `sosabs` versions, which reused them across angles, must be rebuilt).

## Gaseous absorption over 364-2500 nm

```{figure} figures/coupling_full_absorption_350_2500nm.png
:width: 100%

Direct-beam transmittance of each gas for an air mass of 2 (RadCalNet-OC gaseous table, 2 g cm$^{-2}$ of
water vapour), with the windows of the study (grey: the 684-990 nm windows of the tutorial) and the
simulated wavelengths.
```

## Effective path lengths

As in the tutorial, the effective transmittance of each term (path reflectance $R_{atm}$, water-leaving term,
downwelling irradiance $E_d$) is converted into an effective air mass $m^*$, and $\eta = m^*/m$ measures the
coupling ($\eta = 1$: decoupled model exact).

```{figure} figures/coupling_full_eta_350_2500nm.png
:width: 100%

$\eta$ of the three terms per channel (median and spread over the conditions), coloured by the class of
the absorber vertical distribution.
```

- The result of 684-990 nm holds over the whole range. The **path reflectance** sees a fraction of the
  geometric air mass: $\eta_{atm}$ = 0.41 (940 nm) to 0.87 (2150 nm) in median in the H$_2$O bands, 0.48
  (O$_2$-A R branch) to 0.89 (CH$_4$ 2300 nm) for the well-mixed gases, with an IQR of 0.2-0.3.
  $\eta_{atm}$ increases with wavelength beyond 1000 nm (less scattering below the absorbers).
- Two exceptions: O$_3$ (stratospheric, above the scatterers, $\eta_{atm}$ = 1.00) and the weak CO$_2$
  1600 nm band ($\eta_{atm}$ = 1.04), for which the decoupled model is already close to exact.
- The **water-leaving term and the irradiance** see 1-5 % more than the direct beam
  ($\eta$ = 1.04-1.05 below 1000 nm, 1.01-1.03 beyond 1100 nm).

## Error of the decoupled model

- **Channel level** (5 cm$^{-1}$), path reflectance: 5-16 % of the signal in median (95th percentile
  22-38 %) in the strong H$_2$O bands (720, 820, 900, 940, 1130, 1380, 1880, 2450 nm), 1.5-4.4 % in the
  O$_2$ and CO$_2$ 2010 nm bands, below 1 % for O$_3$, CO$_2$ 1600, CH$_4$ 1660 and the 1000 nm H$_2$O wing.
  Water-leaving term and irradiance: 0.1-0.6 % (95th percentile < 1.7 %).
- **Band level** (error of the absorption depth of the band): 45-48 % for OLCI Oa19 / Oa20 and MSI B9,
  29-34 % for OCI 1038 and the cirrus bands (MSI B10, OLI B9, OCI 1378), 21-29 % for Oa13, Oa15 and Oa21,
  6-13 % in the SWIR (1240-2260 nm), 3-6 % for MSI B3 and B4 (O$_3$). Weighted by the absorption depth, the
  band path reflectance is underestimated by 31-34 % (Oa20, B9, cirrus bands), 15 % (Oa13, Oa19), 4.7 %
  (Oa15) and below 2 % elsewhere; the irradiance by 0.2-3 % of the absorption, i.e. < 0.5 % of the signal.

## Photon path-length model

The coupling of each term is described by a Gamma distribution of the absorber paths of the photons, of
mean $\eta m$ and relative variance $s^2$ (`ppld` model of RadCalNet-OC, 8 coefficients per window and term,
57 models fitted in 5 min).

```{figure} figures/coupling_full_ppld_coefficients.png
:width: 100%

Fitted $\eta$ (for $m$ = 2, AOT(550) 0.1, 2 g cm$^{-2}$ of water vapour) and $s^2$ of each window and term.
```

- The model reduces the channel error of the path reflectance to 0.4-3.2 % in median (95th percentile
  1-9 %), and that of the water-leaving term and of the irradiance to < 0.07 % (95th percentile < 0.25 %).
- Band level, decoupled → `ppld`, median [95th percentile], % of the absorption: Oa13 20.9 → 2.2 [6.2],
  Oa15 21.7 → 4.8 [13.8], Oa19 47.8 → 7.9 [19.6], Oa20 44.8 → 6.2 [18.0], MSI B9 44.7 → 6.1 [17.9], cirrus
  bands 29-33 → 5.4-6.1 [19], OCI 1038 34.3 → 9.0, MSI B12 9.0 → 5.3, OCI 2130 10.2 → 6.8, MSI B3 / B4
  3.1 / 6.4 → 1.6 / 3.2; irradiance ≤ 3 → ≤ 0.5 %.
- The residual band TOA error of the path reflectance is 4.6-6.1 % in the 940 nm and cirrus bands,
  1.0-2.6 % in Oa13 / Oa15 / Oa19, < 1 % elsewhere. The gain is small where the absorption is weak and spread
  over a window (Oa21, VIIRS M8, the 1600 nm bands, all with an absorption depth < 5 %).

```{figure} figures/coupling_full_band_error_10nm.png
:width: 100%

Absorption-depth error (median and 95th percentile over the conditions) of Gaussian bands of 10 nm FWHM
across the range, decoupled and photon path-length models.
```

**Robustness.** Refitted without the water vapour column 4 g cm$^{-2}$ (resp. the solar zenith angle 70°)
and evaluated on it, the channel error of the path reflectance increases by at most 0.1 point (resp. 0.5
point, H$_2$O 2450 nm): the model extrapolates in water vapour and geometry. A single model per absorber
class is **not** sufficient (it doubles the error in the SWIR H$_2$O bands): one model per window is needed.

**O$_3$ Chappuis band.** The fitted $\eta_{atm}$ (1.15) and $s^2$ (3.4) are at the limit of the model; the
10 nm bands of 470-520 nm are worse than with the decoupled model. The O$_3$ window should keep the decoupled
model for the path reflectance.

## Absorption tables: SOS-ABS CKD against REPTRAN

The coupling corrections are relative and transfer to another code only as far as the absorption tables
agree.

```{figure} figures/coupling_full_absorption_tables.png
:width: 100%

Absorption depth $1 - T_g$ for $m = 2$ and 2 g cm$^{-2}$ of water vapour: RadCalNet-OC (REPTRAN) against the
CKD of SOS-ABS on the sampled 5 cm$^{-1}$ intervals.
```

The two agree within 10 % in the O$_3$, O$_2$-A, H$_2$O 820 / 940 / 1000, O$_2$ 1270, CO$_2$, CH$_4$ 2300 and
H$_2$O 2450 nm windows. REPTRAN is 25-40 % larger in the H$_2$O 720, 900, 1880, 2150 nm and CH$_4$ 1660 nm
windows, 70 % larger in O$_2$-B (a line-by-line O$_2$ agreed within 2 % in the 684-990 nm study), and
**one third smaller in the H$_2$O 1380 nm window** (cirrus bands). These differences exceed the coupling
residual in several windows.

## Checks

**OPAC aerosols, SOS-ABS against OSOAA** (v4 atmosphere). SOS-ABS with a converted MOPSMAP file
(continental average, RH 70 %, 0.8 µm, SZA 40°) against the OSOAA run of the same TOA LUT: aerosol term
$R(\tau = 0.1) - R(0)$ within -0.6 to +3.9 % (vector) and -1.7 to +2.5 % (scalar), total TOA reflectance
within -0.55 to +1.04 %, Rayleigh -0.35 %; same result with the angles in either order.

**End to end on Sentinel-2.** The candidate full-range table (photon path-length model in the H$_2$O 1000,
1130, 1380, 1880, 2150, 2450 nm, O$_2$ 1270 and CO$_2$ 2010 nm windows) was applied to the AERONET-OC
matchups of 30 cloud-free Sentinel-2 scenes at Venise (AAOT). In the cirrus band B10 the clear-sky L1C
reflectance is 0.50 × 10$^{-3}$ (0.55 × 10$^{-3}$ for single-scattering Rayleigh above the water vapour); the
decoupled model gives 0 and the candidate table 0.20 × 10$^{-3}$ (correlation 0.73): the 1380 nm model,
fitted on the band wings, recovers 40 % of the signal in the saturated core. The other bands are unchanged
(2190 nm: +0.015 × 10$^{-3}$). The set v5 simulates these saturated cores (88 CKD intervals of 1344-1471 and
1801-1964 nm, 16 in MSI B10).

## Recommendations and limits

- One photon path-length model per window beyond 990 nm (H$_2$O 1000, 1130, 1380, 1880, 2150, 2450 nm;
  O$_2$ 1270; CO$_2$ 2010 nm), where the gain is large and robust.
- The decoupled model for O$_3$ (at least for the path reflectance) and the CO$_2$ 1600 nm window
  ($\eta \approx 1$, error < 1 %).
- In the CH$_4$ windows the fit is made on the total transmittance: the code applying it must couple CH$_4$
  too, or drop these windows.
- Limits: maritime aerosols, US62 profiles, 1013 hPa, scalar radiative transfer, Cox-Munk 5 m s$^{-1}$,
  H$_2$O 0.5-4 g cm$^{-2}$, solar angles 20-70°, AOT(550) 0.05-0.4. The 95th percentile of the band residual
  remains 18-20 % of the absorption in the saturated H$_2$O bands (940, 1380 nm). The variants of v4 (other
  aerosols, profiles, pressure, polarisation, wind, glint) will test the sensitivity to these assumptions.
