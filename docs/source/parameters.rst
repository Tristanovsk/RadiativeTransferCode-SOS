Simulation parameters
=====================

Each parameter is given on the command line as ``-Keyword value``. The
second column gives the name of the corresponding argument of
``sos.sos_proc`` in the :doc:`python_binding`. Status: **R** required,
**D** default value from ``SOS.h``, **O** optional, **C** conditional
(see :doc:`running`). File names "without path" are written in the
``SOS`` or ``LOG`` sub-directory of the working directory.

General parameters
------------------

.. list-table::
   :header-rows: 1
   :widths: 25 18 6 51

   * - Keyword
     - Python
     - Status
     - Definition
   * - ``-SOS_Main.ResRoot``
     - ``resroot``
     - R
     - Working directory (complete path).
   * - ``-SOS_Main.Wa``
     - ``wa_simu``
     - R
     - Wavelength of the radiance simulation (µm), between 0.364 and 4.0.
   * - ``-ANG.Thetas``
     - ``tetas``
     - R
     - Solar zenith angle (degrees, 0 < θs < 90).
   * - ``-SOS_Main.Log``
     - ``ficmain_log``
     - O
     - Main log file (without path).

Atmospheric profile
-------------------

The profile is not discretised in altitude but in optical thickness,
from TOA (120 km) to the surface. The molecular optical thickness
decreases exponentially with altitude:

.. math::

   \tau_{mol}(z) = \tau_{mol}^{surf} \exp(-z / h_{mol})

where :math:`\tau_{mol}^{surf}` is given by ``-AP.MOT`` or computed from
the surface pressure :math:`P` (``-AP.Psurf``, hPa) and the wavelength
:math:`\lambda` (µm):

.. math::

   \tau_{mol}^{surf} = \frac{P}{P_0}
   \left( \frac{84.35}{\lambda^4} - \frac{1.225}{\lambda^5}
   + \frac{1.4}{\lambda^6} \right) 10^{-4}

with :math:`P_0` = 1013 hPa (``CTE_HT_STD_PSURF``). With
``-AP.AerProfile.Type 1``, the aerosol optical thickness follows the
same law with the scale height :math:`h_{aer}`:

.. math::

   \tau_{aer}(\lambda, z) = \tau_{aer}^{surf}(\lambda) \exp(-z / h_{aer}),
   \qquad
   \tau_{aer}^{surf}(\lambda) = \tau_{aer}^{surf}(\lambda_{ref})
   \frac{\tilde\sigma_{ext}(\lambda)}{\tilde\sigma_{ext}(\lambda_{ref})}

where :math:`\tau_{aer}^{surf}(\lambda_{ref})` is ``-AER.AOTref`` at
``-AER.Waref`` and :math:`\tilde\sigma_{ext}` the extinction cross
section of the aerosol model. With ``-AP.AerProfile.Type 2``, aerosols
are homogeneously mixed with molecules between ``-AP.AerLayer.Zmin`` and
``-AP.AerLayer.Zmax``, with pure molecular layers above and below (this
case is not available with gaseous absorption).

.. list-table::
   :header-rows: 1
   :widths: 25 18 6 51

   * - Keyword
     - Python
     - Status
     - Definition
   * - ``-AP.Log``
     - ``ficprofil_log``
     - O
     - Log file of the profile computation.
   * - ``-AP.AerProfile.Type``
     - ``iprofil``
     - R
     - 1: exponential profiles (scale heights); 2: aerosol layer between
       two molecular layers (no gaseous absorption).
   * - ``-AP.Psurf``
     - ``psurf``
     - D
     - Surface pressure (hPa). Default ``CTE_HT_STD_PSURF`` = 1013.
   * - ``-AP.MOT``
     - ``tr``
     - O
     - Molecular optical thickness at the simulation wavelength. If not
       given, computed from ``-AP.Psurf``.
   * - ``-AP.HR``
     - ``hr``
     - C
     - Scale height of the molecular profile (km), typically 8. Required
       if the molecular optical thickness is not null.
   * - ``-AER.AOTref``
     - ``aot_ref``
     - R
     - Aerosol optical thickness at the reference wavelength.
   * - ``-AER.Waref``
     - ``waref_aot``
     - C
     - Reference wavelength of ``-AER.AOTref`` (µm). Required if
       ``-AER.AOTref`` > 0.
   * - ``-AP.AerHS.HA``
     - ``ha``
     - C
     - Scale height of the aerosol profile (km). Required for
       ``-AP.AerProfile.Type 1`` with aerosols.
   * - ``-AP.AerLayer.Zmin``, ``-AP.AerLayer.Zmax``
     - ``zmin``, ``zmax``
     - C
     - Bottom and top of the aerosol layer (km), for
       ``-AP.AerProfile.Type 2``.

Gaseous absorption
------------------

Gases taken into account: H\ :sub:`2`\ O, CO\ :sub:`2`, O\ :sub:`3`,
N\ :sub:`2`\ O, CO, CH\ :sub:`4`, O\ :sub:`2` (lines and continuum) and
NO\ :sub:`2` (continuum only). A gas can be removed by setting its
``CTE_ABS_*`` constant to 0 in ``SOS.h`` and recompiling.

.. list-table:: Standard gas profiles (``-AP.AbsProfile.Type``) and their default amounts
   :header-rows: 1
   :widths: 6 34 15 15 15 15

   * - Type
     - Profile
     - H\ :sub:`2`\ O (g cm\ :sup:`-2`)
     - O\ :sub:`3` (DU)
     - CO\ :sub:`2` (ppmv)
     - CH\ :sub:`4` (ppmv)
   * - 0
     - user profile (``-AP.AbsProfile.UserFile``)
     -
     -
     -
     -
   * - 1
     - Tropical
     - 4.068
     - 280.95
     - 330
     - 1.7
   * - 2
     - Mid-latitude summer
     - 2.906
     - 333.12
     - 330
     - 1.7
   * - 3
     - Mid-latitude winter
     - 0.853
     - 376.85
     - 330
     - 1.7
   * - 4
     - Sub-arctic summer
     - 2.080
     - 346.57
     - 330
     - 1.7
   * - 5
     - Sub-arctic winter
     - 0.417
     - 374.78
     - 330
     - 1.7
   * - 6
     - US Standard 1962
     - 1.417
     - 343.03
     - 330
     - 1.7
   * - 7
     - no absorption
     - 0
     - 0
     - 0
     - 0

The total columns of O\ :sub:`3` and H\ :sub:`2`\ O and the surface
concentrations of CO\ :sub:`2` and CH\ :sub:`4` can be changed with
``-AP.O3``, ``-AP.H2O``, ``-AP.CO2`` and ``-AP.CH4``: the standard
profile is then scaled to these values.

.. list-table::
   :header-rows: 1
   :widths: 25 18 6 51

   * - Keyword
     - Python
     - Status
     - Definition
   * - ``-AP.AbsProfile.Type``
     - ``absprofil``
     - R
     - Gas profile, 0 to 7 (table above).
   * - ``-AP.AbsProfile.UserFile``
     - ``ficabsprofil``
     - C
     - User gas profile (complete path), for type 0. Format in
       :doc:`input_files`.
   * - ``-AP.SpectralResol``
     - ``nustep``
     - C
     - Spectral resolution of the CKD coefficients: 1, 5 or
       10 cm\ :sup:`-1`. Required with absorption.
   * - ``-SOS.AbsModeCKD``
     - ``imode_ckd_calcul``
     - C
     - 1: fine modelling of the absorption/scattering coupling;
       2: simplified modelling (see :doc:`running`).
   * - ``-AP.H2O``
     - ``h2o``
     - O
     - Water vapour column (g cm\ :sup:`-2`).
   * - ``-AP.O3``
     - ``o3``
     - O
     - Ozone column (Dobson units).
   * - ``-AP.CO2``
     - ``co2``
     - O
     - CO\ :sub:`2` concentration at the surface (ppmv).
   * - ``-AP.CH4``
     - ``ch4``
     - O
     - CH\ :sub:`4` concentration at the surface (ppmv).

Aerosols
--------

These parameters are needed when ``-AER.AOTref`` is not null.

.. list-table::
   :header-rows: 1
   :widths: 25 18 6 51

   * - Keyword
     - Python
     - Status
     - Definition
   * - ``-AER.Model``
     - ``imod_aer``
     - C
     - 0: mono-modal; 1: WMO; 2: Shettle & Fenn; 3: bimodal log-normal;
       4: phase matrix from an external source; 5: user mixture of modes.
   * - ``-AER.DirMie``
     - ``dir_mie``
     - C
     - Directory of the Mie files (complete path).
   * - ``-AER.Tronca``
     - ``itronc_aer``
     - D
     - Truncation of the forward peak of the phase function: 0 no, 1 yes
       (default).
   * - ``-AER.ResFile``
     - ``ficgranu``
     - D
     - Output file of the aerosol radiative properties (without path).
       Default ``Aerosols.txt``.
   * - ``-AER.UserFile``
     - ``ficuser_aer``
     - O
     - Pre-computed aerosol properties (complete path), same format as
       ``-AER.ResFile``. ``-AER.Waref`` must equal ``-SOS_Main.Wa``.
   * - ``-AER.Log``
     - ``ficgranu_log``
     - O
     - Log file of the aerosol computations.
   * - ``-AER.MieLog``
     - ``ficmie_log``
     - O
     - Log file of the Mie computations.

Mono-modal size distribution (``-AER.Model 0``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Log-normal distribution (``-AER.MMD.SDtype 1``):

.. math::

   N(r) = \frac{1}{r \sigma \sqrt{2\pi}}
   \exp\left( -\frac{\ln^2(r / r_m)}{2\sigma^2} \right)

Junge distribution (``-AER.MMD.SDtype 2``): :math:`N(r) = r_{min}^{-\nu}`
for :math:`r < r_{min}`, :math:`N(r) = r^{-\nu}` for
:math:`r_{min} \le r \le r_{max}` and 0 above :math:`r_{max}`. The slope
:math:`\nu = 3` is a singular value.

.. list-table::
   :header-rows: 1
   :widths: 25 18 6 51

   * - Keyword
     - Python
     - Status
     - Definition
   * - ``-AER.MMD.SDtype``
     - ``igranu``
     - C
     - 1: log-normal; 2: Junge.
   * - ``-AER.MMD.MRwa``, ``-AER.MMD.MIwa``
     - ``rn_wa``, ``in_wa``
     - C
     - Real and imaginary (negative) parts of the refractive index at the
       simulation wavelength.
   * - ``-AER.MMD.MRwaref``, ``-AER.MMD.MIwaref``
     - ``rn_waref``, ``in_waref``
     - C
     - Refractive index at the reference wavelength, if
       ``-AER.Waref`` ≠ ``-SOS_Main.Wa``.
   * - ``-AER.MMD.LNDradius``
     - ``lnd_radius_mmd_aer``
     - C
     - Modal radius :math:`r_m` (µm).
   * - ``-AER.MMD.LNDvar``
     - ``lnd_lnvar_mmd_aer``
     - C
     - Standard deviation :math:`\sigma`.
   * - ``-AER.MMD.JD.slope``
     - ``jd_slope_mmd_aer``
     - C
     - Junge slope :math:`\nu`.
   * - ``-AER.MMD.JD.rmin``
     - ``jd_rmin_mmd_aer``
     - C
     - Minimum radius (µm).
   * - ``-AER.MMD.JD.rmax``
     - ``jd_rmax_mmd_aer``
     - D
     - Maximum radius (µm). Default and recommended: 50.

WMO models (``-AER.Model 1``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Mixtures of the WMO components dust-like, water soluble, oceanic and
soot [WMO1986]_.

.. list-table::
   :header-rows: 1
   :widths: 25 18 6 51

   * - Keyword
     - Python
     - Status
     - Definition
   * - ``-AER.WMO.Model``
     - ``imodele_wmo``
     - C
     - 1: continental; 2: maritime; 3: urban; 4: user mixture.
   * - ``-AER.WMO.DL``, ``-AER.WMO.WS``, ``-AER.WMO.OC``, ``-AER.WMO.SO``
     - ``c_wmo_dl``, ``c_wmo_ws``, ``c_wmo_oc``, ``c_wmo_so``
     - C
     - Volume concentrations (0 to 1, sum = 1) of the dust-like, water
       soluble, oceanic and soot components, for ``-AER.WMO.Model 4``.

Shettle & Fenn models (``-AER.Model 2``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Size distributions and refractive indices depend on the relative
humidity [ShettleFenn1979]_.

.. list-table::
   :header-rows: 1
   :widths: 25 18 6 51

   * - Keyword
     - Python
     - Status
     - Definition
   * - ``-AER.SF.Model``
     - ``imodele_sf``
     - C
     - 1: tropospheric; 2: urban; 3: maritime; 4: coastal.
   * - ``-AER.SF.RH``
     - ``rh``
     - C
     - Relative humidity (0 to 99 %).

Bimodal log-normal (``-AER.Model 3``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fine (FM) and coarse (CM) log-normal modes:

.. math::

   \frac{dV(r)}{d\ln r} = \sum_{i=1}^{2} \frac{C_{V,i}}{\sigma_i \sqrt{2\pi}}
   \exp\left( -\frac{\ln^2(r / r_i)}{2\sigma_i^2} \right)

The mixture is given either by the volume concentrations
:math:`C_{V,i}` (``-AER.BMD.VCdef 1``) or by the ratio of the coarse mode
AOT to the total AOT at the reference wavelength (``-AER.BMD.VCdef 2``).

.. list-table::
   :header-rows: 1
   :widths: 25 18 6 51

   * - Keyword
     - Python
     - Status
     - Definition
   * - ``-AER.BMD.VCdef``
     - ``mode_param_bilnd``
     - C
     - 1: volume concentrations; 2: AOT ratio.
   * - ``-AER.BMD.CoarseVC``, ``-AER.BMD.FineVC``
     - ``user_cv_coarse``, ``user_cv_fine``
     - C
     - Volume concentrations of the coarse and fine modes (VCdef = 1).
   * - ``-AER.BMD.RAOT``
     - ``rtauct_waref``
     - C
     - AOT\ :sub:`coarse` / AOT\ :sub:`total` at the reference wavelength
       (VCdef = 2).
   * - ``-AER.BMD.CM.SDradius``, ``-AER.BMD.CM.SDvar``
     - ``bmd_cm_rmodal``, ``bmd_cm_var``
     - C
     - Modal radius (µm) and standard deviation of the coarse mode.
   * - ``-AER.BMD.CM.MRwa``, ``-AER.BMD.CM.MIwa``
     - ``bmd_cm_mrwa``, ``bmd_cm_miwa``
     - C
     - Refractive index of the coarse mode at the simulation wavelength.
   * - ``-AER.BMD.CM.MRwaref``, ``-AER.BMD.CM.MIwaref``
     - ``bmd_cm_mrwaref``, ``bmd_cm_miwaref``
     - C
     - Refractive index of the coarse mode at the reference wavelength
       (if different).
   * - ``-AER.BMD.FM.SDradius``, ``-AER.BMD.FM.SDvar``
     - ``bmd_fm_rmodal``, ``bmd_fm_var``
     - C
     - Modal radius (µm) and standard deviation of the fine mode.
   * - ``-AER.BMD.FM.MRwa``, ``-AER.BMD.FM.MIwa``
     - ``bmd_fm_mrwa``, ``bmd_fm_miwa``
     - C
     - Refractive index of the fine mode at the simulation wavelength.
   * - ``-AER.BMD.FM.MRwaref``, ``-AER.BMD.FM.MIwaref``
     - ``bmd_fm_mrwaref``, ``bmd_fm_miwaref``
     - C
     - Refractive index of the fine mode at the reference wavelength
       (if different).

External data and user mixtures (``-AER.Model 4`` and ``5``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 18 6 51

   * - Keyword
     - Python
     - Status
     - Definition
   * - ``-AER.ExtData``
     - ``ficextdata_aer``
     - C
     - Phase matrix and cross sections from an external source, e.g.
       non-spherical particles [Dubovik2006]_ (complete path), for model 4.
       ``-AER.Waref`` must equal ``-SOS_Main.Wa``.
   * - ``-AER.DefMixture``
     - ``ficmixture_aer``
     - C
     - Definition of a mixture of log-normal and Junge modes (complete
       path), for model 5; suited to CAMS aerosols [Bozzo2017]_.

The file formats are described in :doc:`input_files`.

Surface
-------

.. list-table:: Surface types (``-SURF.Type``)
   :header-rows: 1
   :widths: 10 90

   * - Type
     - Model
   * - 0
     - Lambertian surface of albedo ``-SURF.Alb``
   * - 1
     - Rough sea, sun glint (Cox and Munk) [CoxMunk1954]_
   * - 2
     - Flat sea (Fresnel reflection)
   * - 3
     - Roujean BRDF [Roujean1992]_
   * - 4
     - Roujean BRDF + Rondeaux BPDF (vegetation) [Rondeaux1991]_
   * - 5
     - Roujean BRDF + Bréon BPDF (bare soils) [Breon1995]_
   * - 6
     - Roujean BRDF + Nadal BPDF (no longer available)
   * - 7
     - Roujean BRDF + Maignan BPDF [Maignan2009]_

The Lambertian term ``-SURF.Alb`` is added to every surface type (for
the sea, it represents the foam and water-leaving reflectance).

The Roujean BRDF is

.. math::

   \rho(\theta_1, \theta_2, \phi) = k_0 + k_1 \rho_1(\theta_1, \theta_2, \phi)
   + k_2 \rho_2(\theta_1, \theta_2, \phi)

.. warning::

   Roujean et al. (1992) give :math:`k_0, k_1, k_2` in percent of
   reflectance. SOS-ABS expects reflectances, i.e. values 100 times
   smaller: for a ploughed field (k0 = 24.3, k1 = 7.3, k2 = 64.2), use
   ``-SURF.Roujean.K0 0.243 -SURF.Roujean.K1 0.073 -SURF.Roujean.K2 0.642``.

.. list-table::
   :header-rows: 1
   :widths: 25 18 6 51

   * - Keyword
     - Python
     - Status
     - Definition
   * - ``-SURF.Type``
     - ``isurf``
     - R
     - Surface type (table above).
   * - ``-SURF.Alb``
     - ``rho``
     - R
     - Lambertian albedo.
   * - ``-SURF.Dir``
     - ``dir_surf``
     - C
     - Directory of the BRDF/BPDF matrix files (complete path), for
       types other than 0.
   * - ``-SURF.File``
     - ``ficsurf``
     - D
     - Pre-computed BRDF/BPDF matrix (complete path). By default (and
       recommended), the matrix is computed by the code.
   * - ``-SURF.Ind``
     - ``surf_ind``
     - C
     - Refractive index of the surface relative to air, at the
       simulation wavelength (types 1, 2, 4, 5, 7).
   * - ``-SURF.Glitter.Wind``
     - ``wind``
     - C
     - Wind speed at the sea surface (m s\ :sup:`-1`), type 1.
   * - ``-SURF.Roujean.K0``, ``K1``, ``K2``
     - ``k0_roujean``, ``k1_roujean``, ``k2_roujean``
     - C
     - Roujean parameters (reflectance units), types ≥ 3.
   * - ``-SURF.Maignan.C``
     - ``coef_c_maignan``
     - C
     - Coefficient :math:`C \exp(-NDVI)` of the Maignan BPDF, type 7.
   * - ``-SURF.Nadal.Alpha``, ``-SURF.Nadal.Beta``
     - ``alpha_nadal``, ``beta_nadal``
     - C
     - Nadal BPDF parameters (type 6, no longer available).
   * - ``-SURF.Log``
     - ``ficsurf_log``
     - O
     - Log file of the surface computations.

Angles
------

The radiance field and the BRDF/BPDF matrices are computed on Gauss
angles, to which the user can add angles. The phase functions are
computed on a separate set of Gauss angles. These numbers set the limit
orders of the series expansions.

.. list-table::
   :header-rows: 1
   :widths: 25 18 6 51

   * - Keyword
     - Python
     - Status
     - Definition
   * - ``-ANG.Rad.NbGauss``
     - ``nbmu_gauss_lum``
     - D
     - Number of Gauss angles for the radiance. Default 24.
   * - ``-ANG.Rad.UserAngFile``
     - ``ficangles_user_lum``
     - O
     - List of user viewing zenith angles added to the Gauss angles
       (complete path).
   * - ``-ANG.Rad.ResFile``
     - ``ficangles_res_lum``
     - D
     - Output file of the angles used for the radiance (without path).
       Default ``SOS_UsedAngles.txt``.
   * - ``-ANG.Aer.NbGauss``
     - ``nbmu_gauss_mie``
     - D
     - Number of Gauss angles for the phase matrix. Default 40.
   * - ``-ANG.Aer.UserAngFile``
     - ``ficangles_user_mie``
     - O
     - List of user scattering angles (≤ 90°, complete path).
   * - ``-ANG.Aer.ResFile``
     - ``ficangles_res_mie``
     - D
     - Output file of the angles used for the phase matrix (without
       path). Default ``Aer_UsedAngles.txt``.
   * - ``-ANG.Log``
     - ``ficanglog``
     - O
     - Log file of the angle computations.

Outputs
-------

.. list-table::
   :header-rows: 1
   :widths: 25 18 6 51

   * - Keyword
     - Python
     - Status
     - Definition
   * - ``-SOS.View``
     - ``itrphi``
     - R
     - 1: radiance in one azimuth plane; 2: polar diagram.
   * - ``-SOS.View.Phi``
     - ``phios``
     - C
     - Relative azimuth (degrees), for ``-SOS.View 1``.
   * - ``-SOS.View.Dphi``
     - ``pas_phi``
     - C
     - Azimuth step (integer degrees), for ``-SOS.View 2``.
   * - ``-SOS.OutputAlt``
     - ``zout``
     - D
     - Output altitude (km, 0 to 120). If not given, upward radiance at
       TOA and downward radiance at the surface.
   * - ``-SOS.ResFileUp``, ``-SOS.ResFileDown``
     - (not in Python)
     - D
     - Upward and downward radiance files (without path). Default
       ``SOS_Up.txt``, ``SOS_Down.txt``.
   * - ``-SOS.ResFileUp.UserAng``, ``-SOS.ResFileDown.UserAng``
     - (not in Python)
     - O
     - Same, restricted to the user angles.
   * - ``-SOS.ResBin``
     - ``ficsos_res_bin``
     - D
     - Binary file of (I, Q, U) as Fourier series (without path).
       Default ``SOS_Result.bin``.
   * - ``-SOS.Trans``
     - ``fictrans``
     - O
     - Direct and diffuse transmittances (time consuming).
   * - ``-SOS.Flux``
     - ``ficflux``
     - O
     - Downward fluxes at the surface and upward flux at TOA.
   * - ``-SOS.Log``
     - ``ficsos_log``
     - O
     - Log file of the successive orders computation.

Expert parameters
~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 18 6 51

   * - Keyword
     - Python
     - Status
     - Definition
   * - ``-SOS.IGmax``
     - ``igmax``
     - D
     - Maximum order of scattering and surface reflection. Default 100.
   * - ``-SOS.Ipolar``
     - ``ipolar``
     - D
     - 0: no polarisation. If not given, polarisation is computed
       (recommended).
   * - ``-SOS.DirectSurf``
     - ``idirsurf``
     - D
     - Direct reflection of the solar beam by the surface (sun glint,
       Fresnel flat sea, BRDF/BPDF), added analytically to the diffuse
       field. 1 (default): complete field. 0: diffuse field only, the
       direct surface terms being added separately (as in the OSOAA LUTs,
       switch ``CTE_LUT_DIRECT_SURF``). Upward field only; the downward
       field and the fluxes are not changed.
