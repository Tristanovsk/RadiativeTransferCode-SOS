Output files
============

Results are written in the ``SOS`` sub-directory, and log files in the
``LOG`` sub-directory, of the working directory ``-SOS_Main.ResRoot``.

Angles used
-----------

``-ANG.Rad.ResFile`` (default ``SOS_UsedAngles.txt``) lists the angles
of the radiance computation and the limit orders of the expansions:

.. code-block:: text

   NB_TOTAL_ANGLES : 41
   NB_GAUSS_ANGLES : 40
   ANGLES_USERFILE :NO_USER_ANGLES
   SOLAR ZENITH ANGLE : 35.000
   INTERNAL_IMUS : 16
   INTERNAL_OS_NB : 80
   INTERNAL_OS_NS : 80
   INTERNAL_OS_NM : 160
   INDEX  COS_ANGLE              WEIGHT                USER_ANGLE
      1   0.99955382265163D+00   0.11449500031887D-02     0
      ...

``INTERNAL_IMUS`` is the index of the solar angle; ``OS_NB`` is the
limit order of the Legendre expansion of the phase functions, ``OS_NS``
that of the Fresnel matrix and of the Fourier series of the radiance,
``OS_NM`` that of the Fourier series of the wave slope distribution.
The last column flags user angles (1).

``-ANG.Aer.ResFile`` (default ``Aer_UsedAngles.txt``) gives the same
information (without solar angle, ``OS_NS`` and ``OS_NM``) for the phase
matrix computation.

Aerosol properties
------------------

``-AER.ResFile`` (default ``Aerosols.txt``):

.. code-block:: text

   EXTINCTION CROSS SECTION (mic^2)    : 0.55281E-03
   SCATTERING CROSS SECTION (mic^2)    : 0.35911E-03
   ASYMMETRY FACTOR (no truncation)    : 0.59091E+00
   TRUNCATION COEFFICIENT              : 0.00000
   SINGLE SCATTERING ALBEDO (truncation): 0.64961
   ---------------------------------
   PHASE MATRIX COEFFICIENTS FOR K=0 TO 80
     ALPHA(K)        BETA11(K)       GAMMA12(K)      ZETA(K)
    0.00000000E+00  0.10000000E+01  0.00000000E+00  0.00000000E+00
    ...

The cross sections are in µm². The single scattering albedo and the
Greek coefficients :math:`\alpha_k, \beta_k, \gamma_k, \zeta_k` of the
Legendre expansion of the phase matrix (:math:`\beta_k` for the
intensity, the others for polarisation), for :math:`k = 0` to
``OS_NB``, are adjusted to the truncation of the phase function.

Radiance files
--------------

``-SOS.ResFileUp`` (default ``SOS_Up.txt``) and ``-SOS.ResFileDown``
(default ``SOS_Down.txt``) give the upward and downward radiance
fields, by default at TOA and at the surface, or at the altitude
``-SOS.OutputAlt`` (linear interpolation in the profile). The Stokes
parameters are normalised to the extraterrestrial solar irradiance:

.. math::

   I = \frac{\pi L(z)}{E_{sun}} \quad (\mathrm{sr}^{-1})

and similarly for Q, U and the polarised intensity
:math:`I_{pol} = \sqrt{Q^2 + U^2}`.

``-SOS.View 1`` (one azimuth plane, ``-SOS.View.Phi``): 31 header lines
(starting with ``#``), then the columns

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Column
     - Content
   * - ``VZA``
     - viewing zenith angle (degrees), signed (see below)
   * - ``SCA_ANG``
     - scattering angle (degrees)
   * - ``I``, ``Q``, ``U``
     - Stokes parameters (sr\ :sup:`-1`)
   * - ``POL_ANG``
     - polarisation angle (degrees; −999 if undefined)
   * - ``POL_RATE``
     - degree of polarisation (%)
   * - ``IPOL``
     - polarised intensity (sr\ :sup:`-1`)

with the Fortran format ``2(2X,F7.2),2X,3(E13.6,2X),2(F7.2,2X),E13.6``.
The negative angles correspond to the azimuth ``Phi + 180°``, the
positive angles to ``Phi``.

``-SOS.View 2`` (polar diagram, step ``-SOS.View.Dphi``): 26 header
lines, then the same columns preceded by the relative azimuth ``PHI``
(degrees), format ``3(2X,F7.2),1X,3(E13.6,2X),2(1X,F7.2),E13.6``.

Reading the files in Python:

.. code-block:: python

   import numpy as np

   names = ['vza', 'sca_ang', 'I', 'Q', 'U', 'pol_ang', 'pol_rate', 'Ipol']
   data = np.loadtxt('SOS/SOS_Up.txt', comments='#')
   up = dict(zip(names, data.T))          # -SOS.View 1
   # -SOS.View 2: names = ['phi'] + names

Angle conventions
-----------------

- The relative azimuth is 0° when the sun and the observer are in
  opposite half-planes with respect to the zenith, and 180° when they
  are in the same half-plane (opposite to the geographical convention).
- The viewing zenith angle is positive in the half-plane of azimuth
  :math:`\phi` and negative in the half-plane :math:`\phi + 180°`. In the
  solar principal plane, the backscattering direction (sun behind the
  observer) has a negative VZA.
- For the upward radiance, VZA = 0 is a propagation towards the zenith
  (a sensor looking at nadir); for the downward radiance, VZA = 0 is a
  propagation towards the nadir (a sensor looking at the zenith).
  VZA = ±90° is the horizon.
- Internally, the solar direction is a downward direction and is
  redefined as 180° − θ\ :sub:`s`.

Polarisation angle
------------------

With :math:`\vec\ell` in the meridian plane and :math:`\vec r` such that
:math:`(\vec\ell, \vec r, \vec k)` is direct,
:math:`Q = I_p \cos 2\chi` and :math:`U = I_p \sin 2\chi`, with
:math:`-\pi/2 < \chi \le \pi/2`:

- :math:`Q > 0`: :math:`\chi = \frac{1}{2}\arctan(U/Q)`;
- :math:`Q < 0`: :math:`\chi = \pm\frac{\pi}{2} + \frac{1}{2}\arctan(U/Q)`
  (+ if :math:`U > 0`, − otherwise);
- :math:`Q = 0`: :math:`\chi = \pm\pi/4` according to the sign of U,
  undefined if :math:`U = 0`.

Binary file
-----------

``-SOS.ResBin`` (default ``SOS_Result.bin``) stores (I, Q, U) as Fourier
series in azimuth. With ``CTE_REMOVE_PREVIOUS_BIN_FILE = 1`` (default),
an existing file is overwritten.

Transmittances
--------------

``-SOS.Trans`` (optional, time consuming) gives the direct and diffuse
transmittances from TOA to the surface for the solar angle, and the
diffuse transmittance from the surface to TOA for each viewing angle, for
the true atmosphere (not the one adjusted to the truncation):

.. code-block:: text

   Solar Zenith Angle : 35.000
   Direct transmission TOA -> surface :  0.6176
   Diffuse transmittance : TOA -> surface
       thetas = 35.000   td(thetas) = 0.1711
   Diffuse transmittance : surface -> TOA
       thetav =  1.712   td(thetav) = 0.1532
       ...

Fluxes
------

``-SOS.Flux`` (optional) gives the direct, diffuse and total downward
fluxes at the surface and the upward diffuse flux at TOA, normalised to
the solar flux at TOA (so the direct flux is the direct transmittance),
followed by the profile of molecular (MOT), aerosol (AOT), gas (GOT) and
total optical thickness:

.. code-block:: text

   Solar Zenith Angle : 40.000
    Downward fluxes at BOA (normalized by TOA solar flux)
      - Downward direct flux at BOA :  0.1562
      - Downward diffuse flux at BOA:  0.0582
      ==> Downward total flux at BOA:  0.2144
    Upward diffuse flux at TOA (normalized by TOA solar flux):  2.1411568E-002
    According to the following profile
    Z(km)    MOT     AOT     GOT    TOTAL
    120.00  0.0000  0.0000 -0.0000  0.0000
    ...
      0.00  0.0085  0.3821  1.0317  1.4222

The direct transmission and flux are only valid without gaseous
absorption or with ``-SOS.AbsModeCKD 2``.

Log files
---------

All log files are optional; the Mie, surface and SOS logs are large and
slow to produce.

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Keyword
     - Content
   * - ``-SOS_Main.Log``
     - execution of the routines: parameters, warnings, errors
   * - ``-ANG.Log``
     - Gauss angles for phase functions and radiance, solar and user
       angles
   * - ``-AP.Log``
     - atmospheric profile
   * - ``-AER.Log``
     - aerosol properties: phase functions of the components and of the
       mixture, truncation, cross sections
   * - ``-AER.MieLog``
     - phase matrices versus size parameter (Mie)
   * - ``-SURF.Log``
     - BRDF/BPDF matrices: Fresnel matrices, wave slope distribution
   * - ``-SOS.Log``
     - successive orders of scattering computation
