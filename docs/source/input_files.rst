Input files
===========

User angles
-----------

``-ANG.Rad.UserAngFile`` lists viewing zenith angles (degrees) at which
the radiance is computed in addition to the Gauss angles, one value per
line, without blank line at the end (at most ``CTE_NBMAX_USER_ANGLES`` =
20 angles):

.. code-block:: text

   20
   25
   30
   35
   40

The radiance files ``-SOS.ResFileUp`` / ``-SOS.ResFileDown`` contain
both Gauss and user angles; ``-SOS.ResFileUp.UserAng`` /
``-SOS.ResFileDown.UserAng`` contain only the user angles.

``-ANG.Aer.UserAngFile`` lists scattering angles (≤ 90°) added to the
Gauss angles of the phase matrix computation; the code adds the
complementary angles 180° − θ. They only appear in the aerosol log file
(``-AER.Log``).

.. important::

   Mie and surface files are named after the user angle *file name*.
   Use a new file name when you change the list of angles (see
   :doc:`running`).

Examples: ``binding/TestFiles/ficAngRad.txt`` and ``ficAngAer.txt``.

User gas profile
----------------

With ``-AP.AbsProfile.Type 0``, ``-AP.AbsProfile.UserFile`` gives the
gas profile on ``CTE_ABS_NBLEV`` = 50 levels, one line per level, after
one header line, with 14 columns:

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - Column
     - Content
   * - 0
     - level index (not used)
   * - 1
     - altitude (km)
   * - 2
     - pressure (hPa)
   * - 3
     - temperature (K)
   * - 4 to 10
     - H\ :sub:`2`\ O, CO\ :sub:`2`, O\ :sub:`3`, N\ :sub:`2`\ O, CO,
       CH\ :sub:`4`, O\ :sub:`2` concentrations (ppmv)
   * - 11
     - air density (not used)
   * - 12
     - NO\ :sub:`2` concentration (ppmv)
   * - 13
     - SO\ :sub:`2` concentration (ppmv, no longer used)

.. code-block:: text

   LEVEL ALT(km) P(hPa) T(K)  + GASES (ppmv), except column 12: air density
   01   0.0 1013.00000 299.7 0.25930E+05 330. 0.28690E-01 0.32000E+00 0.15000E+00 0.17000E+01 0.20900E+06 0.24500E+20 0.23000E-04 0.30000E-03
   02   1.0  904.00000 293.7 0.19490E+05 330. 0.31500E-01 0.32000E+00 0.14500E+00 0.17000E+01 0.20900E+06 0.22310E+20 0.23000E-04 0.27400E-03
   ...
   50 120.0    0.00002 380.0 0.20000E+00  35. 0.50000E-03 0.18510E-03 0.50000E+02 0.30000E-01 0.72500E+05 0.42250E+12 0.17000E-03 0.35800E-06

Aerosol properties from an external source
------------------------------------------

With ``-AER.Model 4``, ``-AER.ExtData`` gives the phase matrix and
cross sections, in a format derived from the DLS tool of O. Dubovik
(spheroids). The values of the first three lines are read after the
``:``.

.. code-block:: text

   EXTINCTION_COEF : 0.19941E+00
   SCATTERING_COEF : 0.19695E+00
   NB_LINES : 163
   Angle    F11        -F12/F11    F22/F11    F33/F11
   180.00   1.810E-01   0.000E+00  1.000E+00  -1.000E+00
   179.14   1.810E-01  -4.686E-07  1.000E+00  -1.000E+00
   ...
     0.00   7.109E+00   0.000E+00  1.000E+00   1.000E+00

- line 1: extinction cross section; line 2: scattering cross section;
- line 3: number of angles (≤ ``CTE_MAXNB_ANG_EXT`` = 200);
- line 4: comment;
- then, for each scattering angle (degrees): :math:`F_{11}`,
  :math:`-F_{12}/F_{11}`, :math:`F_{22}/F_{11}`, :math:`F_{33}/F_{11}`.
  :math:`F_{12}` is negative for Rayleigh scattering, as in SOS-ABS.

Following lines are not read. ``-AER.Waref`` must equal
``-SOS_Main.Wa``.

Pre-computed aerosol properties
-------------------------------

``-AER.UserFile`` skips the aerosol computation. It has the format of the
``-AER.ResFile`` output (see :doc:`output_files`): it can be a file
produced by a previous run or a user file with Greek coefficients. Only
the Greek coefficients, the truncation coefficient and the single
scattering albedo are used; the other values are informative.

Mixture of aerosol modes
------------------------

With ``-AER.Model 5``, ``-AER.DefMixture`` defines up to
``CTE_MAX_NB_MODE_MIXTURE`` = 20 modes. Values are read after the
``:`` separator; the text before it is free (but must not contain
``:``).

.. code-block:: text

   NUMBER OF AEROSOLS MODES: 3
   SIZE DISTRIBUTION MODE 1: LND
      MODAL RADIUS (microns): 0.2
      STANDARD DEVIATION (microns): 0.1
      REFRACTIVE INDEX at WA_SIMU - Real part: 1.20
      REFRACTIVE INDEX at WA_SIMU - Imaginary part (<0): -0.0120
      REFRACTIVE INDEX at WA_REF - Real part: 1.4
      REFRACTIVE INDEX at WA_REF - Imaginary part (<0): -0.014
      PROPORTION TO THE TOTAL AOT at WA_REF: 0.5
   SIZE DISTRIBUTION MODE 2: LND
      ...
      PROPORTION TO THE TOTAL AOT at WA_REF: 0.1
   SIZE DISTRIBUTION MODE 3: JUNGE
      SLOPE: 4.00
      MIN RADIUS (microns): 0.01
      MAX RADIUS (microns): 50.0
      REFRACTIVE INDEX at WA_SIMU - Real part: 1.4
      REFRACTIVE INDEX at WA_SIMU - Imaginary part (<0): -0.001
      REFRACTIVE INDEX at WA_REF - Real part: 1.45
      REFRACTIVE INDEX at WA_REF - Imaginary part (<0): -0.002
      PROPORTION TO THE TOTAL AOT at WA_REF: 0.4

For each mode, in this order:

- ``LND`` (modal radius, standard deviation) or ``JUNGE`` (slope,
  minimum radius, maximum radius), in capital letters;
- refractive index at the simulation wavelength (real, imaginary < 0);
- refractive index at the reference wavelength, required even if it is
  equal to the simulation wavelength;
- proportion of the mode in the total AOT at the reference wavelength.
  The proportions must sum to 1 (tolerance ``CTE_GAP_TOLER_SUM_RATES``).
