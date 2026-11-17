Constants (SOS.h)
=================

``inc/SOS.h`` defines array sizes, default values, data file names and
numerical thresholds shared by all routines and by the Python binding.
They are compile-time constants: **recompile after any change**.
Changing them requires a good knowledge of the numerical method; the
most useful ones are listed below.

Common
------

.. list-table::
   :header-rows: 1
   :widths: 40 15 45

   * - Constant
     - Value
     - Meaning
   * - ``CTE_LENDIR`` / ``CTE_LENFIC1`` / ``CTE_LENFIC2``
     - 350 / 150 / 500
     - length of directory names, file names, complete file names
   * - ``CTE_WAMIN``, ``CTE_WAMAX``
     - 0.364, 4.0
     - spectral range (µm)
   * - ``CTE_NOT_DEFINED_VALUE_INT`` / ``_DBLE``
     - −999 / −999.0
     - values of undefined parameters

Aerosols
--------

.. list-table::
   :header-rows: 1
   :widths: 40 15 45

   * - Constant
     - Value
     - Meaning
   * - ``CTE_MIE_DIM``
     - 10000
     - size of the Mie arrays (> 2 × max size parameter + 20)
   * - ``CTE_MAXNB_ANG_EXT``
     - 3600
     - max number of angles of an external phase matrix (376 for the
       MOPSMAP/OPAC files, see :doc:`tutorials/opac_aerosol_luts`)
   * - ``CTE_DEFAULT_AER_JUNGE_RMAX``
     - 50
     - default maximum radius of the Junge law (µm)
   * - ``CTE_AER_MU1_TRONCA``, ``CTE_AER_MU2_TRONCA``
     - 0.8, 0.94
     - cosines bounding the linearisation of the truncated phase
       function
   * - ``CTE_PH_SEUIL_TRONCA``
     - 0.1
     - the truncation is cancelled below this truncation coefficient
   * - ``CTE_MAX_NB_MODE_MIXTURE``
     - 20
     - max number of modes of a user mixture
   * - ``CTE_AER_DATAWMO``, ``CTE_AER_DATASF``, ``CTE_AER_*_SF``
     - files of ``fic/``
     - WMO and Shettle & Fenn data

Atmospheric profile and absorption
----------------------------------

.. list-table::
   :header-rows: 1
   :widths: 40 15 45

   * - Constant
     - Value
     - Meaning
   * - ``CTE_HT_STD_PSURF``
     - 1013
     - standard pressure (hPa)
   * - ``CTE_TOA_ALT``
     - 120
     - altitude of TOA (km)
   * - ``CTE_OS_NT``, ``CTE_OS_NT_MIN``
     - 600, 100
     - max and min number of layers
   * - ``CTE_TCOUCHE``
     - 0.005
     - optical thickness of a layer
   * - ``CTE_ABS_NBLEV``
     - 50
     - number of levels of the gas profiles
   * - ``CTE_ABS_H2O``, ``_CO2``, ``_O3``, ``_N2O``, ``_CO``, ``_CH4``,
       ``_O2``, ``CTE_CONT_NO2``
     - 1
     - gases taken into account (0 to remove)
   * - ``CTE_CKD_NUMIN``, ``CTE_CKD_NUMAX``
     - 2500, 27500
     - range of the CKD files (cm\ :sup:`-1`)
   * - ``CTE_TAUABS_MAX``
     - 999
     - max absorption optical thickness

Surface
-------

.. list-table::
   :header-rows: 1
   :widths: 40 15 45

   * - Constant
     - Value
     - Meaning
   * - ``CTE_TETAS_LIM_ROUJEAN``, ``CTE_TETAV_LIM_ROUJEAN``
     - 60, 60
     - max solar and viewing zenith angles for the Roujean BRDF (degrees)
   * - ``CTE_PH_NU``
     - 1024
     - azimuths for the wave slope probability

Radiative transfer
------------------

.. list-table::
   :header-rows: 1
   :widths: 40 15 45

   * - Constant
     - Value
     - Meaning
   * - ``CTE_MDF``
     - 0.0279
     - molecular depolarisation factor
   * - ``CTE_DEFAULT_IGMAX``
     - 100
     - default max order of scattering
   * - ``CTE_PH_SEUIL_SUMDIF``, ``CTE_PH_SEUIL_SF``
     - 1e-5, 1e-4
     - convergence thresholds of the orders of scattering and of the
       Fourier series
   * - ``CTE_SOLAR_DISC_SOLID_ANGLE``
     - 6.8e-5
     - solid angle of the solar disc (sr)
   * - ``CTE_REMOVE_PREVIOUS_BIN_FILE``
     - 1
     - overwrite an existing ``-SOS.ResBin`` (0: stop)

Angles and expansion orders
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 40 15 45

   * - Constant
     - Value
     - Meaning
   * - ``CTE_DEFAULT_NBMU_LUM``
     - 24
     - default number of Gauss angles for the radiance
   * - ``CTE_DEFAULT_NBMU_MIE``
     - 40
     - default number of Gauss angles for the phase matrix
   * - ``CTE_OS_NBMU_MAX``
     - 80
     - max angles for the radiance (≥ NbGauss + user angles)
   * - ``CTE_MIE_NBMU_MAX``
     - 100
     - max angles for the phase matrix
   * - ``CTE_NBMAX_USER_ANGLES``
     - 20
     - max number of user angles
   * - ``CTE_OS_NB_MAX``, ``CTE_OS_NS_MAX``, ``CTE_OS_NM_MAX``
     - 200, 136, 336
     - max orders of the expansions (see the advice in ``SOS.h``)

.. tip::

   To use more Gauss angles than ``CTE_OS_NBMU_MAX`` allows, increase
   ``CTE_OS_NBMU_MAX``, ``CTE_NBANGLES_MAX`` and the expansion orders
   consistently (``CTE_OS_NS_MAX = 2 × CTE_OS_NBMU_MAX``,
   ``CTE_OS_NM_MAX = CTE_OS_NB_MAX + CTE_OS_NS_MAX``), then recompile.
   Larger arrays slow down the code.
