Python binding
==============

The routine ``SOS_PROC`` is exposed to Python by ``f2py`` as
``sos.sos_proc``. Build the module as described in :doc:`installation`
(``binding/compilation_f2py``); a complete example is given in
``binding/run_sos.py`` (edit ``REP_RES`` and ``REP_DATA`` to your
directories).

Calling ``sos_proc``
--------------------

All arguments are keyword arguments. Their names and meanings are given
in the *Python* column of :doc:`parameters`. Every argument must be
given: unused ones are set to the "not defined" values of ``SOS.h``,
``-999`` (integers) and ``-999.0`` (floats), or to the placeholder
strings used by the code (``"NO_USER_ANGLES"``, ``"NO_LOG_FILE"``,
``"NO_OUTPUT"``, ...).

Two arguments have no command-line equivalent:

``ier``
   status, set to 0 on input; the routine sets it to −1 on error.
``trace``
   ``True`` to write the main log file (``ficmain_log``).

.. code-block:: python

   import sos

   (nblum, ind_angout, phi, vza,
    sca_ang_up, i_up, q_up, u_up, pol_ang_up, pol_rate_up, l_pol_up,
    sca_ang_down, i_down, q_down, u_down, pol_ang_down, pol_rate_down, l_pol_down,
    flux_dir_down, flux_diff_down, flux_tot_down, flux_diff_up,
    coef_tronca) = sos.sos_proc(resroot=..., wa_simu=0.440, tetas=40.0,
                                ...,  # all other arguments
                                ier=0, trace=True)

``run_sos.py`` gathers the defaults in a dictionary indexed by the
command-line keywords (``dict_sos_all_parameters``), updates it with the
user values (``dict_sos_user``) and converts it into the argument list
with ``set_sos_params``. ``gen_sos_output`` writes ``SOS_Up.txt`` and
``SOS_Down.txt`` files with the same format as the executable.

Outputs
-------

The radiance arrays are indexed ``[azimuth, angle]`` with the relative
azimuth from 0 to 360° by steps of 1° and the angle index from 0 to
``CTE_OS_NBMU_MAX``; only the first ``nblum`` angles are meaningful.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Output
     - Content
   * - ``nblum``
     - number of angles in the output arrays
   * - ``ind_angout``
     - 1 for the angles of the user file ``ficangles_user_lum``, else 0
   * - ``phi``
     - relative azimuths (degrees), 0 to 360 by 1°
   * - ``vza``
     - viewing zenith angles (degrees)
   * - ``sca_ang_up``, ``sca_ang_down``
     - scattering angles (degrees), 2D
   * - ``i_up``, ``q_up``, ``u_up``, ``i_down``, ``q_down``, ``u_down``
     - Stokes parameters (sr\ :sup:`-1`, normalised as in
       :doc:`output_files`), 2D
   * - ``pol_ang_up``, ``pol_ang_down``
     - polarisation angle (degrees)
   * - ``pol_rate_up``, ``pol_rate_down``
     - degree of polarisation (%)
   * - ``l_pol_up``, ``l_pol_down``
     - polarised intensity (sr\ :sup:`-1`), 2D
   * - ``flux_dir_down``
     - direct downward flux at the surface = direct transmittance (true
       atmosphere; valid without absorption or with
       ``imode_ckd_calcul = 2``)
   * - ``flux_diff_down``, ``flux_tot_down``
     - diffuse and total downward fluxes at the surface
   * - ``flux_diff_up``
     - upward diffuse flux at TOA
   * - ``coef_tronca``
     - truncation coefficient of the aerosol phase function

All fluxes are normalised to the solar flux at TOA.
