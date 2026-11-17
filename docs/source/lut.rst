Look-up tables in NetCDF (``sosabs``)
=====================================

The Python package ``sosabs`` (directory ``python/``) runs SOS-ABS from
Python and builds look-up tables (LUT): the simulations for all the
combinations of a few parameters, gathered in one self-documented NetCDF
file. See :doc:`installation` and the tutorial
:doc:`tutorials/water_vapour_900_1000nm`.

One simulation
--------------

.. code-block:: python

   import sosabs

   ds = sosabs.run({
       "-SOS_Main.Wa": 0.94, "-ANG.Thetas": 40., "-ANG.Rad.NbGauss": 24,
       "-AP.AbsProfile.Type": 2, "-AP.H2O": 2.0, "-AP.SpectralResol": 10.,
       "-SOS.AbsModeCKD": 1, "-AP.AerProfile.Type": 1, "-AP.HR": 8.,
       "-AP.AerHS.HA": 2., "-AER.Model": 1, "-AER.WMO.Model": 2,
       "-AER.Waref": 0.55, "-AER.AOTref": 0.1,
       "-SURF.Type": 1, "-SURF.Ind": 1.33, "-SURF.Glitter.Wind": 5.,
       "-SOS.OutputAlt": 0., "-SOS.View.Dphi": 5,
   })

The parameters are given with the command-line keywords or the Python
names of :doc:`parameters`; the others take the default values of the
executable. The radiation field is always computed as a polar diagram
(``-SOS.View 2``, azimuth step ``-SOS.View.Dphi``, 5° by default).

``run`` returns an ``xarray.Dataset``:

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Variable
     - Dimensions
     - Content
   * - ``I_up``, ``Q_up``, ``U_up``
     - ``vza``, ``azi``
     - Stokes parameters of the upward field, :math:`\pi L / E_0`
   * - ``I_down``, ``Q_down``, ``U_down``
     - ``vza``, ``azi``
     - Stokes parameters of the downward field
   * - ``scat_angle_up``, ``scat_angle_down``
     - ``vza``, ``azi``
     - scattering angles (degrees)
   * - ``flux_dir_down``, ``flux_diff_down``, ``flux_tot_down``
     -
     - direct, diffuse and total downward fluxes at the surface,
       normalised to :math:`E_0 \cos\theta_s`
   * - ``flux_diff_up``
     -
     - diffuse upward flux at TOA
   * - ``coef_tronca``
     -
     - truncation coefficient of the aerosol phase function (f2py backend)
   * - ``gauss_weight``, ``user_angle`` (coordinates)
     - ``vza``
     - Gauss weights (0 for the added angles) and user-angle flags

The input parameters are stored as attributes ``sos_<keyword>`` (e.g.
``sos_AP_H2O``), with the normalisation and the azimuth convention.

Backends
~~~~~~~~

``backend="f2py"``
   in-memory call of ``SOS_PROC`` through the module ``sos``
   (``binding/compilation_f2py``): double precision, no result file. By
   default the call is made in a child process (``isolate=True``), so that
   each run starts from a fresh Fortran state and a Fortran ``STOP`` does not
   stop the calling Python session (e.g. a Jupyter kernel).
``backend="exe"``
   run of ``exe/SOS_ABS_MAIN.exe`` and reading of the ASCII files
   (6 significant digits).
``backend="auto"`` (default)
   f2py if the module can be imported, else the executable.

The console output of the Fortran code is written to ``sos_console.log`` in
the working directory of the run; its last lines are given in the error
message when a run fails.

Building a LUT
--------------

.. code-block:: python

   import numpy as np
   import sosabs

   lut = sosabs.build_lut(
       grid={"-SOS_Main.Wa": np.arange(0.900, 0.999, 0.002),
             "-AP.H2O": [0.5, 1.5, 3.0, 4.5]},
       fixed=config,                 # dict of the other parameters
       outfile="lut.nc",
       n_workers=32,                 # default: number of CPUs (max. 32)
       cache_dir="sos_cache")        # optional, reused between LUTs

The order of ``grid`` gives the order of the dimensions; any SOS-ABS
parameter can be a dimension. Usual parameters get short dimension names:

.. list-table::
   :header-rows: 1
   :widths: 35 20 45

   * - Keyword
     - Dimension
     - Units
   * - ``-SOS_Main.Wa``
     - ``wl``
     - µm
   * - ``-ANG.Thetas``
     - ``sza``
     - degree
   * - ``-AP.H2O``, ``-AP.O3``
     - ``h2o``, ``o3``
     - g cm\ :sup:`-2`, Dobson
   * - ``-AER.AOTref``
     - ``aot_ref``
     -
   * - ``-SURF.Glitter.Wind``
     - ``wind``
     - m s\ :sup:`-1`
   * - ``-AP.Psurf``
     - ``pressure``
     - hPa
   * - ``-SOS.OutputAlt``
     - ``z_out``
     - km
   * - ``-SURF.Ind``, ``-SURF.Alb``
     - ``n_surf``, ``albedo``
     -
   * - ``-SOS.AbsModeCKD``, ``-AP.SpectralResol``
     - ``ckd_mode``, ``ckd_resolution``
     - , cm\ :sup:`-1`

Other parameters keep their Python name as dimension name.

Viewing angles
   SOS-ABS adds the solar zenith angle to the output angles. When ``sza`` is
   a LUT dimension, only the angles common to all the runs (Gauss and user
   angles) are kept (``common_vza=True``).

Speed
   * all the azimuths are recombined from one reading of the Fourier series
     (``SOS_TRPHI_MULTI``, see :doc:`algorithm`);
   * the f2py backend avoids the ASCII result files;
   * the runs are parallel and share one cache of Mie and surface matrix
     files (``cache_dir/MIE``, ``cache_dir/SURF``): the code builds each file
     in the private ``TMP`` directory of the run and moves it to the cache
     (atomic rename) only if it does not exist yet, so each file is stored
     once and never read half-written. With the f2py backend each run is
     made in a new process.

   Put ``cache_dir`` on a local disk, not in a synchronised folder
   (Dropbox...): it can reach several GB and changes at each LUT.

Failed runs
   A failed run gives a warning and NaN values in the LUT; the global
   attributes ``n_failed_runs`` and ``failed_runs`` (indices) record it.

NetCDF file
~~~~~~~~~~~

``outfile`` (or ``sosabs.write_netcdf(lut, path)``) writes a compressed
NetCDF4 file (zlib level 4, Stokes parameters stored as ``float32`` by
default, ``dtype="float64"`` to keep the double precision). It contains:

* the coordinates with ``units``, ``long_name`` and ``sos_keyword``;
* the variables of :func:`sosabs.run` with the LUT dimensions in front;
* global attributes: ``Conventions`` (CF-1.8), ``title``, ``history``,
  ``date_created``, ``lut_dimensions``, ``lut_keywords``,
  ``sos_code_version`` (version and git commit of the code), ``sos_backend``,
  ``normalisation``, ``azimuth_convention`` and all the fixed parameters
  (``sos_<keyword>``).

.. code-block:: python

   import xarray as xr

   lut = xr.open_dataset("lut.nc")
   lut.interp(h2o=2.2, vza=40.).sel(azi=135.).I_down   # sky radiance

Conventions
~~~~~~~~~~~

* Normalised radiance :math:`X = \pi L / E_0`, with :math:`E_0` the
  extraterrestrial solar irradiance normal to the beam; the reflectance is
  :math:`X / \cos\theta_s`.
* ``vza`` is the zenith angle of the line of sight: from nadir for the
  upward field, from zenith for the downward field.
* ``azi`` is the relative azimuth between the line of sight and the Sun:
  0° when the sensor looks towards the Sun's azimuth (sun glint in the
  upward field, solar aureole in the downward field), 180° when it looks
  away from the Sun (backscattering).
