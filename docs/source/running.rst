Running simulations
===================

Command line
------------

A simulation is run by calling ``SOS_ABS_MAIN.exe`` with a list of
``-Keyword value`` pairs:

.. code-block:: bash

   $SOS_ABS_ROOT/exe/SOS_ABS_MAIN.exe -Keyword1 value1 -Keyword2 value2 ...

Only the parameters needed by the options selected for the simulation
have to be given. A parameter that is not used by these options is
ignored. If a needed parameter is missing, the run stops and the name of
the missing parameter is printed. The keywords are listed in
:doc:`parameters`.

Each parameter has one of the following statuses:

``R`` (required)
   must be given.
``D`` (default)
   a value is needed; if it is not given, the default of ``inc/SOS.h``
   is used.
``O`` (optional)
   may be given.
``C`` (conditional)
   required by the value of another parameter.

Demo scripts
------------

``exe/runSOS-ABS_demo.ksh`` simulates the radiance over the sea at
910 nm (water vapour band):

- surface pressure 1013 hPa, exponential profiles of molecules and
  aerosols with scale heights of 8 and 2 km,
- WMO maritime aerosols, AOT of 0.3 at 550 nm,
- tropical gas profile, CKD coefficients at 10 cm\ :sup:`-1`, fine
  modelling of the absorption/scattering coupling,
- rough sea (wind 2 m s\ :sup:`-1`), refractive index 1.34, null albedo,
- solar zenith angle of 35°, output in the solar principal plane with
  40 Gauss angles.

.. code-block:: bash

   export dirRESULT=${SOS_ABS_ROOT}/DEMO_RESULTS
   export dirMIE=$dirRESULT/MIE
   export dirSURF=$dirRESULT/BRDF_BPDF

   $SOS_ABS_ROOT/exe/SOS_ABS_MAIN.exe -SOS_Main.Wa 0.910 \
       -SOS_Main.ResRoot ${dirRESULT} -SOS_Main.Log SOS_Main_Demo.Log \
       -ANG.Rad.NbGauss 40 -ANG.Thetas 35. -SOS.View 1 -SOS.View.Phi 0. \
       -SOS.ResFileUp SOS_Up_Demo.txt -SOS.ResFileDown SOS_Down_Demo.txt \
       -AP.Log Profile_Demo.Log \
       -AP.Psurf 1013 -AP.AerProfile.Type 1 -AP.HR 8.0 -AP.AerHS.HA 2.0 \
       -AP.AbsProfile.Type 1 -AP.SpectralResol 10. -SOS.AbsModeCKD 1 \
       -AER.DirMie ${dirMIE} -AER.Model 1 -AER.WMO.Model 2 \
       -AER.Waref 0.550 -AER.AOTref 0.3 \
       -AER.ResFile Aerosols_Demo.txt -AER.Log Aerosols_Demo.Log \
       -SURF.Dir ${dirSURF} -SURF.Type 1 -SURF.Alb 0.00 -SURF.Ind 1.34 \
       -SURF.Glitter.Wind 2.0

``exe/runSOS-ABS_demoPolar.ksh`` runs the same case with
``-SOS.View 2 -SOS.View.Dphi 30``, i.e. a polar diagram of the radiance
for relative azimuths from 0 to 360° by steps of 30°.

The demo takes about 15 s on a recent laptop.

.. warning::

   In a command file, do not leave any blank after the ``\`` at the end
   of a line.

Working directory
-----------------

All outputs are written below the working directory given by
``-SOS_Main.ResRoot``. Its sub-directories are created automatically:

.. code-block:: text

   $ResRoot/
   ├── SOS/   result files: radiances (-SOS.ResFileUp / -SOS.ResFileDown),
   │          binary Fourier series (-SOS.ResBin), aerosol properties
   │          (-AER.ResFile), angles used (-ANG.Rad.ResFile,
   │          -ANG.Aer.ResFile), transmissions and fluxes if requested
   ├── LOG/   log files, if requested
   └── TMP/   temporary files

File names given for outputs are therefore names *without* directory
(see :doc:`output_files`).

Mie and surface databases
-------------------------

The Mie computations (``-AER.DirMie``) and the surface reflection
matrices (``-SURF.Dir``) are stored in files whose names encode their
parameters, for instance:

.. code-block:: text

   MIE/MIE1.381-0.00000-0.0001-00800.00-MU40     # nr, ni, alpha min/max, Gauss angles
   BRDF_BPDF/GLITTER/GLITTER-1.340-02.0-MU40     # index, wind, Gauss angles

These files are computed once and reused by all later simulations with
the same parameters, which saves most of the computation time. Point
``-AER.DirMie`` and ``-SURF.Dir`` to a shared database directory rather
than to a per-simulation directory.

.. important::

   The names of these files include the name of the user angle files
   (``-ANG.Rad.UserAngFile``, ``-ANG.Aer.UserAngFile``), not their
   content. Change the file name when you change the list of angles,
   otherwise the old Mie and interface files are reused.

.. note::

   Concurrent runs that need the same Mie or surface file would all try
   to compute it. For a new set of parameters, run one job first, then
   launch the batch.

Gaseous absorption: fine or simplified CKD
------------------------------------------

The CKD method expands the gas transmission over a spectral interval as
a weighted sum of exponentials. Two options are available
(``-SOS.AbsModeCKD``):

``1`` -- fine modelling
   a full radiative transfer computation is made for each term of the
   CKD expansion and the radiances are summed with the CKD weights. It
   accounts for the absorption/scattering coupling and is recommended in
   strong absorption bands, but its cost grows with the number of terms.
``2`` -- simplified modelling
   the CKD method only gives the absorption optical thickness profile,
   then a single radiative transfer computation is made. Faster, and
   close to the fine method in moderate absorption bands, but not suited
   to strong bands.

The direct transmission and the direct downward flux (``-SOS.Trans``,
``-SOS.Flux``) are only valid without gaseous absorption or with
``-SOS.AbsModeCKD 2``.
