Installation
============

Repository layout
-----------------

.. code-block:: text

   RadiativeTransferCode-SOS/
   ├── src/       Fortran 77 sources (fixed form, .F, C preprocessor)
   ├── inc/       SOS.h: array dimensions, default values, thresholds
   ├── gen/       Makefile.gfortran
   ├── obj/       object files (created by the Makefile)
   ├── exe/       SOS_ABS_MAIN.exe and the demo scripts
   ├── fic/       tabulated data: WMO and Shettle & Fenn aerosol models,
   │              SO2/NO2 profile, CKD coefficients (COEFF_CKD/)
   ├── binding/   f2py compilation script and Python example (run_sos.py)
   ├── doc/       user manual (PDF)
   └── docs/      this documentation

Source files
------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - File
     - Role
   * - ``SOS_ABS_MAIN.F``
     - Main program: reads the ``-Keyword value`` arguments of the
       command line and calls ``SOS_PROC``.
   * - ``SOS_PROC.F``
     - Checks the simulation parameters and runs the sequence of
       modules. This is the routine exposed to Python by the binding.
   * - ``SOS_ANGLES.F``
     - Gauss and user angles, limit orders of the series expansions.
   * - ``SOS_AEROSOLS.F``
     - Aerosol radiative properties (phase matrix, cross sections).
   * - ``SOS_MIE.F``
     - Mie computations versus the size parameter.
   * - ``SOS_SURFACE.F``, ``SOS_NOM_FIC_SURFACE.F``
     - Surface BRDF/BPDF reflection matrices and their file names.
   * - ``SOS_ROUJEAN.F``, ``SOS_GLITTER.F``, ``SOS_SURFACE_BPDF.F``
     - Roujean BRDF, sun glint on a rough sea, Rondeaux / Bréon /
       Maignan BPDF.
   * - ``SOS_PREPA_ABSPROFILE.F``
     - CKD coefficients and gas concentration profile (50 levels).
   * - ``SOS_SUB_TRS.F``
     - CKD absorption coefficients for a given pressure and temperature.
   * - ``SOS_ABSPROFILE.F``
     - Gaseous optical thickness profile for one term of the CKD
       expansion.
   * - ``SOS_PREPA_OS.F``
     - Prepares the parameters of a single radiative transfer
       computation.
   * - ``SOS_PROFIL.F``
     - Optical thickness profile combining molecular and aerosol
       scattering and gaseous absorption, at SOS resolution.
   * - ``SOS.F``
     - Single radiative transfer computation; adjusts the profile to the
       truncation of the aerosol phase function.
   * - ``SOS_OS.F``
     - Core of the successive orders of scattering; results as a Fourier
       series in azimuth.
   * - ``SOS_AGGREGATE.F``
     - Weights and sums the computations of each CKD term.
   * - ``SOS_TRPHI.F``
     - Recombines the radiance for given azimuth angles.

Requirements
------------

- ``gfortran`` (the reference compilation used GCC 10.2.0; GCC 15.2 also
  works)
- ``ksh`` or ``bash`` to run the demo scripts
- for the Python binding: Python 3 with ``numpy`` (which provides
  ``f2py``); tested with Python 3.10 (user manual) and Python 3.9 /
  numpy 1.21. Recent numpy versions (Python ≥ 3.12) build ``f2py``
  modules with ``meson`` and may need ``pip install meson ninja``.

Environment variable
--------------------

``SOS_ABS_ROOT`` must point to the root of the repository. It is used by
the Makefile, by the demo scripts and by the code itself to locate the
data of ``$SOS_ABS_ROOT/fic``. Define it in your shell configuration
file (``.bashrc``, ``.profile``, ...):

.. code-block:: bash

   export SOS_ABS_ROOT=/path/to/RadiativeTransferCode-SOS

Compiling the executable
------------------------

.. code-block:: bash

   cd $SOS_ABS_ROOT/gen
   make -f Makefile.gfortran

This compiles the sources with ``gfortran -O -ffpe-summary=none -DSYSV
-I$SOS_ABS_ROOT/inc``, writes the objects in ``$SOS_ABS_ROOT/obj`` and
links ``$SOS_ABS_ROOT/exe/SOS_ABS_MAIN.exe``. A copy of ``inc/SOS.h`` is
saved in ``exe/config.txt`` to keep track of the constants used for the
compilation.

.. note::

   The constants of ``inc/SOS.h`` (array sizes, default values, gases
   taken into account, ...) are fixed at compile time: recompile after
   any change (see :doc:`constants`). For debugging, uncomment
   ``DEBUG = -g`` in the Makefile and remove the ``-O`` option.

Compiling the Python binding
----------------------------

The script ``binding/compilation_f2py`` compiles the sources with
``-O -fPIC``, preprocesses ``SOS_PROC.F`` and builds the ``sos`` extension
module with ``f2py``. Run it from the ``binding`` directory (it creates the
``outdir`` sub-directory):

.. code-block:: bash

   cd $SOS_ABS_ROOT/binding
   . ./compilation_f2py

.. warning::

   Do not compile the binding with ``-O3``: with gfortran 15 the radiances
   returned by ``sos.sos_proc`` are then NaN (the executable is not
   affected).

The result is a shared library ``sos.<python tag>.so`` in ``binding``;
add this directory to ``PYTHONPATH`` to ``import sos``. See
:doc:`python_binding`.

Python package ``sosabs``
-------------------------

The package ``python/sosabs`` (runs returned as ``xarray`` datasets, NetCDF
look-up tables, see :doc:`lut`) is installed with:

.. code-block:: bash

   pip install -e $SOS_ABS_ROOT/python

It needs ``numpy``, ``xarray`` and ``netCDF4``, and the executable and/or the
f2py module above.

String length limits
--------------------

Directory and file names are stored in fixed-length strings defined in
``inc/SOS.h``: ``CTE_LENDIR = 350`` characters for directories,
``CTE_LENFIC1 = 150`` for file names without path and
``CTE_LENFIC2 = 500`` for complete file names.
