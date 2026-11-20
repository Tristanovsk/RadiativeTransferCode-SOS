SOS-ABS documentation
=====================

**SOS-ABS** (Successive Orders of Scattering method coupled with the gaseous
ABSorptions) is a 1D plane-parallel vector radiative transfer code for the
Earth surface-atmosphere system, from 0.364 to 4 µm. It can be run from the
command line or from Python.

.. include:: _readme.md
   :parser: myst_parser.sphinx_

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   running
   python_binding
   lut
   algorithm
   parameters
   input_files
   ckd
   output_files

.. toctree::
   :maxdepth: 1
   :caption: Tutorials

   tutorials/water_vapour_900_1000nm
   tutorials/gaseous_coupling
   tutorials/gaseous_coupling_full_range
   tutorials/opac_aerosol_luts

.. toctree::
   :maxdepth: 2
   :caption: History

   history/radiative_transfer
   history/sos_loa

.. toctree::
   :maxdepth: 2
   :caption: Roadmap

   roadmap/stokes_v

.. toctree::
   :maxdepth: 1
   :caption: Reference

   constants
   references
