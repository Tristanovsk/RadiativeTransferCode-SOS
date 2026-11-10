"""sosabs: Python tools for the SOS-ABS radiative transfer code.

* :func:`run` runs one simulation and returns an :class:`xarray.Dataset`;
* :func:`build_lut` runs a parameter grid in parallel and writes a NetCDF
  look-up table with its metadata;
* :mod:`sosabs.ckd` converts the CKD absorption coefficients
  (``fic/COEFF_CKD``) to NetCDF and computes band transmittances.
"""

__version__ = "0.1.0"

from .params import PARAMETERS, KEYWORD, NAME, DEFAULTS  # noqa: F401
from .runner import run, sos_root, available_backend  # noqa: F401
from .lut import build_lut, write_netcdf  # noqa: F401
from . import ckd  # noqa: F401
