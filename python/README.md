# sosabs

Python tools for the SOS-ABS radiative transfer code:

- `sosabs.run(params)`: one simulation, returned as an `xarray.Dataset`
  (Stokes parameters I, Q, U of the upward and downward fields versus the
  viewing zenith angle and the relative azimuth, fluxes, input parameters
  as attributes);
- `sosabs.build_lut(grid, fixed, outfile=...)`: parallel simulations on a
  parameter grid, gathered in a NetCDF look-up table with CF-style
  attributes.

Installation (from the repository root):

```bash
export SOS_ABS_ROOT=$PWD
make -f gen/Makefile.gfortran                     # executable (exe backend)
(cd binding && mkdir -p outdir && bash compilation_f2py)   # f2py backend
pip install -e python
```

See the documentation (`docs/`) and the tutorial notebook
`docs/source/tutorials/water_vapour_900_1000nm.ipynb`.
