CKD absorption coefficients
===========================

Gaseous absorption is described by the correlated-k distribution (CKD)
method (see :doc:`algorithm`): over each spectral interval, the
transmission of a gas layer with a column amount :math:`u` (molecules
cm\ :sup:`-2`) is

.. math::

   T(u) = \sum_{i=1}^{n} a_i \exp\left(-k_i(P, T)\, u\right)

with at most :math:`n = 5` terms. For H\ :sub:`2`\ O, :math:`k_i` also
depends on the H\ :sub:`2`\ O mole fraction.

ASCII files
-----------

The coefficients are stored in ``fic/COEFF_CKD/<res>cmm1`` for the spectral
resolutions 1, 5 and 10 cm\ :sup:`-1` (``-AP.SpectralResol``), one file
per gas and 500 cm\ :sup:`-1` chunk between 2500 and 27500 cm\ :sup:`-1`:
``coef_<GAS>_<numax>_<numin>_<res>cmm1``. Gases: H\ :sub:`2`\ O,
CO\ :sub:`2`, O\ :sub:`3`, N\ :sub:`2`\ O, CO, CH\ :sub:`4`, O\ :sub:`2`,
NO\ :sub:`2` (N\ :sub:`2` files are provided but not used).

Each file has a header describing its content, the grids of temperature
(9 values, 160 to 320 K), pressure (31 values, 0.007 to 1099 hPa) and, for
H\ :sub:`2`\ O, mole fraction (12 values), followed by one block per
interval: wavelengths and wavenumbers of the interval and number of terms
:math:`n`, the weights :math:`a_i`, then the :math:`k_i` tables. A block
with :math:`n = 0` means no absorption. These files are read by the
routine ``READ_CKD_COEFF`` (``src/SOS_SUB_TRS.F``).

.. note::

   On the interval lines, the larger wavenumber is given first.

NetCDF version
--------------

For analyses in Python, :mod:`sosabs.ckd` gathers all the files of one
resolution into a single NetCDF file, ``fic/COEFF_CKD/ckd_<res>cmm1.nc``
(not versioned; created on demand):

.. code-block:: python

   import sosabs

   sosabs.ckd.to_netcdf(10)      # about 30 s
   ds = sosabs.ckd.open_ckd(10)

   # weights and coefficients of O2 around 762 nm
   ds.a.sel(gas="O2").sel(wavenumber=13105, method="nearest")
   ds.k.sel(gas="O2").sel(wavenumber=13105, method="nearest")

   # band transmission of a homogeneous layer
   t = sosabs.ckd.transmittance(ds, "O2", amount=4.5e24,
                                pressure=500., temperature=250.)
   t.plot(x="wavelength")

Content of the dataset:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Variable
     - Content
   * - ``n_term(gas, wavenumber)``
     - number of terms (0: no absorption)
   * - ``a(gas, wavenumber, term)``
     - weights :math:`a_i` (NaN for unused terms)
   * - ``k(gas, wavenumber, term, pressure, temperature)``
     - :math:`k_i(P, T)` in cm\ :sup:`2` molecule\ :sup:`-1`, all gases
       but H\ :sub:`2`\ O
   * - ``k_h2o(wavenumber, term, h2o_fraction, h2o_pressure, h2o_temperature)``
     - :math:`k_i(C, P, T)` of H\ :sub:`2`\ O

Coordinates give the central wavenumber, the interval bounds and the
corresponding wavelengths (µm).

For the 10 cm\ :sup:`-1` resolution, the 450 ASCII files (309 MB) become a
59 MB file. Values are stored in single precision, which loses nothing:
the ASCII values are themselves single-precision numbers (the smallest one,
1.4·10\ :sup:`-45`, is the float32 denormal minimum); the H\ :sub:`2`\ O
files, written with 6 digits, are reproduced to 6·10\ :sup:`-8`. Selecting
a gas band takes about 0.1 s against about 0.4 s to parse the ASCII files
of the gas.

.. note::

   :func:`sosabs.ckd.transmittance` interpolates linearly in pressure,
   temperature and mole fraction, whereas SOS-ABS uses a spline in
   temperature: values can differ slightly from those of the Fortran code.
