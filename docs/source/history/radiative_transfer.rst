A short history of radiative transfer
=====================================

Polarised light
---------------

In 1852, George Gabriel Stokes described any beam of light, polarised or not,
by four measurable quantities, :math:`I, Q, U, V`: the total intensity, two
linear polarisation components and the circular component [Stokes1852]_.
These **Stokes parameters** add for incoherent beams, which makes them the
natural unknowns of radiative transfer. The sky light was already known to be
polarised (Arago, 1809), and the neutral points of the sky polarisation
discovered by Arago, Babinet and Brewster were one of the first puzzles that a
theory of multiple scattering had to explain.

In 1871, Lord Rayleigh showed that scattering by particles much smaller than
the wavelength varies as :math:`\lambda^{-4}`, which explains the blue sky,
and that the light scattered at 90° is completely polarised. In 1908, Gustav
Mie gave the exact solution for the scattering of a plane wave by a
homogeneous sphere of any size [Mie1908]_, still used to compute the phase
matrix of aerosols and hydrosols.

The transfer equation
---------------------

The equation of radiative transfer was written at the beginning of the 20th
century for stellar atmospheres: Arthur Schuster (1905) studied the radiation
through a "foggy" atmosphere with two streams, and Karl Schwarzschild (1906)
the radiative equilibrium of the solar atmosphere [Schuster1905]_.
Subrahmanyan Chandrasekhar gave the theory its modern form in *Radiative
Transfer* (1950) [Chandrasekhar1950]_: he wrote the transfer equation for the
full Stokes vector and solved exactly the problem of a Rayleigh scattering
atmosphere, including the polarisation of the sky light and its neutral
points.

Numerical methods
-----------------

With computers, the transfer equation could be solved for realistic
atmospheres with aerosols. Several methods were developed in the 1960s and
1970s and are still used: the discrete ordinates method, the adding-doubling
method, the Monte Carlo method and the **successive orders of scattering**
method, in which the radiation field is built as the sum of the photons
scattered once, twice, and so on. Hansen and Travis (1974) reviewed these
methods and the polarisation of planetary atmospheres [HansenTravis1974]_;
Lenoble (1985) gathered the standard computational procedures
[Lenoble1985]_.

The successive orders of scattering method is simple, it gives the radiance at
any level and its physical meaning is clear. Its convergence is slow for
optically thick media (clouds), but it is well suited to clear atmospheres
and to the ocean, as developed at the Laboratoire d'Optique Atmosphérique
(see :doc:`sos_loa`).

Remote sensing
--------------

From the 1970s, satellite remote sensing of the ocean colour (CZCS, 1978) and
of aerosols required accurate radiative transfer codes to correct the
atmospheric signal, including the reflection on the rough sea surface
described by Cox and Munk (1954) [CoxMunk1954]_. The POLDER instrument (1996)
and the PARASOL mission (2004-2013) measured the polarised reflectance from
space, which made the vector (polarised) codes indispensable. Most of them
solve the transfer equation for :math:`(I, Q, U)`: the circular polarisation
:math:`V` is very small for the Earth atmosphere and is usually neglected (see
:doc:`../roadmap/stokes_v`).

See :doc:`/references` for the full references.
