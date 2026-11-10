The SOS code, from Lille to SOS-ABS
===================================

The Lille school
----------------

The successive orders of scattering code was developed at the Laboratoire
d'Optique Atmosphérique (LOA, University of Lille), where Jean-Luc Deuzé
developed the method during his thesis [Deuze1974]_, in the group of
Maurice Herman and Jacqueline Lenoble.

The rough ocean
---------------

In 1989, Deuzé, Herman and Santer extended the code to an atmosphere bounded
by a rough ocean [Deuze1989]_. Until then, the Fourier series expansion of the
radiation field in azimuth, which makes the problem tractable, was only
compatible with a Lambertian ground or a flat sea. Their idea was to expand
separately

- the Fresnel reflection matrix, as a series of generalised Legendre
  functions of the reflection angle, like a phase matrix (only once, the sea
  water index varying little), and
- the Cox and Munk wave slope distribution, as a Fourier series of the
  azimuth (depending only on the wind speed),

the reflection matrix being the product of the two series. The sunglint of
the direct solar beam, which would need thousands of terms, is computed
exactly and replaces its truncated series. The rough ocean then costs only a
few percent more than a black surface. The paper keeps the full Stokes vector
:math:`(I, Q, U, V)` in its equations; the code, as most of its successors,
solves for :math:`(I, Q, U)` only (see :doc:`../roadmap/stokes_v`).

POLDER and PARASOL
------------------

The code was used for the analysis of the POLDER and PARASOL polarised
measurements, and was published as a reference code in 2007
[Lenoble2007]_. It was then maintained and extended
with CNES by CS (now CS GROUP France), B. Lafrance and X. Lenot:
land BRDF (Roujean) and BPDF models (Rondeaux and Herman, Bréon, Nadal,
Maignan), WMO and Shettle & Fenn aerosol models, bimodal and user aerosols,
user angles, transmissions and fluxes.

The OS method was also extended to the coupled atmosphere-ocean system,
including the radiation field below the sea surface, by Chami et al. with the
OSOA code [Chami2001]_, which became OSOAA [Chami2015]_.

SOS-ABS
-------

SOS-ABS (version 5.1, 2023) adds the gaseous absorption coupled with
scattering. The approach was established with the LOA (P. Dubuisson,
J. Riedi); the CKD coefficients were computed by HYGEOS (M. Compiègne) for
H\ :sub:`2`\ O, CO\ :sub:`2`, O\ :sub:`3`, N\ :sub:`2`\ O, CO, CH\ :sub:`4`,
O\ :sub:`2` and NO\ :sub:`2` from 2500 to 27500 cm\ :sup:`-1`, and the
coupling was implemented by CS GROUP France. The code also provides a Python
binding and multi-mode aerosol mixtures compatible with CAMS. It is
distributed by CNES under the GPL v3 licence:
https://github.com/CNES/RadiativeTransferCode-SOS.

See :doc:`/references` for the full references.
