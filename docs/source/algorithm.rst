Algorithm
=========

Successive orders of scattering
-------------------------------

The code solves the vector radiative transfer equation for the Stokes vector
:math:`(I, Q, U)` in a plane-parallel atmosphere above a reflecting surface.
The radiation field is the sum of the contributions of the photons scattered
(or reflected by the surface) once, twice, ..., :math:`n` times:

.. math::
   :label: sos_sum

   \mathbf{I}(\tau, \mu, \phi) = \sum_{n \ge 1} \mathbf{I}^{(n)}(\tau, \mu, \phi)

where :math:`\mathbf{I}^{(n)}` is computed from :math:`\mathbf{I}^{(n-1)}`
with the integral form of the transfer equation, :math:`\mathbf{I}^{(1)}`
being the single scattering of the direct solar beam. A reflection on the
surface counts as one order. The series is stopped when the contribution of
an order is negligible (``CTE_PH_SEUIL_SUMDIF``), or replaced by a geometric
series once the ratio of two successive orders is stable
(``CTE_PH_SEUIL_CV_SG``), or at the order ``-SOS.IGmax``.

Fourier series in azimuth
-------------------------

The radiation field is expanded in Fourier series of the relative azimuth
[Deuze1989]_:

.. math::
   :label: fourier

   \begin{pmatrix} I \\ Q \\ U \end{pmatrix}(\tau, \mu, \phi)
   = \sum_{s=0}^{S} (2 - \delta_{0s})
   \begin{pmatrix} I^s(\tau, \mu) \cos s\phi \\
                   Q^s(\tau, \mu) \cos s\phi \\
                   U^s(\tau, \mu) \sin s\phi \end{pmatrix}

so that the transfer equation splits into independent equations for each
term :math:`s`, solved on Gauss angles. This requires the expansion of the
phase matrix in generalised spherical functions (Greek coefficients
:math:`\alpha_l, \beta_l, \gamma_l, \zeta_l`, see :doc:`output_files`) and of
the surface reflection matrices in Fourier series. For a rough sea, Deuzé et
al. (1989) obtained it by expanding separately the Fresnel matrix (in
Legendre functions of the reflection angle) and the wave slope distribution
of Cox and Munk (in Fourier series), the sunglint of the direct beam being
computed exactly. The maximum order :math:`S` is 2 for a molecular
atmosphere and ``OS_NB`` with aerosols.

Recombination in azimuth
~~~~~~~~~~~~~~~~~~~~~~~~

The terms :math:`I^s, Q^s, U^s` of all the orders are stored in the binary
file ``-SOS.ResBin``. ``SOS_TRPHI_MULTI`` (``SOS_TRPHI.F``) reads this file
**once**, keeps all the orders in memory and applies Eq. :eq:`fourier` to the
whole list of requested azimuths (two half-planes for ``-SOS.View 1``,
0 to 360° by ``-SOS.View.Dphi`` for ``-SOS.View 2``), with
:math:`\cos s\phi` and :math:`\sin s\phi` computed once per order and
azimuth. Up to version 8.2 the whole file was read again for each azimuth;
the results are identical.

The reflection of the direct solar beam by the surface (sun glint,
Fresnel reflection of a flat sea, Roujean BRDF, Rondeaux, Bréon, Nadal and
Maignan BPDF), which is removed from the successive orders, is then added
analytically for each azimuth (``SOS_TRPHI_DIRECT``). It is attenuated by
the direct transmission of the solar beam from TOA to the surface and from
the surface to the output level, computed by ``SOS_PROC``. With
``-SOS.AbsModeCKD 1`` this transmission is the CKD-weighted sum of the
transmissions of each k-term,
:math:`\sum_k a_k \exp[-\tau_k (1/\mu_s + 1/\mu_v)]` (and likewise for the
direct transmission at :math:`\mu_s` used for the fluxes), instead of the
Beer law applied to the aggregated optical thickness, which is exact only
for an air mass of 1.

Truncation of the aerosol phase function
----------------------------------------

The forward peak of the aerosol phase function is truncated (``-AER.Tronca``)
so that it can be expanded with a moderate number of terms. The truncated
energy is considered as not scattered: the aerosol optical thickness and
single scattering albedo are scaled with the truncation coefficient, and the
direct transmissions are computed with the true atmosphere.

Gaseous absorption
------------------

The gaseous absorption is modelled by the correlated-k distribution (CKD)
method [LacisOinas1991]_: over a spectral interval, the transmission of
each gas is a weighted sum of exponentials, with coefficients tabulated
versus pressure, temperature (and concentration for H\ :sub:`2`\ O) in
``fic/COEFF_CKD``. With ``-SOS.AbsModeCKD 1``, the radiative transfer is
solved for each combination of terms of the 8 gases and the results are
summed with their weights; with ``-SOS.AbsModeCKD 2``, a single computation is
made with the mean absorption optical thickness (see :doc:`running`).

Sequence of the modules
-----------------------

.. The flowchart is drawn with mermaid on the web pages and with graphviz in the
   PDF and ePub, where the mermaid JavaScript cannot run. Keep both in sync.

.. only:: format_html and not builder_epub

   .. mermaid::

      flowchart TD
          A["SOS_ANGLES<br/>Gauss and user angles"] --> B["SOS_AEROSOLS / SOS_MIE<br/>phase matrix, Greek coefficients"]
          B --> C["SOS_SURFACE<br/>BRDF / BPDF Fourier matrices"]
          C --> D["SOS_PREPA_ABSPROFILE<br/>gas profile, CKD coefficients"]
          D --> E["SOS_PREPA_OS"]
          E --> F{"CKD mode"}
          F -- "1: fine" --> G["for each CKD term:<br/>SOS_ABSPROFILE, SOS_PROFILE,<br/>SOS (SOS_OS), SOS_AGGREGATE"]
          F -- "2: simplified" --> H["mean absorption profile,<br/>SOS_PROFILE, SOS (SOS_OS)"]
          G --> I["SOS_TRPHI<br/>radiance at the requested azimuths"]
          H --> I

.. only:: not format_html or builder_epub

   .. graphviz::
      :align: center

      digraph algorithm {
          node [shape=box, fontname="DejaVu Sans", fontsize=10];
          A [label="SOS_ANGLES\nGauss and user angles"];
          B [label="SOS_AEROSOLS / SOS_MIE\nphase matrix, Greek coefficients"];
          C [label="SOS_SURFACE\nBRDF / BPDF Fourier matrices"];
          D [label="SOS_PREPA_ABSPROFILE\ngas profile, CKD coefficients"];
          E [label="SOS_PREPA_OS"];
          F [label="CKD mode", shape=diamond];
          G [label="for each CKD term:\nSOS_ABSPROFILE, SOS_PROFILE,\nSOS (SOS_OS), SOS_AGGREGATE"];
          H [label="mean absorption profile,\nSOS_PROFILE, SOS (SOS_OS)"];
          I [label="SOS_TRPHI\nradiance at the requested azimuths"];
          A -> B -> C -> D -> E -> F;
          F -> G [label="1: fine"];
          F -> H [label="2: simplified"];
          G -> I; H -> I;
      }

The Mie files and the surface matrices are stored in databases and reused
(see :doc:`running`). The source files are described in
:doc:`installation`.
