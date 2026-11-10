The fourth Stokes parameter V
=============================

SOS-ABS solves the vector transfer equation for the first three Stokes
parameters :math:`(I, Q, U)`. The circular polarisation :math:`V` is not
computed. This page explains why it is neglected, when it matters, and how the
code could be extended to the full Stokes vector by following the formalism
of Deuzé et al. (1989) [Deuze1989]_, which is already written for
:math:`(I, Q, U, V)`.

What is neglected
-----------------

For randomly oriented particles with a plane of symmetry, the phase matrix in
the scattering plane has six independent elements:

.. math::
   :label: phase_matrix_4x4

   \mathbf{F}(\Theta) =
   \begin{pmatrix}
   F_{11} & F_{12} & 0 & 0 \\
   F_{12} & F_{22} & 0 & 0 \\
   0 & 0 & F_{33} & F_{34} \\
   0 & 0 & -F_{34} & F_{44}
   \end{pmatrix}

with :math:`F_{22} = F_{11}` and :math:`F_{44} = F_{33}` for spheres. SOS-ABS
keeps the upper-left :math:`3 \times 3` block: :math:`F_{34}` and
:math:`F_{44}` are dropped. Accordingly:

- ``SOS_MIE.F`` computes :math:`F_{11}`, :math:`F_{12}` and :math:`F_{33}`
  (``IMIE``, ``QMIE``, ``UMIE``) but not :math:`F_{34}`, which is the
  imaginary part of :math:`S_2 S_1^*` (the real part gives :math:`F_{33}`);
- the aerosol file (``-AER.ResFile``) holds four Greek coefficients
  :math:`\alpha_l, \beta_l, \gamma_l, \zeta_l`; the external phase matrix
  (``-AER.ExtData``) has no :math:`F_{34}` or :math:`F_{44}` column;
- the surface matrices are :math:`3 \times 3` (``P11`` to ``P33``);
- the Fourier components, the binary file ``-SOS.ResBin`` and the outputs
  hold :math:`I, Q, U`.

Why it is usually small
-----------------------

The sunlight is not polarised. A first scattering produces linear
polarisation only (:math:`F \cdot (1, 0, 0, 0)^T = (F_{11}, F_{12}, 0, 0)^T`).
Circular polarisation appears when a beam already polarised, with
:math:`U \neq 0` in the new scattering plane, is scattered again by a
particle with :math:`F_{34} \neq 0`:

.. math::

   V' = -F_{34}\, U + F_{44}\, V .

:math:`V` is therefore a **multiple scattering** effect that needs particles
with a non-zero :math:`F_{34}`. For Rayleigh scattering, :math:`F_{34} = 0`
and a molecular atmosphere produces no circular polarisation; :math:`F_{34}`
is significant for large or absorbing particles. At the surface, the
:math:`(3,4)` element of the Fresnel matrix is
:math:`\mathrm{Im}(r_l r_r^*)`, zero when the refractive index is real:
Deuzé et al. note that the imaginary part of the sea water index is negligible
in the visible and near-infrared, so that :math:`R_{34} = R_{43} = 0`.

As a result, :math:`V` is generally much smaller than the linear polarisation
in clear atmospheres [HansenTravis1974]_, and its feedback on
:math:`I, Q, U` is negligible: this is why most vector codes, including
SOS-ABS and OSOAA, solve the :math:`3 \times 3` problem.

Why compute it
--------------

- **Completeness and validation**: comparisons with codes that solve for the
  full Stokes vector, such as adding-doubling codes [deHaan1987]_, and with
  vector benchmarks [Kokhanovsky2010]_.
- **Absorbing and non-spherical aerosols** (dust, smoke), for which
  :math:`F_{34}` and :math:`F_{44} \neq F_{33}` are larger, with the
  spheroid phase matrices already accepted by ``-AER.ExtData``.
- **Surfaces with a complex index**, and total internal reflection at the
  water-air interface (in coupled atmosphere-ocean codes), where the Fresnel
  coefficients are complex and the reflection converts linear into circular
  polarisation.
- **Circular polarisation as a biosignature**: homochiral molecules of
  photosynthetic organisms produce a circular polarisation signal
  [Sparks2009]_.

The way forward: extending Deuzé et al. (1989)
----------------------------------------------

The paper of Deuzé, Herman and Santer gives every ingredient for the four
components; the code kept three of them. The steps below follow its
sections.

1. Fourier series of the Stokes vector
   (Eq. 4 of [Deuze1989]_). :math:`V` has the parity of :math:`U`:

   .. math::

      \begin{pmatrix} I \\ Q \\ U \\ V \end{pmatrix}(\tau, \mu, \phi)
      = \sum_{s} (2 - \delta_{0s})
      \begin{pmatrix} I^s \cos s(\phi - \phi_s) \\
                      Q^s \cos s(\phi - \phi_s) \\
                      U^s \sin s(\phi - \phi_s) \\
                      V^s \sin s(\phi - \phi_s) \end{pmatrix}

   so the equations for each term :math:`s` stay independent (Eq. 5), with a
   :math:`4 \times 4` kernel :math:`\mathbf{P}^s(\mu, \mu')` (Eqs. 7-10). The
   rotation matrices :math:`\mathbf{L}(\chi)` (Eq. 3) leave :math:`V`
   unchanged.

2. Phase matrix expansion (Eq. 6 and Appendix, Eqs. A3-A6). Add the
   coefficients of :math:`F_{34}` (:math:`\epsilon_l`, expanded with the same
   generalised spherical functions :math:`P^l_{0,2}` as :math:`\gamma_l` for
   :math:`F_{12}`) and, for non-spherical particles, of :math:`F_{44}`
   (:math:`\delta_l`, Legendre polynomials). In the code:

   - ``SOS_MIE.F``: store :math:`F_{34}` with the other elements (one more
     array, from ``RES1, IMS1, RES2, IMS2``), which changes the format of the
     Mie files (rename them to avoid mixing old and new databases);
   - ``SOS_AEROSOLS.F``: average :math:`F_{34}` (and :math:`F_{44}`) over the
     size distributions and the mixtures, apply the truncation, compute
     :math:`\epsilon_l` (and :math:`\delta_l`), and add the columns to the
     ``-AER.ResFile``, ``-AER.UserFile`` and ``-AER.ExtData`` formats;
   - ``SOS.h``: no new dimension is needed, the arrays of order
     ``CTE_OS_NB_MAX`` are reused.

3. Surface reflection (Eqs. 13-34). Use the full Fresnel matrix of Eq. (14),
   with :math:`R_{33} = R_{44} = \mathrm{Re}(r_l r_r^*)` and
   :math:`R_{34} = -R_{43} = \mathrm{Im}(r_l r_r^*)`; expand it in Legendre
   functions of :math:`\Omega = \pi - 2\omega` like the phase matrix (Eq. 24),
   then combine with the Fourier series of the Cox-Munk distribution
   (Eqs. 26-32). The matrices :math:`\mathbf{M}^c_p, \mathbf{M}^s_p` keep the
   block structure of Eq. (34) with :math:`\mathbf{D} = \mathrm{diag}(1, 1,
   -1, -1)`. For the land BPDF models (Rondeaux, Bréon, Maignan), which are
   built on the Fresnel matrix with a real index, only the :math:`(4,4)`
   element is added. The sunglint correction of the direct beam (Eqs. 35-36)
   is unchanged except for the fourth component. Files: ``SOS_SURFACE.F``,
   ``SOS_GLITTER.F``, ``SOS_SURFACE_BPDF.F``, ``SOS_ROUJEAN.F``,
   ``SOS_NOM_FIC_SURFACE.F`` (new file names for the 4x4 matrices).

4. Successive orders (Eqs. 37-40). In ``SOS_OS.F``, carry
   :math:`V^s_{(n)}` with :math:`U^s_{(n)}` in the source functions, the
   integrations over :math:`\mu'` and the surface reflection; add the
   convergence tests on :math:`V` and write :math:`V^s` in the binary file
   (``SOS.F``, ``SOS_AGGREGATE.F`` for the CKD weighting).

5. Outputs. ``SOS_TRPHI.F``: recombine
   :math:`V(\phi) = \sum_s (2 - \delta_{0s}) V^s \sin s\phi`, add a ``V``
   column and the degree of circular polarisation :math:`V / I` to the
   radiance files, and distinguish the degree of linear polarisation
   :math:`\sqrt{Q^2 + U^2} / I` from the total degree of polarisation
   :math:`\sqrt{Q^2 + U^2 + V^2} / I`. ``SOS_PROC.F`` and the Python binding:
   return ``v_up`` and ``v_down``.

6. Validation. Check that :math:`V = 0` for a Rayleigh atmosphere and
   that :math:`(I, Q, U)` are unchanged when :math:`F_{34} = 0`; then compare
   with the vector benchmarks and an independent 4-component code for
   aerosol cases.

The extra cost is moderate: the matrix products grow from :math:`3 \times 3`
to :math:`4 \times 4` and :math:`V^s` shares the sine terms with
:math:`U^s`. A run-time option (as ``-SOS.Ipolar`` for the polarisation)
would keep the :math:`3 \times 3` mode and the current file formats.

.. todo::

   Implement the full Stokes vector (I, Q, U, V) following the steps above.

See :doc:`/references` for the full references.
