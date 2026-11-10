C******************************************************************************
C* Copyright 2023, Centre National d'Etudes Spatiales (CNES)
C* 
C* This file is part of the SOS-ABS radiative transfer code.
C* 
C* SOS-ABS is free software: you can redistribute it and/or modify
C* it under the terms of the GNU General Public License as published by
C* the Free Software Foundation, either version 3 of the License, or
C* (at your option) any later version.
C* 
C* SOS-ABS is distributed in the hope that it will be useful,
C* but WITHOUT ANY WARRANTY; without even the implied warranty of
C* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
C* GNU General Public License for more details.
C* 
C* You should have received a copy of the GNU General Public License
C* along with SOS-ABS. If not, see <http://www.gnu.org/licenses/>.
C******************************************************************************


C****************************************************************************
C* FILE: SOS_OS.F
C* PROJECT: Successive orders of scattering
C* ROLE: This module manages the procedures for solving the
C*       Radiative Transfer Equation, with polarisation, by successive
C*       orders of scattering.
C*
C* AUTHOR:
C*   Initial code: Laboratoire d'Optique Atmospherique (LOA, Villeneuve d'Ascq).
C*   Code rework and development: CS GROUP France, B.Lafrance
C* DATE: 30/04/02
C*
C* MOD:VERSION:1.0: Upgrade of the LOA OS code
C* MOD:VERSION:2.0: Correction of an error on the surface polarisation
C* MOD:VERSION:2.1: 02/07/2008 Adjustment of the display format of the profile in the
C*                             trace file
C*                             ==> Level number on 4 characters instead of 3.
C* MOD:VERSION:3.0: 10/01/2010 
C*   - Change of the size of the arrays of radiance fields:
C*     ==> Changed from CTE_OS_NBMU to CTE_OS_NBMU_MAX.
C*
C*   - Change of the size of the arrays of the expansion orders:
C*     ==> Changed from CTE_OS_NB to CTE_OS_NB_MAX.
C*
C*   - Added input parameter: NBMU (nb of actual angles for the simulations)
C*   - Added input parameter: OS_NB (max order of the expansions of the phase
C		                    functions in Legendre polynomials to be used).
C*
C*   - Change of the record line of the radiance arrays
C*     at order IS of the Fourier series decomposition
C*     ==> limited to the useful bounds -NBMU:NBMU
C*                  
C*   - Change of the reading line of the arrays of the reflection
C*     matrix at order IS of the Fourier series decomposition:
C*     ==> limited to the useful bounds -NBMU:NBMU
C* 
C*   - Adaptation of routines SOS_MAT_FRESNEL_PLAN_REFL, SOS_NOYAUX,
C*     SOS_INTEGR_EPOPT, SOS_FSOURCE_ORDRE1, SOS_FSOURCE_ORDREIG,
C*     SOS_FSOURCE_DIFF_FRESNEL1, SOS_PARAM_CONV, SOS_ARRET_DIFFUS, 
C*     SOS_ARRET_FOURIER and SOS_AJOUT_QUEUE to the passing of the parameter
C*     LUM_NBMU (nb of useful angles).
C*  
C*   - Adaptation of routine SOS_NOYAUX to the passing of OS_NB
C*     (to loop on the expansion orders).
C*
C*   - The computation of the incident diffuse flux (normalised by the solar irradiance) EMOINS
C*     is passed as output argument
C*     (to provide the diffuse transmission if the ground is black)
C*
C*   - Correction of a loop dimension error for the cancellation of the
C*     surface polarisation (IPOLAR=0 & IMAT_SURF=1): loop start 0 --> 1
C*
C* MOD:VERSION:3.1: 06/06/2013
C*
C*   -  Change of the value of Pi, declared as a constant in the file
C*      (INCTE_PI replaces CTE_PI)
C*
C*   -  Change of the constant SOS_NT replaced by the parameter NT
C*
C*   -  The identifier of the log, previously hard-coded, is replaced by the parameter
C*      IDLOG given as input
C*
C*   -  Addition of SOS_ARRET_DIFFUS_1 and SOS_ARRET_DIFFUS changed to SOS_ARRET_DIFFUS_2
C*
C* MOD:VERSION:3.2: 24/02/2016 
C*    - Minor adjustments to strictly respect the 72 columns
C*
C* MOD:VERSION:3.3: 31/08/2022
C*    - Constants renamed from SOS_* to CTE_*
C*    - Traces converted to English.
C*
C* MOD:VERSION:3.4: 26/09/2023
C*    - Removal of a comma after the command READ(...) when reading the surface matrices.
C*      No impact on the simulation results, but raised a warning at compilation.
C*    - Change of the syntax of multiple loops in routines SOS_NOYAUX, SOS_INTEGR_EPOPT,
C*      SOS_FSOURCE_ORDREIG and SOS_FSOURCE_DIFF_FRESNEL1
C*      for compatibility with Fortran 2018 for the gcc/10.2.0 compiler
C*    - Change of the syntax of obsolete IF tests in routines SOS_FSOURCE_DIFF_FRESNEL1,
C*      SOS_FSOURCE_ORDREIG, SOS_PARAM_CONV, SOS_ARRET_DIFFUS_1 and SOS_ARRET_DIFFUS_2
C*      for compatibility with Fortran 2018 for the gcc/10.2.0 compiler
C*
C*
C****************************************************************************

C----------------------------------------------------------------------------
C Definition of constants
C---------------------------------------------------------------------------- 
C Constants used:
C    CTE_LENDIR : Length of directory names.
C     CTE_LENFIC2 : Length of file names with directory tree.
C    INCTE_PI : Value of PI
C    CTE_OS_NB_MAX :  Maximum value of the limit order of the expansion in
C                     Legendre polynomials (array sizing).
C    CTE_OS_NT  : Number of layers of the atmospheric profile.
C    CTE_OS_NBMU_MAX : Maximum number of positive angles for solving the RTE
C    CTE_OS_IBOR : Minimum order of the Fourier series expansion
C    CTE_PH_SEUIL_CV_SG : Convergence threshold of the geometric series.
C    CTE_PH_SEUIL_SUMDIF : Stop threshold of the multiple scattering.
C    CTE_PH_SEUIL_SF : Threshold to stop the Fourier series
C                      decomposition.
C----------------------------------------------------------------------------
#include "SOS.h"
#define INCTE_PI DACOS(-1.D+00)

C*==========================================================================
C PROCEDURE: SOS_OS
C ==========
C      This procedure computes the diffuse radiation field (I, Q, U) in
C      a scattering atmosphere.
C      The result field is decomposed in Fourier series.
C
C Description of parameters
C ---------------------------
C     NBMU (I4) : (I) Number of actually useful (positive) angles
C
C	 RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (DOUBLE) : (I) Cosines of the angles.
C     GA(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (DOUBLE)  : (I) Weights of the angles.
C
C     OS_NB (I4) : (I) Max order of the expansions of the phase functions
C				       in Legendre polynomials to be used.
C
C             NT (I4) : (I) Number of layers in the atmospheric profile
C
C     FICSURF (CHARACTER*CTE_LENFIC2)  (I) :
C                           name of the file containing the reflection matrix
C                   		            files (Glitter file or other BRDF/BPDF).
C                                     (directory + file name + extension)
C
C     FICOS   (CHARACTER*CTE_LENFIC2)  (I) :
C                           name of the OS result file
C                                     (directory + file name + extension)
C
C     N0  (I4)   (I) : Number of the Gauss mu for the solar zenith angle.
C			      (complementary angle: mus = -rmu(N0))
C			      (N0 > 0 => to take it into account).
C
C     TETAS (double)   (I) : Solar zenith angle.
C      		               (complementary angle: mus= -cos(tetas))
C                           (N0 < 0 => to take it into account).
C
C     RO (double)   (I) : Surface albedo.
C
C     IMAT_SURF  (I4)   (I) : Option of surface BRDF / BPDF
C      		        0 :Lambertian surface of albedo RHO.
C			1 :Lambertian surface + glitter.
C			        2 :Lambertian surface + Fresnel reflection.
C		       3 :Lambertian surface + Roujean BRDF.
C			4 :Lambertian surface + Roujean BRDF + Rondeaux BPDF.
C			5 :Lambertian surface + Roujean BRDF + Breon BPDF.
C			6 :Lambertian surface + Roujean BRDF + Nadal BPDF.
C			7 :Lambertian surface + Roujean BRDF + Maignan BPDF.
C
C     IFRESNEL  (I4)  (I) : Option for the computation of the Fresnel reflection on a flat sea
C                            0 : no reflection on a flat sea.
C                            1 : computation of the Fresnel reflection on a flat sea.
C
C     IND_SURF (double) (I) : Refractive index of the interface surface / atmosphere.
C
C     H(0:CTE_OS_NT) (double)   (I) : Optical thickness at each level.
C
C     XDEL(0:CTE_OS_NT) (double)   (I) : Percentage of aerosols per level.
C     YDEL(0:CTE_OS_NT) (double)   (I) : Percentage of molecules per level.
C
C     ZPROF(0:CTE_OS_NT) (double)   (I) : Altitudes of the levels.
C
C     RON (double)   (I) : (Molecular) depolarisation factor.
C
C     ALPHA(0:CTE_OS_NB_MAX) (double)   (I) : Parameters alpha(k).
C     BETA(0:CTE_OS_NB_MAX)  (double)   (I) : Parameters beta(k).
C     GAMMA(0:CTE_OS_NB_MAX) (double)   (I) : Parameters gamma(k).
C     ZETA(0:CTE_OS_NB_MAX)  (double)   (I) : Parameters zeta(k).
C
C      	ZOUT (F9) : (I) Altitude of the output level
C			         -1 =>  Standard output:
C				 Radiance emerging at TOA
C				         + downward radiance at ground level.
C			          z =>  Upward and downward radiance
C				         for altitude z (0 <= z <= CTE_TOA_ALT).
C		
C     IGMAX  (I4)   (I) : Maximum order of scattering.
C
C     IBORM  (I4)   (I) : Maximum order of the Fourier series expansion
C                         being processed:
C                            2 for a purely molecular atmosphere,
C                            OS_NB in presence of aerosols.
C
C     IPOLAR  (I4)   (I) : Option for simulation without polarisation.
C      			   0 : computations without polarisation.
C			          1 : normal computations with polarisation.
C
C     TRACE (logical)   (I) : = true if writing in the trace file
C              (logical unit number of the trace file =99)
C     IDLOG (I4) : (I) Identifier of the log file
C
C     EMOINS (double)   (O) : Diffuse flux illuminating the ground, normalised by
C                   			 the incident solar flux at TOA.
C
C     EPLUS (double)   (O) : Upward diffuse flux at TOA, normalised by
C                   			 the incident solar flux at TOA.
C
C       IER  (I4)   (O) : error code =0 if no error, =-1 otherwise
C
C
C Other result provided by the procedure
C -------------------------------------
C   Binary OS result file
C      Information in the trace file (optional)
C         (logical unit number of the trace file =99,
C    file created/opened by the calling program)
C
C 
C Environment variable
C ------------------------
C   None
C
C
C Description of the files used
C ---------------------------------
C   -->Content of the Wind file
C       Unformatted binary file, created by SOS_SURFACE
C       It contains the surface reflection matrix
C       (Glitter file or other BRDF/BPDF):
C       i.e. the elements Rij of the reflection matrix
C              R11(NBMU,NBMU) (real),
C              R12(NBMU,NBMU) (real),
C              R13(NBMU,NBMU) (real),
C              R21(NBMU,NBMU) (real),
C              R22(NBMU,NBMU) (real),
C              R23(NBMU,NBMU) (real),
C              R31(NBMU,NBMU) (real),
C              R32(NBMU,NBMU) (real),
C              R33(NBMU,NBMU) (real)
C
C        --> Records on the domain limited to the useful angles
C                                   --------------------------------
C
C   -->Content of the OS result file
C       Unformatted binary file, created by SOS_OS
C       It contains the Stokes parameters Q, U and I for an order S
C       of the Fourier series expansion, versus the angle mu
C       (double precision values).
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C     -Opening of a file
C       Error if the surface file does not exist
C     -Reading/writing in a file
C     -Inconsistency of the output level in the atmospheric profile
C     -Error in a called routine
C
C     Message displayed on screen, the routine stops its computations and
C     returns status -1 to the calling program
C
C     Remark: If the OS result file already exists, it is overwritten
C                by the new file created by this routine
C                (no message and no error)
C
C      Warning case
C     (message without program stop and without error return)
C            Stop of the OS computations for an incidence at the limb.
C
C
C  Routines used by the procedure:
C ------------------------------------
C  This program calls the routines:
C      - SOS_MAT_FRESNEL_PLAN_REFL
C      - SOS_NOYAUX
C      - SOS_INTEGR_EPOPT
C      - SOS_FSOURCE_ORDRE1  
C      - SOS_FSOURCE_ORDREIG
C      - SOS_FSOURCE_DIFF_FRESNEL1 
C      - SOS_PARAM_CONV 
C      - SOS_ARRET_DIFFUS_1
C      - SOS_ARRET_DIFFUS_2
C      - SOS_ARRET_FOURIER
C      - SOS_AJOUT_QUEUE
C
C===========================================================================
      SUBROUTINE SOS_OS(NBMU,RMU,GA,OS_NB,NT,
     &                  FICSURF,FICOS,
     &                  N0,TETAS,RO,IMAT_SURF,IFRESNEL,IND_SURF,
     &                  H,XDEL,YDEL,ZPROF,RON,
     &                  ALPHA,BETA,GAMMA,ZETA,ZOUT,IGMAX,IBORM,
     &                  IPOLAR,TRACE,IDLOG,EMOINS,EPLUS,IER)

      IMPLICIT NONE


C* Definition of variables
C*------------------------------------------------------
				
      DOUBLE PRECISION TETAS ! Solar zenith angle in degrees (0 < tetas < 90).
      				 ! (complementary angle: mus= -cos(tetas))
      				 ! (N0 < 0 => to take it into account).
	
      DOUBLE PRECISION RO	! Surface albedo.
	
      DOUBLE PRECISION IND_SURF  ! Refractive index surface / atmosphere.
      				
      DOUBLE PRECISION H(0:CTE_OS_NT)	!Optical thickness for each level.
      DOUBLE PRECISION XDEL(0:CTE_OS_NT)  !Percentage of aerosols per level.
      DOUBLE PRECISION YDEL(0:CTE_OS_NT)  !Percentage of molecules per level.
      DOUBLE PRECISION ZPROF(0:CTE_OS_NT) !Altitudes of the profile levels

      DOUBLE PRECISION RON	! (Molecular) depolarisation factor.

      DOUBLE PRECISION ALPHA(0:CTE_OS_NB_MAX)	! Parameters alpha(k).
      DOUBLE PRECISION BETA(0:CTE_OS_NB_MAX)	! Parameters beta(k).
      DOUBLE PRECISION GAMMA(0:CTE_OS_NB_MAX)	! Parameters gamma(k).
      DOUBLE PRECISION ZETA(0:CTE_OS_NB_MAX)	! Parameters zeta(k).

        				
      DOUBLE PRECISION I3ZOUT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Stokes parameter I for an order S
      				! of the Fourier series expansion,
				! versus the angle mu.
				! At the output altitudes
      DOUBLE PRECISION Q3ZOUT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Stokes parameter Q.
      DOUBLE PRECISION U3ZOUT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Stokes parameter U.

      DOUBLE PRECISION I3OUT(0:CTE_OS_NT,
     & 			     -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Stokes parameter I for an order S
      				! of the Fourier series expansion,
				! versus the angle mu.
				! For each level of the atmospheric profile
      DOUBLE PRECISION Q3OUT(0:CTE_OS_NT,
     & 			     -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Stokes parameter Q.
      DOUBLE PRECISION U3OUT(0:CTE_OS_NT,
     & 			     -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Stokes parameter U.

      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
                                               ! Array of the Gauss mu.
      DOUBLE PRECISION GA(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   
                                               ! Array of the Gauss weights.

      DOUBLE PRECISION RRMU	! Double inversion of mu.

      DOUBLE PRECISION TAB      ! Value of the cosine of the solar zenith
                                ! angle (mus).

      DOUBLE PRECISION  ZZ      ! Linear interpolation coefficient between 2 altitudes

C*  - Lambertian reflection
      DOUBLE PRECISION XR(0:CTE_OS_NBMU_MAX) ! radiance field of the Lambertian
                                             ! reflection at the ground.

      DOUBLE PRECISION LSOL   ! Radiance of the Lambertian reflection at the ground.

C*  - Contribution of the direct surface reflection (BRDF/BPDF)
C*    to the signal at TOA and at the chosen output level
      DOUBLE PRECISION RII(CTE_OS_NBMU_MAX)	     ! Contribution to I TOA.
      DOUBLE PRECISION RQQ(CTE_OS_NBMU_MAX)	     ! Contribution to Q TOA.
      DOUBLE PRECISION RUU(CTE_OS_NBMU_MAX)	     ! Contribution to U TOA.
      DOUBLE PRECISION RIIOUT(0:CTE_OS_NT,CTE_OS_NBMU_MAX)   ! Contribution to I
                                             		     ! at each level of the profile.
      DOUBLE PRECISION RQQOUT(0:CTE_OS_NT,CTE_OS_NBMU_MAX)   ! Contribution to Q
                                             		     ! at each level of the profile
      DOUBLE PRECISION RUUOUT(0:CTE_OS_NT,CTE_OS_NBMU_MAX)   ! Contribution to U
                                             		     ! at each level of the profile

      DOUBLE PRECISION A	! Transmission in the viewing direction
      				! between the ground and the output level.

      DOUBLE PRECISION RR	! Attenuation of the direct solar beam at the ground
      				! divided by mus.

      DOUBLE PRECISION F11(0:CTE_OS_NBMU_MAX)  ! Element of the Fresnel matrix.
      DOUBLE PRECISION F12(0:CTE_OS_NBMU_MAX)  ! Element of the Fresnel matrix.
      DOUBLE PRECISION F33(0:CTE_OS_NBMU_MAX)  ! Element of the Fresnel matrix.

C*  - For the molecular scattering:
      DOUBLE PRECISION BETA0	! Coefficients of the expansion of the molecular phase
      DOUBLE PRECISION BETA2	! matrix in Legendre functions
      DOUBLE PRECISION GAMMA2	
      DOUBLE PRECISION ALPHA2	

      DOUBLE PRECISION AAA	! Correction factor of the molecular depolarisation.

      DOUBLE PRECISION XPL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      			! Legendre function PSL at order L=2 and S=IS
      			! versus the angle mu.
      DOUBLE PRECISION XRL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      			! Function RSL at order L=2 and S=IS versus mu.
      DOUBLE PRECISION XTL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      			! Function TSL at order L=2 and S=IS versus mu.

C*  - For the aerosol scattering:
	! Element Pij of the aerosol phase matrix for an order
	! IS of the Fourier series decomposition.
	! With mup the incidence direction (-1 < mup < 1)
	! and mu the scattering one (-1 < mu < 1),
	! the scattering is expressed by Pij(mu,mup).
      DOUBLE PRECISION 
     &    BP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &       -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  ! Element P11.
      DOUBLE PRECISION 
     &    GR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &       -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  ! P12 or P21.
      DOUBLE PRECISION 
     &    GT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &       -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  ! -P13 or -P31.
      DOUBLE PRECISION 
     &    ARR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &        -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! P22.
      DOUBLE PRECISION 
     &    ART(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &        -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! -P23 or -P32.
      DOUBLE PRECISION 
     &    ATT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &        -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! P33.

      DOUBLE PRECISION CH(0:CTE_OS_NT)  ! Attenuation of the direct solar beam
                                        ! / 4 for the aerosols + molecules mixture

C* Variables associated with the spatial and vertical integration:
      	! Radiance field of order IG-1 or IG versus the level of the profile
	! and the direction.
      DOUBLE PRECISION I1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	! Field I.
      DOUBLE PRECISION Q1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	! Field Q.
      DOUBLE PRECISION U1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	! Field U.

      DOUBLE PRECISION 
     &    I1fresnel(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      			! Field I for the primary scattering of the
			! solar radiation going up
			! after a direct reflection on a flat sea,
			! according to the Fresnel laws.
      DOUBLE PRECISION 
     &    Q1fresnel(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      			! Field Q.
      DOUBLE PRECISION 
     &    U1fresnel(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      			! Field U.
      
	! Source function vector versus the level of the profile
	! and the direction.
      DOUBLE PRECISION I2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      						! Source function for I.
      DOUBLE PRECISION Q2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      						! Source function for Q.
      DOUBLE PRECISION U2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      						! Source function for U.

C* Variable associated with the loop exit tests:
      DOUBLE PRECISION Z1	! Threshold variable for the stop tests.
     	
C*  - Loop on the order of the Fourier series decomposition:
	 ! Arrays summing the orders S of the decomposition
        ! in Fourier series, versus the angle mu.
      DOUBLE PRECISION I4(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      DOUBLE PRECISION Q4(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      DOUBLE PRECISION U4(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)

        ! Arrays summing the orders S of the decomposition
        ! in Fourier series, with permutation of the sign of the orders S,
	 ! versus the angle mu.
      DOUBLE PRECISION I5(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      DOUBLE PRECISION Q5(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      DOUBLE PRECISION U5(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	

      DOUBLE PRECISION SIGN	! Value 1 or -1.
      DOUBLE PRECISION COEF	! Coefficient 1 or 2.

C*  - Loop on the scatterings:
	! Functions for the convergence test in geometric series,
	! field emerging at TOA (mu > 0) and downward at the ground (mu < 0).
	! The functions are established for an order S of the expansion
	! in Fourier series, versus the angle mu.
      DOUBLE PRECISION I3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)		
                            ! Stokes parameter I for the sum of the
      				! previous scattering orders up to the current IG.
      DOUBLE PRECISION Q3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
      				! Stokes parameter Q.
      DOUBLE PRECISION U3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
      				! Stokes parameter U.
            			
      DOUBLE PRECISION A1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
       			! Stokes parameter I for the order IG-2
				! of the field emerging at TOA (mu>0),
				! or downward at ground level (mu<0).
      DOUBLE PRECISION B1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Stokes parameter Q for the order IG-2.
      DOUBLE PRECISION C1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Stokes parameter U for the order IG-2.

      DOUBLE PRECISION D1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Stokes parameter I for the order IG-1.
      DOUBLE PRECISION E1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Stokes parameter Q for the order IG-1.
      DOUBLE PRECISION F1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      			  	! Stokes parameter U for the order IG-1.

      DOUBLE PRECISION G1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Stokes parameter I for the order IG.
      DOUBLE PRECISION H1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Stokes parameter Q for the order IG.
      DOUBLE PRECISION P1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      			  	! Stokes parameter U for the order IG.

	! Same functions for the output level of the results:
      DOUBLE PRECISION D1OUT(0:CTE_OS_NT,
     & 			 -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Stokes parameter I for the order IG-1.
      DOUBLE PRECISION E1OUT(0:CTE_OS_NT,
     &			        -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Stokes parameter Q for the order IG-1.
      DOUBLE PRECISION F1OUT(0:CTE_OS_NT,
     &			        -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      			  	! Stokes parameter U for the order IG-1.

      DOUBLE PRECISION G1OUT(0:CTE_OS_NT,
     &			        -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Stokes parameter I for the order IG.
      DOUBLE PRECISION H1OUT(0:CTE_OS_NT,
     &			        -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Stokes parameter Q for the order IG.
      DOUBLE PRECISION P1OUT(0:CTE_OS_NT,
     &  		        -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      			  	! Stokes parameter U for the order IG.

      DOUBLE PRECISION EMOINS	! Diffuse flux illuminating the ground, normalised by
      					! the incident solar flux at TOA.

      DOUBLE PRECISION EPLUS	! Diffuse flux emerging at TOA, normalised by
      					! the incident solar flux at TOA.

      DOUBLE PRECISION II2 	! Intermediate value for the computation of I2
                                ! (mu > 0).
      DOUBLE PRECISION QQ2 	! Intermediate value for the computation of Q2
                                ! (mu > 0).
      DOUBLE PRECISION UU2 	! Intermediate value for the computation of U2
                                ! (mu > 0).
      DOUBLE PRECISION XI1	! Value of I of the upward field of order IG-1.
      DOUBLE PRECISION XQ1	! Value of Q of the upward field of order IG-1.
      DOUBLE PRECISION XU1	! Value of U of the upward field of order IG-1.

      DOUBLE PRECISION Z	! Single value of the array of Gauss weights.


C*  - Elements Rij of the surface reflection matrix
      REAL*4 R11(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL*4 R12(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL*4 R13(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL*4 R21(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL*4 R22(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL*4 R23(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL*4 R31(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL*4 R32(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL*4 R33(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)

								
      CHARACTER*CTE_LENFIC2 FICSURF  ! File containing the surface reflection
      				         ! matrix (Glitter file or other BRDF/BPDF).

      CHARACTER*CTE_LENFIC2 FICOS    ! Name of the OS result file.
      

      LOGICAL TRACE     ! = true if writing in the trace file

      INTEGER*4 NT	   ! Number of layers of the atmospheric profile
      
      INTEGER*4 OS_NB 	! Max order of the expansions of the phase functions
				! in Legendre polynomials.

      INTEGER*4 N0		! Number of the Gauss mu for the angle
      				! zenith angle.
				! (complementary angle: mus = -rmu(N0))
				! (N0 > 0 => to take it into account).
               
      INTEGER*4 IMAT_SURF ! Option of surface BRDF / BPDF:
      			     !   0 : Lambertian surface of albedo RHO.
			     !   1 : Lambertian surface + Glitter.
			     !   2 : Lambertian surface + Fresnel reflection.
		            !   3 : Lambertian surface + Roujean BRDF.
			     !   4 : Lambertian surface + Roujean BRDF + Rondeaux BPDF.
			     !   5 : Lambertian surface + Roujean BRDF + Breon BPDF.
			     !   6 : Lambertian surface + Roujean BRDF + Nadal BPDF.
			     !   7 : Lambertian surface + Roujean BRDF + Maignan BPDF.
			  
      INTEGER*4 IFRESNEL ! Index specifying whether the simulations
      			    ! include the Fresnel reflection on a flat sea:
			    !    1 => for a computation of the Fresnel reflection.
			    !    Otherwise => not taken into account.
			 
      DOUBLE PRECISION ZOUT	! Output type:
				!      -1 =>  Standard output:
				!             Radiance emerging at TOA
				!           + downward radiance at ground level.
				!       z =>  Upward and downward radiance
				!             for altitude z (0 <= z <= CTE_TOA_ALT).
			
 	
      INTEGER*4 IGMAX	! Maximum order of scattering.
      
      INTEGER*4 IBORM	! Maximum order of the expansion in Fourier series
      			       ! being processed:
			       !   _ 2 for a purely molecular atmosphere,
			       !   _ OS_NB in presence of aerosols.

      INTEGER*4 IPOLAR	! Option for simulation without polarisation.
      			       ! 0 : computations without polarisation.
			       ! 1 : normal computations with polarisation.
			
      INTEGER*4 JK	! Storage number of mus in the array RMU
      			! (0 by default).

      INTEGER*4 NUP_OUT    ! Output level of I, Q and U for the
                           ! upward field.
      INTEGER*4 NDOWN_OUT  ! Output level of I, Q and U for the
                           ! downward field.

      INTEGER*4 IS	      ! Current order of the Fourier series
                           ! expansion associated with the phase matrix
			   
      INTEGER*4 IG	      ! Current order of interaction
                           ! (scattering + reflection).
								
      INTEGER*4 IER        !error code =0 if no error, =-1 otherwise
      INTEGER*4 IFIN       !nb of significant characters of a string
      INTEGER*4 I,K,J      !array indices
      
      INTEGER*4 NBMU       ! number of useful angles

      INTEGER*4 IDLOG      ! Identifier number of the trace file

C* Initialisation of the error code
C------------------------------------------------------
      IER=0
      
       
C* Opening of the files
C------------------------------------------------------
C     Surface BRDF / BPDF file

      IF (IMAT_SURF.EQ.1) THEN
          OPEN(UNIT=12,FILE=FICSURF,FORM='UNFORMATTED',
     &         STATUS='OLD',ERR=910)      
      ENDIF  
     
C     OS result file
C     If this file already exists, it is overwritten by the new file
      IF (FICOS.NE.'NO_OUTPUT') THEN
          OPEN(UNIT=10,FILE=FICOS,FORM='UNFORMATTED',ERR=913)
      ENDIF

C* Account for the depolarisation factor in the parameters
C* of the molecular phase function.
C---------------------------------------------------------------------
      AAA=RON/(2-RON)
      AAA=(1-AAA)/(1+2*AAA)

      BETA0=1.
      BETA2=0.5*AAA
      GAMMA2=-AAA*DSQRT(1.5D+00)
      ALPHA2=3.*AAA


C* Cut of the atmospheric polarisation.
C---------------------------------------------------------------------
      IF (IPOLAR.EQ.0) THEN
          GAMMA2 = 0.D+00
          ALPHA2 = 0.D+00
	  
	  DO K=0,OS_NB
	     ALPHA(K) = 0.D+00
	     GAMMA(K) = 0.D+00
             ZETA(K)  = 0.D+00
   	  ENDDO
      
      ENDIF !Cut of the polarisation
      

C* Definition of mus (cosine of the solar zenith angle)
C* depending on whether the simulation is made for a Gauss angle
C* or for any angle.
C----------------------------------------------------------
      IF (N0.GT.0) THEN
         TAB=-RMU(N0)
      ELSE
         TAB=-DCOS(INCTE_PI*TETAS/180.)
      ENDIF

C* Storage of mus at index 0 of the array of the Gauss mu (RMU).
C------------------------------------------------------------------
      JK=0
      RMU(JK)=TAB

C* Trace file
C------------------------------------------------------  
      IF (TRACE) THEN
        WRITE(IDLOG,*,err=921) ' '  
        WRITE(IDLOG,*,err=921) ' '	   
        WRITE(IDLOG,*,err=921) 'SOS_OS SUBROUTINE'
	 WRITE(IDLOG,*,err=921) '================='
        WRITE(IDLOG,*,err=921) ' '
        WRITE(IDLOG,*,err=921) 'SIMULATION PARAMETERS'
	 WRITE(IDLOG,*,err=921) '---------------------'
	 WRITE(IDLOG,*,err=921)' '   
        WRITE(IDLOG,*,err=921)' Solar incidence'
	WRITE(IDLOG,*,err=921)'    N0 : ',N0
	IF (N0.LE.0) THEN
	   WRITE(IDLOG,*,err=921)'    Solar zenith angle : ',TETAS 
	ELSE
	   WRITE(IDLOG,*,err=921)
     &       '    Solar zenith angle : ', 
     &            DACOS(RMU(N0))*180.D+00/INCTE_PI
	ENDIF
	WRITE(IDLOG,*,err=921)'    Gauss angle mus : ' ,TAB
	WRITE(IDLOG,*,err=921)
	WRITE(IDLOG,*,err=921)
     &	'    Direct transmittance TOA -> surface : '
	WRITE(IDLOG,*,err=921)
     &'    (AOT adjusted to the truncated phase function, if any) ',
     &                          DEXP(H(NT)/(TAB))
        WRITE(IDLOG,*,err=921) 
        IF (IPOLAR.EQ.1) THEN
	   WRITE(IDLOG,*,err=921)'    Simulation including polarization'
	ELSE
	   WRITE(IDLOG,*,err=921)'    Simulation WITHOUT polarization' 
	ENDIF
	WRITE(IDLOG,*,err=921)' ' 
	WRITE(IDLOG,*,err=921)' Surface' 
	WRITE(IDLOG,*,err=921)'    Surface albedo: ',RO 
        IF (IMAT_SURF.EQ.1) THEN
	  WRITE(IDLOG,*,err=921)
     &   '    Simulation of a surface with BRDF / BPDF'
	  WRITE(IDLOG,*,err=921)'    Reflection matrix file : ', FICSURF
	ELSE
	  WRITE(IDLOG,*,err=921)'    No surface BRDF / BPDF ' 
	ENDIF
	IF (IFRESNEL.EQ.1) THEN
	  WRITE(IDLOG,*,err=921)
     &	   '   Simulation of the reflections on a flat sea'
	  WRITE(IDLOG,*,err=921)
     &	   '    Surface / atmosphere refractive index', IND_SURF  
	ELSE
	  WRITE(IDLOG,*,err=921)'    No reflection on a flat sea' 
	ENDIF				
	WRITE(IDLOG,*,err=921)' '
	WRITE(IDLOG,*,err=921)' Atmospheric profile ' 
	WRITE(IDLOG,*,err=921)
     &	   '  (adjusted to the truncated phase function, if applied)' 
	WRITE(IDLOG,*,err=921)'    Number of layers : ',NT 
	WRITE(IDLOG,*,err=921)
     &	     '      Layer   % aerosols    % molecules   Opt. Thickness'
       DO 3 J=0,NT
	   WRITE(IDLOG,'(8X,I4,5X,F10.5,8X,F10.5,8X,F10.5)',err=921)
     &	         J,XDEL(J),YDEL(J),H(J)
   3   CONTINUE
	WRITE(IDLOG,*,err=921)' '	   
	WRITE(IDLOG,*,err=921)
     &  '   Coefficients of the aerosol phase matrix:'
       WRITE(IDLOG,*,err=921) ' '  
	WRITE(IDLOG,35,err=921)'I','ALPHA(I)','BETA(I)','GAMMA(I)',
     &			    'ZETA(I)'  
       DO K=0,OS_NB
	   WRITE(IDLOG,25,err=921) K,ALPHA(K),BETA(K),GAMMA(K),ZETA(K)
   	ENDDO
	WRITE(IDLOG,*,err=921)' '	
	WRITE(IDLOG,*,err=921)
     &	     '    Coefficients of the molecule phase matrix'
	WRITE(IDLOG,*,err=921)
     &	     '         Depolarization factor : ',RON
	WRITE(IDLOG,*,err=921)
     &       '         Coef F for the depolarization correction : ',AAA
	WRITE(IDLOG,*,err=921)'         Beta0 = ',BETA0
	WRITE(IDLOG,*,err=921)'         Beta2 = ',BETA2	
	WRITE(IDLOG,*,err=921)'         Gamma2 = ',GAMMA2
	WRITE(IDLOG,*,err=921)'         Alpha2 = ',ALPHA2	
	WRITE(IDLOG,*,err=921)' '	
	WRITE(IDLOG,*,err=921)'RESULTS'
	WRITE(IDLOG,*,err=921)'-------'
      END IF


C* Stop of the computations for an incidence at the limb.
C------------------------------------------------------
      IF(TAB.EQ.0.D+00) GOTO 995

C* Consistency test of the output level in the atmospheric profile
C----------------------------------------------------------------
      IF (((ZOUT.LT.0).AND.(ZOUT.NE.-1.0)).OR.(ZOUT.GT.CTE_TOA_ALT))
     &	  GOTO 996
     

C* Computation of the reflection coefficients of the Fresnel matrix
C---------------------------------------------------------------
      IF (IFRESNEL.EQ.1) THEN
          CALL SOS_MAT_FRESNEL_PLAN_REFL(NBMU,RMU,IND_SURF,IPOLAR,
     &	                                 F11,F12,F33)
	  
	  IF (TRACE) THEN
	     WRITE(IDLOG,*,err=921) 
     &       ' Matrix of the Fresnel reflection: '
             DO J= 0,NBMU
                WRITE(IDLOG,40,err=921) 
     &		      'Reflection Angle ',DACOS(RMU(J))*180./INCTE_PI,
     &	              ' F11 ',F11(J),' F12 ',F12(J),' F33 ',F33(J)   
             ENDDO
             WRITE(IDLOG,*,err=921) ' '	   
	  ENDIF
	  
      ENDIF !Fresnel reflection on a flat sea


C* Attenuation of the direct solar beam (divided by 4).
C------------------------------------------------------
      DO 84 I=0,NT
         CH(I)=DEXP(-H(I)/(-TAB))/4.
  84  CONTINUE


C* Initialisation of the arrays for the stop test of the Fourier series
C* decomposition.
C-----------------------------------------------------------------------
      DO 22 J=-NBMU,NBMU
         Q4(J)=0.
         U4(J)=0.
         I4(J)=0.

         Q5(J)=0.
         U5(J)=0.
         I5(J)=0.
 22   CONTINUE


C--------------------------------------------------------------
C
C       ******  FOURIER SERIES DECOMPOSITION   ******
C
C--------------------------------------------------------------
      IF (TRACE) THEN
         WRITE (IDLOG,*,ERR=921) " " 
	  WRITE (IDLOG,*,ERR=921) "FOURIER SERIES EXPANSION"
	  WRITE (IDLOG,*,ERR=921) "************************"
	  WRITE (IDLOG,*,ERR=921) " " 
      ENDIF !Trace file

      SIGN=-1.
	
C* Loop on the orders of the Fourier series decomposition
C------------------------------------------------------
      DO 24 IS=CTE_OS_IBOR,IBORM
 	 
         SIGN=-SIGN
 
C* Trace file
C------------------------------------------------------  
         IF (TRACE) THEN
            WRITE(IDLOG,*,err=921) 
     &       ' Fourier series expansion order: IS = ',IS
            WRITE(IDLOG,*,err=921) 
     &       ' ------------------------------------------'
            WRITE(IDLOG,*,err=921) ' '	   
         ENDIF

C* Removal of the molecular coefficient BETA0 for the computation
C* of the molecular scattering phase matrix for the
C* orders IS > 0.
C--------------------------------------------------------------
         IF(IS.GT.0) BETA0=0.

C* Initialisation of the arrays giving the parameters I, Q, U
C* for the order IS of the Fourier series decomposition.
C------------------------------------------------------
         DO 23 J=-NBMU,NBMU
	 
            I3(J)=0.
            Q3(J)=0.
            U3(J)=0.
	    
	    DO I=0,NT
            	I3OUT(I,J)=0.
            	Q3OUT(I,J)=0.
            	U3OUT(I,J)=0.
	    ENDDO
	    
   23    CONTINUE

C* Reading of the elements of the reflection matrix for
C* the current order IS.
C--------------------------------------------------------------
         IF (IMAT_SURF.EQ.1) THEN
C*      Reading of the files that were written with
C*    array bounds limited to the useful
C*    dimension -NBMU:NBMU
	     READ(12,err=920)
     &	       ((R11(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((R12(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((R13(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((R21(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((R22(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((R23(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((R31(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((R32(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((R33(I,J),I=1,NBMU),J=1,NBMU)  
     	     
C*    Cut of the surface polarisation.
             IF (IPOLAR.EQ.0) THEN 	  
	         DO K=1,NBMU
		    DO J=1,NBMU
	               R12(J,K) = 0.
		       R13(J,K) = 0.
		       R21(J,K) = 0.
		       R22(J,K) = 0.
		       R23(J,K) = 0.
		       R31(J,K) = 0.
		       R32(J,K) = 0.
		       R33(J,K) = 0.	     
   	            ENDDO
		 ENDDO        
             ENDIF !IPOLAR=0
      
	 ENDIF  !IMAT_SURF=1
	   

C* Call of the routine computing the scattering kernels for
C* the order IS of the Fourier series decomposition.
C--------------------------------------------------------------
         CALL SOS_NOYAUX(IS,NBMU,RMU,OS_NB,ALPHA,BETA,GAMMA,ZETA,
     &                   XPL,XRL,XTL,BP,GR,GT,ARR,ART,ATT)



C--------------------------------------------------------------
C
C   ******  PRIMARY SCATTERING AND REFLECTION COMPUTATIONS   ******
C
C--------------------------------------------------------------

C* Computation of the primary scattering source function: I2, Q2, U2
C----------------------------------------------------
         CALL SOS_FSOURCE_ORDRE1(IS,NBMU,NT,JK,XDEL,YDEL,
     &                           BETA0,BETA2,GAMMA2,
     &                           XPL,XRL,XTL,BP,GR,GT,CH,
     &                           I2,Q2,U2)       


C* Boundary condition at the ground
C-------------------------
         DO 108 K=1,NBMU
            I1(NT,K)=0.
            Q1(NT,K)=0.
            U1(NT,K)=0.
            XR(K)=0.

C* Lambertian reflection at the ground. Stored in XR.
C------------------------------------------------
            IF (RO.EQ.0..OR.IS.NE.0) GOTO 1022
            I1(NT,K)=-RO*TAB*DEXP(H(NT)/TAB)
            XR(K)=I1(NT,K)

C* Surface reflection with BRDF / BPDF
C---------------------------------------
 1022       IF (IMAT_SURF.EQ.1) THEN
               RR = DEXP(H(NT)/TAB) / RMU(K)

               I1(NT,K)=  I1(NT,K) + R11(N0,K)*RR
               Q1(NT,K)=  R21(N0,K)*RR
               U1(NT,K)=  R31(N0,K)*RR    	
            ENDIF
	    
  108    CONTINUE


C* Integration of the upward radiance field (mu > 0) from the ground to TOA
C* and of the downward radiance field (mu < 0) from TOA to the ground.
C--------------------------------------------------------------
         CALL SOS_INTEGR_EPOPT(NBMU,RMU,NT,H,I2,Q2,U2,I1,Q1,U1)



C* Fresnel reflection:
C* --------------------
C* For a direct solar illumination, the upward radiance (L1fresnel)
C* obtained by Fresnel reflection is reduced to the direction mu = -mus.
C* It cannot be expanded in Fourier series.
C* On the other hand, the radiance of first scattering associated
C* with L1fresnel is expanded in Fourier series.
C--------------------------------------------------------------
         IF (IFRESNEL.EQ.1) THEN

C*   Computation of the source function for the scattering of L1fresnel
	     CALL SOS_FSOURCE_DIFF_FRESNEL1(
     &   	         F11(0),F12(0),XDEL,YDEL,NT,
     &                 BETA0,BETA2,GAMMA2,ALPHA2,BP,GR,GT,ARR,ART,
     &                 XPL,XRL,XTL,IS,NBMU,TAB,H,I2,Q2,U2)
     

C*   Integration over the atmospheric profile.
C*   Computation of the radiance field for the primary scattering of L1fresnel.
	 
             DO 111 K=1,NBMU
                 I1fresnel(NT,K)=0.	!Initialisation at the ground.
                 Q1fresnel(NT,K)=0.
                 U1fresnel(NT,K)=0. 
111	     CONTINUE
     
 	     CALL SOS_INTEGR_EPOPT(NBMU,RMU,NT,H,I2,Q2,U2,I1fresnel,
     &	                           Q1fresnel,U1fresnel)
     

C* Sum of the primary scattering terms:
C*   - of the direct solar beam + Lambertian reflection
C*   - of the solar beam reflected on a flat sea, according to the Fresnel laws.
            DO 112 I=0,NT
               DO 113 K=-NBMU,NBMU
                     I1(I,K)=I1(I,K)+I1fresnel(I,K)
                     Q1(I,K)=Q1(I,K)+Q1fresnel(I,K)
                     U1(I,K)=U1(I,K)+U1fresnel(I,K) 
113	       CONTINUE		     
112	    CONTINUE

	 ENDIF	!End of the computations for IFRESNEL=1
	 

	 
C* Storage in RII, RQQ and RUU of the contribution of the direct
C* surface reflection (BRDF / BPDF) to the signal for the requested
C* output level (TOA: Nup_out=0, ground: Nup_out = NT).
C--------------------------------------------------------------
         DO 4325 K=1,NBMU
            RII(K)=0. 	
            RQQ(K)=0.
            RUU(K)=0.
	     DO I=0,NT
            	RIIOUT(I,K)=0. 	
            	RQQOUT(I,K)=0.
            	RUUOUT(I,K)=0.
	     ENDDO
 4325    CONTINUE
		
         IF (IMAT_SURF.EQ.1) THEN

            DO 4321 K=1 , NBMU
	    
	       DO I=0,NT
	       
                  A = -( H(NT) - H(I) ) / RMU(K)
                  A= DEXP(A)
                  RIIOUT(I,K) = A*(I1(NT,K)-XR(K))
                  RQQOUT(I,K) = A*Q1(NT,K)
                  RUUOUT(I,K )= A*U1(NT,K)
	
	       ENDDO
	       
	       A = -H(NT)  / RMU(K)
              A= DEXP(A)
	       RII(K) = A*(I1(NT,K)-XR(K))
              RQQ(K) = A*Q1(NT,K)
              RUU(K )= A*U1(NT,K)
	       
 4321       CONTINUE

         ENDIF


C------------------------------------------------------------
C* Storage of the results of the computation of the interaction of order 1
C* (primary scattering and direct reflection).
C-------------------------------------------------------------

C* Storage of the downward field of order IG = 1.
C------------------------------------------------
         DO K=-NBMU,-1

C* At the ground for the convergence test:
            I3(K)=I1(NT,K)
	     Q3(K)=Q1(NT,K)
            U3(K)=U1(NT,K)
            D1(K)=I1(NT,K)
            E1(K)=Q1(NT,K)
            F1(K)=U1(NT,K)
	
C* At the output level of the results:
	   DO I=0,NT
	   
           	I3OUT(I,K)=I1(I,K)
	   	Q3OUT(I,K)=Q1(I,K)
           	U3OUT(I,K)=U1(I,K)
	   
	   ENDDO
	   
         ENDDO


C* Storage of the upward field of order IG = 1.
C-------------------------------------------
         DO K=1,NBMU

C* At TOA for the convergence test:
            I3(K)=I1(0,K)
	     Q3(K)=Q1(0,K)
            U3(K)=U1(0,K)
	     D1(K)=I1(0,K)
            E1(K)=Q1(0,K)
            F1(K)=U1(0,K)
	 
C* At the output level of the results:
	    DO I=0,NT	
	    
            	I3OUT(I,K)=I1(I,K)
	    	Q3OUT(I,K)=Q1(I,K)
            	U3OUT(I,K)=U1(I,K)
	    
	    ENDDO
	    
         ENDDO


C--------------------------------------------------------------
C
C   ******  LOOP ON THE SUCCESSIVE INTERACTIONS   ******
C
C--------------------------------------------------------------

         IG=1

  503    IG=IG+1
  
C* Stop test of the computations if the maximum order of scattering is reached
C-----------------------------------------------------------------------
         IF (IG.GT.IGMAX) GOTO 505


C* Computation of the source function of order IG: I2, Q2, U2
C----------------------------------------------------
         CALL SOS_FSOURCE_ORDREIG(IS,NBMU,NT,XDEL,YDEL,
     &                            BETA0,BETA2,GAMMA2,ALPHA2,
     &                            XPL,XRL,XTL,I1,Q1,U1,
     &                            BP,GR,GT,ARR,ART,ATT,GA,
     &                            I2,Q2,U2)
     
	
C* Boundary condition at the ground
C-------------------------
         DO 1081 K=1,NBMU
            I1(NT,K)=0.
            Q1(NT,K)=0.
            U1(NT,K)=0.
            XR(K)=0.
 1081    CONTINUE	    

	    
C* Integration of the reflection on a Lambertian ground of
C* the downward radiance of order IG-1.
C-------------------------------------------------------
         LSOL=0.
         DO 9876 J=1,NBMU
            LSOL=LSOL+GA(J)*I1(NT,-J)*RMU(J)            
 9876    CONTINUE
         LSOL=2*LSOL*RO

C* Lambertian reflection at the ground. Stored in XR.
C--------------------------------------------------
         IF (RO.EQ.0..OR.IS.NE.0) GOTO 2121
  
         DO 2122 J=1,NBMU
            I1(NT,J) = LSOL  
            XR(J) = LSOL
 2122    CONTINUE

C* Surface reflection (BRDF / BPDF)
C------------------------------------		
 2121   IF (IMAT_SURF.EQ.1) THEN
 
            DO 960 K=1,NBMU
               II2=0.
               QQ2=0.
               UU2=0.
               RRMU=2/RMU(K)
	 
               DO 961 J=1,NBMU
                  Z= GA(J)
                  XI1 = I1(NT,-J)
                  XQ1 = Q1(NT,-J)
                  XU1 = U1(NT,-J)

                  II2 = II2 +Z*( XI1*R11(J,K)+XQ1*R12(J,K)+XU1*R13(J,K))
                  QQ2 = QQ2 +Z*( XI1*R21(J,K)+XQ1*R22(J,K)+XU1*R23(J,K))
                  UU2 = UU2 +Z*( XI1*R31(J,K)+XQ1*R32(J,K)+XU1*R33(J,K))
		  
  961          CONTINUE
		
C* Upward radiance at ground level (order IG)
               I1(NT,K)=II2*RRMU +XR(K)
               Q1(NT,K)=QQ2*RRMU
               U1(NT,K)=UU2*RRMU
  960       CONTINUE

         ENDIF


C* Fresnel reflection on a flat sea
C------------------------------------
      IF (IFRESNEL.EQ.1) THEN
 
          DO 1510 K=1,NBMU
	  
	     I1(NT,K) = I1(NT,K) 
     &            + F11(K)*I1(NT,-K) + F12(K)*Q1(NT,-K)
     
             Q1(NT,K) = Q1(NT,K) 
     &            + F12(K)*I1(NT,-K) + F11(K)*Q1(NT,-K)
     
             U1(NT,K) = U1(NT,K) + F33(K)*U1(NT,-K) 
	     
1510 	  CONTINUE 

      ENDIF

C* Integration of the upward radiance field (mu > 0) from the ground to TOA
C* and of the downward radiance field (mu < 0) from TOA to the ground.
C-------------------------------------------------------------------
         CALL SOS_INTEGR_EPOPT(NBMU,RMU,NT,H,I2,Q2,U2,I1,Q1,U1)

C* Storage of the field of order IG.
C------------------------------
         DO K=-NBMU,-1   
C* At the ground for the convergence test:
            G1(K)=I1(NT,K)
            H1(K)=Q1(NT,K)
            P1(K)=U1(NT,K)  	
	
C* At the output level of the results:
	    DO I=0,NT
	    
	        G1OUT(I,K)=I1(I,K)
               H1OUT(I,K)=Q1(I,K)
               P1OUT(I,K)=U1(I,K)
	    
	    ENDDO
	    
         ENDDO

         DO K=1,NBMU     
C* At TOA for the convergence test:
	     G1(K)=I1(0,K)
            H1(K)=Q1(0,K)
            P1(K)=U1(0,K)
	
C* At the output level of the results:
	    DO I=0,NT

	       G1OUT(I,K)=I1(I,K)
	       H1OUT(I,K)=Q1(I,K)
              P1OUT(I,K)=U1(I,K)
	    
	    ENDDO
	    
         ENDDO
 

C*  Convergence test in geometric series and storage of the order ig-2
C-------------------------------------------------------------------------
         IF (IG.EQ.2) GO TO 506
      
         CALL SOS_PARAM_CONV(NBMU,A1,B1,C1,D1,E1,F1,G1,H1,P1,
     &                       I3,Q3,U3,Z1)

C* If the threshold CTE_PH_SEUIL_CV_SG is exceeded for at least one angle mu,
C* the convergence condition is not reached. The next scattering order
C* is computed.
         IF (Z1.GT.CTE_PH_SEUIL_CV_SG) GO TO 506


C* If the convergence regime is reached, the tail of the geometric
C* series is added to the sum of the previous scatterings.
C-------------------------------------------------------
         CALL SOS_AJOUT_QUEUE(NT,NBMU,D1,E1,F1,G1,H1,P1,
     &                        D1OUT,E1OUT,F1OUT,
     &                        G1OUT,H1OUT,P1OUT,
     &                        I3,Q3,U3,I3OUT,Q3OUT,U3OUT)

C* Trace file
C------------------------------------------------------  
         IF (TRACE) THEN
            WRITE(IDLOG,*,err=921) 
     &      '    Convergence in geometric series from IG = ',IG
            WRITE(IDLOG,*,err=921) ' '	   
         ENDIF

C* Exit of the loop on the scatterings.
C* Next: test on the order of the Fourier series decomposition.
C---------------------------------------------------------------------
         GO TO 505


C* Still in the loop on the scatterings.
C* Storage of the order IG-1 in A1, B1 and C1.
C* Storage of the order IG in D1, E1 and F1.
C* for the test when going to the next order IG+1.
C-------------------------------------------------------
  506    DO 86 K=-NBMU,NBMU
            A1(K)=D1(K)
            B1(K)=E1(K)
            C1(K)=F1(K)
            D1(K)=G1(K)
            E1(K)=H1(K)
            F1(K)=P1(K)
	    
	    DO I=0,NT
	    
               D1OUT(I,K)=G1OUT(I,K)
               E1OUT(I,K)=H1OUT(I,K)
               F1OUT(I,K)=P1OUT(I,K)
	    
	    ENDDO
	    
   86    CONTINUE

C* Computation of the sum of the scatterings
C------------------------------------
17       DO 400 J=1,NBMU

            I3(J)=I3(J)+I1(0,J)
            Q3(J)=Q3(J)+Q1(0,J)
            U3(J)=U3(J)+U1(0,J)
            I3(-J)=I3(-J)+I1(NT,-J)
            Q3(-J)=Q3(-J)+Q1(NT,-J)
            U3(-J)=U3(-J)+U1(NT,-J)

	    DO I=0,NT
	    
               I3OUT(I,J)=I3OUT(I,J)+I1(I,J)
               Q3OUT(I,J)=Q3OUT(I,J)+Q1(I,J)
               U3OUT(I,J)=U3OUT(I,J)+U1(I,J)
               I3OUT(I,-J)=I3OUT(I,-J)+I1(I,-J)
               Q3OUT(I,-J)=Q3OUT(I,-J)+Q1(I,-J)
               U3OUT(I,-J)=U3OUT(I,-J)+U1(I,-J)
	    
	    ENDDO
	    
  400    CONTINUE

C* Stop test on the scatterings when the terms
C* I, Q and U are too small
C-------------------------------------------------------
         CALL SOS_ARRET_DIFFUS_1(NBMU,NT,I1,Q1,U1,Z1)

         IF (Z1.GT.CTE_PH_SEUIL_VALDIF) goto 507

C* Trace file
C------------------------------------------------------  
         IF (TRACE) THEN
            WRITE(IDLOG,*,err=921) 
     &        '    End of the scattering loop due to the test',
     &        ' on I, Q, U values',
     &        ' (too low to continue the calculations) for IG = ',IG
            WRITE(IDLOG,*,err=921) ' '	   
         ENDIF

         GO TO 505

C* Stop test on the scatterings when the next terms
C* become very small
C-------------------------------------------------------
  507    CALL SOS_ARRET_DIFFUS_2(NBMU,NT,I1,Q1,U1,I3,Q3,U3,Z1)
   
         IF (Z1.GT.CTE_PH_SEUIL_SUMDIF) goto 508
      
C* Trace file
C------------------------------------------------------  
         IF (TRACE) THEN
            WRITE(IDLOG,*,err=921) 
     &        '    End of the scattering loop due to the test',
     &        ' comparing radiance of order IG to the accumulation',
     &        ' of previous orders for IG = ',IG
            WRITE(IDLOG,*,err=921) ' '	   
         ENDIF
	 
	 
         GO TO 505

C* Stop test on the maximum number of scatterings.
C---------------------------------------------------
  508    IF(IG.LT.IGMAX) GO TO 503
  
C* Trace file
C------------------------------------------------------  
         IF (TRACE) THEN
            WRITE(IDLOG,*,err=921) 
     &       '    End of the scattering loop because of reaching',
     &       '    the maximal order : IGMAX = ',IGMAX
            WRITE(IDLOG,*,err=921) ' '	   
         ENDIF
	 
  505    CONTINUE

C* Removal of the direct reflection term on the surface (BRDF / BPDF).
C-----------------------------------------------------------------------	
        IF (IMAT_SURF.EQ.1) THEN
	 
            DO 5200 J=1,NBMU
	    
               I3(J)=I3(J)-RII(J)
               Q3(J)=Q3(J)-RQQ(J)
               U3(J)=U3(J)-RUU(J)
	       
	       DO I=0,NT
	       
                  I3OUT(I,J)=I3OUT(I,J)-RIIOUT(I,J)
                  Q3OUT(I,J)=Q3OUT(I,J)-RQQOUT(I,J)
                  U3OUT(I,J)=U3OUT(I,J)-RUUOUT(I,J)
	       
	       ENDDO
	       
 5200       CONTINUE

         ENDIF

C* Computation of the diffuse irradiance at the top and bottom of the atmosphere,
C* normalised by the incident solar flux at TOA.
C*  EMOINS ==> diffuse transmission td(thetas)
C*        if the surface is not reflecting.
C*  EPLUS ==> includes the different surface - atmosphere reflections.
C-------------------------------------------------------------------
         IF(IS.EQ.0) THEN
           EMOINS=0.
           EPLUS=0.
           DO J=1,NBMU
             EMOINS = EMOINS + RMU(J)*GA(J)*I3(-J)
             EPLUS  = EPLUS  + RMU(J)*GA(J)*I3(J)
           ENDDO
	    EMOINS = -EMOINS*2/TAB
	    EPLUS  = -EPLUS *2/TAB
         ENDIF	

C* Computation of the functions of mu I4, Q4, U4, I5, Q5 and U5:
C-------------------------------------------------------
         COEF=2.
         IF (IS.EQ.0) COEF=1.

         DO 67 J=-NBMU,NBMU
            IF (J.EQ.0) GOTO 67
	
	     I4(J)=I4(J)+COEF*I3(J)
            Q4(J)=Q4(J)+COEF*Q3(J)
            U4(J)=U4(J)+COEF*U3(J)

	     I5(J)=I5(J)+COEF*I3(J)*SIGN
            Q5(J)=Q5(J)+COEF*Q3(J)*SIGN
            U5(J)=U5(J)+COEF*U3(J)*SIGN
67       CONTINUE


C* Definition of the profile levels for the upward and downward fields.
C  Default case: ZOUT=-1 => TOA and surface levels
C  Otherwise computation by interpolation of the fluxes at the altitude ZOUT
C-----------------------------------------------------------------------

C* Default case.
C---------------------------

      IF (ZOUT.EQ.-1) THEN
      
C* For the downward field.
C
	 NDOWN_OUT = NT
         DO K=-NBMU,-1

            I3ZOUT(K)=I3OUT(NDOWN_OUT,K)
	     Q3ZOUT(K)=Q3OUT(NDOWN_OUT,K)
            U3ZOUT(K)=U3OUT(NDOWN_OUT,K)
	   
         ENDDO
      
C* For the upward field.
C
	 NUP_OUT = 0
         DO K=1,NBMU

            I3ZOUT(K)=I3OUT(NUP_OUT,K)
	     Q3ZOUT(K)=Q3OUT(NUP_OUT,K)
            U3ZOUT(K)=U3OUT(NUP_OUT,K)
	   
         ENDDO

C* Computation at ZOUT by interpolation.
C-----------------------------------
 
      ELSE
      
c Search for the altitude level below ZOUT
	    J=1
	    DO WHILE(ZOUT.LT.ZPROF(J))
		   J=J+1
	    ENDDO

C Position of the altitude with respect to the previous and next levels
	    ZZ = (ZOUT-ZPROF(J-1))/(ZPROF(J)-ZPROF(J-1))

	    
C* Computation of the downward and upward fields
C

         DO K=-NBMU,NBMU

            I3ZOUT(K)=(1-ZZ)*I3OUT(J-1,K)+ZZ*I3OUT(J,K)
	     Q3ZOUT(K)=(1-ZZ)*Q3OUT(J-1,K)+ZZ*Q3OUT(J,K)
            U3ZOUT(K)=(1-ZZ)*U3OUT(J-1,K)+ZZ*U3OUT(J,K)
	   
         ENDDO
	 	          
      ENDIF

C* Trace file
C------------------------------------------------------  
         IF (TRACE) THEN
            WRITE(IDLOG,*,err=921)	    
     &        '    Scattering accumulation : IS = ',IS
	     IF (ZOUT.EQ.-1) THEN
	       WRITE(IDLOG,*,err=921)	    
     &        '       Upward radiance at level   : ',NUP_OUT
              WRITE(IDLOG,*,err=921) 
     &        '       Downward radiance at level : ',NDOWN_OUT
     	     ELSE
	       WRITE(IDLOG,*,err=921)	    
     &        '       Altitude of output radiance : ',ZOUT	    
	     ENDIF   
            WRITE(IDLOG,*,err=921) '  '	   	    
            DO K=1,NBMU
              WRITE(IDLOG,20,err=921) RMU(K),
     &                             '  ID ',I3ZOUT(-K),
     &                             '  QD ',Q3ZOUT(-K),
     &                             '  UD ',U3ZOUT(-K),	      
     &                             '  IM ',I3ZOUT(K),
     &                             '  QM ',Q3ZOUT(K),
     &                             '  UM ',U3ZOUT(K)     
 	     ENDDO
            WRITE(IDLOG,*,err=921) '  '
	     WRITE(IDLOG,*,err=921) ' '	   	    	   
         ENDIF


C* Recording of the radiance arrays at order IS of
C* the Fourier series decomposition, for the profile level
C* requested by the user.
C* Record bounds of the arrays limited to the useful dimension
C*    dimension -NBMU:NBMU
C--------------------------------------------------------------
      IF (FICOS.NE.'NO_OUTPUT') THEN
         WRITE(10,err=923)(Q3ZOUT(K),K=-NBMU,NBMU),
     & 		     (U3ZOUT(K),K=-NBMU,NBMU),
     &			     (I3ZOUT(K),K=-NBMU,NBMU)
      ENDIF      


C* Stop test on the Fourier series decomposition
C-------------------------------------------------------
         CALL SOS_ARRET_FOURIER(NBMU,I3,Q3,U3,I4,Q4,U4,I5,Q5,U5,Z1)
 
C* If at least one term is above the threshold, the next
C* order IS is computed.
C-------------------------------------------------------
         IF (Z1.GT.CTE_PH_SEUIL_SF) GO TO 24
	 
C* Exit of the loop on IS
C---------------------------
         GOTO 243

C* END OF LOOP ON IS
C----------------------
   24 CONTINUE
   
  243 CONTINUE   
   

C* Trace file
C------------------------------------------------------  
         IF (TRACE) THEN
            WRITE(IDLOG,*,err=921) "--> End of the Fourier series"//
     & 				" expansion for the order : ", IS
            WRITE(IDLOG,*,err=921) ' '
	     WRITE(IDLOG,*,err=921) ' '
	     WRITE(IDLOG,*,err=921)  
     &      ' FLUX (normalized by TOA solar irradiance):'
            WRITE(IDLOG,*,err=921) 
     &      '      Downward diffuse flux at ground level ', EMOINS
            WRITE(IDLOG,*,err=921) 
     &      '        --> including multiple surface / atm. reflections'
            WRITE(IDLOG,*,err=921)
     &      '        --> = Diffuse transmittance td(thetas) only'//
     &      ' for a black surface'
            WRITE(IDLOG,*,err=921) 
     &            '      Upward diffuse flux a at TOA ',EPLUS   
            WRITE(IDLOG,*,err=921) 
     &      '        --> including multiple surface / atm interactions'
            WRITE(IDLOG,*,err=921) ' '	   
         ENDIF

C* Nominal end
C-------------------
      goto 9999
      
C* Error case: return of status -1 to the calling program
C----------------------------------------------------------
  910 WRITE(6,*) '  ERROR on SURFACE file opening for SOS_OS'
      IER=-1  
      goto 9999 
  913 WRITE(6,*) '  ERROR on SOS binary result file opening for SOS_OS'
      IER=-1
      goto 9999    
  920 WRITE(6,*) '  ERROR on SURFACE file reading for SOS_OS'
      IER=-1
      goto 9999
  921 WRITE(6,*) '  ERROR on logfile writing for SOS_OS'
      IER=-1
      goto 9999
  922 WRITE(6,*) '  ERROR on angles file reading for SOS_OS'
      IER=-1
      goto 9999
  923 WRITE(6,*) '  ERROR on SOS binary result file writing for SOS_OS'
      IER=-1
      goto 9999
  996 WRITE(6,*) '  Inconsistency output altitude in the ',
     &           '  atmospheric profile'
      IER=-1
      goto 9999

C* Special case
C--------------------     
  995 WRITE(6,*) '  SOS computation are stopped ',
     &           '  because of a limb incidence'
      goto 9999

      
 9999 CONTINUE

C* File closing
C-------------------
      CLOSE(21)
      IF (IMAT_SURF.EQ.1) CLOSE(12)
      IF (FICOS.NE.'NO_OUTPUT') CLOSE(10)

C* Format 
C----------  
   15 FORMAT(2(E21.14))
   10 FORMAT(F9.6,A5,E12.5,A5,E12.5,A5,E12.5) 
   20 FORMAT(F9.6,6(A5,E12.5)) 
   25 FORMAT(1X,I3,4(2X,E13.5))
   35 FORMAT(3X,A1,3X,A6,11X,A7,10X,A8,9X,A7)
   40 FORMAT(A17,F6.2,3(A5,E13.7))          
      RETURN
      END  	!END OF PROCEDURE SOS_OS




C*===========================================================================
C PROCEDURE: SOS_MAT_FRESNEL_PLAN_REFL
C ==========
C      This procedure computes the elements of the Fresnel reflection
C      matrix for reflection angles equivalent to the Gauss
C      angles.
C
C      For an incidence -mu, the Fresnel reflection matrix is:
C           F11(mu)    F12(mu)     0
C           F12(mu)    F11(mu)     0
C              0         0     F33(mu)
C
C
C Description of parameters
C -------------------------- 
C
C     NBMU (I4) : (I) Number of actually useful (positive) angles
C
C      RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)	 (I) Array of the Gauss mu.
C
C     IND_SURF (double)  (I) Refractive index surface / atmosphere.
C
C     IPOLAR  (I4)   (I) : Option for simulation without polarisation.
C      			   0 : computations without polarisation.
C			          1 : normal computations with polarisation.
C	
C     F11(0:CTE_OS_NBMU_MAX)  (double)  (O)
C     F12(0:CTE_OS_NBMU_MAX)  (double)  (O)
C     F33(0:CTE_OS_NBMU_MAX)  (double)  (O)
C           Elements of the Fresnel reflection matrix.
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C None
C
C==============================================================================
      SUBROUTINE SOS_MAT_FRESNEL_PLAN_REFL(NBMU,RMU,IND_SURF,IPOLAR,
     &                                     F11,F12,F33)
     
      IMPLICIT NONE


C* Definition of variables
C*-----------------------------------------------------------------------

      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
                                               !  Array of the Gauss mu.

      DOUBLE PRECISION IND_SURF  ! Refractive index surface / atmosphere.
      
      DOUBLE PRECISION F11(0:CTE_OS_NBMU_MAX)  ! Element of the Fresnel matrix.
      DOUBLE PRECISION F12(0:CTE_OS_NBMU_MAX)  ! Element of the Fresnel matrix.
      DOUBLE PRECISION F33(0:CTE_OS_NBMU_MAX)  ! Element of the Fresnel matrix.


      DOUBLE PRECISION MU   ! Value of a Gauss mu (reflection direction).
      DOUBLE PRECISION MU2  ! Square of mu.
      DOUBLE PRECISION IND2 ! Square of IND_SURF.
      DOUBLE PRECISION X    
      DOUBLE PRECISION RL   ! Parallel reflection coefficient.
      DOUBLE PRECISION RR   ! Perpendicular reflection coefficient.

      INTEGER*4 NBMU    ! Number of actual angles of the simulations
      INTEGER*4 IPOLAR  ! Option for simulation with / without polarisation.
			      
      INTEGER J  ! Index

C* Processing
C*------------

      DO 10 J=0,NBMU

C* Case J=0 : reflection of the direct solar beam of incidence rmu(0) < 0.
C* Case J>0 : reflection of a ray of incidence -rmu(j) < 0.
	 IF(J.EQ.0) THEN
	    MU = -RMU(0)
	 ELSE
	    MU = RMU(J)
	 ENDIF   
      
         IND2 = IND_SURF*IND_SURF
	 MU2  = MU*MU	 
	 X = DSQRT(IND2-1.D+00+MU2)
	 RL = (IND2*MU-X)/(IND2*MU+X)
	 RR = (MU-X)/(MU+X)
	 
	 F11(J) = (RL*RL+RR*RR)/2.
	 
	 IF(IPOLAR.EQ.1) THEN
	    F12(J) = (RL*RL-RR*RR)/2.
	    F33(J) =  RL*RR
	 ELSE
            F12(J) = 0.D+00
	    F33(J) = 0.D+00
	 ENDIF
	  
	    	     
10    CONTINUE

      END      !END OF PROCEDURE SOS_MAT_FRESNEL_PLAN_REFL
       
       
       
       
C*===========================================================================
C PROCEDURE: SOS_NOYAUX
C ==========
C      This procedure computes the elements of the aerosol phase
C      matrix at order IS of the Fourier series expansion.
C      The computations use the Legendre polynomials and functions.
C      These functions are output for an order L=2. They are used
C      to compute the elements of the molecular phase matrix.
C
C      The following matrix P(mu,mup) is:
C        BP(mu,mup)    GR(mu,mup)    -GT(mu,mup)
C        GR(mup,mu)    ARR(mu,mup)   -ART(mup,mu)
C       -GT(mup,mu)   -ART(mu,mup)    ATT(mu,mup)
C
C      BP(mu,mup), GR(mu,mup), GT(mu,mup), ARR(mu,mup),
C      ART(mu,mup) and ATT(mu,mup) are computed for -1 < mu < 1 and for -1 < mup < 1.
C
C
C Description of parameters
C -------------------------- 
C
C     IS  (I4)  (I)   Order of the Fourier series decomposition.
C     NBMU (I4) : (I) Number of actually useful (positive) angles
C      RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)	 (I) Array of the Gauss mu.
C     OS_NB (I4) : (I) Max order of the expansions of the phase functions
C				       in Legendre polynomials to be used.
C
C     ALPHA(0:CTE_OS_NB_MAX) (double)  (I)    Parameters alpha(k).
C     BETA(0:CTE_OS_NB_MAX)  (double)  (I)    Parameters beta(k).
C     GAMMA(0:CTE_OS_NB_MAX) (double)  (I)    Parameters gamma(k).
C     ZETA(0:CTE_OS_NB_MAX)  (double)  (I)    Parameters zeta(k).
C
C     XPL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (double)  (O)
C      		        Legendre function PSL at order L=2
C   		        and S=IS versus the angle mu.
C     XRL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (double)  (O)
C      		        Function RSL at order L=2 and S=IS.
C     XTL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (double)  (O)
C      		        Function TSL at order L=2 and S=IS
C	             Legendre functions and derived functions for the computation
C	             of the phase matrix.
C
C	
C     BP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (O)
C                                                        Element P11.
C     GR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (O)
C                                                        Element P12 or P21.
C     GT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (O)
C                                                        Element -P13 or -P31.
C     ARR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (double)  (O)
C                                                        Element P22.
C     ART(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (double)  (O)
C                                                        Element -P23 or -P32.
C     ATT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (double)  (O)
C                                                        Element P33.
C             Element Pij of the aerosol phase matrix for an order
C             IS of the Fourier series decomposition.
C             With mup the incidence direction ( -1 < mup < 1)
C             and mu the scattering one ( -1 < mu < 1),
C             the scattering is expressed by Pij(mu,mup).
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C None
C
C==============================================================================
      SUBROUTINE SOS_NOYAUX(IS,NBMU,RMU,OS_NB,ALPHA,BETA,GAMMA,ZETA,
     &                  XPL,XRL,XTL,BP,GR,GT,ARR,ART,ATT)
     
      IMPLICIT NONE


C* Definition of variables
C*-----------------------------------------------------------------------

      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				!  Array of the Gauss mu.

      DOUBLE PRECISION ALPHA(0:CTE_OS_NB_MAX)	! Parameters alpha(k).
      DOUBLE PRECISION BETA(0:CTE_OS_NB_MAX)	! Parameters beta(k).
      DOUBLE PRECISION GAMMA(0:CTE_OS_NB_MAX)	! Parameters gamma(k).
      DOUBLE PRECISION ZETA(0:CTE_OS_NB_MAX)	! Parameters zeta(k).

	! Legendre functions and derived functions for the computation
	! of the phase matrix.
      DOUBLE PRECISION XPL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				!   Legendre function PSL at order L=2
      				!   and S=IS versus the angle mu.
      DOUBLE PRECISION XRL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				!   Function RSL at order L=2 and S=IS.
      DOUBLE PRECISION XTL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				!   Function TSL at order L=2 and S=IS
							
	! Element Pij of the aerosol phase matrix for an order
	! IS of the Fourier series decomposition.
	! With mup the incidence direction (-1 < mup < 1)
	! and mu the scattering one (-1 < mu < 1),
	! the scattering is expressed by Pij(mu,mup).
      DOUBLE PRECISION 
     &    BP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &       -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  ! Element P11.
      DOUBLE PRECISION 
     &    GR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &       -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  ! P12 or P21.
      DOUBLE PRECISION 
     &    GT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &       -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  ! -P13 or -P31.
      DOUBLE PRECISION 
     &    ARR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &        -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! P22.
      DOUBLE PRECISION 
     &    ART(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &        -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! -P23 or -P32.
      DOUBLE PRECISION 
     &    ATT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &        -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! P33.



      DOUBLE PRECISION 
     &    PSL(-1:CTE_OS_NB_MAX,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
      				! Legendre function PSL for the order S=IS
      				! versus the order L and the angle mu.
      DOUBLE PRECISION 
     &    RSL(-1:CTE_OS_NB_MAX,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Function RSL.
      DOUBLE PRECISION 
     &    TSL(-1:CTE_OS_NB_MAX,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				  ! Function TSL.

      DOUBLE PRECISION C	! Single value of the Gauss mu.

      DOUBLE PRECISION A,B,D,E,F     ! Coefficients for the recurrence
      				     ! on PSL, RSL and TSL.

      DOUBLE PRECISION X	! Intermediate computation variable.
      DOUBLE PRECISION XX	! Intermediate computation variable.
      DOUBLE PRECISION YY	! Intermediate computation variable.

      DOUBLE PRECISION SBP	! Summation variable from L = IS to OS_NB for BP.
      DOUBLE PRECISION SATT	! Summation variable from L = IS to OS_NB for SATT.
      DOUBLE PRECISION SARR	! Summation variable from L = IS to OS_NB for SARR.
      DOUBLE PRECISION SGR	! Summation variable from L = IS to OS_NB for SGR.
      DOUBLE PRECISION SGT	! Summation variable from L = IS to OS_NB for SGT.
      DOUBLE PRECISION SART	! Summation variable from L = IS to OS_NB for SART.
      DOUBLE PRECISION R1		
      DOUBLE PRECISION R2	

      DOUBLE PRECISION RAC3	! Square root of 3.
      DOUBLE PRECISION X26	! 2*sqrt(6)

      INTEGER*4 OS_NB 	! Max order of the expansions of the phase functions
			       ! in Legendre polynomials.
      INTEGER IS	   ! Order of the Fourier series decomposition.

      INTEGER*4 NBMU    ! Number of actual angles of the simulations
      INTEGER*4 I
      INTEGER*4 K
      INTEGER*4 IG	! Variable equal to 1 or -1.
      INTEGER*4 L
      INTEGER*4 LP	! Value L+1
      INTEGER*4 LM	! Value L-1
      INTEGER*4 J	! Index of the array RMU.


C* Initialisation
C----------------
      RAC3=DSQRT(3.D+00)
      X26=2.*DSQRT(6.D+00)

C-----------------------------------------------------
C* Initialisation of the functions of order IS for the
C* orders L = 0, 1 and 2.
C-----------------------------------------------------

C* Initialisation for the order IS = 0
C-------------------------------------
      IF (IS.NE.0) GO TO 700
      
      DO 25 J=0,NBMU
         C=RMU(J)
         PSL(0,-J)=1.
         PSL(0,J)=1.
         PSL(1,J)=C
         PSL(1,-J)=-C
         X=(3.*C*C-1.)*0.5
         PSL(2,-J)=X
         PSL(2,J)=X
         RSL(1,J)=0.
         RSL(1,-J)=0.
         X=3.*(1.-C*C)/X26
         RSL(2,-J)=X
         RSL(2,J)=X
         TSL(1,J)=0.
         TSL(1,-J)=0.
         TSL(2,J)=0.
         TSL(2,-J)=0.
   25 CONTINUE
   
      PSL(1,0)=RMU(0)
      RSL(1,0)=0.
      
      GOTO 501

C* Initialisation for the order IS = 1
C-------------------------------------
  700 IF (IS.NE.1) GO TO 701
  
      DO 26 J=0,NBMU
         C=RMU(J)
         X=1.-C*C
         PSL(0,J)=0.
         PSL(0,-J)=0.
         PSL(1,-J)=DSQRT(X*0.5)
         PSL(1,J)=DSQRT(X*0.5)
         PSL(2,J)=C*PSL(1,J)*RAC3
         PSL(2,-J)=-PSL(2,J)
         RSL(1,-J)=0.
         RSL(1,J)=0.
         RSL(2,J)=-C*DSQRT(X)*0.5
         RSL(2,-J)=-RSL(2,J)
         TSL(1,-J)=0.
         TSL(1,J)=0.
         TSL(2,J)=-DSQRT(X)*0.5
         TSL(2,-J)=-DSQRT(X)*0.5
   26 CONTINUE
   
      PSL(2,0)=-PSL(2,0)
      RSL(2,0)=-RSL(2,0)
      RSL(1,0)=0.
      TSL(1,0)=0.
      
      GOTO 501

C* Initialisation for the order IS >= 2
C----------------------------------------
  701 A=1.
  
      DO 27 I=1,IS
        X=I
        A=A*DSQRT((I+IS)/X)*0.5
 27   CONTINUE
 
      B=A*DSQRT(IS/(IS+1.D+00))*DSQRT((IS-1.D+00)/(IS+2.)) 
      
      DO 28 J=0,NBMU
         C=RMU(J)
         XX=1.-C*C
         YY=IS*0.5-1.
         PSL(IS-1,J)=0.
         RSL(IS-1,J)=0.
         TSL(IS-1,J)=0.
         X=A*XX**(IS*0.5)
         PSL(IS,-J)=X
         PSL(IS,J)=X
         X=B*(1.+C*C)*XX**YY
         RSL(IS,-J)=X
         RSL(IS,J)=X
         X=2.*B*C*XX**YY
         TSL(IS,-J)=-X
         TSL(IS,J)=X
   28 CONTINUE


C* Computation of the functions PSL(L,mu), RSL(L,mu) and TSL(L,mu)
C* by recurrence on the order L for fixed IS.
C------------------------------------------------------
  501 K=2

      IF (IS.GT.2) K=IS
      
      IF (K.EQ.OS_NB) GOTO 502
      
      IG=-1
      IF (IS.EQ.1) IG=1

      DO 30 L=K,(OS_NB-1)
         LP=L+1
         LM=L-1

C* Coefficients for the recurrence on PSL
         A=(2*L+1.)/DSQRT((L+IS+1.D+00)*(L-IS+1.))
         B=DSQRT(DFLOAT((L+IS)*(L-IS)))/(2.*L+1.)

C* Coefficients for the recurrence on RSL and TSL
         D=(L+1.)*(2*L+1.)/DSQRT((L+3.D+00)*(L-1.)*(L+IS+1.)  
     &     *(L-IS+1.))
         E=DSQRT((L+2.D+00)*(L-2.)*(L+IS)*(L-IS))/(L*(2.*L+1.))
         F=2.*IS/(L*(L+1.))

C* Application of the recurrence for each Gauss mu
C------------------------------------------------------
         DO 31 J=0,NBMU
            C=RMU(J)
            X=A*(C*PSL(L,J)-B*PSL(LM,J))
            PSL(LP,J)=X
            X=D*(C*RSL(L,J)-F*TSL(L,J)-E*RSL(LM,J))
            RSL(LP,J)=X
            X=D*(C*TSL(L,J)-F*RSL(L,J)-E*TSL(LM,J))
            TSL(LP,J)=X
	    
            IF(J.EQ.0) GO TO 31
	    
            PSL(LP,-J)=IG*PSL(LP,J)
            RSL(LP,-J)=IG*RSL(LP,J)
            TSL(LP,-J)=-IG*TSL(LP,J)
   31    CONTINUE

         IG=-IG
   30 CONTINUE

  502 CONTINUE


C* Storage of PSL(L=2,mu), RSL(L=2,mu) and TSL(L=2,mu)
C------------------------------------------------------
      DO 1005 J=-NBMU,NBMU
         XPL(J)=PSL(2,J)
         XRL(J)=RSL(2,J)
         XTL(J)=TSL(2,J)
 1005 CONTINUE


C------------------------------------------------------
C* Computation of the elements of the aerosol phase matrix
C------------------------------------------------------

C* Loops on the incidence directions (K = mup,
C* with -1 < mup < 1) and emergence directions (J = mu > 0)
C------------------------------------------------------
      DO 321 J=-NBMU,NBMU
         DO 32 K=-NBMU,NBMU
            SBP=0.
            SATT=0.
            SARR=0.
            SGR=0.
            SGT=0.
            SART=0.
	    
            IF(IS.GT.OS_NB) GOTO 1

C* Summation loop on L = IS to OS_NB.
C-----------------------------------------------
            DO 33 L=IS,OS_NB
               R1=TSL(L,J)*TSL(L,K)
               R2=RSL(L,J)*RSL(L,K)
               SBP=SBP+BETA(L)*PSL(L,J)*PSL(L,K)
               SATT=SATT+ALPHA(L)*R1+ZETA(L)*R2
               SARR=SARR+ZETA(L)*R1+ALPHA(L)*R2
               SGR=SGR+GAMMA(L)*PSL(L,J)*RSL(L,K)
               SGT=SGT+GAMMA(L)*PSL(L,J)*TSL(L,K)
               SART=SART+ALPHA(L)*RSL(L,K)*TSL(L,J)+
     &              ZETA(L)*RSL(L,J)*TSL(L,K)
  33        CONTINUE !Loop on L
  
 1          CONTINUE

            BP(J,K)=SBP
            ATT(J,K)=SATT
            ARR(J,K)=SARR
            GR(J,K)=SGR
            GT(J,K)=SGT
            ART(J,K)=SART
   32    CONTINUE !Loop on K
  321 CONTINUE !Loop on J
   
      RETURN
      END  	!END OF PROCEDURE SOS_NOYAUX









C*===========================================================================
C PROCEDURE: SOS_INTEGR_EPOPT
C ==========
C      This procedure performs the integration over the optical thickness
C      of the radiation field:
C        - upward field    : Integration from the ground to TOA.
C        - downward field  : Integration from TOA to the ground.
C
C      The integration is made step by step from the starting level
C      to the final level.
C      This implies, in solving the Radiative Transfer Equation,
C      integrations over the optical thickness between two consecutive
C      levels. These integrations involve the scattering source
C      functions (for I, Q and U - versus the optical thickness level
C      and the propagation direction) for the current orders
C      of scattering and of Fourier series decomposition.
C      The integration between two levels is made by linearisation of the
C      source function over the layer.
C
C       The atmospheric layers are assumed homogeneous.
C
C
C
C Description of parameters
C -------------------------- 
C
C       NBMU (I4)    (I)  Number of actually useful (positive) angles
C      RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)	 (I) Array of the Gauss mu.
C     H(0:CTE_OS_NT) (double) (I)  Optical thickness at each level.
C
C     I2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double) (I)
C                                                  Source function for I
C     Q2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double) (I)
C                                                  Source function for Q
C     U2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double) (I)
C                                                  Source function for U
C                          Source function vector versus the level of the
C                          profile and the direction.
C    	
C     I1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double) (O) Field I.
C     Q1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double) (O) Field Q.
C     U1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double) (O) Field U.
C                          Radiance field versus the level
C                          of the profile and the direction.
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C None
C
C==============================================================================
      SUBROUTINE SOS_INTEGR_EPOPT(NBMU,RMU,NT,H,I2,Q2,U2,I1,Q1,U1)
      
      IMPLICIT NONE
      

C* Definition of variables
C*------------------------------------------------------------------

C* Propagation directions:
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				!  Array of the Gauss mu.

C* Optical thickness profile:
      DOUBLE PRECISION H(0:CTE_OS_NT)	!  Optical thickness for each level.
     	
C* Source function vector versus the level of the profile and the direction.
      DOUBLE PRECISION I2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      					!  Source function for I.
      DOUBLE PRECISION Q2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      					!  Source function for Q.
      DOUBLE PRECISION U2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      					!  Source function for U.

C* Radiance field versus the level of the profile and the direction.
      DOUBLE PRECISION I1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	!  Field I.
      DOUBLE PRECISION Q1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	!  Field Q.
      DOUBLE PRECISION U1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	!  Field U.

      DOUBLE PRECISION RMUK 	! Value K of the array of the Gauss mu

      DOUBLE PRECISION ZI1	! Intermediate value for the computation of I1
      				! at order IG.
      DOUBLE PRECISION ZQ1	! Intermediate value for the computation of Q1.
      DOUBLE PRECISION ZU1	! Intermediate value for the computation of U1.

      DOUBLE PRECISION DTAU   	! Optical thickness of a layer.

      DOUBLE PRECISION ATT	! Attenuation between two consecutive levels.

      DOUBLE PRECISION MATT	       ! 1 - att.

      DOUBLE PRECISION ATTDTAU	! att * dtau

      DOUBLE PRECISION A,B	! Linear interpolation coefficients

      INTEGER*4 NT		! Number of layers of the atmospheric profile
      
      INTEGER*4 I,JJ		! Level indices.
      INTEGER*4 K		! Direction index.
      INTEGER*4 NBMU        ! Number of actual angles of the simulations
		
C------------------------------------------------------	
C* Integration of the upward radiance field (mu > 0) from the ground to TOA.
C*
C* Step-by-step integration from the ground to TOA by linearisation
C* of the source function between two successive levels.
C------------------------------------------------------
      DO 107 K=1,NBMU

         RMUK=RMU(K)

C* Initialisation with the value at the ground.
         ZI1=I1(NT,K)
         ZQ1=Q1(NT,K)
         ZU1=U1(NT,K)
		
         DO 108 I=NT-1,0,-1
            JJ=I+1
            DTAU = H(JJ)-H(I)
            ATT = DEXP(-DTAU/RMUK)
            MATT = 1.D+00 - ATT
            ATTDTAU =  ATT * DTAU

            B = I2(I,K)
            A = (I2(JJ,K) - B) / DTAU
            ZI1 = ZI1*ATT + MATT*(A*RMUK+B) - A*ATTDTAU
            I1(I,K)=ZI1

            B = Q2(I,K)
            A = (Q2(JJ,K) - B) / DTAU
            ZQ1 = ZQ1*ATT + MATT*(A*RMUK+B) - A*ATTDTAU
            Q1(I,K)=ZQ1

            B = U2(I,K)
            A = (U2(JJ,K) - B) / DTAU
            ZU1 = ZU1*ATT + MATT*(A*RMUK+B) - A*ATTDTAU
            U1(I,K)=ZU1
108      CONTINUE     !Loops on I
107   CONTINUE        !Loops on K


C------------------------------------------------------
C* Integration of the downward radiance field (mu < 0) from TOA to the ground
C*
C* Step-by-step integration from TOA to the ground by linearisation
C* of the source function between two successive levels.
C* Initialisation with a null downward diffuse field at TOA.
C------------------------------------------------------
      DO 109 K=-NBMU,-1

         RMUK=RMU(K)

C* Initialisation with the value at TOA.
         I1(0,K)=0.
         Q1(0,K)=0.
         U1(0,K)=0.
         ZI1=0.D+00
         ZQ1=0.D+00
         ZU1=0.D+00

         DO 110 I=1,NT
            JJ=I-1
            DTAU = H(I)-H(JJ)
            ATT = DEXP(DTAU/RMUK)
            MATT = 1.D+00 - ATT
            ATTDTAU =  ATT * DTAU

            B = I2(I,K)
            A = (B - I2(JJ,K)) / DTAU
            ZI1 = ZI1*ATT + MATT*(A*RMUK+B) + A*ATTDTAU
            I1(I,K)=ZI1

            B = Q2(I,K)
            A = (B - Q2(JJ,K)) / DTAU
            ZQ1 = ZQ1*ATT + MATT*(A*RMUK+B) + A*ATTDTAU
            Q1(I,K)=ZQ1

            B = U2(I,K)
            A = (B- U2(JJ,K)) / DTAU
            ZU1 = ZU1*ATT + MATT*(A*RMUK+B) + A*ATTDTAU
            U1(I,K)=ZU1
  110    CONTINUE  !Loops on I
  109 CONTINUE  !Loops on K

      RETURN
      END  	!END OF PROCEDURE SOS_INTEGR_EPOPT






C*============================================================================
C PROCEDURE: SOS_FSOURCE_ORDRE1
C ==========
C      This procedure computes the source function vector of the
C      primary scattering, from the incident solar radiation field and the
C      scattering properties of the atmosphere (phase function).
C      
C      Computation for each level of the profile and for each propagation
C      direction.
C
C
C Description of parameters
C -------------------------- 
C
C       IS   (I4)    (I)  Current order of the Fourier series expansion.
C       NBMU (I4)    (I)  Number of actually useful (positive) angles
C       JK   (I4)    (I)  Storage number of mus in the array RMU
C                        (0 by default).
C       XDEL(0:CTE_OS_NT)  (double)  (I)  Percentage of aerosols per level.
C       YDEL(0:CTE_OS_NT)  (double)  (I)  Percentage of molecules per level.
C
C       BETA0	  (double)  (I)
C       BETA2   (double)  (I) 	.
C       GAMMA2  (double)  (I)
C                   Coefficients of the expansion of the molecular
C                   phase matrix in Legendre functions
C
C       XPL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C      		     Legendre function PSL at order L=2 and S=IS
C      		     versus the angle mu.
C       XRL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C     		     Function RSL at order L=2 and S=IS versus mu.
C       XTL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C     		     Function TSL at order L=2 and S=IS versus mu.
C	
C       BP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C                                                        Element P11.
C       GR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C                                                        Element P12 or P21.
C       GT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C                                                        Element -P13 or -P31.
C             Element Pij of the aerosol phase matrix for an order
C             IS of the Fourier series decomposition.
C             With mup the incidence direction (-1 < mup < 1)
C             and mu the scattering one (-1 < mu < 1),
C             the scattering is expressed by Pij(mu,mup).
C
C       CH(0:CTE_OS_NT)  (double)  (I) 	Attenuation of the direct solar
C                                       beam/4.
C
C       I2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (O)
C       Q2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (O)
C       U2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (O)
C	         Source function vector I,Q,U versus the level of the profile
C	         and the direction.
C 
C
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C None
C
C==============================================================================
      SUBROUTINE SOS_FSOURCE_ORDRE1(IS,NBMU,NT,JK,XDEL,YDEL,
     &                              BETA0,BETA2,GAMMA2,
     &                              XPL,XRL,XTL,BP,GR,GT,CH,I2,Q2,U2)
      
      IMPLICIT NONE
      

C* Definition of variables
C*------------------------------------------------------------------

      DOUBLE PRECISION XDEL(0:CTE_OS_NT)   !Percentage of aerosols per level.
      DOUBLE PRECISION YDEL(0:CTE_OS_NT)   !Percentage of molecules per level.

      DOUBLE PRECISION BETA0	! Coefficients of the expansion of the molecular
      DOUBLE PRECISION BETA2	! phase matrix in Legendre functions.
      DOUBLE PRECISION GAMMA2

      DOUBLE PRECISION XPL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      			      ! Legendre function PSL at order L=2 and S=IS
      			      ! versus the angle mu.
      		
      DOUBLE PRECISION XRL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      			! Function RSL at order L=2 and S=IS versus mu
      DOUBLE PRECISION XTL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      			! Function TSL at order L=2 and S=IS versus mu

	! Element Pij of the aerosol phase matrix for an order
	! IS of the Fourier series decomposition.
	! With mup the incidence direction (-1 < mup < 1)
	! and mu the scattering one (-1 < mu < 1),
	! the scattering is expressed by Pij(mu,mup).
      DOUBLE PRECISION 
     &    BP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &       -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  ! Element P11.
      DOUBLE PRECISION 
     &    GR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &       -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  ! P12 or P21.
      DOUBLE PRECISION 
     &    GT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &       -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  ! -P13 or -P31.
    

      DOUBLE PRECISION CH(0:CTE_OS_NT)  ! Attenuation of the direct solar
                                        ! beam / 4.

	! Source function vector versus the level of the profile
	! and the direction.
      DOUBLE PRECISION I2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      						!  Source function for I.
      DOUBLE PRECISION Q2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      						!  Source function for Q.
      DOUBLE PRECISION U2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      						!  Source function for U.

      DOUBLE PRECISION SPL   ! Value of the function XPL for mu = mus.
      DOUBLE PRECISION SA1   ! Element P11 of the molecular phase matrix.
      DOUBLE PRECISION SB1   ! Element P21 of the molecular phase matrix.
      DOUBLE PRECISION SC1   ! Element P31 of the molecular phase matrix.
      DOUBLE PRECISION SA2   ! Element P11 of the aerosol phase matrix.
      DOUBLE PRECISION SB2   ! Element P21 of the aerosol phase matrix.
      DOUBLE PRECISION SC2   ! Element P31 of the aerosol phase matrix.
 
      DOUBLE PRECISION PCRAY	! Value of the percentage of molecules.
      DOUBLE PRECISION PCAER	! Value of the percentage of aerosols.

      DOUBLE PRECISION ATTDIR	! Single value of the array CH.

      INTEGER*4 NT      ! Number of layers of the atmospheric profile
          
      INTEGER*4 IS	! Current order of the Fourier series expansion.
      INTEGER*4 NBMU ! Number of actual angles of the simulations
      
      INTEGER JK	! Storage number of mus in the array RMU
      			! (0 by default).

      INTEGER*4 J,K


C-------------------------------------------------------------      
C* Computation of the elements P11 = SAi , P21 = SBi  and P31 = SCi
C* of the phase matrix, according to the order IS of the Fourier
C* series expansion (i = 1 : molecules, i = 2 : aerosols).
C* The incidence direction is the solar illumination
C* direction (mus <=> JK = 0).
C-------------------------------------------------------------      


C* Loop on the scattering directions.
C------------------------------------------
      DO 100 J=-NBMU,NBMU

        IF ((IS-2).GT.0) GOTO 201
      
C* Computations for IS <= 2.
C-----------------------
 200    SPL=XPL(JK)
        SA1=BETA0+BETA2*XPL(J)*SPL
        SA2=BP(JK,J)
        SB1=GAMMA2*XRL(J)*SPL
        SB2=GR(JK,J)
        SC1=GAMMA2*XTL(J)*SPL
        SC2=GT(JK,J)

        GOTO 202

C* Computations for IS > 2. The molecular phase
C* matrix is null for IS > 2
C------------------------------------------------
 201    SA2=BP(JK,J)
        SA1=0.
        SB2=GR(JK,J)
        SB1=0.
        SC2=GT(JK,J)
        SC1=0.

202     CONTINUE


C------------------------------------------------------------------------
C* Computation of the primary scattering source function for each level
C* of the profile and each direction.
C------------------------------------------------------------------------
        DO 101 K=0,NT
           ATTDIR=CH(K)
           PCRAY=YDEL(K)
           PCAER=XDEL(K)
           I2(K,J)=ATTDIR*(SA2*PCAER+SA1*PCRAY)
           Q2(K,J)=ATTDIR*(SB2*PCAER+SB1*PCRAY)
           U2(K,J)=-ATTDIR*(SC2*PCAER+SC1*PCRAY)
  101    CONTINUE
  
  100 CONTINUE
  

      END  	!END OF PROCEDURE SOS_FSOURCE_ORDRE1







C*============================================================================
C PROCEDURE: SOS_FSOURCE_ORDREIG
C ==========
C      This procedure computes the source function vector of the scattering
C      of order IG>1, from the radiation field of order IG-1 and the
C      scattering properties of the atmosphere (phase function).
C
C      It computes the elements of the phase matrix,
C      according to the order IS of the Fourier series expansion.
C      The following matrix P(mu,mup) is built indirectly:
C         BP(mu,mup)    GR(mu,mup)    -GT(mu,mup)
C         GR(mup,mu)    ARR(mu,mup)   -ART(mup,mu)
C        -GT(mup,mu)   -ART(mu,mup)    ATT(mu,mup)
C 
C      The terms actually computed are equivalent following the rules
C      of index permutation.
C	
C      A spatial integration is applied on the incidence directions.
C
C      Computation for each level of the profile and for each propagation
C      direction.
C
C
C
C Description of parameters
C -------------------------- 
C   
C       IS   (I4)    (I)  Current order of the Fourier series expansion.
C       NBMU (I4)    (I)  Number of actually useful (positive) angles
C       XDEL(0:CTE_OS_NT)  (double)  (I)  Percentage of aerosols per level.
C       YDEL(0:CTE_OS_NT)  (double)  (I)  Percentage of molecules per level.
C
C       BETA0	  (double)  (I)
C       BETA2   (double)  (I) 	.
C       GAMMA2  (double)  (I)
C       ALPHA2  (double)  (I)
C                   Coefficients of the expansion of the molecular
C                   phase matrix in Legendre functions
C
C       XPL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C      		     Legendre function PSL at order L=2 and S=IS
C      		     versus the angle mu.
C       XRL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C     		     Function RSL at order L=2 and S=IS versus mu.
C       XTL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C     		     Function TSL at order L=2 and S=IS versus mu.
C
C      	
C       I1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)    Field I.
C       Q1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)    Field Q.
C       U1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)    Field U.
C         	     Radiance field of order IG versus the
C                  level of the profile and the direction.
C	
C       BP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C                                                        Element P11.
C       GR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C                                                        Element P12 or P21.
C       GT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C                                                        Element -P13 or -P31.
C       ARR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (double)  (I)
C                                                        Element P22.
C       ART(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (double)  (I)
C                                                        Element -P23 or -P32.
C       ATT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (double)  (I)
C                                                        Element P33.
C             Element Pij of the aerosol phase matrix for an order
C             IS of the Fourier series decomposition.
C             With mup the incidence direction ( -1 < mup < 1)
C             and mu the scattering one ( -1 < mu < 1),
C             the scattering is expressed by Pij(mu,mup).
C
C       GA(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I) Array of the Gauss weights.
C
C       I2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (O)
C       Q2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (O)
C       U2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (O)
C	         Source function vector I,Q,U versus the level of the profile
C	         and the direction.
C
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C None
C
C==============================================================================
      SUBROUTINE SOS_FSOURCE_ORDREIG(IS,NBMU,NT,XDEL,YDEL,
     &                               BETA0,BETA2,GAMMA2,ALPHA2,
     &                               XPL,XRL,XTL,I1,Q1,U1,
     &                               BP,GR,GT,ARR,ART,ATT,GA,
     &                               I2,Q2,U2)
      
      IMPLICIT NONE
      

C* Definition of variables
C*------------------------------------------------------------------

      DOUBLE PRECISION XDEL(0:CTE_OS_NT)  ! Percentage of aerosols per level.
      DOUBLE PRECISION YDEL(0:CTE_OS_NT)  ! Percentage of molecules per level.

      DOUBLE PRECISION BETA0	! Coefficients of the expansion of the molecular
      DOUBLE PRECISION BETA2	! phase matrix in Legendre functions.
      DOUBLE PRECISION GAMMA2	
      DOUBLE PRECISION ALPHA2	

      DOUBLE PRECISION XPL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      			      ! Legendre function PSL at order L=2 and S=IS
      			      ! versus the angle mu.
      DOUBLE PRECISION XRL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      			! Function RSL at order L=2 and S=IS versus mu.
      DOUBLE PRECISION XTL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      			! Function TSL at order L=2 and S=IS versus mu.


      	! Radiance field of order IG-1 versus the level of the profile
	! and the direction.
      DOUBLE PRECISION I1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	!  Field I.
      DOUBLE PRECISION Q1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	!  Field Q.
      DOUBLE PRECISION U1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	!  Field U.

	! Source function vector versus the level of the profile
	! and the direction.
      DOUBLE PRECISION I2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      						!  Source function for I.
      DOUBLE PRECISION Q2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      						!  Source function for Q.
      DOUBLE PRECISION U2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      						!  Source function for U.


	! Element Pij of the aerosol phase matrix for an order
	! IS of the Fourier series decomposition.
	! With mup the incidence direction (-1 < mup < 1)
	! and mu the scattering one (-1 < mu < 1),
	! the scattering is expressed by Pij(mu,mup).
      DOUBLE PRECISION 
     &    BP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &       -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  ! Element P11.
      DOUBLE PRECISION 
     &    GR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &       -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  ! P12 or P21.
      DOUBLE PRECISION 
     &    GT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &       -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  ! -P13 or -P31.
      DOUBLE PRECISION 
     &    ARR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &        -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! P22.
      DOUBLE PRECISION 
     &    ART(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &        -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! -P23 or -P32.
      DOUBLE PRECISION 
     &    ATT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &        -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! P33.


      DOUBLE PRECISION GA(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! Array of the Gauss weights

      DOUBLE PRECISION XPK	! Value of XPL for the scattering
      				! directions mu > 0 (index K).
      DOUBLE PRECISION XRK	! Value of XRL for mu > 0.
      DOUBLE PRECISION XTK	! Value of XTL for mu > 0.
      DOUBLE PRECISION YPK	! Value of XPL for -mu.
      DOUBLE PRECISION YRK	! Value of XRL for -mu.
      DOUBLE PRECISION YTK	! Value of XTL for -mu.
      	
      DOUBLE PRECISION XPJ	! Value of XPL for the incident directions
      				! mu = mup > 0 (index J).
      DOUBLE PRECISION XRJ	! Value of XRL for mu = mup > 0.
      DOUBLE PRECISION XTJ	! Value of XTL for mu = mup > 0.
      DOUBLE PRECISION YPJ	! Value of XPL for mu = -mup.
      DOUBLE PRECISION YRJ	! Value of XRL for mu = -mup.
      DOUBLE PRECISION YTJ	! Value of XTL for mu = -mup.

      DOUBLE PRECISION XI1	! Value of I of the upward field of order IG-1.
      DOUBLE PRECISION XQ1	! Value of Q of the upward field of order IG-1.
      DOUBLE PRECISION XU1	! Value of U of the upward field of order IG-1.
      DOUBLE PRECISION XI2	! Value of I of the downward field of order IG-1.
      DOUBLE PRECISION XQ2	! Value of Q of the downward field of order IG-1.
      DOUBLE PRECISION XU2	! Value of U of the downward field of order IG-1.

      DOUBLE PRECISION II1	! Intermediate value for the computation of I2
      				! (mu < 0).
      DOUBLE PRECISION QQ1	! Intermediate value for the computation of Q2
      				! (mu < 0).
      DOUBLE PRECISION UU1	! Intermediate value for the computation of U2
      				! (mu < 0).
      DOUBLE PRECISION II2 	! Intermediate value for the computation of I2
      				! (mu > 0).
      DOUBLE PRECISION QQ2 	! Intermediate value for the computation of Q2
      				! (mu > 0).
      DOUBLE PRECISION UU2 	! Intermediate value for the computation of U2
      				! (mu > 0).

      DOUBLE PRECISION Z   ! Single value of the array of Gauss weights.

      DOUBLE PRECISION PCRAY	! Value of the percentage of molecules.
      DOUBLE PRECISION PCAER	! Value of the percentage of aerosols.
 
      DOUBLE PRECISION BPJK	    ! Value of the elements of the phase
      DOUBLE PRECISION BPJMK    ! matrix of the molecules + aerosols mixture.
      DOUBLE PRECISION GTJMK
      DOUBLE PRECISION GTJK
      DOUBLE PRECISION GTKJ
      DOUBLE PRECISION GTKMJ
      DOUBLE PRECISION GRJK
      DOUBLE PRECISION GRKJ
      DOUBLE PRECISION GRKMJ
      DOUBLE PRECISION GRJMK
      DOUBLE PRECISION ARRJK
      DOUBLE PRECISION ARRJMK
      DOUBLE PRECISION ARTJK
      DOUBLE PRECISION ARTJMK
      DOUBLE PRECISION ARTKJ
      DOUBLE PRECISION ARTKMJ
      DOUBLE PRECISION ATTJK
      DOUBLE PRECISION ATTJMK

      INTEGER*4 NT	! Number of layers of the atmospheric profile
      INTEGER*4 IS	! Current order of the Fourier series expansion.
      INTEGER*4 NBMU    ! Number of actual angles of the simulations

      INTEGER*4 I
      INTEGER*4 K
      INTEGER*4 J
      


C------------------------------------------------------
C* Computation of the elements of the phase matrix, according to
C* the order IS of the Fourier series expansion.
C------------------------------------------------------

      IF ((IS-2).GT.0) GOTO 211
	
      
C--------------------------------
C* Computations for IS <= 2.
C-------------------------------

C* Loop on the scattering directions mu (mu > 0).
C------------------------------------------------------
 210  DO 455 K=1,NBMU

C* Value of the functions PSL, RSL and TSL at order L=2 and S=IS
C* versus the angle mup and -mup
C* (for the computation of the molecular phase matrix
C*  by the calling program).
C------------------------------------------------------
         XPK=XPL(K)
         XRK=XRL(K)
         XTK=XTL(K)
         YPK=XPL(-K)
         YRK=XRL(-K)
         YTK=XTL(-K)

C* Loop on the profile.
C-----------------------
         DO 466 I=0,NT
            II1=0.
            II2=0.
            QQ1=0.
            QQ2=0.
            UU1=0.
            UU2=0.
            PCAER=XDEL(I)
            PCRAY=YDEL(I)

C* Loop on the incidence directions mup.
C------------------------------------------------------
            DO 477 J=1,NBMU

C* Computation of the elements of the phase matrix P(mu,mup)
C* at order IS for the molecules + aerosols mixture.
C------------------------------------------------------
               XPJ=XPL(J)
               XRJ=XRL(J)
               XTJ=XTL(J)
               YPJ=XPL(-J)
               YRJ=XRL(-J)
               YTJ=XTL(-J)

               BPJK=BP(J,K)*PCAER+PCRAY*(BETA0+BETA2*XPJ*XPK)
               BPJMK=BP(J,-K)*PCAER+PCRAY*(BETA0+BETA2*XPJ*YPK)
               GTJMK=GT(J,-K)*PCAER+PCRAY*(GAMMA2*XPJ*YTK)
               GTJK=GT(J,K)*PCAER+PCRAY*(GAMMA2*XPJ*XTK)
               GTKMJ=GT(K,-J)*PCAER+PCRAY*(GAMMA2*XPK*YTJ)
               GTKJ=GT(K,J)*PCAER+PCRAY*(GAMMA2*XPK*XTJ)
               GRJK=GR(J,K)*PCAER+PCRAY*(GAMMA2*XPJ*XRK)
               GRJMK=GR(J,-K)*PCAER+PCRAY*(GAMMA2*XPJ*YRK)
               GRKJ=GR(K,J)*PCAER+PCRAY*(GAMMA2*XPK*XRJ)
               GRKMJ=GR(K,-J)*PCAER+PCRAY*(GAMMA2*XPK*YRJ)
               ARRJK=ARR(J,K)*PCAER+PCRAY*(ALPHA2*XRJ*XRK)
               ARRJMK=ARR(J,-K)*PCAER+PCRAY*(ALPHA2*XRJ*YRK)
               ARTJK=ART(J,K)*PCAER+PCRAY*(ALPHA2*XTJ*XRK)
               ARTJMK=ART(J,-K)*PCAER+PCRAY*(ALPHA2*XTJ*YRK)
               ARTKJ=ART(K,J)*PCAER+PCRAY*(ALPHA2*XTK*XRJ)
               ARTKMJ=ART(K,-J)*PCAER+PCRAY*(ALPHA2*XTK*YRJ)
               ATTJMK=ATT(J,-K)*PCAER+PCRAY*(ALPHA2*XTJ*YTK)
               ATTJK=ATT(J,K)*PCAER+PCRAY*(ALPHA2*XTJ*XTK)

C* Computation of the source function of order IG, for each
C* level of the profile and each direction.
C------------------------------------------------------
               Z=GA(J)

C* Retrieval of the radiance values of order IG-1.
C------------------------------------------------------
               XI1=I1(I,J)
               XI2=I1(I,-J)
               XQ1=Q1(I,J)
               XQ2=Q1(I,-J)
               XU1=U1(I,J)
               XU2=U1(I,-J)

C* Spatial integration over mup.
C---------------------------------
               II2=II2+Z*(XI1*BPJK+XI2*BPJMK+XQ1*GRKJ
     &             +XQ2*GRKMJ-XU1*GTKJ-XU2*GTKMJ)
               II1=II1+Z*(XI1*BPJMK+XI2*BPJK+XQ1*GRKMJ
     &             +XQ2*GRKJ+XU1*GTKMJ+XU2*GTKJ)
               QQ2=QQ2+Z*(XI1*GRJK+XI2*GRJMK+XQ1*ARRJK
     &             +XQ2*ARRJMK+XU2*ARTJMK-XU1*ARTJK)
               QQ1=QQ1+Z*(XI1*GRJMK+XI2*GRJK+XQ1*ARRJMK
     &             +XQ2*ARRJK-XU1*ARTJMK+XU2*ARTJK)
               UU2=UU2-Z*(XI1*GTJK-XI2*GTJMK+XQ1*ARTKJ
     &             +XQ2*ARTKMJ-XU1*ATTJK-XU2*ATTJMK)
               UU1=UU1-Z*(XI1*GTJMK-XI2*GTJK-XQ1*ARTKMJ
     &             -XQ2*ARTKJ-XU1*ATTJMK-XU2*ATTJK)
 477        CONTINUE

C* Storage of the results
C--------------------------
            I2(I,K)=II2 * 0.5
            I2(I,-K)=II1 * 0.5
            Q2(I,K)=QQ2 * 0.5
            Q2(I,-K)=QQ1 * 0.5
            U2(I,K)=UU2 * 0.5
            U2(I,-K)=UU1 * 0.5
	    
	    
 466     CONTINUE   ! Loop on I
	    
 455  CONTINUE   ! Loop on K
 
      GOTO 212
      
      
C-----------------------------
C* Computations for IS > 2.
C-----------------------------

C* Loop on the scattering directions mu (mu > 0).
C------------------------------------------------------
 211  DO 45 K=1,NBMU

C* Loop on the profile.
C------------------------
         DO 46 I=0,NT
            II1=0.
            II2=0.
            QQ1=0.
            QQ2=0.
            UU1=0.
            UU2=0.
            PCAER=XDEL(I)

C* Loop on the incidence directions mup
C-------------------------------------------
            DO 47 J=1,NBMU

C* Computation of the elements of the aerosol phase matrix P(mu,mup)
C* at order IS, weighted by the aerosol proportion.
C------------------------------------------------------
               BPJK=BP(J,K)*PCAER
               BPJMK=BP(J,-K)*PCAER
               GTJMK=GT(J,-K)*PCAER
               GTJK=GT(J,K)*PCAER
               GTKMJ=GT(K,-J)*PCAER
               GTKJ=GT(K,J)*PCAER
               GRJK=GR(J,K)*PCAER
               GRJMK=GR(J,-K)*PCAER
               GRKJ=GR(K,J)*PCAER
               GRKMJ=GR(K,-J)*PCAER
               ARRJK=ARR(J,K)*PCAER
               ARRJMK=ARR(J,-K)*PCAER
               ARTJK=ART(J,K)*PCAER
               ARTJMK=ART(J,-K)*PCAER
               ARTKJ=ART(K,J)*PCAER
               ARTKMJ=ART(K,-J)*PCAER
               ATTJMK=ATT(J,-K)*PCAER
               ATTJK=ATT(J,K)*PCAER
	       
C* Computation of the source function of order IG, for each
C* level of the profile and each direction.
C------------------------------------------------------
               Z=GA(J)

C* Retrieval of the radiance values of order IG-1.
C------------------------------------------------------
               XI1=I1(I,J)
               XI2=I1(I,-J)
               XQ1=Q1(I,J)
               XQ2=Q1(I,-J)
               XU1=U1(I,J)
               XU2=U1(I,-J)

C* Spatial integration over mup.
C-------------------------------
               II2=II2+Z*(XI1*BPJK+XI2*BPJMK+XQ1*GRKJ
     &             +XQ2*GRKMJ-XU1*GTKJ-XU2*GTKMJ)
               II1=II1+Z*(XI1*BPJMK+XI2*BPJK+XQ1*GRKMJ
     &             +XQ2*GRKJ+XU1*GTKMJ+XU2*GTKJ)
               QQ2=QQ2+Z*(XI1*GRJK+XI2*GRJMK+XQ1*ARRJK
     &             +XQ2*ARRJMK+XU2*ARTJMK-XU1*ARTJK)
               QQ1=QQ1+Z*(XI1*GRJMK+XI2*GRJK+XQ1*ARRJMK
     &             +XQ2*ARRJK-XU1*ARTJMK+XU2*ARTJK)
               UU2=UU2-Z*(XI1*GTJK-XI2*GTJMK+XQ1*ARTKJ
     &             +XQ2*ARTKMJ-XU1*ATTJK-XU2*ATTJMK)
               UU1=UU1-Z*(XI1*GTJMK-XI2*GTJK-XQ1*ARTKMJ
     &             -XQ2*ARTKJ-XU1*ATTJMK-XU2*ATTJK)
     
   47       CONTINUE


C* Storage of the results.
C--------------------------
            I2(I,K)=II2 * 0.5
            I2(I,-K)=II1 * 0.5
            Q2(I,K)=QQ2  * 0.5
            Q2(I,-K)=QQ1 * 0.5
            U2(I,K)=UU2  * 0.5
            U2(I,-K)=UU1 * 0.5

   46    CONTINUE ! Loop on I
   45 CONTINUE    ! Loop on K

  212 CONTINUE

      RETURN
      END  	!END OF PROCEDURE SOS_FSOURCE_ORDREIG







C*===========================================================================
C PROCEDURE: SOS_FSOURCE_DIFF_FRESNEL1
C ==========
C      This procedure computes the source function of the scattering
C      of the solar radiation going up after a direct reflection
C      on a flat sea, according to the Fresnel laws (scattering of the
C      reflection field of order 1).
C
C
C Description of parameters
C -------------------------- 
C
C      F11sun  (double)  (I) Element F11 of the Fresnel reflection
C				 matrix for the solar incidence.
C      F12sun  (double)  (I) Same for the element F12.
C
C       XDEL(0:CTE_OS_NT)  (double)  (I)  Percentage of aerosols per level.
C       YDEL(0:CTE_OS_NT)  (double)  (I)  Percentage of molecules per level.
C
C       BETA0	  (double)  (I)
C      BETA2   (double)  (I)
C       GAMMA2  (double)  (I)
C       ALPHA2  (double)  (I)
C             Coefficients of the expansion of the molecular
C             phase matrix in Legendre functions.
C
C       BP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C                                                        Element P11.
C       GR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C                                                        Element P12 or P21.
C       GT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C                                                        Element -P13 or -P31.
C       ARR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (double)  (I)
C                                                        Element P22.
C       ART(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (double)  (I)
C                                                        Element -P23 or -P32.
C             Element Pij of the aerosol phase matrix for an order
C             IS of the Fourier series decomposition.
C             With mup the incidence direction ( -1 < mup < 1)
C             and mu the scattering one ( -1 < mu < 1),
C             the scattering is expressed by Pij(mu,mup).
C
C       XPL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C      		        Legendre function PSL at order L=2
C   		        and S=IS versus the angle mu.
C       XRL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C      		        Function RSL at order L=2 and S=IS.
C       XTL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C      		        Function TSL at order L=2 and S=IS
C	             Legendre functions and derived functions for the computation
C	             of the phase matrix.
C
C     IS  (I4)  (I)   Order of the Fourier series decomposition.
C
C       NBMU (I4)    (I)  Number of actually useful (positive) angles
C
C      MUS (double)  (I)   Cosine of the solar zenith angle.
C
C     H(0:CTE_OS_NT) (double) (I)  Optical thickness at each level.
C
C       I2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (O)
C                                                  Source function for I
C       Q2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (O)
C                                                  Source function for Q
C       U2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (O)
C                                                  Source function for U
C             Source function vector, versus the level of the
C             profile and the direction, for the scattering of
C             the light obtained by Fresnel reflection of the direct
C             solar beam.
C                          
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C None
C
C==============================================================================
      SUBROUTINE SOS_FSOURCE_DIFF_FRESNEL1(F11sun,F12sun,XDEL,YDEL,NT,
     &                                     BETA0,BETA2,GAMMA2,ALPHA2,
     &                                     BP,GR,GT,ARR,ART,
     &                                     XPL,XRL,XTL,IS,NBMU,
     &					        MUS,H,I2,Q2,U2)
     
      IMPLICIT NONE


C* Definition of variables
C*-----------------------------------------------------------------------     
	! Source function vector versus the level of the profile
	! and the direction.
      DOUBLE PRECISION I2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      						!  Source function for I.
      DOUBLE PRECISION Q2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      						!  Source function for Q.
      DOUBLE PRECISION U2(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      						!  Source function for U.
						
      DOUBLE PRECISION F11sun   ! Element of the Fresnel matrix
      			           ! with as incidence the direction
			           ! of solar illumination.
      DOUBLE PRECISION F12sun   ! Element of the Fresnel matrix.

      DOUBLE PRECISION XDEL(0:CTE_OS_NT)  !Percentage of aerosols per level.
      DOUBLE PRECISION YDEL(0:CTE_OS_NT)  !Percentage of molecules per level.

      DOUBLE PRECISION MUS               ! Cosine of the solar zenith angle.
      DOUBLE PRECISION H(0:CTE_OS_NT)    !  Optical thickness for each level.
           
      DOUBLE PRECISION YR   ! Percentage of molecules at level K of the profile.
      DOUBLE PRECISION YYR  ! Percentage of molecules at level K+1 of the profile.
      DOUBLE PRECISION XP   ! Percentage of aerosols at level K of the profile.
      DOUBLE PRECISION XXP  ! Percentage of aerosols at level K+1 of the profile.


      DOUBLE PRECISION XPL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      			! Legendre function PSL at order L=2 and S=IS
      			! versus the angle mu.
      DOUBLE PRECISION XRL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      			! Function RSL at order L=2 and S=IS versus mu.
      DOUBLE PRECISION XTL(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      			! Function TSL at order L=2 and S=IS versus mu.

      DOUBLE PRECISION SPL   ! Value of the function XPL for mu = mus
           	
      DOUBLE PRECISION BETA0 ! Coefficients of the expansion of the molecular
      DOUBLE PRECISION BETA2 ! phase matrix in Legendre functions
      DOUBLE PRECISION GAMMA2	
      DOUBLE PRECISION ALPHA2	

	! Element Pij of the aerosol phase matrix for an order
	! IS of the Fourier series decomposition.
	! With mup the incidence direction (-1 < mup < 1)
	! and mu the scattering one (-1 < mu < 1),
	! the scattering is expressed by Pij(mu,mup).
      DOUBLE PRECISION 
     &    BP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &       -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  ! Element P11.
      DOUBLE PRECISION 
     &    GR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &       -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  ! P12 or P21.
      DOUBLE PRECISION 
     &    GT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &       -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  ! -P13 or -P31.
      DOUBLE PRECISION 
     &    ARR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &        -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! P22.
      DOUBLE PRECISION 
     &    ART(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX,
     &        -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! -P23 or -P32.
     
           
      DOUBLE PRECISION BP0MJ
      DOUBLE PRECISION BP0J
      DOUBLE PRECISION GRJ0
      DOUBLE PRECISION GR0J
      DOUBLE PRECISION GR0MJ
      DOUBLE PRECISION GRMJ0
      DOUBLE PRECISION GT0MJ
      DOUBLE PRECISION GT0J
      DOUBLE PRECISION ARR0MJ
      DOUBLE PRECISION ARR0J
      DOUBLE PRECISION ARTJ0
      DOUBLE PRECISION ARTMJ0
       
      
      DOUBLE PRECISION COEFNT
      DOUBLE PRECISION COEFK
      DOUBLE PRECISION COEFKP1 
      
      INTEGER*4 NT   ! Number of layers of the atmospheric profile
      INTEGER IS	! Order of the Fourier series decomposition.
      INTEGER K	! Index of the profile level.
      INTEGER J	! Index of the incidence angle.
      INTEGER*4 NBMU ! Number of actual angles of the simulations
      
       
C* Processing
C*------------

C* Initialisation

      DO J=-NBMU,NBMU
         DO K=0,NT
	      I2(K,J)=0.
	      Q2(K,J)=0.
	      U2(K,J)=0.
         ENDDO    
      ENDDO 


      COEFNT  = DEXP(2.*H(NT)/MUS) / 4.  

      SPL=XPL(0)

      
      DO 1500 K= 0 , NT-1
      
         YR=YDEL(K)
         XP=XDEL(K)
         
         YYR=YDEL(K+1)
         XXP=XDEL(K+1)
	 
         DO 1502 J=1,NBMU

C* Phase matrix elements for the molecules + aerosols mixture
            IF (IS.LE.2) THEN	 
	    
                BP0MJ=BP(0,-J)*XP+(BETA0+BETA2*XPL(-J)*SPL)*YR
                BP0J=BP(0,J)*XXP+(BETA0+BETA2*XPL(J)*SPL)*YYR
            
                GRJ0=GR(J,0)*XXP+YYR*XRL(0)*XPL(J)*GAMMA2
	         GR0J=GR(0,J)*XXP+YYR*XRL(J)*XPL(0)*GAMMA2
                GR0MJ=GR(0,-J)*XP+YR*XRL(-J)*SPL*GAMMA2
	         GRMJ0 =GR(-J,0)*XP+YR*GAMMA2*XRL(0)*XPL(-J)	
	    
                GT0MJ=GT(0,-J)*XP+YR*GAMMA2*SPL*XTL(-J)
                GT0J=GT(0,J)*XXP+YYR*GAMMA2*SPL*XTL(J)
	                
                ARR0MJ=ARR(0,-J)*XP+ALPHA2*YR*XRL(0)*XRL(-J)
                ARR0J=ARR(0,J)*XXP+ALPHA2*YYR*XRL(0)*XRL(J)            
            
                ARTJ0=ART(J,0)*XXP+YYR*ALPHA2*XTL(J)*XRL(0)
                ARTMJ0=ART(-J,0)*XP+YR*ALPHA2*XTL(-J)*XRL(0)	
		
            ELSE
	    
                BP0MJ=BP(0,-J)*XP
                BP0J=BP(0,J)*XXP
            
                GRJ0=GR(J,0)*XXP
	         GR0J=GR(0,J)*XXP
                GR0MJ=GR(0,-J)*XP
	         GRMJ0 =GR(-J,0)*XP
	    
                GT0MJ=GT(0,-J)*XP
                GT0J=GT(0,J)*XXP
	                
                ARR0MJ=ARR(0,-J)*XP
                ARR0J=ARR(0,J)*XXP          
            
                ARTJ0=ART(J,0)*XXP
                ARTMJ0=ART(-J,0)*XP	  
		  
	    ENDIF
   
	     
	     
C* Source function for the upward field
            COEFK   = COEFNT * DEXP(-H(K)/MUS)
	    
            I2(K,J)= COEFK * ( F11sun*BP0MJ + F12sun*GRMJ0 )		
            Q2(K,J)= COEFK * (  F11sun*GR0MJ + F12sun*ARR0MJ )
            U2(K,J)= COEFK * (  F11sun*GT0MJ + F12sun*ARTMJ0 )		
       
C* Source function for the downward field
            COEFKP1 = COEFNT * DEXP(-H(K+1)/MUS) 
	    
            I2(K+1,-J)=  COEFKP1 * ( F11sun*BP0J + F12sun*GRJ0 )
            Q2(K+1,-J)=  COEFKP1 * ( F11sun*GR0J + F12sun*ARR0J )
            U2(K+1,-J)=  COEFKP1 * ( F11sun*GT0J + F12sun*ARTJ0 )
	     
 1502    CONTINUE
 1500 CONTINUE
      

      END      !END OF PROCEDURE SOS_FSOURCE_DIFF_FRESNEL1
       
       
       
       
       
       

C*===========================================================================
C PROCEDURE: SOS_PARAM_CONV
C ==========
C      This procedure computes the convergence parameter for the test
C      of convergence in geometric series.
C
C      Principle:
C      --------
C      The convergence in geometric series of the successive scatterings
C      is analysed from the orders IG-2, IG-1 and IG.
C      The test is made for the radiance field emerging at TOA and for
C      the field illuminating the ground. If there is convergence at these limits, it is assumed
C      that there is convergence in the whole atmosphere.
C
C      Assume a behaviour of convergent geometric series from
C      the scattering order IG-2.
C      The ratio of the series is estimated between the terms IG-2 and IG-1, as well as
C      between the terms IG-1 and IG:
C 		_ q  = D1 / A1  (ratio of the orders IG-1 / IG-2 for I)
C		_ q' = G1 / D1  (ratio of the orders IG / IG-1 for I)
C     It is tested that q and q' are close (condition for a geometric series
C     approximation).
C     At the same time, the ratio of the sum of the terms of the convergent geometric
C     series from the order IG (order IG / (1-q) )
C     to the sum of the previous terms is tested
C     (from 1 to IG-1, stored in I3, Q3 or U3).
C
C     The routine outputs the convergence parameter Z1, which is tested by the
C     calling program.
C
C Description of parameters
C -------------------------- 
C
C       NBMU (I4)    (I)  Number of actually useful (positive) angles
C
C    I3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)		
C                         Stokes parameter I for the sum of the
C      		          previous scattering orders up to the current IG.
C    Q3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)  Stokes parameter Q.
C    U3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)  Stokes parameter U.
C            			
C    A1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C       	          Stokes parameter I for the order IG-2.
C    B1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C      			  Stokes parameter Q for the order IG-2.
C    C1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C      			  Stokes parameter U for the order IG-2.
C
C    D1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C      			  Stokes parameter I for the order IG-1.
C    E1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C      			  Stokes parameter Q for the order IG-1.
C    F1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C      			  Stokes parameter U for the order IG-1.
C
C    G1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C      		          Stokes parameter I for the order IG.
C    H1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C      			  Stokes parameter Q for the order IG.
C    P1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C      			  Stokes parameter U for the order IG.
C
C    Z1	   (double)   (O)    Convergence parameter.
C
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C None
C
C==============================================================================
      SUBROUTINE SOS_PARAM_CONV(NBMU,A1,B1,C1,D1,E1,F1,G1,H1,P1,
     &                          I3,Q3,U3,Z1)
      
      IMPLICIT NONE
      

C* Definition of variables
C*------------------------------------------------------------------

      DOUBLE PRECISION Z1	! Convergence parameter.

	! Functions for the convergence test in geometric series,
	! field emerging at TOA (mu > 0) and downward at the ground (mu < 0).
	! The functions are established for an order S of the expansion
	! in Fourier series, versus the angle mu.
      DOUBLE PRECISION I3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)		
                            ! Stokes parameter I for the sum of the
      				! previous scattering orders up to the current IG.
      DOUBLE PRECISION Q3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
      				! Stokes parameter Q.
      DOUBLE PRECISION U3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
      				! Stokes parameter U.
            			
      DOUBLE PRECISION A1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
       			! Stokes parameter I for the order IG-2.
      DOUBLE PRECISION B1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Stokes parameter Q for the order IG-2.
      DOUBLE PRECISION C1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Stokes parameter U for the order IG-2.

      DOUBLE PRECISION D1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Stokes parameter I for the order IG-1.
      DOUBLE PRECISION E1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Stokes parameter Q for the order IG-1.
      DOUBLE PRECISION F1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      			  	! Stokes parameter U for the order IG-1.

      DOUBLE PRECISION G1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Stokes parameter I for the order IG.
      DOUBLE PRECISION H1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Stokes parameter Q for the order IG.
      DOUBLE PRECISION P1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      			  	! Stokes parameter U for the order IG.

      DOUBLE PRECISION Y	! Variable for the stop tests.

      INTEGER*4 K
      INTEGER*4 NBMU    ! Number of actual angles of the simulations

C* Determination of the maximum value of Z1 for all the angles mu.
C----------------------------------------------------------------------------
      Z1=0.D+00

      DO 99 K=-NBMU,NBMU

         IF (K.NE.0) THEN 
 
            IF ((A1(K).NE.0.D+00).AND.(D1(K).NE.0.D+00)
     &          .AND.(I3(K).NE.0.D+00)) THEN
                Y=((G1(K)/D1(K)-D1(K)/A1(K))
     &            /((1-G1(K)/D1(K))**2)*(G1(K)/I3(K)))
                 Z1=DMAX1(Z1,DABS(Y))
            ENDIF

	    IF ((B1(K).NE.0.D+00).AND.(E1(K).NE.0.D+00)
     &          .AND.(Q3(K).NE.0.D+00)) THEN
                Y=((H1(K)/E1(K)-E1(K)/B1(K))
     &            /((1-H1(K)/E1(K))**2)*(H1(K)/Q3(K)))
                Z1=DMAX1(Z1,DABS(Y))
            ENDIF

	    IF ((C1(K).NE.0.D+00).AND.(F1(K).NE.0.D+00)
     &          .AND.(U3(K).NE.0.D+00)) THEN
                Y=((P1(K)/F1(K)-F1(K)/C1(K))
     &            /((1-P1(K)/F1(K))**2)*(P1(K)/U3(K)))
                Z1=DMAX1(Z1,DABS(Y))
            ENDIF
 
         ENDIF !End of test K!=0
 
   
   99 CONTINUE

      RETURN
      END  	!END OF PROCEDURE SOS_PARAM_CONV




C*============================================================================
C PROCEDURE: SOS_ARRET_DIFFUS_1
C ==========
C      This procedure prepares the stop test on the scatterings when the
c      terms I, Q, U are very small
C
C
C Description of parameters
C -------------------------- 
C    
C       NBMU (I4)    (I)  Number of actually useful (positive) angles
C
C    I1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)  (I)
C    Q1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)  (I)
C    U1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)  (I)
C      	 Radiance field of order IG versus the level of the
C                          profile and the direction.
C
C
C    Z1   (double)   (O)    Threshold variable for the stop tests.
C
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C None
C
C==============================================================================
      SUBROUTINE SOS_ARRET_DIFFUS_1(NBMU,NT,I1,Q1,U1,Z1)
      
      IMPLICIT NONE
      

C* Definition of variables
C*------------------------------------------------------------------

      DOUBLE PRECISION Z1	! Threshold variable for the stop tests.

      	! Radiance field of order IG versus the level of the profile
	! and the direction.
      DOUBLE PRECISION I1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	!  Field I.
      DOUBLE PRECISION Q1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	!  Field Q.
      DOUBLE PRECISION U1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	!  Field U.
 
      INTEGER*4 NT	! Number of layers of the atmospheric profile
      INTEGER*4 IND	! Layer number for mu > 0 or mu < 0
      INTEGER*4 K 
      INTEGER*4 NBMU ! Number of actual angles of the simulations

C------------------------------------------------------ 
C* The value of the term of order IG is compared with the sum of the
C* previous terms.
C* Determination of the maximum value of their ratio
C* for the upward field at TOA and the downward field at the ground.
C------------------------------------------------------
      Z1=0.D+00
      DO 88 K=-NBMU,NBMU

	 IF (K.NE.0) THEN
	 
	    IF (K.LT.0) THEN
	       IND=NT		! At the ground
	    ELSE
	       IND=0		! At TOA
	    ENDIF
    
            Z1=DMAX1(Z1,DABS(I1(IND,K)))
	    
            Z1=DMAX1(Z1,DABS(Q1(IND,K)))
	    
            Z1=DMAX1(Z1,DABS(U1(IND,K)))

	 ENDIF	 
	 
   88 CONTINUE
   
      RETURN

      END  	!END OF PROCEDURE SOS_ARRET_DIFFUS_1

C*============================================================================
C PROCEDURE: SOS_ARRET_DIFFUS_2
C ==========
C      This procedure prepares the stop test on the scatterings when the
c      next terms become very small
C
C
C Description of parameters
C -------------------------- 
C    
C       NBMU (I4)    (I)  Number of actually useful (positive) angles
C
C    I1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)  (I)
C    Q1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)  (I)
C    U1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)  (I)
C      	 Radiance field of order IG versus the level of the
C                          profile and the direction.
C
C    I3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)	
C            Stokes parameter I for the sum of the scattering orders
C            up to the current IG.
C    Q3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)  (I) Stokes parameter Q.
C    U3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)  (I) Stokes parameter U
C
C    Z1   (double)   (O)    Threshold variable for the stop tests.
C
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C None
C
C==============================================================================
      SUBROUTINE SOS_ARRET_DIFFUS_2(NBMU,NT,I1,Q1,U1,I3,Q3,U3,Z1)
      
      IMPLICIT NONE
      

C* Definition of variables
C*------------------------------------------------------------------

      DOUBLE PRECISION Z1	! Threshold variable for the stop tests.

      	! Radiance field of order IG versus the level of the profile
	! and the direction.
      DOUBLE PRECISION I1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	!  Field I.
      DOUBLE PRECISION Q1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	!  Field Q.
      DOUBLE PRECISION U1(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	!  Field U.

	! Functions for the convergence test in geometric series,
	! field emerging at TOA (mu > 0) and downward at the ground (mu < 0).
	! The functions are established for an order S of the expansion
	! in Fourier series, versus the angle mu.
      DOUBLE PRECISION I3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)		
                            ! Stokes parameter I for the sum of the
      				! previous scattering orders up to the current IG.
      DOUBLE PRECISION Q3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
      				! Stokes parameter Q.
      DOUBLE PRECISION U3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
      				! Stokes parameter U.
 
      DOUBLE PRECISION Y	! Variable for the stop tests.

      INTEGER*4 NT	! Number of layers of the atmospheric profile
      INTEGER*4 IND	! Layer number for mu > 0 or mu < 0
      INTEGER*4 K 
      INTEGER*4 NBMU ! Number of actual angles of the simulations

C------------------------------------------------------ 
C* The value of the term of order IG is compared with the sum of the
C* previous terms.
C* Determination of the maximum value of their ratio
C* for the upward field at TOA and the downward field at the ground.
C------------------------------------------------------
      Z1=0.D+00
      DO 88 K=-NBMU,NBMU
      	 
	 IF (K.NE.0) THEN
	 
	    IF (K.LT.0) THEN
	       IND=NT		! At the ground
	    ELSE
	       IND=0		! At TOA
	    ENDIF
	

	    IF (I3(K).NE.0.D+00) THEN
                Y=I1(IND,K)/I3(K)
                Z1=DMAX1(Z1,DABS(Y))
	    ENDIF
	    	 
            IF (Q3(K).NE.0.D+00) THEN
                Y=Q1(IND,K)/Q3(K)
                Z1=DMAX1(Z1,DABS(Y))
	    ENDIF
	 
            IF (U3(K).NE.0.D+00) THEN
                Y=U1(IND,K)/U3(K)
                Z1=DMAX1(Z1,DABS(Y))
	    ENDIF
	 
	 ENDIF !End of test K!=0
	 
	 
   88 CONTINUE
   
      RETURN
      END  	!END OF PROCEDURE SOS_ARRET_DIFFUS_2





C*============================================================================
C PROCEDURE: SOS_ARRET_FOURIER
C ==========
C      This procedure prepares the stop test on the Fourier series
c      decomposition.
c
C
C Description of parameters
C -------------------------- 
C
C       NBMU (I4)    (I)  Number of actually useful (positive) angles
C
C    I4(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)  (I) 	
C    Q4(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)  (I) 
C    U4(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)  (I) 
C        Arrays summing the orders S of the decomposition
C        in Fourier series, versus the angle mu.
C
C    I5(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)  (I) 	
C    Q5(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)  (I) 
C    U5(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)  (I) 	
C        Arrays summing the orders S of the decomposition
C        in Fourier series, with permutation of the sign of the orders S,
C      		     versus the angle mu.
C
C    I3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)	
C            Stokes parameter I for the sum of the scattering orders
C            up to the current IG.
C    Q3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)  (I) Stokes parameter Q.
C    U3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)  (I) Stokes parameter U
C
C    Z1   (double)   (O)    Threshold variable for the stop tests.
C
C
C   
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C None
C
C==============================================================================
      SUBROUTINE SOS_ARRET_FOURIER(NBMU,I3,Q3,U3,I4,Q4,U4,I5,Q5,U5,Z1)
      
      IMPLICIT NONE
      

C* Definition of variables
C*------------------------------------------------------------------

      DOUBLE PRECISION Z1	! Threshold variable for the stop tests.

	! Arrays summing the orders S of the decomposition
        ! in Fourier series, versus the angle mu.
      DOUBLE PRECISION I4(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      DOUBLE PRECISION Q4(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      DOUBLE PRECISION U4(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)

        ! Arrays summing the orders S of the decomposition
        ! in Fourier series, with permutation of the sign of the orders S,
	! versus the angle mu.
      DOUBLE PRECISION I5(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      DOUBLE PRECISION Q5(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      DOUBLE PRECISION U5(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	

	! Field emerging at TOA (mu > 0) and downward at the ground (mu < 0).
	! The functions are established for an order S of the expansion
	! in Fourier series, versus the angle mu.
      DOUBLE PRECISION I3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)		
                            ! Stokes parameter I for the sum of the
      				! previous scattering orders up to the current IG.
      DOUBLE PRECISION Q3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
      				! Stokes parameter Q.
      DOUBLE PRECISION U3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
      				! Stokes parameter U.
 

      DOUBLE PRECISION Y	! Variable for the stop tests.
      
      INTEGER*4 J
      INTEGER*4 NBMU    ! Number of actual angles of the simulations

C------------------------------------------------------------------------- 
C* The stop of the Fourier series decomposition is tested by
C* comparing the order IS (for each direction mu) with
C* the sum of the previous contributions.
C* The maximum value Z1 of this ratio is output for the field emerging
C* at TOA and downward at ground level.
C* If the relative contribution (in absolute value) is below
C* the threshold CTE_PH_SEUIL_SF, the Fourier series decomposition is stopped
C* in the calling program.
C----------------------------------------------------------------------------
      Z1=0.D+00
      DO 169 J=-NBMU,NBMU
         IF (J.EQ.0) GOTO 169

         IF (Q4(J).EQ.0.D+00) GO TO 71
	 
         Y=DABS(Q3(J)/Q4(J))
         Z1=DMAX1(Z1,Y)
	 
   71    IF (I4(J).EQ.0.D+00) GO TO 73
   
         Y=DABS(I3(J)/I4(J))
         Z1=DMAX1(Z1,Y)
	 
   73    IF (U4(J).EQ.0.D+00) GO TO 74
   
         Y=DABS(U3(J)/U4(J))
         Z1=DMAX1(Z1,Y)

   74    IF (Q5(J).EQ.0.D+00) GO TO 771
   
         Y=DABS(Q3(J)/Q5(J))
         Z1=DMAX1(Z1,Y)
	 
  771    IF (U5(J).EQ.0.D+00) GO TO 772
  
         Y=DABS(U3(J)/U5(J))
         Z1=DMAX1(Z1,Y)
	 
  772    IF (I5(J).EQ.0.D+00) GO TO 169
  
         Y=DABS(I3(J)/I5(J))
         Z1=DMAX1(Z1,Y)

  169 CONTINUE

      RETURN
      END  	!END OF PROCEDURE SOS_ARRET_FOURIER





C*========================================================================
C PROCEDURE: SOS_AJOUT_QUEUE
C ==========
C      This procedure adds, if the convergence regime is reached,
C      the tail of the geometric series to the sum of the previous
C      scatterings.
C
C
C Description of parameters
C -------------------------- 
C
C             NT (I4) : (I) Number of layers in the atmospheric profile
C
C       NBMU (I4)    (I)  Number of actually useful (positive) angles
C
C    D1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I) 
C      			  Stokes parameter I for the order IG-1.
C    E1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C      			  Stokes parameter Q for the order IG-1.
C    F1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C      			  Stokes parameter U for the order IG-1.
C
C    G1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C      		          Stokes parameter I for the order IG.
C    H1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C      			  Stokes parameter Q for the order IG.
C    P1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C      			  Stokes parameter U for the order IG.
C
C    D1OUT(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I) 
C                            Stokes parameter I for the order IG-1
C                            (functions for the output profile of the results)
C    E1OUT(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C      			  Stokes parameter Q for the order IG-1.
C    F1OUT(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C      			  Stokes parameter U for the order IG-1.
C
C    G1OUT(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C      		          Stokes parameter I for the order IG.
C                            (functions for the output profile of the results)
C    H1OUT(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C      			  Stokes parameter Q for the order IG.
C    P1OUT(-0:CTE_OS_NT,CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I)
C      			  Stokes parameter U for the order IG.
C
C
C    I3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I/O)	
C            Stokes parameter I for the sum of the scattering orders
C            up to the current IG.
C    Q3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)  (I/O) Stokes parameter Q.
C    U3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)  (I/O) Stokes parameter U
C
C    I3OUT(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I/O)	
C            Stokes parameter I for an order S of the expansion in
C            Fourier series, versus the angle mu (output profile
C            of the results).
C    Q3OUT(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I/O) Stokes parameter Q.
C    U3OUT(0:CTE_OS_NT,-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  (double)  (I/O) Stokes parameter U.
C
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C None
C
C==============================================================================
      SUBROUTINE SOS_AJOUT_QUEUE(NT,NBMU,D1,E1,F1,G1,H1,P1,
     &                           D1OUT,E1OUT,F1OUT,
     &                           G1OUT,H1OUT,P1OUT,
     &                           I3,Q3,U3,I3OUT,Q3OUT,U3OUT)
      
      IMPLICIT NONE
      

C* Definition of variables
C*------------------------------------------------------------------

	! Field emerging at TOA (mu > 0) and downward at the ground (mu < 0).
	! The functions are established for an order S of the expansion
	! in Fourier series, versus the angle mu.
      DOUBLE PRECISION I3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)		
                            ! Stokes parameter I for the sum of the
      				! previous scattering orders up to the current IG.
      DOUBLE PRECISION Q3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
      				! Stokes parameter Q.
      DOUBLE PRECISION U3(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
      				! Stokes parameter U
				
      DOUBLE PRECISION D1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Stokes parameter I for the order IG-1.
      DOUBLE PRECISION E1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Stokes parameter Q for the order IG-1.
      DOUBLE PRECISION F1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      			  	! Stokes parameter U for the order IG-1.

      DOUBLE PRECISION G1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Stokes parameter I for the order IG.
      DOUBLE PRECISION H1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Stokes parameter Q for the order IG.
      DOUBLE PRECISION P1(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      			  	! Stokes parameter U for the order IG.

	! Same functions for the output level of the results:
      DOUBLE PRECISION D1OUT(0:CTE_OS_NT,
     &			     -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Stokes parameter I for the order IG-1.
      DOUBLE PRECISION E1OUT(0:CTE_OS_NT,
     &			     -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Stokes parameter Q for the order IG-1.
      DOUBLE PRECISION F1OUT(0:CTE_OS_NT,
     &			     -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      			  	! Stokes parameter U for the order IG-1.

      DOUBLE PRECISION G1OUT(0:CTE_OS_NT,
     &			     -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Stokes parameter I for the order IG.
      DOUBLE PRECISION H1OUT(0:CTE_OS_NT,
     &			     -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      				! Stokes parameter Q for the order IG.
      DOUBLE PRECISION P1OUT(0:CTE_OS_NT,
     &			     -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)
      			  	! Stokes parameter U for the order IG.


      DOUBLE PRECISION I3OUT(0:CTE_OS_NT,
     &			     -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				!(O)   Stokes parameter I for an order S
      				! of the Fourier series expansion,
				! versus the angle mu.
      DOUBLE PRECISION Q3OUT(0:CTE_OS_NT,
     &			     -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				!(O)   Stokes parameter Q.
      DOUBLE PRECISION U3OUT(0:CTE_OS_NT,
     &			     -CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				!(O)   Stokes parameter U.


      DOUBLE PRECISION I1QUEUE	
      		! Estimated value of the tail of the geometric series for I.
      DOUBLE PRECISION Q1QUEUE
      		! Estimated value of the tail of the geometric series for Q.
      DOUBLE PRECISION U1QUEUE
      		! Estimated value of the tail of the geometric series for U.

      INTEGER*4 I,J
      INTEGER*4 NBMU    ! Number of actual angles of the simulations
      INTEGER*4 NT	! Number of layers of the atmospheric profile

C* For the field emerging at TOA and downward at the ground
C------------------------------------------------------
      DO 610 J=-NBMU,NBMU
      
         IF (J.EQ.0) GOTO 610
	
         IF (D1(J).EQ.0.) THEN
            I1QUEUE = 0.
         ELSE	
            I1QUEUE = G1(J) / ( 1-G1(J)/D1(J) )
         ENDIF
	 	
         IF (E1(J).EQ.0.) THEN
            Q1QUEUE = 0.
         ELSE	
            Q1QUEUE = H1(J) / ( 1-H1(J)/E1(J) )
         ENDIF	
	
         IF (F1(J).EQ.0.) THEN
            U1QUEUE = 0.
         ELSE	
            U1QUEUE = P1(J) / ( 1-P1(J)/F1(J) )
         ENDIF	

         I3(J)=I3(J)+I1QUEUE
         Q3(J)=Q3(J)+Q1QUEUE
         U3(J)=U3(J)+U1QUEUE
	
  610 CONTINUE

C* For the upward and downward fields at the output level of the results
C------------------------------------------------------------------------
      DO 612 J=-NBMU,NBMU
      
         IF (J.EQ.0) GOTO 612
	
	  DO I=0,NT
	 
         	IF (D1OUT(I,J).EQ.0.) THEN
            	     I1QUEUE = 0.
         	ELSE	
            	     I1QUEUE = G1OUT(I,J) / ( 1-G1OUT(I,J)/D1OUT(I,J) )
         	ENDIF
	 	
         	IF (E1OUT(I,J).EQ.0.) THEN
            	     Q1QUEUE = 0.
         	ELSE	
            	     Q1QUEUE = H1OUT(I,J) / ( 1-H1OUT(I,J)/E1OUT(I,J) )
         	ENDIF	
	
         	IF (F1OUT(I,J).EQ.0.) THEN
            	     U1QUEUE = 0.
         	ELSE	
            	     U1QUEUE = P1OUT(I,J) / ( 1-P1OUT(I,J)/F1OUT(I,J) )
         	ENDIF	

         	I3OUT(I,J)=I3OUT(I,J)+I1QUEUE
         	Q3OUT(I,J)=Q3OUT(I,J)+Q1QUEUE
         	U3OUT(I,J)=U3OUT(I,J)+U1QUEUE
	   
	  ENDDO  
	 	
  612 CONTINUE

      RETURN
      END  	!END OF PROCEDURE SOS_AJOUT_QUEUE
