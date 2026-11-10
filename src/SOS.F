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


C******************************************************************************
C* FILE: SOS.F
C* PROJECT: Successive orders of scattering
C* ROLE: Main module of the OS code
C*           - reads the atmospheric profile and adapts it to a possible
C*             truncation of the aerosol phase function,
C*           - solves the radiative transfer equation by successive
C*             orders of interaction (scattering/reflection),
C*           - recombines the Stokes parameters from their
C*             decomposition in Fourier series.
C*
C* AUTHOR:
C*   Initial code: Laboratoire d'Optique Atmospherique (LOA, Villeneuve d'Ascq).
C*   Code rework and development: CS GROUP France, B.Lafrance, X.Lenot
C* DATE: 30/04/2002
C*
C* MOD:VERSION:1.0: Upgrade of the LOA OS code
C*             2.0: Reading of the aerosol and profile parameters,
C*                  Update of the comments.
C*             2.1: 07/07/2008 
C*                  Change of the reading format of the parameters WIND and TETAS
C*                  to correct a decoding error for integer values.
C*			WIND changed from F4.1 to F4
C*                   TETAS changed from F9.5 to F9
C*
C* MOD:VERSION:3.0: 24/11/2009 
C*       - Change of the argument passing:
C*       - passing by pairs (keyword, value).
C*       - adaptation of the error case management.
C*       - Error messages converted to English.
C*
C*
C* MOD:VERSION:3.1: 22/12/2009 
C*                  Change of the reading format of the parameters Alpha and Beta of the Nadal model
C*                  to correct a decoding error for integer values.
C*			ALPHA changed from F6.4 to F6
C*			BETA  changed from F4.1 to F4
C*
C*                  
C* MOD:VERSION:4.0: 10/01/2010 
C*
C*   - Removal of the constant CTE_OS_FICANGLE (from SOS.h)
C*     in favour of passing the name of the angle file (produced by SOS_ANGLES) as argument.
C*
C*   - Change of the size of the radiance arrays:
C*     ==> Changed from CTE_OS_NBMU to CTE_OS_NBMU_MAX.
C*     ==> Loop on the angles up to LUM_NBMU, equal to the actual number of angles used.
C*
C*   - Removal of the constant CTE_OS_NB (from SOS.h)
C*     in favour of passing the value OS_NB (read in the angle file) as argument
C*
C*   - Change of the size of the arrays of the decomposition in Legendre polynomials:
C*     ==> Changed from CTE_OS_NB to CTE_OS_NB_MAX.
C*     ==> Loop on the orders up to OS_NB.
C* 
C*   - Angle management:
C*         * Removal of the reading of a Gauss angle file
C*         * Added call to routine SOS_ANGLES to compute the Gauss angles and their weights
C*           and to add user angles (with null weights): useful angles = LUM_NBMU.
C*         * Removal of the reading of SOS.Imus and SOS.thetas: these parameters
C*           are obtained by reading the content of the angle file.
C*         * Added optional input parameters SOS.ResFileUp.UserAng and SOS.ResFileDown.UserAng
C*           to store the output radiances only for the "user" angles
C*           --> Adaptation of routine SOS_TRPHI_OPTION.
C*
C*   - Change of the size of the radiance arrays:
C*     ==> Changed from CTE_OS_NBMU to CTE_OS_NBMU_MAX.
C*     ==> Loop on the angles up to LUM_NBMU, equal to the actual number of angles used.
C*
C*   - Routine SOS_TRPHI_OPTION:
C*         * Change of the reading line of the arrays of
C*           the reflection matrix at order IS of the Fourier series decomposition:
C*           ==> limited to the useful bounds -LUM_NBMU:LUM_NBMU
C*         * Passing of the array IND_ANGOUT defining the angles for which an output
C*           is required in the optional files SOS.ResFileUp.UserAng
C*           and SOS.ResFileDown.UserAng
C*
C*   - Adaptation of routines SOS_OS & SOS_TRPHI_OPTION to the passing of the parameter
C*     LUM_NBMU (nb of useful angles) and OS_NB (max order of the FS expansions).
C*
C*   - Adaptation to the passing of the parameter FICTRANS (-SOS.Trans) to output the
C*     transmission values.
C*
C* MOD:VERSION:4.1: 24/06/2010 
C*   -  FORMAT labels ,X, changed to ,1X, to avoid compilation errors.
C*   -  Fix to delete the file FICOS_TRANS_TMP at the end of the transmission computations.
C*
C* MOD:VERSION:6.0: 06/06/2013
C*   -  The SOS.F function of the previous version (V 5.1) has been split into 3 in
C*	 this new version: SOS_PREPAOS.F, SOS.F and SOS_FINOS.F. The function SOS.F
C*	 is run for each CKD exponential, while SOS_PREPAOS.F and SOS_FIN_OS.F
C*	 are run before and after the loop on the exponentials, respectively (see SOS_MAIN.F)
C*
C*   -  The input parameters are given as a list of parameters
C*      and no longer as pairs with associated keyword => all
C*	 parameters are therefore mandatory
C*
C*   -  Change of the value of Pi, declared as a constant in the file
C*      (INCTE_PI replaces CTE_PI)
C*
C* MOD:VERSION:6.1: 29/05/2020
C*   -  Adjustment of the format of the atmospheric profile file to the improved
C*      precision of the profile (format label 70)
C*
C* MOD:VERSION:6.2: 07/07/2022
C*    - Constants renamed from SOS_* to CTE_*
C*
C******************************************************************************

C----------------------------------------------------------------------------
C Definition of constants
C---------------------------------------------------------------------------- 
C Constants used:
C    CTE_LENDIR : Length of directory names.
C    CTE_LENFIC2 : Length of file names with directory tree.
C    CTE_OS_NBMU_MAX : Maximum number of positive angles for solving the RTE
C    CTE_OS_NB_MAX :  Maximum value of the limit order of the expansion in
C                     Legendre polynomials.
C    CTE_OS_NT  : Number of layers of the atmospheric profile.
C----------------------------------------------------------------------------
#include "SOS.h"
#define INCTE_PI DACOS(-1.D+00)

C==============================================================================
C PROGRAM: SOS
C ==========      
C      Main program of the OS code
C
C Input data of the program
C ------------------------------
C    -->Parameters/arguments of the program
C
C             FICOS (CHARACTER) : (I) Complete name of the OS result file.
C                                 (directory + file name + extension)
C
C             FICTRANS (CHARACTER): (I) complete name of the file of transmission values
C                                   (directory + file name + extension)
C
C             FICPROFIL(CHARACTER) : (I) complete name of the PROFIL result file
C                                    (directory + file name + extension)
C
C             NT (I4) : (I) Number of layers in the atmospheric profile
C
C      	ZOUT (F9) : (I) Altitude of the output level
C			         -1 =>  Standard output:
C				         Radiance emerging at TOA
C				         + downward radiance at ground level.
C			          z =>  Upward and downward radiance
C				         for altitude z (0 <= z <= CTE_TOA_ALT).
C			
C             IGMAX (I4) : (I) Maximum order of interaction.
C		 
C             IPOLAR (I4) : (I) Option to cut the polarisation.
C                               0 : simulation without polarisation.
C			           1 : normal simulation with polarisation.
C
C             RON  (DOUBLE) : (I) (Molecular) depolarisation factor
C                             (decoding format: real F9.5)
C
C             IND_SURF (DOUBLE) : (I) Refractive index surface / air.
C
C             RHO (DOUBLE) : (I) Surface albedo
C
C             IMAT_SURF (I4) :(I) Index specifying whether the simulations
C      			      use a BRDF or BPDF reflection model:
C			             1 : for BRDF or BPDF computation (N0 > 0 must be used).
C			             otherwise : no BRDF or BPDF simulation.
C
C             IFRESNEL (I4) : (I) Index specifying whether the simulations
C      			      include the Fresnel reflection on a flat sea:
C			             1 : for a computation of the Fresnel reflection.
C			             Otherwise : not taken into account.
C        
C             FICSURF (CHARACTER) : ==> Associated with: ISURF = 1 3 4 5 or 6 (see SOS_PREPA_OS)
C					complete name of the file containing the
C                   			surface reflection matrix (Glitter file or other BRDF/BPDF).
C                   			(directory + file name + extension)
C                    		> DEFAULT if a file of the BRDF / BPDF database is used.
C
C      	N0 (I4) : (I) Number of the Gauss mu for the solar zenith angle.
C			      (complementary angle: mus = -rmu(N0))
C			      (N0 > 0 => to take it into account).
C
C      	PIZ (DOUBLE) : (I) Single scattering albedo.
C
C      	PIZTR (DOUBLE) : (I) Single scattering albedo after truncation
C
C      	A (DOUBLE) : (I) Truncation coefficient
C
C      	RMU (DOUBLE) : (I) 1D array of 2*CTE_OS_NBMU_MAX + 1 elements
C			            containing the cosine of the Gauss angles and of the
C			            angles provided by the user
C
C      	GA (DOUBLE) : (I) 1D array of 2*CTE_OS_NBMU_MAX + 1 elements
C			           containing the weight associated with each Gauss angle and each
C			           angle provided by the user
C 
C             TETAS (DOUBLE) : (I) Solar zenith angle in degrees (0 < tetas < 90).
C      		               (complementary angle: mus= -cos(tetas))
C
C             OS_NB (I4) : (I) Order of the expansions of the phase functions
C			          in Legendre polynomials.
C
C             LUM_NBMU (I4) : (I) Number of actual angles for the radiance computations
C
C             ALPHA (DOUBLE) : (I) 1D array of CTE_OS_NB_MAX + 1 elements
C				       containing the Alpha coefficients of the decomposition in Legendre
C				       functions of the elements Pij of the phase function
C
C             BETA (DOUBLE) : (I) 1D array of CTE_OS_NB_MAX + 1 elements
C				      containing the Beta coefficients of the decomposition in Legendre
C				      functions of the elements Pij of the phase function
C
C             GAMMA (DOUBLE) : (I) 1D array of CTE_OS_NB_MAX + 1 elements
C				       containing the Gamma coefficients of the decomposition in Legendre
C				       functions of the elements Pij of the phase function
C
C             ZETA (DOUBLE) : (I) 1D array of CTE_OS_NB_MAX + 1 elements
C				      containing the Zeta coefficients of the decomposition in Legendre
C				      functions of the elements Pij of the phase function
C
C             TTOT_TRONC (DOUBLE) : (O) Total optical thickness at ground (truncated)
C
C             TTOT_VRAI (DOUBLE) : (O) Total optical thickness at ground (not truncated)
C
C             TAUOUT (DOUBLE) : (O) Optical thickness at the output level for the upward field
C			               (truncated)
C
C             EMOINS  (DOUBLE): (O) Diffuse flux illuminating the ground, normalised by
C                   			 the incident solar flux at TOA.
C             
C             EPLUS  (DOUBLE) : (O) Diffuse flux emerging at TOA, normalised by
C                   			 the incident solar flux at TOA.
C
C             TDIFMUS (DOUBLE) : (O) Diffuse transmission for the solar incidence direction
C      				 of the equivalent atmosphere
C
C             TDIFMUG (DOUBLE) : (O) 1D array of 2*CTE_OS_NBMU_MAX + 1 elements
C			     	         containing the diffuse transmission of the equivalent atmosphere
C				         for the Gauss angles and the angles provided by the user
C
C	       TRACE (LOGICAL) : (I) If true => trace files are written
C
C	       IDLOG (I4) : (I) Identifier of the log file
C
C             IER (I4) : (O) error code =0 if no error, =-1 otherwise
C
C    -->Files
C       The surface reflection matrix file
C       The Profile file
C       The transmission file
C
C
C Results provided by the program
C ----------------------------------
C   OS result file (Fourier series decomposition)
C   Diffuse transmissions
C   Optical thicknesses (ground and output level)
C   Trace file (optional) (logical unit number of the trace file =IDLOG)
C
C   If one or several files with the same names already exist, they are overwritten
C
C
C 
C Description of the files used
C ---------------------------------
C
C   -->Content of the BRDF/BPDF file
C       Unformatted binary file, created by SOS_SURFACE
C       It contains the surface reflection matrix
C       (Glitter file or other BRDF/BPDF):
C       i.e. the elements Pij of the reflection matrix:
C              P11(LUM_NBMU,LUM_NBMU) (real),
C              P12(LUM_NBMU,LUM_NBMU) (real),
C              P13(LUM_NBMU,LUM_NBMU) (real),
C              P21(LUM_NBMU,LUM_NBMU) (real),
C              P22(LUM_NBMU,LUM_NBMU) (real),
C              P23(LUM_NBMU,LUM_NBMU) (real),
C              P31(LUM_NBMU,LUM_NBMU) (real),
C              P32(LUM_NBMU,LUM_NBMU) (real),
C              P33(LUM_NBMU,LUM_NBMU) (real)
C
C        --> Records on the domain limited to the useful angles
C                                   --------------------------------
C
C   -->Content of the OS result file (FICOS)
C       Unformatted binary file, created by SOS_OS
C       It contains the Stokes parameters Q, U and I for an order S
C       of the Fourier series expansion, versus the angle mu
C       (double precision values).
C
C   -->Content of the profile file
C	ASCII file containing, for each level of the atmospheric profile:
C	The level index, the molecules + aerosols optical thickness,
C	the absorption optical thickness, the total optical thickness,
C	the percentage of aerosols and the percentage of molecules
C
C Common used:
C --------------
C  None
C
C
C Error cases:
C ------------        
C     -Incorrect decoding of parameters
C     -Error in a subroutine
C     -Error when opening the trace file
C     -Invalid parameters
C
C     Message displayed on screen, program stop and
C     return of status 1 to the shell
C   
C
C  Subroutines used:
C --------------------------
C  This program calls the routines:
C     - SOS_OS
C
C============================================================================
       SUBROUTINE SOS(FICOS,FICTRANS,FICPROFIL,
     &		      NT,ZOUT,IGMAX,IPOLAR,RON,           
     &	             IND_SURF,RHO,IMAT_SURF,IFRESNEL,FICSURF,N0,PIZ,
     &		      PIZTR,A,RMU,GA,TETAS,OS_NB,LUM_NBMU,ALPHA,BETA,
     &		      GAMMA,ZETA,TTOT_TRONC,TTOT_VRAI,TAUOUT,TDIFMUS,
     &		      TDIFMUG,EMOINS,EPLUS,TRACE,IDLOG,IER)      
      
      IMPLICIT NONE

C* Definition of variables
C*-----------------------------------------------------------------
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
                         !Cosines of the Gauss (and user) angles.
      DOUBLE PRECISION GA(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
                         !Weights of the Gauss (and user) angles.
 
 
      DOUBLE PRECISION TETAS ! Solar zenith angle in degrees (0 < tetas < 90).
      				 ! (complementary angle: mus= -cos(tetas))

      
      DOUBLE PRECISION IND_SURF  ! Refractive index surface / air.
           
      DOUBLE PRECISION RON	! (Molecular) depolarisation factor.

      DOUBLE PRECISION RHO	! Surface albedo.
      
      DOUBLE PRECISION A	! Truncation coefficient.
      DOUBLE PRECISION PIZ	! Single scattering albedo.
      DOUBLE PRECISION PIZTR	! Single scattering albedo after truncation.

      DOUBLE PRECISION ALPHA(0:CTE_OS_NB_MAX)	! Parameters alpha(k).
      DOUBLE PRECISION BETA(0:CTE_OS_NB_MAX)	! Parameters beta(k).
      DOUBLE PRECISION GAMMA(0:CTE_OS_NB_MAX)	! Parameters gamma(k).
      DOUBLE PRECISION ZETA(0:CTE_OS_NB_MAX)	! Parameters zeta(k).

      DOUBLE PRECISION VA,VATR
      DOUBLE PRECISION VR
      DOUBLE PRECISION VG
      
      DOUBLE PRECISION H(0:CTE_OS_NT)	! Optical thickness for each level
      DOUBLE PRECISION Hma(0:CTE_OS_NT)	! Molecules + aerosols optical thickness for each level
      DOUBLE PRECISION Ha(0:CTE_OS_NT)	! Absorption optical thickness for each level
      DOUBLE PRECISION HTR(0:CTE_OS_NT)	! Truncated optical thickness for each level
      DOUBLE PRECISION XDEL(0:CTE_OS_NT)  ! Percentage of aerosols per level
      DOUBLE PRECISION YDEL(0:CTE_OS_NT)  ! Percentage of molecules per level
      DOUBLE PRECISION ZPROF(0:CTE_OS_NT) ! Altitudes of the profile levels (in km)

      DOUBLE PRECISION TTOT_TRONC ! Total optical thickness at ground (truncated).
      DOUBLE PRECISION TTOT_VRAI  ! Total optical thickness at ground (not truncated).
      DOUBLE PRECISION TAUOUT     ! Optical thickness at the
      				      ! output level of the results for the upward field (truncated).

      DOUBLE PRECISION EMOINS   ! Diffuse flux illuminating the ground, normalised by
      				    ! the incident solar flux at TOA.

      DOUBLE PRECISION EPLUS    ! Diffuse flux emerging at TOA, normalised by
      				    ! the incident solar flux at TOA.

      DOUBLE PRECISION EPLUS_SOLNOIR    ! Diffuse flux emerging at TOA, normalised by
      				            ! the incident solar flux at TOA, for a black ground.
								
      DOUBLE PRECISION TDIFMUS  !Diffuse transmission for the solar incidence direction
      				    !of the equivalent atmosphere
				               
      DOUBLE PRECISION TDIFMUG(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
      				!Diffuse transmission for the Gauss angles of the equivalent atmosphere

      DOUBLE PRECISION TDIFMUS_TMP  !Diffuse transmission for each Gauss angle
      				        !of the equivalent atmosphere

      DOUBLE PRECISION  ZZ      ! Linear interpolation coefficient between 2 altitudes

      INTEGER*4 NT      ! Number of layers in the profile
     
      INTEGER*4 N0	! Number of the Gauss mu for the solar
      			! zenith angle.
			! (complementary angle: mus = -rmu(N0))
			! (N0 > 0 => to take it into account).

				
      DOUBLE PRECISION ZOUT	! Output type:
				!      -1 =>  Standard output:
				!             Radiance emerging at TOA
				!           + downward radiance at ground level.
				!       z =>  Upward and downward radiance
				!             for altitude z (0 <= z <= CTE_TOA_ALT).
			

      INTEGER*4 IGMAX	! Maximum order of scattering.
      
      INTEGER*4 IPOLAR	! Option for simulation without polarisation.
      			       ! 0 : computations without polarisation.
			       ! 1 : normal computations with polarisation.

      INTEGER*4 ITRPHI	! Index of output type.
			       !   1 =>  Viewing plane for one azimuth value.
			       !   2 =>  Polar diagram versus the azimuth
			       !         and viewing directions.
			
      INTEGER*4 PAS_PHI	! Azimuth step in degrees (for ITRPHI=2)

      INTEGER*4 IGLI	! Index specifying whether the simulations
      			! include the Glitter:
			!    1 => for Glitter computation (N0 > 0 must then be used).
			!    Otherwise => no Glitter simulation.
	
      INTEGER*4 IFRESNEL ! Index specifying whether the simulations
      			    ! include the Fresnel reflection on a flat sea:
			    !    1 => for a computation of the Fresnel reflection.
			    !    Otherwise => not taken into account.
	
		 
	
      INTEGER*4 IMAT_SURF ! Index specifying whether the simulations
      			     ! use a BRDF or BPDF reflection model:
			     !    1 => for BRDF or BPDF computation (N0 > 0 must be used).
			     !    Otherwise => no BRDF or BPDF simulation.
			  			  			 											
      INTEGER*4 IBORM ! Maximum order of the Fourier series expansion
      			 ! being processed:
			 !   _ 2 for a purely molecular atmosphere,
			 !   _ NB in presence of aerosols.

      INTEGER*4 IER     ! error code =0 if no error, =-1 otherwise

      INTEGER*4 K       ! Array index
      INTEGER*4 I,J     ! Array index
      
      INTEGER*4 LDIRSURF ! length of the character string REPSURF.
      INTEGER*4 LFICSURF ! length of the character string FICSURF.
      
      INTEGER*4 NBMU_GAUSS      ! Number of Gauss angles used
      INTEGER*4 LUM_NBMU        ! Number of actual angles of the simulations
      INTEGER*4 OS_NB 	 ! Order of the expansions of the phase functions
			        ! in Legendre polynomials.
      
      LOGICAL TRACE     ! = true if writing in the trace file
      LOGICAL LTA       ! = true if the aerosol optical thickness is null
      LOGICAL EX	! Checks the existence of a file

      CHARACTER*CTE_LENFIC2 FICOS       ! Name of the OS result file.
      CHARACTER*CTE_LENFIC2 FICOS_TRANS_TMP   ! Name of the temporary OS file (for the transmissions)
      
      CHARACTER*CTE_LENFIC2 FICSURF ! File containing the surface reflection
      				        ! matrix (Glitter file or other BRDF/BPDF).
      
      CHARACTER*CTE_LENFIC2 FICPROFIL   ! name of the PROFIL input file

      CHARACTER*CTE_LENFIC2 FICTRANS    ! Name of the optional output file of the transmissions.

             
      INTEGER*4 IDLOG	        ! Identifier number of the trace file
     	  
      	        		 

C* Initialisation
C-----------------
      IER=0
      LTA=.TRUE.			!null aerosol optical thickness
      VA=0.
      VATR=0.
      VR=0.


	 
      
C* Reading of the PROFIL file
C  and storage of the non-truncated optical thickness
C  -----------------------------------------------
      
      OPEN(2,FILE=FICPROFIL,STATUS='OLD',ERR=912)
      
      DO I=0,NT
c         READ(2,70,err=922)J,ZPROF(I),Hma(I),Ha(I),H(I),XDEL(I),YDEL(I)
         READ(2,70,err=922)J,ZPROF(I),H(I),XDEL(I),YDEL(I)
      ENDDO

      TTOT_VRAI=H(NT)
      
            
C* Adaptation of the profile to the truncation if A is not null
C  -------------------------------------------------
      HTR(0)=H(0)

      IF(A.NE.0.) THEN
         DO I=1,NT	
	    VA = XDEL(I)*(H(I)-H(I-1))
	    VATR = VA*(1-PIZ*0.5*A)
	    VR = YDEL(I)*(H(I)-H(I-1))
	    VG = (1-XDEL(I)-YDEL(I))*(H(I)-H(I-1))
	    
	    HTR(I) = (VATR+VR+VG) + HTR(I-1)
	    
	    XDEL(I) = VATR/(VATR+VR+VG)
	    YDEL(I) = VR/(VATR+VR+VG)
         ENDDO
      ENDIF

      DO I=0,NT
      	 IF(A.NE.0.) H(I)=HTR(I)	! Case with truncation
         XDEL(I)=XDEL(I)*PIZTR
      	 IF(XDEL(I).NE.0.) LTA=.FALSE.	!TA not null
      ENDDO

      
C* Maximum order of the Fourier series decomposition
C* (2 for pure Rayleigh)
C---------------------------------------------------------
      IBORM=OS_NB    
      IF (LTA) IBORM=2
	
C* Call of the OS
C------------------
      CALL SOS_OS(LUM_NBMU,RMU,GA,OS_NB,NT,
     &            FICSURF,FICOS,
     &            N0,TETAS,RHO,IMAT_SURF,IFRESNEL,IND_SURF,
     &            H,XDEL,YDEL,ZPROF,RON,
     &            ALPHA,BETA,GAMMA,ZETA,
     &            ZOUT,IGMAX,IBORM,IPOLAR,TRACE,IDLOG,
     &            EMOINS,EPLUS,IER)
      IF (IER.NE.0) GOTO 995   ! error case


C* Definition of the optical thickness at the profile level
C* for the output of the upward field
C-------------------------------------------------------
      IF (ZOUT.EQ.-1) THEN
	       TAUOUT=H(0)		!Toa.
      ELSE
      
c Search for the altitude level below ZOUT
	    J=1
	    DO WHILE(ZOUT.LT.ZPROF(J))
		   J=J+1
	    ENDDO

C Position of the altitude with respect to the previous and next levels
	    ZZ = (ZOUT-ZPROF(J-1))/(ZPROF(J)-ZPROF(J-1))

C* Computation of the downward and upward fields
            TAUOUT=(1-ZZ)*H(J-1)+ZZ*H(J)
	    
      ENDIF

C* Optical thickness at ground level (possibly truncated)
      TTOT_TRONC=H(NT)
      
C* Computation of the direct and diffuse transmissions
C-----------------------------------------------

C     No writing of the results (I,Q,U) in a file: FICOS_TRANS_TMP="NO_OUTPUT"
C     Conditions to simulate a diffuse transmission:
C        RHO = 0 : black ground (no Lambertian reflection)
C        IMAT_SURF = 0 : black ground (no specular reflection)
C        IFRESNEL = 0 : black ground (no specular reflection)
C        ZOUT = -1 : output for down field = ground
C        IBORM = 0 : expansion limited to order IS=0
C        TRACE = .FALSE. : the TRACE file is not completed
C
C     The diffuse transmissions output by SOS_OS are given for
C     the equivalent atmosphere (after truncation of the aerosol phase
C     function). They are converted to the value corresponding to
C     the true atmosphere.

      IF (FICTRANS.NE.'NO_OUTPUT') THEN

	  FICOS_TRANS_TMP = "NO_OUTPUT"	!File not created

C     --> Diffuse transmission TOA --> Ground for the solar incidence direction
	  	 	  
	  CALL SOS_OS(LUM_NBMU,RMU,GA,OS_NB,NT,FICSURF,
     &                FICOS_TRANS_TMP,
     &                N0,TETAS,0.D+00,0,0,IND_SURF,
     &                H,XDEL,YDEL,ZPROF,RON,ALPHA,BETA,GAMMA,ZETA,
     &                -1.D+00,IGMAX,0,IPOLAR,.FALSE.,
     &		        IDLOG,TDIFMUS,EPLUS_SOLNOIR,IER)
          IF (IER.NE.0) GOTO 995   ! error case
	
	  
C     --> Diffuse transmission Ground --> TOA:
C         Loop on the values IMUS: N0 = J	for computation TOA --> Ground (reciprocity)
	  DO J=1,LUM_NBMU
	  
	  
             CALL SOS_OS(LUM_NBMU,RMU,GA,OS_NB,NT,FICSURF,
     &                   FICOS_TRANS_TMP,
     &                   J,TETAS,0.D+00,0,0,IND_SURF,
     &                   H,XDEL,YDEL,ZPROF,RON,ALPHA,BETA,GAMMA,ZETA,
     &                   -1.D+00,IGMAX,0,IPOLAR,.FALSE.,
     &			    IDLOG,TDIFMUS_TMP,EPLUS_SOLNOIR,IER)
             IF (IER.NE.0) GOTO 995   ! error case
	      
	     TDIFMUG(J)= TDIFMUS_TMP
	     
          ENDDO	    
	  
      ENDIF !End of test whether the transmission file is to be created
      
C* File closing
C-------------------
      CLOSE(1)
      CLOSE(2)
      CLOSE(20)
      
C* Format 
C----------
  100 FORMAT(17hNB_TOTAL_ANGLES :,I4)
  110 FORMAT(17hNB_GAUSS_ANGLES :,I4)
  120 FORMAT(17hANGLES_USERFILE :,A)  
  200 FORMAT(16hINTERNAL_OS_NB :,I4)
  410 FORMAT(20hSOLAR ZENITH ANGLE :,F7.3)
  420 FORMAT(15hINTERNAL_IMUS :,I4)
  610 FORMAT(I4,1X,2D21.14,1X,I4)

101   FORMAT(12h    thetas =,1x,f6.3,3x,12htd(thetas) =,1x,f7.4)  
201   FORMAT(12h    thetav =,1x,f6.3,3x,12htd(thetav) =,1x,f7.4)     

C* Nominal end
C-------------------
      GOTO 9999
      
C* Error case and return of status 1 to the shell
C----------------------------------------------

            
  912 WRITE(6,*) '  ERROR on PROFILE file opening for SOS'
       IER=-1
      GOTO 9999 
      
  913 WRITE(6,*) '  ERROR on AEROSOLS file opening for SOS'
       IER=-1
      GOTO 9999    
      
  922 WRITE(6,*) '  ERROR on PROFILE file reading for SOS'
       IER=-1
      GOTO 9999
      
  923 WRITE(6,*) '  ERROR on AEROSOLS file reading for SOS'
       IER=-1
      GOTO 9999
      
  995 WRITE(6,*) '  ERROR on subroutine SOS_OS'
       IER=-1
      GOTO 9999
                                                                                         
                               
C* Format
C---------   
   50 FORMAT(2X,F9.5)
   60 FORMAT(4(2X,E16.8))
c   70 FORMAT(2X,I5,F10.5,5(F13.7))
   70 FORMAT(2X,I5,F10.5,3(E15.8))
  223 FORMAT(I4,1X,D21.14,2X,F7.3,2X,D21.14,2X,I2)
    
 9999 RETURN
 
      END	!END OF PROGRAM SOS




