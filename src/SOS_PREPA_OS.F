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
C* FILE: SOS_PREPA_OS.F
C* PROJECT: Successive orders of scattering
C* ROLE: Main module of the OS code
C*           - reads the parameters of the aerosol size distributions,
C*           - reads the angular parameters
C*
C* AUTHOR:
C*   Initial code: Laboratoire d'Optique Atmospherique (LOA, Villeneuve d'Ascq).
C*   Code rework and development: CS GROUP France, B.Lafrance, X.Lenot
C*
C* DATE: 24/06/2013
C* MOD:VERSION:1.0: File created by splitting the SOS.F function into 3
C*		      sub-functions: SOS.F, SOS_PREPA_OS.F and SOS_FIN_OS.F
C*
C* MOD:VERSION:1.1: 24/02/2016 
C*    - Minor adjustments to strictly respect the 72 columns
C*
C*    - Change of the way the formatted files are read (Angles, Aerosols)
C*      to be compatible with the gfortran compiler.
C*
C* MOD:VERSION:1.2: 17/07/2020 
C*    - Adjustments to introduce the Maignan model.
C*
C* MOD:VERSION:1.3: 07/07/2022
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
C
C Specific constants:
C    LENLINE  : Maximum length of a data line read in a file
C----------------------------------------------------------------------------
#include "SOS.h"
#define INCTE_PI DACOS(-1.D+00)
#define LENLINE 1000

C==============================================================================
C PROGRAM: SOS_PREPA_OS
C ==========      
C      Program initialising the variables for the OS computation: reading of the
C      angle and aerosol files. Identification of the computation mode
C      of the surface reflectance and of the associated file.
C
C Input data of the program
C ------------------------------
C    -->Parameters/arguments of the program
C
C            FICANGLES (CHARACTER) : (I) complete name of the file of the angles used
C                                      for the radiance and BRDF computations
C                 		           (directory + file name + extension)
C
C	      FICGRANU (CHARACTER) : (I) complete name of the GRANU result file
C                                      (directory + file name + extension)
C
C            ISURF  (I4)    : (I) Surface type ISURF
C			            (decoding format: integer I4)
C      		0 :Lambertian surface of albedo RHO.
C			1 :Lambertian surface + glitter.
C			2 :Lambertian surface + Fresnel reflection.
C		       3 :Lambertian surface + Roujean BRDF.
C			4 :Lambertian surface + Roujean BRDF + Rondeaux BPDF.
C			5 :Lambertian surface + Roujean BRDF + Breon BPDF.
C			6 :Lambertian surface + Roujean BRDF + Nadal BPDF.
C
C            IND_SURF (I4) : (I) Refractive index surface / air.
C
C            WIND (DOUBLE) : (I) ==> Associated with: ISURF = 1
C 	      		       Wind speed in m/s.
C		
C            K0_ROUJEAN (DOUBLE) : (I) ==> Associated with: ISURF >= 3
C 	       		      Parameter k0 of the Roujean BRDF model.
C		
C            K1_ROUJEAN (DOUBLE) : (I) ==> Associated with: ISURF >= 3
C 	      			     Parameter k1 of the Roujean BRDF model.
C		
C            K2_ROUJEAN (DOUBLE) : (I) ==> Associated with: ISURF >= 3
C 	      			     Parameter k2 of the Roujean BRDF model.
C
C            ALPHA_NADAL (DOUBLE) : (I) ==> Associated with: ISURF = 6
C 	      			     Parameter alpha of the Nadal BPDF model.
C
C            BETA_NADAL (DOUBLE) : (I) ==> Associated with: ISURF = 6
C 	      			     Parameter beta of the Nadal BPDF model.
C
C            COEF_C_MAIGNAN (DOUBLE) : (I) ==> Associated with: ISURF = 7
C 	      		            Parameter C*exp(-NDVI) of the Maignan BPDF model.
C
C
C            REPSURF (CHARACTER) : (I) ==> Associated with: ISURF= 1 3 4 5 or 6 and FICSURF= "DEFAULT"
C             			     Storage directory of the surface reflection matrix
C                   		            files (Glitter file or other BRDF/BPDF).
C
C            FICSURF (CHARACTER) : (I) ==> Associated with: ISURF= 1 3 4 5 or 6
C             			    complete name of the file containing the surface
C                   		           reflection matrix (Glitter file or other BRDF/BPDF).
C                   		           (directory + file name + extension)
C                    		    => DEFAULT if a file of the BRDF / BPDF database is used.
C
C            TRACE (LOGICAL) : (I) If true => trace files are written
C
C	      IDLOG (I4) : (I) Identifier of the log file
C
C            LUM_NBMU (I4) : (O) Number of actual angles for the radiance computations
C
C            N0 (I4) : (O) Number of the Gauss mu for the solar zenith angle.
C			     (complementary angle: mus = -rmu(N0))
C			     (N0 > 0 => to take it into account).
C
C            PIZ (DOUBLE) : (O) Single scattering albedo.
C
C            PIZTR (DOUBLE) : (O) Single scattering albedo after truncation
C
C            A (DOUBLE) : (O) Truncation coefficient
C
C            RMU (DOUBLE) : (O) 1D array of 2*CTE_OS_NBMU_MAX + 1 elements
C			           containing the cosine of the Gauss angles and of the
C			           angles provided by the user
C
C            GA (DOUBLE) : (O) 1D array of 2*CTE_OS_NBMU_MAX + 1 elements
C			          containing the weight associated with each Gauss angle and each
C			          angle provided by the user
C 
C            OS_NB (I4) : (O) Order of the expansions of the phase functions
C			         in Legendre polynomials.
C
C            IMAT_SURF (I4) :(O) Index specifying whether the simulations
C      			     use a BRDF or BPDF reflection model:
C			            1 : for BRDF or BPDF computation (N0 > 0 must be used).
C			            otherwise : no BRDF or BPDF simulation.
C
C            IFRESNEL (I4) : (O) Index specifying whether the simulations
C      			     include the Fresnel reflection on a flat sea:
C			      1 : for a computation of the Fresnel reflection.
C			      Otherwise : not taken into account.
C
C            IGLI (I4) : (O) Index specifying whether the simulations include the Glitter
C			        1 : for a Glitter computation
C			        Otherwise : not taken into account.
C 
C            IROUJEAN (I4) : (O) Index specifying whether the simulations
C      			     include the Roujean reflection model:
C			            1 : for a Roujean BRDF computation (N0 > 0 must be used).
C			            Otherwise : not taken into account.
C
C            IRONDEAUX (I4) : (O) Index specifying whether the simulations
C      			      include the Rondeaux - Herman reflection model:
C			             1 : for a Rondeaux - Herman BRDF computation (N0 > 0 must be used).
C			             Otherwise : not taken into account.
C
C            IBREON (I4) : (O) Index specifying whether the simulations
C      			   include the Breon reflection model:
C			          1 : for a Breon BRDF computation (N0 > 0 must be used).
C			          Otherwise : not taken into account.
C
C            IMAIGNAN (I4) : (O) Index specifying whether the simulations
C      			     include the Maignan reflection model:
C			            1 : for a Maignan BRDF computation (N0 > 0 must be used).
C			            Otherwise : not taken into account.
C
C            INADAL (I4) : (O) Index specifying whether the simulations
C      			   include the Nadal reflection model:
C			          1 : for a Nadal BRDF computation (N0 > 0 must be used).
C			          Otherwise : not taken into account.
C
C            ALPHA (DOUBLE) : (O) 1D array of CTE_OS_NB_MAX + 1 elements
C				      containing the Alpha coefficients of the decomposition in Legendre
C				      functions of the elements Pij of the phase function
C
C            BETA (DOUBLE) : (O) 1D array of CTE_OS_NB_MAX + 1 elements
C				     containing the Beta coefficients of the decomposition in Legendre
C				     functions of the elements Pij of the phase function
C
C            GAMMA (DOUBLE) : (O) 1D array of CTE_OS_NB_MAX + 1 elements
C				      containing the Gamma coefficients of the decomposition in Legendre
C				      functions of the elements Pij of the phase function
C
C            ZETA (DOUBLE) : (O) 1D array of CTE_OS_NB_MAX + 1 elements
C				     containing the Zeta coefficients of the decomposition in Legendre
C				     functions of the elements Pij of the phase function
C
C            IND_ANGOUT_FIN (DOUBLE) : (O) Indices of the angle arrays for which there must be
C                                      an output in the result file of the user angles
C      				    (1: output required)
C                                      Indices starting from 0 for Python.
C				 				
C            IER (I4) : (O) error code =0 if no error, =-1 otherwise
C 
C
C    -->Files
C       The file of angles and internal parameters for the expansion limits.
C       The surface reflection matrix file
C       The Aerosols file
C
C
C Results provided by the program
C ----------------------------------
C   Content of the angle file
C   Content of the aerosol file
C   Type of modelling of the surface reflection
C   Name of the surface file
C   Trace file (optional)
C              (logical unit number of the trace file =IDLOG)
C
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
C
C    -->Content of the angle file:
C       Formatted ascii file, created by SOS_ANGLES.
C       It contains the information on the angles to be used
C       and on the adapted expansion orders:
C		NB_TOTAL_ANGLES :    total nb of angles to be used
C             NB_GAUSS_ANGLES :    nb of Gauss angles
C             ANGLES_USERFILE :    Name of the user angles file
C	       SOLAR ZENITH ANGLE : Solar zenith angle (degrees)
C		INTERNAL_IMUS :      Corresponding array index number
C             INTERNAL_OS_NB :   Order of the expansions of the phase functions
C				     in Legendre polynomials to be used.
C             INTERNAL_OS_NS :   Maximum order to be used for the expansion in
C                                Legendre polynomials of the Fresnel matrix elements
C                                and in Fourier series of the radiance computations.
C                                (not used in this program)
C             INTERNAL_OS_NM :   Maximum order of the expansion in Fourier series
C                                of the G function of the reflection matrix computations.
C                                (not used in this program)
C             INDEX COS_ANGLE WEIGHT OUTPUT: lines giving the angle number, its cosine
C                                            and its weight, and whether it is a "user" output
C                                            with format I4,1X,2D21.14,1X,I4
C         
C   -->Content of the MIE file
C       Unformatted binary file, created by SOS_MIE
C       1 header record with
C                Refractive index (real part)
C                Refractive index (imaginary part)
C                Limit size parameter
C                Number of positive Gauss angles
C
C      1 record for each value of the size parameter with:
C        - Size parameter (real),
C        - Extinction efficiency factor (real),
C        - Scattering efficiency factor (real),
C        - Asymmetry factor (double precision),
C        - Array of the phase function Mie(mu) (real),
C        - Array of the polarised phase function Qmie(mu) (real),
C        - Array of the polarised phase function Umie(mu) (real).
C
C    --> The arrays have a useful size limited to MIE_NBMU:
C        they are created and read with this limitation.
C
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
C   - SOS_NOM_FICGLITTER
C   - SOS_NOM_FICROUJ
C   - SOS_NOM_FICRH
C   - SOS_NOM_FICBREON
C   - SOS_NOM_FICNADAL
C   - SOS_NOM_FICMAIGNAN
C
C============================================================================
       SUBROUTINE SOS_PREPA_OS(FICANGLES, FICGRANU, 
     &	       ISURF, IND_SURF,WIND, 
     &		K0_ROUJEAN, K1_ROUJEAN, K2_ROUJEAN, 
     &		ALPHA_NADAL, BETA_NADAL, COEF_C_MAIGNAN,
     &		REPSURF,FICSURF,TRACE,IDLOG,
     &		LUM_NBMU,N0,PIZ,PIZTR,A,RMU,GA,OS_NB,IMAT_SURF,
     &		IGLI,IFRESNEL,IROUJEAN,IRONDEAUX,IBREON,INADAL,IMAIGNAN,
     &		ALPHA,BETA,GAMMA,ZETA,IND_ANGOUT_FIN,IER)      
      
      IMPLICIT NONE

C* Definition of variables
C*-----------------------------------------------------------------
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
                                   !Cosines of the Gauss (and user) angles.
      DOUBLE PRECISION GA(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
                                   !Weights of the Gauss (and user) angles.
 
 
      DOUBLE PRECISION TETAS	! Solar zenith angle in degrees (0 < tetas < 90).
      				       ! (complementary angle: mus= -cos(tetas))
      
      DOUBLE PRECISION IND_SURF   ! Refractive index surface / air.
      
      DOUBLE PRECISION WIND       ! Wind speed in m/s.
      
      DOUBLE PRECISION K0_ROUJEAN ! Parameter k0 of the Roujean BRDF model.
      DOUBLE PRECISION K1_ROUJEAN ! Parameter k1 of the Roujean BRDF model.
      DOUBLE PRECISION K2_ROUJEAN ! Parameter k2 of the Roujean BRDF model.
      
      DOUBLE PRECISION ALPHA_NADAL ! Parameter alpha of the Nadal BPDF model.
      DOUBLE PRECISION BETA_NADAL  ! Parameter beta of the Nadal BPDF model.

      DOUBLE PRECISION COEF_C_MAIGNAN !Parameter C*exp(-NDVI) of the Maignan model.
            
      DOUBLE PRECISION A	! Truncation coefficient.
      DOUBLE PRECISION PIZ	! Single scattering albedo.
      DOUBLE PRECISION PIZTR	! Single scattering albedo after truncation.

      DOUBLE PRECISION ALPHA(0:CTE_OS_NB_MAX)	! Parameters alpha(k).
      DOUBLE PRECISION BETA(0:CTE_OS_NB_MAX)	! Parameters beta(k).
      DOUBLE PRECISION GAMMA(0:CTE_OS_NB_MAX)	! Parameters gamma(k).
      DOUBLE PRECISION ZETA(0:CTE_OS_NB_MAX)	! Parameters zeta(k).
								               
     
      INTEGER*4 N0	! Number of the Gauss mu for the solar
      			! zenith angle.
			! (complementary angle: mus = -rmu(N0))
			! (N0 > 0 => to take it into account).

      INTEGER*2 ISURF	! Surface type:
			       ! 0 :Lambertian surface of albedo RHO.
			       ! 1 :Lambertian surface + glitter.
			       ! 2 :Lambertian surface + Fresnel reflection.
			       ! 3 :Lambertian surface + Roujean BRDF.
			       ! 4 :Lambertian surface + Roujean BRDF + Rondeaux BPDF.
			       ! 5 :Lambertian surface + Roujean BRDF + Breon BPDF.
			       ! 6 :Lambertian surface + Roujean BRDF + Nadal BPDF.
				
      INTEGER*4 IGLI	! Index specifying whether the simulations
      			! include the Glitter:
			!    1 => for Glitter computation (N0 > 0 must then be used).
			!    Otherwise => no Glitter simulation.
	
      INTEGER*4 IFRESNEL ! Index specifying whether the simulations
      			    ! include the Fresnel reflection on a flat sea:
			    !    1 => for a computation of the Fresnel reflection.
			    !    Otherwise => not taken into account.
	
      INTEGER*4 IROUJEAN ! Index specifying whether the simulations
      			    ! use the Roujean reflection model:
			    !    1 => for Roujean BRDF computation (N0 > 0 must be used).
			    !    Otherwise => no Roujean BRDF simulation.
			 
      INTEGER*4 IRONDEAUX ! Index specifying whether the simulations
      			     ! use the Rondeaux - Herman reflection model:
			     !    1 => for RH BPDF computation (N0 > 0 must be used).
			     !    Otherwise => no RH BPDF simulation.
			
      INTEGER*4 IBREON   ! Index specifying whether the simulations
      			    ! use the Breon reflection model:
			    !    1 => for Breon BPDF computation (N0 > 0 must be used).
			    !    Otherwise => no Breon BPDF simulation.

      INTEGER*4 INADAL   ! Index specifying whether the simulations
      		           ! use the Nadal reflection model:
		           !    1 => for Nadal BPDF computation (N0 > 0 must be used).
			    !    Otherwise => no Nadal BPDF simulation.

      INTEGER*4 IMAIGNAN ! Index specifying whether the simulations
      		           ! use the Maignan reflection model:
		           !    1 => for Maignan BPDF computation (N0 > 0 must be used).
			    !    Otherwise => no Maignan BPDF simulation.
			 
	
      INTEGER*4 IMAT_SURF ! Index specifying whether the simulations
      			     ! use a BRDF or BPDF reflection model:
			     !    1 => for BRDF or BPDF computation (N0 > 0 must be used).
			     !    Otherwise => no BRDF or BPDF simulation.
			  			  			 											

      INTEGER*4 IER     ! error code =0 if no error, =-1 otherwise

      INTEGER*4 K       ! Array index
      INTEGER*4 I,J     ! Array index
      
      INTEGER*4 LDIRSURF ! length of the character string REPSURF.
      INTEGER*4 LFICSURF ! length of the character string FICSURF.
      
      INTEGER*4 NBMU_GAUSS   ! Number of Gauss angles used
      INTEGER*4 LUM_NBMU     ! Number of actual angles of the simulations
      INTEGER*4 OS_NB 	 ! Order of the expansions of the phase functions
			        ! in Legendre polynomials.

      INTEGER*4 IND_ANGOUT(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
                            ! Indices of the angle arrays for which there must be
      				! an output in the result file of the user angles
				! (1: output required)
			
      INTEGER*4 IND_ANGOUT_FIN(0:CTE_OS_NBMU_MAX)  
                            ! Indices of the angle arrays for which there must be an output in the result file of the user angles
                            ! ==> 1: output required

      INTEGER*4 IFIN        ! Nb of significant characters of a string.

			               
      
      LOGICAL TRACE   ! = true if writing in the trace file
      LOGICAL EX	! Checks the existence of a file

      CHARACTER*CTE_LENDIR REPSURF	! Storage directory of the surface reflection
      				       ! matrices (Glitter file or other BRDF/BPDF).
      CHARACTER*CTE_LENFIC2 FICSURF ! File containing the surface reflection
      				        ! matrix (Glitter file or other BRDF/BPDF).
      CHARACTER*CTE_LENFIC2 FICANGLES           ! path and name of the file of the actual angles
      CHARACTER*CTE_LENFIC2 FICANGLES_USER      ! path and name of the "user" angles file
      
      CHARACTER*CTE_LENFIC2 FICGRANU	! name of the AEROSOLS input file

      CHARACTER*CTE_LENFIC1 FICS ! name of the SURFACE file rebuilt from
      				     ! the BRDF / BPDF simulation arguments.

      CHARACTER*LENLINE LIGNE_TEXTE    ! Character line for decoding parameters
      CHARACTER*1000 CVAL	           ! Character string to read unused values
             
      INTEGER*4 IDLOG	     ! Identifier number of the trace file
     	  
      	        		 

C* Initialisation
C-----------------
      IER=0
      IGLI=0
      IFRESNEL=0
      IROUJEAN=0
      IRONDEAUX=0
      IBREON=0
      INADAL=0
      IMAIGNAN=0
      IMAT_SURF=0
      IF (ISURF.EQ.1) IGLI=1
      IF (ISURF.EQ.2) IFRESNEL=1
      IF (ISURF.GE.3) IROUJEAN=1
      IF (ISURF.EQ.4) IRONDEAUX=1
      IF (ISURF.EQ.5) IBREON=1
      IF (ISURF.EQ.6) INADAL=1
      IF (ISURF.EQ.7) IMAIGNAN=1
      IF ((IGLI.EQ.1).OR.(IROUJEAN.EQ.1).OR.(IRONDEAUX.EQ.1).
     &     OR.(IBREON.EQ.1).OR.(INADAL.EQ.1).
     &     OR.(IMAIGNAN.EQ.1)) IMAT_SURF=1


C* Reading of the angle and weight tables used by the simulations
C* -----------------------------------------------------------------     
      OPEN(UNIT=10,FILE=FICANGLES,STATUS='OLD',ERR=950)

      READ(10,'(a)',ERR=951) LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) LUM_NBMU
      
      READ(10,'(a)',ERR=951) LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) NBMU_GAUSS
      
      READ(10,'(a)',ERR=951) LIGNE_TEXTE
      !Only the part after the separator ':' is kept
      IFIN=INDEX(LIGNE_TEXTE,':')	
      LIGNE_TEXTE=LIGNE_TEXTE(IFIN+1:LENLINE)	
      !Search for the end of the useful line ' ' to reduce the string size
      IFIN=INDEX(LIGNE_TEXTE,' ')
      FICANGLES_USER=LIGNE_TEXTE(1:IFIN)

      READ(10,'(a)',ERR=951) LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) TETAS

      READ(10,'(a)',ERR=951) LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) N0

      READ(10,'(a)',ERR=951) LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) OS_NB

      READ(10,*,ERR=951) CVAL
      READ(10,*,ERR=951) CVAL
      READ(10,*,ERR=951) CVAL

      DO J=1,LUM_NBMU
         READ(10,*,ERR=951) CVAL, RMU(J), GA(J), IND_ANGOUT(J)
      ENDDO 

      
      DO J=1,LUM_NBMU
         RMU(-J) = -RMU(J)
	  GA(-J) =  GA(J)
	  IND_ANGOUT(-J) = IND_ANGOUT(J)
	  IND_ANGOUT_FIN(J-1) = IND_ANGOUT(J)	!Indexing from 0 for Python
      ENDDO 
      
      !Applicability check
      IF (ABS(N0).GT.LUM_NBMU) GOTO 1000
      IF ((TETAS.LT.0.D-00).OR.(TETAS.GE.90.D-00)) GOTO 1001
      IF ((ISURF.NE.0).AND.(N0.LT.0)) GOTO 1012
       
      IF (TRACE) THEN
        WRITE (IDLOG,*,ERR=960) "USED ANGLES"
	 WRITE (IDLOG,*,ERR=960) "***********"
	 WRITE (IDLOG,*,ERR=960) " "
	 WRITE (IDLOG,*,ERR=960) "Total number of angles : ",LUM_NBMU
	 WRITE (IDLOG,*,ERR=960) "Number of Gauss angles : ",NBMU_GAUSS
	 WRITE (IDLOG,*,ERR=960) "Solar incidence : "
	 WRITE (IDLOG,*,ERR=960) " --> Solar zenith angle: ", TETAS
	 WRITE (IDLOG,*,ERR=960) " --> Index in the angles table: ", N0 
	 WRITE (IDLOG,*,ERR=960) " "
        IF (FICANGLES_USER.NE.'NO_USER_ANGLES') THEN
	     WRITE(IDLOG,*,ERR=960)  
     &	       "A file of user's angles is used :", FICANGLES_USER
     	     WRITE(IDLOG,*,ERR=960) 
     &       "Number of user's angles + solar angle SZA :",
     &         (LUM_NBMU-NBMU_GAUSS)
	 ELSE
	     WRITE(IDLOG,*,ERR=960)  
     &	       "No file of user's angles"
	 ENDIF !End of test: no user angle file
        WRITE (IDLOG,*,ERR=960) " " 
	 WRITE (IDLOG,*,ERR=960) "Max order NB ", OS_NB
	 WRITE (IDLOG,*,ERR=960) " " 
	 WRITE (IDLOG,*,err=960) "List of used angles :"
	 WRITE (IDLOG,*,err=960)  
     &  "Index        Mu             angle (deg)    Weigth"// 
     &  "         user's angle"
	 DO 1200 J=-LUM_NBMU,LUM_NBMU
	    IF (J.EQ.0) GOTO 1200
	    IF (J.LT.0) THEN 
	        WRITE(IDLOG,223,ERR=960) J, RMU(J),
     &	                     -DACOS(-RMU(J))*180./INCTE_PI, GA(J),
     &                       IND_ANGOUT(-J)
            ELSE
	        WRITE(IDLOG,223,ERR=960) J, RMU(J),
     &	                      DACOS(RMU(J))*180./INCTE_PI, GA(J), 
     &      		        IND_ANGOUT(J) 
	    ENDIF
1200    CONTINUE	
	 WRITE (IDLOG,*,ERR=960) " " 
	 WRITE (IDLOG,*,ERR=960) " "  	
      ENDIF !Trace file

	     	
      
C* Determination of the surface BRDF or BPDF file
C  ---------------------------------------------------

      IF (IMAT_SURF.EQ.1) THEN

C* Retrieval of the name of the surface reflection file
C  ----------------------------------------------------------	
         IF (FICSURF.EQ.'DEFAULT') THEN

	     IF (ISURF.EQ.1) THEN
	         CALL SOS_NOM_FICGLITTER(NBMU_GAUSS,FICANGLES_USER,
     &		                         IND_SURF,WIND,FICS)	        
	     ENDIF
	     
	     IF (ISURF.EQ.3) THEN
	         CALL SOS_NOM_FICROUJ(NBMU_GAUSS,FICANGLES_USER,
     &                                K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &                                FICS)    	        
	     ENDIF
	     
	     IF (ISURF.EQ.4) THEN
	         CALL SOS_NOM_FICRH(NBMU_GAUSS,FICANGLES_USER,IND_SURF,
     &	                            K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &                             FICS)
	     ENDIF
	     
	     IF (ISURF.EQ.5) THEN
	         CALL SOS_NOM_FICBREON(NBMU_GAUSS,FICANGLES_USER,
     &	                               IND_SURF,
     &                                K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &                                FICS)
	     ENDIF
	     
	     IF (ISURF.EQ.6) THEN
	         CALL SOS_NOM_FICNADAL(NBMU_GAUSS,FICANGLES_USER,
     &                                IND_SURF,ALPHA_NADAL,BETA_NADAL,
     &	                               K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &                                FICS)
	     ENDIF

	     IF (ISURF.EQ.7) THEN
	         CALL SOS_NOM_FICMAIGNAN(NBMU_GAUSS,FICANGLES_USER,
     &                                  IND_SURF,COEF_C_MAIGNAN,
     &                                  K0_ROUJEAN,K1_ROUJEAN,
     &	                                 K2_ROUJEAN,FICS) 	
	     ENDIF

	     
	     LDIRSURF = INDEX(REPSURF,' ')
            LDIRSURF = LDIRSURF - 1
            IF(LDIRSURF.LE.0) LDIRSURF = CTE_LENDIR
            FICSURF = REPSURF(1:LDIRSURF)//'/'//FICS
	 
	 ENDIF        
		
C* Check of the existence of the reflection matrix file
C  ----------------------------------------------------------	
	 
        LFICSURF=INDEX(FICSURF,' ')-1 
        INQUIRE(FILE=FICSURF(1:LFICSURF),EXIST=EX)
        IF (EX.EQV..FALSE.) GOTO 1020   !Error case
	       
      ENDIF ! End of the test for the use of the reflection matrix file


      
C* Reading of the AEROSOLS file
C  ---------------------------
      OPEN(UNIT=1,FILE=FICGRANU,STATUS='OLD',ERR=913)
	 
      ! Extinction cross section
      READ(1,'(a)',ERR=923) LIGNE_TEXTE	

      ! Scattering cross section
      READ(1,'(a)',ERR=923) LIGNE_TEXTE

      ! Asymmetry factor (not truncated)
      READ(1,'(a)',ERR=923) LIGNE_TEXTE
          
      ! Truncation coefficient
      READ(1,'(a)',ERR=923) LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) A
	 
      ! Single scattering albedo (possibly truncated)
      READ(1,'(a)',ERR=923) LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*)  PIZTR

      !Reading of the information on the column contents
      READ(1,'(a)',ERR=923) LIGNE_TEXTE
      READ(1,'(a)',ERR=923) LIGNE_TEXTE
      READ(1,'(a)',ERR=923) LIGNE_TEXTE
      
      DO K=0,OS_NB
         READ(1,*,ERR=923) ALPHA(K),BETA(K),GAMMA(K),ZETA(K)
      ENDDO      


C* Computation of the non-truncated single scattering albedo
C* from the truncated value
C  ---------------------------------------------------
      PIZ=PIZTR/(1+0.5*A*(PIZTR-1))
            
     
      
      
C* File closing
C-------------------
      CLOSE(1)
      CLOSE(10)

      
C* Format 
C----------
  223 FORMAT(I4,1X,D21.14,2X,F7.3,2X,D21.14,2X,I2)  


C* Nominal end
C-------------------
      GOTO 9999
      
C* Error case and return of status 1 to the shell
C----------------------------------------------    

  950 WRITE(6,*) '  SOS_PREPA_OS : ERROR_950'
      WRITE(6,*) '  ERROR on angles file opening for SOS'
      GOTO 9998

  951 WRITE(6,*) '  SOS_PREPA_OS : ERROR_951'
      WRITE(6,*) '  ERROR on angles file reading for SOS'
      GOTO 9998

  960 WRITE(6,*) '  SOS_PREPA_OS : ERROR_960'      
      WRITE(6,*) '  ERROR on logfile writing for SOS'
      GOTO 9998
             
  913 WRITE(6,*) '  SOS_PREPA_OS : ERROR_913'      
      WRITE(6,*) '  ERROR on AEROSOLS file opening for SOS'
      GOTO 9998    

  922 WRITE(6,*) '  SOS_PREPA_OS : ERROR_922'      
      WRITE(6,*) '  ERROR on PROFILE file reading for SOS'
      GOTO 9998
 
  923 WRITE(6,*) '  SOS_PREPA_OS : ERROR_923'        
      WRITE(6,*) '  ERROR on AEROSOLS file reading for SOS'
      GOTO 9998 

 1000 WRITE(6,*) '  SOS_PREPA_OS : ERROR_1000'                 
      WRITE(6,*) '  ERROR on angles file value for IMUS' 
      WRITE(6,*) '  Out of available range.'
      WRITE(6,*) '  --> Defined value : ',N0
      WRITE(6,*) '      Have to be define between 1 and ',LUM_NBMU
      GOTO 9998

 1001 WRITE(6,*) '  SOS_PREPA_OS : ERROR_1001'                 
      WRITE(6,*) '  ERROR on TETAS value : out of available range'
      GOTO 9998      

 1012 WRITE(6,*) '  SOS_PREPA_OS : ERROR_1012'                 
      WRITE(6,*) '  ERROR on angles file value for IMUS' 
      WRITE(6,*) '  Out of available range.'
      WRITE(6,*) '  --> ISURF > 0 requires IMUS > 0 '
      WRITE(6,*) '      IMUS corresponding to the number '
      WRITE(6,*) '      of the solar angle in the angles file.'
      GOTO 9998     
                                                                                                       
 1020 WRITE(6,*) '  SOS_PREPA_OS : ERROR_1020'                 
      WRITE(6,*) '  ERROR on FICSURF value : file does not exist.'
      WRITE(6,*) '  ISURF: ',ISURF
      WRITE(6,*) '  IMAT_SURF: ',IMAT_SURF     
      GOTO 9998 
        
 9998 IER=-1
    
 9999 RETURN
 
      END	!END OF PROGRAM SOS




