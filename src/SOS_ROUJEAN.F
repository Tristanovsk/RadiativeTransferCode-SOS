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
C* FILE: SOS_ROUJEAN.F
C* PROJECT: Successive orders of scattering
C* ROLE: Generates the Roujean BRDF files.
C*
C* AUTHOR:
C*    CS GROUP France, B.Lafrance.
C* DATE: 30/04/2002
C*
C* MOD:VERSION:1.0: Creation of the program
C* MOD:VERSION:1.1: Improvement of the comments
C*
C* MOD:VERSION:2.0: 10/01/2010 
C*
C*   - Removal of the constant CTE_OS_FICANGLE (from SOS.h) and of the reading of the corresponding
C*     Gauss angle file, in favour of passing the angle (RMU) and weight (CHR) tables as arguments.
C*
C*   - Change of the size of the radiance arrays:
C*     ==> Changed from CTE_OS_NBMU to CTE_OS_NBMU_MAX.
C*     ==> Loop on the angles up to LUM_NBMU (new input argument)
C*         equal to the actual number of angles used.
C*
C*   - Removal of the constants CTE_OS_NB (from SOS.h)
C*     in favour of passing the value OS_NB as argument.
C*
C*   - Change of the size of the arrays of the decomposition in Legendre polynomials:
C*     ==> Changed from CTE_OS_NB to CTE_OS_NB_MAX.
C*     ==> Loop on the orders up to OS_NB (new input argument)
C*
C*   - Routine SOS_FSF_ROUJEAN: change of the routine call to include the argument LUM_NBMU
C*     (nb of useful angles) and to limit the loops to the useful bounds LUM_NBMU,
C*     and to include the argument OS_NB.
C*
C*   - Routine SOS_MISE_FORMAT_RJ: change of the record line of the
C*     arrays of the reflection matrix at order IS of the Fourier series
C*     decomposition:
C*     ==> limited to the useful bounds -LUM_NBMU:LUM_NBMU
C*       Rem: no change of the record bounds of the temporary file
C*             RES_MAT_REFLEX since it is not a final product.
C*     Use of OS_NB.
C*
C* MOD:VERSION:2.1: 06/06/2013
C*   -  Change of the value of Pi, declared as a constant in the file
C*      (INCTE_PI replaces CTE_PI)
C*
C* MOD:VERSION:2.2: 24/02/2016 
C*    - Minor adjustments to strictly respect the 72 columns
C*
C* MOD:VERSION:2.3: 18/05/2020
C*
C*   - Added test in routine SOS_CALC_F_ROUJEAN
C*     for the computation of the cosine of the angle Zeta, to avoid
C*     numerical rounding giving a value 1+Epsilon or -1-Epsilon.
C*
C*   - Introduction of a new threshold constant: CTE_SEUIL_NUM
C*     (defined in inc/SOS.h)
C* 
C*   - Change of the stored values of the BRDF function F
C*     (by routine SOS_CALC_F_ROUJEAN):
C*     ==> F * cos(Reflection angle) * cos(Incidence angle) is stored
C*         without any threshold on the angles for the computation of the weighting
C*         cosines. Indeed, the computations of the incident flux at the surface
C*         in the core of the OS code (SOS_OS) need the true
C*         angles.
C*
C* MOD:VERSION:2.4: 31/08/2022
C*    - Constants renamed from SOS_* to CTE_*
C*    - Traces converted to English.
C*
C******************************************************************************

C----------------------------------------------------------------------------
C Definition of constants
C---------------------------------------------------------------------------- 
C Constants used:
C    CTE_LENDIR : Length of directory names.
C     CTE_LENFIC2 : Length of file names with directory tree.
C    CTE_OS_NBMU_MAX : Maximum number of positive angles.
C    INCTE_PI : Value of PI
C    CTE_OS_NB_MAX :  Maximum value of the limit order of the expansion in
C                     Legendre polynomials.
C    CTE_PH_NU : Number of azimuth angles (2**NQ)
C    CTE_TETAS_LIM_ROUJEAN : Limit solar zenith angle for the BRDF
C                            computations of the Roujean model.
C    CTE_TETAV_LIM_ROUJEAN : Limit viewing zenith angle for the BRDF
C                            computations of the Roujean model.
C    CTE_SEUIL_SF_ROUJEAN : Threshold to stop the Fourier series
C                           decomposition of the Roujean BRDF.
C    CTE_SEUIL_NUM : Threshold of numerical precision
C
C----------------------------------------------------------------------------
#include "SOS.h"
#define INCTE_PI DACOS(-1.D+00)

C==============================================================================
C PROCEDURE: SOS_ROUJEAN
C ==========      
C      This procedure generates the BRDF files of the Roujean model.
C
C
C Input data of the program
C ------------------------------
C    -->Parameters/arguments of the program
C       LUM_NBMU (I4) (I) : Number of actually useful (positive) angles
C	 RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (DOUBLE) (I) : Cosines of the angles.
C	 CHR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (DOUBLE) (I) : Weights of the angles.
C       OS_NB (I4) : (I) Order of the expansions of the phase functions
C		           in Legendre polynomials to be used.
C 	 K0  (DOUBLE)  (I) : Coefficient of the Roujean model.
C 	 K1  (DOUBLE)  (I) : Coefficient of the Roujean model.
C 	 K2  (DOUBLE)  (I) : Coefficient of the Roujean model.
C       FIC_RES_MAT_REFLEX (CHARACTER*CTE_LENFIC2)  :  
C          (I) Name of the temporary file RES_MAT_REFLEX with complete directory tree.
C       FICROUJ  (CHARACTER*CTE_LENFIC2) (I) : complete name of the Roujean file
C                                              (complete directory tree)
C       TRACE  (LOGICAL) (I) : true if writing in the trace file.
C       IER  (I4)   (O) : error code =0 if no error, =-1 otherwise
C
C    -->File
C       The Gauss angle file
C
C Result provided by the program
C --------------------------------
C   Roujean BRDF file
C   Trace file (optional),
C               if a file with the same name already exists, it is overwritten
C              (logical unit number of the trace file =99)
C
C 
C Environment variable
C ------------------------
C   None
C
C
C Description of the files used
C ---------------------------------
C
C   -->Content of the Roujean BRDF file
C       Unformatted binary file, created by SOS_MISE_FORMAT_RJ
C       It contains the surface reflection matrix:
C       i.e. the elements Pij of the reflection matrix
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
C        --> Record on the domain limited to the useful angles
C                                  --------------------------------
C
C   -->Content of the temporary file RES_MAT_REFLEX
C       Unformatted binary file,
C       File created in the directory where the program is run
C       If a file with the same name already exists, it is overwritten
C       Record for each pair (Teta1, Teta2) of the
C       Fourier series decomposition of the Roujean BRDF:
C       M11 (n real values)
C
C    The temporary file is overwritten if it exists, then
C    deleted at the end of the program
C
C
C Common used:
C --------------
C None
C
C
C Error cases:
C ------------
C     -Opening of a file
C      Error if the Roujean BRDF file already exists
C     -Reading/writing in a file
C
C     Message displayed on screen, program stop and
C     return of status 1 to the shell
C        
C   
C
C  Subroutines used:
C --------------------------
C  This procedure calls the routines:
C      - SOS_FSF_ROUJEAN
C      - SOS_MISE_FORMAT_RJ
C
C
C=============================================================================

      SUBROUTINE SOS_ROUJEAN(LUM_NBMU,RMU,CHR,OS_NB,K0,K1,K2,
     &                       FIC_RES_MAT_REFLEX,FICROUJ,TRACE,IER)

      IMPLICIT NONE


C* Definition of variables
C*--------------------------------------------------------------------------
      CHARACTER*CTE_LENFIC2 FICROUJ  ! (I) Complete name of the result ROUJEAN file.

      CHARACTER*CTE_LENFIC2 FIC_RES_MAT_REFLEX ! Name of the file RES_MAT_REFLEX
                                               ! with complete directory tree

      LOGICAL TRACE         ! (I) true if writing in the trace file

      DOUBLE PRECISION K0	! (I) Parameter of the Roujean model.
      DOUBLE PRECISION K1	! (I) Parameter of the Roujean model.
      DOUBLE PRECISION K2	! (I) Parameter of the Roujean model.
      	
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   ! Array of the Gauss mu
      DOUBLE PRECISION CHR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   ! Gauss weights.
      
      
      INTEGER*4 LUM_NBMU ! Number of actual angles of the simulations
      INTEGER*4 OS_NB    ! Order of the expansions of the phase functions
			    ! in Legendre polynomials.
      INTEGER*4 I	   ! Array index.
      INTEGER*4 IFIN    ! nb of significant characters of a string
      INTEGER*4 IER     ! (O) error code =0 if no error, =-1 otherwise

C* Initialisation of the error code
C------------------------------------------------------
      IER=0

C* Opening of the trace file
C------------------------------------------------------
      IF (TRACE) THEN
         WRITE(99,*,err=921)
	  WRITE(99,*,err=921) '*****************************************'
         WRITE(99,*,err=921) 'Calculation of the ROUJEAN file :' 
	  WRITE(99,*,err=921) FICROUJ 
         WRITE(99,*,err=921) ' '
         WRITE(99,*,err=921)'Simulation parameters : '
         WRITE(99,*,err=921)' k0 =  ',K0
	  WRITE(99,*,err=921)' k1 =  ',K1
	  WRITE(99,*,err=921)' k2 =  ',K2
	  WRITE(99,*,err=921)
	  WRITE(99,*,err=921)' Limit of viewing zenith angle = ',
     &	                     CTE_TETAV_LIM_ROUJEAN
	  WRITE(99,*,err=921)' Limit of solar zenith angle = ',
     &	                     CTE_TETAS_LIM_ROUJEAN
         WRITE(99,*,err=921)
	  WRITE(99,*,err=921)	 	 
      END IF
        


C* Fourier series decomposition of the function F(Teta1,Teta2,phi)
C* and recording of the BRDF expanded in Fourier series.
C* (BRDF converted into normalised radiance and adapted to the OS code)
C* => Creation of the temporary file RES_MAT_REFLEX.
C---------------------------------------------------------------------
      CALL SOS_FSF_ROUJEAN(LUM_NBMU,RMU,OS_NB,K0,K1,K2,
     &                     FIC_RES_MAT_REFLEX,TRACE,IER)
      IF (IER.NE.0) GOTO 995   ! error case
      
      
C* Rearrangement of the results to be compatible with the OS code.
C-----------------------------------------------------------------------
      CALL SOS_MISE_FORMAT_RJ(FIC_RES_MAT_REFLEX,FICROUJ,
     &                        LUM_NBMU,OS_NB,IER)
      IF (IER.NE.0) GOTO 996   ! error case

C* Nominal end
C-------------------
      GOTO 9999
      
C* Error case and return of status 1 to the shell
C----------------------------------------------
  912 WRITE(6,*) 'Erreur a l''ouverture du fichier angles de Gauss'
      IER=-1
      GOTO 9999
  921 WRITE(6,*) 'Erreur d''écriture sur le fichier trace Surface'
      IER=-1
      GOTO 9999
  922 WRITE(6,*) 'Erreur de lecture sur le fichier angles de Gauss'
      IER=-1
      GOTO 9999      
  995 WRITE(6,*) 'Erreur dans la routine SOS_FSF_ROUJEAN'
      IER=-1
      GOTO 9999
  996 WRITE(6,*) 'Erreur dans la routine SOS_MISE_FORMAT_RJ'
      IER=-1
      GOTO 9999  
      
C* Format
C---------   
 2    FORMAT(2(D21.14))    
              
 9999 CONTINUE

C* Deletion of the temporary files
C---------------------------------------   
      OPEN(UNIT=3,FILE=FIC_RES_MAT_REFLEX,FORM='UNFORMATTED')
      CLOSE(3,STATUS='DELETE') 
                                                      
      RETURN
      END    !END OF PROCEDURE SOS_ROUJEAN





C==============================================================================
C PROCEDURE: SOS_FSF_ROUJEAN
C ==========
C      This procedure performs the Fourier series decomposition
C      (on phi) of the function F(Teta1,Teta2,phi):
C
C	 F(Teta1,Teta2,phi) = K0 + K1*F1(Teta1,Teta2,phi) + K2*F2(Teta1,Teta2,phi)
C
C	 with Teta1 : incidence zenith angle: tetas  (0 < Teta1 < PI/2),
C             Teta2 : reflection zenith angle: tetav (0 < Teta2 < PI/2),
C             phi   : relative azimuth angle between the directions
C                     of incidence and reflection with the convention
C                     of the Roujean model (0 <= phi <= PI),
C
C	      K0, K1, K2 : coefficients of the Roujean model,
C			
C             F1 = [ (PI-phi)*cos(phi)+sin(phi) ]*tgs*tgv / 2*PI
C                 -[ tgs + tgv + SQRT{tgs**2+tgv**2-2*tgs*tgv*cos(phi)} ] / PI
C
C		  with tgs = tan(tetas) and tgv = tan(tetav),
C
C            F2 = [ (PI/2-zeta)*cos(zeta)+sin(zeta) ] * (4/3*PI) 
C                 / [ cos(tetas)+cos(tetav) ]   - 1/3
C
C                 with zeta the phase angle defined by:
C                 cos(zeta) = cos(tetas)*cos(tetav)+sin(tetas)*sin(tetav)*cos(phi). 
C
C
C      The BRDF of the Roujean model is converted from reflectance into a
C      quantity useful for the OS: normalised radiance * cos(Teta2),
C      i.e.: BRDF * cos(Teta1) * cos(Teta2).
C
C
C      The procedure stores the results in the file "RES_MAT_REFLEX"
C      as:
C         Decomposed function F for (Teta1,Teta2) : (E(IS),IS=0,CTE_OS_NB)
C
C
C   REM : phiOS = PI - phiRoujean
C   *****************************
C
C Description of parameters
C --------------------------
C      LUM_NBMU (I4) (I) : Number of actually useful (positive) angles.
C      RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)	 (I) : Array of the Gauss mu.
C      OS_NB (I4) : (I) Order of the expansions of the phase functions
C		        in Legendre polynomials to be used.
C      K0	(double)  (I) : Coefficient of the Roujean model.
C      K1	(double)  (I) : Coefficient of the Roujean model.
C      K2	(double)  (I) : Coefficient of the Roujean model.
C      FIC_RES_MAT_REFLEX (CHARACTER)  : (O) Name of the temporary file RES_MAT_REFLEX
C                                            with complete directory tree.
C      TRACE    (logical) (I) : = true if writing in the trace file (unit 99)
C      IER  (I4)   (O) : error code =0 if no error, =-1 otherwise
C
C Other result provided by the procedure
C -------------------------------------
C      Temporary file RES_MAT_REFLEX
C
C 
C Description of the files used
C ---------------------------------
C
C   -->Content of the temporary file RES_MAT_REFLEX
C       Unformatted binary file,
C       File created in the directory where the program is run
C       If a file with the same name already exists, it is overwritten
C       Record for each pair (Teta1, Teta2) of the
C       Fourier series decomposition of the Roujean BRDF:
C       M11 (n real values)
C
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C     -Opening of a file
C     -Reading/writing in a file
C     -BRDF < 0
C
C     Message displayed on screen, the routine stops its computations and
C     returns status -1 to the calling program
C
C     
C  Routines used by the procedure:
C ------------------------------------
C  This program calls the routines:
C        - SOS_CALC_F_ROUJEAN
C
C==============================================================================
      SUBROUTINE SOS_FSF_ROUJEAN(LUM_NBMU,RMU,OS_NB,K0,K1,K2,
     &                           FIC_RES_MAT_REFLEX,TRACE,IER)
 	
      IMPLICIT NONE

C* Definition of variables
C*-------------------------------------------------------------------
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! Array of the Gauss mu

      DOUBLE PRECISION K0	! Coefficient of the Roujean model.
      DOUBLE PRECISION K1	! Coefficient of the Roujean model.
      DOUBLE PRECISION K2	! Coefficient of the Roujean model.

      DOUBLE PRECISION E(0:CTE_OS_NB_MAX) ! Result of the decomposition in
      				              ! Fourier series of the function F.
	
      DOUBLE PRECISION PHIOS	! Relative azimuth for the OS convention (radians).
      DOUBLE PRECISION PHIRJ	! Relative azimuth for the Roujean convention (radians).
      DOUBLE PRECISION C1	! Cos(Teta1).
      DOUBLE PRECISION S1	! Sin(Teta1).
      DOUBLE PRECISION C2	! Cos(Teta2).
      DOUBLE PRECISION S2	! Sin(Teta2).

      DOUBLE PRECISION F	! Value of the function F(phiRJ).

      DOUBLE PRECISION FMIN	! Minimum value kept for the function
      				! F(phiRJ) for fixed Theta1 and Theta2.
      DOUBLE PRECISION FMAX	! Maximum value of the function F(phiRJ).

      DOUBLE PRECISION U(0:CTE_PH_NU)  ! Function F(phiRJ) between 0 and PI
      			                  ! for fixed Theta1 and Theta2.

      DOUBLE PRECISION Q	! Step in phi.

      DOUBLE PRECISION T1	! Recombination of the function F(phi) from
      				! E(IS).

      DOUBLE PRECISION Y
      DOUBLE PRECISION B1	! Test value to stop the decomposition
      				! in Fourier series.
      DOUBLE PRECISION B1_PREC  ! Value of B1 for order IS-1.

      DOUBLE PRECISION PHIFMIN	!Relative azimuth giving Fmin (OS convention).
      DOUBLE PRECISION PHIFMAX	!Relative azimuth giving Fmax (OS convention).
     
	
      DOUBLE PRECISION FRECOMB	! Value of the function F recombined from its
      				       ! Fourier series decomposition.
      DOUBLE PRECISION ECART	! Relative difference between the exact value of F
      				       ! and the recombined one.
      DOUBLE PRECISION ECARTMAX ! Maximum relative difference between the exact value of F
      				    ! and the recombined one.
      
      DOUBLE PRECISION PHIECARTMAX ! Azimuth angle of the configuration giving ECARTMAX.
      DOUBLE PRECISION TVECARTMAX  ! Viewing zenith angle for ECARTMAX.
      DOUBLE PRECISION TSECARTMAX  ! Solar zenith angle for ECARTMAX.
      DOUBLE PRECISION F_ECARTMAX  ! Exact value of F giving ECARTMAX.
      DOUBLE PRECISION FRECOMB_ECARTMAX  ! Recombined value of F giving ECARTMAX.
      DOUBLE PRECISION SOMECART	! Sum of the relative differences.

      DOUBLE PRECISION ECARTMAX_C1C2 ! Maximum relative difference between the exact value of F
      				         ! and the recombined one
				         !=> for the current pair (Teta1, Teta2).
      DOUBLE PRECISION PHIECARTMAX_C1C2 ! Azimuth angle for ECARTMAX_C1C2.
      DOUBLE PRECISION F_ECARTMAX_C1C2  ! Exact value of F giving ECARTMAX_C1C2.
      DOUBLE PRECISION FRECOMB_ECARTMAX_C1C2  ! Recombined value of F giving ECARTMAX_C1C2.
            
      INTEGER*4 LUM_NBMU ! Number of actual angles of the simulations
      INTEGER*4 OS_NB ! Order of the expansions of the phase functions
			 ! in Legendre polynomials.
 	    
      INTEGER*4 IL	! Maximum useful order for the function E.
      INTEGER*4 I1	! Index of the angle number Teta1.
      INTEGER*4 I2	! Index of the angle number Teta2.
      INTEGER*4 IS	! Order IS of the Fourier series decomposition.
      INTEGER*4 IS2	! Order IS of the Fourier series decomposition.
      INTEGER*4 IER  ! error code =0 if no error, =-1 otherwise
      INTEGER*4 I
      INTEGER*4 COMPT	! Number of configurations for the sum of the relative differences.

      	! Element M11 of the reflection matrix versus the order IS
	! of the Fourier series, for an incidence at Teta1 and a
	! reflection at Teta2
      REAL M11(0:CTE_OS_NB_MAX)		
      
      LOGICAL TRACE

      CHARACTER*CTE_LENFIC2 FIC_RES_MAT_REFLEX 
                     ! Name of the file RES_MAT_REFLEX
                     ! with complete directory tree
     
	
C* Initialisation
C---------------------------------------------------------	
      ECARTMAX = -1.D+300
      SOMECART = 0.
      COMPT = 0
      

C* Creation of the result file of SOS_FSF_ROUJEAN
C---------------------------------------------------------
      OPEN(UNIT=4,FILE=FIC_RES_MAT_REFLEX,FORM='UNFORMATTED',ERR=991)

C* Loop on the angle Teta1 > 0.
C---------------------------------
      DO  3 I1=1 ,LUM_NBMU

C* Cosine and sine of Teta1: incidence angle
C---------------------------------
         C1=RMU(I1)
         S1=DSQRT(1-C1*C1)

C*   Loop on the angle teta2 > 0.
C---------------------------------
         DO 4 I2=1 ,LUM_NBMU	

C*   Cosine and sine of Teta2: reflection angle
C---------------------------------
           C2=RMU(I2)
           S2=DSQRT(1-C2*C2)

C*   Determination of F for phi = 0 to PI with step PI/CTE_PH_NU
C-------------------------------------------------------------         
           Q=INCTE_PI/CTE_PH_NU
           DO 6 I=0,CTE_PH_NU
              PHIOS=Q*I
	       PHIRJ = INCTE_PI - PHIOS	!Adaptation to the Roujean convention
              CALL SOS_CALC_F_ROUJEAN(K0,K1,K2,C1,S1,C2,S2,PHIRJ,F)     
              U(I)=F
	      
	      IF(F.LT.0.D+00) GOTO 993	!BRDF < 0 => error case.
 6         CONTINUE

		
C*   Determination of the min and max values of F
C*   for fixed I1 and I2 but varying PHI => for the trace file
C--------------------------------------------------------------
	   IF (TRACE) THEN
	       FMAX=-1.D+300
	       FMIN=1.D+300
	       DO I=0,CTE_PH_NU
	          IF(U(I).GT.FMAX) THEN
                   FMAX = U(I)
		     PHIFMAX = Q*I		!Azimuth convention OS
                 ENDIF   
	          IF(U(I).LT.FMIN) THEN
		     FMIN = U(I)
		     PHIFMIN = Q*I		!Azimuth convention OS
	          ENDIF   	 
	       ENDDO

	       WRITE(99,1000,err=921) ACOS(C1)*180./INCTE_PI, 
     &	                              ACOS(C2)*180./INCTE_PI,
     &                               PHIFMIN*180./INCTE_PI,FMIN,
     &				         PHIFMAX*180./INCTE_PI,FMAX	          	 
           ENDIF !Trace 
	      


	 
C*   Fourier series decomposition.
C*   Loop on the order of the series.
C--------------------------------------------------------------------

C--------- Initialisations
C--------------------------------       
          DO IS=0,OS_NB
	      E(IS) = 0.D+00
	   ENDDO
	   
	   B1_PREC = 1.D+300		


C--------- Loop on the order IS
C-------------------------------- 	   
	   DO 13 IS=0,OS_NB
              

C--------- Determination of E(IS)
C--------------------------------      
	      Y = 0.D+00      
	      DO I=0,CTE_PH_NU
	         PHIOS = I*Q		 		!Azimuth convention OS
	         Y=Y+U(I)*DCOS(IS*PHIOS)
	      ENDDO
	      E(IS) = Y*Q / INCTE_PI
	      
C--------- Determination of the maximum difference between exact and recombined F: B1
C---------------------------------------------------------------------      
	      B1=0.D+00			!initialisation
	      DO I=0,CTE_PH_NU
	         PHIOS=Q*I
	         PHIRJ = INCTE_PI - PHIOS	!Adaptation to the Roujean convention
                CALL SOS_CALC_F_ROUJEAN(K0,K1,K2,C1,S1,C2,S2,PHIRJ,F)     

C-------------- Computation of the recombined F from the function E(IS) for fixed PHI.
                 T1 = E(0)
		 DO IS2 = 1, IS		 
		     T1 = T1+2.*E(IS2)*DCOS(IS2*PHIOS)	!Azimuth convention OS		     
		 ENDDO !End of loop on IS2
		 
		 DO IS2 = 0, IS		 
		     IF(IS2.EQ.0) THEN
	                T1=E(0)			
	             ELSE 
	                T1=T1+2.*E(IS2)*DCOS(IS2*PHIOS)	!Azimuth convention OS
	             ENDIF		     
		 ENDDO !End of loop on IS2
       
                 B1=MAX( DABS((T1-F)/F) , B1)
		 		 
	      ENDDO !End of loop on I (phi = 0 to PI)
	           	     

C--------- Stop test of the decomposition:
C---------    * Stop if the maximum relative difference between exact and recombined F
C---------      is below the threshold CTE_SEUIL_SF_ROUJEAN
C---------      or if it increases with respect to the previous order.
C---------    * Otherwise the loop on IS goes on.
C-------------------------------------------------------------------------------- 		
              IF (B1.GT.CTE_SEUIL_SF_ROUJEAN) GOTO 14
              IL=IS	      
              GOTO 17	!Exit of the loop on IS
	      
 14           IF(B1.LT.B1_PREC) GOTO 141
              IL=IS-1
	      GOTO 17	!Exit of the loop on IS
	      
 141	      B1_PREC=B1	        
 13        CONTINUE	!End of loop on IS
	   IL=OS_NB


	     



C------------------------------------------------------------------------
C* Storage of the Fourier series decomposition of the Roujean
C* reflection matrix for the directions Teta1 and Teta2
C--------------------------------------------------------------------------

C* Type conversion for each order IS: double precision -> real
C---------------------------------------------------------------------
17          DO IS=0,OS_NB
                M11(IS)=E(IS)  		               
            ENDDO

C* Recording for the pair (Teta1, Teta2) of the decomposition
C* in Fourier series of the BRDF (converted into normalised radiance
C* and adapted to the OS code).
C-----------------------------------------------------------------
            WRITE(4,err=992) M11



C* Trace file: Computation of quality analysis parameters
C*                 of the Fourier series decomposition
C*                 and recordings.
C-----------------------------------------------------------------
           IF (TRACE) THEN

	      ECARTMAX_C1C2=0.
	      	      
	      Q=INCTE_PI/CTE_PH_NU
	      DO I=0,CTE_PH_NU
	         PHIOS=I*Q
		  PHIRJ = INCTE_PI - PHIOS		!Adaptation to the Roujean convention
	         CALL SOS_CALC_F_ROUJEAN(K0,K1,K2,C1,S1,C2,S2,PHIRJ,F)
		 
                FRECOMB = E(0)
	         DO IS=1,IL
	            FRECOMB = FRECOMB + 2.*DCOS(IS*PHIOS)*E(IS)
	         ENDDO
            
	         IF(F.NE.0.D+00) THEN
	            ECART = DABS((FRECOMB-F)/F)   

C*       Analysis for all the azimuths of the current pair (Teta1, Teta2)
		     IF(ECART.GT.ECARTMAX_C1C2) THEN
		        ECARTMAX_C1C2 = ECART
		        PHIECARTMAX_C1C2 = PHIOS*180./INCTE_PI	!Azimuth convention OS
	               F_ECARTMAX_C1C2 = F
	               FRECOMB_ECARTMAX_C1C2 = FRECOMB
		     ENDIF 
		    
C*       Analysis for all the azimuths and all the pairs (Teta1, Teta2)
	            IF(ECART.GT.ECARTMAX) THEN
	               ECARTMAX = ECART
	               PHIECARTMAX = PHIOS*180./INCTE_PI	!Azimuth convention OS
	               TVECARTMAX = ACOS(C2)*180./INCTE_PI
	               TSECARTMAX = ACOS(C1)*180./INCTE_PI
	               F_ECARTMAX = F
	               FRECOMB_ECARTMAX = FRECOMB
	            ENDIF   
	    
	            SOMECART=SOMECART + ECART
	            COMPT = COMPT + 1
		 ENDIF
	      ENDDO !End of loop on I (phi)
         	
	      	
	      IF(F_ECARTMAX_C1C2.NE.0.D+00) THEN	 
	         WRITE(99,2000,err=921) IL,PHIECARTMAX_C1C2,
     &		 			    F_ECARTMAX_C1C2,
     &		                         FRECOMB_ECARTMAX_C1C2,
     &					    ECARTMAX_C1C2*100.
	      ELSE
	      	 WRITE(99,2500,err=921) IL,PHIECARTMAX_C1C2,
     &		 			   F_ECARTMAX_C1C2,
     &		                        FRECOMB_ECARTMAX_C1C2
	      ENDIF	 
	      WRITE(99,*,err=921)   	    
                 
	      
           END IF !End of Trace
      
      	
		                            
 4       CONTINUE	!End of loop on Teta2 (I2)
 
 3    CONTINUE	!End of loop on Teta1 (I1)
 
 


C* Trace file: Recording of the quality analysis parameters
C*                 of the Fourier series decomposition.
C----------------------------------------------------------------- 
      IF (TRACE) THEN
         WRITE(99,*,err=921)
	  WRITE(99,*,err=921)
         WRITE(99,*,err=921)  
     &	  "Analysis of the Fourier expansion of the ROUJEAN's BRDF: "
	
	 IF(COMPT.GT.1) THEN
	    WRITE(99,*,err=921)
	    WRITE(99,*,err=921)
	    WRITE(99,*,err=921)  
     &	    'Mean relative difference from the recombination'//
     &	    ' (% excluding null values):', SOMECART*100/COMPT
	    WRITE(99,*,err=921)
	    WRITE(99,*,err=921)  
     &	    'Max relative difference from the recombination (%):',
     &	    ECARTMAX*100
           WRITE(99,*,err=921)   
     &	    'Related delta phi angle (SOS convention) :',PHIECARTMAX
           WRITE(99,*,err=921)    
     &	    'Related viewing angle :',TVECARTMAX
	    WRITE(99,*,err=921)    
     &	    'Related solar angle :',TSECARTMAX
	    WRITE(99,*,err=921)    
     &	    'Related value of BRDF :',F_ECARTMAX
	    WRITE(99,*,err=921)   
     &	    'Related value of the recombined BRDF :',FRECOMB_ECARTMAX 
        ELSE
	    WRITE(99,*,err=921)
	    WRITE(99,*,err=921) 'WARNING : BRDF forced to zero !'
	    WRITE(99,*,err=921) 'Unsuitable simulation parameters'
	 ENDIF 
	 
	 WRITE(99,*,err=921)
	 WRITE(99,*,err=921)
	 
      END IF 
      
      

       
      
      
C* Nominal end
C-------------------
      IER=0
      GOTO 9999
      
C* Error case: return of status -1 to the calling program
C------------------------------------------------------
  921 WRITE(6,*) 'Erreur d''écriture sur le fichier trace Surface'
      IER=-1
      GOTO 9999
  991 WRITE(6,*) 'SOS_FSF_ROUJEAN : Erreur a l''ouverture d''un fichier'
      IER=-1
      GOTO 9999
  992 WRITE(6,*) 'SOS_FSF_ROUJEAN',
     &		 'Erreur de lecture ou d''écriture sur un fichier'
      IER=-1 
      GOTO 9999
  993 WRITE(6,*) 'SOS_FSF_ROUJEAN : ATTENTION : PRESENCE DE BRDF < 0 !',
     &		 'Paramètres de simulation inadaptes.'
      IER=-1 
      GOTO 9999
      
 9999 CONTINUE


C* File closing: 'RES_MAT_REFLEX'
C-------------------------------------
      CLOSE(4)
 
 
C* Format
C--------

1000  FORMAT(6hTeta1:,1X,F5.2,2X,6hTeta2:,1X,F5.2,2X,
     &15h; (PHIos,Fmin):,1X,F6.2,2X,E11.4,2X,
     &15h; (PHIos,Fmax):,1X,F6.2,2X,E11.4)    
 
2000  FORMAT(5X,12hLimit order:,I3,3h ; ,14h Error max => ,1X,4hphi:,
     &1X,F6.2,3h, F:,1X,E11.4,2X,11h, F recomb:,1X,E11.4,2X,
     &12h, Ecart (%):,1X,F7.3) 
     
2500  FORMAT(5X,12hLimit order:,I3,3h ; ,14h Error max => ,1X,4hphi:,
     &1X,F6.2,3h, F:,1X,E11.4,2X,11h, F recomb:,1X,E11.4)                               
      RETURN
      END         !END OF PROCEDURE SOS_FSF_ROUJEAN





C==============================================================================
C PROCEDURE: SOS_CALC_F_ROUJEAN
C ==========
C      This procedure computes the BRDF of the Roujean model:
C	 F(Teta1,Teta2,phi) = K0 + K1*F1(Teta1,Teta2,phi) + K2*F2(Teta1,Teta2,phi)
C
C	 with Teta1 : incidence zenith angle: tetas  (0 < Teta1 < PI/2),
C             Teta2 : reflection zenith angle: tetav(0 < Teta2 < PI/2),
C             phi   : relative azimuth angle between the directions
C                     of incidence and reflection with the convention
C                     of the Roujean model  (0 <= phi <= PI),
C
C	      K0, K1, K2 : coefficients of the Roujean model,
C			
C            F1 = [ (PI-phi)*cos(phi)+sin(phi) ]*tgs*tgv / 2*PI
C                 -[ tgs + tgv + SQRT{tgs**2+tgv**2-2*tgs*tgv*cos(phi)} ] / PI
C
C		  with tgs = tan(tetas) and tgv = tan(tetav),
C
C            F2 = [ (PI/2-zeta)*cos(zeta)+sin(zeta) ] * (4/3*PI) 
C                 / [ cos(tetas)+cos(tetav) ]   - 1/3
C
C                 with zeta the phase angle defined by:
C                 cos(zeta) = cos(tetas)*cos(tetav)+sin(tetas)*sin(tetav)*cos(phi). 
C
C
C      The BRDF of the Roujean model is converted from reflectance into a
C      quantity useful for the OS: normalised radiance * cos(Teta2),
C      i.e.: BRDF * cos(Teta1) * cos(Teta2).
C
C Description of parameters
C --------------------------
C       K0	(double)   (I) Coefficient of the Roujean model.
C       K1	(double)   (I) Coefficient of the Roujean model.
C       K2	(double)   (I) Coefficient of the Roujean model.
C	 C1	(double)   (I) Cos(Teta1).
C	 S1	(double)   (I) Sin(Teta1).
C	 C2	(double)   (I) Cos(Teta2).
C	 S2	(double)   (I) Sin(Teta2).
C       PHI     (double)   (I) Relative azimuth for the convention
C                              of the Roujean model (radians).
C 	 F	(double)   (O) Value of F(Teta1,Teta2,phi).
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
      SUBROUTINE SOS_CALC_F_ROUJEAN(K0,K1,K2,C1,S1,C2,S2,PHI,F)

      IMPLICIT NONE

C* Definition of variables
C*-----------------------------------------------------------------------
      DOUBLE PRECISION K0	! (I) Parameter of the Roujean model.
      DOUBLE PRECISION K1	! (I) Parameter of the Roujean model.
      DOUBLE PRECISION K2	! (I) Parameter of the Roujean model.

      DOUBLE PRECISION C1	! (I)  Cos(Teta1).
      DOUBLE PRECISION S1	! (I)  Sin(Teta1).
      DOUBLE PRECISION C2	! (I)  Cos(Teta2).
      DOUBLE PRECISION S2	! (I)  Sin(Teta2).
      DOUBLE PRECISION PHI	! (I) Relative azimuth for the convention
      				!     of Roujean (radians).
      
      DOUBLE PRECISION xPHI	! Value of PHI actually used:
      				! => value between 0 and PI.
      
      DOUBLE PRECISION xC1	! Value of C1 actually used.
      DOUBLE PRECISION xS1	! Value of S1 actually used.
      DOUBLE PRECISION xC2	! Value of C2 actually used.
      DOUBLE PRECISION xS2	! Value of S2 actually used.
      
      DOUBLE PRECISION F	! (O) Value of the Roujean BRDF
      				! for the observation geometry
				! (converted into normalised radiance * cos(teta2).
      
      DOUBLE PRECISION F1	! Function of the ROUJEAN BRDF.
      DOUBLE PRECISION F2	! Function of the ROUJEAN BRDF.
      
      DOUBLE PRECISION COSPHI	! Cosine of the relative azimuth.
      DOUBLE PRECISION TANTS	! Tangent of the incidence angle.
      DOUBLE PRECISION TANTV	! Tangent of the reflection angle.
            
      DOUBLE PRECISION ZETA	! Phase angle.
      DOUBLE PRECISION COSZETA	! Cosine of the phase angle between the
      				       ! directions of incidence and reflection.

      DOUBLE PRECISION PI	! Value of PI
      

C* Processing
C---------------------------------    
      PI=4.D+00*DATAN(1.D+00)

C*   Check of the angles
C---------------------------------------------------------      

C*   Conversion of PHI between 0 and PI
      xPHI = PHI
      IF(xPHI.LT.0.D+00) xPHI = - xPHI
      IF(xPHI.GT.PI) xPHI = 2.*PI - xPHI

      
C*   Check of the viewing and illumination zenith angles
      xC1=C1
      xS1=S1
      xC2=C2
      xS2=S2  
      
      IF (DACOS(C1)*180./PI.GT.CTE_TETAS_LIM_ROUJEAN) THEN
          xC1 = DCOS(CTE_TETAS_LIM_ROUJEAN*PI/180.D+00)
 	   xS1 = DSIN(CTE_TETAS_LIM_ROUJEAN*PI/180.D+00)
      ENDIF
      IF (DACOS(C2)*180./PI.GT.CTE_TETAV_LIM_ROUJEAN) THEN
          xC2 = DCOS(CTE_TETAV_LIM_ROUJEAN*PI/180.D+00)
 	   xS2 = DSIN(CTE_TETAV_LIM_ROUJEAN*PI/180.D+00)
      ENDIF	
	
	
C*   Computation of cos(phi) and of the tangents of Tetas and Tetav
C---------------------------------------------------------     	
      COSPHI = DCOS(xPHI)
      TANTS =  xS1 / xC1
      TANTV =  xS2 / xC2

C*   Computation of F1
C---------------------------------------------------------      
      F1 = 0.5 * ( (PI-xPHI)*COSPHI + DSIN(xPHI) ) * TANTS * TANTV 
      F1 = F1 - TANTS - TANTV
      F1 = F1 - DSQRT(TANTS*TANTS + TANTV*TANTV -2.*TANTV*TANTS*COSPHI)
      F1 = F1 / PI
      
      
C*   Computation of F2
C---------------------------------------------------------            
      COSZETA = xC1*xC2+xS1*xS2*COSPHI  

      !Management of the case of a small numerical rounding error giving a value beyond 1 or -1
      IF ( DABS(DABS(COSZETA)-1.D+00).LE.CTE_SEUIL_NUM ) THEN
         IF ( (COSZETA.GE.(1.-CTE_SEUIL_NUM))
     &       .AND.(COSZETA.LE.(1.+CTE_SEUIL_NUM)) ) THEN
            COSZETA = 1.D+00
         ELSE
            COSZETA = -1.D+00
         ENDIF
      ENDIF

      ZETA = DACOS(COSZETA)

      F2 = 4.*( (PI/2. - ZETA)*COSZETA + DSIN(ZETA) ) 
     & / ( 3.*PI*(xC1+xC2))
      F2 = F2 - (1.D+00/3.D+00)



C*   Computation of F
C---------------------------------------------------------            
      F = K0 + K1*F1 + K2*F2

C*   Conversion reflectance -> normalised radiance (*C1)
C*   and adaptation to the OS -> *C2.
C*
C*   NB : 
C*    - The adaptation to the OS code by multiplication by
C*      cos(tetav) must use the value corresponding to the
C*      Gauss angle file because this is the value
C*      used in the computations of the OS code
C*      (C2 is used and not xC2).
C*    - Likewise, the conversion from reflectance to normalised radiance
C*      is made for the true value of cos(TetaIncident), not the one
C*      thresholded by CTE_TETAS_LIM_ROUJEAN, because the flux computations
C*      in SOS_OS need the true incidence angle:
C*      C1 is used and not xC1.
                   
C--------------------------------------------------------- 	
      F = F*C2*C1	     

      RETURN
      END         !END OF PROCEDURE SOS_CALC_F_ROUJEAN










C==============================================================================
C PROCEDURE:  SOS_MISE_FORMAT_RJ
C ==========
C      This procedure reorders the elements of the matrix Mij into a matrix
C      Pij usable by the OS code.
C
C      The elements M11, stored in the file "RES_MAT_REFLEX", give
C      for each pair of angles (Teta1, Teta2) the associated Fourier series
C      decomposition. The other (polarisation) elements are not
C      precomputed since they are set to zero.
C
C      The elements Pij, stored in the output file, give
C      for each order of the Fourier series decomposition
C      the reflection associated with all the pairs (Teta1, Teta2).
C
C Description of parameters
C --------------------------
C     FIC_RES_MAT_REFLEX (CHARACTER)  : (I) Name of the temporary file RES_MAT_REFLEX
C                                           with complete directory tree.

C      FIC  (CHARACTER*CTE_LENFIC2)  (I)  Complete name of the result file.
C      NBMU (I4)  : (I)  Number of actually useful (positive) angles.
C      OS_NB (I4) : (I)  Order of the expansions of the phase functions
C		         in Legendre polynomials to be used.
C      IER  (I4)   (O) : error code =0 if no error, =-1 otherwise
C
C
C Other result provided by the procedure
C -------------------------------------
C      Result file: FIC, BRDF file of the Roujean model.
C 
C Description of the files used
C ---------------------------------
C   -->Content of the temporary file RES_MAT_REFLEX
C       Unformatted binary file, created by SOS_FSF_ROUJEAN
C       File created in the directory where the program is run
C       Record for each pair (Teta1, Teta2) of
C       the Fourier series decomposition of: M11 (n real values)
C
C   -->Content of the ROUJEAN BRDF file
C       Unformatted binary file, created by SOS_MISE_FORMAT_RJ
C       It contains the surface reflection matrix
C       i.e. the elements Pij of the reflection matrix
C              P11(CTE_OS_NBMU,CTE_OS_NBMU) (real),
C              P12(CTE_OS_NBMU,CTE_OS_NBMU) (real),
C              P13(CTE_OS_NBMU,CTE_OS_NBMU) (real),
C              P21(CTE_OS_NBMU,CTE_OS_NBMU) (real),
C              P22(CTE_OS_NBMU,CTE_OS_NBMU) (real),
C              P23(CTE_OS_NBMU,CTE_OS_NBMU) (real),
C              P31(CTE_OS_NBMU,CTE_OS_NBMU) (real),
C              P32(CTE_OS_NBMU,CTE_OS_NBMU) (real),
C              P33(CTE_OS_NBMU,CTE_OS_NBMU) (real)
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C     -Opening of a file
C      The file RES_MAT_REFLEX does not exist
C      The Roujean BRDF file exists
C     -Reading/writing in a file
C
C     Message displayed on screen, the routine stops its computations and
C     returns status -1 to the calling program
C
C
C=============================================================================
      SUBROUTINE SOS_MISE_FORMAT_RJ(FIC_RES_MAT_REFLEX,FIC,
     &                              NBMU,OS_NB,IER)

      IMPLICIT NONE

C* Definition of variables
C*--------------------------------------------------------------------
      CHARACTER*CTE_LENFIC2 FIC	! Name of the result file.

      CHARACTER*CTE_LENFIC2 FIC_RES_MAT_REFLEX 
                     ! Name of the file RES_MAT_REFLEX
                     ! with complete directory tree

      	! Element M11 of the reflection matrix versus the order IS
	! of the Fourier series, for an incidence at Teta1 and a
	! reflection at Teta2.
      REAL M11(0:CTE_OS_NB_MAX)		
      

      	! Elements Pij of the reflection matrix versus the directions
	! of incidence and reflection, for an order IS of the expansion
	! in Fourier series
      REAL P11(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P12(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P13(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)	
      REAL P21(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P22(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P23(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P31(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P32(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P33(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)

      INTEGER*4 IS	! Order of the Fourier series decomposition.
      INTEGER*4 I	! Index of the incidence angle Teta1.
      INTEGER*4 J	! Index of the reflection angle Teta2.
      INTEGER*4 IER     ! error code =0 if no error, =-1 otherwise
      INTEGER*4 NBMU    ! Number of actual angles of the simulations
      INTEGER*4 OS_NB ! Order of the expansions of the phase functions
			 ! in Legendre polynomials.
       
C* File opening
C------------------------------------------------------------
      OPEN(UNIT=1,FILE=FIC_RES_MAT_REFLEX,FORM='UNFORMATTED',
     &     STATUS='OLD',ERR=991)
      OPEN(UNIT=2,FILE=FIC,FORM='UNFORMATTED',
     &     STATUS='NEW',ERR=991)

C* Loop on the order of the Fourier series decomposition
C------------------------------------------------------------
      DO 1 IS=0,OS_NB

C* Loop on the incidence directions Teta1
C---------------------------------------------
         DO 2 I=1,NBMU

C* Loop on the reflection directions Teta2
C-----------------------------------------------
            DO 3 J=1,NBMU
               READ(1,err=992) M11
               P11(I,J)=M11(IS)
               P12(I,J)=0. 
               P13(I,J)=0.  
               P21(I,J)=0.
               P22(I,J)=0. 
               P23(I,J)=0.
               P31(I,J)=0.
               P32(J,I)=0.
               P33(J,I)=0.
 3          CONTINUE
 2       CONTINUE

C*       Recording of the arrays of the reflection matrix
C*       at order IS of the Fourier series decomposition.
C*       Record bounds of the arrays limited to the
C*       useful dimension -NBMU:NBMU
C-------------------------------------------------------------- 
	 WRITE(2,err=992) ((P11(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P12(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P13(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P21(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P22(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P23(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P31(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P32(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P33(I,J),I=1,NBMU),J=1,NBMU)

C* Back to the beginning of the file before the next order IS.
C---------------------------------------------------------------
         REWIND (1)
 1    CONTINUE


      
C* Nominal end
C-------------------
      IER=0
      GOTO 9999
      
C* Error case: return of status -1 to the calling program
C------------------------------------------------------
  991 WRITE(6,*) 'SOS_MISE_FORMAT_RJ : '//
     &           'Erreur a l''ouverture d''un fichier'
      IER=-1
      GOTO 9999
  992 WRITE(6,*) 'SOS_MISE_FORMAT_RJ : '//
     &           'Erreur de lecture ou d''écriture sur un fichier'
      IER=-1 
      GOTO 9999

 9999 CONTINUE

C* File closing
C-------------------
      CLOSE(1)
      CLOSE(2)
      
      RETURN
      END         !END OF PROCEDURE SOS_MISE_FORMAT_RJ





