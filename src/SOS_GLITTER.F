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
C* FILE: SOS_GLITTER.F
C* PROJECT: Successive orders of scattering
C* ROLE: Generates the "glitter" files for the glitter computation.
C*
C* AUTHOR:
C*   Initial code: Laboratoire d'Optique Atmospherique (LOA, Villeneuve d'Ascq).
C*   Code rework and development: CS GROUP France, B.Lafrance
C* DATE: 30/04/02
C*
C* MOD:VERSION:1.0: Upgrade of the LOA OS code
C* MOD:VERSION:1.1: Update of the comments
C*                  
C* MOD:VERSION:2.0: 10/01/2010 
C*
C*   - Removal of the constant SOS_OS_FICANGLE (from SOS.h) and of the reading of the corresponding
C*     Gauss angle file, in favour of passing the angle (RMU) and weight (CHR) tables as arguments.
C*
C*   - Change of the size of the radiance arrays:
C*     ==> Changed from CTE_OS_NBMU to CTE_OS_NBMU_MAX.
C*     ==> Loop on the angles up to LUM_NBMU (new input argument)
C*         equal to the actual number of angles used.
C*
C*   - Removal of the constants CTE_OS_NB, CTE_OS_NS and CTE_OS_NM (from SOS.h)
C*     in favour of passing the values OS_NB, OS_NS, OS_NM as arguments.
C*
C*   - Change of the size of the arrays of the decomposition in Legendre polynomials:
C*     ==> Changed from CTE_OS_NB to CTE_OS_NB_MAX.
C*     ==> Changed from CTE_OS_NS to CTE_OS_NS_MAX.
C*     ==> Changed from CTE_OS_NM to CTE_OS_NM_MAX.
C*     ==> Loop on the orders up to OS_NB, OS_NS or OS_NM (new input arguments)
C*
C*   - Routine SOS_GSF: change of the routine call to include the argument LUM_NBMU (nb of useful angles)
C*     and to limit the loops to the useful bounds LUM_NBMU, and to include the argument OS_NM.
C*
C*   - Adaptation of the calls of routines SOS_MAT_FRESNEL, SOS_MAT_REFLEXION
C*     and SOS_MISE_FORMAT to the passing of the parameters LUM_NBMU (nb of useful angles),
C*     RMU (table of the cosines of the angles) and CHR (table of the weights of the angles),
C*     as well as OS_NB, OS_NS and OS_NM.
C*
C* MOD:VERSION:6.0: 06/06/2013
C*
C*   - Change of the value of Pi, declared as a constant in the file
C*     (INCTE_PI replaces CTE_PI).
C*
C* MOD:VERSION:7.0: 16/06/2020
C*
C*   - Names of the temporary files (RES_GSF, RES_FRESNEL, RES_MAT_REFLEX) passed
C*     with complete directory tree.
C*   - Consequent change of the arguments given to routine SOS_GSF.
C*   - Routine SOS_GSF called with the parameter SIG (of Cox&Munk) and routine
C*     SOS_MAT_REFLEXION with the parameter (1./SIG) instead of the wind speed.
C*
C* MOD:VERSION:7.1: 31/08/2022
C*    - Constants renamed from SOS_* to CTE_*
C*    - Traces converted to English.
C*
C*
C******************************************************************************

C----------------------------------------------------------------------------
C Definition of constants
C---------------------------------------------------------------------------- 
C Constants used:
C    CTE_LENFIC2 : Length of file names with directory tree.
C    CTE_OS_NBMU_MAX : Maximum number of positive angles.
C    INCTE_PI : Value of PI.
C    CTE_OS_NB_MAX :  Maximum value of the limit order of the expansion in
C                     Legendre polynomials.
C    CTE_OS_NS_MAX :  Maximum value of the limit order of the expansion in Legendre
C                     polynomials of the Fresnel matrix elements
C                     and in Fourier series of the radiances.
C    CTE_OS_NM_MAX :  Maximum value of the limit order of the expansion
C                     in Fourier series of the G function.
C    CTE_PH_TEST : Factor for the comparison of GMIN and GMAX
C    CTE_PH_NU : Number of azimuth angles (2**NQ)
C    CTE_PH_NQ : Power of 2 of CTE_PH_NU (NU=2**NQ)
C
C----------------------------------------------------------------------------
#include "SOS.h"
#define INCTE_PI DACOS(-1.D+00)

C==============================================================================
C PROCEDURE: SOS_GLITTER
C ==========      
C      This procedure generates the surface reflection files
C      for the glitter computation.
C
C
C Input data of the program
C ------------------------------
C    -->Parameters/arguments of the program
C       LUM_NBMU (I4) : (I) Number of actually useful (positive) angles
C	 RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (DOUBLE) : (I) Cosines of the angles.
C	 CHR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (DOUBLE) : (I) Weights of the angles.
C 	 WIND  (DOUBLE)  : (I) Wind speed in m/s.
C       IND   (DOUBLE)  : (I) Refractive index water / air.
C       OS_NB (I4) : (I) Order of the expansions of the phase functions
C		           in Legendre polynomials to be used.
C       OS_NS (I4) : (I) Maximum order to be used for the expansion in
C                        Legendre polynomials of the Fresnel matrix elements
C                        and in Fourier series of the radiance computations.
C       OS_NM (I4) : (I) Maximum order of the expansion in Fourier series
C                        of the G function of the reflection matrix computations.
C       FIC_RES_FRESNEL (CHARACTER*CTE_LENFIC2)  : 
C                 (I) Name of the temporary file RES_FRESNEL with complete directory tree.
C       FIC_RES_GSF (CHARACTER*CTE_LENFIC2)  :  
C                 (I) Name of the temporary file RES_GSF with complete directory tree.
C       FIC_RES_FRESNEL (CHARACTER*CTE_LENFIC2)  :  
C                 (I) Name of the temporary file RES_FRESNEL with complete directory tree.
C       FIC_RES_MAT_REFLEX (CHARACTER*CTE_LENFIC2)  :  
C                 (I) Name of the temporary file RES_MAT_REFLEX with complete directory tree.
C       FICGLITTER  (CHARACTER*CTE_LENFIC2) : (O) complete name of the result Glitter file.
C       TRACE  (LOGICAL)   ! = true if writing in the trace file.
C       IER  (I4)   (O) : error code =0 if no error, =-1 otherwise
C
C
C Result provided by the program
C --------------------------------
C   Glitter file
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
C   -->Content of the Glitter file
C       Unformatted binary file, created by SOS_MISE_FORMAT
C       It contains the surface reflection matrix,
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
C   -->Content of the temporary file RES_GSF
C       Unformatted binary file, created by SOS_GSF
C       File created in the directory where the program is run
C       If a file with the same name already exists, it is overwritten
C       Each record contains:
C           -number of the angle Teta1 (integer value)
C           -number of the angle Teta2 (integer value)
C           -maximum useful order of the decomposition: IL  (integer value)
C           -Decomposed G function for (Teta1,Teta2) : (E(J),J=0,IL)
C            (IL+1 double precision values)
C
C   -->Content of the temporary file RES_FRESNEL
C      Ascii file, record format 4E15.8,
C      created by SOS_MAT_FRESNEL
C      File created in the directory where the program is run
C      If a file with the same name already exists, it is overwritten
C      Number of records: from k=0 to OS_NS
C      Each record contains the parameters ALPHA(k), BETA(k),
C      GAMMA(k) and ZETA(k) associated with the Fresnel matrix
C
C   -->Content of the temporary file RES_MAT_REFLEX
C       Unformatted binary file, created by SOS_MAT_REFLEXION
C       File created in the directory where the program is run
C       If a file with the same name already exists, it is overwritten
C       Record for each pair (Teta1, Teta2) of
C       the Fourier series decomposition of the Fresnel matrix:
C       M11,M12,M13,M21,M22,M23,M31,M32,M33 
C       (n real values)
C
C    The 3 temporary files are overwritten if they exist, then
C    deleted at the end of the program.
C
C
C Common used:
C --------------
C  None
C
C
C Error cases:
C ------------
C     -Opening of a file
C      Error if the Glitter file already exists
C     -Reading/writing in a file
C
C     Message displayed on screen, program stop and
C     return of status 1 to the shell
C        
C
C  Subroutines used:
C --------------------------
C  This procedure calls the routines:
C      - SOS_GSF
C      - SOS_MAT_FRESNEL 
C      - SOS_MAT_REFLEXION 
C      - SOS_MISE_FORMAT
C
C=============================================================================

      SUBROUTINE SOS_GLITTER(LUM_NBMU,RMU,CHR,WIND,IND,
     &                       OS_NB,OS_NS,OS_NM,
     &                       FIC_RES_GSF,FIC_RES_FRESNEL,
     &                       FIC_RES_MAT_REFLEX,
     &                       FICGLITTER,TRACE,IER)

      IMPLICIT NONE


C* Definition of variables
C*--------------------------------------------------------------------------
      CHARACTER*CTE_LENFIC2 FICGLITTER	! Complete name of the result GLITTER file.

      CHARACTER*CTE_LENFIC2 FIC_RES_GSF ! Name of the file RES_GSF
                                        ! with complete directory tree

      CHARACTER*CTE_LENFIC2 FIC_RES_FRESNEL ! Name of the file RES_FRESNEL
                                            ! with complete directory tree

      CHARACTER*CTE_LENFIC2 FIC_RES_MAT_REFLEX ! Name of the file RES_MAT_REFLEX
                                               ! with complete directory tree
            
      DOUBLE PRECISION WIND	! Wind speed (m/s).

      DOUBLE PRECISION SIG 	! Square of the parameter Sigma of the Cox and Munk theory.
      
      DOUBLE PRECISION IND	! Refractive index water / air
	
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   ! Array of the Gauss mu
      DOUBLE PRECISION CHR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   ! Gauss weights.
      
      LOGICAL TRACE     ! = true if writing in the trace file

      INTEGER*4 OS_NB   ! Order of the expansions of the phase functions
			   ! in Legendre polynomials.
      INTEGER*4 OS_NS   ! Maximum order of the expansion in Legendre polynomials
			   ! of the Fresnel matrix elements.
      INTEGER*4 OS_NM   ! Maximum order of the expansion in Fourier series
			   ! of the G function of the reflection matrix computations.
			
      INTEGER*4 LUM_NBMU ! Number of actual angles of the simulations
      INTEGER*4 I	    ! Array index.
      INTEGER*4 IFIN     ! nb of significant characters of a string
      INTEGER*4 IER      ! error code =0 if no error, =-1 otherwise




C* Initialisation of the error code
C------------------------------------------------------
      IER=0

      
C* Opening of the trace file
C------------------------------------------------------
      IF (TRACE) THEN
         WRITE(99,*,err=921)	 
         WRITE(99,*,err=921) 'Calculation of Sun Glint file'  
         WRITE(99,*,err=921) ' '
         WRITE(99,*,err=921)'Simulation parameters '
         WRITE(99,*,err=921)' Wind velocity  (m/s) ',WIND
         WRITE(99,*,err=921)
     &        ' Refractive index of sea water ',IND
         WRITE(99,*,err=921)' OS_NB ',OS_NB
	  WRITE(99,*,err=921)' OS_NS ',OS_NS
	  WRITE(99,*,err=921)' OS_NM ',OS_NM
	  WRITE(99,*,err=921) ' ' 	 	 
      END IF
        
C* Square of the parameter Sigma of the Cox and Munk theory.
C--------------------------------------------------------
      SIG = .003 + .00512*WIND

C* Fourier series decomposition of the function G(Teta1,Teta2,phi).
C---------------------------------------------------------------------
      CALL SOS_GSF(LUM_NBMU,RMU,SIG,OS_NM,FIC_RES_GSF,IER)
      IF (IER.NE.0) GOTO 995   ! error case
      
C* Expansion of the Fresnel matrix in Legendre functions.
C------------------------------------------------------------------
      CALL SOS_MAT_FRESNEL(LUM_NBMU,RMU,CHR,IND,OS_NS,
     &                     FIC_RES_FRESNEL,TRACE,IER)
      IF (IER.NE.0) GOTO 996   ! error case
      
C* Computation of the reflection matrix expanded in Fourier series.
C------------------------------------------------------------------
      CALL SOS_MAT_REFLEXION((1./SIG),LUM_NBMU,RMU,OS_NB,OS_NS,OS_NM,
     &                       FIC_RES_FRESNEL,FIC_RES_GSF,
     &                       FIC_RES_MAT_REFLEX,IER)
      IF (IER.NE.0) GOTO 997   ! error case
      
C* Rearrangement of the results to be compatible with the OS code.
C-----------------------------------------------------------------------
      CALL SOS_MISE_FORMAT(FIC_RES_MAT_REFLEX, FICGLITTER,
     &                     LUM_NBMU,OS_NB,IER)
      IF (IER.NE.0) GOTO 998   ! error case


C* Nominal end
C-------------------
      GOTO 9999
      
C* Error case and return of status 1 to the shell
C----------------------------------------------
  912 WRITE(6,*) 'Erreur a l''ouverture du fichier angles de Gauss'
      IER=-1
      GOTO 9999
  921 WRITE(6,*) 'Erreur d''écriture sur le fichier trace Glitter'
      IER=-1
      GOTO 9999
  922 WRITE(6,*) 'Erreur de lecture sur le fichier angles de Gauss'
      IER=-1
      GOTO 9999      
  995 WRITE(6,*) 'Erreur dans la routine SOS_GSF'
      IER=-1
      GOTO 9999
  996 WRITE(6,*) 'Erreur dans la routine SOS_MAT_FRESNEL'
      IER=-1
      GOTO 9999
  997 WRITE(6,*) 'Erreur dans la routine SOS_MAT_REFLEXION'
      IER=-1
      GOTO 9999
  998 WRITE(6,*) 'Erreur dans la routine SOS_MISE_FORMAT'
      IER=-1
      GOTO 9999         
           
C* Format
C---------   
 2    FORMAT(2(D21.14))    
              
 9999 CONTINUE

C* Deletion of the temporary files
C---------------------------------------   
      OPEN(UNIT=3,FILE=FIC_RES_GSF,FORM='UNFORMATTED') 
      CLOSE(3,STATUS='DELETE')     
      OPEN(UNIT=3,FILE=FIC_RES_FRESNEL) 
      CLOSE(3,STATUS='DELETE')       
      OPEN(UNIT=3,FILE=FIC_RES_MAT_REFLEX,FORM='UNFORMATTED')
      CLOSE(3,STATUS='DELETE') 
                                                      
      RETURN
      END    !END OF PROGRAM SOS_GLITTER





C==============================================================================
C PROCEDURE: SOS_GSF
C ==========
C      This procedure performs the Fourier series decomposition
C      (on phi) of the function G(Teta1,Teta2,phi):
C
C	 G(Teta1,Teta2,phi) = exp [-tan(TetaN)**2 / SIG] / cos(TetaN)**4
C
C	 with  Teta1 : incidence zenith angle  (0 < Teta1 < PI/2),
C             Teta2 : reflection zenith angle (0 < Teta2 < PI/2),
C             phi   : relative azimuth angle between the directions
C                     of incidence and reflection,
C
C	       TetaN : zenith angle of the normal to the wave facet,
C                     computed from Teta1, Teta2 and phi.
C			
C             SIG   : square of the parameter Sigma of the Cox and Munk theory.
C
C
C      Stores the results in the file "RES_GSF" as:
C	   Number of the angle Teta1 : I1,
C         Number of the angle Teta2 : I2,
C         Maximum useful order of the decomposition: IL,
C         Decomposed G function for (Teta1,Teta2) : (E(J),J=0,IL)
C
C Description of parameters
C --------------------------
C      LUM_NBMU (I4) : (I) Number of actually useful (positive) angles.
C      RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double) : (I) Array of the Gauss mu.
C      SIG	(double) : (I)  Square of the parameter Sigma of the Cox and Munk theory.
C      OS_NM (I4) : (I) Maximum order of the expansion in Fourier series
C                       of the G function of the reflection matrix computations.
C      FIC_RES_GSF  (CHARACTER)  : (I) Name of the temporary file RES_GSF
C                                      with complete directory tree.
C
C      IER  (I4) : (O) error code =0 if no error, =-1 otherwise
C
C Other result provided by the procedure
C -------------------------------------
C      Temporary file RES_GSF
C 
C Description of the files used
C ---------------------------------
C   -->Content of the temporary file RES_GSF
C       Unformatted binary file, created by SOS_GSF
C       File created in the directory where the program is run
C       If a file with the same name already exists, it is overwritten
C       Each record contains:
C           -number of the angle Teta1 (integer value)
C           -number of the angle Teta2 (integer value)
C           -maximum useful order of the decomposition: IL  (integer value)
C           -Decomposed G function for (Teta1,Teta2) : (E(J),J=0,IL)
C            (IL+1 double precision values)
C
C
C Common used:
C --------------
C     None
C
C Error cases:
C ------------
C     -Opening of a file
C     -Reading/writing in a file
C
C     Message displayed on screen, the routine stops its computations and
C     returns status -1 to the calling program
C
C     
C  Routines used by the procedure:
C ------------------------------------
C  This program calls the routine:
C        - SOS_CALCG
C
C==============================================================================
      SUBROUTINE SOS_GSF(LUM_NBMU,RMU,SIG,OS_NM,FIC_RES_GSF,IER)
 	
      IMPLICIT NONE

C* Definition of variables
C*-------------------------------------------------------------------
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! Array of the Gauss mu

      DOUBLE PRECISION E(0:CTE_OS_NM_MAX) ! Result of the decomposition in
      				              ! Fourier series of the function G.
	
      DOUBLE PRECISION SIG 	! Square of the parameter Sigma of the Cox and Munk theory.

      DOUBLE PRECISION PHI	! Relative azimuth (radians).

      DOUBLE PRECISION G	! Value of the function G(phi).

      DOUBLE PRECISION C1	! Cos(Teta1).
      DOUBLE PRECISION S1	! Sin(Teta1).
      DOUBLE PRECISION C2	! Cos(Teta2).
      DOUBLE PRECISION S2	! Sin(Teta2).
      DOUBLE PRECISION C12	! Cos(Teta1) * Cos(Teta2).
      DOUBLE PRECISION S12	! Sin(Teta1) * Sin(Teta2).
      DOUBLE PRECISION CS12	! 0.5 * [ Cos(Teta1) + Cos(Teta2) ]**2 .

      DOUBLE PRECISION GMIN	! Minimum value kept for the function
      				! G(phi) for fixed Theta1 and Theta2.
      DOUBLE PRECISION GMAX	! Maximum value of the function G(phi) = G(0).

      DOUBLE PRECISION PHIB	! Maximum azimuth angle phi for the
      				! description of the function G(phi).
			
      DOUBLE PRECISION PHI1	! Intermediate angle for the computation of PHIB.
      DOUBLE PRECISION PHI2	! Intermediate angle for the computation of PHIB.
						
      DOUBLE PRECISION U(0:CTE_PH_NU)  ! Function G(phi) between 0 and PHIB
      			                  ! for fixed Theta1 and Theta2.

      DOUBLE PRECISION Q	! Step in phi.

      DOUBLE PRECISION T2	! Recombination of the function G(phi) from
      				! E(IS).
      DOUBLE PRECISION T1	! Sum of the E(IS).

      DOUBLE PRECISION X,Y,Z
      DOUBLE PRECISION XT,B1

      INTEGER*4 OS_NM ! Maximum order of the expansion in Fourier series
			 ! of the G function of the reflection matrix computations.
			 
      INTEGER*4 LUM_NBMU ! Number of actual angles of the simulations
      
      INTEGER*4 IL	! Maximum useful order for the function E.
      INTEGER*4 I1	! Index of the angle number Teta1.
      INTEGER*4 I2	! Index of the angle number Teta2.
      INTEGER*4 IS	! Order IS of the Fourier series decomposition.
      INTEGER*4 IER  ! error code =0 if no error, =-1 otherwise

      INTEGER*4 I
      INTEGER*4 IA,IP,K,J
      
      CHARACTER*CTE_LENFIC2 FIC_RES_GSF ! Name of the file RES_GSF
                                        ! with complete directory tree



C* Opening of the result file
C--------------------------------
      OPEN(UNIT=30,FILE=FIC_RES_GSF,FORM='UNFORMATTED',ERR=991)

C* Loop on the angle Teta1 > 0.
C---------------------------------
      DO  3 I1=1 ,LUM_NBMU

C* Cosine and sine of Teta1
C---------------------------------
         C1=RMU(I1)
         S1=DSQRT(1-C1*C1)

C*   Loop on the angle teta2 > 0.
C---------------------------------
         DO 4 I2=1 ,I1

C*   Cosine and sine of Teta2
C---------------------------------
           C2=RMU(I2)
           S2=DSQRT(1-C2*C2)
           C12=C1*C2
           S12=S1*S2
           CS12=(C1+C2)
           CS12=.5*CS12*CS12

C*   Determination of G for phi = 0.
C-------------------------------------
           PHI=0.0
           CALL SOS_CALCG(CS12,C12,S12,SIG,PHI,G)
           U(0)=G
           GMAX=G

C*   Determination of G for phi = INCTE_PI.
C-----------------------------------------
           PHI=INCTE_PI
           CALL SOS_CALCG(CS12,C12,S12,SIG,PHI,G)
           U(CTE_PH_NU)=G
           GMIN=G

C*   Evaluation of the difference between Gmin and Gmax.
C*      1st case: The difference is moderate.
C*               For fixed Teta1 and Teta2, the function
C*               G(phi) is defined between 0 and PI with a step of PI/CTE_PH_NU.
C*      2nd case: The difference is large.
C*               The function G(phi) becomes negligible compared with Gmax
C*		   beyond an angle PHIB to be determined.
C*               For fixed Teta1 and Teta2, the function
C*               G(phi) is defined between 0 and PHIB with a step of PHIB/CTE_PH_NU.
C----------------------------------------------------------------------
           X=CTE_PH_TEST*GMIN
           IF(X.GE.GMAX) THEN
C* 1st case: The difference between Gmin and Gmax is moderate.
C--------------------------------------------------
              PHIB=INCTE_PI

              Q=INCTE_PI/CTE_PH_NU
              DO 6 I=1,CTE_PH_NU
                 PHI=Q*I
                 CALL SOS_CALCG(CS12,C12,S12,SIG,PHI,G)
                 U(I)=G
 6            CONTINUE

           ELSE
C* 2nd case: The difference between Gmin and Gmax is large.
C*          Search for the angle PHIB beyond which G(phi) becomes
C*          negligible compared with Gmax.
C------------------------------------------------------------------
C*      Initialisation
              PHI1=0
              PHI2=INCTE_PI
	
C*      Estimation of PHIB
C* (loop on label 9)
C---------------------------	
 9            PHIB=.5*(PHI1+PHI2)

C*      Computation of G(PHIB)
C---------------------------
              CALL SOS_CALCG(CS12,C12,S12,SIG,PHIB,G)
	
C*      Check of the difference between G(PHIB) and Gmax.
C*      Case 2A: The difference is still large.
C*              PHIB must be re-estimated.
C*              PHI1 or PHI2 is narrowed, according to the sign of the difference.
C*              PHIB is estimated again as the mean of PHI1 and PHI2
C*              and the new value G(PHIB) is computed.
C*      Case 2B: The difference is moderate.
C*              PHIB is a satisfactory estimate of the value of the
C*              angle phi up to which G(phi) remains significant
C*              compared with Gmax. The function G(phi) is described between 0 and
C*              PHIB with a step of PHIB/CTE_PH_NU.
C-------------------------------------------------------------------
              X=CTE_PH_TEST*G
              IF (DABS(X-GMAX).LT..01*GMAX) GOTO 8
	
C*      Case 2A: The difference between G(PHIB) and GMAX is still large.
C*              The value of PHIB is refined.
C-------------------------------------------------------------
              IF (X.LE.GMAX) GOTO 88
              PHI1=PHIB
              GOTO 9
 88           PHI2=PHIB
              GOTO 9
	
C*      Case 2B: The difference between G(PHIB) and GMAX is moderate.
C*              The estimate of PHIB is satisfactory.
C*              G(phi) is defined between 0 and PHIB with a step of
C*              PHIB/CTE_PH_NU.
C*              Gmin becomes G(PHIB).
C--------------------------------------------------------------------
 8            Q=PHIB/CTE_PH_NU

              DO 10 I=1,CTE_PH_NU
                 PHI=Q*I
                 CALL SOS_CALCG(CS12,C12,S12,SIG,PHI,G)
                 U(I)=G
 10           CONTINUE
              GMIN=U(CTE_PH_NU)
              T1=0.
              T2=0.
           ENDIF


C*   Fourier series decomposition.
C*   Loop on the order of the series.
C--------------------------------------------------------------------
           DO 13 IS=0,OS_NM
              Z=.5*(GMAX+GMIN*DCOS(IS*PHIB))

              IA=1
              DO 14 I=1,CTE_PH_NQ
                 IA=2*IA
                 IP=CTE_PH_NU/IA
                 Y=0
                 DO 15 J=1,IA,2
                    K=IP*J
                    Y=Y+U(K)*DCOS(IS*K*Q)
 15              CONTINUE
                 Y=2*Y/IA
	
                 XT=DABS(Z-Y)/Z
                 IF (XT.LT..0001) GOTO 16
                 Z=.5*(Y+Z)
 14           CONTINUE

 16           E(IS)=PHIB*Z/INCTE_PI
              IF (IS.NE.0) GOTO 133
              T1=E(0)
              T2=E(0)
              GOTO 13
	
 133          T1=T1+2*E(IS)
              T2=T2+2*DCOS(IS*PHIB*.5)*E(IS)

C*   Stop test of the decomposition - otherwise the loop on IS goes on.
C---------------------------------------------------------------------------
              B1=DABS(T1-GMAX)/GMAX
              IF (B1.GT..001) GOTO 13
              IL=IS
              GOTO 17	!Exit of the loop
 13        CONTINUE	!End of loop on IS
           IL=OS_NM

C*   Recording of the results
C---------------------------------
 17        WRITE(30,err=992)I1,I2,IL,(E(J),J=0,IL)

 4       CONTINUE	!End of loop on Teta2 (I2)
 
 3    CONTINUE	!End of loop on Teta1 (I1)
      
C* Nominal end
C-------------------
      IER=0
      GOTO 9999
      
C* Error case: return of status -1 to the calling program
C------------------------------------------------------
  991 WRITE(6,*) 'SOS_GSF Erreur a l''ouverture d''un fichier'
      IER=-1
      GOTO 9999
  992 WRITE(6,*) 'SOS_GSF',
     &		 'Erreur de lecture ou d''écriture sur un fichier'
      IER=-1 
      GOTO 9999

 9999 CONTINUE

C* File closing
C-------------------
      CLOSE(30)
      
      RETURN
      END         !END OF PROCEDURE SOS_GSF





C==============================================================================
C PROCEDURE: SOS_CALCG
C ==========
C      This procedure computes:
C	 G(Teta1,Teta2,phi) = exp [-tan(TetaN)**2 / SIG] / cos(TetaN)**4
C
C	 with  Teta1 : incidence zenith angle  (0 < Teta1 < PI/2),
C             Teta2 : reflection zenith angle (0 < Teta2 < PI/2),
C             phi   : relative azimuth angle between the directions
C                     of incidence and reflection,
C
C	       TetaN : zenith angle of the normal to the wave facet,
C                     computed from Teta1, Teta2 and phi.
C			
C             SIG   : square of the parameter Sigma of the Cox and Munk theory.
C
C			
C
C Description of parameters
C --------------------------
C	PHI	(double)   (I) Relative azimuth (radians).
C	C12	(double)   (I) Cos(Teta1) * Cos(Teta2).
C	S12	(double)   (I) Sin(Teta1) * Sin(Teta2).
C	CS12	(double)   (I) 0.5 * [ Cos(Teta1) + Cos(Teta2) ]**2 .
C 	SIG 	(double)   (I) Square of the parameter Sigma of the
C    			        Cox and Munk theory.
C 	G	(double)   (O) Value of G(TetaN).
C
C
C Common used:
C --------------
C     None
C
C Error cases:
C ------------
C     None
C
C==============================================================================
      SUBROUTINE SOS_CALCG(CS12,C12,S12,SIG,PHI,G)

      IMPLICIT NONE

C* Definition of variables
C*-----------------------------------------------------------------------
      DOUBLE PRECISION PHI	! Relative azimuth (radians).

      DOUBLE PRECISION C12	! Cos(Teta1) * Cos(Teta2).
      DOUBLE PRECISION S12	! Sin(Teta1) * Sin(Teta2).
      DOUBLE PRECISION CS12	! 0.5 * [ Cos(Teta1) + Cos(Teta2) ]**2 .

      DOUBLE PRECISION SIG 	! Square of the parameter Sigma of the
      				! Cox and Munk theory.

      DOUBLE PRECISION G	! Value of G(TetaN).
      
      DOUBLE PRECISION COSTETAD	! Cosine of the scattering angle between the
      				       ! directions of incidence and reflection.
      DOUBLE PRECISION X	! 1 / cos(TetaN)**2.


C* Computation
C---------------------------------
      COSTETAD=-C12+S12*DCOS(PHI)
      X=(1-COSTETAD)/CS12
      G=X*X*DEXP (-(X-1)/SIG)
      
      RETURN
      END         !END OF PROCEDURE SOS_CALCG








