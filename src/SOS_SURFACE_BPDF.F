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
C* FILE: SOS_SURFACE_BPDF.F
C* PROJECT: Successive orders of scattering
C* ROLE: Generates the BPDF files of the Rondeaux-Herman, Bréon and Maignan models.
C*
C* AUTHOR:
C*   CS GROUP France, B.Lafrance.
C* DATE: 12/08/2020
C*
C* MOD:VERSION:1.0: Creation of the program.
C*    - Reuse of existing routines: SOS_GSF_RONDEAUX_BREON, SOS_F21SF_NADAL, SOS_CALC_F21_NADAL_SUR_FRESNEL.
C*    - Development of new routines: SOS_GSF_MAIGNAN and SOS_CALCG_MAIGNAN.
C*
C* MOD:VERSION:1.1: 31/08/2022
C*    - Constants renamed from SOS_* to CTE_*
C*    - Traces converted to English.
C*
C******************************************************************************

C----------------------------------------------------------------------------
C Definition of constants
C---------------------------------------------------------------------------- 
C Constants used:
C     CTE_LENFIC2 : Length of file names with directory tree.
C    CTE_OS_NBMU_MAX : Maximum number of positive angles.
C    INCTE_PI : Value of PI.
C    CTE_OS_NB_MAX :  Maximum value of the limit order of the expansion in
C                     Legendre polynomials.
C    CTE_OS_NS_MAX :    Maximum value of the limit order of the expansion in Legendre
C                     polynomials of the Fresnel matrix elements
C                     and in Fourier series of the radiances.
C    CTE_OS_NM_MAX :    Maximum value of the limit order of the expansion
C                     in Fourier series of the G function.
C    CTE_PH_TEST : Factor for the comparison of GMIN and GMAX
C    CTE_PH_NU : Number of azimuth angles (2**NQ)
C    CTE_PH_NQ : Power of 2 of CTE_PH_NU (NU=2**NQ)
C
C----------------------------------------------------------------------------
#include "SOS.h"
#define INCTE_PI DACOS(-1.D+00)

C==============================================================================
C PROCEDURE: SOS_SURFACE_BPDF
C ==========      
C      This procedure generates the surface reflection files of the BPDF models

C          Rondeaux model:
C           => Rp(Teta1, Teta2, phi) =  G(Teta1, Teta2, phi) * Fp(Teta1, Teta2, phi) 
C                with G(Teta1, Teta2, phi) = 1 / [ 4 * (cosTeta1+cosTeta2) ]
C
C          Bréon model:
C           => Rp(Teta1, Teta2, phi) =  G(Teta1, Teta2, phi) * Fp(Teta1, Teta2, phi) 
C                with G(Teta1, Teta2, phi) = 1 / [ 4*cosTeta1*cosTeta2 ]
C
C          Maignan model:
C           => Rp(Teta1, Teta2, phi) =  G(Teta1, Teta2, phi) * Fp(Teta1, Teta2, phi) 
C                with G(Teta1, Teta2, phi) = C*exp[-tan(TetaN)] / [ 4 * (cosTeta1+cosTeta2) ]
C
C	     with Teta1 : incidence zenith angle  (0 < Teta1 < PI/2),
C                Teta2 : reflection zenith angle (0 < Teta2 < PI/2),
C             phi   : relative azimuth angle between the directions
C                        of incidence and reflection,
C                Fp(Teta1, Teta2, phi) : Fresnel reflection matrix,
C                TetaN : incidence angle with respect to the normal of the facet
C                        which, for an incidence TetaN, gives a reflection at Teta2 for the azimuth angle phi,
C                C     : constant coefficient of the Maignan model (includes the term exp(-NDVI) of the model).
C
C
C Input data of the program
C ------------------------------
C    -->Parameters/arguments of the program
C       LUM_NBMU (I4) : (I) Number of actually useful (positive) angles
C	 RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (DOUBLE) : (I) Cosines of the angles.
C	 CHR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) (DOUBLE) : (I) Weights of the angles.
C       IND   (DOUBLE)  : (I) Refractive index water / air.
C	 ISURF  (I2) : (I) Surface type ISURF
C			       (decoding format: integer I2)
C			4 :Lambertian surface + Roujean BRDF + Rondeaux BPDF.
C			5 :Lambertian surface + Roujean BRDF + Breon BPDF.
C			6 :Lambertian surface + Roujean BRDF + Nadal BPDF.
C			7 :Lambertian surface + Roujean BRDF + Maignan BPDF.
C 	 ALPHA   (DOUBLE)  (I) : Parameter alpha of the Nadal model ==> Associated with ISURF = 6
C 	 BETA    (DOUBLE)  (I) : Parameter beta of the Nadal model  ==> Associated with ISURF = 6
C       COEF_C  (DOUBLE)  (I) : Parameter C*exp(-NDVI) of the Maignan BPDF model
C                               ==> Associated with: ISURF = 7
C       OS_NB (I4) : (I) Order of the expansions of the phase functions
C				     in Legendre polynomials to be used.
C       OS_NS (I4) : (I) Maximum order to be used for the expansion in
C                                Legendre polynomials of the Fresnel matrix elements
C                                and in Fourier series of the radiance computations.
C       OS_NM (I4) : (I) Maximum order of the expansion in Fourier series
C                                of the G function of the reflection matrix computations.
C       FIC_RES_FRESNEL (CHARACTER*CTE_LENFIC2)  : 
C                 (I) Name of the temporary file RES_FRESNEL with complete directory tree.
C       FIC_RES_GSF (CHARACTER*CTE_LENFIC2)  :  
C                 (I) Name of the temporary file RES_GSF with complete directory tree.
C       FIC_RES_FRESNEL (CHARACTER*CTE_LENFIC2)  :  
C                 (I) Name of the temporary file RES_FRESNEL with complete directory tree.
C       FIC_RES_MAT_REFLEX (CHARACTER*CTE_LENFIC2)  :  
C          (I) Name of the temporary file RES_MAT_REFLEX with complete directory tree.
C       FICBPDF  (CHARACTER*CTE_LENFIC2) : (O) complete name of the result file.
C       TRACE  (LOGICAL)   ! = true if writing in the trace file.
C       IER  (I4)   (O) : error code =0 if no error, =-1 otherwise
C
C
C Result provided by the program
C --------------------------------
C   BRDF / BPDF file
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
C   -->Content of the Result file
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
C       Unformatted binary file, created by SOS_GSF_RONDEAUX_BREON, SOS_F21SF_NADAL or SOS_GSF_MAIGNAN
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
C None
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
C      - SOS_GSF_RONDEAUX_BREON
C      - SOS_F21SF_NADAL and SOS_CALC_F21_NADAL_SUR_FRESNEL
C      - SOS_GSF_MAIGNAN and SOS_CALCG_MAIGNAN
C      - SOS_MAT_FRESNEL 
C      - SOS_MAT_REFLEXION 
C      - SOS_MISE_FORMAT
C
C=============================================================================

      SUBROUTINE SOS_SURFACE_BPDF(LUM_NBMU,RMU,CHR,IND,
     &                            ISURF,ALPHA,BETA,COEF_C,
     &                            OS_NB,OS_NS,OS_NM,
     &                            FIC_RES_GSF,FIC_RES_FRESNEL,
     &                            FIC_RES_MAT_REFLEX,
     &                            FICBPDF,TRACE,IER)

      IMPLICIT NONE


C* Definition of variables
C*--------------------------------------------------------------------------
      CHARACTER*CTE_LENFIC2 FICBPDF	! Complete name of the result BPDF/BRDF file.

      CHARACTER*CTE_LENFIC2 FIC_RES_GSF ! Name of the file RES_GSF
                                        ! with complete directory tree

      CHARACTER*CTE_LENFIC2 FIC_RES_FRESNEL ! Name of the file RES_FRESNEL
                                            ! with complete directory tree

      CHARACTER*CTE_LENFIC2 FIC_RES_MAT_REFLEX ! Name of the file RES_MAT_REFLEX
                                               ! with complete directory tree

      DOUBLE PRECISION IND	! Refractive index water / air

      DOUBLE PRECISION ALPHA	!Parameter alpha of the Nadal model.
      DOUBLE PRECISION BETA	       !Parameter beta of the Nadal model.
      DOUBLE PRECISION COEF_C      !Parameter C*exp(-NDVI) of the Maignan model.
	
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   ! Array of the Gauss mu
      DOUBLE PRECISION CHR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   ! Gauss weights.
      
      LOGICAL TRACE     ! = true if writing in the trace file

      INTEGER*2 ISURF   ! Surface type

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
         WRITE(99,*,err=921) 'Calculation of the MAIGNAN BPDF file'  
	  WRITE(99,*,err=921) FICBPDF  
         WRITE(99,*,err=921) ' '
         WRITE(99,*,err=921)'Simulation parameters ' 
         WRITE(99,*,err=921)
     &        ' Refractive index surface /air ',IND
         WRITE(99,*,err=921)' OS_NB ',OS_NB
	  WRITE(99,*,err=921)' OS_NS ',OS_NS
	  WRITE(99,*,err=921)' OS_NM ',OS_NM
	  WRITE(99,*,err=921) ' ' 	 	 
      END IF
        

C* Fourier series decomposition of the function G(Teta1,Teta2,phi).
C---------------------------------------------------------------------
      IF (ISURF.EQ.4) THEN	!Rondeaux BPDF model
          CALL SOS_GSF_RONDEAUX_BREON(LUM_NBMU,RMU,1,OS_NM,
     &                                FIC_RES_GSF,TRACE,IER)
      ENDIF

      IF (ISURF.EQ.5) THEN	!Bréon BPDF model
          CALL SOS_GSF_RONDEAUX_BREON(LUM_NBMU,RMU,2,OS_NM,
     &                                FIC_RES_GSF,TRACE,IER)
      ENDIF

C			4 :Lambertian surface + Roujean BRDF + Rondeaux BPDF.
C			5 :Lambertian surface + Roujean BRDF + Breon BPDF.



      IF (ISURF.EQ.6) THEN
          CALL SOS_F21SF_NADAL(LUM_NBMU,RMU,IND,ALPHA,BETA,OS_NB,
     &                         FIC_RES_GSF,TRACE,IER)
      ENDIF

      IF (ISURF.EQ.7) THEN
          CALL SOS_GSF_MAIGNAN(LUM_NBMU,RMU,COEF_C,OS_NM,
     &                         FIC_RES_GSF,IER)
      ENDIF
     
      IF (IER.NE.0) GOTO 995   ! error case
      
      
C* Expansion of the Fresnel matrix in Legendre functions.
C------------------------------------------------------------------
      CALL SOS_MAT_FRESNEL(LUM_NBMU,RMU,CHR,IND,OS_NS,
     &                     FIC_RES_FRESNEL,TRACE,IER)
      IF (IER.NE.0) GOTO 996   ! error case
      
C* Computation of the reflection matrix expanded in Fourier series.
C------------------------------------------------------------------
      CALL SOS_MAT_REFLEXION(1.D+00,LUM_NBMU,RMU,OS_NB,OS_NS,OS_NM,
     &                       FIC_RES_FRESNEL,FIC_RES_GSF,
     &                       FIC_RES_MAT_REFLEX,IER)
      IF (IER.NE.0) GOTO 997   ! error case
      
C* Rearrangement of the results to be compatible with the OS code.
C-----------------------------------------------------------------------
      CALL SOS_MISE_FORMAT(FIC_RES_MAT_REFLEX, FICBPDF,
     &                     LUM_NBMU,OS_NB,IER)
      IF (IER.NE.0) GOTO 998   ! error case


C* Nominal end
C-------------------
      GOTO 9999
      
C* Error case and return of status 1 to the shell
C----------------------------------------------
  912 WRITE(6,*) 'Erreur à l''ouverture du fichier angles de Gauss'
      IER=-1
      GOTO 9999
  921 WRITE(6,*) 'Erreur d''écriture sur le fichier trace Glitter'
      IER=-1
      GOTO 9999
  922 WRITE(6,*) 'Erreur de lecture sur le fichier angles de Gauss'
      IER=-1
      GOTO 9999      
  995 IF(ISURF.EQ.6) THEN
         WRITE(6,*) 'Erreur dans la routine SOS_F21SF_NADAL'
      ENDIF
      IF(ISURF.EQ.7) THEN
         WRITE(6,*) 'Erreur dans la routine SOS_GSF_MAIGNAN'
      ENDIF
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
      END    !END OF PROGRAM SOS_SURFACE_BPDF




C==============================================================================
C PROCEDURE: SOS_GSF_RONDEAUX_BREON
C ==========
C      This procedure performs the Fourier series decomposition
C      (on phi) of the function G(Teta1,Teta2,phi):
C
C         Rondeaux model: G(Teta1, Teta2, phi) = 1. / [ 4 * (cosTeta1+cosTeta2) ]
C
C         Breon model:    G(Teta1, Teta2, phi) = 1. / [ 4*cosTeta1*cosTeta2 ]
C
C	    with  Teta1 : incidence zenith angle  (0 < Teta1 < PI/2),
C                Teta2 : reflection zenith angle (0 < Teta2 < PI/2),
C             phi   : relative azimuth angle between the directions
C                     of incidence and reflection,
C			
C
C      Stores the results in the file "RES_GSF" as:
C	   Number of the angle Teta1 : I1,
C         Number of the angle Teta2 : I2,
C         Maximum useful order of the decomposition: IL
C            (IL = 0 for these azimuth-independent functions),
C         Decomposed G function for (Teta1,Teta2) : (E(J),J=0,IL)
C
C Description of parameters
C --------------------------
C      LUM_NBMU (I4) : (I) Number of actually useful (positive) angles.
C      RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)	 (I) : Array of the Gauss mu.
C      IRONDEAUX (integer*4) (I) : Model index 1 <-> Rondeaux, Breon otherwise.
C      OS_NM (I4) : (I) Maximum order of the expansion in Fourier series
C                           of the G function of the reflection matrix computations.
C      FIC_RES_GSF (CHARACTER*CTE_LENFIC2)  :  
C                 (I) Name of the temporary file RES_GSF with complete directory tree.
C      TRACE (logical) (I) : Trace option.
C      IER  (I4)   (O) : error code =0 if no error, =-1 otherwise.
C
C Other result provided by the procedure
C -------------------------------------
C      Temporary file RES_GSF
C 
C Description of the files used
C ---------------------------------
C   -->Content of the temporary file RES_GSF
C       Unformatted binary file, created by SOS_GSF_RONDEAUX_BREON
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
C None
C
C Error cases:
C ------------
C     -Opening of a file
C     -Reading/writing in a file
C
C     Message displayed on screen, the routine stops its computations and
C     returns status -1 to the calling program
C
C==============================================================================
      SUBROUTINE SOS_GSF_RONDEAUX_BREON(LUM_NBMU,RMU,IRONDEAUX,
     &                                  OS_NM,FIC_RES_GSF,TRACE,IER)
 	
      IMPLICIT NONE

C* Definition of variables
C*-------------------------------------------------------------------
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! Array of the Gauss mu

      DOUBLE PRECISION E(0:CTE_OS_NM_MAX)	! Result of the decomposition in
      				        ! Fourier series of the function G.

      DOUBLE PRECISION C1	! Cos(Teta1).
      DOUBLE PRECISION S1	! Sin(Teta1).
      DOUBLE PRECISION C2	! Cos(Teta2).
      DOUBLE PRECISION S2	! Sin(Teta2).

      LOGICAL TRACE     ! = true if writing in the trace file.

      INTEGER*4 OS_NM   ! Maximum order of the expansion in Fourier series
			   ! of the G function of the reflection matrix computations.
      INTEGER*4 LUM_NBMU  ! Number of actual angles of the simulations
      INTEGER*4 IRONDEAUX ! 1 if simulation of the Rondeaux model.
      		            ! Otherwise simulation of the Breon model.
      INTEGER*4 IL	! Maximum useful order for the function E.
      INTEGER*4 I1	! Index of the angle number Teta1.
      INTEGER*4 I2	! Index of the angle number Teta2.
      INTEGER*4 IER  ! error code =0 if no error, =-1 otherwise
      INTEGER*4 J	! Array index.
      
      CHARACTER*CTE_LENFIC2 FIC_RES_GSF ! Name of the file RES_GSF
                                        ! with complete directory tree


C* Opening of the result file
C--------------------------------
      OPEN(UNIT=30,FILE=FIC_RES_GSF,FORM='UNFORMATTED',ERR=991)

C* Trace file
C--------------------------------
      IF (TRACE) THEN
         IF(IRONDEAUX.EQ.1) THEN
	     WRITE(99,*,err=921) 'Values of the function G Rondeaux'     
	  ELSE
	     WRITE(99,*,err=921) 'Values of the function G Breon'       
         ENDIF
	    	 
      ENDIF  
	   
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
          

C*   G function at order IS = 0 for the Rondeaux model.
C*   => Independent of the azimuth: single term IS = 0.
C---------------------------------
           IF(IRONDEAUX.EQ.1) THEN
	      E(0) = 1. / (1./C1 + 1./C2)
	   ELSE 
	      E(0) = 1.
	   ENDIF   
	   IL=0

C* Writing of the results in the trace file
C--------------------------------
	   IF (TRACE) THEN
	       WRITE(99,1000,err=921) ACOS(C1)*180./INCTE_PI, 
     &	                              ACOS(C2)*180./INCTE_PI,
     &                               E(0)         	 
           ENDIF  
	      
                 

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
  921 WRITE(6,*) 'Erreur d''écriture sur le fichier trace Surface'
      IER=-1
      GOTO 9999
  991 WRITE(6,*) 'SOS_GSF_RONDEAUX_BREON: ',
     &           'Erreur a l''ouverture d''un fichier'
      IER=-1
      GOTO 9999
  992 WRITE(6,*) 'SOS_GSF_RONDEAUX_BREON: ',
     &		 'Erreur de lecture ou d''écriture sur un fichier'
      IER=-1 
      GOTO 9999

 9999 CONTINUE

C* File closing
C-------------------
      CLOSE(30)

C* Format
C--------
1000  FORMAT(6hTeta1:,1X,F5.2,2X,6hTeta2:,1X,F5.2,2X,
     &2hG:,1X,E11.4)    
           
      RETURN
      END         !END OF PROCEDURE SOS_GSF_RONDEAUX_BREON



C==============================================================================
C PROCEDURE: SOS_F21SF_NADAL
C ==========
C      This procedure performs the Fourier series decomposition
C      (on phi) of the function F(Teta1,Teta2,phi):
C
C	 F(Teta1,Teta2,phi) = F12nadal(Teta1,Teta2,phi) 
C                            / F12fresnel(Teta1,Teta2,phi)
C	 with Teta1 : incidence zenith angle: tetas  (0 < Teta1 < PI/2),
C             Teta2 : reflection zenith angle: tetav (0 < Teta2 < PI/2),
C             phi   : relative azimuth angle between the directions
C                     of incidence and reflection with the convention
C                     of the OS (0 <= phi <= PI),
C			
C             F12fresnel = 0.5 * (RR*RR - RL*RL) 
C             with the sign convention of the Nadal model.
C
C                with RL the parallel reflection coefficient,
C                computed from the refractive index of the surface
C                and the reflection angle teta:
C                RL(teta) = [ IND**2 * cos(teta) - SQRT(IND**2 - sin(teta)**2) ]
C                         / [ IND**2 * cos(teta) + SQRT(IND**2 - sin(teta)**2) ]
C
C                RR the perpendicular reflection coefficient:
C                RR(teta) = [ cos(teta) - SQRT(IND**2 - sin(teta)**2) ]
C                         / [ cos(teta) + SQRT(IND**2 - sin(teta)**2) ]
C
C                
C
C             F12nadal = alpha * 
C                  [ 1 -exp(-beta*F12fresnel / (cosTeta1 + cosTeta2) ) ]
C
C                 with alpha and beta the parameters of the Nadal model.
C
C      The BPDF of the Nadal model is converted from polarised reflectance
C      into a quantity useful for the OS: normalised radiance * cos(Teta2),
C      i.e.: BPDF * cos(Teta1) * cos(Teta2).
C
C
C      The procedure stores the results in the file "RES_MAT_REFLEX"
C      as:
C         Decomposed function F for (Teta1,Teta2) : (E(IS),IS=0,OS_NB)
C
C Description of parameters
C --------------------------
C      LUM_NBMU (I4) (I) : Number of actually useful (positive) angles.
C      RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)	 (I) Array of the Gauss mu.
C      IND     (double)  (I) : Refractive index of the surface.
C      ALPHA   (double)  (I) : Parameter alpha of the Nadal model.
C      BETA    (double)  (I) : Parameter beta of the Nadal model.
C      OS_NB (I4) (I) : Order of the expansions of the phase functions
C		           in Legendre polynomials to be used.
C      TRACE    (logical) (I) = true if writing in the trace file
C       IER  (I4)   (O) : error code =0 if no error, =-1 otherwise
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
C       Fourier series decomposition of the function F21 of the Nadal
C       model.
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
C
C     Message displayed on screen, the routine stops its computations and
C     returns status -1 to the calling program
C
C     
C  Routines used by the procedure:
C ------------------------------------
C  This program calls the routines:
C        - SOS_CALC_F21_NADAL_SUR_FRESNEL
C
C==============================================================================
      SUBROUTINE SOS_F21SF_NADAL(LUM_NBMU,RMU,IND,ALPHA,BETA,
     &                           OS_NB,FIC_RES_GSF,TRACE,IER)
 	
      IMPLICIT NONE

C* Definition of variables
C*-------------------------------------------------------------------
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! Array of the Gauss mu

      DOUBLE PRECISION IND	 ! Refractive index surface / air
      DOUBLE PRECISION ALPHA !Parameter alpha of the Nadal model.
      DOUBLE PRECISION BETA	 !Parameter beta of the Nadal model.
      
      DOUBLE PRECISION E(0:CTE_OS_NB_MAX) ! Result of the decomposition in
      				              ! Fourier series of the function F.
	
      DOUBLE PRECISION PHI	! Relative azimuth for the OS convention (radians).

      DOUBLE PRECISION C1	! Cos(Teta1).
      DOUBLE PRECISION S1	! Sin(Teta1).
      DOUBLE PRECISION C2	! Cos(Teta2).
      DOUBLE PRECISION S2	! Sin(Teta2).

      DOUBLE PRECISION F	! Value of the function F(phi).
      DOUBLE PRECISION FMIN	! Minimum value kept for the function
      				! F(phi) for fixed Theta1 and Theta2.
      DOUBLE PRECISION FMAX	! Maximum value of the function F(phi).

      DOUBLE PRECISION U(0:CTE_PH_NU)  ! Function F(phi) between 0 and PI
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
      				    ! and the recombined one
				    !=> for all the pairs (Teta1, Teta2).
      
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
      
      LOGICAL TRACE

      INTEGER*4 OS_NB ! Order of the expansions of the phase functions
			 ! in Legendre polynomials.
				
      INTEGER*4 LUM_NBMU ! Number of actual angles of the simulations
      	    
      INTEGER*4 J	! Array index.
      INTEGER*4 IL	! Maximum useful order for the function E.
      INTEGER*4 I1	! Index of the angle number Teta1.
      INTEGER*4 I2	! Index of the angle number Teta2.
      INTEGER*4 IS	! Order IS of the Fourier series decomposition.
      INTEGER*4 IS2	! Order IS of the Fourier series decomposition.
      INTEGER*4 IER  ! error code =0 if no error, =-1 otherwise
      INTEGER*4 I
      INTEGER*4 COMPT	! Number of configurations for the sum of the relative differences.
     
      CHARACTER*CTE_LENFIC2 FIC_RES_GSF ! Name of the file RES_GSF
                                        ! with complete directory tree
	
C* Initialisation
C---------------------------------------------------------	
      ECARTMAX = -1.D+300
      SOMECART = 0.
      COMPT = 0
      
     
C* Creation of the result file of SOS_MAT_REFLEXION
C---------------------------------------------------------
      OPEN(UNIT=30,FILE=FIC_RES_GSF,FORM='UNFORMATTED',ERR=991)

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
	
              PHI=Q*I
              CALL SOS_CALC_F21_NADAL_SUR_FRESNEL(IND,ALPHA,BETA,
     &	                                           C1,S1,C2,S2,PHI,F)   
              U(I)=F
	      
 6         CONTINUE

		
C*   Determination of the min and max values of F
C*   for fixed I1 and I2 but varying PHI.
C----------------------------------------------
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
           ENDIF  !trace
	      
	 
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
	         PHI = I*Q		 		!Azimuth convention OS
	         Y=Y+U(I)*DCOS(IS*PHI)
	      ENDDO
	      E(IS) = Y*Q / INCTE_PI
	      
C--------- Determination of the maximum difference between exact and recombined F: B1
C---------------------------------------------------------------------      
	      B1=0.D+00			!initialisation
	      DO I=0,CTE_PH_NU
	         PHI = I*Q		 		!Azimuth convention OS
	         CALL SOS_CALC_F21_NADAL_SUR_FRESNEL(IND,ALPHA,BETA,
     &	                                          C1,S1,C2,S2,PHI,F)

C-------------- Computation of the recombined F from the function E(IS) for fixed PHI.
                T1=E(0)	
		  DO IS2 = 1, IS		  
	              T1=T1+2.*E(IS2)*DCOS(IS2*PHI)	!Azimuth convention OS          		     
		  ENDDO !End of loop on IS2
       
                B1=MAX( DABS((T1-F)/F) , B1)
		 		 
	      ENDDO !End of loop on I (phi = 0 to PI)
	           	     

C--------- Stop test of the decomposition:
C---------    * Stop if the maximum relative difference between exact and recombined F
C---------      below the threshold CTE_SEUIL_SF_NADAL
C---------      or if it increases with respect to the previous order.
C---------    * Otherwise the loop on IS goes on.
C-------------------------------------------------------------------------------- 		
              IF (B1.GT.CTE_SEUIL_SF_NADAL) GOTO 14
              IL=IS	      
              GOTO 17	!Exit of the loop on IS
	      
 14           IF(B1.LT.B1_PREC) GOTO 141
              IL=IS-1
	      GOTO 17	!Exit of the loop on IS
	      
 141	      B1_PREC=B1	        
 13        CONTINUE	!End of loop on IS
	   IL=OS_NB



C* Recording for the pair (Teta1, Teta2) of the decomposition
C* in Fourier series of the BPDF (converted into normalised radiance
C* and adapted to the OS code).
C-----------------------------------------------------------------
 17        WRITE(30,err=992)I1,I2,IL,(E(J),J=0,IL)


C* Trace file: Computation of quality analysis parameters
C*                 of the Fourier series decomposition
C*                 and recordings.
C-----------------------------------------------------------------
           IF (TRACE) THEN

	      ECARTMAX_C1C2=0.
	            	      
	      Q=INCTE_PI/CTE_PH_NU
	      DO I=0,CTE_PH_NU
	         PHI=I*Q
	         CALL SOS_CALC_F21_NADAL_SUR_FRESNEL(IND,ALPHA,BETA,
     &	                                             C1,S1,C2,S2,PHI,F)   
		 
                FRECOMB = E(0)
	         DO IS=1,IL
	            FRECOMB = FRECOMB + 2.*DCOS(IS*PHI)*E(IS)
	         ENDDO
            
	         IF(F.NE.0.D+00) THEN
	            ECART = DABS((FRECOMB-F)/F)  
		    
C*       Analysis for all the azimuths of the current pair (Teta1, Teta2)
		     IF(ECART.GT.ECARTMAX_C1C2) THEN
		        ECARTMAX_C1C2 = ECART
		        PHIECARTMAX_C1C2 = PHI*180./INCTE_PI	!Azimuth convention OS
	               F_ECARTMAX_C1C2 = F
	               FRECOMB_ECARTMAX_C1C2 = FRECOMB
		     ENDIF 
		    
C*       Analysis for all the azimuths and all the pairs (Teta1, Teta2)
	            IF(ECART.GT.ECARTMAX) THEN
	               ECARTMAX = ECART
	               PHIECARTMAX = PHI*180./INCTE_PI	!Azimuth convention OS
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
     &					    F_ECARTMAX_C1C2,
     &		                         FRECOMB_ECARTMAX_C1C2,
     &					    ECARTMAX_C1C2*100.
	      ELSE
	      	  WRITE(99,2500,err=921) IL,PHIECARTMAX_C1C2,
     &				 	    F_ECARTMAX_C1C2,
     &		                         FRECOMB_ECARTMAX_C1C2
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
     &	  "Analysis of the Fourier expansion of the NADAL's BPDF" 
	
	  IF(COMPT.GT.1) THEN
	     WRITE(99,*,err=921)
	     WRITE(99,*,err=921)
	     WRITE(99,*,err=921) 
     &	          'Mean relative difference from the recombination'//
     &	          ' (% excluding null values):', SOMECART*100/COMPT
	     WRITE(99,*,err=921)
	     WRITE(99,*,err=921)   
     &	          'Max relative difference from the recombination (%):',
     &	          ECARTMAX*100
            WRITE(99,*,err=921)   
     &	          'Related delta phi angle (SOS convention) :', 
     &	          PHIECARTMAX
            WRITE(99,*,err=921)    
     &	          'Related viewing angle :',TVECARTMAX
	     WRITE(99,*,err=921)    
     &	          'Related solar angle :',TSECARTMAX
	     WRITE(99,*,err=921)   
     &	          'Related value of BRDF :',F_ECARTMAX
	     WRITE(99,*,err=921)
     &	          'Related value of the recombined BRDF :',
     &	          FRECOMB_ECARTMAX
         ELSE
	     WRITE(99,*,err=921)
	     WRITE(99,*,err=921) 'WARNING : nulle BPDF !'
	     WRITE(99,*,err=921) 'Simulation parameters to be checked'
	 ENDIF 
	 
	 WRITE(99,*,err=921)
	 WRITE(99,*,err=921)
	 
      END IF 
      
     
      
C* Nominal end
C-------------------
      IER=0
      goto 9999
      
C* Error case: return of status -1 to the calling program
C------------------------------------------------------
  921 WRITE(6,*) 'Erreur d''écriture sur le fichier trace Surface'
      IER=-1
      goto 9999
  991 WRITE(6,*) 'SOS_F21SF_NADAL : Erreur a l''ouverture d''un fichier'
      IER=-1
      goto 9999
  992 WRITE(6,*) 'SOS_F21SF_NADAL',
     &		 'Erreur de lecture ou d''écriture sur un fichier'
      IER=-1 
      goto 9999
      
 9999 CONTINUE


C* File closing: 'RES_MAT_REFLEX'
C-------------------------------------
      CLOSE(30)
 
 
C* Format
C--------

1000  FORMAT(6hTeta1:,1X,F5.2,2X,6hTeta2:,1X,F5.2,2X,
     &15h; (PHIos,Fmin):,1X,F6.2,2X,E11.4,2X,
     &15h; (PHIos,Fmax):,1X,F6.2,2X,E11.4)    
      
2000  FORMAT(5X,12hLimit Order:,I3,3h ; ,14h Error max => ,1X,4hphi:,
     &1X,F6.2,3h, F:,1X,E11.4,2X,11h, F recomb:,1X,E11.4,2X,
     &12h, Ecart (%):,1X,F7.3)         
     
2500  FORMAT(5X,12hLimit Order:,I3,3h ; ,14h Error max => ,1X,4hphi:,
     &1X,F6.2,3h, F:,1X,E11.4,2X,11h, F recomb:,1X,E11.4)         
                           
      RETURN
      END         !END OF PROCEDURE SOS_F21SF_NADAL



C==============================================================================
C PROCEDURE: SOS_CALC_F21_NADAL_SUR_FRESNEL
C		
C      This procedure computes the function G(Teta1, Teta2, phi) for the
C      Nadal model:
C
C      G(Teta1, Teta2, phi) 
C        =  F21Nadal(Teta1, Teta2, phi) / F21Fresnel(Teta1, Teta2, phi)
C
C        with F21Nadal(Teta1, Teta2, phi)
C             = Alpha * [ 1 -exp(-Beta*F21Fresnel / ( cosTeta1+cosTeta2) ) ]
C
C        and F21Fresnel(Teta1, Teta2, phi) = 0.5 * (RR*RR - RL*RL)
C        with the sign convention of the Nadal model.
C
C            RL is the parallel reflection coefficient
C
C               RL(teta) = [ IND**2 * cos(teta) - SQRT(IND**2 - sin(teta)**2) ]
C                        / [ IND**2 * cos(teta) + SQRT(IND**2 - sin(teta)**2) ]
C
C            RR is the perpendicular reflection coefficient
C
C               RR(teta) = [ cos(teta) - SQRT(IND**2 - sin(teta)**2) ]
C                        / [ cos(teta) + SQRT(IND**2 - sin(teta)**2) ]
C
C            computed from the refractive index of the surface
C            and the reflection angle teta.
C
C
C
C      The BPDF of the Nadal model is converted from polarised reflectance
C      into a quantity useful for the OS: normalised radiance * cos(Teta2),
C      i.e.: BPDF * cos(Teta1) * cos(Teta2).
C
C Description of parameters
C --------------------------
C       IND	(double)   (I) Refractive index surface / air.
C       ALPHA	(double)   (I) Parameter alpha of the Nadal model.
C       BETA	(double)   (I) Parameter beta of the Nadal model.
C	 C1	(double)   (I) Cos(Teta1).
C	 S1	(double)   (I) Sin(Teta1).
C	 C2	(double)   (I) Cos(Teta2).
C	 S2	(double)   (I) Sin(Teta2).
C       PHI   (double)   (I) Relative azimuth for the OS convention (radians).
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
      SUBROUTINE SOS_CALC_F21_NADAL_SUR_FRESNEL(IND,ALPHA,BETA,
     &	                                        C1,S1,C2,S2,PHI,F)   

      IMPLICIT NONE

C* Definition of variables
C*-----------------------------------------------------------------------
      DOUBLE PRECISION IND	! (I)  Refractive index surface / air
      DOUBLE PRECISION ALPHA! (I)  Parameter alpha of the Nadal model.
      DOUBLE PRECISION BETA	! (I)  Parameter beta of the Nadal model.
      DOUBLE PRECISION C1	! (I)  Cos(Teta1).
      DOUBLE PRECISION S1	! (I)  Sin(Teta1).
      DOUBLE PRECISION C2	! (I)  Cos(Teta2).
      DOUBLE PRECISION S2	! (I)  Sin(Teta2).
      DOUBLE PRECISION PHI	! (I) Relative azimuth for the convention
      				!      OS (radians).
            
      DOUBLE PRECISION F	! (O) Value of the function F21nadal/F21fresnel
      				! for the observation geometry
				! (converted into normalised radiance * cos(teta2).
      
      DOUBLE PRECISION F21FRESNEL	! Value of the Fresnel F21 function.
      DOUBLE PRECISION F21NADAL	! Value of the Nadal F21 function.
      
      DOUBLE PRECISION COEF

      DOUBLE PRECISION COSDIF	! Cosine of the scattering angle.
      DOUBLE PRECISION COSW	! Cosine of the reflection angle W
      				! (incidence angle with respect to the
				!  normal of the interface)
      
      DOUBLE PRECISION V	! Square sine of W.
      DOUBLE PRECISION X	! 

      DOUBLE PRECISION IND2	! Square of the refractive index surface / air
      DOUBLE PRECISION RL	! Parallel reflection coefficient.
      DOUBLE PRECISION RR	! Perpendicular reflection coefficient.

      DOUBLE PRECISION PI	! Value of PI
      

C* Processing
C---------------------------------    
      PI=4.D+00*DATAN(1.D+00)

      COEF = 4.*C1*PI
	
C*   Computation of the cosine of the scattering angle
C---------------------------------------------------------     	
      COSDIF=-C1*C2+DSQRT(1.-C1*C1)*DSQRT(1.-C2*C2)*DCOS(PHI)

C*   Computation of the cosine of the reflection angle
C---------------------------------------------------------   
      COSW=DSQRT(.5*(1-COSDIF))
      
C*   Computation of the square sine of the reflection angle
C--------------------------------------------------------- 
      V=.5*(1+COSDIF)  
     
C*   Computation of the parallel and perpendicular reflection
C*   coefficients with the Fresnel theory
C---------------------------------------------------------         
      IND2=IND*IND         
      X=DSQRT(IND2-V)     
      RL=(IND2*COSW-X)/(IND2*COSW+X)
      RR=(COSW-X)/(COSW+X)
    

C*   Computation of the Fresnel F21 (Nadal sign convention)
C--------------------------------------------------------- 
      F21FRESNEL = 0.5D+00 * (RR*RR-RL*RL)      
      
C*   Computation of the Nadal F21
C--------------------------------------------------------- 
      F21NADAL = -BETA * F21FRESNEL / (C1 + C2)     
      F21NADAL = ALPHA * ( 1.D+00 - DEXP(F21NADAL) )

C*   Computation of F
C---------------------------------------------------------            
      IF(F21FRESNEL.LT.1.0D-10) THEN  !Limited expansion in the neighbourhood
      					   !of 0 for F21Fresnel
         F = ALPHA*BETA / (C1 + C2)
      ELSE	 
         F = F21NADAL / F21FRESNEL
      ENDIF


C*   Conversion reflectance -> normalised radiance (*C1)
C*   and adaptation to the OS -> *C2.
C--------------------------------------------------------- 
      F = F*COEF*C2*C1		

      RETURN
      
      END         !END OF PROCEDURE SOS_CALC_F21_NADAL_SUR_FRESNEL







C==============================================================================
C PROCEDURE: SOS_GSF_MAIGNAN
C ==========
C
C      This procedure performs the Fourier series decomposition
C      (on phi) of the function G(Teta1,Teta2,phi):
C
C	 G(Teta1,Teta2,phi) = SIG * exp [-tan(TetaN)] / [ 4 * (cosTeta1+cosTeta2) ]
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
C       OS_NM (I4) : (I) Maximum order of the expansion in Fourier series
C                                of the G function of the reflection matrix computations.
C      FIC_RES_GSF  (CHARACTER)  : (I) Name of the temporary file RES_GSF
C                                            with complete directory tree.
C
C     IER (I4) : (O) error code =0 if no error, =-1 otherwise
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
C None
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
C        - SOS_CALCG_MAIGNAN
C
C==============================================================================
      SUBROUTINE SOS_GSF_MAIGNAN(LUM_NBMU,RMU,COEF_C,OS_NM,
     &                           FIC_RES_GSF,IER)
      IMPLICIT NONE

C* Definition of variables
C*-------------------------------------------------------------------
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX) ! Array of the Gauss mu

      DOUBLE PRECISION COEF_C	! Coefficient C of the Maignan BPDF model

      DOUBLE PRECISION E(0:CTE_OS_NM_MAX) ! Result of the decomposition in
      				              ! Fourier series of the function G.

      DOUBLE PRECISION PHI	! Relative azimuth (radians).

      DOUBLE PRECISION G	! Value of the function G(phi).

      DOUBLE PRECISION C1	! Cos(Teta1).
      DOUBLE PRECISION S1	! Sin(Teta1).
      DOUBLE PRECISION C2	! Cos(Teta2).
      DOUBLE PRECISION S2	! Sin(Teta2).
      DOUBLE PRECISION S12	! Sin(Teta1) * Sin(Teta2).

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

      INTEGER*4 OS_NM  ! Maximum order of the expansion in Fourier series
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
           S12=S1*S2

C*   Determination of G for phi = 0.
C-------------------------------------
           PHI=0.0
           CALL SOS_CALCG_MAIGNAN(C1,C2,S12,PHI,COEF_C,G)
           U(0)=G
           GMAX=G

C*   Determination of G for phi = INCTE_PI.
C-----------------------------------------
           PHI=INCTE_PI
           CALL SOS_CALCG_MAIGNAN(C1,C2,S12,PHI,COEF_C,G)
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
                 CALL SOS_CALCG_MAIGNAN(C1,C2,S12,PHI,COEF_C,G)
                 U(I)=G
 6            CONTINUE

           ELSE
C* 2nd case: The difference between Gmin and Gmax is large.
C*          Search for the angle PHIB beyond which G(phi) becomes
C*          negligible compared with Gmax.
C------------------------------------------------------------------
C* Initialisation
              PHI1=0
              PHI2=INCTE_PI
	
C*      Estimation of PHIB
C* (loop on label 9)
C---------------------------	
 9            PHIB=.5*(PHI1+PHI2)

C*      Computation of G(PHIB)
C---------------------------
              CALL SOS_CALCG_MAIGNAN(C1,C2,S12,PHIB,COEF_C,G)
	
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
                 CALL SOS_CALCG_MAIGNAN(C1,C2,S12,PHI,COEF_C,G)
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
      END         !END OF PROCEDURE SOS_GSF_MAIGNAN





C==============================================================================
C PROCEDURE: SOS_CALCG_MAIGNAN
C ==========
C       The weighting function of the Fresnel matrix for the Maignan
C       BPDF model is:
C	 G(Teta1,Teta2,phi) = C * exp [-tan(TetaI)] / [ 4 * (cosTeta1+cosTeta2) ]
C
C	 with  Teta1 : incidence zenith angle  (0 < Teta1 < PI/2),
C            Teta2 : reflection zenith angle (0 < Teta2 < PI/2),
C             phi   : relative azimuth angle between the directions
C                    of incidence and reflection,
C
C	      TetaI : incidence angle with respect to the normal of the surface facet,
C                    computed from Teta1, Teta2 and phi.
C
C	      C : site coefficient (independent of the observation geometry).
C
C      This procedure computes: G*4*cosTeta1*cosTeta2
C      i.e. C * exp [-tan(TetaI)] / [ 1/cosTeta1 + 1/cosTeta2) ]
C
C Description of parameters
C --------------------------
C	 C1	(double)   (I) Cos(Teta1).
C	 C2	(double)   (I) Cos(Teta2).
C	S12	(double)   (I) Sin(Teta1) * Sin(Teta2).
C	PHI	(double)   (I) Relative azimuth (radians).
C	COEF_C	(double)   (I) Coefficient C of the Maignan model.
C 	G	(double)   (O) Value of G(TetaN).
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
      SUBROUTINE SOS_CALCG_MAIGNAN(C1,C2,S12,PHI,COEF_C,G)

      IMPLICIT NONE

C* Definition of variables
C*-----------------------------------------------------------------------
      DOUBLE PRECISION PHI	! Relative azimuth (radians).

      DOUBLE PRECISION C1	! Cos(Teta1).
      DOUBLE PRECISION C2	! Cos(Teta2).
      DOUBLE PRECISION S12	! Sin(Teta1) * Sin(Teta2).      
      DOUBLE PRECISION COS_2I  ! Cosine of twice the incidence angle.
      DOUBLE PRECISION TAN2_I  ! Square of the tangent of the incidence angle.
      DOUBLE PRECISION COEF_C ! Coefficient C of the Maignan BPDF model
      DOUBLE PRECISION G	! Value of G(Teta1,Teta2,phi).

C* Computation
C---------------------------------
      COS_2I = C1*C2-S12*DCOS(PHI)
      TAN2_I = (1-COS_2I) / (1+COS_2I)

      !Management of the case of numerical rounding errors (negative square tangent impossible).
      IF (TAN2_I.LT.0.D+00) TAN2_I=0.D+00

      G= COEF_C *DEXP(-DSQRT(TAN2_I)) / (1./C1 + 1./C2)

      RETURN
      END         !END OF PROCEDURE SOS_CALCG_MAIGNAN








