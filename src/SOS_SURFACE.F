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
C* FILE: SOS_SURFACE.F
C* PROJECT: Successive orders of scattering
C* ROLE: This module generates the reflection matrices for different types of
C*	 models.
C*       
C*
C* AUTHOR:
C*   Algorithmic basis: Laboratoire d'Optique Atmospherique (LOA, Villeneuve d'Ascq).
C*   Code rework and development: CS GROUP France, B.Lafrance, X.Lenot
C*
C* DATE: creation 1.0 : 30/04/2002
C*
C* MOD:VERSION:1.1: Upgrade of the LOA codes (codes related to the Glitter).
C*                  Developments to adapt the codes to the other models.
C* MOD:VERSION:2.0: Correction of an error on the surface polarisation
C*                  for routine SOS_MAT_REFLEXION.
C* MOD:VERSION:2.1: 07/07/2008 
C*                  Change of the reading format of the parameter WIND,
C*                  to correct a decoding error for integer values.
C*			WIND changed from F4.1 to F4
C* MOD:VERSION:3.0: 17/11/2009 
C*                  * Change of the argument passing:
C*       - passing by pairs (keyword, value).
C*       - adaptation of the error case management.
C*                  * Error messages converted to English.
C* MOD:VERSION:3.1: 22/12/2009 
C*                  Change of the reading format of the parameters Alpha and Beta of the Nadal model
C*                  to correct a decoding error for integer values.
C*			ALPHA changed from F6.4 to F6
C*			BETA  changed from F4.1 to F4
C*                  
C* MOD:VERSION:4.0: 10/01/2010 
C*
C*   - Removal of the constant SOS_OS_FICANGLE (from SOS.h)
C*     in favour of passing the name of the angle file (produced by SOS_ANGLES) as argument.
C*
C*   - Change of the size of the radiance arrays:
C*     ==> Changed from CTE_OS_NBMU to CTE_OS_NBMU_MAX.
C*     ==> Loop on the angles up to LUM_NBMU, equal to the actual number of angles used.
C*
C*   - Removal of the constants CTE_OS_NB, CTE_OS_NS and CTE_OS_NM (from SOS.h)
C*     in favour of passing the values OS_NB, OS_NS, OS_NM (read in the angle file) as arguments
C*
C*   - Change of the size of the arrays of the decomposition in Legendre polynomials:
C*     ==> Changed from CTE_OS_NB to CTE_OS_NB_MAX.
C*     ==> Changed from CTE_OS_NS to CTE_OS_NS_MAX.
C*     ==> Changed from CTE_OS_NM to CTE_OS_NM_MAX.
C*     ==> Loop on the orders up to OS_NB, OS_NS or OS_NM
C*
C*   - Routines SOS_MISE_FORMAT & SOS_BPDF_AJOUT_BRDF: change of the
C*     record or reading line of the arrays of the reflection
C*     matrix at order IS of the Fourier series decomposition:
C*     ==> limited to the useful bounds -LUM_NBMU:LUM_NBMU
C*       Rem: no change of the record bounds of the temporary files
C*             RES_FRESNEL, SOS_MAT_REFLEX since they are not
C*             final products.
C*
C*   - Adaptation of the call of the file naming routines
C*     (SOS_NOM_FICGLITTER, SOS_NOM_FICROUJ, SOS_NOM_FICRH, SOS_NOM_FICBREON, SOS_NOM_FICNADAL)
C*     to the new parameters NBMU_GAUSS and FICANGLES_USER
C*
C*   - Adaptation of routine SOS_BPDF_AJOUT_BRDF to the passing of the parameter
C*     LUM_NBMU (nb of useful angles) and OS_NB
C*
C*   - Adaptation of the call of routines SOS_GLITTER, SOS_ROUJEAN, SOS_RONDEAUX_BREON,
C*     SOS_NADAL and SOS_BPDF_AJOUT_BRDF to the passing of the parameters
C*     LUM_NBMU (nb of useful angles), RMU (table of the cosines of the angles)
C*     and CHR (table of the weights of the angles), as well as OS_NB, OS_NS and OS_NM
C*
C* MOD:VERSION:4.1: 24/06/2010 
C*     FORMAT labels ,X, changed to ,1X, to avoid compilation errors.
C*
C* MOD:VERSION:4.2: 06/06/2013
C*
C*   -  The input parameters are given as a list of parameters
C*      and no longer as pairs with associated keyword => all
C*	 parameters are therefore mandatory
C*
C*   -  Change of the value of Pi, declared as a constant in the file
C*      (INCTE_PI replaces CTE_PI)
C*  
C* MOD:VERSION:4.3: 24/02/2016 
C*    - Minor adjustments to strictly respect the 72 columns
C*
C*    - Change of the way the formatted files are read (Angles, ...)
C*      to be compatible with the gfortran compiler.
C* 
C*  
C* MOD:VERSION:5.0: 12/08/2020
C*   - Added routines for the computation of the Maignan BPDF matrix (with Roujean BRDF):
C*     --> Call of routine SOS_NOM_FICMAIGNAN (defined in the file SOS_NOM_FIC_SURFACE.F)
C*     --> Introduction of the new input parameter COEF_C_MAIGNAN, used in the Maignan model.
C*
C*   - Call of the new routine SOS_SURFACE_BPDF instead of SOS_RONDEAUX_BREON and SOS_NADAL (removed)
C*     to manage the computation of the surface BRDF/BPDF matrices (Rondeaux, Bréon, Nadal and Maignan models).
C*   
C*   - Introduction of a directory DIR_TMP for the location of the temporary files.
C*
C*   -  Names of the temporary files (RES_GSF, RES_FRESNEL, RES_MAT_REFLEX, FICTMP_BPDF) passed,
C*      with complete directory tree, to the subroutines SOS_GLITTER, SOS_ROUJEAN, SOS_SURFACE_BPDF,
C*      SOS_MAT_FRESNEL and SOS_MAT_REFLEXION.
C*
C*   - Change of the thresholds for the tests of limit values of the coefficients KO, K1 and K2
C*     of the Roujean model: changed from 1000 to 10
C*     (some cases with K2 larger than 1 in the paper of ROUJEAN ET AL (A Bidirectional Reflectance Model
C*      of the Earth's Surface for the Correction of Remote Sensing Data, JOURNAL OF GEOPHYSICAL RESEARCH,1992).
C*
C* MOD:VERSION:5.1: 31/08/2022
C*    - Constants renamed from SOS_* to CTE_*
C*    - Traces converted to English.
C*
C* MOD:VERSION:5.2: 26/09/2023
C*    - Removal of a comma after the command READ(...) in two places.
C*      No impact on the simulation results, but raised a warning at compilation.
C*    - Syntax changes for compatibility with Fortran 2018 for the gcc/10.2.0 compiler:
C*         - Zero-completion loop in routine SOS_MAT_REFLEXION
C*         - Initialisation loops for order IS >= 2 in routine SOS_NOYAUX_FRESNEL
C*         - Loops on Teta1 and Teta2 in routine SOS_NOYAUX_FRESNEL
C*
C*
C******************************************************************************

C----------------------------------------------------------------------------
C Definition of constants
C---------------------------------------------------------------------------- 
C Constants used:
C    CTE_LENDIR : Length of directory names.
C     CTE_LENFIC2 : Length of file names with directory tree.
C    CTE_LENKEYWORD : Maximum size of the keywords for argument passing
C    INCTE_PI : Value of PI
C    CTE_OS_NBMU_MAX : Maximum number of positive angles.
C    CTE_OS_NB_MAX :  Maximum value of the limit order of the expansion in
C                     Legendre polynomials.
C    CTE_OS_NS_MAX :    Maximum value of the limit order of the expansion in Legendre
C                     polynomials of the Fresnel matrix elements
C                     and in Fourier series of the radiances.
C    CTE_OS_NM_MAX :    Maximum value of the limit order of the expansion
C                     in Fourier series of the G function.
C
C Specific constants:
C    LENLINE  : Maximum length of a data line read in a file
C----------------------------------------------------------------------------
#include "SOS.h"
#define INCTE_PI DACOS(-1.D+00)
#define LENLINE 1000

C==============================================================================
C PROGRAM: SOS_SURFACE
C ==========  
C      This program computes the reflection matrices (BRDF and/or BPDF)
C      of different models:
C	- Glitter reflection,
C	- Roujean BRDF model to simulate the reflection in intensity
C        of vegetation covers,
C	- Simplified BPDF model of Rondeaux and Herman (RH) to simulate the
C	  reflection in polarisation of vegetation covers,
C	- Simplified BPDF model of Breon, to simulate the reflection
C        in polarisation of bare soils,
C	- Nadal BPDF model, whose interest is to simulate
C	  the polarisation of different types of surfaces, in a
C	  simplified but representative way.
C        ==> PROBLEM OF APPLICABILITY TO THE SOS FORMALISM
C            MODEL NOT USABLE.
C	- Maignan BPDF model, representative of desert surfaces.
C
C       The Roujean BRDF is added to the last 4 BPDF models
C       (Rondeaux, Breon, Nadal, Maignan).
C
C Input data of the program
C ------------------------------
C    -->Parameters/arguments of the program
C
C            ISURF  (I4)    : (I) Surface type ISURF
C			       (decoding format: integer I4)
C      		0 :Lambertian surface of albedo RHO.
C			1 :Lambertian surface + glitter.
C			2 :Lambertian surface + Fresnel reflection.
C		       3 :Lambertian surface + Roujean BRDF.
C			4 :Lambertian surface + Roujean BRDF + Rondeaux BPDF.
C			5 :Lambertian surface + Roujean BRDF + Breon BPDF.
C			6 :Lambertian surface + Roujean BRDF + Nadal BPDF.
C			7 :Lambertian surface + Roujean BRDF + Maignan BPDF.
C
C            FICANGLES (CHARACTER) : (I) complete name of the file of the angles used
C                                         for the radiance and BRDF computations
C				  (directory + file name + extension)
C
C             IND (I4) : (I) Refractive index surface / air.
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
C             FICLOG (CHARACTER): (I) complete name of the trace file (directory + file name + extension)
C                                    or =0 if no trace
C
C             DIRS (CHARACTER) : (I) Storage directory of the reflection matrix files
C                   		            files (Glitter file or other BRDF/BPDF).
C
C             DIRROUJEAN (CHARACTER) : (I) ==> Associated with: ISURF = 3 4 5 6 or 7
C					storage directory of the Roujean BRDF files,
C                                  applied to the BPDF models.
C
C             DIR_TMP  (CHARACTER)  : (I) temporary directory where the files are generated
C
C     IER (I4) : (O) error code =0 if no error, =-1 otherwise
C
C                                 
C
C    -->Files
C       The Trace file.
C       The result file containing the expression, for the OS code, of the
C       Fourier series decomposition of the surface reflection matrix.
C       Intermediate files (not kept):
C           comments and description in the called routines.
C       The file of angles and internal parameters for the expansion limits.
C         
C       
C
C Results provided by the program
C ----------------------------------
C BRDF / BPDF file of the surface reflection matrix
C for a Fourier series decomposition
C
C Environment variable
C ------------------------
C   None
C
C
C Description of the files used
C ---------------------------------
C
C   -->Content of the result BRDF / BPDF file
C       Unformatted binary file, created by SOS_MISE_FORMAT
C       It contains the surface reflection matrix
C       (GLITTER file or other BRDF / BPDF):
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
C
C    -->Content of the angle file:
C       Formatted ascii file, created by SOS_ANGLES.
C       It contains the information on the angles to be used
C       and on the adapted expansion orders:
C		  NB_TOTAL_ANGLES :  Total nb of angles to be used
C               NB_GAUSS_ANGLES :  Nb of Gauss angles
C               ANGLES_USERFILE :  Name of the user angles file
C	         SOLAR ZENITH ANGLE : Solar zenith angle (degrees)
C		  INTERNAL_IMUS :      Corresponding array index number
C               INTERNAL_OS_NB :   Order of the expansions of the phase functions
C		           in Legendre polynomials to be used.
C             INTERNAL_OS_NS :   Maximum order to be used for the expansion in
C                        Legendre polynomials of the Fresnel matrix elements
C                        and in Fourier series of the radiance computations.
C             INTERNAL_OS_NM :   Maximum order of the expansion in Fourier series
C                                of the G function of the reflection matrix computations.
C               INDEX COS_ANGLE WEIGHT OUTPUT: lines giving the angle number, its cosine
C                                            and its weight, and whether it is a "user" output
C                                              with format I4,1X,2D21.14,1X,I4
C         
C
C Common used:
C --------------
C None
C
C
C Error cases:
C ------------        
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
C     - SOS_NOM_FICGLITTER   |
C     - SOS_NOM_FICROUJ      |
C     - SOS_NOM_FICRH        | -> File SOS_NOM_FIC_SURFACE.F
C     - SOS_NOM_FICBREON     |
C     - SOS_NOM_FICNADAL     |
C     - SOS_GLITTER
C     - SOS_ROUJEAN
C     - SOS_RONDEAUX_BREON
C     - SOS_NADAL
C     - SOS_BPDF_AJOUT_BRDF
C
C============================================================================
      SUBROUTINE SOS_SURFACE(ISURF, FICANGLES, IND, WIND,  
     &                       K0_ROUJEAN, K1_ROUJEAN, K2_ROUJEAN,   
     &                       ALPHA_NADAL, BETA_NADAL,   
     &                       COEF_C_MAIGNAN, FICLOG,                
     &                       DIRS, DIRROUJEAN, DIR_TMP, IER)
      
      IMPLICIT NONE
      
C* Definition of variables
C*-------------------------------------------------------------------------- 
         
      CHARACTER*CTE_LENDIR DIRS	     ! Access path to the result SURFACE file.
      CHARACTER*CTE_LENDIR DIRROUJEAN   ! Storage directory of the Roujean BRDF
                                        ! files, applied to the BPDF models.
      CHARACTER*CTE_LENDIR DIR_TMP      ! Access path to the SURFACE files being generated
      CHARACTER*CTE_LENFIC1 FICS        ! Name of the SURFACE file (without directory tree).
      CHARACTER*CTE_LENFIC2 FICLOG      ! Name of the SURFACE trace file.
      CHARACTER*CTE_LENFIC2 FICANGLES   ! Path and name of the file of the actual angles
      CHARACTER*CTE_LENFIC2 FICANGLES_USER  ! Path and name of the "user" angles file
      CHARACTER*CTE_LENFIC2 FICSURF_TMP     ! Path and name of the SURFACE file being generated
      CHARACTER*CTE_LENFIC2 FIC_RES_MAT_REFLEX ! Name of the file RES_MAT_REFLEX
                                               ! with complete directory tree.
      CHARACTER*CTE_LENFIC2 FIC_RES_GSF ! Name of the file RES_GSF with complete directory tree.
      CHARACTER*CTE_LENFIC2 FIC_RES_FRESNEL ! Name of the file RES_FRESNEL with complete directory tree.
      CHARACTER*CTE_LENFIC2 FICGLITTER ! Complete name of the Glitter file.
      CHARACTER*CTE_LENFIC2 FICROUJ    ! Complete name of the ROUJEAN model file.
      CHARACTER*CTE_LENFIC2 FICRH      ! Complete name of the RH model file.
      CHARACTER*CTE_LENFIC2 FICBREON   ! Complete name of the BREON model file.
      CHARACTER*CTE_LENFIC2 FICNADAL   ! Complete name of the NADAL model file.
      CHARACTER*CTE_LENFIC2 FICMAIGNAN ! Complete name of the MAIGNAN model file.
      CHARACTER*CTE_LENFIC2 FICTMP_BPDF ! Complete name of the temporary file for the BPDF without the BRDF complement.
      CHARACTER*CTE_LENCOM  COMMAND     ! Command string for system call.
  
      CHARACTER*LENLINE LIGNE_TEXTE    ! Character line for decoding parameters
      CHARACTER*1000 CVAL	       ! Character string to read unused values
      
      LOGICAL EX		! Checks the existence of a file.
      LOGICAL TRACE     	! = true if writing in the trace file

      INTEGER*4 NBMU_GAUSS  ! Number of Gauss angles used
      INTEGER*4 LUM_NBMU    ! Number of actual angles of the simulations
      INTEGER*4 LDIRS	! Length of the character string DIRS.
      INTEGER*4 LDIRTMP     ! Length of the character string DIRTMP
      INTEGER*4 LDIRROUJEAN	! Length of the character string DIRROUJEAN
      INTEGER*4 LFICS	! Length of the character string FICS.
      INTEGER*4 LFICTMP	! Length of the character string FICSURF_TMP
      INTEGER*4 LCOM        ! Length of the character string COMMAND
      INTEGER*4 IER     	! Error code =0 if no error, =-1 otherwise
      INTEGER*4 J	       ! Array index
      INTEGER*4 OS_NB 	! Order of the expansions of the phase functions
				! in Legendre polynomials.
      INTEGER*4 OS_NS    	! Maximum order of the expansion in Legendre polynomials
				! of the Fresnel matrix elements.
      INTEGER*4 OS_NM    	! Maximum order of the expansion in Fourier series
				! of the G function of the reflection matrix computations.
  
      INTEGER*2 ISURF	! Surface type

      INTEGER*4 IFIN        ! Nb of significant characters of a string.
 

      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
                                             !Cosines of the Gauss (and user) angles.
      DOUBLE PRECISION CHR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
                                             !Weights of the Gauss angles.
          
      DOUBLE PRECISION WIND	! Wind speed (m/s).
      
      DOUBLE PRECISION IND	! Refractive index of the interface.
      
      DOUBLE PRECISION K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN	!Parameters of the Roujean model
      
      DOUBLE PRECISION ALPHA_NADAL	!Parameter alpha of the Nadal model.
      DOUBLE PRECISION BETA_NADAL	!Parameter beta of the Nadal model.
      DOUBLE PRECISION COEF_C_MAIGNAN !Parameter C*exp(-NDVI) of the Maignan model.
   
      
C* Initialisation of the error code
C------------------------------------------------------
      IER=0
     
C --- Naming of the temporary BRDF/BPDF file
C     for its creation in the work space before moving it
C     to the archive space.
      LDIRTMP = INDEX(DIR_TMP,'/TMP') +3
      FICSURF_TMP=DIR_TMP(1:LDIRTMP)//'/'//"SURFACE_FILE_BUILDING_TMP"
      LFICTMP = INDEX(FICSURF_TMP,' ') - 1

C      Temporary file RES_GSF
      FIC_RES_GSF = DIR_TMP(1:LDIRTMP)//'/'//"RES_GSF"
      FIC_RES_FRESNEL = DIR_TMP(1:LDIRTMP)//'/'//"RES_FRESNEL"
      FIC_RES_MAT_REFLEX = DIR_TMP(1:LDIRTMP)//'/'//"RES_MAT_REFLEX"
      FICTMP_BPDF = DIR_TMP(1:LDIRTMP)//'/'//"FICTMP_BPDF"

C     Removal of the possibly existing temporary files
      INQUIRE(FILE=FICSURF_TMP,EXIST=EX)
      IF (EX) THEN 
         OPEN(UNIT=3,FILE=FICSURF_TMP,FORM='UNFORMATTED')
         CLOSE(3,STATUS='DELETE') 
      ENDIF 

      INQUIRE(FILE=FIC_RES_GSF,EXIST=EX)
      IF (EX) THEN 
         OPEN(UNIT=3,FILE=FIC_RES_GSF,FORM='UNFORMATTED')
         CLOSE(3,STATUS='DELETE') 
      ENDIF

      INQUIRE(FILE=FIC_RES_FRESNEL,EXIST=EX)
      IF (EX) THEN 
         OPEN(UNIT=3,FILE=FIC_RES_FRESNEL,FORM='UNFORMATTED')
         CLOSE(3,STATUS='DELETE') 
      ENDIF

      INQUIRE(FILE=FIC_RES_MAT_REFLEX,EXIST=EX)
      IF (EX) THEN 
         OPEN(UNIT=3,FILE=FIC_RES_MAT_REFLEX,FORM='UNFORMATTED')
         CLOSE(3,STATUS='DELETE') 
      ENDIF

      INQUIRE(FILE=FICTMP_BPDF,EXIST=EX)
      IF (EX) THEN 
         OPEN(UNIT=3,FILE=FICTMP_BPDF,FORM='UNFORMATTED')
         CLOSE(3,STATUS='DELETE') 
      ENDIF


C* Opening of the trace file
C------------------------------------------------------
      IF (FICLOG.EQ.'NO_LOG_FILE') THEN
	  TRACE=.FALSE.
      ELSE
         TRACE=.TRUE.
	  OPEN (99,FILE=FICLOG,ERR=911)
      ENDIF
            	 

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
      
      READ(10,*,ERR=951) CVAL 
      READ(10,*,ERR=951) CVAL
      
      READ(10,'(a)',ERR=951) LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) OS_NB
      
      READ(10,'(a)',ERR=951) LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) OS_NS
      
      READ(10,'(a)',ERR=951) LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) OS_NM
      
      READ(10,*,ERR=951) CVAL
      
      DO J=1,LUM_NBMU
         READ(10,*,ERR=951) CVAL, RMU(J), CHR(J), CVAL
      ENDDO 

      DO J=1,LUM_NBMU
         RMU(-J) = -RMU(J)
	 CHR(-J) =  CHR(J)
      ENDDO 

	
      
      IF (TRACE) THEN
        WRITE (99,*,ERR=960) "USED ANGLES"
	 WRITE (99,*,ERR=960) "***********"
	 WRITE (99,*,ERR=960) " "
	 WRITE(99,*,ERR=960) "Total number of angles : ", LUM_NBMU
	 WRITE(99,*,ERR=960) "Number of Gauss angles : ", NBMU_GAUSS
	 WRITE (99,*,ERR=960) " "
        IF (FICANGLES_USER.NE.'NO_USER_ANGLES') THEN
	     WRITE(99,*,ERR=960)  
     &	       "A file of user's angles is used :",
     &         FICANGLES_USER
     	     WRITE(99,*,ERR=960) "Number of user's angles :",
     &         (LUM_NBMU-NBMU_GAUSS)
	 ELSE
	     WRITE(99,*,ERR=960)  
     &	       "No file of user's angles"
	 ENDIF
        WRITE (99,*,ERR=960) " " 
	 WRITE (99,*,ERR=960) "Max order NB ", OS_NB
	 WRITE (99,*,ERR=960) "Max order NS ", OS_NS
	 WRITE (99,*,ERR=960) "Max order NM ", OS_NM
	 WRITE (99,*,ERR=960) " " 
	 WRITE (99,*,err=960) "List of used angles :"
	 WRITE (99,*,err=960) "   Index, mu, angle (degrees), weigth"
	 DO 1000 J=-LUM_NBMU,LUM_NBMU
	    IF (J.EQ.0) GOTO 1000
	    IF (J.LT.0) THEN 
	        WRITE(99,222,ERR=960) J, RMU(J),
     &	                  -DACOS(-RMU(J))*180./INCTE_PI, CHR(J)
            ELSE
	        WRITE(99,222,ERR=960) J, RMU(J),
     &	                   DACOS(RMU(J))*180./INCTE_PI, CHR(J)
	    ENDIF
1000    CONTINUE	
	 WRITE (99,*,ERR=960) " " 
	 WRITE (99,*,ERR=960) " "  	
      ENDIF !Trace file


c* Size of the complete path to the output directory
C* -----------------------------------------------------------------     

      LDIRS = INDEX(DIRS,' ')
      LDIRS = LDIRS - 1
      IF(LDIRS.LE.0) LDIRS = CTE_LENDIR

      LDIRROUJEAN = INDEX(DIRROUJEAN,' ')
      LDIRROUJEAN = LDIRROUJEAN - 1
      IF(LDIRROUJEAN.LE.0) LDIRROUJEAN = CTE_LENDIR


C =====================================================================
C* GLITTER case
C =====================================================================

      IF (ISURF.EQ.1) THEN

	
C* Building of the name of the GLITTER file
C  --------------------------------------
	 CALL SOS_NOM_FICGLITTER(NBMU_GAUSS,FICANGLES_USER,IND,WIND,
     &	                         FICS)      
	 FICGLITTER = DIRS(1:LDIRS)//'/'//FICS	
	 LFICS=INDEX(FICGLITTER,' ')-1	
	 WRITE(6,*) "File : ",FICGLITTER(1:LFICS)
	 
C* If the GLITTER file already exists, the program is not run
C  ---------------------------------------------------------------     
         INQUIRE(FILE=FICGLITTER(1:LFICS),EXIST=EX)
	  
         IF (EX) THEN  
            WRITE(6,*)'==> The file has already been calculated'
	  
	 ELSE   
      
C* Computation of the Glitter file
C  ------------------------- 
           WRITE(6,*)'==> Starting the file computation'
	    CALL SOS_GLITTER(LUM_NBMU,RMU,CHR,WIND,IND,
     &	                     OS_NB,OS_NS,OS_NM,
     &                      FIC_RES_GSF,FIC_RES_FRESNEL,
     &                      FIC_RES_MAT_REFLEX,
     &	                     FICSURF_TMP,TRACE,IER)
            IF (IER.NE.0) GOTO 993   ! error case

C* Move of FICSURF_TMP to FICGLITTER
C---------------------------------------------------------------------
            COMMAND=
     &      'mv '//FICSURF_TMP(1:LFICTMP)//' '//FICGLITTER(1:LFICS)
            LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   !The search starts after the last blank at the beginning of the string
            LCOM = LCOM + (4+LFICTMP)
            IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020

            INQUIRE(FILE=FICGLITTER(1:LFICS),EXIST=EX)
            IF (EX) THEN !If the BRDF/BPDF file exists, the temporary file (generated by another process) is deleted
               OPEN(UNIT=3,FILE=FICSURF_TMP,FORM='UNFORMATTED')
               CLOSE(3,STATUS='DELETE') 
            ELSE !Otherwise it is moved to the storage space of the BRDF/BPDF files
               CALL SYSTEM(COMMAND)  
            ENDIF    
         ENDIF

	  
      ENDIF ! Glitter case
	 



C =====================================================================
C* ROUJEAN BRDF case
C*   NB: Also computed for the BPDF models.
C =====================================================================

      IF (ISURF.GE.3) THEN


C* Consistency check with respect to the format F7.3
C  Values of K0, K1, K2 forced to the format F7.3
C  -----------------------------------------	
        IF((K0_ROUJEAN.LT.0.).OR.(K0_ROUJEAN.GE.10.)) GOTO 9910
	 IF((K1_ROUJEAN.LT.0.).OR.(K1_ROUJEAN.GE.10.)) GOTO 9910
	 IF((K2_ROUJEAN.LT.0.).OR.(K2_ROUJEAN.GE.10.)) GOTO 9910
	
	 K0_ROUJEAN = DNINT(K0_ROUJEAN*1000.D+00) / 1000.D+00
	 K1_ROUJEAN = DNINT(K1_ROUJEAN*1000.D+00) / 1000.D+00
	 K2_ROUJEAN = DNINT(K2_ROUJEAN*1000.D+00) / 1000.D+00	 	 

	
C* Building of the name of the ROUJEAN file
C  -----------------------------------------
	 CALL SOS_NOM_FICROUJ(NBMU_GAUSS,FICANGLES_USER,
     &                       K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,FICS)     
	 FICROUJ = DIRROUJEAN(1:LDIRROUJEAN)//'/'//FICS 
	 LFICS=INDEX(FICROUJ,' ')-1      
	 WRITE(6,*) "File : ",FICROUJ(1:LFICS)
	 
C* If the file already exists, the program is not run
C  --------------------------------------------------------	 
         INQUIRE(FILE=FICROUJ(1:LFICS),EXIST=EX)	 
	 
         IF (EX) THEN  
            WRITE(6,*)'==> The file has already been calculated'
	    
	 ELSE   
      
C* Computation of the ROUJEAN file
C  -----------------------------
            WRITE(6,*)'==> Starting the file computation'
	     CALL SOS_ROUJEAN(LUM_NBMU,RMU,CHR,OS_NB,
     &	                      K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,
     &	                      FIC_RES_MAT_REFLEX,FICSURF_TMP,TRACE,IER)

            IF (IER.NE.0) GOTO 994   ! error case

C* Move of FICSURF_TMP to FICROUJ
C---------------------------------------------------------------------
            COMMAND=
     &      'mv '//FICSURF_TMP(1:LFICTMP)//' '//FICROUJ(1:LFICS)
            LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   !The search starts after the last blank at the beginning of the string
            LCOM = LCOM + (4+LFICTMP)
            IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020

            INQUIRE(FILE=FICROUJ(1:LFICS),EXIST=EX)
            IF (EX) THEN !If the BRDF/BPDF file exists, the temporary file (generated by another process) is deleted
               OPEN(UNIT=3,FILE=FICSURF_TMP,FORM='UNFORMATTED')
               CLOSE(3,STATUS='DELETE') 
            ELSE !Otherwise it is moved to the storage space of the BRDF/BPDF files
               CALL SYSTEM(COMMAND)  
            ENDIF
	       
         ENDIF
	 
      ENDIF ! Case of the Roujean BRDF



	 
C =====================================================================
C* Case of the Rondeaux - Herman BPDF
C =====================================================================

      IF (ISURF.EQ.4) THEN

C* Building of the name of the Rondeaux - Herman file
C  ---------------------------------------------------
	 CALL SOS_NOM_FICRH(NBMU_GAUSS,FICANGLES_USER,IND,
     &	                    K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,FICS)
	 FICRH = DIRS(1:LDIRS)//'/'//FICS	 
	 LFICS=INDEX(FICRH,' ')-1
	 WRITE(6,*) "File : ",FICRH(1:LFICS)
	 
C* If the file already exists, the program is not run
C  -------------------------------------------------------	 
        INQUIRE(FILE=FICRH(1:LFICS),EXIST=EX)
	 	 
        IF (EX) THEN  
            WRITE(6,*)'==> The file has already been calculated'
	  
	 ELSE   
  
            
C* Computation of the Rondeaux - Herman file
C  --------------------------------------
           WRITE(6,*)'==> Starting the file computation'
	    CALL SOS_SURFACE_BPDF(LUM_NBMU,RMU,CHR,IND,
     &                           ISURF,ALPHA_NADAL,BETA_NADAL,  
     &                           COEF_C_MAIGNAN,
     &	                          OS_NB,OS_NS,OS_NM,
     &                           FIC_RES_GSF,FIC_RES_FRESNEL,
     &                           FIC_RES_MAT_REFLEX,
     &	                          FICTMP_BPDF,TRACE,IER)     
            IF (IER.NE.0) GOTO 995   ! error case


		
C* Addition of the Roujean BRDF
C---------------------------------------------------------------------
            CALL SOS_BPDF_AJOUT_BRDF(FICTMP_BPDF,FICROUJ,LUM_NBMU,OS_NB,
     &	                              FICSURF_TMP,IER)
	     IF (IER.NE.0) GOTO 998   ! error case

C* Move of FICSURF_TMP to FICRH
C---------------------------------------------------------------------
            COMMAND=
     &      'mv '//FICSURF_TMP(1:LFICTMP)//' '//FICRH(1:LFICS)
            LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   !The search starts after the last blank at the beginning of the string
            LCOM = LCOM + (4+LFICTMP)
            IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020

            INQUIRE(FILE=FICRH(1:LFICS),EXIST=EX)
            IF (EX) THEN !If the BRDF/BPDF file exists, the temporary file (generated by another process) is deleted
               OPEN(UNIT=3,FILE=FICSURF_TMP,FORM='UNFORMATTED')
               CLOSE(3,STATUS='DELETE') 
            ELSE !Otherwise it is moved to the storage space of the BRDF/BPDF files
               CALL SYSTEM(COMMAND)  
            ENDIF


C* Deletion of the temporary file
C---------------------------------------   
	     OPEN(UNIT=3,FILE=FICTMP_BPDF,FORM='UNFORMATTED')
            CLOSE(3,STATUS='DELETE') 
	    	                   
	        	       
        ENDIF ! If the file is to be created.

      ENDIF ! Case of the Rondeaux - Herman BPDF



C =====================================================================
C* Case of the BREON BPDF
C =====================================================================

      IF (ISURF.EQ.5) THEN

C* Building of the name of the BREON file
C  ---------------------------------------
	 CALL SOS_NOM_FICBREON(NBMU_GAUSS,FICANGLES_USER,IND,
     &	                       K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,FICS)   
	 FICBREON = DIRS(1:LDIRS)//'/'//FICS
	 LFICS=INDEX(FICBREON,' ')-1
	 WRITE(6,*) FICBREON(1:LFICS)
	 
C* If the file already exists, the program is not run
C  --------------------------------------------------------	 
        INQUIRE(FILE=FICBREON(1:LFICS),EXIST=EX)	 
	 
        IF (EX) THEN  
            WRITE(6,*)'==> The file has already been calculated'
	  
	 ELSE        
            
C* Computation of the Breon file
C  --------------------------------------
            WRITE(6,*)'==> Starting the file computation'
	     CALL SOS_SURFACE_BPDF(LUM_NBMU,RMU,CHR,IND,
     &                            ISURF,ALPHA_NADAL,BETA_NADAL,  
     &                            COEF_C_MAIGNAN,
     &	                           OS_NB,OS_NS,OS_NM,
     &                            FIC_RES_GSF,FIC_RES_FRESNEL,
     &                            FIC_RES_MAT_REFLEX,
     &	                           FICTMP_BPDF,TRACE,IER)  
            IF (IER.NE.0) GOTO 995   ! error case

C* Addition of the Roujean BRDF
C---------------------------------------------------------------------
            CALL SOS_BPDF_AJOUT_BRDF(FICTMP_BPDF,FICROUJ,LUM_NBMU,OS_NB,
     &                               FICSURF_TMP,IER)
	     IF (IER.NE.0) GOTO 998   ! error case

C* Move of FICSURF_TMP to FICBREON
C---------------------------------------------------------------------
            COMMAND=
     &      'mv '//FICSURF_TMP(1:LFICTMP)//' '//FICBREON(1:LFICS)
            LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   !The search starts after the last blank at the beginning of the string
            LCOM = LCOM + (4+LFICTMP)
            IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020

            INQUIRE(FILE=FICBREON(1:LFICS),EXIST=EX)
            IF (EX) THEN !If the BRDF/BPDF file exists, the temporary file (generated by another process) is deleted
               OPEN(UNIT=3,FILE=FICSURF_TMP,FORM='UNFORMATTED')
               CLOSE(3,STATUS='DELETE') 
            ELSE !Otherwise it is moved to the storage space of the BRDF/BPDF files
               CALL SYSTEM(COMMAND)  
            ENDIF
	    
C* Deletion of the temporary file
C---------------------------------------   
	     OPEN(UNIT=3,FILE=FICTMP_BPDF,FORM='UNFORMATTED')
            CLOSE(3,STATUS='DELETE') 
	       
        ENDIF ! If the file is to be created.

      ENDIF ! Case of the Breon BPDF



C =====================================================================
C* Case of the NADAL BPDF
C =====================================================================

      IF (ISURF.EQ.6) THEN	

		
C* Consistency check with respect to the formats F6.4 and F4.1
C  Values of alpha and beta forced to these formats
C  -----------------------------------------	
        IF((ALPHA_NADAL.LT.0.).OR.(ALPHA_NADAL.GE.10.)) GOTO 9915
	 IF((BETA_NADAL.LT.0.).OR.(BETA_NADAL.GE.100.)) GOTO 9915

	 ALPHA_NADAL = DNINT(ALPHA_NADAL*10000.D+00) / 10000.D+00
	 BETA_NADAL = DNINT(BETA_NADAL*10.D+00) / 10.D+00	
		
C* Building of the name of the NADAL file
C  ---------------------------------------
	 CALL SOS_NOM_FICNADAL(NBMU_GAUSS,FICANGLES_USER,IND,
     &	                       ALPHA_NADAL,BETA_NADAL,
     &	                       K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,FICS) 	        
	 FICNADAL = DIRS(1:LDIRS)//'/'//FICS	 
	 LFICS=INDEX(FICNADAL,' ')-1
	 WRITE(6,*) "File : ",FICNADAL(1:LFICS)
	 
C* If the file already exists, the program is not run
C  --------------------------------------------------------	 
        INQUIRE(FILE=FICNADAL(1:LFICS),EXIST=EX)
	 	 
        IF (EX) THEN  
            WRITE(6,*)'==> The file has already been calculated'
	  
	 ELSE   
            
C* Computation of the NADAL file
C  ---------------------------
            WRITE(6,*)'==> Starting the file computation'
	     CALL SOS_SURFACE_BPDF(LUM_NBMU,RMU,CHR,IND,
     &                            ISURF,ALPHA_NADAL,BETA_NADAL,  
     &                            COEF_C_MAIGNAN,
     &	                           OS_NB,OS_NS,OS_NM,
     &                            FIC_RES_GSF,FIC_RES_FRESNEL,
     &                            FIC_RES_MAT_REFLEX,
     &	                           FICTMP_BPDF,TRACE,IER)  
            IF (IER.NE.0) GOTO 997   ! error case

C* Addition of the Roujean BRDF
C---------------------------------------------------------------------
            CALL SOS_BPDF_AJOUT_BRDF(FICTMP_BPDF,FICROUJ,LUM_NBMU,OS_NB,
     &                               FICSURF_TMP,IER)
	     IF (IER.NE.0) GOTO 998   ! error case

C* Move of FICSURF_TMP to FICNADAL
C---------------------------------------------------------------------
            COMMAND=
     &      'mv '//FICSURF_TMP(1:LFICTMP)//' '//FICNADAL(1:LFICS)
            LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   !The search starts after the last blank at the beginning of the string
            LCOM = LCOM + (4+LFICTMP)
            IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020

            INQUIRE(FILE=FICNADAL(1:LFICS),EXIST=EX)
            IF (EX) THEN !If the BRDF/BPDF file exists, the temporary file (generated by another process) is deleted
               OPEN(UNIT=3,FILE=FICSURF_TMP,FORM='UNFORMATTED')
               CLOSE(3,STATUS='DELETE') 
            ELSE !Otherwise it is moved to the storage space of the BRDF/BPDF files
               CALL SYSTEM(COMMAND)  
            ENDIF  

C* Deletion of the temporary file
C---------------------------------------   
	     OPEN(UNIT=3,FILE=FICTMP_BPDF,FORM='UNFORMATTED')
            CLOSE(3,STATUS='DELETE') 
	       
        ENDIF ! If the file is to be created.
	 
      ENDIF ! Case of the Nadal BPDF




C =====================================================================
C* Case of the MAIGNAN BPDF
C =====================================================================

      IF (ISURF.EQ.7) THEN

C* Building of the name of the MAIGNAN file
C  ---------------------------------------
	 CALL SOS_NOM_FICMAIGNAN(NBMU_GAUSS,FICANGLES_USER,IND,
     &	                         COEF_C_MAIGNAN,
     &	                         K0_ROUJEAN,K1_ROUJEAN,K2_ROUJEAN,FICS)	        
	 FICMAIGNAN = DIRS(1:LDIRS)//'/'//FICS
	 LFICS=INDEX(FICMAIGNAN,' ')-1  
	 WRITE(6,*) "File : ",FICMAIGNAN(1:LFICS)
	 
C* If the file already exists, the program is not run
C  --------------------------------------------------------	 
        INQUIRE(FILE=FICMAIGNAN(1:LFICS),EXIST=EX)	 
	 
        IF (EX) THEN  
            WRITE(6,*)'==> The file has already been calculated'
	  
	 ELSE        
            
C* Computation of the Maignan file
C  --------------------------------------
            WRITE(6,*)'==> Starting the file computation'

	     CALL SOS_SURFACE_BPDF(LUM_NBMU,RMU,CHR,IND,
     &                            ISURF,ALPHA_NADAL,BETA_NADAL,  
     &                            COEF_C_MAIGNAN,
     &	                           OS_NB,OS_NS,OS_NM,
     &                            FIC_RES_GSF,FIC_RES_FRESNEL,
     &                            FIC_RES_MAT_REFLEX,
     &	                           FICTMP_BPDF,TRACE,IER)    
            IF (IER.NE.0) GOTO 996   ! error case


C* Addition of the Roujean BRDF
C---------------------------------------------------------------------
            CALL SOS_BPDF_AJOUT_BRDF(FICTMP_BPDF,FICROUJ,LUM_NBMU,OS_NB,
     &                               FICSURF_TMP,IER)
	     IF (IER.NE.0) GOTO 998   ! error case

C* Move of FICSURF_TMP to FICMAIGNAN
C---------------------------------------------------------------------
            COMMAND=
     &      'mv '//FICSURF_TMP(1:LFICTMP)//' '//FICMAIGNAN(1:LFICS)
            LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   !The search starts after the last blank at the beginning of the string
            LCOM = LCOM + (4+LFICTMP)
            IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020

            INQUIRE(FILE=FICMAIGNAN(1:LFICS),EXIST=EX)
            IF (EX) THEN !If the BRDF/BPDF file exists, the temporary file (generated by another process) is deleted
               OPEN(UNIT=3,FILE=FICSURF_TMP,FORM='UNFORMATTED')
               CLOSE(3,STATUS='DELETE') 
            ELSE !Otherwise it is moved to the storage space of the BRDF/BPDF files
               CALL SYSTEM(COMMAND)  
            ENDIF
	    
C* Deletion of the temporary file
C---------------------------------------   
	     OPEN(UNIT=3,FILE=FICTMP_BPDF,FORM='UNFORMATTED')
            CLOSE(3,STATUS='DELETE') 
	       
        ENDIF ! If the file is to be created.

      ENDIF ! Case of the Maignan BPDF




C* File closing
C-------------------

      IF (TRACE) CLOSE(99)      
      CLOSE(10) !Angle file
       
C* Nominal end
C-------------------
      GOTO 9999
            
C* Error case and return of status 1 to the shell
C----------------------------------------------
  800 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.AngFile Value'
      GOTO 9998 
      
  801 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Type Value (format : I2)'
      GOTO 9998  
      
  802 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Log Value'
      GOTO 9998       
            
  803 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Dir Value'
      GOTO 9998        
                   
  804 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.DirBRDFtoBPDFcalculation Value'
      GOTO 9998        
                   
  805 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.ind Value'
      GOTO 9998        
             
  810 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Glitter.Wind Value (format : F4)'
      GOTO 9998        
           
  811 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Roujean.K0 Value (format : F7.3)',
     &           '     -SURF.Roujean.K1 Value (format : F7.3)',
     &           '     -SURF.Roujean.K2 Value (format : F7.3)'
      GOTO 9998        
           
  812 WRITE(6,*) '  ERROR on parameters for SOS_SURFACE'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SURF.Nadal.Alpha Value (format : F6.4)',
     &           '     -SURF.Nadal.Beta Value (format : F4.1)'
      GOTO 9998        
                    
  990 WRITE(6,*) '  ERROR on argument loading for SOS_SURFACE'
      WRITE(6,*) '  --> Check required formats.'
      GOTO 9998  
      
 9901 WRITE(6,*) '  ERROR on -SURF.Type value : out of available range.'
      CALL EXIT(1)
      GOTO 9998  
      
 9910 WRITE(6,*) '  ERROR on -SURF.Roujean.K* values : '
      WRITE(6,*) '  out of available range : ]0, 1000.]'
      GOTO 9998 
      
 9915 WRITE(6,*) '  ERROR on -SURF.Nadal.Alpha or .Beta values : '
      WRITE(6,*) '  out of available range : 0 < Alpha <= 10',
     &           ' & 0 < Beta <= 100. '
      GOTO 9998      
              
  911 WRITE(6,*) '  ERROR on logfile opening for SOS_SURFACE'
      GOTO 9998   
  
  950 WRITE(6,*) '  ERROR on angles file opening for SOS_SURFACE'
      GOTO 9998

  951 WRITE(6,*) '  ERROR on angles file reading for SOS_SURFACE'
      GOTO 9998
      
  960 WRITE(6,*) '  ERROR on logfile writing for SOS_SURFACE'
      GOTO 9998
                     
  993 WRITE(6,*) '  ERROR on subroutine SOS_GLITTER'
      GOTO 9998
      
  994 WRITE(6,*) '  ERROR on subroutine SOS_ROUJEAN'
      GOTO 9998
      
 9951 WRITE(6,*) '  ERROR while deleting temp. file MAT_BPDF_FORMATOS'
      GOTO 9998      
      
  995 WRITE(6,*) '  ERROR on subroutine SOS_RONDEAUX_BREON'
      GOTO 9998
      
  996 WRITE(6,*) '  ERROR on subroutine SOS_MAIGNAN'
      GOTO 9998
      
  997 WRITE(6,*) '  ERROR on subroutine SOS_NADAL'
      GOTO 9998
      
  998 WRITE(6,*) '  ERROR on subroutine SOS_BPDF_AJOUT_BRDF'
      GOTO 9998


 1020 WRITE(6,*) '  ERROR : mv command'
      WRITE(6,*) '  => The command to move the temporary SURFACE file',
     &           ' to the final depository is maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998

9998  IER=-1
      GOTO 9999

C* Format
C---------   
  222 FORMAT(I4,1X,D21.14,2X,F7.3,2X,D21.14,2X)

    
 9999 RETURN
 	 

      END 	!END OF PROGRAM SOS_SURFACE
      
      
      


C==============================================================================
C PROCEDURE: SOS_MAT_FRESNEL
C ==========
C      This procedure computes the expansion of the Fresnel matrix
C      in the basis of the Legendre polynomials.
C			
C      The Fresnel matrix, expressed in a frame tied to the reflection plane,
C      is:
C                  | R11(teta)  R12(teta)    0      |
C         F(teta)= | R12(teta)  R11(teta)    0      |
C                  |    0          0      R33(teta) |
C
C         with:
C           _ teta : reflection angle with respect to the normal of the interface,
C           _ R11(teta) = 0.5 * (RL*RL + RR*RR),
C           _ R12(teta) = 0.5 * (RL*RL - RR*RR),
C           _ R33(teta) = RL*RR,
C
C         The relative refractive index of the interface, between the reflection
C         medium and the transmission medium, is IND.
C
C         RL and RR are computed from IND and teta by:
C           _ RL : the parallel reflection coefficient
C
C               RL(teta) = [ IND**2 * cos(teta) - SQRT(IND**2 - sin(teta)**2) ]
C                        / [ IND**2 * cos(teta) + SQRT(IND**2 - sin(teta)**2) ]
C
C           _ RR : the perpendicular reflection coefficient
C
C               RR(teta) = [ cos(teta) - SQRT(IND**2 - sin(teta)**2) ]
C                        / [ cos(teta) + SQRT(IND**2 - sin(teta)**2) ]
C
C
C
C      The expansion gives the parameters ALPHA(k), BETA(k), GAMMA(k) and
C      ZETA(k) associated with the Fresnel matrix.
C      They are stored from k=0 to OS_NS in the file "RES_FRESNEL".
C
C
C      NB: IMPLEMENTATION APPROACH:
C      *******************************
C          The implementation of the code is adapted to a computation of the
C          reflection matrix (routine SOS_MAT_REFLEXION) using only positive
C           values of cos(Gauss angle). To achieve this:
C          - The complementary angle of the true scattering angle is used
C            (2*reflection angle): DEM 1
C            (routine SOS_MAT_FRESNEL).
C          - The opposite of the true parallel reflection coefficient
C            is used: DEM 2 (routine SOS_MAT_FRESNEL).
C
C          Downstream, the routine computing the reflection matrix introduces
C          a coefficient (-1)**N in the summation on the terms of the
C          Fourier series to make the whole consistent: DEM 3.
C
C
C Description of parameters
C --------------------------
C      LUM_NBMU (I4) : (I) Number of actually useful (positive) angles.
C	RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double) : (I) Array of the Gauss mu
C                               (complementary angles of the scattering angles).
C	CHR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double) : (I) Array of the Gauss weights.
C 	IND   (double) :   (I) Refractive index water / air.
C       OS_NS (I4) : (I) Maximum order to be used for the expansion in
C                        Legendre polynomials of the Fresnel matrix elements
C                        and in Fourier series of the radiance computations.
C      FIC_RES_FRESNEL (CHARACTER)  : (O) Name of the temporary file RES_FRESNEL
C                                      with complete directory tree.
C      TRACE (logical)  : (I) = true if writing in the trace file
C              (logical unit number of the trace file =99)
C     IER (I4) : (O) error code =0 if no error, =-1 otherwise
C
C Other result provided by the procedure
C -------------------------------------
C      Temporary file RES_FRESNEL
C      Information in the trace file (optional)
C         (logical unit number of the trace file =99,
C          file created/opened by the calling program)
C
C 
C Description of the files used
C ---------------------------------
C   -->Content of the temporary file RES_FRESNEL
C      Ascii file, record format 4E15.8,
C      created by SOS_MAT_FRESNEL
C      File created in the directory where the program is run
C      If a file with the same name already exists, it is overwritten
C      Number of records: from k=0 to OS_NS
C      Each record contains the parameters ALPHA(k), BETA(k),
C      GAMMA(k) and ZETA(k) associated with the Fresnel matrix
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
C     -Reading/writing in a file
C
C     Message displayed on screen, the routine stops its computations and
C     returns status -1 to the calling program
C
C
C==============================================================================
      SUBROUTINE SOS_MAT_FRESNEL(LUM_NBMU,RMU,CHR,IND,OS_NS,
     &                           FIC_RES_FRESNEL,TRACE,IER)

      IMPLICIT NONE

C* Definition of variables
C*----------------------------------------------------------------------
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Array of the Gauss mu
      				! (complementary angles of the scattering angles).
      DOUBLE PRECISION CHR(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Array of the Gauss weights.

      DOUBLE PRECISION IND	! Refractive index water / air.
	
      DOUBLE PRECISION ALPHA(0:CTE_OS_NS_MAX)    ! Parameters alpha(k).
      DOUBLE PRECISION BETA(0:CTE_OS_NS_MAX)	! Parameters beta(k).
      DOUBLE PRECISION GAMMA(0:CTE_OS_NS_MAX)    ! Parameters gamma(k).
      DOUBLE PRECISION ZETA(0:CTE_OS_NS_MAX)	! Parameters zeta(k).

      DOUBLE PRECISION DELTA(0:CTE_OS_NS_MAX)	! Parameters delta(k).

      DOUBLE PRECISION A	!
      DOUBLE PRECISION B	!
      DOUBLE PRECISION C	!

      DOUBLE PRECISION RL	! Parallel reflection coefficient.
      DOUBLE PRECISION RR	! Perpendicular reflection coefficient.

      DOUBLE PRECISION R11(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
                                ! Element R11 of the Fresnel matrix
      				! versus the reflection angle.
      DOUBLE PRECISION R12(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)     
      				! Element R12 of the Fresnel matrix.
      DOUBLE PRECISION R33(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)    
      				! Element R33 of the Fresnel matrix

      DOUBLE PRECISION IT	! Value of R11 recomputed from its
      				! decomposition in BETA(k).
      DOUBLE PRECISION QT	! Value of R12 recomputed from its
      				! decomposition in GAMMA(k).
      DOUBLE PRECISION UT	! Value of R33 recomputed from its
      				! decomposition in DELTA(k).

      DOUBLE PRECISION PL(-1:CTE_OS_NS_MAX+1)
                      !Legendre polynomials of order L computed for one angle
      DOUBLE PRECISION POL(0:CTE_OS_NS_MAX+1)
                      !Legendre function of order S=2 and variable order L
		        !(computed for one angle)
		       			
      DOUBLE PRECISION D	!Variable for the computation of POL.
      DOUBLE PRECISION E	!Variable for the computation of POL.

      DOUBLE PRECISION X
      DOUBLE PRECISION XRMU	!Value used from the array of Gauss angles

      DOUBLE PRECISION X2				
      DOUBLE PRECISION XX	
      DOUBLE PRECISION XXX

      DOUBLE PRECISION CO1     !Variable for the computation of ALPHA(k) and ZETA(k)
      DOUBLE PRECISION CO2     !Variable for the computation of ALPHA(k) and ZETA(k)
      DOUBLE PRECISION CO3     !Variable for the computation of ALPHA(k) and ZETA(k)

      DOUBLE PRECISION SOM1    !Variable for the computation of ALPHA(k) and ZETA(k)
      DOUBLE PRECISION SOM2    !Variable for the computation of ALPHA(k) and ZETA(k)
      DOUBLE PRECISION SOM3    !Variable for the computation of ALPHA(k) and ZETA(k)
      DOUBLE PRECISION SOM4    !Variable for the computation of ALPHA(k) and ZETA(k)

      LOGICAL TRACE        ! = true if writing in the trace file

      INTEGER*4 LUM_NBMU ! Number of actual angles of the simulations

      INTEGER*4 OS_NS   ! Maximum order of the expansion in Legendre polynomials
			   ! of the Fresnel matrix elements.
      INTEGER*4 IER     ! error code =0 if no error, =-1 otherwise
      INTEGER*4 J	! Array index.
      INTEGER*4 K	! Array index.
      INTEGER*4	NN
      INTEGER*4	MM
      INTEGER*4 I	

      CHARACTER*CTE_LENFIC2 FIC_RES_FRESNEL ! Name of the file RES_FRESNEL
                                            ! with complete directory tree


C* Opening of the result file
C---------------------------------
      OPEN(UNIT=3,FILE=FIC_RES_FRESNEL,ERR=991)

C* Initialisation
C---------------------------------	
      DO 101 K=0,OS_NS
         BETA(K)=0.D+00
         GAMMA(K)=0.D+00
         ALPHA(K)=0.D+00
         ZETA(K)=0.D+00
         DELTA(K)=0.D+00
  101 CONTINUE
  
C* Trace file
C-----------------
      IF (TRACE) THEN
         WRITE(99,*,err=921)
	 WRITE(99,*,err=921) 
     &	    'Calculation of coefficients of the Fresnel matrix ' //
     &	    'expansion into Legendre functions' 
      END IF
			
C* Computation of the elements of the Fresnel reflection matrix
C--------------------------------------------------------------
      DO 104 J=-LUM_NBMU,LUM_NBMU
         IF(J.EQ.0) GOTO 104

C* Cosine of the complementary angle of the scattering angle:
C* scattering angle = PI - 2* reflection angle.
C-------------------------------------------------------------	 	
         C=RMU(J)				! DEM 1
	
C* Cosine of the reflection angle:
C------------------------------------------------------
         C=DSQRT(.5*(1+C))


C* Parallel reflection coefficients (opposite of the true value)
C* and perpendicular
C-----------------------------------------------------------------	
         A=DSQRT(IND*IND-1.D+00+C*C)
	 B=IND*IND*C
         RL=-(B-A)/(B+A)			! DEM 2  	
         RR=(C-A)/(C+A)

C* Elements of the Fresnel matrix
C-----------------------------------	  	   		
         R11(J)=.5*(RL*RL+RR*RR)
         R12(J)=.5*(RL*RL-RR*RR)
         R33(J)=RL*RR
	 
C* Trace file
C-----------------
         IF (TRACE) THEN
	   WRITE(99,200,err=921) 
     &           ' Refracted angle:',DACOS(C)*90./INCTE_PI,
     &           ' R11 ',R11(J),' R12 ',R12(J),' R33 ',R33(J)
         END IF
	 
 104  CONTINUE


C* Computation of the parameters BETA(k) of the decomposition of the function
C* R11(teta) in Legendre polynomials.
C------------------------------------------------------------------
      DO 4 J=-LUM_NBMU,LUM_NBMU
         IF(J.EQ.0) GO TO 4
         X=R11(J)*CHR(J)
         XRMU=RMU(J)

C* Computation of the Legendre polynomials of order K, for the fixed value of mu.
C------------------------------------------------------
         PL(-1)=0.D+00
         PL(0)=1.D+00
         DO 5 K=0,OS_NS
            PL(K+1)=((2*K+1.)*XRMU*PL(K)-K*PL(K-1))/(K+1.)
            BETA(K)=BETA(K)+X*PL(K)
    5    CONTINUE
    4 CONTINUE

      DO 6 K=0,OS_NS
           BETA(K)=(2*K+1)*BETA(K)*.5
    6 CONTINUE
    
C* Trace file
C-----------------
      IF (TRACE) THEN
          WRITE(99,*,err=921) '  '
	   WRITE(99,*,err=921) ' Recomposition of the R11 function'
      
C* Check of the recombination of R11
C--------------------------------------
           DO 8 J=-LUM_NBMU,LUM_NBMU
              IT=0.D+00
              IF(J.EQ.0) GO TO 8
              XRMU=RMU(J)
              PL(-1)=0.D+00
              PL(0)=1.D+00
              DO 9 K=0,OS_NS
                PL(K+1)=((2*K+1.)*XRMU*PL(K)-K*PL(K-1))/(K+1.)
                IT=IT+BETA(K)*PL(K)
    9         CONTINUE
    
	       WRITE(99,203,err=921) ' RMU ',RMU(J),' R11 ',R11(J),
     &	                            ' R11 recomp  ',IT         
    8      CONTINUE
      END IF !Trace file
      
C* Computation of the parameters GAMMA(k) and DELTA(k) of the decomposition of the functions
C* R12 and R33 in Legendre polynomials or functions.
C-----------------------------------------------------------------------------
      DO 11 J=-LUM_NBMU,LUM_NBMU
         IF(J.EQ.0)  GO TO 11
         XXX=CHR(J)*R12(J)
         XX=CHR(J)*R33(J)
         POL(0)=0.D+00
         POL(1)=0.D+00
         XRMU=RMU(J)
         PL(-1)=0.D+00
         PL(0)=1.D+00
         POL(2)=3.*(1.-XRMU**2)/2./DSQRT(6.0D+00)
	
         DO 12 K=2,OS_NS
            D=(2.*K+1.)/DSQRT(1.D+00*(K+3.)*(K-1.))
            E=DSQRT(1.D+00*(K+2.)*(K-2.))/(2.*K+1.)
            POL(K+1)=D*(XRMU*POL(K)-E*POL(K-1))
            GAMMA(K)=GAMMA(K)+XXX*POL(K)
  12     CONTINUE

         DO 14 K=0,OS_NS
            PL(K+1)=((2.*K+1.)*XRMU*PL(K)-K*PL(K-1))/(K+1.)
            DELTA(K)=DELTA(K)+XX*PL(K)
  14     CONTINUE

   11 CONTINUE

      DO 15 K=0,OS_NS
         DELTA(K)=DELTA(K)*(2.*K+1.)*.5
         GAMMA(K)=GAMMA(K)*(2.*K+1.)*.5
   15 CONTINUE



C* Trace file
C-----------------
      IF (TRACE) THEN
          WRITE(99,*,err=921) '  '
	   WRITE(99,*,err=921) ' Recomposition of the R12 function'
      
C* Check of the recombination of R12
C--------------------------------------
           DO 81 J=-LUM_NBMU,LUM_NBMU
                 QT=0.D+00
	      
	   	   IF(J.EQ.0)  GO TO 81
                 POL(0)=0.D+00
                 POL(1)=0.D+00
                 XRMU=RMU(J)
                 POL(2)=3.*(1.-XRMU**2)/2./DSQRT(6.0D+00)
   
                 DO 91 K=2,OS_NS
		      D=(2.*K+1.)/DSQRT(1.D+00*(K+3.)*(K-1.))
                    E=DSQRT(1.D+00*(K+2.)*(K-2.))/(2.*K+1.)
		      POL(K+1)=D*(XRMU*POL(K)-E*POL(K-1))
		      QT=QT+GAMMA(K)*POL(K)
  91             CONTINUE
    
	          WRITE(99,203,err=921) ' RMU ',RMU(J),' R12 ',R12(J),
     &	                            ' R12 recomp  ',QT 
  81       CONTINUE
  
  
           WRITE(99,*,err=921) '  '
	   WRITE(99,*,err=921) ' Recomposition of the R33 function'
      
C* Check of the recombination of R33
C--------------------------------------
           DO 82 J=-LUM_NBMU,LUM_NBMU
	      UT=0.D+00
	      
	      IF(J.EQ.0)  GO TO 82
                 PL(-1)=0.D+00
                 PL(0)=1.D+00
		 XRMU=RMU(J)

                 DO 92 K=0,OS_NS
                    PL(K+1)=((2*K+1.)*XRMU*PL(K)-K*PL(K-1))/(K+1.)
                    UT=UT+DELTA(K)*PL(K)
  92             CONTINUE
    
	         WRITE(99,203,err=921) ' RMU ',RMU(J),' R33 ',R33(J),
     &                              ' R33 recomp  ',UT   
  82       CONTINUE
      END IF !Trace file
      
      
      
C* Computation of the parameters ALPHA(k) and ZETA(k) from BETA(k) and DELTA(k).
C-------------------------------------------------------------------------------
      DO 16 I=2,OS_NS
         CO1=4*(2*I+1.)/FLOAT(I)/(I-1.)/(I+1.)/(I+2.)
         CO2=I*(I-1.)/((I+1.)*(I+2.))
         CO3=CO2*DELTA(I)
         CO2=CO2*BETA(I)
         NN=I*.5
         MM=(I-1)*.5
         SOM1=0.D+00
         SOM2=0.D+00
         SOM3=0.D+00
         SOM4=0.D+00
         DO 17 J=1,NN
            X2=(I-1.)*(I-1.)-3.*(2*J-1.)*(I-J)
            SOM1=SOM1+X2*BETA(I-2*J)
            SOM2=SOM2+X2*DELTA(I-2*J)
  17     CONTINUE

         DO 18 J=0,MM
            X2=(I-1.)*(I-1.)-3.*J*(2*I-2*J-1.)
            SOM3=SOM3+X2*BETA(I-2*J-1)
            SOM4=SOM4+X2*DELTA(I-2*J-1)
  18     CONTINUE

         ZETA(I)=CO3-CO1*(SOM2-SOM3)
         ALPHA(I)=CO2-CO1*(SOM1-SOM4)
  16  CONTINUE


C*   Recording of the results
C---------------------------------
      DO 20 K=0,OS_NS
         WRITE(3,207,err=992)ALPHA(K),BETA(K),GAMMA(K),ZETA(K)
   20 CONTINUE



C* Trace file
C-----------------
      IF (TRACE) THEN
           WRITE(99,*,err=921) '  '
	   WRITE(99,*,err=921) 'Limit order : ',OS_NS
	   WRITE(99,215,err=921)'K','ALPHA(K)','BETA(K)','GAMMA(K)',
     &                      'DELTA(I)','ZETA(K)'  
           DO 45 K=0,OS_NS
	     WRITE(99,205,err=921) K,ALPHA(K),BETA(K),GAMMA(K),DELTA(K),
     &	                           ZETA(K)
   45      CONTINUE
      END IF 
   
   
       
C* Nominal end
C-------------------
      IER=0
      goto 9999
      
C* Error case: return of status -1 to the calling program
C------------------------------------------------------
  991 WRITE(6,*) '  ERROR on a file opening for SOS_MAT_FRESNEL'
      IER=-1
      goto 9999
  921 WRITE(6,*) '  ERROR on logfile writing for SOS_MAT_FRESNEL'
      IER=-1
      goto 9999
  992 WRITE(6,*) '  ERROR on a file reading/writing for SOS_MAT_FRESNEL'
      IER=-1 
      goto 9999

 9999 CONTINUE

C* File closing
C-------------------
      CLOSE(3)
   
C* Format
C---------   
  207 FORMAT(4(E15.8))
  200 FORMAT(A17,F6.2,3(A5,E13.7))
  203 FORMAT(2(A5,E13.7),A12,E13.7)
  205 FORMAT(1X,I3,5(2X,E13.5))
  215 FORMAT(3X,A1,7X,A8,7X,A7,8X,A8,7X,A8,7X,A7)    
      RETURN
      END         !END OF PROCEDURE SOS_MAT_FRESNEL




C==============================================================================
C PROCEDURE: SOS_MAT_REFLEXION
C ==========
C      This procedure computes the surface reflection matrix
C      for the glitter.
C      The incidence and reflection angles are expressed as being
C      below 90 degrees (Gauss mu > 0).
C
C
C
C      NB: IMPLEMENTATION APPROACH:
C      *******************************
C          The implementation of the code is adapted to a computation of the
C          reflection matrix (routine SOS_MAT_REFLEXION) using only positive
C           values of cos(Gauss angle). To achieve this:
C          - The complementary angle of the true scattering angle is used
C            (2*reflection angle): DEM 1
C            (routine SOS_MAT_FRESNEL).
C          - The opposite of the true parallel reflection coefficient
C            is used: DEM 2 (routine SOS_MAT_FRESNEL).
C
C          This routine computing the reflection matrix introduces
C          a coefficient (-1)**K in the computation of the terms of the
C          Fourier series to make the whole consistent: DEM 3.
C
C
C Description of parameters
C --------------------------
C     COEF  (double): (I) Multiplicative coefficient of the reflection matrix.
C      LUM_NBMU (I4) : (I) Number of actually useful (positive) angles.
C     RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)   (double)    (I) Array of the Gauss mu
C                             (complementary angles of the scattering angles).
C             OS_NB (I4) : (I) Order of the expansions of the phase functions
C		           in Legendre polynomials to be used.
C       OS_NS (I4) : (I) Maximum order to be used for the expansion in
C                        Legendre polynomials of the Fresnel matrix elements
C                        and in Fourier series of the radiance computations.
C       OS_NM (I4) : (I) Maximum order of the expansion in Fourier series
C                                of the G function of the reflection matrix computations.
C      FIC_RES_GSF  (CHARACTER)  : (I) Name of the temporary file RES_GSF
C                                      with complete directory tree.
C     FIC_RES_FRESNEL (CHARACTER)  : (I) Name of the temporary file RES_FRESNEL
C                                      with complete directory tree.
C      FIC_RES_MAT_REFLEX (CHARACTER)  : (O) Name of the temporary file RES_MAT_REFLEX
C                                      with complete directory tree.
C       IER  (I4)   (O) : error code =0 if no error, =-1 otherwise
C
C
C 
C Description of the files used
C ---------------------------------
C   -->Content of the temporary file RES_GSF
C       Unformatted binary file, created by SOS_GSF
C       File created in the directory where the program is run
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
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C     -Opening of a file
C      The file RES_GSF or the file RES_FRESNEL does not exist
C     -Reading/writing in a file
C
C     Message displayed on screen, the routine stops its computations and
C     returns status -1 to the calling program
C
C
C     
C  Routines used by the procedure:
C -------------------------------------
C  This program calls the routines:
C         -SOS_NOYAUX_FRESNEL
C
C=============================================================================
      SUBROUTINE SOS_MAT_REFLEXION(COEF,LUM_NBMU,RMU,OS_NB,OS_NS,OS_NM,
     &                             FIC_RES_FRESNEL, FIC_RES_GSF,
     &                             FIC_RES_MAT_REFLEX, IER)

      IMPLICIT NONE

C* Definition of variables
C*--------------------------------------------------------------------
      DOUBLE PRECISION COEF !Multiplicative coefficient of the reflection matrix.
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Array of the Gauss mu
				!(complementary angles of the scattering angles)

      DOUBLE PRECISION ALPHA(0:CTE_OS_NS_MAX)   ! Parameters alpha(k).
      DOUBLE PRECISION BETA(0:CTE_OS_NS_MAX)    ! Parameters beta(k).
      DOUBLE PRECISION GAMMA(0:CTE_OS_NS_MAX)   ! Parameters gamma(k).
      DOUBLE PRECISION ZETA(0:CTE_OS_NS_MAX)    ! Parameters zeta(k).

      DOUBLE PRECISION G(0:CTE_OS_NM_MAX)	! G function decomposed in Fourier
      					        ! series.

       ! Element Pij of the Fresnel reflection matrix decomposed in Fourier
	! series for an order IS = 0 to OS_NS.
	! Computation for the incidence direction Teta1 (RMU1)
	! and the reflection direction Teta2 (RMU2).
      DOUBLE PRECISION BP(0:CTE_OS_NS_MAX,2)    !   Element P11.
      DOUBLE PRECISION GR(0:CTE_OS_NS_MAX,2)    !   Element P12 or P21.
      DOUBLE PRECISION GT(0:CTE_OS_NS_MAX,2)    !   Element -P13 or -P31.
      DOUBLE PRECISION ARR(0:CTE_OS_NS_MAX,2)   !   Element P22.
      DOUBLE PRECISION ART(0:CTE_OS_NS_MAX,2)   !   Element -P23 or -P32.
      DOUBLE PRECISION ATT(0:CTE_OS_NS_MAX,2)   !   Element P33.

      DOUBLE PRECISION R111    ! Variable for the computation of the term M11(IS,1).
      DOUBLE PRECISION R121    ! Variable for the computation of the term M12(IS,1).
      DOUBLE PRECISION R122    ! Variable for the computation of the term M12(IS,2).
      DOUBLE PRECISION R131    ! Variable for the computation of the term M13(IS,1).
      DOUBLE PRECISION R132    ! Variable for the computation of the term M13(IS,2).
      DOUBLE PRECISION R211    ! Variable for the computation of the term M21(IS,1).
      DOUBLE PRECISION R212    ! Variable for the computation of the term M21(IS,2).
      DOUBLE PRECISION R221    ! Variable for the computation of the term M22(IS,1).
      DOUBLE PRECISION R222    ! Variable for the computation of the term M22(IS,2).
      DOUBLE PRECISION R231    ! Variable for the computation of the term M23(IS,1).
      DOUBLE PRECISION R232    ! Variable for the computation of the term M23(IS,2).
      DOUBLE PRECISION R311    ! Variable for the computation of the term M31(IS,1).
      DOUBLE PRECISION R312    ! Variable for the computation of the term M31(IS,2).
      DOUBLE PRECISION R321    ! Variable for the computation of the term M32(IS,1).
      DOUBLE PRECISION R322    ! Variable for the computation of the term M32(IS,2).
      DOUBLE PRECISION R331    ! Variable for the computation of the term M33(IS,1).
      DOUBLE PRECISION R332    ! Variable for the computation of the term M33(IS,2).

      DOUBLE PRECISION X,Y	! Multiplicative coefficients of the reflection
      				! matrix.

      	! Elements Mij of the reflection matrix versus the order IS
	! of the Fourier series, for 2 directions.
	! The direction index is:
	!     _ 1 for an incidence at Teta1 and a reflection at Teta2,
	!     _ 2 for an incidence at Teta2 and a reflection at Teta1.
      REAL M11(0:CTE_OS_NB_MAX,2)		
      REAL M12(0:CTE_OS_NB_MAX,2)	
      REAL M13(0:CTE_OS_NB_MAX,2)		
      REAL M21(0:CTE_OS_NB_MAX,2)	
      REAL M22(0:CTE_OS_NB_MAX,2)		
      REAL M23(0:CTE_OS_NB_MAX,2)		
      REAL M31(0:CTE_OS_NB_MAX,2)	
      REAL M32(0:CTE_OS_NB_MAX,2)		
      REAL M33(0:CTE_OS_NB_MAX,2)	

      INTEGER*4 LUM_NBMU ! Number of actual angles of the simulations
      
      INTEGER*4 OS_NB ! Order of the expansions of the phase functions
			 ! in Legendre polynomials.
      INTEGER*4 OS_NS ! Maximum order of the expansion in Legendre polynomials
			 ! of the Fresnel matrix elements.
      INTEGER*4 OS_NM ! Maximum order of the expansion in Fourier series
			 ! of the G function of the reflection matrix computations.
 
      INTEGER*4 I	! Index of the angle Teta1.
      INTEGER*4 J	! Index of the angle Teta2.

      INTEGER*4 L	! Order of the decomposition in Legendre polynomials.
      INTEGER*4 LIM	! Maximum order of the decomposition in Legendre
      			! polynomials.
      INTEGER*4 II,II2	! Indices stored with the function G(IS).
      INTEGER*4 IJ	
      INTEGER*4 IS	! Order of the Fourier series decomposition.
      INTEGER*4 K
      INTEGER*4 I1	! Value K+IS
      INTEGER*4 I2	! Value ABS(K-IS)
      INTEGER*4 IM	! Value 1 or -1.
      INTEGER*4 IER  ! error code =0 if no error, =-1 otherwise
  
      CHARACTER*CTE_LENFIC2 FIC_RES_MAT_REFLEX 
                     ! Name of the file RES_MAT_REFLEX
                     ! with complete directory tree

      CHARACTER*CTE_LENFIC2 FIC_RES_GSF ! Name of the file RES_GSF
                                        ! with complete directory tree

      CHARACTER*CTE_LENFIC2 FIC_RES_FRESNEL ! Name of the file RES_FRESNEL
                                            ! with complete directory tree

      
C* Creation of the result file of SOS_MAT_REFLEXION
C---------------------------------------------------------
      OPEN(UNIT=4,FILE=FIC_RES_MAT_REFLEX,FORM='UNFORMATTED',ERR=991)


C* Reading of the parameters of the decomposition of the Fresnel matrix
C* in Legendre polynomials or functions.
C----------------------------------------------------------------------
      OPEN(UNIT=2,FILE=FIC_RES_FRESNEL,STATUS='OLD',ERR=991)
      
      DO 951 L=0,OS_NS
         READ(2,1300,err=992) ALPHA(L),BETA(L),GAMMA(L),ZETA(L)
  951 CONTINUE
  
      CLOSE(2)
      
C* Opening of the file RES_GSF
C---------------------------------------------------------
      OPEN(UNIT=3,FILE=FIC_RES_GSF,FORM='UNFORMATTED',
     &     STATUS='OLD',ERR=991)

C* Loop on the incidence and reflection directions.
C--------------------------------------------------------
      DO 25 I=1,LUM_NBMU   !Teta1 : incidence zenith angle
                           !(0 < Teta1 < PI/2),
         DO 2 J=1,I        !Teta2 : reflection zenith angle
	                    !(0 < Teta2 < PI/2),

C* Reading of the Fourier series decomposition of the function G
C* for the angles Teta1 and Teta2, with a maximum useful order L.
C-------------------------------------------------------------------
            READ(3,err=992)II,II2,LIM,(G(K),K=0,LIM)

C* The array is completed with 0.
C---------------------------------
            DO IJ=LIM+1,OS_NM
               G(IJ)=0.
            ENDDO
	      
C------------------------------------------------------------------------
C* Fourier series decomposition of the Fresnel reflection matrix
C* for the directions Teta1 and Teta2
C--------------------------------------------------------------------------

C*  Computation of the functions BP,GR,GT,ARR,ART and ATT versus the order IS
C*  for the directions Teta1 and Teta2
C---------------------------------------------------------------------------
            CALL SOS_NOYAUX_FRESNEL(RMU(I),RMU(J),OS_NS,
     &                              ALPHA,BETA,GAMMA,ZETA,
     &                              BP,GR,GT,ARR,ART,ATT)

C* Computation of the elements of the reflection matrix for each order IS.
C---------------------------------------------------------------------
            DO 24 IS=0,OS_NB
               X=COEF*G(IS)/4.
	       
C*   Elements for IS=0.
               R111     =X*BP(0,1)
               R121     =X*GR(0,1)
               R122     =X*GR(0,2)
               R131     =0.
               R132     =0.
               R231     =0.
               R232     =0.
               R211     =X*GR(0,2)
               R212     =X*GR(0,1)
               R221     =X*ARR(0,2)
               R222     =X*ARR(0,1)
               R311     =0.
               R312     =0.
               R321     =0.
               R322     =0.
               R331     =X*ATT(0,2)
               R332     =X*ATT(0,1)

C*   Elements for the summation from K=1 to OS_NS
               IM=1
               DO 3 K=1,OS_NS
                 IM=-IM					! DEM 3
                 I1=K+IS
                 I2=ABS(K-IS)
                 IF ((I1.GT.LIM).AND.(I2.GT.LIM)) GOTO 3	
                 X=COEF*IM*(G(I1)+G(I2))/4.		! DEM 3
                 Y=COEF*IM*(G(I2)-G(I1))/4.		! DEM 3
                 R111     =R111     +BP(K,1)*X
                 R121     =R121     +GR(K,1)*X
                 R122     =R122     +GR(K,2)*X
                 R131     =R131     +GT(K,1)*Y
                 R132     =R132     +GT(K,2)*Y
                 R211     =R211     +GR(K,2)*X
                 R212     =R212     +GR(K,1)*X
                 R221     =R221     +ARR(K,2)*X
                 R222     =R222     +ARR(K,1)*X
                 R231     =R231     +ART(K,2)*Y
                 R232     =R232     +ART(K,1)*Y
                 R311     =R311     +GT(K,2)*Y
                 R312     =R312     +GT(K,1)*Y
                 R321     =R321     +ART(K,1)*Y
                 R322     =R322     +ART(K,2)*Y
                 R331     =R331     +ATT(K,2)*X
                 R332     =R332     +ATT(K,1)*X
  3            CONTINUE

C* Assignment of the results to the result arrays Rij(IS)
 244           M11(IS,1)=R111
               M11(IS,2)=R111
               M12(IS,1)=R121
               M12(IS,2)=R122
               M13(IS,1)=R131
               M13(IS,2)=R132
               M21(IS,1)=R211
               M21(IS,2)=R212
               M22(IS,1)=R221
               M22(IS,2)=R222
               M23(IS,1)=R231
               M23(IS,2)=R232
               M31(IS,1)=-R311
               M31(IS,2)=-R312
               M32(IS,1)=-R321
               M32(IS,2)=-R322
               M33(IS,1)=-R331
               M33(IS,2)=-R332
 24         CONTINUE	!End of loop on IS

C* Recording of the order for the pair (Teta1, Teta2) of
C* the Fourier series decomposition of the Fresnel matrix.
C-----------------------------------------------------------------
            WRITE(4,err=992) M11,M12,M13,M21,M22,M23,M31,M32,M33

 2       CONTINUE	!End of loop on Teta2 (J)

 25   CONTINUE	!End of loop on Teta1 (I)
 
 
C* Nominal end
C-------------------
      IER=0
      goto 9999
      
C* Error case: return of status -1 to the calling program
C----------------------------------------------------------    
  991 WRITE(6,*) '  ERROR on a file opening for SOS_MAT_REFLEXION'
      IER=-1
      goto 9999
  992 WRITE(6,*) '  ERROR on a file reading/writing for ',
     &           'SOS_MAT_REFLEXION '
      IER=-1 
      goto 9999

 9999 CONTINUE

C* File closing
C-------------------
      CLOSE(2)
      CLOSE(3)
      CLOSE(4)
 
C* Format
C---------   
 1300 FORMAT(4(E15.8))  
   
      RETURN
      END         !END OF PROCEDURE SOS_MAT_REFLEXION





C==============================================================================
C PROCEDURE: SOS_NOYAUX_FRESNEL
C ==========
C      This procedure computes the Fourier series decomposition of
C      the Fresnel matrix versus the azimuth, for the directions Teta1 and Teta2.
C
C      The following matrix elements are computed:
C         BP(IS,mu)    GR(IS,mu)    GT(IS,mu)
C         ARR(IS,mu)   ART(IS,mu)   ATT(IS,mu)
C      for the orders IS of the Fourier series decomposition
C      and for mu = cos(Teta1) or cos(Teta2).
C
C
C Description of parameters
C --------------------------
C	RMU1   (double)	 (I) Cos(Teta1).
C	RMU2   (double)	 (I) Cos(Teta2).
C       OS_NS (I4)           (I) Maximum order to be used for the expansion in
C                                Legendre polynomials of the elements of the Fresnel
C                                matrix and in Fourier series of the radiance
C                                computations.
C
C	ALPHA(0:CTE_OS_NS_MAX)  (double)	(I) Parameters alpha(k).
C	BETA(0:CTE_OS_NS_MAX)   (double)	(I) Parameters beta(k).
C	GAMMA(0:CTE_OS_NS_MAX)  (double)	(I) Parameters gamma(k).
C	ZETA(0:CTE_OS_NS_MAX)   (double)	(I) Parameters zeta(k).
C
C       Element Pij of the Fresnel reflection matrix decomposed in Fourier
C       series for an order IS = 0 to OS_NS.
C       Computation for the incidence direction Teta1 (RMU1)
C       and the reflection direction Teta2 (RMU2).
C	BP(0:CTE_OS_NS_MAX,2)   (double)	(O) Element P11.
C	GR(0:CTE_OS_NS_MAX,2)   (double)	(O) Element P12 or P21.
C	GT(0:CTE_OS_NS_MAX,2)   (double)	(O) Element -P13 or -P31.
C	ARR(0:CTE_OS_NS_MAX,2)  (double)	(O) Element P22.
C	ART(0:CTE_OS_NS_MAX,2)  (double)	(O) Element -P23 or -P32.
C	ATT(0:CTE_OS_NS_MAX,2)  (double)	(O) Element P33.
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
C
C=============================================================================
      SUBROUTINE SOS_NOYAUX_FRESNEL(RMU1,RMU2,OS_NS,
     &                              ALPHA,BETA,GAMMA,ZETA,
     &                              BP,GR,GT,ARR,ART,ATT)

      IMPLICIT NONE

C* Definition of variables
C*----------------------------------------------------------------------
      DOUBLE PRECISION RMU1	! Cos(Teta1).
      DOUBLE PRECISION RMU2	! Cos(Teta2).

      DOUBLE PRECISION ALPHA(0:CTE_OS_NS_MAX)    ! Parameters alpha(k).
      DOUBLE PRECISION BETA(0:CTE_OS_NS_MAX)	! Parameters beta(k).
      DOUBLE PRECISION GAMMA(0:CTE_OS_NS_MAX)    ! Parameters gamma(k).
      DOUBLE PRECISION ZETA(0:CTE_OS_NS_MAX)	! Parameters zeta(k).

      DOUBLE PRECISION R(2)	! array containing cos(Teta1) and cos(Teta2).

      	! Element Pij of the Fresnel reflection matrix decomposed in Fourier
	! series for an order IS = 0 to OS_NS.
	! Computation for the incidence direction Teta1 (RMU1)
	! and the reflection direction Teta2 (RMU2).
      DOUBLE PRECISION BP(0:CTE_OS_NS_MAX,2)    ! Element P11.
      DOUBLE PRECISION GR(0:CTE_OS_NS_MAX,2)    ! Element P12 or P21.
      DOUBLE PRECISION GT(0:CTE_OS_NS_MAX,2)    ! Element -P13 or -P31.
      DOUBLE PRECISION ARR(0:CTE_OS_NS_MAX,2)   ! Element P22.
      DOUBLE PRECISION ART(0:CTE_OS_NS_MAX,2)   ! Element -P23 or -P32.
      DOUBLE PRECISION ATT(0:CTE_OS_NS_MAX,2)   ! Element P33.

      DOUBLE PRECISION PSL(-1:CTE_OS_NS_MAX,2)  
                            ! Legendre function PSL for the order S=IS
      				! versus the order L and the angle mu.
      DOUBLE PRECISION RSL(-1:CTE_OS_NS_MAX,2)  ! Function RSL.
      DOUBLE PRECISION TSL(-1:CTE_OS_NS_MAX,2)  ! Function TSL.

      DOUBLE PRECISION C	     ! Single value of the Gauss mu.
      DOUBLE PRECISION A,B,D,E,F ! Coefficients for the recurrence
      				     ! on PSL, RSL and TSL.

      DOUBLE PRECISION X	! Intermediate computation variable.
      DOUBLE PRECISION XX	! Intermediate computation variable.
      DOUBLE PRECISION YY	! Intermediate computation variable.

      DOUBLE PRECISION RAC3	! Square root of 3.
      DOUBLE PRECISION X26	! 2*sqrt(6)

      INTEGER*4 OS_NS   ! Maximum order of the expansion in Legendre polynomials
			   ! of the Fresnel matrix elements.

      INTEGER*4 IS	! Order of the Fourier series decomposition.
      INTEGER*4 I
      INTEGER*4 J
      INTEGER*4 K
      INTEGER*4 L


C* Initialisation
C-------------------------
      RAC3=DSQRT(3.D+00)
      X26=2.*DSQRT(6.D+00)

C* Storage of cos(Teta1) and cos(Teta2) in the array R.
C------------------------------------------------------
      R(1)=RMU1
      R(2)=RMU2

C* Fourier series decomposition from IS = 0 to OS_NS
C----------------------------------------------------------
      DO 1 IS=0,OS_NS

C-----------------------------------------------------
C* Initialisation of the functions of order IS for the
C* orders L = 0, 1 and 2.
C-----------------------------------------------------

C* Initialisation for the order IS = 0
C-------------------------------------
         IF (IS.NE.0) GOTO 700
         DO 25 J=1,2
            C=R(J)
            PSL(0,J)=1
            PSL(1,J)=C
            X=(3*C*C-1)*0.5
            PSL(2,J)=X
            RSL(1,J)=0
            X=3*(1-C*C)/X26
            RSL(2,J)=X
            TSL(1,J)=0.
            TSL(2,J)=0.
   25    CONTINUE
         GOTO 501

C* Initialisation for the order IS = 1
C-------------------------------------
  700    IF (IS.NE.1) GOTO 701
         DO 26 J=1,2
            C=R(J)
            X=1-C*C
            PSL(0,J)=0
            PSL(1,J)=DSQRT(X*0.5)
            PSL(2,J)=C*PSL(1,J)*RAC3
            TSL(1,J)=0.
            RSL(1,J)=0
            RSL(2,J)=-C*DSQRT(X)*0.5
            TSL(2,J)=-DSQRT(X)*0.5
   26    CONTINUE
         GOTO 501

C* Initialisation for the order IS >= 2
C----------------------------------------
  701    A=1
         DO I=1,IS
            X=I
            A=A*DSQRT((I+IS)/X)*0.5
         ENDDO
	 
         B=A*DSQRT(IS/(IS+1.D+00))*DSQRT((IS-1.D+00)/(IS+2.)) 
         DO J=1,2
             C=R(J)
             XX=1-C*C
             YY=IS*0.5
             PSL(IS-1,J)=0.
             RSL(IS-1,J)=0.
             TSL(IS-1,J)=0.
             X=XX**YY
             PSL(IS,J)=A*X
             YY=YY-1
             X=XX**YY
             RSL(IS,J) = B*(1+C*C)*X
             TSL(IS,J)= 2*B*C*X
         ENDDO

C* Computation of the functions PSL(L,mu), RSL(L,mu) and TSL(L,mu)
C* by recurrence on the order L for fixed IS.
C------------------------------------------------------
  501    K=2
         IF (IS.GT.2) K=IS

         DO 30 L=K,OS_NS-1
C* Coefficients for the recurrence on PSL
             A=(2*L+1.)/DSQRT((L+IS+1.D+00)*(L-IS+1.))
             B=DSQRT(DFLOAT((L+IS)*(L-IS)))/(2.*L+1.) 

C* Coefficients for the recurrence on RSL and TSL
             D=(L+1.)*(2*L+1.)/DSQRT((L+3.D+00)*(L-1.)*(L+IS+1.)
     &         *(L-IS+1.))
             E=DSQRT((L+2.D+00)*(L-2.)*(L+IS)*(L-IS))/(L*(2.*L+1.))
             F=2.*IS/(L*(L+1.))

C* Application of the recurrence for each Gauss mu
C------------------------------------------------------
             DO 31 J=1,2
                C=R(J)
                X=A*(C*PSL(L,J)-B*PSL(L-1,J))
                PSL(L+1,J)=X
                X=D*(C*RSL(L,J)-F*TSL(L,J)-E*RSL(L-1,J))
                RSL(L+1,J)=X
                X=D*(C*TSL(L,J)-F*RSL(L,J)-E*TSL(L-1,J))
                TSL(L+1,J)=X
   31        CONTINUE
   30    CONTINUE


C------------------------------------------------------
C* Computation of the elements of the reflection matrix
C------------------------------------------------------

C* Loops on the directions Teta1 and Teta2
C------------------------------------------------------
         DO K=1,2
             J=3-K

C* Initialisation of arrays
             BP(IS,K)=0.
             ARR(IS,K)=0.
             ATT(IS,K)=0
             GR(IS,K)=0.
             GT(IS,K)=0.
             ART(IS,K)=0

C* Summation loop on L = IS to OS_NS
C-----------------------------------------------
             DO L=IS,OS_NS
                BP(IS,K)=BP(IS,K)+BETA(L)*PSL(L,J)*PSL(L,K)
	         GR(IS,K)=GR(IS,K)+GAMMA(L)*PSL(L,J)*RSL(L,K)
                GT(IS,K)=GT(IS,K)+GAMMA(L)*PSL(L,J)*TSL(L,K)
                ATT(IS,K)=ATT(IS,K)+ALPHA(L)*TSL(L,J)*TSL(L,K)+
     &                    ZETA(L)*RSL(L,J)*RSL(L,K)
                ARR(IS,K)=ARR(IS,K)+ZETA(L)*TSL(L,J)*TSL(L,K)+
     &		          ALPHA(L)*RSL(L,J)*RSL(L,K)
                ART(IS,K)=ART(IS,K)+ALPHA(L)*RSL(L,K)*TSL(L,J)+
     &                    ZETA(L)*RSL(L,J)*TSL(L,K)
             ENDDO  ! loop on L
	     
         ENDDO  ! loop on K

   1  CONTINUE
      RETURN
      END         !END OF PROCEDURE SOS_NOYAUX_FRESNEL






C==============================================================================
C PROCEDURE:  SOS_MISE_FORMAT
C ==========
C      This procedure reorders the elements of the matrix Mij into a matrix
C      Pij usable by the OS code.
C
C      The elements Mij, stored in the file "RES_MAT_REFLEX", give
C      for each pair of angles (Teta1, Teta2) the associated Fourier series
C      decomposition.
C
C      The elements Pij, stored in the output file, give
C      for each order of the Fourier series decomposition
C      the reflection associated with all the pairs (Teta1, Teta2).
C
C Description of parameters
C --------------------------
C     FIC_RES_MAT_REFLEX (CHARACTER)  : (I) Name of the temporary file RES_MAT_REFLEX
C                                      with complete directory tree.
C      FICGLITTER  (CHARACTER*CTE_LENFIC2) : (O)  Complete name of the result file.
C      NBMU (I4)  : (I)  Number of actually useful (positive) angles.
C             OS_NB (I4) : (I) Order of the expansions of the phase functions
C		           in Legendre polynomials to be used.
C     IER (I4) : (O) error code =0 if no error, =-1 otherwise
C
C
C Other result provided by the procedure
C -------------------------------------
C      Result file: FICGLITTER, GLITTER file
C 
C Description of the files used
C ---------------------------------
C   -->Content of the temporary file RES_MAT_REFLEX
C       Unformatted binary file, created by SOS_MAT_REFLEXION
C       File created in the directory where the program is run
C       Record for each pair (Teta1, Teta2) of
C       the Fourier series decomposition of the Fresnel matrix:
C       M11,M12,M13,M21,M22,M23,M31,M32,M33 
C       (n real values)
C
C   -->Content of the GLITTER file
C       Unformatted binary file, created by SOS_MISE_FORMAT
C       It contains the surface reflection matrix
C       (GLITTER file or other BRDF):
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
C        --> Records on the domain limited to the useful angles
C                                   --------------------------------
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C     -Opening of a file
C      The file RES_MAT_REFLEX does not exist
C      The Glitter file exists
C     -Reading/writing in a file
C
C     Message displayed on screen, the routine stops its computations and
C     returns status -1 to the calling program
C
C
C=============================================================================
      SUBROUTINE SOS_MISE_FORMAT(FIC_RES_MAT_REFLEX, FICGLITTER,
     &                           NBMU,OS_NB,IER)

      IMPLICIT NONE

C* Definition of variables
C*--------------------------------------------------------------------
      CHARACTER*CTE_LENFIC2 FICGLITTER	! Name of the result file.

      CHARACTER*CTE_LENFIC2 FIC_RES_MAT_REFLEX ! Name of the file RES_MAT_REFLEX
                                               ! with complete directory tree

      	! Elements Mij of the reflection matrix versus the order IS
	! of the Fourier series, for 2 directions.
	! The direction index is:
	!     _ 1 for an incidence at Teta1 and a reflection at Teta2,
	!     _ 2 for an incidence at Teta2 and a reflection at Teta1.
      REAL M11(0:CTE_OS_NB_MAX,2)		
      REAL M12(0:CTE_OS_NB_MAX,2)		
      REAL M13(0:CTE_OS_NB_MAX,2)			
      REAL M21(0:CTE_OS_NB_MAX,2)	
      REAL M22(0:CTE_OS_NB_MAX,2)		
      REAL M23(0:CTE_OS_NB_MAX,2)		
      REAL M31(0:CTE_OS_NB_MAX,2)	
      REAL M32(0:CTE_OS_NB_MAX,2)	
      REAL M33(0:CTE_OS_NB_MAX,2)

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

      INTEGER*4 NBMU     ! Number of actual angles of the simulations
      INTEGER*4 OS_NB    ! Order of the expansions of the phase functions
			    ! in Legendre polynomials.
      INTEGER*4 IS	! Order of the Fourier series decomposition.
      INTEGER*4 I	! Index of the incidence angle Teta1.
      INTEGER*4 J	! Index of the reflection angle Teta2.
      INTEGER*4 IER  ! error code =0 if no error, =-1 otherwise
	
C* File opening
C------------------------------------------------------------	
	
      OPEN(UNIT=1,FILE=FIC_RES_MAT_REFLEX,FORM='UNFORMATTED',
     &     STATUS='OLD',ERR=991)
      OPEN(UNIT=2,FILE=FICGLITTER,FORM='UNFORMATTED',
     &     STATUS='NEW',ERR=991)

C* Loop on the order of the Fourier series decomposition
C------------------------------------------------------------
      DO 1 IS=0,OS_NB

C* Loop on the incidence directions Teta1
C---------------------------------------------
         DO 2 I=1,NBMU

C* Loop on the reflection directions Teta2
C-----------------------------------------------

C*       Reading of the file RES_MAT_REFLEX stored for the max size
C*       of the arrays: CTE_OS_NBMU_MAX
            DO 3 J=1,I
               READ(1,err=992) M11,M12,M13,M21,M22,M23,M31,M32,M33
               P11(I,J)=M11(IS,1)
               P11(J,I)=M11(IS,2)
               P12(I,J)=M12(IS,1)
               P12(J,I)=M12(IS,2)
               P13(I,J)=M13(IS,1)
               P13(J,I)=M13(IS,2)
               P21(I,J)=M21(IS,1)
               P21(J,I)=M21(IS,2)
               P22(I,J)=M22(IS,1)
               P22(J,I)=M22(IS,2)
               P23(I,J)=M23(IS,1)
               P23(J,I)=M23(IS,2)
               P31(I,J)=M31(IS,1)
               P31(J,I)=M31(IS,2)
               P32(I,J)=M32(IS,1)
               P32(J,I)=M32(IS,2)
               P33(I,J)=M33(IS,1)
               P33(J,I)=M33(IS,2)
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
      goto 9999
      
C* Error case: return of status -1 to the calling program
C------------------------------------------------------
  991 WRITE(6,*) '  ERROR on a file opening for SOS_MISE_FORMAT'
      IER=-1
      goto 9999
  992 WRITE(6,*) '  ERROR on a file reading/writing for SOS_MISE_FORMAT'
      IER=-1 
      goto 9999

 9999 CONTINUE

C* File closing
C-------------------
      CLOSE(1)
      CLOSE(2)
      
      RETURN
      END         !END OF PROCEDURE SOS_MISE_FORMAT

      
      
      
C==============================================================================
C PROCEDURE: SOS_BPDF_AJOUT_BRDF
C ==========
C      This procedure reads two BRDF or BPDF files formatted for the OS
C      and sums their terms to create a final file, sum of the first
C      two.
C
C Description of parameters
C --------------------------
C       FICSURF1 (CHARACTER*CTE_LENFIC2) (I) : name of the 1st surface reflection file,
C                                              in the OS format.
C                                             (complete name with directory).
C       FICSURF2 (CHARACTER*CTE_LENFIC2) (I) : name of the 2nd surface reflection file,
C                                              in the OS format.
C                                             (complete name with directory).
C       NBMU (I4)  (I) : number of actually useful (positive) angles.
C      OS_NB (I4) (I) : Order of the expansions of the phase functions
C		           in Legendre polynomials to be used.
C       FICS   (CHARACTER*CTE_LENFIC2)  (I)  : name of the result SURFACE file.
C                                             (complete name with directory).
C       IER  (I4)   (O) : error code =0 if no error, =-1 otherwise
C
C Description of the files used
C ---------------------------------
C   -->Content of the BRDF or BPDF file
C       Unformatted binary file.
C       It contains the surface reflection matrix (BRDF / BPDF file):
C       i.e. the elements Rij of the reflection matrix
C       for OS_NB orders of the Fourier series decomposition.
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
      SUBROUTINE SOS_BPDF_AJOUT_BRDF(FICSURF1,FICSURF2,NBMU,OS_NB,
     &                               FICS,IER)

      IMPLICIT NONE

C* Definition of variables
C*--------------------------------------------------------------------------
      CHARACTER*CTE_LENFIC2 FICSURF1 ! name of the 1st surface reflection
      					  ! file, in the OS format.
      CHARACTER*CTE_LENFIC2 FICSURF2 ! name of the 2nd surface reflection
      					  ! file, in the OS format.
      CHARACTER*CTE_LENFIC2 FICS     ! name of the result surface reflection
      					  ! file, in the OS format.
					  
      	! Elements Pij of the reflection matrix versus the directions
	! of incidence and reflection, for an order IS of the expansion
	! in Fourier series. OS input format for the 1st file.
      REAL P11_SURF1(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P12_SURF1(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P13_SURF1(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)	
      REAL P21_SURF1(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P22_SURF1(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P23_SURF1(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P31_SURF1(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P32_SURF1(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P33_SURF1(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
           
        ! Elements Pij of the reflection matrix versus the directions
	! of incidence and reflection, for an order IS of the expansion
	! in Fourier series. OS input format for the 2nd file.
      REAL P11_SURF2(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P12_SURF2(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P13_SURF2(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)	
      REAL P21_SURF2(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P22_SURF2(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P23_SURF2(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P31_SURF2(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P32_SURF2(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P33_SURF2(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)   
          
        ! Elements Pij of the reflection matrix versus the directions
	! of incidence and reflection, for an order IS of the expansion
	! in Fourier series. OS input format for the result file.
      REAL P11(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P12(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P13(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)	
      REAL P21(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P22(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P23(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P31(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P32(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)
      REAL P33(CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX)

      INTEGER*4 NBMU    ! Number of actual angles of the simulations
      INTEGER*4 OS_NB   ! Order of the expansions of the phase functions
			   ! in Legendre polynomials.

      INTEGER*4 IS	! Order of the Fourier series decomposition.
      INTEGER*4 I	! Index of the incidence angle Teta1.
      INTEGER*4 J	! Index of the reflection angle Teta2.
      INTEGER*4 IER  ! error code =0 if no error, =-1 otherwise
      	      
      

C* Opening of the files
C------------------------------------------------------------		    
      OPEN(UNIT=1,FILE=FICSURF1,FORM='UNFORMATTED',
     &     STATUS='OLD',ERR=991)
      OPEN(UNIT=2,FILE=FICSURF2,FORM='UNFORMATTED',
     &     STATUS='OLD',ERR=991)
      OPEN(UNIT=3,FILE=FICS,FORM='UNFORMATTED',
     &     STATUS='NEW',ERR=991)

	   
C* Loop on the order of the Fourier series decomposition
C------------------------------------------------------------      
      DO 1 IS=0,OS_NB

C*      Reading of the files that were written with
C*    array bounds limited to the useful
C*    dimension -NBMU:NBMU
	   READ(1,err=992)
     &	       ((P11_SURF1(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P12_SURF1(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P13_SURF1(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P21_SURF1(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P22_SURF1(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P23_SURF1(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P31_SURF1(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P32_SURF1(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P33_SURF1(I,J),I=1,NBMU),J=1,NBMU)  
                 
	   READ(2,err=992)
     &	       ((P11_SURF2(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P12_SURF2(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P13_SURF2(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P21_SURF2(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P22_SURF2(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P23_SURF2(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P31_SURF2(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P32_SURF2(I,J),I=1,NBMU),J=1,NBMU),
     &	       ((P33_SURF2(I,J),I=1,NBMU),J=1,NBMU) 

           
C* Loop on the incidence directions Teta1
C---------------------------------------------
         DO 2 I=1,NBMU

C* Loop on the reflection directions Teta2
C-----------------------------------------------
            DO 3 J=1,NBMU
               P11(I,J) = P11_SURF1(I,J) + P11_SURF2(I,J)  
               P12(I,J) = P12_SURF1(I,J) + P12_SURF2(I,J)  
               P13(I,J) = P13_SURF1(I,J) + P13_SURF2(I,J)  
               P21(I,J) = P21_SURF1(I,J) + P21_SURF2(I,J)
               P22(I,J) = P22_SURF1(I,J) + P22_SURF2(I,J)
               P23(I,J) = P23_SURF1(I,J) + P23_SURF2(I,J)  
               P31(I,J) = P31_SURF1(I,J) + P31_SURF2(I,J)  
               P32(I,J) = P32_SURF1(I,J) + P32_SURF2(I,J)  
               P33(I,J) = P33_SURF1(I,J) + P33_SURF2(I,J)  
 3          CONTINUE
 2       CONTINUE
 
C*       Recording of the arrays of the reflection matrix
C*       at order IS of the Fourier series decomposition.
C*       Record bounds of the arrays limited to the
C*       useful dimension -NBMU:NBMU
C-------------------------------------------------------------- 
	 WRITE(3,err=992) ((P11(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P12(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P13(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P21(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P22(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P23(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P31(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P32(I,J),I=1,NBMU),J=1,NBMU),
     &	                  ((P33(I,J),I=1,NBMU),J=1,NBMU)
     
 1    CONTINUE


C* Nominal end
C-------------------
      IER=0
      GOTO 9999
 
 
C* Error case: return of status -1 to the calling program
C----------------------------------------------------------
  991 WRITE(6,*) '  ERROR on a file opening for SOS_BPDF_AJOUT_BRDF'
      IER=-1
      goto 9999
  992 WRITE(6,*) '  ERROR on a file reading/writing for ',
     &           'SOS_BPDF_AJOUT_BRDF'
      IER=-1 
      goto 9999

 9999 CONTINUE

C* File closing
C-------------------
      CLOSE(1)
      CLOSE(2)
      CLOSE(3)
      
      RETURN            
      END      !END OF PROCEDURE SOS_BPDF_AJOUT_BRDF
      
