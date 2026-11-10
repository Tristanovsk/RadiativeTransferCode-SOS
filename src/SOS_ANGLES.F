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
C* FILE: SOS_ANGLES.F
C* PROJECT: Successive orders of scattering
C* ROLE: Definition of the angles used to discretise the simulated
C*       radiative properties: phase functions, BRDF/BPDF, radiance fields.
C*       - For the phase function computations.
C*       - For the BRDF/BPDF matrix and radiance field computations.
C*       ==> Estimates the Gauss angles and associated weights for spatial integrations.
C*       ==> Adds user angles with null weights (no impact on the integrated values).
C*       ==> Adds the solar zenith angle (with null weight).
C*
C* AUTHOR:
C*   CS GROUP France, B.Lafrance, 

C* MOD:VERSION:1.0: 11/01/2010 Creation
C*
C* MOD:VERSION:1.1: 24/06/2010 
C*    - Constant CTE_DEFAULT_NBMU_LUM (from SOS.h) added as default value
C*      of the number of Gauss angles NBMU_GAUSS when it is not given
C*      by the user
C*    - FORMAT labels ,X, changed to ,1X, to avoid compilation errors.
C*
C* MOD:VERSION:6.0: 06/06/2013
C*
C*   -  The input parameters are given as a list of parameters
C*      and no longer as pairs with associated keyword => all
C*	 parameters are therefore mandatory
C*
C*   -  Change of the value of Pi, declared as a constant in the file
C*      (INCTE_PI replaces SOS_PI)
C*
C* MOD:VERSION:6.1: 31/08/2022
C*    - Constants renamed from SOS_* to CTE_*
C*    - Clarification of the result file: column OUTPUT renamed USER_ANGLE
C*    - Traces converted to English.
C*
C* MOD:VERSION:6.2: 26/09/2023
C*    - Change of the syntax of the test on XDD in routine SOS_GAUSS for
C*      compatibility with Fortran 2018 for the gcc/10.2.0 compiler
C*    
C*    
C******************************************************************************

C----------------------------------------------------------------------------
C Definition of constants
C---------------------------------------------------------------------------- 
C Constants used:
C    CTE_LENFIC2 : Length of file names with directory tree.
C    CTE_LENKEYWORD : Maximum size of the keywords for argument passing.
C    CTE_NBMAX_USER_ANGLES : Maximum number of user angles
C    INCTE_PI : Value of PI
C    CTE_MIE_NBMU_MAX : Maximum number of positive Gauss angles
C                       for MIE computations
C    CTE_OS_NBMU_MAX  : Maximum number of positive Gauss angles
C                       for radiance computations
C    CTE_OS_NB_MAX :    Maximum value of the limit order of the expansion in
C                       Legendre polynomials.
C    CTE_OS_NS_MAX :    Maximum value of the limit order of the expansion in Legendre
C                       polynomials of the Fresnel matrix elements
C                       and in Fourier series of the radiances.
C    CTE_OS_NM_MAX :    Maximum value of the limit order of the expansion
C                       in Fourier series of the G function.
C    CTE_NBANGLES_MAX : Maximum number of positive angles for the
C                       declaration of the arrays XMU and XHR of the routine
C    CTE_DEFAULT_NBMU_MIE : Default number of positive Gauss angles for the
C                           phase function computations (if not given by the user).
C    CTE_DEFAULT_OS_NB : Default order of the expansions of the phase functions
C                        in Legendre polynomials.
C    CTE_DEFAULT_NBMU_LUM : Default number of positive Gauss angles for the
C                           radiance computations (if not given by the user).
C    CTE_DEFAULT_OS_NS    : Default maximum order of the expansion in Legendre polynomials
C                           of the Fresnel matrix elements.
C    CTE_DEFAULT_OS_NM    : Default maximum order of the expansion in Fourier series
C                           of the G function of the reflection matrix computations.
C    CTE_SEUIL_ECART_MUS  : Threshold for the test of the difference between cos(teta(j)) and cos(tetas).
C----------------------------------------------------------------------------
#include "SOS.h"
#define INCTE_PI DACOS(-1.D+00)


C==============================================================================
C PROGRAM: SOS_ANGLES
C ==========
C      Program generating the angle files used by the codes SOS_AEROSOLS,
C      SOS_SURFACE and SOS.
C
C Input data of the program
C ------------------------------
C    -->Parameters/arguments of the program
C
C            NBMU_GAUSS_LUM (I4) : (I) Number of Gauss angles for the radiance computations.
C				           if it is CTE_NOT_DEFINED_VALUE_INT, then
C				           the default value CTE_DEFAULT_NBMU_LUM is used
C
C 	     TETAS   (DOUBLE) : (I) Solar zenith angle in degrees (0 < tetas < 90).
C      		               (complementary angle: mus= -cos(tetas))
C
C            FICANGLES_USER_LUM (CHARACTER) : (I) complete name of the "user" angles file
C                                                 for the radiance computations.
C                 			               (directory + file name + extension)
C
C            NBMU_GAUSS_MIE (I4) : (I) Number of Gauss angles for the phase function computations
C				           if it is CTE_NOT_DEFINED_VALUE_INT, then
C				           the default value CTE_DEFAULT_NBMU_MIE is used
C
C            FICANGLES_USER_MIE (CHARACTER) : (I) complete name of the "user" angles file
C                                                 for the phase function computations.
C                 			            (directory + file name + extension)
C
C            FICLOG (CHARACTER): (I) complete name of the trace file (directory + file name + extension)
C                                    or =0 if no trace
C
C            FICRES_LUM (CHARACTER) : (I) complete name of the file of the angles used
C                                         for the radiance and BRDF computations
C                 		          (directory + file name + extension)
C
C            FICRES_MIE (CHARACTER) : (I) complete name of the file of the angles used
C                                         for the phase function computations.
C                 		          (directory + file name + extension)
C
C      	     IER (I4) : (O) error code =0 if no error, =-1 otherwise
C
C Result provided by the procedure
C --------------------------------
C   The files describing the angles to be used for the phase function computations
C   and for the radiance + reflection matrix computations.
C   These files also contain the parameters giving the expansion orders
C   of the various functions: OS_NB (phase functions),
C                            OS_NS (max order of the expansions of the
C                                   Fresnel functions and of the radiances)
C   			        OS_NM (maximum order of the expansion in Fourier series
C                                   of the G function of the reflection matrix computations
C 
C Description of the files used
C ---------------------------------
C
C    User file FICANGLES_USER:
C           contains a list of angles (in degrees) to be added to the Gauss angles,
C           with null weights (=> no impact on the
C           spatial integrations)
C
C        Format: One angle value per record line.
C                 Angles must be between 0 and 90 degrees.
C                 Maximum number of records: CTE_NBMAX_USER_ANGLES
C
C    Trace file: identifier = 99
C
C    Result file FICRES_MIE:
C           contains the following information
C		  NB_TOTAL_ANGLES :  total nb of angles to be used
C               NB_GAUSS_ANGLES :  nb of Gauss angles
C               ANGLES_USERFILE :  Name of the user angles file
C               INTERNAL_OS_NB :   Order of the expansions of the phase functions
C				       in Legendre polynomials to be used.
C               INDEX COS_ANGLE WEIGHT : lines giving the angle number, its cosine and its weight
C                                        with format I4,1X,2D21.14
C
C    Result file FICRES_LUM:
C           contains the following information
C		  NB_TOTAL_ANGLES :    total nb of angles to be used
C               NB_GAUSS_ANGLES :    nb of Gauss angles
C               ANGLES_USERFILE :    Name of the user angles file
C	         SOLAR ZENITH ANGLE : Solar zenith angle (degrees)
C		  INTERNAL_IMUS :      Corresponding array index number
C               INTERNAL_OS_NB :     Order of the expansions of the phase functions
C				         in Legendre polynomials to be used.
C               INTERNAL_OS_NS :     Maximum order to be used for the expansion in
C                                    Legendre polynomials of the Fresnel matrix elements
C                                    and in Fourier series of the radiance computations.
C               INTERNAL_OS_NM :     Maximum order of the expansion in Fourier series
C                                    of the G function of the reflection matrix computations.
C               INDEX COS_ANGLE WEIGHT OUTPUT: lines giving the angle number, its cosine
C                                              and its weight, and whether it is a "user" output
C                                              with format I4,1X,2D21.14,1X,I4
C         
C
C Common used:
C --------------
C  None
C
C Error cases:
C ------------
C     - Incorrect decoding of parameters
C     - Error in a subroutine
C     - Error when opening, reading and writing the trace file
C       and the I/O files: user angles file and result files
C     - Invalid parameters
C     - Inconsistency of the dimension CTE_NBANGLES_MAX of the angle arrays of the routine
C       with respect to that of the arrays used in the SOS code:
C          MIE computations: CTE_MIE_NBMU_MAX and radiance: CTE_OS_NBMU_MAX.
C     - Inconsistency of a user angle with the expected range (0 to 90 degrees)
C     - Number of user angles too large compared with CTE_NBMAX_USER_ANGLES
C     - Number of useful angles (Gauss + user) too large compared with CTE_OS_NBMU_MAX
C
C     Message displayed on screen, program stop and
C     return of status 1 to the shell
C
C
C  Subroutines used:
C --------------------------
C  This program calls the routines:
C      - SOS_ANGLES_GAUSS_USER
C      - SOS_GAUSS  
C
C==============================================================================

      SUBROUTINE SOS_ANGLES(NBMU_GAUSS_LUM, TETAS,  
     &                      FICANGLES_USER_LUM,                 
     &			       NBMU_GAUSS_MIE, FICANGLES_USER_MIE, 
     &                      FICLOG,FICRES_LUM, FICRES_MIE, IER) 
      
      IMPLICIT NONE


C ANGLES constants
C -----------------
      DOUBLE PRECISION CONVDEGRAD       !Conversion coefficient from degrees to radians.

C* Definition of variables
C--------------------------------------------
      
      
      DOUBLE PRECISION TETAS	! Solar zenith angle in degrees (0 < tetas < 90).
      				       ! (complementary angle: mus= -cos(tetas))
      DOUBLE PRECISION XMUS     ! Cosine of the solar zenith angle
      DOUBLE PRECISION XMU(-CTE_NBANGLES_MAX:CTE_NBANGLES_MAX) !Cosines of the angles (Gauss and user)
      DOUBLE PRECISION XHR(-CTE_NBANGLES_MAX:CTE_NBANGLES_MAX) !Weights of the angles (Gauss and user)
      INTEGER*2 IND_ANGOUT(-CTE_NBANGLES_MAX:CTE_NBANGLES_MAX) !Indices of the angle arrays for which
      							       !there must be an output in the
							       !result file of the user angles
							       !(1: output required)
      DOUBLE PRECISION XMU_copy(-CTE_NBANGLES_MAX:CTE_NBANGLES_MAX) !copy of XMU before adding TETAS
      DOUBLE PRECISION XHR_copy(-CTE_NBANGLES_MAX:CTE_NBANGLES_MAX) !copy of XHR before adding TETAS
      INTEGER*2 IND_ANGOUT_copy(-CTE_NBANGLES_MAX:CTE_NBANGLES_MAX) !copy of IND_ANGOUT before adding TETAS
      
      INTEGER*2 NBARGS           ! Number of input arguments
     
      INTEGER*4 IARGC     	     ! Type of the function IARGC
      INTEGER*4 NPARAM 	     ! Number of the input parameter (type INTEGER*4 mandatory for f77)
      INTEGER*4 NBMU_GAUSS_LUM   ! Number of Gauss angles used for the radiance computations
      INTEGER*4 NBMU_LUM         ! Number of actual angles of the radiance simulations
      INTEGER*4 NBMU_GAUSS_MIE   ! Number of Gauss angles used for the Mie computations
      INTEGER*4 NBMU_MIE         ! Number of actual angles of the Mie computations
      INTEGER*4 IER     	     ! error code =0 if no error, =-1 otherwise
      INTEGER*4 J          	     ! Indices of the angle arrays
      INTEGER*4 IMUS             ! Index of the solar zenith angle in the angle array
      INTEGER*4 OS_NB 	     ! Order of the expansions of the phase functions
				     ! in Legendre polynomials.
      INTEGER*4 OS_NS    	     ! Maximum order of the expansion in Legendre polynomials
				     ! of the Fresnel matrix elements.
      INTEGER*4 OS_NM    	     ! Maximum order of the expansion in Fourier series
				     ! of the G function of the reflection matrix computations.
               
      CHARACTER*CTE_LENFIC2 FICLOG              ! name of the trace file
      			!Radiance computations:
      CHARACTER*CTE_LENFIC2 FICANGLES_USER_LUM  ! path and name of the "user" angles file
      CHARACTER*CTE_LENFIC2 FICRES_LUM          ! path and name of the file of the angles used
      			!Phase function computations:
      CHARACTER*CTE_LENFIC2 FICANGLES_USER_MIE  ! path and name of the "user" angles file
      CHARACTER*CTE_LENFIC2 FICRES_MIE          ! path and name of the file of the angles used
      


      LOGICAL TRACE     ! = true if writing in the trace file
      


	  	

	 
C* Initialisation
C*-----------------------------------------------------------------
      IER=0
      
      !Conversion coefficient from degrees to radians.
      CONVDEGRAD=INCTE_PI/180.D+00
      
C* Check of the parameters
C*-----------------------------------------------------------------

C     Optional arguments
C     --> Default value assigned if NBMU_GAUSS_MIE is not given
      IF (NBMU_GAUSS_MIE.EQ.CTE_NOT_DEFINED_VALUE_INT) THEN
          NBMU_GAUSS_MIE = CTE_DEFAULT_NBMU_MIE
	  OS_NB = CTE_DEFAULT_OS_NB
      ELSE
          !For a phase function discretised on NBMU_GAUSS_MIE (positive) Gauss angles,
	  !its decomposition in Legendre polynomials can only be correct if it is carried
	  !up to an order of at least 2 times the nb of Gauss angles used for the
	  !discretisation of the function.
          OS_NB = 2*NBMU_GAUSS_MIE	  
      ENDIF 

C     --> Default value assigned if NBMU_GAUSS_LUM is not given
      IF (NBMU_GAUSS_LUM.EQ.CTE_NOT_DEFINED_VALUE_INT) THEN
          NBMU_GAUSS_LUM = CTE_DEFAULT_NBMU_LUM
	  OS_NS = CTE_DEFAULT_OS_NS
	  OS_NM = CTE_DEFAULT_OS_NM
      ELSE
          !For a radiance field discretised on NBMU_GAUSS_LUM (positive) Gauss angles,
	  !the decomposition in Legendre polynomials of the Fresnel reflection functions
	  !and of the atmospheric radiances is only correct if it is carried up to
	  !an order of at least 2 times the nb of Gauss angles used to discretise
	  !the radiance.
          OS_NS = 2*NBMU_GAUSS_LUM
	  
	  !The computation of the G function of the BRDF computations (glitter, ...) requires a maximum order
	  !of the Fourier series expansions such that OS_NM >= OS_NB + OS_NS
	  OS_NM = OS_NB + OS_NS
      ENDIF 
      

C     --> Check of the limit orders of the expansions for the computation
C         of the coefficients of the phase and reflection matrices
      IF (OS_NB.LT.2*NBMU_GAUSS_MIE) GOTO 945
      IF (OS_NS.LT.2*NBMU_GAUSS_LUM) GOTO 946
      IF (OS_NM.LT.(OS_NB+OS_NS))    GOTO 947
       
      IF (OS_NB.GT.CTE_OS_NB_MAX) GOTO 9451
      IF (OS_NS.GT.CTE_OS_NS_MAX) GOTO 9461
      IF (OS_NM.GT.CTE_OS_NM_MAX) GOTO 9471
                  


C* Opening of the trace file
C------------------------------------------------------
      IF (FICLOG.EQ.'NO_LOG_FILE') THEN
	 TRACE=.FALSE.
      ELSE
         TRACE=.TRUE.
	 OPEN (99,FILE=FICLOG,ERR=990)
      ENDIF
  

C* Preparation of the angles for the phase function computations
C--------------------------------------------------------------      
      CALL SOS_ANGLES_GAUSS_USER("MIE",NBMU_GAUSS_MIE,
     &				 FICANGLES_USER_MIE,TRACE,
     &                       NBMU_MIE, XMU, XHR, IER) 
      
      IF (IER.EQ.-1) GOTO 955
      

C* Writing of the result angle file for the phase function
C* computations
C-------------------------------------------------------------------         
      OPEN(UNIT=10,FILE=FICRES_MIE,FORM='FORMATTED',ERR=960)
      WRITE(10,100,ERR=961) NBMU_MIE      
      WRITE(10,110,ERR=961) NBMU_GAUSS_MIE 
      WRITE(10,120,ERR=961) FICANGLES_USER_MIE
      WRITE(10,200,ERR=961) OS_NB
      WRITE(10,250,ERR=961) 
      DO J=1,NBMU_MIE
         WRITE(10,300,ERR=961) J, XMU(J), XHR(J)
      ENDDO
      
C* Preparation of the Gauss and user angles
C* for the radiance and BRDF/BPDF computations
C-------------------------------------------------------------------      
      CALL SOS_ANGLES_GAUSS_USER("LUM",NBMU_GAUSS_LUM,
     &                           FICANGLES_USER_LUM,TRACE,
     &                           NBMU_LUM, XMU, XHR, IER) 
      
      IF (IER.EQ.-1) GOTO 955      

C* Preparation of the array of indices to keep for an output of the
C* user angles
C-------------------------------------------------------------------         
      DO J=1,NBMU_LUM
      
         IF (XHR(J).EQ.0.D+00) THEN
	    IND_ANGOUT(J)=1
	 ELSE
	    IND_ANGOUT(J)=0
         ENDIF
	 
      ENDDO

C* Addition of the solar incidence angle
C-------------------------------------------------------------------         
      IMUS=-1
      XMUS=DCOS(TETAS*CONVDEGRAD)
      
      !Search for a match with an already defined angle
      !-->Only the value of IMUS is kept, the arrays are not modified
      J=1
      DO J=1,NBMU_LUM 
         IF (DABS(XMUS-XMU(J)).LT.CTE_SEUIL_ECART_MUS) IMUS=J
      ENDDO
      
      !If TETAS is not an already defined angle
      !-->Added to the angle array, without output required for the user angles
      IF (IMUS.EQ.-1) THEN 
         !The index IMUS is estimated in the new array of NBMU_LUM+1 elements.
	 
	 !-->Storage of the arrays without TETAS
         DO J=1,NBMU_LUM
	    XMU_copy(J) = XMU(J)
	    XHR_copy(J) = XHR(J)
	    IND_ANGOUT_copy(J) = IND_ANGOUT(J)
	 ENDDO
	 
	 !The array XMU_copy is ordered from the largest to the smallest element
	 !-->Check whether XMUS must be the first element of the new array
	 IF (XMUS.GT.XMU_copy(1)) THEN 
	     IMUS=1
	 ELSE
	     !-->Check whether XMUS must be the last element of the new array
	     IF (XMUS.LT.XMU_copy(NBMU_LUM)) THEN 
	        IMUS=NBMU_LUM+1
	     ELSE
	        !-->Search for the element J of the initial array to be replaced by XMUS
	        DO J=1,NBMU_LUM-1
		   IF ( (XMUS.GT.XMU_copy(J+1)).
     &                   AND.(XMUS.LT.XMU_copy(J)) ) IMUS=J+1
                ENDDO 
	     ENDIF
	 ENDIF
	 
	 !-->Update of the array size
	 NBMU_LUM=NBMU_LUM+1
	 IF (NBMU_LUM.GT.CTE_OS_NBMU_MAX) GOTO 980
	 
	 !-->Positioning of TETAS in the arrays
	 XMU(IMUS) = XMUS
	 XHR(IMUS) = 0.D+00
	 IND_ANGOUT(IMUS)= 0
	 
	 !-->Completion of the arrays with index shift
	 DO J=IMUS+1,NBMU_LUM  
	    XMU(J) = XMU_copy(J-1)
	    XHR(J) = XHR_copy(J-1)
	    IND_ANGOUT(J) = IND_ANGOUT_copy(J-1)
	 ENDDO
	    
        !-->Completion of the angles defined on J=-NBMU_UTILES,0
         DO J=1,NBMU_LUM
	        XMU(-J) = -XMU(J)
	        XHR(-J) =  XHR(J)
		IND_ANGOUT(-J) = IND_ANGOUT(J)
         ENDDO 
         XMU(0) = 0.D+00
         XHR(0) = 0.D+00
	 IND_ANGOUT(0) = 0
	 	 
      ENDIF !Test whether IMUS is undefined
 
      
C* Trace information
C-------------------------------------------------------------------         
      IF (TRACE) THEN
        WRITE (99,*,ERR=990) "Addition of the solar zenith angle"
	 WRITE (99,*,ERR=990) "---------------------------------"
	 WRITE (99,*,ERR=990) " "
	 WRITE(99,*,ERR=990) "TETAS : ", TETAS
	 WRITE(99,*,ERR=990) "IMUS  : ", IMUS
	 WRITE(99,*,ERR=990) "Final number of angles : ", NBMU_LUM 
        WRITE (99,*,ERR=990) " " 
	 WRITE (99,*,err=990) "List of used angles :"
	 WRITE (99,*,err=990) 
     &	 "   Index, mu, angle (deg.), weight, user angles flag"
	 DO J=1,NBMU_LUM
	     WRITE(99,223,ERR=990) J, XMU(J),DACOS(XMU(J))/CONVDEGRAD,
     &	                           XHR(J),IND_ANGOUT(J)  
	 ENDDO
	 WRITE (99,*,ERR=990) " " 
	 WRITE (99,*,ERR=990) " "  	
      ENDIF !Trace file
		

C* Writing of the result angle file for the phase function
C* computations
C-------------------------------------------------------------------         
      OPEN(UNIT=20,FILE=FICRES_LUM,FORM='FORMATTED',ERR=960)
      WRITE(20,100,ERR=961) NBMU_LUM
      WRITE(20,110,ERR=961) NBMU_GAUSS_LUM
      WRITE(20,120,ERR=961) FICANGLES_USER_LUM
      WRITE(20,410,ERR=961) TETAS
      WRITE(20,420,ERR=961) IMUS
      WRITE(20,200,ERR=961) OS_NB
      WRITE(20,500,ERR=961) OS_NS
      WRITE(20,510,ERR=961) OS_NM
      WRITE(20,550,ERR=961) 
      DO J=1,NBMU_LUM
         WRITE(20,600,ERR=961) J, XMU(J), XHR(J), IND_ANGOUT(J)
      ENDDO
     
     
C* File closing
C-------------------
      CLOSE(10)
      CLOSE(20)
      		
C* Nominal end
C-------------------
      GOTO 9999     

C* Error case and return of status 1 to the shell
C----------------------------------------------
  900 WRITE(6,*) '  ERROR on argument loading for SOS_ANGLES'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -ANG.Rad.NbGauss value : Number of Gauss'
      WRITE(6,*) '       angles for radiance computations (I4).' 
       IER=-1
      GOTO 9999     

  905 WRITE(6,*) '  ERROR on argument loading for SOS_ANGLES'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -ANG.Aer.NbGauss value : Number of Gauss'
      WRITE(6,*) '       angles for radiance computations (I4).' 
       IER=-1
      GOTO 9999     
      
  910 WRITE(6,*) '  ERROR on argument loading for SOS_ANGLES'
      WRITE(6,*) '  --> Check required formats for -ANG.Thetas value :'
      WRITE(6,*) '      Solar zenithal angle (F9, degrees).'
       IER=-1
      GOTO 9999  

  920 WRITE(6,*) '  ERROR on parameters for SOS_ANGLES'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -ANG.Rad.Thetas Value'
       IER=-1
      GOTO 9999 

  925 WRITE(6,*) '  ERROR on -ANG.Rad.Thetas value : ',
     &           'out of available range'
       IER=-1
      GOTO 9999 
      
      
  930 WRITE(6,*) '  ERROR on parameters for SOS_ANGLES'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -ANG.Rad.ResFile Value'
       IER=-1
      GOTO 9999 

  935 WRITE(6,*) '  ERROR on parameters for SOS_ANGLES'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -ANG.Aer.ResFile Value'
       IER=-1
      GOTO 9999        
                 
  940 WRITE(6,*) '  ERROR on parameters for SOS_ANGLES'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -ANG.Log Value'
       IER=-1
      GOTO 9999 

  945 WRITE(6,*) '  ERROR : CTE_DEFAULT_OS_NB <= 2*CTE_DEFAULT_NBMU_MIE'
      WRITE(6,*) '  --> check value in SOS.h and compile again'
       IER=-1
      GOTO 9999 

  946 WRITE(6,*) '  ERROR : CTE_DEFAULT_OS_NS <= 2*CTE_DEFAULT_NBMU_LUM'
      WRITE(6,*) '  --> check value in SOS.h and compile again'
       IER=-1
      GOTO 9999   

  947 WRITE(6,*) '  ERROR : ',
     &    'CTE_DEFAULT_OS_NM < CTE_DEFAULT_OS_NB + CTE_DEFAULT_OS_NS'
      WRITE(6,*) '  --> check value in SOS.h and compile again'
       IER=-1
      GOTO 9999                

 9451 WRITE(6,*) '  ERROR : OS_NB > CTE_OS_NB_MAX'
      WRITE(6,*) '  --> check value in SOS.h and compile again'
       IER=-1
      GOTO 9999 

 9461 WRITE(6,*) '  ERROR : OS_NS > CTE_OS_NS_MAX'
      WRITE(6,*) '  --> check value in SOS.h and compile again'
       IER=-1
      GOTO 9999   

 9471 WRITE(6,*) '  ERROR : OS_NM > CTE_OS_NM_MAX'
      WRITE(6,*) '  --> check value in SOS.h and compile again'
       IER=-1
      GOTO 9999                
            
  950 WRITE(6,*) '  ERROR on logfile opening for SOS_ANGLES'
       IER=-1
      GOTO 9999  

  955 WRITE(6,*) '  ERROR on subroutine SOS_ANGLES_GAUSS_USER'
       IER=-1
      GOTO 9999

  960 WRITE(6,*) '  ERROR on ANGLES file opening for SOS_ANGLES'
       IER=-1
      GOTO 9999  

  961 WRITE(6,*) '  ERROR on ANGLES file writing for SOS_ANGLES'
       IER=-1
      GOTO 9999  
                   
  980 WRITE(6,*) "Error on the number of used radiance angles "
      WRITE(6,*) "(gauss angles + possibly solar angle)  : "
      WRITE(6,*) "  --> greater than the maximal number allowed "
      WRITE(6,*) "       CTE_OS_NBMU_MAX =",CTE_OS_NBMU_MAX
       IER=-1
      GOTO 9999 
      
  990 WRITE(6,*) '  ERROR on logfile writing for SOS_ANGLES'
       IER=-1
      GOTO 9999
                                              
 9999 RETURN
 
C* Format
C---------   

      
  100 FORMAT(17hNB_TOTAL_ANGLES :,I4)
  110 FORMAT(17hNB_GAUSS_ANGLES :,I4)
  120 FORMAT(17hANGLES_USERFILE :,A)
  200 FORMAT(16hINTERNAL_OS_NB :,I4)
  250 FORMAT(35hINDEX   COS_ANGLE            WEIGHT)
  300 FORMAT(I4,1X,2D21.14)
  
  410 FORMAT(20hSOLAR ZENITH ANGLE :,F7.3)
  420 FORMAT(15hINTERNAL_IMUS :,I4)
  500 FORMAT(16hINTERNAL_OS_NS :,I4)
  510 FORMAT(16hINTERNAL_OS_NM :,I4)
  550 FORMAT
     &(58hINDEX   COS_ANGLE            WEIGHT             USER_ANGLE)
  600 FORMAT(I4,1X,2D21.14,1X,I4)
  223 FORMAT(I4,1X,D21.14,2X,F7.3,2X,D21.14,2X,I2)
  
      END	!END OF PROGRAM SOS_ANGLES


      

C==============================================================================
C PROCEDURE: SOS_ANGLES_GAUSS_USER
C ==========
C
C  The routine computes the Gauss angles and weights and adds, if required,
C  the angles of the user file.
C
C  It sorts the angles
C
C Description of parameters
C ----------------------------		  
C
C       USAGE (CHAR*3)  (I) : indicator of the use of the angle file:
C                             MIE computations or simulation of the radiation field LUM
C	 NBMU_GAUSS (I4) (I) : number of requested Gauss angles
C       FICANGLES_USER (CHAR*CTE_LENFIC2) (I) : name of the "user" angles file
C                                               (complete path)
C       TRACE (LOGICAL) (I)  : true if writing in the trace file
C       NBMU_UTILES (I4) (O) : number of angles after concatenation of the Gauss angles
C                              and of the "user" angles
C	 XMU (double)     (O) : array of the cosines of the computed angles.
C	 XHR (double)     (O) : array of the weights associated with the computed angles.
C       IER (I4)         (O) : error index: 0 if no error, -1 otherwise.
C     
C Result provided by the procedure
C --------------------------------
C   The arrays XMU and XHR of the cosines and weights of the useful angles
C   The nb of angles NBMU_UTILES
C   An error return IER
C   
C Description of the files used
C ---------------------------------
C
C    User file FICANGLES_USER:
C           contains a list of angles (in degrees) to be added to the Gauss angles,
C           with null weights (=> no impact on the
C           spatial integrations)
C
C        Format: One angle value per record line.
C                 Angles must be between 0 and 90 degrees.
C                 Maximum number of records: CTE_NBMAX_USER_ANGLES
C
C    Trace file: identifier = 99 (consistent with the opening by the calling code)
C
C Common used:
C --------------
C  None
C
C Error cases:
C ------------
C  None
C
C  Subroutines used:
C --------------------------
C  SOS_GAUSS
C
C==============================================================================
      
      SUBROUTINE  SOS_ANGLES_GAUSS_USER(USAGE,NBMU_GAUSS,FICANGLES_USER,
     &                                  TRACE,NBMU_UTILES,XMU,XHR,IER) 

      IMPLICIT NONE


C ANGLES constants
C -----------------
      DOUBLE PRECISION CONVDEGRAD       !Conversion coefficient from degrees to radians.

C* Definition of variables
C--------------------------------------------
      INTEGER*4 NBMU_GAUSS      ! Number of Gauss angles used
      INTEGER*4 NBMU_USER       ! Number of angles in the user file
      INTEGER*4 NBMU_UTILES     ! Number of actual angles of the simulations
      INTEGER*4 IER             ! error code =0 if no error, =-1 otherwise
      INTEGER*4 J, K            ! Indices of the angle arrays
			    
      DOUBLE PRECISION XMU(-CTE_NBANGLES_MAX:CTE_NBANGLES_MAX) !Cosines of the angles (Gauss and user)
      DOUBLE PRECISION XHR(-CTE_NBANGLES_MAX:CTE_NBANGLES_MAX) !Weights of the angles (Gauss and user)
      DOUBLE PRECISION VAL	!Angle value read in the user file
      DOUBLE PRECISION MU_COMP(CTE_NBMAX_USER_ANGLES)  
                                        !Cosines of the angles (Gauss and user)
 
      CHARACTER*CTE_LENFIC2 FICANGLES_USER  ! path and name of the "user" angles file
      CHARACTER*3 USAGE     ! Indicator of the use of the angle file:
      			    ! MIE computations or simulation of the radiation field LUM


      LOGICAL TRACE     ! = true if writing in the trace file

C* Initialisation
C-----------------
      IER=0
      NBMU_USER=0 
      
      !Conversion coefficient from degrees to radians.
      CONVDEGRAD=INCTE_PI/180.D+00
 
C* Validity check of the parameter USAGE
C------------------------------------------------------      
      IF ((USAGE.NE.'MIE').AND.(USAGE.NE.'LUM')) GOTO 950
      
C* Consistency check of the dimensions
C------------------------------------------------------      
      IF (USAGE.EQ.'MIE') THEN
          IF (CTE_NBANGLES_MAX.LT.CTE_MIE_NBMU_MAX) GOTO 960
      ENDIF
      
      IF (USAGE.EQ.'LUM') THEN
          IF (CTE_NBANGLES_MAX.LT.CTE_OS_NBMU_MAX) GOTO 961
      ENDIF
      
      IF (NBMU_GAUSS.GE.CTE_NBANGLES_MAX) GOTO 962
           
      
C* Computation of the Gauss angles and weights
C  -----------------------------------
      CALL SOS_GAUSS(NBMU_GAUSS+1, XMU, XHR)
	
C* Reading of the user angles
C------------------------------------------------------ 	
      IF (FICANGLES_USER.NE.'NO_USER_ANGLES') THEN
      	 
         OPEN(1,FILE=FICANGLES_USER,STATUS='OLD',ERR=970)      
      
         J=0  
   20	  READ(1,*,ERR=971,END=30) VAL
              IF ((VAL.LT.0.).OR.(VAL.GT.90.)) GOTO 972              
	       J=J+1
	       IF (J.GT.CTE_NBMAX_USER_ANGLES) GOTO 973
	       MU_COMP(J)=DCOS(VAL*CONVDEGRAD)
	  GOTO 20
   30	  CLOSE(1)
         NBMU_USER=J
    
      ENDIF

C* Determination of the number of useful angles
C------------------------------------------------------      
      NBMU_UTILES = NBMU_GAUSS + NBMU_USER

      IF (USAGE.EQ.'MIE') THEN
          IF (NBMU_UTILES.GT.CTE_MIE_NBMU_MAX) GOTO 980
      ENDIF 
      
      IF (USAGE.EQ.'LUM') THEN
          IF (NBMU_UTILES.GT.CTE_OS_NBMU_MAX) GOTO 981
      ENDIF 

C*  Addition of the user angles to the Gauss angles
C------------------------------------------------------           
      IF (NBMU_USER.NE.0) THEN
      
          DO J=1,NBMU_USER
	     XMU(NBMU_GAUSS+J) = MU_COMP(J)
	     XHR(NBMU_GAUSS+J) = 0.D+00
	  ENDDO 

      ENDIF
      
C*  Sorting of the angles: case USAGE="MIE" --> sorted from MUmin to MUmax
C*                         case USAGE="LUM" --> sorted from MUmax to MUmin
C*    Sorting of the angles defined on J=1,NBMU_UTILES
C------------------------------------------------------           
      IF (USAGE.EQ.'MIE') THEN
 
          DO J=1,NBMU_UTILES	     
	     DO K=J+1,NBMU_UTILES
	        IF (XMU(J).GT.XMU(K)) THEN		    		    
		    !Swap of the values XMU(J) and XMU(K)
		    VAL=XMU(J)
		    XMU(J)=XMU(K)
		    XMU(K)=VAL
		    !Swap of the values XMU(J) and XHR(K)
		    VAL=XHR(J)
		    XHR(J)=XHR(K)
		    XHR(K)=VAL   	    	    
		ENDIF
	     ENDDO 
	  ENDDO 
	  
      ENDIF !USAGE=MIE	  

      IF (USAGE.EQ.'LUM') THEN
 
          DO J=1,NBMU_UTILES	     
	     DO K=J+1,NBMU_UTILES
	        IF (XMU(J).LT.XMU(K)) THEN		    		    
		    !Swap of the values XMU(J) and XMU(K)
		    VAL=XMU(J)
		    XMU(J)=XMU(K)
		    XMU(K)=VAL
		    !Swap of the values XMU(J) and XHR(K)
		    VAL=XHR(J)
		    XHR(J)=XHR(K)
		    XHR(K)=VAL       
		ENDIF
	     ENDDO 
	  ENDDO 
	  
      ENDIF !USAGE=LUM	  
      	  	
	  	  
C*  Completion of the angles defined on J=-NBMU_UTILES,0
C------------------------------------------------------      	  
      DO J=1,NBMU_UTILES
	     XMU(-J) = -XMU(J)
	     XHR(-J) = XHR(J)
      ENDDO 
      XMU(0) = 0.D+00
      XHR(0) = 0.D+00
	
C*  Information to the TRACE file
C------------------------------------------------------   
      IF (TRACE) THEN
         IF (USAGE.EQ.'MIE') THEN
	     WRITE (99,*,ERR=990) "CALCULS DE FONCTIONS DE PHASE"
	     WRITE (99,*,ERR=990) "*****************************"
	 ENDIF
	 IF (USAGE.EQ.'LUM') THEN
	     WRITE (99,*,ERR=990) "CALCULS DE LUMINANCE"
	     WRITE (99,*,ERR=990) "********************"
	 ENDIF
	 WRITE (99,*,ERR=990) " "
	 WRITE(99,*,ERR=990) "Nombre d'angles de Gauss : ", NBMU_GAUSS
	 WRITE (99,*,ERR=990) " "
         IF (FICANGLES_USER.NE.'NO_USER_ANGLES') THEN
	     WRITE(99,*,ERR=990)  
     &	       "Usage du fichier d'angles complementaires :",
     &         FICANGLES_USER
     	     WRITE(99,*,ERR=990) "Nb d'angles complementaires :",
     &         (NBMU_UTILES-NBMU_GAUSS)
	 ELSE
	     WRITE(99,*,ERR=990)  
     &	       "Pas de fichier d'angles complementaires"
	 ENDIF
        WRITE (99,*,ERR=990) " " 
	 WRITE (99,*,err=990) "Liste des angles utilises :"
	 WRITE (99,*,err=990) "   Indice, mu, angle (degrés), poids"
	 DO J=1,NBMU_UTILES
	     WRITE(99,222,ERR=990) J, XMU(J),
     &	                           DACOS(XMU(J))/CONVDEGRAD, XHR(J)  
	 ENDDO
	 WRITE (99,*,ERR=990) " " 
	 WRITE (99,*,ERR=990) " "  	
      ENDIF !Trace file
C* Nominal end
C-------------------
      GOTO 9999


C* Error case and return of status -1 to the calling program
C----------------------------------------------------------            
  950  WRITE(6,*) "Error on USAGE parameters : has to be 'MIE' or 'LUM'"
       IER=-1
       GOTO 9999        
       
  960  WRITE(6,*) "Inconsistence on size of tables in SOS.h."
       WRITE(6,*) " --> required : CTE_NBANGLES_MAX >= CTE_MIE_NBMU_MAX"
       IER=-1
       GOTO 9999   

  961  WRITE(6,*) "Inconsistence on size of tables in SOS.h."
       WRITE(6,*) " --> required : CTE_NBANGLES_MAX >= CTE_OS_NBMU_MAX"
       IER=-1
       GOTO 9999        

  962  WRITE(6,*) "Inconsistence on size of tables in SOS.h."
       WRITE(6,*) " --> Number of required Gauss angles is ", NBMU_GAUSS
       WRITE(6,*) "     greater than the limit CTE_NBANGLES_MAX :", 
     &            CTE_NBANGLES_MAX
       IER=-1
       GOTO 9999        
          
  970  WRITE(6,*) "Error while opening the user's angles file : "
       WRITE(6,*)  FICANGLES_USER
       WRITE(6,*) " --> check its localization"
       IER=-1
       GOTO 9999      

  971  WRITE(6,*) "Error while reading the user's angles file :"
       WRITE(6,*)  FICANGLES_USER
       IER=-1
       GOTO 9999   

  972  WRITE(6,*) "Error on an user's angles value in file: "
       WRITE(6,*)  FICANGLES_USER
       WRITE(6,*) "  --> have to defined between 0 and 90°"
       IER=-1
       GOTO 9999   

  973  WRITE(6,*) "Error on the user's angles file : "
       WRITE(6,*)  FICANGLES_USER
       WRITE(6,*) "  --> The number of angles has to be < or = to"
       WRITE(6,*) "       CTE_NBMAX_USER_ANGLES =",CTE_NBMAX_USER_ANGLES
       IER=-1
       GOTO 9999

  980  WRITE(6,*) "Error on the number of used phaze function angles : "
       WRITE(6,*) "  --> greater than the maximal number allowed "
       WRITE(6,*) "       CTE_MIE_NBMU_MAX =",CTE_MIE_NBMU_MAX
       IER=-1
       GOTO 9999 
       
  981  WRITE(6,*) "Error on the number of used radiance angles : "
       WRITE(6,*) "  --> greater than the maximal number allowed "
       WRITE(6,*) "       CTE_OS_NBMU_MAX =",CTE_OS_NBMU_MAX
       IER=-1
       GOTO 9999 
      
  990 WRITE(6,*) '  ERROR on logfile writing for SOS_ANGLES'
      CALL EXIT(1)
      GOTO 9999
              
9999  RETURN


C* Format
C---------   
  222 FORMAT(I4,1X,D21.14,2X,F7.3,2X,D21.14,2X)
   
      END !END OF PROCEDURE SOS_ANGLES_GAUSS_USER
             
	     
	     
	     
C==============================================================================
C PROCEDURE: SOS_GAUSS
C ==========
C
C  Original routine provided by the LOA (M. Herman, October 2009).
C  No documentation.
C
C  The routine computes the Gauss angles and weights for MM-1 values.
C
C  The generated arrays cover the values XMU(-MM:MM) and XHR(-MM:MM)
C  with XMU(-MM)=-1 and XMU(MM)=1 having null weights.
C
C  A constraint is to have CTE_NBANGLES_MAX >= MM
C
C Description of parameters
C ----------------------------
C	MM (I4) (I) : number of requested Gauss angles
C	XMU (double) (O) : array of the cosines of the computed angles.
C	XHR (double) (O) : array of the weights associated with the computed angles.
C
C Result provided by the procedure
C --------------------------------
C   The arrays XMU and XHR
C   
C Description of the files used
C ---------------------------------
C  None
C
C Common used:
C --------------
C  None
C
C Error cases:
C ------------
C  None
C
C  Subroutines used:
C --------------------------
C  None
C
C==============================================================================
      
      SUBROUTINE  SOS_GAUSS(MM,AMU,PMU)
C     N : ORDER OF THE QUADRATURE ON XL TO XU

C* Definition of variables
C*-----------------------------------------------------------------
      PARAMETER (IX=2*CTE_NBANGLES_MAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Z(IX),PA(IX),W(IX),R(IX)
      DOUBLE PRECISION AMU(-CTE_NBANGLES_MAX:CTE_NBANGLES_MAX)  !Gauss angles
      DOUBLE PRECISION PMU(-CTE_NBANGLES_MAX:CTE_NBANGLES_MAX)  !Gauss weights

 
      TOL = 1.0D-15
      PI = INCTE_PI
      N=2*MM-2
      XL=-1
      XU=+1
      AA = 2.0D+00/PI**2
      AB = -62.0D+00/(3.0D+00*PI**4)
      AC = 15116.0D+00/(15.0D+00*PI**6)
      AD = -12554474.D+00/(105.0D+00*PI**8)
      PA(1) = 1.D0
      EN = N
      NP1 = N+1
      U= 1.0D+00-(2.0D+00/PI)**2
      D = 1.0D+00/DSQRT((EN+0.5D+00)**2+U/4.0D+00)

      DO 100 I= 1,N
         SM = I
         AZ = 4.0D+00*SM-1.0D+00
         AE = AA/AZ
         AF = AB/AZ**3
         AG = AC/AZ**5
         AH = AD/AZ**7
         Z(I) = 0.25D+00*PI*(AZ+AE+AF+AG+AH)
	 
100   CONTINUE
	
      DO 200 K = 1,N
         X = DCOS(Z(K)*D)
    1    PA(2) = X
    	
         DO 201  NN = 3,NP1
         ENN = NN-1
         PA(NN) =
     &   ((2.0D+00*ENN-1.0D+00)*X*PA(NN-1)-(ENN-1.0D+00)*PA(NN-2))/ENN
  201    CONTINUE
  
         PNP = EN*(PA(N)-X*PA(NP1))/(1.0D+00-X*X)
         XI = X-PA(NP1)/PNP
         XD = DABS(XI-X)
         XDD = XD-TOL
	 
         IF (XDD.LE.0.D+00) THEN	 
	     R(K) = X
             W(K) = 2.0D+00*(1.0D+00-X*X)/(EN*PA(N))**2
	  ELSE
	     X = XI
             GO TO 1
	  ENDIF	    
	    
  200 CONTINUE
	  
      AP = (XU-XL)/2.D0
      BP = (XU+XL)/2.D0
      DO I=1,MM-1
         K=MM-I
         AMU(K)=BP+AP*R(I)
         PMU(K)=AP*W(I)
         AMU(-K)=-AMU(K)
         PMU(-K)=PMU(K)
      ENDDO
    	
      AMU(-MM)=-1.
      AMU(MM)=1.
      AMU(0)=0.
      PMU(0)=0.
      PMU(-MM)=0.
      PMU(MM)=0.
C
      RETURN
      END !END OF PROCEDURE SOS_GAUSS
      
