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
C* FILE: SOS_NOM_FIC_SURFACE.F
C* PROJECT: Successive orders of scattering
C* ROLE: This module manages the naming of the BRDF / BPDF files for the OS.
C*
C* AUTHOR:
C*    CS GROUP France, B.Lafrance.
C*
C* DATE: 30/04/2002
C* MOD:VERSION:1.0: Creation
C*
C* MOD:VERSION:1.1: Correction of the field for the limit solar angle of the
C*                  computations of the Roujean model (VS becomes SL: Solar Limit angle).
C* MOD:VERSION:1.2: Update of the comments.
C*
C* MOD:VERSION:2.0: 10/01/2010 
C*     - Creation of the procedure SOS_COMPLETE_NOM_FIC_SURF
C*     - Addition and management of the input arguments NBMU_GAUSS and FICANGLES_USER
C*       and call of the procedure SOS_COMPLETE_NOM_FIC_SURF by the procedures
C*          SOS_NOM_FICGLITTER, SOS_NOM_FICROUJ, SOS_NOM_FICRH, SOS_NOM_FICBREON, 
C*          SOS_NOM_FICNADAL
C*
C* MOD:VERSION:2.1: 24/02/2016 
C*    - Minor adjustments to strictly respect the 72 columns
C*
C* MOD:VERSION:2.2: 17/07/2020 
C*    - Identifier "ROUJ" changed to "ROUJEAN"
C*      in the name of the produced files
C*    - Addition of the routine SOS_NOM_FICMAIGNAN
C*
C* MOD:VERSION:2.3: 07/07/2022
C*    - Constants renamed from SOS_* to CTE_*
C*
C*
C******************************************************************************
#include "SOS.h"

C==============================================================================
C PROCEDURE: SOS_NOM_FICGLITTER
C ==========
C      This procedure generates the name of the Glitter file from:
C	- the refractive index water / air: IND
C      - the wind speed: WIND
C      - the number of Gauss angles used NBMU_GAUSS
C      - the name of the user angle file FICANGLES_USER
C          (without its directory tree)
C
C      The following formats are imposed:
C         IND with format F5.3.
C         WIND with format F4.1,
C
C      The name of the Glitter file must be:
C
C      GLITTER-IND-WIND-MUnb-NomFicUser with
C      - IND with format F5.3 (e.g. 1.500)
C      - WIND with format F4.1 (e.g. 12.5)
C      - MUnb the number of Gauss angles (e.g. MU24)
C      - NomFicUser the name of the user file (without directory tree) (e.g. MyAngles)
C
C      WIND is converted into a character string CWIND(1:4)
C      For this, Format F4.1 * 10 => A4
C      The last character CRN(4) is the decimal of WIND,
C           etc.
C
C Description of input / output parameters
C ------------------------------------------
C	NBMU_GAUSS (I4)	: (I) Number of Gauss angles used
C       FICANGLES_USER (CHARACTER*CTE_LENFIC2)	: (I) Name of the "user"
C                                                 angles file (complete path)
C	WIND (double)	: Wind speed m/s
C	IND (double)	: Water index
C	FICS (CHARACTER*CTE_LENFIC1) : (O) name of the GLITTER file (without directory tree)
C
C 
C Constants used:
C --------------------
C     CTE_LENFIC1 : Length of file names without directory tree.
C     CTE_LENFIC2 : Length of file names with directory tree.
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C None
C
C  Subroutines used:
C --------------------------
C  This program calls the routine:
C     - SOS_COMPLETE_NOM_FIC_SURF
C
C==============================================================================
      SUBROUTINE SOS_NOM_FICGLITTER(NBMU_GAUSS,FICANGLES_USER,
     &                              IND,WIND,FICS) 


      IMPLICIT NONE

C* Declaration of variables
C  -------------------------         
      DOUBLE PRECISION WIND 	  ! (I) Wind speed m/s
      DOUBLE PRECISION IND	  ! (I) Water index.
				
      CHARACTER*CTE_LENFIC2 FICANGLES_USER  ! (I) path and name of the angle file
                                            ! of the "user" angles
      
      CHARACTER*CTE_LENFIC1 FICS  ! (O) name of the GLITTER file
          
      CHARACTER*3 CWIND	! Character strings
      CHARACTER*4 CIND
      
      INTEGER*2 I   
      INTEGER*4 NBMU_GAUSS  ! (I) Number of Gauss angles used
      			      	 
	    
C --------- Conversion of the input parameters into character strings
      WRITE(CIND,'(I4)') int(IND*1000)		
      WRITE(CWIND,'(I3)') int(WIND*10)  		

C --------- Blanks are replaced by 0
      DO I=1,4
         IF( CIND(I:I).EQ.' ') CIND(I:I) = '0'
      ENDDO   

      DO I=1,3
	 IF( CWIND(I:I).EQ.' ') CWIND(I:I) = '0'
      ENDDO
      
     
      FICS = 'GLITTER-'
           
C -------- Building of the character string for the index IND
      FICS(9:9) = CIND(1:1)     ! Units digit
      FICS(10:10) = '.'	    ! Decimal point
      FICS(11:13) = CIND(2:4)   ! First three decimals
      FICS(14:14) = '-'

C -------- Building of the character string for WIND
      FICS(15:16) = CWIND(1:2)     ! Units digits
      FICS(17:17) = '.'	       ! Decimal point
      FICS(18:18) = CWIND(3:3)     ! first decimal
      FICS(19:19) = '-'

C -------- Completion of the end of the file name with MUnb-NomFicUser
      CALL SOS_COMPLETE_NOM_FIC_SURF(NBMU_GAUSS,FICANGLES_USER,FICS)
      	  
     	      
      END !SOS_NOM_FICGLITTER
 



C==============================================================================
C PROCEDURE: SOS_NOM_FICROUJ
C ==========
C      This procedure generates the name of the file of the ROUJEAN model
C      from the parameters K0,K1 and K2
C
C      The following formats are imposed:
C	  kx with format F7.3.
C
C      The name of the Roujean file must be:
C
C      ROUJEAN-K0-K1-K2-VL{CTE_TETAV_LIM_ROUJEAN}-
C          -SL{CTE_TETAS_LIM_ROUJEAN}-MUnb-NomFicUser with
C      - kx with format F5.3 (e.g. 000.500)
C      - CTE_TETAV_LIM_ROUJEAN and CTE_TETAS_LIM_ROUJEAN respectively the
C        limit viewing and solar zenith angles for the BRDF computations
C        (in degrees).
C      - MUnb the number of Gauss angles (e.g. MU24)
C      - NomFicUser the name of the user file (without directory tree) (e.g. MyAngles)
C
C      Kx is converted into a character string CKx(1:6)
C      For this, Format F7.3 * 1000 => A6
C        The first character CKx(1) is the hundreds digit of Kx,
C        The second character CKx(2) is the tens digit of Kx,
C        The third character CKx(3) is the units digit of Kx,
C           etc.
C
C
C Description of input / output parameters
C ------------------------------------------
C	NBMU_GAUSS (I4)	: (I) Number of Gauss angles used
C       FICANGLES_USER (CHARACTER*CTE_LENFIC2)	: (I) Name of the "user"
C                                                 angles file (complete path)
C	K0  (double)	: Parameter of the Roujean model
C	K1  (double)	: Parameter of the Roujean model
C	K2  (double)	: Parameter of the Roujean model
C	FICS (CHARACTER*CTE_LENFIC1) : (O) name of the Roujean file (without directory tree)
C
C 
C Constants used:
C --------------------
C     CTE_TETAV_LIM_ROUJEAN : Maximum viewing zenith angle for the BRDF
C                             computations.
C     CTE_TETAS_LIM_ROUJEAN : Maximum solar zenith angle for the BRDF
C                             computations.
C     CTE_LENFIC1 : Length of file names without directory tree.
C     CTE_LENFIC2 : Length of file names with directory tree.
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C None
C
C  Subroutines used:
C --------------------------
C  This program calls the routine:
C     - SOS_COMPLETE_NOM_FIC_SURF
C
C==============================================================================
      SUBROUTINE SOS_NOM_FICROUJ(NBMU_GAUSS,FICANGLES_USER,K0,K1,K2,
     &                           FICS) 


      IMPLICIT NONE

C* Declaration of variables
C  -------------------------              
      DOUBLE PRECISION K0,K1,K2    ! (I) Parameters of the Roujean model

      CHARACTER*CTE_LENFIC2 FICANGLES_USER  ! (I) path and name of the angle file
                                            ! of the "user" angles
      CHARACTER*CTE_LENFIC1 FICS   ! (O) Name of the ROUJEAN file
          
      CHARACTER*6 CK0,CK1,CK2	   ! Character strings
      CHARACTER*4 CVLIM,CSLIM	   ! Character strings
      
      INTEGER*2 I
      INTEGER*4 NBMU_GAUSS  ! (I) Number of Gauss angles used

	    
C --------- Conversion of the input parameters into character strings
      WRITE(CK0,'(I6)') int(K0*1000)		
      WRITE(CK1,'(I6)') int(K1*1000)
      WRITE(CK2,'(I6)') int(K2*1000)
      WRITE(CVLIM,'(I2)') int(CTE_TETAV_LIM_ROUJEAN)
      WRITE(CSLIM,'(I2)') int(CTE_TETAS_LIM_ROUJEAN)
      
C --------- Blanks are replaced by 0
      DO I=1,6
         IF( CK0(I:I).EQ.' ') CK0(I:I) = '0'
         IF( CK1(I:I).EQ.' ') CK1(I:I) = '0'   
         IF( CK2(I:I).EQ.' ') CK2(I:I) = '0'
      ENDDO   
      DO I=1,2
         IF( CVLIM(I:I).EQ.' ') CVLIM(I:I) = '0'
         IF( CSLIM(I:I).EQ.' ') CSLIM(I:I) = '0'
      ENDDO 
      
      FICS = 'ROUJEAN-'

C -------- Building of the character string for K0
      FICS(9:11) = CK0(1:3)      ! Units digits
      FICS(12:12) = '.'	     ! Decimal point
      FICS(13:15) = CK0(4:6)     ! First three decimals
      FICS(16:16) = '-'

C -------- Building of the character string for K1
      FICS(17:19) = CK1(1:3)     ! Units digits
      FICS(20:20) = '.'	     ! Decimal point
      FICS(21:23) = CK1(4:6)     ! First three decimals
      FICS(24:24) = '-'

C -------- Building of the character string for K2
      FICS(25:27) = CK2(1:3)     ! Units digits
      FICS(28:28) = '.'	     ! Decimal point
      FICS(29:31) = CK2(4:6)     ! First three decimals
      FICS(32:32) = '-'

C -------- Building of the character string for the limit
C -------- on the viewing angle.
      FICS(33:34) = 'VL'	    
      FICS(35:36) = CVLIM(1:2)     ! Units digit
      FICS(37:37) = '-'	        
 
C -------- Building of the character string for the limit
C -------- on the solar angle.
      FICS(38:39) = 'SL'	    
      FICS(40:41) = CSLIM(1:2)     ! Units digit
      FICS(42:42) = '-'	        
                        
C -------- Completion of the end of the file name with MUnb-NomFicUser
      CALL SOS_COMPLETE_NOM_FIC_SURF(NBMU_GAUSS,FICANGLES_USER,FICS)
      
      END !SOS_NOM_FICROUJ
      


      
C==============================================================================
C PROCEDURE: SOS_NOM_FICRH
C ==========
C      This procedure generates the name of the Rondeaux - Herman file
C      from:
C	- the refractive index vegetation / air: IND
C       - The associated BRDF model for the intensity component:
C           ROUJEAN model with the parameters K0,K1 and K2.
C
C      The following formats are imposed:
C         IND with format F5.3.
C	  kx with format F7.3.
C
C      The name of the RH file must be:
C
C      RH-IND-ROUJEAN-K0-K1-K2-VL{CTE_TETAV_LIM_ROUJEAN}-
C        -SL{CTE_TETAS_LIM_ROUJEAN}-MUnb-NomFicUser 
C      with
C      - IND with format F5.3 (e.g. 1.500)
C      - kx with format F7.3 (e.g. 000.500)
C      - CTE_TETAV_LIM_ROUJEAN and CTE_TETAS_LIM_ROUJEAN respectively the
C        limit viewing and solar zenith angles for the BRDF computations
C        (in degrees).
C      - MUnb the number of Gauss angles (e.g. MU24)
C      - NomFicUser the name of the user file (without directory tree) (e.g. MyAngles)
C
C      IND is converted into a character string IND(1:4)
C      For this, Format 5.3 * 1000 => A4
C           The first character CIND(1) is the units digit of IND,
C           etc.
C      Kx is converted into a character string CKx(1:6)
C      For this, Format F7.3 * 1000 => A6
C        The first character CKx(1) is the hundreds digit of Kx,
C        The second character CKx(2) is the tens digit of Kx,
C        The third character CKx(3) is the units digit of Kx,
C           etc.
C
C
C Description of input / output parameters
C ------------------------------------------
C	NBMU_GAUSS (I4)	: (I) Number of Gauss angles used
C       FICANGLES_USER (CHARACTER*CTE_LENFIC2)	: (I) Name of the "user"
C                                                 angles file (complete path)
C	IND (double)	: Vegetation index
C	K0  (double)	: Parameter of the Roujean model
C	K1  (double)	: Parameter of the Roujean model
C	K2  (double)	: Parameter of the Roujean model
C	FICS (CHARACTER*CTE_LENFIC1) : (O) name of the RH file (without directory tree)
C
C 
C Constants used:
C --------------------
C     CTE_TETAV_LIM_ROUJEAN : Maximum viewing zenith angle for the BRDF
C                             computations.
C     CTE_TETAS_LIM_ROUJEAN : Maximum solar zenith angle for the BRDF
C                             computations.
C     CTE_LENFIC1 : Length of file names without directory tree.
C     CTE_LENFIC2 : Length of file names with directory tree.
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C None
C
C  Subroutines used:
C --------------------------
C  This program calls the routine:
C     - SOS_COMPLETE_NOM_FIC_SURF
C
C==============================================================================
      SUBROUTINE SOS_NOM_FICRH(NBMU_GAUSS,FICANGLES_USER,IND,K0,K1,K2,
     &                         FICS) 


      IMPLICIT NONE

C* Declaration of variables
C  -------------------------         
      DOUBLE PRECISION IND		! (I) Vegetation index.
      DOUBLE PRECISION K0,K1,K2 	! (I) Parameters of the Roujean model.
      CHARACTER*CTE_LENFIC2 FICANGLES_USER  ! (I) path and name of the angle file
                                            ! of the "user" angles
      CHARACTER*CTE_LENFIC1 FICS	! (O) name of the RH file.

      CHARACTER*6 CK0,CK1,CK2	! Character strings
      CHARACTER*4 CVLIM,CSLIM	! Character strings
      CHARACTER*4 CIND 		! Character string of the vegetation index.
      
      INTEGER*2 I
      INTEGER*4 NBMU_GAUSS  ! (I) Number of Gauss angles used
      			      	 

	    
C --------- Conversion of the input parameters into character strings
      WRITE(CIND,'(I4)') int(IND*1000)		
      WRITE(CK0,'(I6)') int(K0*1000)		
      WRITE(CK1,'(I6)') int(K1*1000)
      WRITE(CK2,'(I6)') int(K2*1000)
      WRITE(CVLIM,'(I2)') int(CTE_TETAV_LIM_ROUJEAN)
      WRITE(CSLIM,'(I2)') int(CTE_TETAS_LIM_ROUJEAN)
      
C --------- Blanks are replaced by 0
      DO I=1,4
	 IF( CIND(I:I).EQ.' ') CIND(I:I) = '0'
      ENDDO   
      DO I=1,6
         IF( CK0(I:I).EQ.' ') CK0(I:I) = '0'
         IF( CK1(I:I).EQ.' ') CK1(I:I) = '0'   
         IF( CK2(I:I).EQ.' ') CK2(I:I) = '0'
      ENDDO   
      DO I=1,2
         IF( CVLIM(I:I).EQ.' ') CVLIM(I:I) = '0'
         IF( CSLIM(I:I).EQ.' ') CSLIM(I:I) = '0'
      ENDDO 
  
      
      FICS = 'RH-'
      
C -------- Building of the character string for the index IND
      FICS(4:4) = CIND(1:1)     ! Units digit
      FICS(5:5) = '.'	    ! Decimal point
      FICS(6:8) = CIND(2:4)     ! First three decimals

C -------- Addition of the name of the associated BRDF file
      FICS(9:17) = '-ROUJEAN-'

C -------- Building of the character string for K0
      FICS(18:20) = CK0(1:3)     ! Units digits
      FICS(21:21) = '.'	     ! Decimal point
      FICS(22:24) = CK0(4:6)     ! First three decimals
      FICS(25:25) = '-'

C -------- Building of the character string for K1
      FICS(26:28) = CK1(1:3)     ! Units digits
      FICS(29:29) = '.'	     ! Decimal point
      FICS(30:32) = CK1(4:6)     ! First three decimals
      FICS(33:33) = '-'

C -------- Building of the character string for K2
      FICS(34:36) = CK2(1:3)     ! Units digits
      FICS(37:37) = '.'	     ! Decimal point
      FICS(38:40) = CK2(4:6)     ! First three decimals
      FICS(41:41) = '-'

C -------- Building of the character string for the limit
C -------- on the viewing angle.
      FICS(42:43) = 'VL'	    
      FICS(44:45) = CVLIM(1:2)     ! Units digit
      FICS(46:46) = '-'	        
 
C -------- Building of the character string for the limit
C -------- on the solar angle.
      FICS(47:48) = 'SL'	    
      FICS(49:50) = CSLIM(1:2)     ! Units digit
      FICS(51:51) = '-'	 
                   
C -------- Completion of the end of the file name with MUnb-NomFicUser
      CALL SOS_COMPLETE_NOM_FIC_SURF(NBMU_GAUSS,FICANGLES_USER,FICS)

      END !SOS_NOM_FICRH      
      


      
      
C==============================================================================
C PROCEDURE: SOS_NOM_FICBREON
C ==========
C      This procedure generates the name of the BREON file from:
C	- the refractive index soil / air: IND
C       - The associated BRDF model for the intensity component:
C           ROUJEAN model with the parameters K0,K1 and K2.
C
C      The following formats are imposed:
C         IND with format F5.3.
C	  kx with format F7.3.
C
C      The name of the Breon file must be:
C
C      BREON-IND-ROUJEAN-K0-K1-K2-VL{CTE_TETAV_LIM_ROUJEAN}-
C        -SL{CTE_TETAS_LIM_ROUJEAN}-MUnb-NomFicUser
C      with
C      - IND with format F5.3 (e.g. 1.500)
C      - kx with format F7.3 (e.g. 0.500)
C      - CTE_TETAV_LIM_ROUJEAN and CTE_TETAS_LIM_ROUJEAN respectively the
C        limit viewing and solar zenith angles for the BRDF computations
C        (in degrees).
C      - MUnb the number of Gauss angles (e.g. MU24)
C      - NomFicUser the name of the user file (without directory tree) (e.g. MyAngles)
C
C      IND is converted into a character string IND(1:4)
C      For this, Format 5.3 * 1000 => A4
C           The first character CIND(1) is the units digit of IND,
C           etc.
C      Kx is converted into a character string CKx(1:6)
C      For this, Format F7.3 * 1000 => A6
C        The first character CKx(1) is the hundreds digit of Kx,
C        The second character CKx(2) is the tens digit of Kx,
C        The third character CKx(3) is the units digit of Kx,
C           etc.
C
C Description of input / output parameters
C ------------------------------------------
C	NBMU_GAUSS (I4)	: (I) Number of Gauss angles used
C       FICANGLES_USER (CHARACTER*CTE_LENFIC2)	: (I) Name of the "user"
C                                                 angles file (complete path)
C	IND (double)	: Soil index
C	K0  (double)	: Parameter of the Roujean model
C	K1  (double)	: Parameter of the Roujean model
C	K2  (double)	: Parameter of the Roujean model
C	FICS (CHARACTER*CTE_LENFIC1) : (O) name of the Breon file (without directory tree)
C
C 
C Constants used:
C --------------------
C     CTE_TETAV_LIM_ROUJEAN : Maximum viewing zenith angle for the BRDF
C                             computations.
C     CTE_TETAS_LIM_ROUJEAN : Maximum solar zenith angle for the BRDF
C                             computations.
C     CTE_LENFIC1 : Length of file names without directory tree.
C     CTE_LENFIC2 : Length of file names with directory tree.
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C None
C
C  Subroutines used:
C --------------------------
C  This program calls the routine:
C     - SOS_COMPLETE_NOM_FIC_SURF
C
C==============================================================================
      SUBROUTINE SOS_NOM_FICBREON(NBMU_GAUSS,FICANGLES_USER,IND,
     &                            K0,K1,K2,FICS) 


      IMPLICIT NONE

C* Declaration of variables
C  -------------------------         
      DOUBLE PRECISION IND		! (I) Soil index.
      DOUBLE PRECISION K0,K1,K2 	! (I) Parameters of the Roujean model.
     
      CHARACTER*CTE_LENFIC2 FICANGLES_USER  ! (I) path and name of the "user" angles file
      CHARACTER*CTE_LENFIC1 FICS	! (O) name of the Breon file.

      CHARACTER*6 CK0,CK1,CK2	! Character strings
      CHARACTER*4 CVLIM,CSLIM	! Character strings
      CHARACTER*4 CIND 		! Character string of the vegetation index.
      
      INTEGER*2 I
      INTEGER*4 NBMU_GAUSS  ! (I) Number of Gauss angles used
	    
C --------- Conversion of the input parameters into character strings
      WRITE(CIND,'(I4)') int(IND*1000)		
      WRITE(CK0,'(I6)') int(K0*1000)		
      WRITE(CK1,'(I6)') int(K1*1000)
      WRITE(CK2,'(I6)') int(K2*1000)
      WRITE(CVLIM,'(I2)') int(CTE_TETAV_LIM_ROUJEAN)
      WRITE(CSLIM,'(I2)') int(CTE_TETAS_LIM_ROUJEAN)
      
C --------- Blanks are replaced by 0
      DO I=1,4
	 IF( CIND(I:I).EQ.' ') CIND(I:I) = '0'
      ENDDO      
      DO I=1,6
         IF( CK0(I:I).EQ.' ') CK0(I:I) = '0'
         IF( CK1(I:I).EQ.' ') CK1(I:I) = '0'   
         IF( CK2(I:I).EQ.' ') CK2(I:I) = '0'
      ENDDO   
      DO I=1,2
         IF( CVLIM(I:I).EQ.' ') CVLIM(I:I) = '0'
         IF( CSLIM(I:I).EQ.' ') CSLIM(I:I) = '0'
      ENDDO 
  
      
      FICS = 'BREON-'
      
C -------- Building of the character string for the index IND
      FICS(7:7) = CIND(1:1)      ! Units digit
      FICS(8:8) = '.'	     ! Decimal point
      FICS(9:11) = CIND(2:4)     ! First three decimals

C -------- Addition of the name of the associated BRDF file
      FICS(12:20) = '-ROUJEAN-'

C -------- Building of the character string for K0
      FICS(21:23) = CK0(1:3)     ! Units digits
      FICS(24:24) = '.'	     ! Decimal point
      FICS(25:27) = CK0(4:6)     ! First three decimals
      FICS(28:28) = '-'

C -------- Building of the character string for K1
      FICS(29:31) = CK1(1:3)     ! Units digits
      FICS(32:32) = '.'	     ! Decimal point
      FICS(33:35) = CK1(4:6)     ! First three decimals
      FICS(36:36) = '-'

C -------- Building of the character string for K2
      FICS(37:39) = CK2(1:3)     ! Units digits
      FICS(40:40) = '.'	     ! Decimal point
      FICS(41:43) = CK2(4:6)     ! First three decimals
      FICS(44:44) = '-'

C -------- Building of the character string for the limit
C -------- on the viewing angle.
      FICS(45:46) = 'VL'	    
      FICS(47:48) = CVLIM(1:2)     ! Units digit
      FICS(49:49) = '-'	        
 
C -------- Building of the character string for the limit
C -------- on the solar angle.
      FICS(50:51) = 'SL'	    
      FICS(52:53) = CSLIM(1:2)     ! Units digit
      FICS(54:54) = '-'	 
                   
C -------- Completion of the end of the file name with MUnb-NomFicUser
      CALL SOS_COMPLETE_NOM_FIC_SURF(NBMU_GAUSS,FICANGLES_USER,FICS)


      END !SOS_NOM_FICBREON     
      


 
 
      
C==============================================================================
C PROCEDURE: SOS_NOM_FICNADAL
C ==========
C      This procedure generates the name of the file of the NADAL model
C      from:
C	- the saturation value alpha: ALPHA
C	- the value of beta: BETA
C       - the refractive index surface / atmosphere: IND
C
C      The following formats are imposed:
C         IND with format F5.3.
C         ALPHA with format F6.4,
C	  BETA with format F4.1.
C
C      The name of the Nadal file must be:
C
C      NADAL-IND-K0-K1-K2-MUnb-NomFicUser with
C      - ALPHA with format F6.4 (e.g. 0.0229)
C      - BETA with format F3.1 (e.g. 38)
C      - MUnb the number of Gauss angles (e.g. MU24)
C      - NomFicUser the name of the user file (without directory tree) (e.g. MyAngles)
C
C      ALPHA is converted into a character string CALPHA(1:6)
C      For this, Format F6.4 * 10000 => A6
C      The first character CALPHA(1) is the units digit of ALPHA,
C           etc.
C
C      BETA is converted into a character string CBETA(1:4)
C      For this, Format F4.1 * 10 => A4
C      The first character CBETA(1) is the tens digit of BETA,
C           etc.
C
C Description of input / output parameters
C ------------------------------------------
C	NBMU_GAUSS (I4)	: (I) Number of Gauss angles used
C       FICANGLES_USER (CHARACTER*CTE_LENFIC2)	: (I) Name of the "user"
C                                                 angles file (complete path)
C       IND (double)	: Refractive index surface / atmosphere
C	ALPHA (double)	: Parameter alpha of the Nadal model
C	BETA  (double)	: Parameter beta of the Nadal model
C	K0  (double)	: Parameter of the Roujean model
C	K1  (double)	: Parameter of the Roujean model
C	K2  (double)	: Parameter of the Roujean model
C	FICS (CHARACTER*CTE_LENFIC1) : (O) name of the Nadal file (without directory tree)
C
C 
C Constants used:
C --------------------
C     CTE_TETAV_LIM_ROUJEAN : Maximum viewing zenith angle for the BRDF
C                             computations.
C     CTE_TETAS_LIM_ROUJEAN : Maximum solar zenith angle for the BRDF
C                             computations.
C     CTE_LENFIC1 : Length of file names without directory tree.
C     CTE_LENFIC2 : Length of file names with directory tree.
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C None
C
C  Subroutines used:
C --------------------------
C  This program calls the routine:
C     - SOS_COMPLETE_NOM_FIC_SURF
C
C==============================================================================
      SUBROUTINE SOS_NOM_FICNADAL(NBMU_GAUSS,FICANGLES_USER,IND,
     &                            ALPHA,BETA,K0,K1,K2,FICS) 


      IMPLICIT NONE

C* Declaration of variables
C  -------------------------  
      DOUBLE PRECISION IND	 ! (I) Soil index
      DOUBLE PRECISION ALPHA ! (I) Parameter alpha of the Nadal model
      DOUBLE PRECISION BETA	 ! (I) Parameter beta of the Nadal model
      DOUBLE PRECISION K0,K1,K2  ! (I) Parameters of the Roujean model.
      
      CHARACTER*CTE_LENFIC2 FICANGLES_USER  ! (I) path and name of the "user" angles file
      CHARACTER*CTE_LENFIC1 FICS ! (O) name of the NADAL file
          
      CHARACTER*5 CALPHA	! Character strings
      CHARACTER*3 CBETA 
      CHARACTER*4 CIND 
      CHARACTER*6 CK0,CK1,CK2	
      CHARACTER*4 CVLIM,CSLIM	
             
      INTEGER*2 I
      INTEGER*4 NBMU_GAUSS  ! (I) Number of Gauss angles used
      


C --------- Conversion of the input parameters into character strings
      WRITE(CIND,'(I4)') int(IND*1000)				
      WRITE(CALPHA,'(I5)') int(ALPHA*10000)
      WRITE(CBETA,'(I3)') int(BETA*10)	
      WRITE(CK0,'(I6)') int(K0*1000)		
      WRITE(CK1,'(I6)') int(K1*1000)
      WRITE(CK2,'(I6)') int(K2*1000)
      WRITE(CVLIM,'(I2)') int(CTE_TETAV_LIM_ROUJEAN)
      WRITE(CSLIM,'(I2)') int(CTE_TETAS_LIM_ROUJEAN)
            
    
C --------- Blanks are replaced by 0
      DO I=1,4
         IF( CIND(I:I).EQ.' ') CIND(I:I) = '0'
      ENDDO   
      	     
      DO I=1,5
         IF( CALPHA(I:I).EQ.' ') CALPHA(I:I) = '0'
      ENDDO         
      DO I=1,3
	 IF( CBETA(I:I).EQ.' ') CBETA(I:I) = '0'
      ENDDO   
      
      DO I=1,6
         IF( CK0(I:I).EQ.' ') CK0(I:I) = '0'
         IF( CK1(I:I).EQ.' ') CK1(I:I) = '0'   
         IF( CK2(I:I).EQ.' ') CK2(I:I) = '0'
      ENDDO   
      DO I=1,2
         IF( CVLIM(I:I).EQ.' ') CVLIM(I:I) = '0'
         IF( CSLIM(I:I).EQ.' ') CSLIM(I:I) = '0'
      ENDDO 



      
      FICS = 'NADAL-'
  
C -------- Building of the character string for the index IND
      FICS(7:7) = CIND(1:1)       ! Units digit
      FICS(8:8) = '.'	      ! Decimal point
      FICS(9:11) = CIND(2:4)      ! First three decimals
      FICS(12:12) = '-'
            
C -------- Building of the character string for alpha
      FICS(13:13) = CALPHA(1:1)     ! Units digit
      FICS(14:14) = '.'	        ! Decimal point
      FICS(15:18) = CALPHA(2:5)     ! First four decimals
      FICS(19:19) = '-'

C -------- Building of the character string for beta
      FICS(20:21) = CBETA(1:2)      ! Units digit
      FICS(22:22) = '.'
      FICS(23:23) = CBETA(3:3)      ! Units digit

C -------- Addition of the name of the associated BRDF file
      FICS(24:32) = '-ROUJEAN-'

C -------- Building of the character string for K0
      FICS(33:35) = CK0(1:3)     ! Units digits
      FICS(36:36) = '.'	     ! Decimal point
      FICS(37:39) = CK0(4:6)     ! First three decimals
      FICS(40:40) = '-'
		
C -------- Building of the character string for K1
      FICS(41:43) = CK1(1:3)     ! Units digits
      FICS(44:44) = '.'	     ! Decimal point
      FICS(45:47) = CK1(4:6)     ! First three decimals
      FICS(48:48) = '-'
		
C -------- Building of the character string for K2
      FICS(49:51) = CK2(1:3)     ! Units digits
      FICS(52:52) = '.'	     ! Decimal point
      FICS(53:55) = CK2(4:6)     ! First three decimals
      FICS(56:56) = '-'
		
C -------- Building of the character string for the limit
C -------- on the viewing angle.
      FICS(57:58) = 'VL'	    
      FICS(59:60) = CVLIM(1:2)     ! Units digit
      FICS(61:61) = '-'	        
 
C -------- Building of the character string for the limit
C -------- on the solar angle.
      FICS(62:63) = 'SL'	    
      FICS(64:65) = CSLIM(1:2)     ! Units digit
      FICS(66:66) = '-'	 
      
    			                 
C -------- Completion of the end of the file name with MUnb-NomFicUser
      CALL SOS_COMPLETE_NOM_FIC_SURF(NBMU_GAUSS,FICANGLES_USER,FICS)
   
      END !SOS_NOM_FICNADAL      


 

      
C==============================================================================
C PROCEDURE: SOS_NOM_FICMAIGNAN
C ==========
C      This procedure generates the name of the MAIGNAN file from:
C	- The refractive index soil / air: IND
C      - The coefficient C of the Maignan model
C      - The associated BRDF model for the intensity component:
C         ROUJEAN model with the parameters K0,K1 and K2.
C
C      The following formats are imposed:
C        IND with format F5.3.
C        C with format F7.3
C	  kx with format F7.3.
C
C      The name of the Maignan file must be:
C
C      MAIGNAN-IND-C-ROUJEAN-K0-K1-K2-VL{CTE_TETAV_LIM_ROUJEAN}-
C        -SL{CTE_TETAS_LIM_ROUJEAN}-MUnb-NomFicUser
C      with
C      - IND with format F5.3 (e.g. 1.500)
C      - C with format F7.3 (e.g. 007.290)
C      - kx with format F7.3 (e.g. 000.500)
C      - CTE_TETAV_LIM_ROUJEAN and CTE_TETAS_LIM_ROUJEAN respectively the
C        limit viewing and solar zenith angles for the BRDF computations
C        (in degrees).
C      - MUnb the number of Gauss angles (e.g. MU24)
C      - NomFicUser the name of the user file (without directory tree) (e.g. MyAngles)
C
C      IND is converted into a character string IND(1:4)
C      For this, Format 5.3 * 1000 => A4
C           The first character CIND(1) is the units digit of IND,
C           etc.
C      Kx is converted into a character string CKx(1:6)
C      For this, Format F7.3 * 1000 => A6
C        The first character CKx(1) is the hundreds digit of Kx,
C        The second character CKx(2) is the tens digit of Kx,
C        The third character CKx(3) is the units digit of Kx,
C           etc.
C      Same for the coefficient C.
C
C Description of input / output parameters
C ------------------------------------------
C	NBMU_GAUSS (I4)	: (I) Number of Gauss angles used
C       FICANGLES_USER (CHARACTER*CTE_LENFIC2)	: (I) Name of the "user"
C                                                 angles file (complete path)
C	IND (double)	  : (I) Refractive index surface / air
C      COEF_C (double) : (I) Parameter of the Maignan model
C	K0  (double)	: (I) Parameter of the Roujean model
C	K1  (double)	: (I) Parameter of the Roujean model
C	K2  (double)	: (I) Parameter of the Roujean model
C	FICS (CHARACTER*CTE_LENFIC1) : (O) name of the Maignan file (without directory tree)
C
C 
C Constants used:
C --------------------
C     CTE_TETAV_LIM_ROUJEAN : Maximum viewing zenith angle for the BRDF
C                             computations.
C     CTE_TETAS_LIM_ROUJEAN : Maximum solar zenith angle for the BRDF
C                             computations.
C     CTE_LENFIC1 : Length of file names without directory tree.
C     CTE_LENFIC2 : Length of file names with directory tree.
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C None
C
C  Subroutines used:
C --------------------------
C  This program calls the routine:
C     - SOS_COMPLETE_NOM_FIC_SURF
C
C==============================================================================
      SUBROUTINE SOS_NOM_FICMAIGNAN(NBMU_GAUSS,FICANGLES_USER,IND,
     &                              COEF_C,K0,K1,K2,FICS) 


      IMPLICIT NONE

C* Declaration of variables
C  -------------------------         
      DOUBLE PRECISION IND		! (I) Soil index.
      DOUBLE PRECISION COEF_C      ! (I) Parameter C*exp(-NDVI) of the Maignan model.
      DOUBLE PRECISION K0,K1,K2 	! (I) Parameters of the Roujean model.
     
      CHARACTER*CTE_LENFIC2 FICANGLES_USER  ! (I) path and name of the "user" angles file
      CHARACTER*CTE_LENFIC1 FICS	! (O) name of the Maignan file.

      CHARACTER*6 CK0,CK1,CK2	! Character strings
      CHARACTER*4 CVLIM,CSLIM	! Character strings
      CHARACTER*4 CIND 		! Character string of the vegetation index.
      CHARACTER*6 CC               ! Character strings for the coefficient C
      
      INTEGER*2 I
      INTEGER*4 NBMU_GAUSS  ! (I) Number of Gauss angles used
	    
C --------- Conversion of the input parameters into character strings
      WRITE(CIND,'(I4)') int(IND*1000)		
      WRITE(CC,'(I6)')  int(COEF_C*1000)		
      WRITE(CK0,'(I6)') int(K0*1000)		
      WRITE(CK1,'(I6)') int(K1*1000)
      WRITE(CK2,'(I6)') int(K2*1000)
      WRITE(CVLIM,'(I2)') int(CTE_TETAV_LIM_ROUJEAN)
      WRITE(CSLIM,'(I2)') int(CTE_TETAS_LIM_ROUJEAN)
      
C --------- Blanks are replaced by 0
      DO I=1,4
	 IF( CIND(I:I).EQ.' ') CIND(I:I) = '0'
      ENDDO      
      DO I=1,6
         IF( CC(I:I).EQ.' ') CC(I:I) = '0'
         IF( CK0(I:I).EQ.' ') CK0(I:I) = '0'
         IF( CK1(I:I).EQ.' ') CK1(I:I) = '0'   
         IF( CK2(I:I).EQ.' ') CK2(I:I) = '0'
      ENDDO   
      DO I=1,2
         IF( CVLIM(I:I).EQ.' ') CVLIM(I:I) = '0'
         IF( CSLIM(I:I).EQ.' ') CSLIM(I:I) = '0'
      ENDDO 
  
      
      FICS = 'MAIGNAN-'
      
C -------- Building of the character string for the index IND
      FICS(9:9) = CIND(1:1)       ! Units digit
      FICS(10:10) = '.'	      ! Decimal point
      FICS(11:13) = CIND(2:4)     ! First three decimals

C -------- Building of the character string for the coefficient C
      FICS(14:15) = '-C'
      FICS(16:18) = CC(1:3)      ! Units digits
      FICS(19:19) = '.'	     ! Decimal point
      FICS(20:22) = CC(4:6)      ! First three decimals


C -------- Addition of the name of the associated BRDF file
      FICS(23:31) = '-ROUJEAN-'

C -------- Building of the character string for K0
      FICS(32:34) = CK0(1:3)     ! Units digits
      FICS(35:35) = '.'	     ! Decimal point
      FICS(36:38) = CK0(4:6)     ! First three decimals
      FICS(39:39) = '-'

C -------- Building of the character string for K1
      FICS(40:42) = CK1(1:3)     ! Units digits
      FICS(43:43) = '.'	     ! Decimal point
      FICS(44:46) = CK1(4:6)     ! First three decimals
      FICS(47:47) = '-'

C -------- Building of the character string for K2
      FICS(48:50) = CK2(1:3)     ! Units digits
      FICS(51:51) = '.'	     ! Decimal point
      FICS(52:54) = CK2(4:6)     ! First three decimals
      FICS(55:55) = '-'

C -------- Building of the character string for the limit
C -------- on the viewing angle.
      FICS(56:57) = 'VL'	    
      FICS(58:59) = CVLIM(1:2)     ! Units digit
      FICS(60:60) = '-'	        
 
C -------- Building of the character string for the limit
C -------- on the solar angle.
      FICS(61:62) = 'SL'	    
      FICS(63:64) = CSLIM(1:2)     ! Units digit
      FICS(65:65) = '-'	 
                   
C -------- Completion of the end of the file name with MUnb-NomFicUser
      CALL SOS_COMPLETE_NOM_FIC_SURF(NBMU_GAUSS,FICANGLES_USER,FICS)


      END !SOS_NOM_FICMAIGNAN     
      



      
C==============================================================================
C PROCEDURE: SOS_COMPLETE_NOM_FIC_SURF
C ==========

C
C Description of input / output parameters
C ------------------------------------------
C	NBMU_GAUSS (I4)	: (I) Number of Gauss angles used
C       FICANGLES_USER (CHARACTER*CTE_LENFIC2)	: (I) Name of the "user"
C                                                 angles file (complete path)
C	FICS (CHARACTER*CTE_LENFIC1) : (I/O) name of the BRDF/BPDF file (without directory tree)
C 
C Constants used:
C --------------------
C     CTE_LENFIC1 : Length of file names without directory tree.
C     CTE_LENFIC2 : Length of file names with directory tree.
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C None
C==============================================================================
      SUBROUTINE SOS_COMPLETE_NOM_FIC_SURF(NBMU_GAUSS,FICANGLES_USER,
     &                                     FICS)

      IMPLICIT NONE

C* Declaration of variables
C  -------------------------  
      CHARACTER*CTE_LENFIC2 FICANGLES_USER  ! (I) path and name of the "user" angles file
      CHARACTER*CTE_LENFIC1 FICS            ! (I/O) name of the GLITTER file

      CHARACTER*4 CNBMU     ! Character string associated with NBMU_GAUSS
      
      INTEGER*2 I
      
      INTEGER*4 NBMU_GAUSS  ! (I) Number of Gauss angles used
      INTEGER*4 IDEB, IDEB_CH, IFIN, IFIN_CH  ! start and end positions of a character string
      
      			      	 

	    
C --------- Conversion of NBMU_GAUSS into a character string
      WRITE(CNBMU,'(I4)') NBMU_GAUSS      
      

C -------- Building of the character string for MUgauss
C --------Structured to have NBMU_GAUSS on the left (without "blanks")
      IFIN=INDEX(FICS,' ')
                  
      FICS(IFIN:IFIN+1) = 'MU'
      
      IF (NBMU_GAUSS.LT.10) THEN
          FICS(IFIN+2:IFIN+2) = CNBMU(4:4)	!MUn
      ELSE 
          IF (NBMU_GAUSS.LT.100) THEN
              FICS(IFIN+2:IFIN+3) = CNBMU(3:4)  !MUnn
	  ELSE 
	     IF (NBMU_GAUSS.LT.1000) THEN
                 FICS(IFIN+2:IFIN+4) = CNBMU(2:4) !MUnnn
	     ELSE !NBMU_GAUSS < 10000
	     	 FICS(IFIN+2:IFIN+5) = CNBMU(1:4) !MUnnnn
	     ENDIF !NBMU_GAUSS < 1000
	  ENDIF !NBMU_GAUSS < 100
      ENDIF !NBMU_GAUSS < 10
  
       
C -------- Building of the character string for NomFicUser
C -------- Structured to have NBMU_GAUSS on the left (without "blanks")


      IF (FICANGLES_USER.NE.'NO_USER_ANGLES') THEN

C       Determination of the position of the file name without its directory tree
	  IFIN_CH=INDEX(FICANGLES_USER,' ')
          IFIN_CH=IFIN_CH-1
	  IF (IFIN_CH.LE.0) IFIN_CH=CTE_LENDIR
	  
	  IDEB_CH=0
	  	  
10	  IDEB=INDEX(FICANGLES_USER(IDEB_CH:IFIN_CH),'/')
	  IDEB_CH=IDEB_CH+IDEB
          IF (IDEB.NE.0) GOTO 10	  

C       Addition of the name of the "user" angle file to the file name
          IFIN=INDEX(FICS,' ')  
	  FICS(IFIN:IFIN) = '-'
	  
	  IFIN=IFIN+1
	  FICS(IFIN:IFIN+IFIN_CH-IDEB_CH+1) 
     &	  = FICANGLES_USER(IDEB_CH:IFIN_CH)
	  
      ENDIF 
   
      END !SOS_COMPLETE_NOM_FIC_SURF   
