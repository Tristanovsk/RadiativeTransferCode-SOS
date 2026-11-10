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
C* FILE: SOS_SUB_TRS.F
C* PROJECT: Successive orders of scattering
C* ROLE: Library of functions for the atmospheric absorption computations
C*	 This file comes from the LOA TRS code, which computes the atmospheric
C*	 transmission due to absorption. It has been adapted to be included
C*	 in the OS
C*
C* AUTHOR:
C*   Initial code: Laboratoire d'Optique Atmospherique (LOA, Villeneuve d'Ascq).
C*   Code rework and development: CS GROUP France, X.Lenot, B.Lafrance
C* D
C*
C* MOD:VERSION:1.0: DATE: 24/06/2013 : File creation
C*
C* MOD:VERSION:1.1: DATE: 24/02/2016
C*    - Array NEXP changed from Double Precision to Integer
C*
C* MOD:VERSION:2.0: DATE: 08/06/2020
C*    - Adjustments for the use of new files of CKD coefficients,
C*      including the continuums together with the absorption lines
C*      ==> Removal of the routines related to the continuum data
C*          CONTINUUMNEW, CONT_H2O, CONT_O3, CONT_N2, CONT_O2, CONT_NO2, 
C*          CONT_SO2, O3CHAPW, ABSCNT, CONTHLBL, CNTCO2LBL and FUDGE.
C*      ==> Removal of the associated DATA BLOCKs.
C*      ==> Complete rework of the routines READ_CKD_COEFF and COEFF_ABS_CKD.
C*      ==> Change of the dimension of the array RO of the absorber amounts
C*          in routine DATATM to include NO2
C*
C*    - Removal of the obsolete routine ZEROIT.
C*
C* MOD:VERSION:2.1: DATE: 26/10/2021
C*    - Correction of the management of a surface pressure given by the user in routine DATATM.
C*
C*
C* MOD:VERSION:2.2: DATE: 07/07/2022
C*    - Removal of the use of the environment variable SOS_ABS_FIC.
C*    - Constants renamed from SOS_* to CTE_*
C*
C* MOD:VERSION:2.3: DATE: 16/06/2023
C*    - Change of the prefixes of the CKD files 'coef_xxx_' to give the gas in capital letters.
C*    - Parameters KDIS_KI_H2O, TAB_CONC_H2O and NB_CONC_H2O added as arguments of routine READ_CKD_COEFF.
C*    - Adjustments of routine READ_CKD_COEFF to read the H2O CKD data
C*    - Correction of routine READ_CKD_COEFF: the test on NUMAX_FILE and NUMIN_FILE changed from > 10000 to >= 10000
C*                                                to include 10000 in the case of a 5-digit encoding
C*
C*
C*
C* INCLUDED FUNCTIONS:
C*  - DATATM : reading of the absorbing gas profiles
C*  - COEFF_ABS_CKD : computation of the CKD coefficients for one exponential of an
C*                    absorber and a given spectral interval
C*  - READ_CKD_COEFF  : reading of the tables of CKD coefficients
C*  - TROPICA : reading of the "tropical" profile
C*  - MIDLASU : reading of the "Middle Latitude Summer" profile
C*  - MIDLAWI : reading of the "Middle Latitude Winter" profile
C*  - SUBSUMM : reading of the Sub-Arctic Summer profile
C*  - SUBWINT : reading of the Sub-Arctic Winter profile
C*  - USTAD62 : reading of the USST62 profile (U.S. Standard Atmosphere)
C* 
C*----------------------------------------------------------------------------
C* Definition of constants
C*---------------------------------------------------------------------------- 
C* Constants used:
C*    CTE_NBABS : Number of absorbers
C*    CTE_CKD_NUMAX : Maximum wavenumber covered by the CKD coeff files (in cm-1)
C*    CTE_CKD_NB_NU_PER_FILE : Number of wavelengths per CKD file
C*    CTE_CKD_NWVL_MAX : Maximum number of wavelengths in the CKD coeff tables
C*    CTE_CKD_NAI_MAX  : Maximum number of exponentials of the expansions
C*    CTE_CKD_NT_MAX   : Maximum number of temperatures
C*    CTE_CKD_NP_MAX   : Maximum number of pressures



C******************************************************************************

#include "SOS.h"
#define LENLINE 1000


C==============================================================================
C PROCEDURE: COEFF_ABS_CKD      
C ==========
C      Interpolates the coefficients ki(P,T) of one exponential, one gas and one wavelength
C      at the pressure and temperature given as input.
C
C      Linear interpolation on pressure.
C      Then Splint interpolation on temperature.
C
C Description of parameters
C --------------------------
C
C     NABS  (I2) : (I) Index of the absorber
C                        ==> for KDIS_KI(IT, IP_INF, IK, NABS, LAMB)
C
C     LAMB1 (I4) : (I) Index of the CKD interval
C			   
C     IK  (I4) : (I) Index of the exponential for the absorption of the gas considered
C
C     PRS (double) : (I) Mean pressure of the layer (mbar) for which ki is computed
C
C     TMP (double) : (I) Mean temperature of the layer (K) for which ki is computed
C
C     CONC_H2O (double) : (I) Mean H2O concentration (ppmv) of the layer for which ki is computed
C
C     TAB_PRES (double): (I) Array of the pressure values associated with the coefficients ki(P,T)
C                            of CTE_CKD_NP_MAX elements.
C     NB_PRES (integer): (I) Number of actual pressure values.
C
C     TAB_TEMP (double): (I) Array of the temperature values associated with the coefficients ki(P,T)
C                            of CTE_CKD_NT_MAX elements.
C     NB_TEMP (integer): (I) Number of actual temperature values.
C
C     TAB_CONC_H2O (double): (I) Array of the H2O concentration values associated with the coefficients ki(P,T,C)
C                                of water vapour, of CTE_CKD_NC_MAX elements (in ppmv)
C     NB_CONC_H2O (integer): (I) Number of actual concentration values.
C
C     KDIS_KI(double) : (I) Array of the coefficients ki of the exponentials of the decomposition, containing
C                           (CTE_CKD_NT_MAX * CTE_CKD_NP_MAX * CTE_CKD_NAI_MAX * CTE_NBABS * CTE_CKD_NWVL_MAX)
C                           elements covering wavelengths, species, orders i of the decomposition in
C                           exponentials, pressures and temperatures
C
C     KDIS_KI_H2O(double) : (I) Array of the coefficients ki of the exponentials of the water vapour decomposition,
C                               containing (CTE_CKD_NT_MAX * CTE_CKD_NP_MAX * CTE_CKD_NC_MAX * CTE_CKD_NAI_MAX *  CTE_CKD_NWVL_MAX)
C                               elements covering wavelengths, orders i of the decomposition in exponentials,
C                               H2O concentrations, pressures and temperatures
C
C  
C   Output:
C      XK (double) : (O) Value of ki(PRS,TMP) interpolated at the pressure PRS and temperature TMP
C                        from the tabulated values ki(P,T).
C
C      IER (I4) : (O) error code =0 if no error, =-1 otherwise
C
C Environment variable
C ------------------------
C   SOS_ABS_ROOT : defines the access path to the main directory of the SOS_ABS code
C   This environment variable is used to locate the files under $SOS_ABS_ROOT/fic
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C     Temperature or pressure out of bounds
C     Splint interpolation error
C     Negative absorption coefficient
C
C==============================================================================
      SUBROUTINE COEFF_ABS_CKD(NABS, LAMB, IK, PRS, TMP, CONC_H2O,
     &                         TAB_PRES, NB_PRES, TAB_TEMP, NB_TEMP,
     &                         TAB_CONC_H2O, NB_CONC_H2O,
     &                         KDIS_KI, KDIS_KI_H2O, XK, IER)
c
      IMPLICIT NONE
            
      INTEGER*2 NABS    ! (I) Index of the absorber
      INTEGER*4 LAMB	   ! (I) Band index of the table KDIS_KI of the CKD coef
      INTEGER*4 IK	   ! (I) exponential index
      
      DOUBLE PRECISION  PRS ! (I) Mean pressure of the layer (mbar)
      DOUBLE PRECISION  TMP ! (I) Mean temperature of the layer (K)
      DOUBLE PRECISION  CONC_H2O ! (I) Mean H2O concentration of a layer of the profile (ppmv)
      DOUBLE PRECISION  XK

      DOUBLE PRECISION  KDIS_KI(CTE_CKD_NT_MAX, CTE_CKD_NP_MAX, 
     &                          CTE_CKD_NAI_MAX, CTE_NBABS,  
     &                          CTE_CKD_NWVL_MAX)
                        ! (I) Values of the coefficients ki
                        !     per wavelength, species,
                        !     order i of the decomposition in exponentials,
                        !     pressure and temperature

      DOUBLE PRECISION  KDIS_KI_H2O(CTE_CKD_NT_MAX, CTE_CKD_NP_MAX, 
     &                              CTE_CKD_NC_MAX, CTE_CKD_NAI_MAX, 
     &                              CTE_CKD_NWVL_MAX )
                        ! (I) Values of the coefficients ki
                        !     per wavelength,
                        !     order i of the decomposition in exponentials,
                        !     concentration, pressure and temperature

      DOUBLE PRECISION  XKI(CTE_CKD_NT_MAX)     ! Values KDIS_KI(T) after linear interpolation on the ki values
                                                ! of the tabulated pressures around the layer pressure
                                                ! (for a fixed species NABS, a fixed order IK of the decomposition in exponentials
                                                !  and a fixed wavelength LAMB)

      DOUBLE PRECISION  XKI_H2O(CTE_CKD_NT_MAX, CTE_CKD_NP_MAX) 
                                                ! Values KDIS_KI(T,P) after linear interpolation on the ki values
                                                ! of the tabulated H2O concentrations around the layer concentration
                                                ! (for the species H2O, a fixed order IK of the decomposition in exponentials
                                                !  and a fixed wavelength LAMB)
                        
      INTEGER NB_TEMP, NB_PRES  ! (I) Number of temperature and pressure values in the CKD file
      DOUBLE PRECISION TAB_TEMP(CTE_CKD_NT_MAX) ! (I) Array of the temperature values read
      DOUBLE PRECISION TAB_PRES(CTE_CKD_NP_MAX) ! (I) Array of the pressure values read

      INTEGER NB_CONC_H2O  ! (I) Number of H2O concentration values in the H2O CKD file
      DOUBLE PRECISION TAB_CONC_H2O(CTE_CKD_NC_MAX) ! (I) Array of the H2O concentration values read

      DOUBLE PRECISION X1, X2, Y1, Y2 ! Pairs of values for linear interpolation
      
      INTEGER*4 IT_INF, IP_INF, IC_INF    ! Lower index of the bounds of the pressure, temperature
                                          ! and H2O concentration values in the tables of the CKD files

      INTEGER*4 IT,IP   ! Temperature and pressure indices
                                              
      INTEGER*4 IER     ! (O) Error code =0 if no error, =-1 otherwise
      

C------------------------------------------------------
C* Error code
      IER=0
      
C     Management of the extreme temperature cases:
C     Assumes that TAB_TEMP sorts the temperatures in increasing order
      IF(TAB_TEMP(NB_TEMP).LT.TAB_TEMP(1)) GOTO 900   ! error case
      IF(TMP.GT.TAB_TEMP(NB_TEMP)) TMP = TAB_TEMP(NB_TEMP)
      IF(TMP.LT.TAB_TEMP(1)) TMP = TAB_TEMP(1)

C     Management of the extreme pressure cases:
C     Assumes that TAB_PRES sorts the pressures in increasing order
C     ==> NEW WITH RESPECT TO THE INITIAL ROUTINE
      IF(TAB_PRES(NB_PRES).LT.TAB_PRES(1)) GOTO 910   ! error case
      IF(PRS.LE.TAB_PRES(1))THEN
         XK=0.
         GOTO 9999
      ENDIF
      IF(PRS.GT.TAB_PRES(NB_PRES)) PRS = TAB_PRES(NB_PRES)

C     Management of the extreme water vapour concentration cases:
C     Assumes that TAB_CONC_H2O sorts the concentrations in increasing order
      IF(TAB_CONC_H2O(NB_CONC_H2O).LT.TAB_CONC_H2O(1)) GOTO 915   ! error case
      IF(CONC_H2O.GT.TAB_CONC_H2O(NB_CONC_H2O)) THEN
         CONC_H2O = TAB_CONC_H2O(NB_CONC_H2O)
      ENDIF
      IF(CONC_H2O.LT.TAB_CONC_H2O(1)) CONC_H2O = TAB_CONC_H2O(1)


C     Search for the table indices around the value P
      IP_INF = 1
      DO WHILE((TAB_PRES(IP_INF).LE.PRS).AND.(IP_INF.LT.NB_PRES))
         IP_INF = IP_INF + 1
      ENDDO
      IP_INF = IP_INF - 1
      IF (IP_INF.EQ.NB_PRES) GOTO 920   ! error case

C     Estimation of the value XK (CKD coefficient interpolated for the conditions P,T,Concentration)
	
      IF (NABS.EQ.1) THEN     ! Specific case of H2O

         !Search for the table indices around the H2O concentration value
         IC_INF = 1
         DO WHILE((TAB_CONC_H2O(IC_INF).LE.CONC_H2O)
     &             .AND.(IC_INF.LT.NB_CONC_H2O))
            IC_INF = IC_INF + 1
         ENDDO
         IC_INF = IC_INF - 1
         IF (IC_INF.EQ.NB_CONC_H2O) GOTO 916   ! error case

         !Linear interpolation on the H2O concentration
         X1 = TAB_CONC_H2O(IC_INF)
         X2 = TAB_CONC_H2O(IC_INF+1) 
         DO IP=1,NB_PRES
            DO IT=1,NB_TEMP
               Y1 = KDIS_KI_H2O(IT, IP, IC_INF,   IK, LAMB)
               Y2 = KDIS_KI_H2O(IT, IP, IC_INF+1, IK, LAMB)       
               CALL SOS_INTERPOL(Y1,Y2,X1,X2,CONC_H2O,XKI_H2O(IT,IP))			
            ENDDO 
         ENDDO   

         !Linear interpolation on pressure
         X1 = TAB_PRES(IP_INF)
         X2 = TAB_PRES(IP_INF+1) 
         DO IT=1,NB_TEMP
            Y1 = XKI_H2O(IT,IP_INF)
            Y2 = XKI_H2O(IT,IP_INF+1)      
            CALL SOS_INTERPOL(Y1,Y2,X1,X2,PRS,XKI(IT))
         ENDDO 
             
      ELSE   ! Case of the gases other than H2O

         !Linear interpolation on pressure
         X1 = TAB_PRES(IP_INF)
         X2 = TAB_PRES(IP_INF+1) 
         DO IT=1,NB_TEMP
            Y1 = KDIS_KI(IT, IP_INF,   IK, NABS, LAMB)
            Y2 = KDIS_KI(IT, IP_INF+1, IK, NABS, LAMB)       
            CALL SOS_INTERPOL(Y1,Y2,X1,X2,PRS,XKI(IT))
         ENDDO     

      ENDIF ! End of block: case of the gases other than H2O
          
      !Spline interpolation on temperature ==> coef XK = searched ki
      CALL SOS_INTERPO_SPLINT(NB_TEMP,TAB_TEMP,XKI,1,TMP,XK,IER)
      IF (IER.NE.0) GOTO 922   ! error case

      !Physical consistency test of the value XK found: ki(P,T) >= 0 is required
      !Otherwise the interpolation on T is redone with a linear interpolator
      IF (XK.LT.0.) THEN

         IT_INF = 1
         DO WHILE((TAB_TEMP(IT_INF).LE.TMP).AND.(IT_INF.LT.NB_TEMP))
            IT_INF = IT_INF + 1
         ENDDO
         IT_INF = IT_INF - 1
         IF (IT_INF.EQ.NB_TEMP) GOTO 921   ! error case

         CALL SOS_INTERPOL(XKI(IT_INF),XKI(IT_INF+1),
     &                     TAB_TEMP(IT_INF),TAB_TEMP(IT_INF+1),TMP,XK)

         IF (XK.LT.0.) GOTO 923 ! Error case
      ENDIF ! End of processing if XK > 0.

     

C* Nominal end
C----------------------------------------------      
      GOTO 9999
      

C* Error case and return of status 1 to the shell
C----------------------------------------------

  900 WRITE(6,*) '  COEFF_ABS_CKD : ERROR_900 : '
      WRITE(6,*) '  --> Temperatures are not ascending values'
      WRITE(6,*) '  --> Unconsistent content of CKD datafiles'
      GOTO 9998

  910 WRITE(6,*) '  COEFF_ABS_CKD : ERROR_910 : '
      WRITE(6,*) '  --> Pressures are not ascending values'
      WRITE(6,*) '  --> Unconsistent content of CKD datafiles'
      GOTO 9998      

  915 WRITE(6,*) '  COEFF_ABS_CKD : ERROR_915 : '
      WRITE(6,*) '  --> H2O concentrations are not ascending values'
      WRITE(6,*) '  --> Unconsistent content of CKD datafiles'
      GOTO 9998  

  916 WRITE(6,*) '  COEFF_ABS_CKD : ERROR_916 : '
      WRITE(6,*) '  --> H2O concentration value in the atmospheric '
      WRITE(6,*) '      profile unconsistent with the range of '
      WRITE(6,*) '      concentration covered by the H2O CKD datafiles'
      GOTO 9998    
     
  920 WRITE(6,*) '  COEFF_ABS_CKD : ERROR_920 : '
      WRITE(6,*) '  --> Pressure value in the atmospheric profile'
      WRITE(6,*) '      unconsistent with the range of pressures '
      WRITE(6,*) '      covered by CKD datafiles'
      GOTO 9998      

  921 WRITE(6,*) '  COEFF_ABS_CKD : ERROR_921 : '
      WRITE(6,*) '  --> Temperature value in the atmospheric profile'
      WRITE(6,*) '      unconsistent with the range of temperatures '
      WRITE(6,*) '      covered by CKD datafiles'
      GOTO 9998      

  922 WRITE(6,*) '  COEFF_ABS_CKD : ERROR_922 : ' 
      WRITE(6,*) 
     &   '  --> ERROR while interpolating ki coeff on Temperature'
      GOTO 9998     

  923 WRITE(6,*) '  COEFF_ABS_CKD : ERROR_923 : ' 
      WRITE(6,*) 
     &   '  --> Calculations give ki < 0 : uncorrect value!'
      GOTO 9998      

 9998 IER=-1
      GOTO 9999
      
      
9999  RETURN
      END
c-----------------------------------------------------------------------
C



C==============================================================================
C PROCEDURE: READ_CKD_COEFF      
C ==========
C      Reading of the number of exponentials of the CKD decomposition of each absorber
C      for the wavelength range covered.
C      
C      The coefficients Ki and Ai are read for each exponential.
C
C      The CKD data files are stored under fic/
C
C      For a given wavelength:
C          If the nb of exponentials is 0: there is no absorption
C          If the nb of exponentials is 1: only a continuum is given
C          If the nb of exponentials > 1: the data of an absorption line are given,
C                                          possibly coupled with a continuum
C               The values ai and ki(P,T) are read.
C      
C
C Description of parameters
C --------------------------
C  
C   Inputs:
C   -------
C     NABS  (I2) : (I) Array index identifying the absorber.
C     JABS  (I2) : (I) Index specifying whether the gas is to be processed (1 to process it, 0 otherwise)
C                      If it is not to be processed:
C                         - the nb of exponentials = 1,
C                         - the value of the distribution function ai =1
C                         - the coefficients ki(P,T) = 0
C                      If it is to be processed, these coefficients are read.
C
C     NU (double) : (I) Wavenumber (cm-1) of the simulation.
C
C     NUSTEP (double) : (I) Spectral resolution in cm-1
C                          of the CKD coefficient files to read.
C  
C   Outputs:
C   -------
C     NEXP (I4) :   (O)  Array of (CTE_NBABS * CTE_CKD_NWVL_MAX) elements containing, for each absorber
C		           and each spectral interval, the number of exponentials of the CKD method/
C
C     NUMAX, NUMIN (double) : (O) Max and min wavenumber (cm-1) of the spectral coverage of the CKD file.
C
C     KDIS_AI (double) : (O) Array of (CTE_CKD_NAI_MAX * CTE_NBABS * CTE_CKD_NWVL_MAX) elements containing,
C			        for each wavelength, each species and each exponential of the decomposition,
C			        the value of the weighting coefficient ai of the exponential.
C
C     KDIS_KI(double) : (O) Array of the coefficients ki of the exponentials of the decomposition, containing
C                           (CTE_CKD_NT_MAX * CTE_CKD_NP_MAX * CTE_CKD_NAI_MAX * CTE_NBABS * CTE_CKD_NWVL_MAX)
C                           elements covering wavelengths, species, orders i of the decomposition in
C
C     KDIS_KI_H2O(double) : (O) Array of the coefficients ki of the exponentials of the water vapour decomposition,
C                               containing (CTE_CKD_NT_MAX * CTE_CKD_NP_MAX * CTE_CKD_NC_MAX * CTE_CKD_NAI_MAX *  CTE_CKD_NWVL_MAX)
C                               elements covering wavelengths, orders i of the decomposition in exponentials,
C                               H2O concentrations, pressures and temperatures
C
C     TAB_PRES (double): (O) Array of the pressure values associated with the coefficients ki(P,T)
C                            of CTE_CKD_NP_MAX elements.
C     NB_PRES (integer): (O) Number of actual pressure values
C
C     TAB_TEMP (double): (O) Array of the temperature values associated with the coefficients ki(P,T)
C                            of CTE_CKD_NT_MAX elements.
C     NB_TEMP (integer): (O) Number of actual temperature values
C
C     TAB_CONC_H2O (double): (O) Array of the H2O concentration values associated with the coefficients ki(P,T,C)
C                                of water vapour, of CTE_CKD_NC_MAX elements.
C
C     NB_CONC_H2O (integer): (O) Number of actual concentration values
C
C     IER (I4) : (O) error code =0 if no error, =-1 otherwise
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C     Spectral step of the files larger than 100.
C     Gas species number out of bounds.
C     Inconsistency of the file contents on the spectral step and the min and max spectral bounds.
C
C==============================================================================
      SUBROUTINE READ_CKD_COEFF (NABS,JABS,NU,NUSTEP,
     &                           NEXP,KDIS_AI,KDIS_KI,KDIS_KI_H2O,
     &                           NUMAX,NUMIN,
     &                           TAB_PRES,NB_PRES,TAB_TEMP,NB_TEMP,
     &                           TAB_CONC_H2O,NB_CONC_H2O,IER)
c


      IMPLICIT NONE

C* Definition of parameters

      INTEGER*2 NABS    ! (I) Index of the absorber
      INTEGER*2 JABS    ! (I) Status of the absorber (1 if processed, 0 otherwise)


      DOUBLE PRECISION NU  ! (I) Wavenumber (cm-1)
      
      DOUBLE PRECISION NUSTEP     
                        ! (I) Step in cm-1 of the CKD coefficient files
         
      INTEGER*4  NEXP(CTE_NBABS,CTE_CKD_NWVL_MAX)   
                        ! (O) Nb of exponentials of the CKD expansion
                        !     per wavelength and species
                        
      DOUBLE PRECISION  KDIS_AI(CTE_CKD_NAI_MAX, CTE_NBABS,  
     &                          CTE_CKD_NWVL_MAX) 
                        ! (O) Values of the coefficients ai
                        !     per wavelength, species
                        !     and order i of the decomposition in exponentials
                                          
      DOUBLE PRECISION  KDIS_KI(CTE_CKD_NT_MAX, CTE_CKD_NP_MAX, 
     &                          CTE_CKD_NAI_MAX, CTE_NBABS,  
     &                          CTE_CKD_NWVL_MAX)
                        ! (O) Values of the coefficients ki
                        !     per wavelength, species,
                        !     order i of the decomposition in exponentials,
                        !     pressure and temperature

      DOUBLE PRECISION  KDIS_KI_H2O(CTE_CKD_NT_MAX, CTE_CKD_NP_MAX, 
     &                              CTE_CKD_NC_MAX, CTE_CKD_NAI_MAX, 
     &                              CTE_CKD_NWVL_MAX )
                        ! (O) Values of the coefficients ki
                        !     per wavelength,
                        !     order i of the decomposition in exponentials,
                        !     concentration, pressure and temperature
      
      DOUBLE PRECISION NUMAX       ! (O) Max wavenumber (cm-1) covered by the CKD file
      DOUBLE PRECISION NUMIN       ! (O) Min wavenumber covered by the CKD file
      DOUBLE PRECISION RESOLUTION  ! Spectral resolution (cm-1)

      INTEGER*4 NUMAX_FILE  ! Max wavenumber bound (cm-1) of the file to read
      INTEGER*4 NUMIN_FILE  ! Min wavenumber bound (cm-1) of the file to read

      CHARACTER*CTE_LENDIR DIRE	    ! Length of the directory names
      CHARACTER*CTE_LENFIC1 FICDATA ! File name without directory tree
            
      INTEGER*4 IFIN       		! nb of significant characters of a string
      INTEGER*4 STRLEN       	! nb of characters of a string
      
      CHARACTER*50 STR_ESPECE	! Character strings of the gas type
      
      INTEGER*4 IER     ! error code =0 if no error, =-1 otherwise
     
      INTEGER IWA, IP, IT, IC ! Array index on wavelength, pressure, temperature and H2O concentration
      INTEGER LIG             !Line index for file reading
      INTEGER I
      
      CHARACTER*LENLINE LIGNE_TEXTE     !Character line for decoding parameters
  
      INTEGER NB_TEMP, NB_PRES  ! (O) Number of temperature and pressure values in the CKD file
      DOUBLE PRECISION TAB_TEMP(CTE_CKD_NT_MAX) ! (O) Array of the temperature values read
      DOUBLE PRECISION TAB_PRES(CTE_CKD_NP_MAX) ! (O) Array of the pressure values read

      INTEGER NB_CONC_H2O  ! (O) Number of H2O concentration values in the H2O CKD file
      DOUBLE PRECISION TAB_CONC_H2O(CTE_CKD_NC_MAX) ! (O) Array of the H2O concentration values read
        
      DOUBLE PRECISION WAMIN_INTER, WAMOY_INTER, WAMAX_INTER  ! Min, mean and max wavelengths (microns) of the spectral interval read
      DOUBLE PRECISION NUMIN_INTER, NUMAX_INTER               ! Min and max wavenumbers (cm-1) of the spectral interval read
      INTEGER NMAXAI_INTER             ! Number of exponentials for the spectral interval read
      INTEGER NAI                      ! Order of the decomposition in exponentials for the spectral interval read

      INTEGER NB_WA  ! Number of wavelengths in the CKD file
      INTEGER NAI_LU, NIC_LU, NIP_LU !Values read in the CKD files (exponential index, H2O concentration, pressure)
      INTEGER NB_LIG_HEADER ! Number of header lines of the CKD file
      
C* Initialisation
C------------------------------------------------------

C* Error code
      IER=0
      
      IF((NUSTEP.NE.1).AND.(NUSTEP.NE.5).AND.(NUSTEP.NE.10)) GOTO 905 !Error case

C* Check of the gas species index
      IF ((NABS.LT.1).OR.(NABS.GT.8)) GOTO 900

C* Test whether the processed species belongs to the gases selected for the simulation
      IF(JABS.EQ.0) THEN
         !Case where the species is not to be processed: parameters adjusted for a null weight
         !-----------------------------------

         IF (NABS.EQ.1) THEN  !Specific case of H2O
            DO IWA=1,CTE_CKD_NWVL_MAX
               NEXP(NABS,IWA)=1
               KDIS_AI(1,NABS,IWA) = 1.0
               DO IC=1,CTE_CKD_NC_MAX
                  DO IP=1,CTE_CKD_NP_MAX
                     DO IT=1,CTE_CKD_NT_MAX
                        KDIS_KI_H2O(IT, IP, IC, 1, IWA) = 0.0
                     ENDDO
                  ENDDO
               ENDDO
             ENDDO   

         ELSE  !Case of the gases other than H2O
            DO IWA=1,CTE_CKD_NWVL_MAX
               NEXP(NABS,IWA)=1
               KDIS_AI(1,NABS,IWA) = 1.0
               DO IP=1,CTE_CKD_NP_MAX
                  DO IT=1,CTE_CKD_NT_MAX
                     KDIS_KI(IT, IP, 1, NABS, IWA) = 0.0
                  ENDDO
               ENDDO
            ENDDO   
         ENDIF !End of test NABS=1
    
        
      ELSE ! End of the block with marker FIN_1
         !General case where the species is to be processed: reading of the CKD coefficients
         !-------------------------------------
         
         
C  Retrieval of the name of the data directory

         CALL GETENV('SOS_ABS_ROOT',DIRE)
         IFIN=INDEX(DIRE,' ')
         IFIN=IFIN-1
         IF (IFIN.LE.0) GOTO 925
      
         IF ((IFIN+4).GT.CTE_LENDIR) GOTO 926
      
         DIRE(IFIN+1:IFIN+4)='/fic'
         IFIN=IFIN+4
         IF (IFIN.LE.0) IFIN=CTE_LENDIR


         WRITE(DIRE(IFIN+1:IFIN+11),'(A11)') "/COEFF_CKD/"

         IF(NUSTEP.LT.10) THEN !Coded on one character
            WRITE(DIRE(IFIN+12:IFIN+12),'(I1)') int(NUSTEP)
            IFIN = IFIN + 12
         ELSE !Coded on two characters
            WRITE(DIRE(IFIN+12:IFIN+13),'(I2)') int(NUSTEP)
            IFIN = IFIN + 13  
         ENDIF

         WRITE(DIRE(IFIN+1:IFIN+5),'(A5)') "cmm1/"
         IFIN = IFIN + 5


C  Name of the file to read, as 'coef_H2O_11500_9250_10cmm1' or 'coef_H2O_16000_13750_5cmm1'
         IF(NABS.EQ.1) STR_ESPECE='coef_H2O_'
         IF(NABS.EQ.2) STR_ESPECE='coef_CO2_'
         IF(NABS.EQ.3) STR_ESPECE='coef_O3_'
         IF(NABS.EQ.4) STR_ESPECE='coef_N2O_'
         IF(NABS.EQ.5) STR_ESPECE='coef_CO_'
         IF(NABS.EQ.6) STR_ESPECE='coef_CH4_'
         IF(NABS.EQ.7) STR_ESPECE='coef_O2_'
         IF(NABS.EQ.8) STR_ESPECE='coef_NO2_'

         STRLEN=9
         IF((NABS.EQ.3).OR.(NABS.EQ.5).OR.(NABS.EQ.7)) STRLEN=8        
         FICDATA(1:STRLEN)=STR_ESPECE(1:STRLEN)

         !Determination of the spectral bounds of the file to read
         NUMIN_FILE = CTE_CKD_NUMAX - CTE_CKD_NB_NU_PER_FILE*NUSTEP
         DO WHILE (NUMIN_FILE.GT.NU)  
            NUMIN_FILE = NUMIN_FILE - CTE_CKD_NB_NU_PER_FILE*NUSTEP
         ENDDO
         NUMAX_FILE = NUMIN_FILE + CTE_CKD_NB_NU_PER_FILE*NUSTEP

         IF (NUMAX_FILE.GE.10000) THEN  !Assumed below 100000
             WRITE(FICDATA(STRLEN+1:STRLEN+5),'(I5)') int(NUMAX_FILE)
             STRLEN = STRLEN+5
         ELSE !Assumed above 1000
             WRITE(FICDATA(STRLEN+1:STRLEN+4),'(I4)') int(NUMAX_FILE)
             STRLEN = STRLEN+4
         ENDIF

         WRITE(FICDATA(STRLEN+1:STRLEN+1),'(A1)') "_"
         STRLEN = STRLEN+1

         IF (NUMIN_FILE.GE.10000) THEN  !Assumed below 100000
             WRITE(FICDATA(STRLEN+1:STRLEN+5),'(I5)') int(NUMIN_FILE)
             STRLEN = STRLEN+5
         ELSE !Assumed above 1000
             WRITE(FICDATA(STRLEN+1:STRLEN+4),'(I4)') int(NUMIN_FILE)
             STRLEN = STRLEN+4
         ENDIF

         WRITE(FICDATA(STRLEN+1:STRLEN+1),'(A1)') "_"
         STRLEN = STRLEN+1

         IF(NUSTEP.LT.10) THEN !Coded on one character
            WRITE(FICDATA(STRLEN+1:STRLEN+1),'(I1)') int(NUSTEP)
            STRLEN = STRLEN + 1
         ELSE !Coded on two characters
            WRITE(FICDATA(STRLEN+1:STRLEN+2),'(I2)') int(NUSTEP)
            STRLEN = STRLEN + 2  
         ENDIF
               
         FICDATA(STRLEN+1:STRLEN+4)='cmm1'
         STRLEN=STRLEN+4       

C
C  Reading of the coefficients
         OPEN(1,FILE=DIRE(1:IFIN)//FICDATA(1:STRLEN), 
     &        STATUS='OLD', ERR = 910)

         !Header reading
         IF (NABS.EQ.1) THEN
            NB_LIG_HEADER=21
         ELSE
            NB_LIG_HEADER=18
         ENDIF
         DO LIG=1,NB_LIG_HEADER
            READ(1,'(a)',ERR=911) LIGNE_TEXTE
         ENDDO

         !Reading of the information on the spectral domain
         READ(1,*,ERR=911) NUMAX, NUMIN, RESOLUTION
         IF (RESOLUTION.NE.NUSTEP) GOTO 9120
         IF (NUMAX.NE.NUMAX_FILE)  GOTO 9121
         IF (NUMIN.NE.NUMIN_FILE)  GOTO 9122

         NB_WA = int((NUMAX-NUMIN)/RESOLUTION)
         IF(NB_WA.GT.CTE_CKD_NWVL_MAX) GOTO 9123

         !Reading of the information on the temperature sampling
         READ(1,*,ERR=913) NB_TEMP         
         READ(1,*,ERR=913) (TAB_TEMP(I), I=1,NB_TEMP )      

         !Reading of the information on the pressure sampling
         READ(1,*,ERR=914) NB_PRES         
         READ(1,*,ERR=914) (TAB_PRES(I), I=1,NB_PRES )



         IF (NABS.EQ.1) THEN  !Specific case of H2O

            !Reading of the information on the H2O concentration sampling
            READ(1,*,ERR=914) NB_CONC_H2O     
            READ(1,*,ERR=914) (TAB_CONC_H2O(I), I=1,NB_CONC_H2O )

            !Reading of the information on the ai and ki data: H2O case
            DO IWA=1,NB_WA
               READ(1,*,ERR=915) WAMIN_INTER, WAMOY_INTER, WAMAX_INTER,
     &                           NUMIN_INTER, NUMAX_INTER, NMAXAI_INTER

               IF (NMAXAI_INTER.EQ.0) THEN    !Case without absorption
                  NEXP(NABS,IWA) = 1          !--> Equivalent to one exponential with
                  KDIS_AI(1,NABS,IWA) = 1.0   !--> weight ai = 1 and value = 1, hence ki = 0 (transmission=1)
                  
                  DO IC=1,NB_CONC_H2O 
                     DO IP=1,NB_PRES 
                        DO IT=1,NB_TEMP
                           KDIS_KI_H2O(IT, IP, IC, 1, IWA) = 0.0
                        ENDDO
                     ENDDO
                  ENDDO
               
               ELSE        !Case with absorption
            
                  NEXP(NABS,IWA) = NMAXAI_INTER
                  READ(1,*,ERR=916) 
     &                (KDIS_AI(NAI,NABS,IWA),NAI=1,NMAXAI_INTER)
            
                  DO NAI=1,NMAXAI_INTER
                     DO IC=1, NB_CONC_H2O 
                        DO IP=1, NB_PRES
			      READ(1,*,ERR=917)  NAI_LU, NIC_LU, NIP_LU,
     &                     (KDIS_KI_H2O(IT,IP,IC,NAI,IWA),IT=1,NB_TEMP)
                        ENDDO
                     ENDDO
                  ENDDO                        
                                                
               ENDIF !End of case without/with absorption
         
            ENDDO

         ELSE   ! Case of the gases other than H2O (end of the block with marker FIN_2):
        
            !Reading of the information on the ai and ki data
            DO IWA=1,NB_WA
               READ(1,*,ERR=915) WAMIN_INTER, WAMOY_INTER, WAMAX_INTER,
     &                           NUMIN_INTER, NUMAX_INTER, NMAXAI_INTER

               IF (NMAXAI_INTER.EQ.0) THEN    !Case without absorption
                  NEXP(NABS,IWA) = 1          !--> Equivalent to one exponential with
                  KDIS_AI(1,NABS,IWA) = 1.0   !--> weight ai = 1 and value = 1, hence ki = 0 (transmission=1)
                  DO IP=1,NB_PRES
                     DO IT=1,NB_TEMP
                        KDIS_KI(IT, IP, 1, NABS, IWA) = 0.0
                     ENDDO
                  ENDDO
               
               ELSE        !Case with absorption
            
                  NEXP(NABS,IWA) = NMAXAI_INTER
                  READ(1,*,ERR=916) 
     &                (KDIS_AI(NAI,NABS,IWA),NAI=1,NMAXAI_INTER)
            
                  DO NAI=1,NMAXAI_INTER
                     DO IP=1, NB_PRES
			   READ(1,*,ERR=917)  NAI_LU, NIP_LU, 
     &                  (KDIS_KI(IT, IP, NAI, NABS, IWA), IT=1,NB_TEMP)
                     ENDDO
                  ENDDO                        
                                                
               ENDIF !End of case without/with absorption
         
            ENDDO

         ENDIF !Marker FIN_2: end of test NABS=1
                   
 2       CLOSE(1)
      ENDIF !Marker FIN_1: end of the test whether the processed species belongs to the selected gases
      
      
C* Nominal end
C-------------------
      GOTO 9999
      
C* Error case and return of status 1 to the shell
C----------------------------------------------

  900 WRITE(6,*) '  SOS_SUB_TRS : ERROR_900: '
      WRITE(6,*) '  Index of gas not in the range [1,8]'
      WRITE(6,*) '  --> Error in the code'
      GOTO 9998

  905 WRITE(6,*) '  SOS_SUB_TRS : ERROR_905: '
      WRITE(6,*) '  The required spectral resolution is not supported :'
      WRITE(6,*) NUSTEP, 'cm-1'
      WRITE(6,*) '  --> See available files in repertory fic/'
      GOTO 9998

 925  WRITE(6,*) '  SOS_SUB_TRS : ERROR_925: '
      WRITE(6,*) '  => Error while getting SOS_ABS_ROOT variable'
      GOTO 9998

 926  WRITE(6,*) '  SOS_SUB_TRS : ERROR_926: '
      WRITE(6,*) '  => Length of $SOS_ABS_ROOT/fic is too long' 
      WRITE(6,*) '  => CTE_LENDIR has to be adjusted.' 
      GOTO 9998
      
  910 WRITE(6,*) 'Error while opening the file of CKD coefficients'
      WRITE(6,*) 'File :',DIRE(1:IFIN)//FICDATA(1:STRLEN)
      GOTO 9998

  911 WRITE(6,*) 'Error while reading the file of CKD coefficients'
      WRITE(6,*) 'File :',DIRE(1:IFIN)//FICDATA(1:STRLEN)
      GOTO 9998

 9120 WRITE(6,*) 'Not consistent value of spectral resolution in'
      WRITE(6,*) 'file :',DIRE(1:IFIN)//FICDATA(1:STRLEN)
      WRITE(6,*) 'Value in file is (cm-1) : ', RESOLUTION
      WRITE(6,*) 'Requested value is (cm-1) : ', NUSTEP
      GOTO 9998

 9121 WRITE(6,*) 'Not consistent value of maximal wave number in'
      WRITE(6,*) 'file :',DIRE(1:IFIN)//FICDATA(1:STRLEN)
      WRITE(6,*) 'Value in file is (cm-1) : ', NUMAX
      WRITE(6,*) 'Expected value is (cm-1) : ', NUMAX_FILE
      GOTO 9998

 9122 WRITE(6,*) 'Not consistent value of minimal wave number in'
      WRITE(6,*) 'file :',DIRE(1:IFIN)//FICDATA(1:STRLEN)
      WRITE(6,*) 'Value in file is (cm-1) : ', NUMIN
      WRITE(6,*) 'Expected value is (cm-1) : ', NUMIN_FILE
      GOTO 9998

 9123 WRITE(6,*) 'Not consistent value of CTE_CKD_NWVL_MAX in fic/SOS.h'
      WRITE(6,*) 'with respect to the number of wavelengths in '
      WRITE(6,*) 'absorption spectral data files:'
      WRITE(6,*) '- Max number of wavelengths (CTE_CKD_NWVL_MAX)  : ', 
     &            CTE_CKD_NWVL_MAX 
      WRITE(6,*) '- Number of wavelengths in files : ', NB_WA
      GOTO 9998

  913 WRITE(6,*) 'Error while reading the file of CKD coefficients'
      WRITE(6,*) 'file :',DIRE(1:IFIN)//FICDATA(1:STRLEN)
      WRITE(6,*) '--> Block of temperature data'
      GOTO 9998
 
  914 WRITE(6,*) 'Error while reading the file of CKD coefficients'
      WRITE(6,*) 'file :',DIRE(1:IFIN)//FICDATA(1:STRLEN)
      WRITE(6,*) '--> Block of pressure data'
      GOTO 9998
 
  915 WRITE(6,*) 'Error while reading the file of CKD coefficients'
      WRITE(6,*) 'file :',DIRE(1:IFIN)//FICDATA(1:STRLEN)
      WRITE(6,*) '--> Spectral range information'
      GOTO 9998
 
  916 WRITE(6,*) 'Error while reading the file of CKD coefficients'
      WRITE(6,*) 'file :',DIRE(1:IFIN)//FICDATA(1:STRLEN)
      WRITE(6,*) '--> Coeff ai'
      GOTO 9998
 
  917 WRITE(6,*) 'Error while reading the file of CKD coefficients'
      WRITE(6,*) 'file :',DIRE(1:IFIN)//FICDATA(1:STRLEN)
      WRITE(6,*) '--> Coeff ki'
      GOTO 9998        
          
 9998 IER=-1
                           
9999  RETURN
      END
C* End of the routine READ_CKD_COEFF


     

 
c----------
      SUBROUTINE DATATM(RO,P,T,ALT,IATM,DENS,DONUSER,NLEVEL,PSURF)
C
C  *********************************************************************
C  * SUBPROGRAM FOR READING FROM 'CARDS' OR DISK THE PRESSURES, TEMPERATURES*
C  *                      AND CONCENTRATIONS                            *
C  *********************************************************************
C
      PARAMETER(NCORPS=CTE_NBABS,MAXCLY=50)
      
      DOUBLE PRECISION DONUSER(MAXCLY,13)
      DOUBLE PRECISION ALT(MAXCLY), RO(NCORPS,MAXCLY)
      DOUBLE PRECISION P(MAXCLY),T(MAXCLY)
      DOUBLE PRECISION ROH2O(MAXCLY),ROO3(MAXCLY),RON2O(MAXCLY)
      DOUBLE PRECISION ROCH4(MAXCLY),ROCO(MAXCLY),ROCO2(MAXCLY)
      DOUBLE PRECISION ROO2(MAXCLY),DENS(MAXCLY)
      
      INTEGER*2 IATM 
      
      DOUBLE PRECISION PSURF  ! (I/O) Atmospheric pressure at the surface (mbar)

      DOUBLE PRECISION COEF	  ! Multiplicative factor of the pressure profile
      
C----------------------------------------------------------------------
      
      
      IF(IATM.LE.0)THEN  !User profile data

        COEF = 1.0
        IF(PSURF.GT.0.) COEF = PSURF / DONUSER(1,2)

        DO J=1,NLEVEL
           ALT(J)=DONUSER(J,1)
           P(J)=DONUSER(J,2)*COEF
           T(J)=DONUSER(J,3)
           ROH2O(J)=DONUSER(J,4)
           ROCO2(J)=DONUSER(J,5)
           ROO3(J)=DONUSER(J,6)
           RON2O(J)=DONUSER(J,7)
           ROCO(J)=DONUSER(J,8)
           ROCH4(J)=DONUSER(J,9)
           ROO2(J)=DONUSER(J,10)
           DENS(J)=DONUSER(J,11)
        ENDDO
      ENDIF 
      !Data of the predefined profiles for H2O CO2 O3 N2O CO CH4 O2
      IF(IATM.EQ.1)CALL TROPICA(P,T,ALT,ROH2O,ROCO2,ROO3,RON2O,ROCO,
     $                                ROCH4,ROO2,DENS)
      IF(IATM.EQ.2)CALL MIDLASU(P,T,ALT,ROH2O,ROCO2,ROO3,RON2O,ROCO,
     $                                ROCH4,ROO2,DENS)
      IF(IATM.EQ.3)CALL MIDLAWI(P,T,ALT,ROH2O,ROCO2,ROO3,RON2O,ROCO,
     $                                ROCH4,ROO2,DENS)
      IF(IATM.EQ.4)CALL SUBSUMM(P,T,ALT,ROH2O,ROCO2,ROO3,RON2O,ROCO,
     $                                ROCH4,ROO2,DENS)
      IF(IATM.EQ.5)CALL SUBWINT(P,T,ALT,ROH2O,ROCO2,ROO3,RON2O,ROCO,
     $                                ROCH4,ROO2,DENS)
      IF(IATM.EQ.6)CALL USTAD62(P,T,ALT,ROH2O,ROCO2,ROO3,RON2O,ROCO,
     $                                ROCH4,ROO2,DENS)

      IF(IATM.GT.0)THEN

        COEF = 1.0
        IF(PSURF.GT.0.) COEF = PSURF / P(1)

        DO J=1,NLEVEL
           DONUSER(J,1)=ALT(J)
           DONUSER(J,2)=P(J)*COEF
           DONUSER(J,3)=T(J)
           DONUSER(J,4)=ROH2O(J)
           DONUSER(J,5)=ROCO2(J)
           DONUSER(J,6)=ROO3(J)
           DONUSER(J,7)=RON2O(J)
           DONUSER(J,8)=ROCO(J)
           DONUSER(J,9)=ROCH4(J)
           DONUSER(J,10)=ROO2(J)
           DONUSER(J,11)=DENS(J)
        ENDDO
      ENDIF

      DO IK=1,NLEVEL
         RO(2,IK)=ROCO2(IK)*1.0E-06*44.0/28.97
         RO(4,IK)=RON2O(IK)*1.0E-06*44.0/28.97
         RO(5,IK)=ROCO(IK)*1.0E-06*28.0/28.97
         RO(6,IK)=ROCH4(IK)*1.0E-06*16.0/28.97
         RO(7,IK)=ROO2(IK)*1.0E-06*32.0/28.97
      ENDDO
      DO IK=1,NLEVEL
         ROH2O(IK)=ROH2O(IK)*1.0E-06*18.0/28.97
         RO(1,IK)=ROH2O(IK)/(1+ROH2O(IK))
         RO(3,IK)=ROO3(IK)*1.0E-06*48.0/28.97
      ENDDO
          
C
      RETURN
      END




      SUBROUTINE TROPICA(P,T,ALT,ROH2O,ROCO2,ROO3,RON2O,ROCO,ROCH4,ROO2,
     $                  DENS)

      PARAMETER(NLEVEL=50)
      DOUBLE PRECISION  P(NLEVEL),T(NLEVEL),ALT(NLEVEL)
      DOUBLE PRECISION  ROH2O(NLEVEL),ROO3(NLEVEL)
      DOUBLE PRECISION  RON2O(NLEVEL),ROO2(NLEVEL)
      DOUBLE PRECISION  ROCH4(NLEVEL),ROCO(NLEVEL),ROCO2(NLEVEL)
      DOUBLE PRECISION  P1(NLEVEL),T1(NLEVEL),ALT1(NLEVEL)
      DOUBLE PRECISION  ROH2O1(NLEVEL),ROO31(NLEVEL),RON2O1(NLEVEL)
      DOUBLE PRECISION  ROCH41(NLEVEL),ROCO1(NLEVEL),ROCO21(NLEVEL)
      DOUBLE PRECISION  ROO21(NLEVEL),DENS(NLEVEL),DENS1(NLEVEL)


C----------------------------------------------------------------------
C     DATA ALT (KM)  /
      DATA ALT1/
     C       0.0,       1.0,       2.0,       3.0,       4.0,
     C       5.0,       6.0,       7.0,       8.0,       9.0,
     C      10.0,      11.0,      12.0,      13.0,      14.0,
     C      15.0,      16.0,      17.0,      18.0,      19.0,
     C      20.0,      21.0,      22.0,      23.0,      24.0,
     C      25.0,      27.5,      30.0,      32.5,      35.0,
     C      37.5,      40.0,      42.5,      45.0,      47.5,
     C      50.0,      55.0,      60.0,      65.0,      70.0,
     C      75.0,      80.0,      85.0,      90.0,      95.0,
     C     100.0,     105.0,     110.0,     115.0,     120.0/
      DATA P1/
     C 1.013E+03, 9.040E+02, 8.050E+02, 7.150E+02, 6.330E+02,
     C 5.590E+02, 4.920E+02, 4.320E+02, 3.780E+02, 3.290E+02,
     C 2.860E+02, 2.470E+02, 2.130E+02, 1.820E+02, 1.560E+02,
     C 1.320E+02, 1.110E+02, 9.370E+01, 7.890E+01, 6.660E+01,
     C 5.650E+01, 4.800E+01, 4.090E+01, 3.500E+01, 3.000E+01,
     C 2.570E+01, 1.763E+01, 1.220E+01, 8.520E+00, 6.000E+00,
     C 4.260E+00, 3.050E+00, 2.200E+00, 1.590E+00, 1.160E+00,
     C 8.540E-01, 4.560E-01, 2.390E-01, 1.210E-01, 5.800E-02,
     C 2.600E-02, 1.100E-02, 4.400E-03, 1.720E-03, 6.880E-04,
     C 2.890E-04, 1.300E-04, 6.470E-05, 3.600E-05, 2.250E-05/
C----------------------------------------------------------------------
      DATA T1/
     C    299.70,    293.70,    287.70,    283.70,    277.00,
     C    270.30,    263.60,    257.00,    250.30,    243.60,
     C    237.00,    230.10,    223.60,    217.00,    210.30,
     C    203.70,    197.00,    194.80,    198.80,    202.70,
     C    206.70,    210.70,    214.60,    217.00,    219.20,
     C    221.40,    227.00,    232.30,    237.70,    243.10,
     C    248.50,    254.00,    259.40,    264.80,    269.60,
     C    270.20,    263.40,    253.10,    236.00,    218.90,
     C    201.80,    184.80,    177.10,    177.00,    184.30,
     C    190.70,    212.00,    241.60,    299.70,    380.00/
C----------------------------------------------------------------------
      DATA ROH2O1/
     C 2.593E+04, 1.949E+04, 1.534E+04, 8.600E+03, 4.441E+03,
     C 3.346E+03, 2.101E+03, 1.289E+03, 7.637E+02, 4.098E+02,
     C 1.912E+02, 7.306E+01, 2.905E+01, 9.900E+00, 6.220E+00,
     C 4.000E+00, 3.000E+00, 2.900E+00, 2.750E+00, 2.600E+00,
     C 2.600E+00, 2.650E+00, 2.800E+00, 2.900E+00, 3.200E+00,
     C 3.250E+00, 3.600E+00, 4.000E+00, 4.300E+00, 4.600E+00,
     C 4.900E+00, 5.200E+00, 5.500E+00, 5.700E+00, 5.900E+00,
     C 6.000E+00, 6.000E+00, 6.000E+00, 5.400E+00, 4.500E+00,
     C 3.300E+00, 2.100E+00, 1.300E+00, 8.500E-01, 5.400E-01,
     C 4.000E-01, 3.400E-01, 2.800E-01, 2.400E-01, 2.000E-01/
C     DATA CO2       /
      DATA ROCO21/
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.280E+02, 3.200E+02, 3.100E+02, 2.700E+02,
     C 1.950E+02, 1.100E+02, 6.000E+01, 4.000E+01, 3.500E+01/
C     DATA OZONE     /
      DATA ROO31/
     C 2.869E-02, 3.150E-02, 3.342E-02, 3.504E-02, 3.561E-02,
     C 3.767E-02, 3.989E-02, 4.223E-02, 4.471E-02, 5.000E-02,
     C 5.595E-02, 6.613E-02, 7.815E-02, 9.289E-02, 1.050E-01,
     C 1.256E-01, 1.444E-01, 2.500E-01, 5.000E-01, 9.500E-01,
     C 1.400E+00, 1.800E+00, 2.400E+00, 3.400E+00, 4.300E+00,
     C 5.400E+00, 7.800E+00, 9.300E+00, 9.850E+00, 9.700E+00,
     C 8.800E+00, 7.500E+00, 5.900E+00, 4.500E+00, 3.450E+00,
     C 2.800E+00, 1.800E+00, 1.100E+00, 6.500E-01, 3.000E-01,
     C 1.800E-01, 3.300E-01, 5.000E-01, 5.200E-01, 5.000E-01,
     C 4.000E-01, 2.000E-01, 5.000E-02, 5.000E-03, 5.000E-04/
C     DATA  N2O      /
      DATA RON2O1/
     C 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01,
     C 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.195E-01,
     C 3.179E-01, 3.140E-01, 3.095E-01, 3.048E-01, 2.999E-01,
     C 2.944E-01, 2.877E-01, 2.783E-01, 2.671E-01, 2.527E-01,
     C 2.365E-01, 2.194E-01, 2.051E-01, 1.967E-01, 1.875E-01,
     C 1.756E-01, 1.588E-01, 1.416E-01, 1.165E-01, 9.275E-02,
     C 6.693E-02, 4.513E-02, 2.751E-02, 1.591E-02, 9.378E-03,
     C 4.752E-03, 3.000E-03, 2.065E-03, 1.507E-03, 1.149E-03,
     C 8.890E-04, 7.056E-04, 5.716E-04, 4.708E-04, 3.932E-04,
     C 3.323E-04, 2.837E-04, 2.443E-04, 2.120E-04, 1.851E-04/
C     DATA CO        /
      DATA ROCO1/
     C 1.500E-01, 1.450E-01, 1.399E-01, 1.349E-01, 1.312E-01,
     C 1.303E-01, 1.288E-01, 1.247E-01, 1.185E-01, 1.094E-01,
     C 9.962E-02, 8.964E-02, 7.814E-02, 6.374E-02, 5.025E-02,
     C 3.941E-02, 3.069E-02, 2.489E-02, 1.966E-02, 1.549E-02,
     C 1.331E-02, 1.232E-02, 1.232E-02, 1.307E-02, 1.400E-02,
     C 1.521E-02, 1.722E-02, 1.995E-02, 2.266E-02, 2.487E-02,
     C 2.738E-02, 3.098E-02, 3.510E-02, 3.987E-02, 4.482E-02,
     C 5.092E-02, 5.985E-02, 6.960E-02, 9.188E-02, 1.938E-01,
     C 5.688E-01, 1.549E+00, 3.849E+00, 6.590E+00, 1.044E+01,
     C 1.705E+01, 2.471E+01, 3.358E+01, 4.148E+01, 5.000E+01/
C     DATA  CH4      /
      DATA ROCH41/
     C 1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00,
     C 1.700E+00, 1.700E+00, 1.699E+00, 1.697E+00, 1.693E+00,
     C 1.685E+00, 1.675E+00, 1.662E+00, 1.645E+00, 1.626E+00,
     C 1.605E+00, 1.582E+00, 1.553E+00, 1.521E+00, 1.480E+00,
     C 1.424E+00, 1.355E+00, 1.272E+00, 1.191E+00, 1.118E+00,
     C 1.055E+00, 9.870E-01, 9.136E-01, 8.300E-01, 7.460E-01,
     C 6.618E-01, 5.638E-01, 4.614E-01, 3.631E-01, 2.773E-01,
     C 2.100E-01, 1.651E-01, 1.500E-01, 1.500E-01, 1.500E-01,
     C 1.500E-01, 1.500E-01, 1.500E-01, 1.400E-01, 1.300E-01,
     C 1.200E-01, 1.100E-01, 9.500E-02, 6.000E-02, 3.000E-02/
C     DATA O2        /
      DATA ROO21/
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.000E+05, 1.900E+05, 1.800E+05,
     C 1.600E+05, 1.400E+05, 1.200E+05, 9.400E+04, 7.250E+04/
C     DATA DENSITY   /
      DATA DENS1/
     C 2.450E+19, 2.231E+19, 2.028E+19, 1.827E+19, 1.656E+19,
     C 1.499E+19, 1.353E+19, 1.218E+19, 1.095E+19, 9.789E+18,
     C 8.747E+18, 7.780E+18, 6.904E+18, 6.079E+18, 5.377E+18,
     C 4.697E+18, 4.084E+18, 3.486E+18, 2.877E+18, 2.381E+18,
     C 1.981E+18, 1.651E+18, 1.381E+18, 1.169E+18, 9.920E+17,
     C 8.413E+17, 5.629E+17, 3.807E+17, 2.598E+17, 1.789E+17,
     C 1.243E+17, 8.703E+16, 6.147E+16, 4.352E+16, 3.119E+16,
     C 2.291E+16, 1.255E+16, 6.844E+15, 3.716E+15, 1.920E+15,
     C 9.338E+14, 4.314E+14, 1.801E+14, 7.043E+13, 2.706E+13,
     C 1.098E+13, 4.445E+12, 1.941E+12, 8.706E+11, 4.225E+11/
C----------------------------------------------------------------------

      DO J=1,NLEVEL
         P(J)=P1(J)
         T(J)=T1(J)
         ALT(J)=ALT1(J)
         ROH2O(J)=ROH2O1(J)
         ROCO2(J)=ROCO21(J)
         ROO3(J)=ROO31(J)
         RON2O(J)=RON2O1(J)
         ROCO(J)=ROCO1(J)
         ROCH4(J)=ROCH41(J)
         ROO2(J)=ROO21(J)
         DENS(J)=DENS1(J)
      ENDDO

      RETURN
      END

      SUBROUTINE MIDLASU(P,T,ALT,ROH2O,ROCO2,ROO3,RON2O,ROCO,ROCH4,ROO2,
     $                  DENS)

      PARAMETER(NLEVEL=50)
      DOUBLE PRECISION  P(NLEVEL),T(NLEVEL),ALT(NLEVEL)
      DOUBLE PRECISION  ROH2O(NLEVEL),ROO3(NLEVEL)
      DOUBLE PRECISION  RON2O(NLEVEL),ROO2(NLEVEL)
      DOUBLE PRECISION  ROCH4(NLEVEL),ROCO(NLEVEL),ROCO2(NLEVEL)
      DOUBLE PRECISION P1(NLEVEL),T1(NLEVEL),ALT1(NLEVEL)
      DOUBLE PRECISION  ROH2O1(NLEVEL),ROO31(NLEVEL),RON2O1(NLEVEL)
      DOUBLE PRECISION  ROCH41(NLEVEL),ROCO1(NLEVEL),ROCO21(NLEVEL)
      DOUBLE PRECISION  ROO21(NLEVEL),DENS(NLEVEL),DENS1(NLEVEL)

C     DATA ALT (KM)  /
      DATA ALT1/
     C       0.0,       1.0,       2.0,       3.0,       4.0,
     C       5.0,       6.0,       7.0,       8.0,       9.0,
     C      10.0,      11.0,      12.0,      13.0,      14.0,
     C      15.0,      16.0,      17.0,      18.0,      19.0,
     C      20.0,      21.0,      22.0,      23.0,      24.0,
     C      25.0,      27.5,      30.0,      32.5,      35.0,
     C      37.5,      40.0,      42.5,      45.0,      47.5,
     C      50.0,      55.0,      60.0,      65.0,      70.0,
     C      75.0,      80.0,      85.0,      90.0,      95.0,
     C     100.0,     105.0,     110.0,     115.0,     120.0/
C     DATA PRESSURE  /
      DATA P1/
     C 1.013E+03, 9.020E+02, 8.020E+02, 7.100E+02, 6.280E+02,
     C 5.540E+02, 4.870E+02, 4.260E+02, 3.720E+02, 3.240E+02,
     C 2.810E+02, 2.430E+02, 2.090E+02, 1.790E+02, 1.530E+02,
     C 1.300E+02, 1.110E+02, 9.500E+01, 8.120E+01, 6.950E+01,
     C 5.950E+01, 5.100E+01, 4.370E+01, 3.760E+01, 3.220E+01,
     C 2.770E+01, 1.907E+01, 1.320E+01, 9.300E+00, 6.520E+00,
     C 4.640E+00, 3.330E+00, 2.410E+00, 1.760E+00, 1.290E+00,
     C 9.510E-01, 5.150E-01, 2.720E-01, 1.390E-01, 6.700E-02,
     C 3.000E-02, 1.200E-02, 4.480E-03, 1.640E-03, 6.250E-04,
     C 2.580E-04, 1.170E-04, 6.110E-05, 3.560E-05, 2.270E-05/
C     DATA TEMPERATUR/
      DATA T1/
     C    294.20,    289.70,    285.20,    279.20,    273.20,
     C    267.20,    261.20,    254.70,    248.20,    241.70,
     C    235.30,    228.80,    222.30,    215.80,    215.70,
     C    215.70,    215.70,    215.70,    216.80,    217.90,
     C    219.20,    220.40,    221.60,    222.80,    223.90,
     C    225.10,    228.45,    233.70,    239.00,    245.20,
     C    251.30,    257.50,    263.70,    269.90,    275.20,
     C    275.70,    269.30,    257.10,    240.10,    218.10,
     C    196.10,    174.10,    165.10,    165.00,    178.30,
     C    190.50,    222.20,    262.40,    316.80,    380.00/
C     DATA  H2O      /
      DATA ROH2O1/
     C 1.876E+04, 1.378E+04, 9.680E+03, 5.984E+03, 3.813E+03,
     C 2.225E+03, 1.510E+03, 1.020E+03, 6.464E+02, 4.129E+02,
     C 2.472E+02, 9.556E+01, 2.944E+01, 8.000E+00, 5.000E+00,
     C 3.400E+00, 3.300E+00, 3.200E+00, 3.150E+00, 3.200E+00,
     C 3.300E+00, 3.450E+00, 3.600E+00, 3.850E+00, 4.000E+00,
     C 4.200E+00, 4.450E+00, 4.700E+00, 4.850E+00, 4.950E+00,
     C 5.000E+00, 5.100E+00, 5.300E+00, 5.450E+00, 5.500E+00,
     C 5.500E+00, 5.350E+00, 5.000E+00, 4.400E+00, 3.700E+00,
     C 2.950E+00, 2.100E+00, 1.330E+00, 8.500E-01, 5.400E-01,
     C 4.000E-01, 3.400E-01, 2.800E-01, 2.400E-01, 2.000E-01/
C     DATA CO2       /
      DATA ROCO21/
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.280E+02, 3.200E+02, 3.100E+02, 2.700E+02,
     C 1.950E+02, 1.100E+02, 6.000E+01, 4.000E+01, 3.500E+01/
C     DATA OZONE     /
      DATA ROO31/
     C 3.017E-02, 3.337E-02, 3.694E-02, 4.222E-02, 4.821E-02,
     C 5.512E-02, 6.408E-02, 7.764E-02, 9.126E-02, 1.111E-01,
     C 1.304E-01, 1.793E-01, 2.230E-01, 3.000E-01, 4.400E-01,
     C 5.000E-01, 6.000E-01, 7.000E-01, 1.000E+00, 1.500E+00,
     C 2.000E+00, 2.400E+00, 2.900E+00, 3.400E+00, 4.000E+00,
     C 4.800E+00, 6.000E+00, 7.000E+00, 8.100E+00, 8.900E+00,
     C 8.700E+00, 7.550E+00, 5.900E+00, 4.500E+00, 3.500E+00,
     C 2.800E+00, 1.800E+00, 1.300E+00, 8.000E-01, 4.000E-01,
     C 1.900E-01, 2.000E-01, 5.700E-01, 7.500E-01, 7.000E-01,
     C 4.000E-01, 2.000E-01, 5.000E-02, 5.000E-03, 5.000E-04/
C     DATA  N2O      /
      DATA RON2O1/
     C 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01,
     C 3.200E-01, 3.200E-01, 3.200E-01, 3.195E-01, 3.163E-01,
     C 3.096E-01, 2.989E-01, 2.936E-01, 2.860E-01, 2.800E-01,
     C 2.724E-01, 2.611E-01, 2.421E-01, 2.174E-01, 1.843E-01,
     C 1.607E-01, 1.323E-01, 1.146E-01, 1.035E-01, 9.622E-02,
     C 8.958E-02, 8.006E-02, 6.698E-02, 4.958E-02, 3.695E-02,
     C 2.519E-02, 1.736E-02, 1.158E-02, 7.665E-03, 5.321E-03,
     C 3.215E-03, 2.030E-03, 1.397E-03, 1.020E-03, 7.772E-04,
     C 6.257E-04, 5.166E-04, 4.352E-04, 3.727E-04, 3.237E-04,
     C 2.844E-04, 2.524E-04, 2.260E-04, 2.039E-04, 1.851E-04/
C     DATA CO        /
      DATA ROCO1/
     C 1.500E-01, 1.450E-01, 1.399E-01, 1.349E-01, 1.312E-01,
     C 1.303E-01, 1.288E-01, 1.247E-01, 1.185E-01, 1.094E-01,
     C 9.962E-02, 8.964E-02, 7.814E-02, 6.374E-02, 5.025E-02,
     C 3.941E-02, 3.069E-02, 2.489E-02, 1.966E-02, 1.549E-02,
     C 1.331E-02, 1.232E-02, 1.232E-02, 1.307E-02, 1.400E-02,
     C 1.521E-02, 1.722E-02, 1.995E-02, 2.266E-02, 2.487E-02,
     C 2.716E-02, 2.962E-02, 3.138E-02, 3.307E-02, 3.487E-02,
     C 3.645E-02, 3.923E-02, 4.673E-02, 6.404E-02, 1.177E-01,
     C 2.935E-01, 6.815E-01, 1.465E+00, 2.849E+00, 5.166E+00,
     C 1.008E+01, 1.865E+01, 2.863E+01, 3.890E+01, 5.000E+01/
C     DATA  CH4      /
      DATA ROCH41/
     C 1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00, 1.697E+00,
     C 1.687E+00, 1.672E+00, 1.649E+00, 1.629E+00, 1.615E+00,
     C 1.579E+00, 1.542E+00, 1.508E+00, 1.479E+00, 1.451E+00,
     C 1.422E+00, 1.390E+00, 1.356E+00, 1.323E+00, 1.281E+00,
     C 1.224E+00, 1.154E+00, 1.066E+00, 9.730E-01, 8.800E-01,
     C 7.888E-01, 7.046E-01, 6.315E-01, 5.592E-01, 5.008E-01,
     C 4.453E-01, 3.916E-01, 3.389E-01, 2.873E-01, 2.384E-01,
     C 1.944E-01, 1.574E-01, 1.500E-01, 1.500E-01, 1.500E-01,
     C 1.500E-01, 1.500E-01, 1.500E-01, 1.400E-01, 1.300E-01,
     C 1.200E-01, 1.100E-01, 9.500E-02, 6.000E-02, 3.000E-02/
C     DATA O2        /
      DATA ROO21/
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.000E+05, 1.900E+05, 1.800E+05,
     C 1.600E+05, 1.400E+05, 1.200E+05, 9.400E+04, 7.250E+04/
C     DATA DENSITY   /
      DATA DENS1/
     C 2.496E+19, 2.257E+19, 2.038E+19, 1.843E+19, 1.666E+19,
     C 1.503E+19, 1.351E+19, 1.212E+19, 1.086E+19, 9.716E+18,
     C 8.656E+18, 7.698E+18, 6.814E+18, 6.012E+18, 5.141E+18,
     C 4.368E+18, 3.730E+18, 3.192E+18, 2.715E+18, 2.312E+18,
     C 1.967E+18, 1.677E+18, 1.429E+18, 1.223E+18, 1.042E+18,
     C 8.919E+17, 6.050E+17, 4.094E+17, 2.820E+17, 1.927E+17,
     C 1.338E+17, 9.373E+16, 6.624E+16, 4.726E+16, 3.398E+16,
     C 2.500E+16, 1.386E+16, 7.668E+15, 4.196E+15, 2.227E+15,
     C 1.109E+15, 4.996E+14, 1.967E+14, 7.204E+13, 2.541E+13,
     C 9.816E+12, 3.816E+12, 1.688E+12, 8.145E+11, 4.330E+11/


      DO J=1,NLEVEL
         P(J)=P1(J)
         T(J)=T1(J)
         ALT(J)=ALT1(J)
         ROH2O(J)=ROH2O1(J)
         ROCO2(J)=ROCO21(J)
         ROO3(J)=ROO31(J)
         RON2O(J)=RON2O1(J)
         ROCO(J)=ROCO1(J)
         ROCH4(J)=ROCH41(J)
         ROO2(J)=ROO21(J)
         DENS(J)=DENS1(J)
      ENDDO

      RETURN
      END

      SUBROUTINE MIDLAWI(P,T,ALT,ROH2O,ROCO2,ROO3,RON2O,ROCO,ROCH4,ROO2,
     $                  DENS)

      PARAMETER(NLEVEL=50)
      DOUBLE PRECISION  P(NLEVEL),T(NLEVEL),ALT(NLEVEL)
      DOUBLE PRECISION  ROH2O(NLEVEL),ROO3(NLEVEL)
      DOUBLE PRECISION  RON2O(NLEVEL),ROO2(NLEVEL)
      DOUBLE PRECISION  ROCH4(NLEVEL),ROCO(NLEVEL),ROCO2(NLEVEL)
      DOUBLE PRECISION  P1(NLEVEL),T1(NLEVEL),ALT1(NLEVEL)
      DOUBLE PRECISION  ROH2O1(NLEVEL),ROO31(NLEVEL),RON2O1(NLEVEL)
      DOUBLE PRECISION  ROCH41(NLEVEL),ROCO1(NLEVEL),ROCO21(NLEVEL)
      DOUBLE PRECISION  ROO21(NLEVEL),DENS(NLEVEL),DENS1(NLEVEL)


C     DATA ALT (KM)  /
      DATA ALT1/
     C       0.0,       1.0,       2.0,       3.0,       4.0,
     C       5.0,       6.0,       7.0,       8.0,       9.0,
     C      10.0,      11.0,      12.0,      13.0,      14.0,
     C      15.0,      16.0,      17.0,      18.0,      19.0,
     C      20.0,      21.0,      22.0,      23.0,      24.0,
     C      25.0,      27.5,      30.0,      32.5,      35.0,
     C      37.5,      40.0,      42.5,      45.0,      47.5,
     C      50.0,      55.0,      60.0,      65.0,      70.0,
     C      75.0,      80.0,      85.0,      90.0,      95.0,
     C     100.0,     105.0,     110.0,     115.0,     120.0/
C     DATA PRESSURE  /
      DATA P1/
     C 1.018E+03, 8.973E+02, 7.897E+02, 6.938E+02, 6.081E+02,
     C 5.313E+02, 4.627E+02, 4.016E+02, 3.473E+02, 2.993E+02,
     C 2.568E+02, 2.199E+02, 1.882E+02, 1.611E+02, 1.378E+02,
     C 1.178E+02, 1.007E+02, 8.610E+01, 7.360E+01, 6.280E+01,
     C 5.370E+01, 4.580E+01, 3.910E+01, 3.340E+01, 2.860E+01,
     C 2.440E+01, 1.646E+01, 1.110E+01, 7.560E+00, 5.180E+00,
     C 3.600E+00, 2.530E+00, 1.800E+00, 1.290E+00, 9.400E-01,
     C 6.830E-01, 3.620E-01, 1.880E-01, 9.500E-02, 4.700E-02,
     C 2.220E-02, 1.030E-02, 4.560E-03, 1.980E-03, 8.770E-04,
     C 4.074E-04, 2.000E-04, 1.057E-04, 5.980E-05, 3.600E-05/
C     DATA TEMPERATUR/
      DATA T1/
     C    272.20,    268.70,    265.20,    261.70,    255.70,
     C    249.70,    243.70,    237.70,    231.70,    225.70,
     C    219.70,    219.20,    218.70,    218.20,    217.70,
     C    217.20,    216.70,    216.20,    215.70,    215.20,
     C    215.20,    215.20,    215.20,    215.20,    215.20,
     C    215.20,    215.50,    217.40,    220.40,    227.90,
     C    235.50,    243.20,    250.80,    258.50,    265.10,
     C    265.70,    260.60,    250.80,    240.90,    230.70,
     C    220.40,    210.10,    199.80,    199.50,    208.30,
     C    218.60,    237.10,    259.50,    293.00,    333.00/
C     DATA  H2O      /
      DATA ROH2O1/
     C 4.316E+03, 3.454E+03, 2.788E+03, 2.088E+03, 1.280E+03,
     C 8.241E+02, 5.103E+02, 2.321E+02, 1.077E+02, 5.566E+01,
     C 2.960E+01, 1.000E+01, 6.000E+00, 5.000E+00, 4.800E+00,
     C 4.700E+00, 4.600E+00, 4.500E+00, 4.500E+00, 4.500E+00,
     C 4.500E+00, 4.500E+00, 4.530E+00, 4.550E+00, 4.600E+00,
     C 4.650E+00, 4.700E+00, 4.750E+00, 4.800E+00, 4.850E+00,
     C 4.900E+00, 4.950E+00, 5.000E+00, 5.000E+00, 5.000E+00,
     C 4.950E+00, 4.850E+00, 4.500E+00, 4.000E+00, 3.300E+00,
     C 2.700E+00, 2.000E+00, 1.330E+00, 8.500E-01, 5.400E-01,
     C 4.000E-01, 3.400E-01, 2.800E-01, 2.400E-01, 2.000E-01/
C     DATA CO2       /
      DATA ROCO21/
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.280E+02, 3.200E+02, 3.100E+02, 2.700E+02,
     C 1.950E+02, 1.100E+02, 6.000E+01, 4.000E+01, 3.500E+01/
C     DATA OZONE     /
      DATA ROO31/
     C 2.778E-02, 2.800E-02, 2.849E-02, 3.200E-02, 3.567E-02,
     C 4.720E-02, 5.837E-02, 7.891E-02, 1.039E-01, 1.567E-01,
     C 2.370E-01, 3.624E-01, 5.232E-01, 7.036E-01, 8.000E-01,
     C 9.000E-01, 1.100E+00, 1.400E+00, 1.800E+00, 2.300E+00,
     C 2.900E+00, 3.500E+00, 3.900E+00, 4.300E+00, 4.700E+00,
     C 5.100E+00, 5.600E+00, 6.100E+00, 6.800E+00, 7.100E+00,
     C 7.200E+00, 6.900E+00, 5.900E+00, 4.600E+00, 3.700E+00,
     C 2.750E+00, 1.700E+00, 1.000E-00, 5.500E-01, 3.200E-01,
     C 2.500E-01, 2.300E-01, 5.500E-01, 8.000E-01, 8.000E-01,
     C 4.000E-01, 2.000E-01, 5.000E-02, 5.000E-03, 5.000E-04/
C     DATA  N2O      /
      DATA RON2O1/
     C 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01,
     C 3.200E-01, 3.200E-01, 3.200E-01, 3.195E-01, 3.163E-01,
     C 3.096E-01, 2.989E-01, 2.936E-01, 2.860E-01, 2.800E-01,
     C 2.724E-01, 2.611E-01, 2.421E-01, 2.174E-01, 1.843E-01,
     C 1.621E-01, 1.362E-01, 1.230E-01, 1.124E-01, 1.048E-01,
     C 9.661E-02, 8.693E-02, 7.524E-02, 6.126E-02, 5.116E-02,
     C 3.968E-02, 2.995E-02, 2.080E-02, 1.311E-02, 8.071E-03,
     C 4.164E-03, 2.629E-03, 1.809E-03, 1.321E-03, 1.007E-03,
     C 7.883E-04, 6.333E-04, 5.194E-04, 4.333E-04, 3.666E-04,
     C 3.140E-04, 2.717E-04, 2.373E-04, 2.089E-04, 1.851E-04/
C     DATA CO        /
      DATA ROCO1/
     C 1.500E-01, 1.450E-01, 1.399E-01, 1.349E-01, 1.312E-01,
     C 1.303E-01, 1.288E-01, 1.247E-01, 1.185E-01, 1.094E-01,
     C 9.962E-02, 8.964E-02, 7.814E-02, 6.374E-02, 5.025E-02,
     C 3.941E-02, 3.069E-02, 2.489E-02, 1.966E-02, 1.549E-02,
     C 1.331E-02, 1.232E-02, 1.232E-02, 1.307E-02, 1.400E-02,
     C 1.498E-02, 1.598E-02, 1.710E-02, 1.850E-02, 1.997E-02,
     C 2.147E-02, 2.331E-02, 2.622E-02, 3.057E-02, 3.803E-02,
     C 6.245E-02, 1.480E-01, 2.926E-01, 5.586E-01, 1.078E+00,
     C 1.897E+00, 2.960E+00, 4.526E+00, 6.862E+00, 1.054E+01,
     C 1.709E+01, 2.473E+01, 3.359E+01, 4.149E+01, 5.000E+01/
C     DATA  CH4      /
      DATA ROCH41/
     C 1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00, 1.697E+00,
     C 1.687E+00, 1.672E+00, 1.649E+00, 1.629E+00, 1.615E+00,
     C 1.579E+00, 1.542E+00, 1.508E+00, 1.479E+00, 1.451E+00,
     C 1.422E+00, 1.390E+00, 1.356E+00, 1.323E+00, 1.281E+00,
     C 1.224E+00, 1.154E+00, 1.066E+00, 9.730E-01, 8.800E-01,
     C 7.931E-01, 7.130E-01, 6.438E-01, 5.746E-01, 5.050E-01,
     C 4.481E-01, 3.931E-01, 3.395E-01, 2.876E-01, 2.386E-01,
     C 1.944E-01, 1.574E-01, 1.500E-01, 1.500E-01, 1.500E-01,
     C 1.500E-01, 1.500E-01, 1.500E-01, 1.400E-01, 1.300E-01,
     C 1.200E-01, 1.100E-01, 9.500E-02, 6.000E-02, 3.000E-02/
C     DATA O2        /
      DATA ROO21/
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.000E+05, 1.900E+05, 1.800E+05,
     C 1.600E+05, 1.400E+05, 1.200E+05, 9.400E+04, 7.250E+04/
C     DATA DENS1/
      DATA DENS1/
     C 2.711E+19, 2.420E+19, 2.158E+19, 1.922E+19, 1.724E+19,
     C 1.542E+19, 1.376E+19, 1.225E+19, 1.086E+19, 9.612E+18,
     C 8.472E+18, 7.271E+18, 6.237E+18, 5.351E+18, 4.588E+18,
     C 3.931E+18, 3.368E+18, 2.886E+18, 2.473E+18, 2.115E+18,
     C 1.809E+18, 1.543E+18, 1.317E+18, 1.125E+18, 9.633E+17,
     C 8.218E+17, 5.536E+17, 3.701E+17, 2.486E+17, 1.647E+17,
     C 1.108E+17, 7.540E+16, 5.202E+16, 3.617E+16, 2.570E+16,
     C 1.863E+16, 1.007E+16, 5.433E+15, 2.858E+15, 1.477E+15,
     C 7.301E+14, 3.553E+14, 1.654E+14, 7.194E+13, 3.052E+13,
     C 1.351E+13, 6.114E+12, 2.952E+12, 1.479E+12, 7.836E+11/

      DO J=1,NLEVEL
         P(J)=P1(J)
         T(J)=T1(J)
         ALT(J)=ALT1(J)
         ROH2O(J)=ROH2O1(J)
         ROCO2(J)=ROCO21(J)
         ROO3(J)=ROO31(J)
         RON2O(J)=RON2O1(J)
         ROCO(J)=ROCO1(J)
         ROCH4(J)=ROCH41(J)
         ROO2(J)=ROO21(J)
         DENS(J)=DENS1(J)
      ENDDO

      RETURN
      END

      SUBROUTINE SUBSUMM(P,T,ALT,ROH2O,ROCO2,ROO3,RON2O,ROCO,ROCH4,ROO2,
     $                  DENS)

      PARAMETER(NLEVEL=50)
      DOUBLE PRECISION  P(NLEVEL),T(NLEVEL),ALT(NLEVEL)
      DOUBLE PRECISION  ROH2O(NLEVEL),ROO3(NLEVEL)
      DOUBLE PRECISION  RON2O(NLEVEL),ROO2(NLEVEL)
      DOUBLE PRECISION  ROCH4(NLEVEL),ROCO(NLEVEL),ROCO2(NLEVEL)
      DOUBLE PRECISION  P1(NLEVEL),T1(NLEVEL),ALT1(NLEVEL)
      DOUBLE PRECISION  ROH2O1(NLEVEL),ROO31(NLEVEL),RON2O1(NLEVEL)
      DOUBLE PRECISION  ROCH41(NLEVEL),ROCO1(NLEVEL),ROCO21(NLEVEL)
      DOUBLE PRECISION  ROO21(NLEVEL),DENS(NLEVEL),DENS1(NLEVEL)

C     DATA ALT (KM)  /
      DATA ALT1/
     C       0.0,       1.0,       2.0,       3.0,       4.0,
     C       5.0,       6.0,       7.0,       8.0,       9.0,
     C      10.0,      11.0,      12.0,      13.0,      14.0,
     C      15.0,      16.0,      17.0,      18.0,      19.0,
     C      20.0,      21.0,      22.0,      23.0,      24.0,
     C      25.0,      27.5,      30.0,      32.5,      35.0,
     C      37.5,      40.0,      42.5,      45.0,      47.5,
     C      50.0,      55.0,      60.0,      65.0,      70.0,
     C      75.0,      80.0,      85.0,      90.0,      95.0,
     C     100.0,     105.0,     110.0,     115.0,     120.0/
C     DATA PRESSURE  /
      DATA P1/
     C 1.010E+03, 8.960E+02, 7.929E+02, 7.000E+02, 6.160E+02,
     C 5.410E+02, 4.740E+02, 4.130E+02, 3.590E+02, 3.108E+02,
     C 2.677E+02, 2.300E+02, 1.977E+02, 1.700E+02, 1.460E+02,
     C 1.260E+02, 1.080E+02, 9.280E+01, 7.980E+01, 6.860E+01,
     C 5.900E+01, 5.070E+01, 4.360E+01, 3.750E+01, 3.228E+01,
     C 2.780E+01, 1.923E+01, 1.340E+01, 9.400E+00, 6.610E+00,
     C 4.720E+00, 3.400E+00, 2.480E+00, 1.820E+00, 1.340E+00,
     C 9.870E-01, 5.370E-01, 2.880E-01, 1.470E-01, 7.100E-02,
     C 3.200E-02, 1.250E-02, 4.510E-03, 1.610E-03, 6.060E-04,
     C 2.480E-04, 1.130E-04, 6.000E-05, 3.540E-05, 2.260E-05/
C     DATA TEMPERATUR/
      DATA T1/
     C    287.20,    281.70,    276.30,    270.90,    265.50,
     C    260.10,    253.10,    246.10,    239.20,    232.20,
     C    225.20,    225.20,    225.20,    225.20,    225.20,
     C    225.20,    225.20,    225.20,    225.20,    225.20,
     C    225.20,    225.20,    225.20,    225.20,    226.60,
     C    228.10,    231.00,    235.10,    240.00,    247.20,
     C    254.60,    262.10,    269.50,    273.60,    276.20,
     C    277.20,    274.00,    262.70,    239.70,    216.60,
     C    193.60,    170.60,    161.70,    161.60,    176.80,
     C    190.40,    226.00,    270.10,    322.70,    380.00/
C     DATA  H2O      /
      DATA ROH2O1/
     C 1.194E+04, 8.701E+03, 6.750E+03, 4.820E+03, 3.380E+03,
     C 2.218E+03, 1.330E+03, 7.971E+02, 3.996E+02, 1.300E+02,
     C 4.240E+01, 1.330E+01, 6.000E+00, 4.450E+00, 4.000E+00,
     C 4.000E+00, 4.000E+00, 4.050E+00, 4.300E+00, 4.500E+00,
     C 4.600E+00, 4.700E+00, 4.800E+00, 4.830E+00, 4.850E+00,
     C 4.900E+00, 4.950E+00, 5.000E+00, 5.000E+00, 5.000E+00,
     C 5.000E+00, 5.000E+00, 5.000E+00, 5.000E+00, 5.000E+00,
     C 4.950E+00, 4.850E+00, 4.500E+00, 4.000E+00, 3.300E+00,
     C 2.700E+00, 2.000E+00, 1.330E+00, 8.500E-01, 5.400E-01,
     C 4.000E-01, 3.400E-01, 2.800E-01, 2.400E-01, 2.000E-01/
C     DATA CO2       /
      DATA ROCO21/
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.280E+02, 3.200E+02, 3.100E+02, 2.700E+02,
     C 1.950E+02, 1.100E+02, 6.000E+01, 4.000E+01, 3.500E+01/
C     DATA OZONE     /
      DATA ROO31/
     C 2.412E-02, 2.940E-02, 3.379E-02, 3.887E-02, 4.478E-02,
     C 5.328E-02, 6.564E-02, 7.738E-02, 9.114E-02, 1.420E-01,
     C 1.890E-01, 3.050E-01, 4.100E-01, 5.000E-01, 6.000E-01,
     C 7.000E-01, 8.500E-01, 1.000E+00, 1.300E+00, 1.700E+00,
     C 2.100E+00, 2.700E+00, 3.300E+00, 3.700E+00, 4.200E+00,
     C 4.500E+00, 5.300E+00, 5.700E+00, 6.900E+00, 7.700E+00,
     C 7.800E+00, 7.000E+00, 5.400E+00, 4.200E+00, 3.200E+00,
     C 2.500E+00, 1.700E+00, 1.200E+00, 8.000E-01, 4.000E-01,
     C 2.000E-01, 1.800E-01, 6.500E-01, 9.000E-01, 8.000E-01,
     C 4.000E-01, 2.000E-01, 5.000E-02, 5.000E-03, 5.000E-04/
C     DATA  N2O      /
      DATA RON2O1/
     C 3.100E-01, 3.100E-01, 3.100E-01, 3.100E-01, 3.079E-01,
     C 3.024E-01, 2.906E-01, 2.822E-01, 2.759E-01, 2.703E-01,
     C 2.651E-01, 2.600E-01, 2.549E-01, 2.494E-01, 2.433E-01,
     C 2.355E-01, 2.282E-01, 2.179E-01, 2.035E-01, 1.817E-01,
     C 1.567E-01, 1.350E-01, 1.218E-01, 1.102E-01, 9.893E-02,
     C 8.775E-02, 7.327E-02, 5.941E-02, 4.154E-02, 3.032E-02,
     C 1.949E-02, 1.274E-02, 9.001E-03, 6.286E-03, 4.558E-03,
     C 2.795E-03, 1.765E-03, 1.214E-03, 8.866E-04, 6.756E-04,
     C 5.538E-04, 4.649E-04, 3.979E-04, 3.459E-04, 3.047E-04,
     C 2.713E-04, 2.439E-04, 2.210E-04, 2.017E-04, 1.851E-04/
C     DATA CO        /
      DATA ROCO1/
     C 1.500E-01, 1.450E-01, 1.399E-01, 1.349E-01, 1.312E-01,
     C 1.303E-01, 1.288E-01, 1.247E-01, 1.185E-01, 1.094E-01,
     C 9.962E-02, 8.964E-02, 7.814E-02, 6.374E-02, 5.025E-02,
     C 3.941E-02, 3.069E-02, 2.489E-02, 1.966E-02, 1.549E-02,
     C 1.331E-02, 1.232E-02, 1.232E-02, 1.307E-02, 1.400E-02,
     C 1.510E-02, 1.649E-02, 1.808E-02, 1.997E-02, 2.183E-02,
     C 2.343E-02, 2.496E-02, 2.647E-02, 2.809E-02, 2.999E-02,
     C 3.220E-02, 3.650E-02, 4.589E-02, 6.375E-02, 1.176E-01,
     C 3.033E-01, 7.894E-01, 1.823E+00, 3.402E+00, 5.916E+00,
     C 1.043E+01, 1.881E+01, 2.869E+01, 3.892E+01, 5.000E+01/
C     DATA  CH4      /
      DATA ROCH41/
     C 1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00, 1.697E+00,
     C 1.687E+00, 1.672E+00, 1.649E+00, 1.629E+00, 1.615E+00,
     C 1.579E+00, 1.542E+00, 1.506E+00, 1.471E+00, 1.434E+00,
     C 1.389E+00, 1.342E+00, 1.290E+00, 1.230E+00, 1.157E+00,
     C 1.072E+00, 9.903E-01, 9.170E-01, 8.574E-01, 8.013E-01,
     C 7.477E-01, 6.956E-01, 6.442E-01, 5.888E-01, 5.240E-01,
     C 4.506E-01, 3.708E-01, 2.992E-01, 2.445E-01, 2.000E-01,
     C 1.660E-01, 1.500E-01, 1.500E-01, 1.500E-01, 1.500E-01,
     C 1.500E-01, 1.500E-01, 1.500E-01, 1.400E-01, 1.300E-01,
     C 1.200E-01, 1.100E-01, 9.500E-02, 6.000E-02, 3.000E-02/
C     DATA O2        /
      DATA ROO21/
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.000E+05, 1.900E+05, 1.800E+05,
     C 1.600E+05, 1.400E+05, 1.200E+05, 9.400E+04, 7.250E+04/
C     DATA DENSITY   /
      DATA DENS1/
     C 2.549E+19, 2.305E+19, 2.080E+19, 1.873E+19, 1.682E+19,
     C 1.508E+19, 1.357E+19, 1.216E+19, 1.088E+19, 9.701E+18,
     C 8.616E+18, 7.402E+18, 6.363E+18, 5.471E+18, 4.699E+18,
     C 4.055E+18, 3.476E+18, 2.987E+18, 2.568E+18, 2.208E+18,
     C 1.899E+18, 1.632E+18, 1.403E+18, 1.207E+18, 1.033E+18,
     C 8.834E+17, 6.034E+17, 4.131E+17, 2.839E+17, 1.938E+17,
     C 1.344E+17, 9.402E+16, 6.670E+16, 4.821E+16, 3.516E+16,
     C 2.581E+16, 1.421E+16, 7.946E+15, 4.445E+15, 2.376E+15,
     C 1.198E+15, 5.311E+14, 2.022E+14, 7.221E+13, 2.484E+13,
     C 9.441E+12, 3.624E+12, 1.610E+12, 7.951E+11, 4.311E+11/

      DO J=1,NLEVEL
         P(J)=P1(J)
         T(J)=T1(J)
         ALT(J)=ALT1(J)
         ROH2O(J)=ROH2O1(J)
         ROCO2(J)=ROCO21(J)
         ROO3(J)=ROO31(J)
         RON2O(J)=RON2O1(J)
         ROCO(J)=ROCO1(J)
         ROCH4(J)=ROCH41(J)
         ROO2(J)=ROO21(J)
         DENS(J)=DENS1(J)
      ENDDO

      RETURN
      END

      SUBROUTINE SUBWINT(P,T,ALT,ROH2O,ROCO2,ROO3,RON2O,ROCO,ROCH4,ROO2,
     $                  DENS)

      PARAMETER(NLEVEL=50)
      DOUBLE PRECISION  P(NLEVEL),T(NLEVEL),ALT(NLEVEL)
      DOUBLE PRECISION  ROH2O(NLEVEL),ROO3(NLEVEL)
      DOUBLE PRECISION  RON2O(NLEVEL),ROO2(NLEVEL)
      DOUBLE PRECISION  ROCH4(NLEVEL),ROCO(NLEVEL),ROCO2(NLEVEL)
      DOUBLE PRECISION  P1(NLEVEL),T1(NLEVEL),ALT1(NLEVEL)
      DOUBLE PRECISION  ROH2O1(NLEVEL),ROO31(NLEVEL),RON2O1(NLEVEL)
      DOUBLE PRECISION  ROCH41(NLEVEL),ROCO1(NLEVEL),ROCO21(NLEVEL)
      DOUBLE PRECISION  ROO21(NLEVEL),DENS(NLEVEL),DENS1(NLEVEL)


C     DATA ALT (KM)  /
      DATA ALT1/
     C       0.0,       1.0,       2.0,       3.0,       4.0,
     C       5.0,       6.0,       7.0,       8.0,       9.0,
     C      10.0,      11.0,      12.0,      13.0,      14.0,
     C      15.0,      16.0,      17.0,      18.0,      19.0,
     C      20.0,      21.0,      22.0,      23.0,      24.0,
     C      25.0,      27.5,      30.0,      32.5,      35.0,
     C      37.5,      40.0,      42.5,      45.0,      47.5,
     C      50.0,      55.0,      60.0,      65.0,      70.0,
     C      75.0,      80.0,      85.0,      90.0,      95.0,
     C     100.0,     105.0,     110.0,     115.0,     120.0/
C     DATA PRESSURE  /
      DATA P1/
     C 1.013E+03, 8.878E+02, 7.775E+02, 6.798E+02, 5.932E+02,
     C 5.158E+02, 4.467E+02, 3.853E+02, 3.308E+02, 2.829E+02,
     C 2.418E+02, 2.067E+02, 1.766E+02, 1.510E+02, 1.291E+02,
     C 1.103E+02, 9.431E+01, 8.058E+01, 6.882E+01, 5.875E+01,
     C 5.014E+01, 4.277E+01, 3.647E+01, 3.109E+01, 2.649E+01,
     C 2.256E+01, 1.513E+01, 1.020E+01, 6.910E+00, 4.701E+00,
     C 3.230E+00, 2.243E+00, 1.570E+00, 1.113E+00, 7.900E-01,
     C 5.719E-01, 2.990E-01, 1.550E-01, 7.900E-02, 4.000E-02,
     C 2.000E-02, 9.660E-03, 4.500E-03, 2.022E-03, 9.070E-04,
     C 4.230E-04, 2.070E-04, 1.080E-04, 6.000E-05, 3.590E-05/
C     DATA TEMPERATUR/
      DATA T1/
     C    257.20,    259.10,    255.90,    252.70,    247.70,
     C    240.90,    234.10,    227.30,    220.60,    217.20,
     C    217.20,    217.20,    217.20,    217.20,    217.20,
     C    217.20,    216.60,    216.00,    215.40,    214.80,
     C    214.20,    213.60,    213.00,    212.40,    211.80,
     C    211.20,    213.60,    216.00,    218.50,    222.30,
     C    228.50,    234.70,    240.80,    247.00,    253.20,
     C    259.30,    259.10,    250.90,    248.40,    245.40,
     C    234.70,    223.90,    213.10,    202.30,    211.00,
     C    218.50,    234.00,    252.60,    288.50,    333.00/
C     DATA  H2O      /
      DATA ROH2O1/
     C 1.405E+03, 1.615E+03, 1.427E+03, 1.166E+03, 7.898E+02,
     C 4.309E+02, 2.369E+02, 1.470E+02, 3.384E+01, 2.976E+01,
     C 2.000E+01, 1.000E+01, 6.000E+00, 4.450E+00, 4.500E+00,
     C 4.550E+00, 4.600E+00, 4.650E+00, 4.700E+00, 4.750E+00,
     C 4.800E+00, 4.850E+00, 4.900E+00, 4.950E+00, 5.000E+00,
     C 5.000E+00, 5.000E+00, 5.000E+00, 5.000E+00, 5.000E+00,
     C 5.000E+00, 5.000E+00, 5.000E+00, 5.000E+00, 5.000E+00,
     C 4.950E+00, 4.850E+00, 4.500E+00, 4.000E+00, 3.300E+00,
     C 2.700E+00, 2.000E+00, 1.330E+00, 8.500E-01, 5.400E-01,
     C 4.000E-01, 3.400E-01, 2.800E-01, 2.400E-01, 2.000E-01/
C     DATA CO2       /
      DATA ROCO21/
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.280E+02, 3.200E+02, 3.100E+02, 2.700E+02,
     C 1.950E+02, 1.100E+02, 6.000E+01, 4.000E+01, 3.500E+01/
C     DATA OZONE     /
      DATA ROO31/
     C 1.802E-02, 2.072E-02, 2.336E-02, 2.767E-02, 3.253E-02,
     C 3.801E-02, 4.446E-02, 7.252E-02, 1.040E-01, 2.100E-01,
     C 3.000E-01, 3.500E-01, 4.000E-01, 6.500E-01, 9.000E-01,
     C 1.200E+00, 1.500E+00, 1.900E+00, 2.450E+00, 3.100E+00,
     C 3.700E+00, 4.000E+00, 4.200E+00, 4.500E+00, 4.600E+00,
     C 4.700E+00, 4.900E+00, 5.400E+00, 5.900E+00, 6.200E+00,
     C 6.250E+00, 5.900E+00, 5.100E+00, 4.100E+00, 3.000E+00,
     C 2.600E+00, 1.600E+00, 9.500E-01, 6.500E-01, 5.000E-01,
     C 3.300E-01, 1.300E-01, 7.500E-01, 8.000E-01, 8.000E-01,
     C 4.000E-01, 2.000E-01, 5.000E-02, 5.000E-03, 5.000E-04/
C     DATA  N2O      /
      DATA RON2O1/
     C 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01,
     C 3.200E-01, 3.200E-01, 3.200E-01, 3.195E-01, 3.163E-01,
     C 3.096E-01, 2.989E-01, 2.936E-01, 2.860E-01, 2.800E-01,
     C 2.724E-01, 2.611E-01, 2.421E-01, 2.174E-01, 1.843E-01,
     C 1.621E-01, 1.362E-01, 1.230E-01, 1.122E-01, 1.043E-01,
     C 9.570E-02, 8.598E-02, 7.314E-02, 5.710E-02, 4.670E-02,
     C 3.439E-02, 2.471E-02, 1.631E-02, 1.066E-02, 7.064E-03,
     C 3.972E-03, 2.508E-03, 1.726E-03, 1.260E-03, 9.602E-04,
     C 7.554E-04, 6.097E-04, 5.024E-04, 4.210E-04, 3.579E-04,
     C 3.080E-04, 2.678E-04, 2.350E-04, 2.079E-04, 1.851E-04/
C     DATA CO        /
C     DATA CO        /
      DATA ROCO1/
     C 1.500E-01, 1.450E-01, 1.399E-01, 1.349E-01, 1.312E-01,
     C 1.303E-01, 1.288E-01, 1.247E-01, 1.185E-01, 1.094E-01,
     C 9.962E-02, 8.964E-02, 7.814E-02, 6.374E-02, 5.025E-02,
     C 3.941E-02, 3.069E-02, 2.489E-02, 1.966E-02, 1.549E-02,
     C 1.331E-02, 1.232E-02, 1.232E-02, 1.307E-02, 1.400E-02,
     C 1.521E-02, 1.722E-02, 2.037E-02, 2.486E-02, 3.168E-02,
     C 4.429E-02, 6.472E-02, 1.041E-01, 1.507E-01, 2.163E-01,
     C 3.141E-01, 4.842E-01, 7.147E-01, 1.067E+00, 1.516E+00,
     C 2.166E+00, 3.060E+00, 4.564E+00, 6.877E+00, 1.055E+01,
     C 1.710E+01, 2.473E+01, 3.359E+01, 4.149E+01, 5.000E+01/
C     DATA  CH4      /
      DATA ROCH41/
     C 1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00, 1.697E+00,
     C 1.687E+00, 1.672E+00, 1.649E+00, 1.629E+00, 1.615E+00,
     C 1.579E+00, 1.542E+00, 1.506E+00, 1.471E+00, 1.434E+00,
     C 1.389E+00, 1.342E+00, 1.290E+00, 1.230E+00, 1.161E+00,
     C 1.084E+00, 1.014E+00, 9.561E-01, 9.009E-01, 8.479E-01,
     C 7.961E-01, 7.449E-01, 6.941E-01, 6.434E-01, 5.883E-01,
     C 5.238E-01, 4.505E-01, 3.708E-01, 3.004E-01, 2.453E-01,
     C 1.980E-01, 1.590E-01, 1.500E-01, 1.500E-01, 1.500E-01,
     C 1.500E-01, 1.500E-01, 1.500E-01, 1.400E-01, 1.300E-01,
     C 1.200E-01, 1.100E-01, 9.500E-02, 6.000E-02, 3.000E-02/
C     DATA O2        /
      DATA ROO21/
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.000E+05, 1.900E+05, 1.800E+05,
     C 1.600E+05, 1.400E+05, 1.200E+05, 9.400E+04, 7.250E+04/
C     DATA DENSITY   /
      DATA DENS1/
     C 2.855E+19, 2.484E+19, 2.202E+19, 1.950E+19, 1.736E+19,
     C 1.552E+19, 1.383E+19, 1.229E+19, 1.087E+19, 9.440E+18,
     C 8.069E+18, 6.898E+18, 5.893E+18, 5.039E+18, 4.308E+18,
     C 3.681E+18, 3.156E+18, 2.704E+18, 2.316E+18, 1.982E+18,
     C 1.697E+18, 1.451E+18, 1.241E+18, 1.061E+18, 9.065E+17,
     C 7.742E+17, 5.134E+17, 3.423E+17, 2.292E+17, 1.533E+17,
     C 1.025E+17, 6.927E+16, 4.726E+16, 3.266E+16, 2.261E+16,
     C 1.599E+16, 8.364E+15, 4.478E+15, 2.305E+15, 1.181E+15,
     C 6.176E+14, 3.127E+14, 1.531E+14, 7.244E+13, 3.116E+13,
     C 1.403E+13, 6.412E+12, 3.099E+12, 1.507E+12, 7.814E+11/
 
      DO J=1,NLEVEL
         P(J)=P1(J)
         T(J)=T1(J)
         ALT(J)=ALT1(J)
         ROH2O(J)=ROH2O1(J)
         ROCO2(J)=ROCO21(J)
         ROO3(J)=ROO31(J)
         RON2O(J)=RON2O1(J)
         ROCO(J)=ROCO1(J)
         ROCH4(J)=ROCH41(J)
         ROO2(J)=ROO21(J)
         DENS(J)=DENS1(J)
      ENDDO

      RETURN
      END

      SUBROUTINE USTAD62(P,T,ALT,ROH2O,ROCO2,ROO3,RON2O,ROCO,ROCH4,ROO2,
     $                  DENS)

      PARAMETER(NLEVEL=50)
      DOUBLE PRECISION  P(NLEVEL),T(NLEVEL),ALT(NLEVEL)
      DOUBLE PRECISION  ROH2O(NLEVEL),ROO3(NLEVEL)
      DOUBLE PRECISION  RON2O(NLEVEL),ROO2(NLEVEL)
      DOUBLE PRECISION  ROCH4(NLEVEL),ROCO(NLEVEL),ROCO2(NLEVEL)
      DOUBLE PRECISION  P1(NLEVEL),T1(NLEVEL),ALT1(NLEVEL)
      DOUBLE PRECISION  ROH2O1(NLEVEL),ROO31(NLEVEL),RON2O1(NLEVEL)
      DOUBLE PRECISION  ROCH41(NLEVEL),ROCO1(NLEVEL),ROCO21(NLEVEL)
      DOUBLE PRECISION  ROO21(NLEVEL),DENS(NLEVEL),DENS1(NLEVEL)


C     DATA ALT (KM)  /
      DATA ALT1/
     C       0.0,       1.0,       2.0,       3.0,       4.0,
     C       5.0,       6.0,       7.0,       8.0,       9.0,
     C      10.0,      11.0,      12.0,      13.0,      14.0,
     C      15.0,      16.0,      17.0,      18.0,      19.0,
     C      20.0,      21.0,      22.0,      23.0,      24.0,
     C      25.0,      27.5,      30.0,      32.5,      35.0,
     C      37.5,      40.0,      42.5,      45.0,      47.5,
     C      50.0,      55.0,      60.0,      65.0,      70.0,
     C      75.0,      80.0,      85.0,      90.0,      95.0,
     C     100.0,     105.0,     110.0,     115.0,     120.0/
C     DATA PRESSURE  /
      DATA P1/
     C 1.013E+03, 8.988E+02, 7.950E+02, 7.012E+02, 6.166E+02,
     C 5.405E+02, 4.722E+02, 4.111E+02, 3.565E+02, 3.080E+02,
     C 2.650E+02, 2.270E+02, 1.940E+02, 1.658E+02, 1.417E+02,
     C 1.211E+02, 1.035E+02, 8.850E+01, 7.565E+01, 6.467E+01,
     C 5.529E+01, 4.729E+01, 4.047E+01, 3.467E+01, 2.972E+01,
     C 2.549E+01, 1.743E+01, 1.197E+01, 8.010E+00, 5.746E+00,
     C 4.150E+00, 2.871E+00, 2.060E+00, 1.491E+00, 1.090E+00,
     C 7.978E-01, 4.250E-01, 2.190E-01, 1.090E-01, 5.220E-02,
     C 2.400E-02, 1.050E-02, 4.460E-03, 1.840E-03, 7.600E-04,
     C 3.200E-04, 1.450E-04, 7.100E-05, 4.010E-05, 2.540E-05/
C     DATA TEMPERATUR/
      DATA T1/
     C    288.20,    281.70,    275.20,    268.70,    262.20,
     C    255.70,    249.20,    242.70,    236.20,    229.70,
     C    223.30,    216.80,    216.70,    216.70,    216.70,
     C    216.70,    216.70,    216.70,    216.70,    216.70,
     C    216.70,    217.60,    218.60,    219.60,    220.60,
     C    221.60,    224.00,    226.50,    230.00,    236.50,
     C    242.90,    250.40,    257.30,    264.20,    270.60,
     C    270.70,    260.80,    247.00,    233.30,    219.60,
     C    208.40,    198.60,    188.90,    186.90,    188.40,
     C    195.10,    208.80,    240.00,    300.00,    360.00/
C     DATA  H2O      /
      DATA ROH2O1/
     C 7.745E+03, 6.071E+03, 4.631E+03, 3.182E+03, 2.158E+03,
     C 1.397E+03, 9.254E+02, 5.720E+02, 3.667E+02, 1.583E+02,
     C 6.996E+01, 3.613E+01, 1.906E+01, 1.085E+01, 5.927E+00,
     C 5.000E+00, 3.950E+00, 3.850E+00, 3.825E+00, 3.850E+00,
     C 3.900E+00, 3.975E+00, 4.065E+00, 4.200E+00, 4.300E+00,
     C 4.425E+00, 4.575E+00, 4.725E+00, 4.825E+00, 4.900E+00,
     C 4.950E+00, 5.025E+00, 5.150E+00, 5.225E+00, 5.250E+00,
     C 5.225E+00, 5.100E+00, 4.750E+00, 4.200E+00, 3.500E+00,
     C 2.825E+00, 2.050E+00, 1.330E+00, 8.500E-01, 5.400E-01,
     C 4.000E-01, 3.400E-01, 2.800E-01, 2.400E-01, 2.000E-01/
C     DATA CO2       /
      DATA ROCO21/
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02, 3.300E+02,
     C 3.300E+02, 3.280E+02, 3.200E+02, 3.100E+02, 2.700E+02,
     C 1.950E+02, 1.100E+02, 6.000E+01, 4.000E+01, 3.500E+01/
C     DATA OZONE     /
      DATA ROO31/
     C 2.660E-02, 2.931E-02, 3.237E-02, 3.318E-02, 3.387E-02,
     C 3.768E-02, 4.112E-02, 5.009E-02, 5.966E-02, 9.168E-02,
     C 1.313E-01, 2.149E-01, 3.095E-01, 3.846E-01, 5.030E-01,
     C 6.505E-01, 8.701E-01, 1.187E+00, 1.587E+00, 2.030E+00,
     C 2.579E+00, 3.028E+00, 3.647E+00, 4.168E+00, 4.627E+00,
     C 5.118E+00, 5.803E+00, 6.553E+00, 7.373E+00, 7.837E+00,
     C 7.800E+00, 7.300E+00, 6.200E+00, 5.250E+00, 4.100E+00,
     C 3.100E+00, 1.800E+00, 1.100E+00, 7.000E-01, 3.000E-01,
     C 2.500E-01, 3.000E-01, 5.000E-01, 7.000E-01, 7.000E-01,
     C 4.000E-01, 2.000E-01, 5.000E-02, 5.000E-03, 5.000E-04/
C     DATA  N2O      /
      DATA RON2O1/
     C 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01,
     C 3.200E-01, 3.200E-01, 3.200E-01, 3.200E-01, 3.195E-01,
     C 3.179E-01, 3.140E-01, 3.095E-01, 3.048E-01, 2.999E-01,
     C 2.944E-01, 2.877E-01, 2.783E-01, 2.671E-01, 2.527E-01,
     C 2.365E-01, 2.194E-01, 2.051E-01, 1.967E-01, 1.875E-01,
     C 1.756E-01, 1.588E-01, 1.416E-01, 1.165E-01, 9.275E-02,
     C 6.693E-02, 4.513E-02, 2.751E-02, 1.591E-02, 9.378E-03,
     C 4.752E-03, 3.000E-03, 2.065E-03, 1.507E-03, 1.149E-03,
     C 8.890E-04, 7.056E-04, 5.716E-04, 4.708E-04, 3.932E-04,
     C 3.323E-04, 2.837E-04, 2.443E-04, 2.120E-04, 1.851E-04/
C     DATA CO        /
      DATA ROCO1/
     C 1.500E-01, 1.450E-01, 1.399E-01, 1.349E-01, 1.312E-01,
     C 1.303E-01, 1.288E-01, 1.247E-01, 1.185E-01, 1.094E-01,
     C 9.962E-02, 8.964E-02, 7.814E-02, 6.374E-02, 5.025E-02,
     C 3.941E-02, 3.069E-02, 2.489E-02, 1.966E-02, 1.549E-02,
     C 1.331E-02, 1.232E-02, 1.232E-02, 1.307E-02, 1.400E-02,
     C 1.498E-02, 1.598E-02, 1.710E-02, 1.850E-02, 2.009E-02,
     C 2.220E-02, 2.497E-02, 2.824E-02, 3.241E-02, 3.717E-02,
     C 4.597E-02, 6.639E-02, 1.073E-01, 1.862E-01, 3.059E-01,
     C 6.375E-01, 1.497E+00, 3.239E+00, 5.843E+00, 1.013E+01,
     C 1.692E+01, 2.467E+01, 3.356E+01, 4.148E+01, 5.000E+01/
C     DATA  CH4      /
      DATA ROCH41/
     C 1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00, 1.700E+00,
     C 1.700E+00, 1.700E+00, 1.699E+00, 1.697E+00, 1.693E+00,
     C 1.685E+00, 1.675E+00, 1.662E+00, 1.645E+00, 1.626E+00,
     C 1.605E+00, 1.582E+00, 1.553E+00, 1.521E+00, 1.480E+00,
     C 1.424E+00, 1.355E+00, 1.272E+00, 1.191E+00, 1.118E+00,
     C 1.055E+00, 9.870E-01, 9.136E-01, 8.300E-01, 7.460E-01,
     C 6.618E-01, 5.638E-01, 4.614E-01, 3.631E-01, 2.773E-01,
     C 2.100E-01, 1.650E-01, 1.500E-01, 1.500E-01, 1.500E-01,
     C 1.500E-01, 1.500E-01, 1.500E-01, 1.400E-01, 1.300E-01,
     C 1.200E-01, 1.100E-01, 9.500E-02, 6.000E-02, 3.000E-02/
C     DATA O2        /
      DATA ROO21/
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05, 2.090E+05,
     C 2.090E+05, 2.090E+05, 2.000E+05, 1.900E+05, 1.800E+05,
     C 1.600E+05, 1.400E+05, 1.200E+05, 9.400E+04, 7.250E+04/
C     DATA DENSITY   /
      DATA DENS1/
     C 2.548E+19, 2.313E+19, 2.094E+19, 1.891E+19, 1.704E+19,
     C 1.532E+19, 1.373E+19, 1.228E+19, 1.094E+19, 9.719E+18,
     C 8.602E+18, 7.589E+18, 6.489E+18, 5.546E+18, 4.739E+18,
     C 4.050E+18, 3.462E+18, 2.960E+18, 2.530E+18, 2.163E+18,
     C 1.849E+18, 1.575E+18, 1.342E+18, 1.144E+18, 9.765E+17,
     C 8.337E+17, 5.640E+17, 3.830E+17, 2.524E+17, 1.761E+17,
     C 1.238E+17, 8.310E+16, 5.803E+16, 4.090E+16, 2.920E+16,
     C 2.136E+16, 1.181E+16, 6.426E+15, 3.386E+15, 1.723E+15,
     C 8.347E+14, 3.832E+14, 1.711E+14, 7.136E+13, 2.924E+13,
     C 1.189E+13, 5.033E+12, 2.144E+12, 9.688E+11, 5.114E+11/

      DO J=1,NLEVEL
         P(J)=P1(J)
         T(J)=T1(J)
         ALT(J)=ALT1(J)
         ROH2O(J)=ROH2O1(J)
         ROCO2(J)=ROCO21(J)
         ROO3(J)=ROO31(J)
         RON2O(J)=RON2O1(J)
         ROCO(J)=ROCO1(J)
         ROCH4(J)=ROCH41(J)
         ROO2(J)=ROO21(J)
         DENS(J)=DENS1(J)
      ENDDO

      RETURN
      END


