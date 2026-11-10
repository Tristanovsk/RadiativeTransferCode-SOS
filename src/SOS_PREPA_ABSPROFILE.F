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
C* FILE: SOS_PREPA_ABSPROFILE.F
C* PROJECT: Successive orders of scattering
C* ROLE: Preparation of the computation of the atmospheric absorption profile
C*
C* AUTHOR:
C*   Initial code: Laboratoire d'Optique Atmospherique (LOA, Villeneuve d'Ascq).
C*   Code rework and development: CS GROUP France, X.Lenot, B.Lafrance
C* 
C*
C* MOD:VERSION:1.0: DATE: 06/06/2013
C*    File creation to adapt the OS to the CKD method (atmospheric absorption)
C*
C* MOD:VERSION:1.1: DATE: 24/02/2016
C*    - Minor adjustments to strictly respect the 72 columns
C*    - Array NEXP changed from Double Precision to Integer
C*
C* MOD:VERSION:1.2: DATE: 06/04/2019
C*    - Correction of the CH4 information in the table ATMOCM
C       ==> ATMOCM(6) changed from the wrong value 13.837E+22 to 3.837E+22
C*
C* MOD:VERSION:2.0: DATE: 02/09/2021
C*    - Adjustments for the removal of the distinction between absorption lines
C*      and continuums, with new CKD data mixing both.
C*      ==> Removal of the array TAUCONT that gave the optical thicknesses
C*          of the absorption continuum.
C*      ==> Extension of the species table IABS to include the gases
C*          absorbing by continuums only (NO2 only)
C*          and removal of the table ICONT.
C*      ==> Removal of the tables AA and BB correcting the temperature
C*          dependence of the CKD coefficients.
C*      ==> Array XKI is replaced by KDIS_KI with the additional Temperature variable.
C*      ==> Array AL is replaced by KDIS_AI
C*      ==> Adjustment of the I/O parameters.
C*      ==> Computation of the wavelength index of the CKD coefficient arrays
C*          from the chosen spectral step and the spectral bounds of the CKD data
C*
C*    - Correction of the weighting of the profiles by the user amounts of H2O and O3.
C*
C*    - Possibility to weight the amounts of CH4 and CO2
C*      ==> Parameters CH4 and CO2 added.
C*      ==> CH4 and CO2 profiles weighted by the ratio between the surface amount
C*          defined by the user and that of the default profile.
C*
C* MOD:VERSION:2.1: 07/07/2022
C*    - Removal of the use of the environment variable SOS_ABS_FIC.
C*    - Constants renamed from SOS_* to CTE_*
C*    - Traces converted to English.
C*
C* MOD:VERSION:2.2: 16/06/2023
C*    - Parameters KDIS_KI_H2O, TAB_CONC_H2O and NB_CONC_H2O added to the call of routine READ_CKD_COEFF
C******************************************************************************

C----------------------------------------------------------------------------
C Definition of constants
C---------------------------------------------------------------------------- 
C Constants used:
C    CTE_LENFIC2 : Length of file names with directory tree.
C    CTE_ABS_NBLEV : number of levels of the absorption profile
c    CTE_NBABS : number of absorbers
C    CTE_ABS_NBCOL : number of columns in an absorption profile given
C		       as input
C    CTE_ABS_H2O  : Flag to take H2O into account (0 no / 1 yes)
C    CTE_ABS_CO2  : Flag to take CO2 into account (0 no / 1 yes)
C    CTE_ABS_O3   : Flag to take O3 into account  (0 no / 1 yes)
C    CTE_ABS_N2O  : Flag to take N2O into account (0 no / 1 yes)
C    CTE_ABS_CO   : Flag to take CO into account  (0 no / 1 yes)
C    CTE_ABS_CH4  : Flag to take CH4 into account (0 no / 1 yes)
C    CTE_ABS_O2   : Flag to take O2 into account  (0 no / 1 yes)
C    CTE_CONT_NO2 : Flag to take NO2 into account (0 no / 1 yes)
C
C

C----------------------------------------------------------------------------
#include "SOS.h"



C==============================================================================
C PROGRAM: SOS_PREPA_ABSPROFILE
C ==========
C      This program prepares the computations for the absorption profile:
C
C        - It reads the profile file: user file or predefined files
C          NB: the predefined files separate the NO2 continuum profile from the other gases
C              --> specific SO2-NO2 file, independent of the profile type.
C
C        - It converts the profile data into numbers of particles/cm2 per layer,
C          consistent with the dimension of the coefficients ki in cm2/particle.
C
C        - It overrides the H2O and O3 profiles with the integrated contents given by the user.
C
C        - It scales the CH4 and CO2 profiles by the ratio between the surface amount
C          defined by the user and that of the default profile.
C
C        - It reads, for all the wavelengths of the domain covered by the CKD files,
C          the coefficients of the CKD files:
C             * Number of exponentials of the CKD decomposition of each absorber
C               for the wavelength range covered.
C             * coefficients Ki and Ai for each exponential.
C
C        - It provides the array index corresponding to the processed wavelength
C
C
C Input data of the program
C ------------------------------
C     -->Parameters/arguments of the program
C
C
C     WA (DOUBLE) : (I) Wavelength (microns).
C
C     PSURF (DOUBLE)  : (I) Atmospheric pressure at the surface (mbar).
C
C     H2O   (DOUBLE)  : (I) Integrated water vapour content (g/cm2).
C
C     O3   (DOUBLE)  : (I) Integrated ozone content (dobson).
C
C     CO2  (DOUBLE)  : (I) Carbon dioxide concentration at the surface (ppmv).
C
C     CH4  (DOUBLE)  : (I) Methane concentration at the surface (ppmv).
C
C     ABSPROFIL (I2) : (I) Type of atmospheric profile for absorption
C		         0:User, 1:TROPICAL, 2:MLS, 3:MLW, 4:SAS, 5:SAW,
C                      6:USST62, 7:No absorption
C
C     FICABSPROFIL (CHARACTER) : (I) associated with ABSPROFIL = 0
C				  complete name of the atmospheric absorption profile file
C				  (directory + file name + extension)
C
C     TRACE (LOGICAL) : (I) If true => trace files are written
C
C     IDLOG (I4) : (I) Identifier of the log file
C
C     NU  (double) : (O)   Wavenumber in cm-1
C
C     LAMB1  (I4) : (O)   Index of the CKD interval
C
C     IABS  (I2) : (O) Array of CTE_NBABS elements. Each element corresponds to an
C			  absorbing gas (Absorbers: H2O CO2 O3 N2O CO CH4 O2 NO2).
C	 		  For each absorber, 1 if processed with the CKD method, 0 otherwise
C			   
C
C     USERPROFIL (double) (O)  Array of (CTE_ABS_NBLEV * CTE_ABS_NBCOL) elements
C			    containing for each absorbing gas the profile expressed in ppmv
C			    It contains, in this order:
C                  		USERPROFIL(J, 1) : altitude (km)
C                  		USERPROFIL(J, 2) : pressure (mbar)
C                  		USERPROFIL(J, 3) : temperature (K)
C                  		USERPROFIL(J, 4) : H2O (ppmv)
C                  		USERPROFIL(J, 5) : CO2 (ppmv)
C                  		USERPROFIL(J, 6) : O3  (ppmv)
C                  		USERPROFIL(J, 7) : N2O (ppmv)
C                  		USERPROFIL(J, 8) : CO  (ppmv)
C                  		USERPROFIL(J, 9) : CH4 (ppmv)
C                  		USERPROFIL(J, 10): O2  (ppmv)
C                 		USERPROFIL(J, 11): density 		(not used)
C                  		USERPROFIL(J, 12): NO2 (ppmv)
C                  		USERPROFIL(J, 13): SO2 (ppmv)	(not used)
C
C
C     ALTABS (double) : (O)  Array of (CTE_ABS_NBLEV) containing the altitude levels of the profile
C
C     RO (double) (O) : Array of (CTE_NBABS * CTE_ABS_NBLEV) elements
C		          containing, for each absorbing gas of the CKD method, the profile
C			   expressed in particles/cm2 for each layer.
C			   It contains, in this order: H2O, CO2, 03, N2O, CO, CH4, O2, NO2
C                       --> Gases with fine line absorption (not only continuums)
C
C     NEXP (I4) :   (O)  Array of (CTE_NBABS * CTE_CKD_NWVL_MAX) elements containing, for each absorber
C		           and each spectral interval, the number of exponentials of the CKD method
C
C
C     KDIS_AI (double) : (O) Array of (CTE_CKD_NAI_MAX * CTE_NBABS * CTE_CKD_NWVL_MAX) elements containing,
C			        for each wavelength, each species and each exponential of the decomposition,
C			        the value of the weighting coefficient ai of the exponential.
C
C     KDIS_KI(double) : (O) Array of the coefficients ki of the exponentials of the decomposition, containing
C                           (CTE_CKD_NT_MAX * CTE_CKD_NP_MAX * CTE_CKD_NAI_MAX * CTE_NBABS * CTE_CKD_NWVL_MAX)
C                           elements covering wavelengths, species, orders i of the decomposition in
C                           exponentials, pressures and temperatures
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
C                                of water vapour, of CTE_CKD_NC_MAX elements (in ppmv)
C     NB_CONC_H2O (integer): (O) Number of actual concentration values
C
C     IER (I4) : (O) error code =0 if no error, =-1 otherwise
C
C
C
C Results provided by the program
C ----------------------------------
C	- Profile of absorbing gases
C	- CKD coefficients
C
C Environment variable
C ------------------------
C   SOS_ABS_ROOT : defines the access path to the main directory of the SOS_ABS code
C   This environment variable is used to locate the files under $SOS_ABS_ROOT/fic
C
C Common used:
C --------------
C     None
C
C Error cases:
C ------------
C     -Writing in a log file
C     -Reading of the user profile file
C
C  Subroutines used:
C --------------------------
C    Functions included in SOS_SUB_TRS.F: DATATM, READ_CKD_COEFF
C
C==============================================================================

      SUBROUTINE SOS_PREPA_ABSPROFILE(
     &			      WA, NUSTEP, PSURF, H2O, O3, CO2, CH4,
     &			      ABSPROFIL, FICABSPROFIL, TRACE, IDLOG,
     &		 	      NU, LAMB1, IABS, USERPROFIL, ALTABS,
     &			      RO, NEXP, KDIS_AI, KDIS_KI, KDIS_KI_H2O,
     &                     TAB_PRES, NB_PRES, TAB_TEMP, NB_TEMP,
     &                     TAB_CONC_H2O,NB_CONC_H2O, IER)
      			 
      IMPLICIT NONE

C* Definition of parameters
      INTEGER*4  NLEVEL  ! nb of altitude levels of the profile
      PARAMETER(NLEVEL=CTE_ABS_NBLEV)

      INTEGER    MAXKD   ! Max number of exponentials per absorbing gas (CKD method)
      PARAMETER (MAXKD=CTE_CKD_NAI_MAX)
      
      INTEGER    MAXP    ! Max number of atmospheric pressure levels (CKD method)
      PARAMETER (MAXP=CTE_CKD_NP_MAX)

      INTEGER    MAXT    ! Max number of atmospheric temperature levels (CKD method)
      PARAMETER (MAXT=CTE_CKD_NT_MAX)




      
C* Definition of variables
C--------------------------------------------

C*   Simu
      DOUBLE PRECISION WA       ! (I) Wavelength (microns).
      DOUBLE PRECISION NU       ! Wavenumber (cm-1)
      DOUBLE PRECISION NUSTEP   ! (I) Wavenumber step (cm-1)
      DOUBLE PRECISION NUMIN, NUMAX ! Min and max wavenumber (cm-1) of the spectral data

C*   Profiles
      INTEGER*2 ABSPROFIL    ! (I) Type of the absorption profile
      CHARACTER*CTE_LENFIC2 FICABSPROFIL ! (I) Absorption profile file provided by the user
      DOUBLE PRECISION ALTC(1:CTE_ABS_NBLEV)	! Altitudes of the absorption profile (increasing order)
      DOUBLE PRECISION ALTABS(1:CTE_ABS_NBLEV)   ! Altitudes of the absorption profile (decreasing order)
      DOUBLE PRECISION USERPROFIL(CTE_ABS_NBLEV,CTE_ABS_NBCOL)  ! absorption profile provided by the user
      DOUBLE PRECISION RO(CTE_NBABS,NLEVEL),DENS(NLEVEL)	 ! Nb of particles/cm2 per species and layer
      DOUBLE PRECISION DP(NLEVEL),P(NLEVEL),T(NLEVEL)

C*   CKD coefficients
      INTEGER*4  NEXP(CTE_NBABS, CTE_CKD_NWVL_MAX) 
                        ! (O) Number of exponentials per wavelength and species

      DOUBLE PRECISION  KDIS_AI(CTE_CKD_NAI_MAX, CTE_NBABS,  
     &                          CTE_CKD_NWVL_MAX) 
                        ! (O) Values of the coefficients ai
                        !     per wavelength, species
                        !     and order i of the decomposition in exponentials
                                          
      DOUBLE PRECISION  KDIS_KI(CTE_CKD_NT_MAX, CTE_CKD_NP_MAX, 
     &                          CTE_CKD_NAI_MAX, CTE_NBABS,  
     &                          CTE_CKD_NWVL_MAX )
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

      INTEGER NB_TEMP, NB_PRES  ! (O) Number of temperature and pressure values in the CKD file
      DOUBLE PRECISION TAB_TEMP(CTE_CKD_NT_MAX) ! (O) Array of the temperature values read
      DOUBLE PRECISION TAB_PRES(CTE_CKD_NP_MAX) ! (O) Array of the pressure values read
 
      INTEGER NB_CONC_H2O  ! (O) Number of H2O concentration values in the H2O CKD file
      DOUBLE PRECISION TAB_CONC_H2O(CTE_CKD_NC_MAX) ! (O) Array of the H2O concentration values read (ppmv)

          
C*   Continuum
      DOUBLE PRECISION WVCONT   ! Wavelength
      
C*   Gaseous absorption
      DOUBLE PRECISION PSURF  ! (I) Atmospheric pressure at the surface (mbar).
      DOUBLE PRECISION H2O    ! (I) Integrated water vapour content (g/cm2).
      DOUBLE PRECISION O3     ! (I) Integrated ozone content (Dobson = 1000 cm.atm).
      DOUBLE PRECISION CO2    ! (I) Carbon dioxide concentration at the surface (ppmv).
      DOUBLE PRECISION CH4    ! (I) Methane concentration at the surface (ppmv).

      INTEGER*2 IABS(1:CTE_NBABS)   ! Array identifying the simulated absorption bands
                                    ! (H2O, CO2, O3, N2O, CO, CH4, O2, NO2)

      DOUBLE PRECISION QO3      ! adjusted ozone content
      DOUBLE PRECISION QH2O     ! adjusted water vapour content

      DOUBLE PRECISION CO2_DEFAULT_PPMV_SURF  ! CO2 concentration at the surface of the default profile (in ppmv)
      DOUBLE PRECISION CH4_DEFAULT_PPMV_SURF  ! CH4 concentration at the surface of the default profile (in ppmv)

      INTEGER*2 NABS            ! index of the absorber
      INTEGER*2 JABS		    ! Status of the absorber (1 if processed, 0 otherwise)
      
C*    Other variables
      INTEGER*4 I,J,II 	! Level.
      INTEGER*4 K

      INTEGER*4 LAMB1

      INTEGER*4 IER     	! error code =0 if no error, =-1 otherwise

      LOGICAL TRACE        	! = true if writing in the trace file

      INTEGER*4 IDLOG	        ! Identifier number of the trace file

      CHARACTER*CTE_LENDIR DIRFIC	! access path to the SO2-NO2 profile file
      INTEGER*4 IFIN       		! nb of significant characters of a string
      
C*    Atmospheric data
      DOUBLE PRECISION ATMOCM(CTE_NBABS) !Coef per gas species: H2O, CO2, O3, N2O, CO, CH4, O2, NO2
                                         ! = (10*Avogadro number) /(Molar mass * gravity acceleration)
      DATA (ATMOCM(I),I=1,CTE_NBABS)/3.410E+22,1.395E+22,1.279E+22,
     &1.395E+22,2.192E+22,3.837E+22,1.918E+22,1.3340E+22/
     
      DOUBLE PRECISION PDSMOL(CTE_NBABS)
      DATA (PDSMOL(I),I=1,CTE_NBABS)/18.,44.,48.,44.,28.,16.,32.,46./



C* Initialisation of the error code
C------------------------------------------------------
      IER=0
      

 	
C* Initialisation of the arrays for the absorbers taken into account
C--------------------------------------------------------------------

C
C-------- ABSORBERS       --------------------------------------
C LINES AND/OR CONTINUUM
C 1-H2O  2-CO2  3-O3  4-N2O  5-CO  6-CH4  7-O2  8-NO2 
      
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C Absorbers taken into account (defined in SOS.h)
      IABS(1)=CTE_ABS_H2O
      IABS(2)=CTE_ABS_CO2
      IABS(3)=CTE_ABS_O3
      IABS(4)=CTE_ABS_N2O
      IABS(5)=CTE_ABS_CO
      IABS(6)=CTE_ABS_CH4
      IABS(7)=CTE_ABS_O2
      IABS(8)=CTE_CONT_NO2
      
C Initialisation of the absorption profile versus altitude
      DO I=1,CTE_ABS_NBLEV
         ALTABS(I)=0.
      ENDDO

C If Psurf, H2O, O3, CO2 or CH4 are not defined, the default
C values of the profiles are used
      
      IF (TRACE) THEN      

         IF (H2O.NE.CTE_NOT_DEFINED_VALUE_DBLE) THEN
            WRITE(IDLOG,*,err=921) 'H2O (g/cm2) :', H2O
         ENDIF            
         
         IF (O3.NE.CTE_NOT_DEFINED_VALUE_DBLE) THEN
            WRITE(IDLOG,*,err=921) 'O3 (cm.atm) :', O3
         ENDIF
       
         IF (CO2.NE.CTE_NOT_DEFINED_VALUE_DBLE) THEN
            WRITE(IDLOG,*,err=921) 'CO2 at surface (ppmv) :', CO2
         ENDIF
         
         IF (CH4.NE.CTE_NOT_DEFINED_VALUE_DBLE) THEN
            WRITE(IDLOG,*,err=921) 'CH4 at surface (ppmv) :', CH4
         ENDIF
         
      ENDIF         
   
   
C Reading of the user profile if it is given as input
c USERPROFIL(J,K), for each level J, enter (K=1,CTE_ABS_NBCOL):
C Z(km),P(hPa),T(K), H2O, CO2, O3, N2O, CO, CH4, O2, density (not used), NO2, SO2
c with gaz concentrations in ppmv.
C
C Or reading of the predefined files.

      IF(ABSPROFIL.EQ.0) THEN !User file case
         write(6,*) ' Reading user file for gas profile:',FICABSPROFIL
         OPEN(2,FILE=FICABSPROFIL,err=1012)
         DO I=1,CTE_ABS_NBLEV
            READ(2,*,err=1021)II,(USERPROFIL(I,K),K=1,CTE_ABS_NBCOL)
         ENDDO
         CLOSE(2)

      ELSE !Case of predefined profile files

         !Name of the data directory
         CALL GETENV('SOS_ABS_ROOT',DIRFIC)
         IFIN=INDEX(DIRFIC,' ')
         IFIN=IFIN-1
         IF (IFIN.LE.0) GOTO 925
      
         IF ((IFIN+4).GT.CTE_LENDIR) GOTO 926
      
         DIRFIC(IFIN+1:IFIN+4)='/fic'
         IFIN=IFIN+4
         IF (IFIN.LE.0) IFIN=CTE_LENDIR
      
         !Reading of the SO2,NO2 profile and computation of the NO2 content in ppm (mg/kg) on the profile
         !--> Fills the information on the NO2 content in RO(8,*)
         WRITE(6,*) '  SOS_PREPA_ABSPROFILE: read profile SO2,NO2'
         OPEN(1,FILE=DIRFIC(1:IFIN)//'/SO2-NO2')
         DO J=1,NLEVEL
            READ(1,*,ERR=927)USERPROFIL(J,13),USERPROFIL(J,12)	!The file gives SO2 in the 1st column
								       !NO2 stored at index 12, SO2 at index 13
         ENDDO
         CLOSE(1)

      ENDIF

C     Retrieval of the gas amounts per level in ppm (mg/kg) on the profile (array RO)
C     for gases 1 to 7: 1-H2O  2-CO2  3-O3  4-N2O  5-CO  6-CH4  7-O2
C     --> Use of the profile in ppmv*10^6 of the array USERPROFIL,
C         filled beforehand by the user file option (ABSPROFIL.EQ.0)
C         or filled directly by the routine DATATM if ABSPROFIL > 0.
C     --> At this stage, the array RO gives the profile per level (not per atmospheric layer).
      CALL DATATM(RO,P,T,ALTC,ABSPROFIL,DENS,USERPROFIL,NLEVEL,PSURF)

C     Conversion of the NO2 amount, not handled by the routine DATATM,
C     to ppm (mg/kg) on the profile ==> Completes the array RO
      DO J=1,NLEVEL           
         RO(8,J)=USERPROFIL(J,12)*1.0E-06*46/28.9
      ENDDO

      IF(CO2.GE.0.) CO2_DEFAULT_PPMV_SURF = RO(2,1)*28.97/(44.0E-06)  ! Conversion back to ppmv
      IF(CH4.GE.0.) CH4_DEFAULT_PPMV_SURF = RO(6,1)*28.97/(16.0E-06)


      
c     Inversion of the altitudes for SOS => from the highest to the lowest
      DO J=1,NLEVEL
	     ALTABS(J)=ALTC(NLEVEL-J+1)
      ENDDO
C
C     Conversion of the profile of contents from ppm per level to particles/cm2 per layer
C                                             ----------                   ----------
      DO J=1,NLEVEL-1
         DP(J)=P(J)-P(J+1)
         DO K=1,CTE_NBABS
            RO(K,J)=DP(J)*(RO(K,J)+RO(K,J+1))/2.0*ATMOCM(K)	!Expressed in particles/cm2 per layer
c            RO(K,J)=RO(K,J)/6.022E+23*PDSMOL(K)			!Expressed in g/cm2
         ENDDO
      ENDDO
      
C     H2O correction
      IF(H2O.GE.0.)THEN
         WRITE(6,*) 
     &   "  Correction of the profile with the user's H2O amount"
         QH2O=0.
         DO J=1,NLEVEL
            QH2O=QH2O+RO(1,J)	! Integrated content of the profile in particles/cm2
         ENDDO

         QH2O=QH2O/6.022E+23*PDSMOL(1)	!Conversion to g/cm2
         DO J=1,NLEVEL
            RO(1,J)=RO(1,J)*H2O/QH2O	! Kept in particles/cm2
            USERPROFIL(J,4)=USERPROFIL(J,4)*H2O/QH2O
         ENDDO
      ENDIF

C
C     O3 correction
      IF(O3.GE.0.)THEN
         WRITE(6,*) 
     &   "  Correction of the profile with the user's O3 amount"
         O3=O3/1000.		! Conversion Dobson ==> cm.atm
         QO3=0.         
         DO J=1,NLEVEL
            QO3=QO3+RO(3,J)	! Integrated content of the profile in particles/cm2
         ENDDO

         QO3=QO3/6.022E+23*PDSMOL(3)	!Conversion particles/cm2 ==> g/cm2
         QO3=QO3*466.23			!Conversion g/cm2 ==> cm.atm
         DO J=1,NLEVEL
            RO(3,J)=RO(3,J)*O3/QO3				! Kept in particles/cm2
            USERPROFIL(J,6)=USERPROFIL(J,6)*O3/QO3
         ENDDO
      ENDIF
   
C     CO2 correction
      IF(CO2.GE.0.)THEN
         DO J=1,NLEVEL
            RO(2,J)=RO(2,J)*CO2/CO2_DEFAULT_PPMV_SURF		! Kept in particles/cm2
            USERPROFIL(J,5)=USERPROFIL(J,5)*CO2/CO2_DEFAULT_PPMV_SURF
         ENDDO         
      ENDIF

C     CH4 correction
      IF(CH4.GE.0.)THEN
         DO J=1,NLEVEL
            RO(6,J)=RO(6,J)*CH4/CH4_DEFAULT_PPMV_SURF		! Kept in particles/cm2
            USERPROFIL(J,9)=USERPROFIL(J,9)*CH4/CH4_DEFAULT_PPMV_SURF
         ENDDO         
      ENDIF

C*******************************************************
C*******************************************************
C*
C* Reading of the CKD coefficients
C* Estimation of the index LAMB1 of the tables
C*
C*******************************************************
C*******************************************************
C     Determination of the wavenumber for the absorption computation
C     (nu = 1 / Wa in cm-1)
      NU=1.0E+4/WA   
      IF ((NU.GT.CTE_CKD_NUMAX).OR.(NU.LT.CTE_CKD_NUMIN)) GOTO 905

c     Reading of the absorption coefficients except H2O
      !WRITE(6,*) '  Read Absorption Coefficients'
      !WRITE(6,*) '  Wave number (cm-1) : ', NU
      DO NABS=1,CTE_NBABS               !Covers the 8 gas species
         JABS=IABS(NABS)
         !print*,"    ===> NABS : ", NABS
         CALL READ_CKD_COEFF(NABS, JABS, NU, NUSTEP,
     &                       NEXP, KDIS_AI, KDIS_KI, KDIS_KI_H2O,
     &                       NUMAX, NUMIN,
     &                       TAB_PRES, NB_PRES, TAB_TEMP, NB_TEMP,
     &				 TAB_CONC_H2O, NB_CONC_H2O, IER)

         IF (IER.NE.0) GOTO 900   ! error case
      ENDDO
       
C     Index of the CKD arrays corresponding to the processed wavenumber NU
      LAMB1 = 1+INT((NUMAX-NU)/NUSTEP)  


C* Writing of the trace file
C-------------------------------------------
      IF (TRACE) THEN      

         WRITE(IDLOG,*,err=921) "Wavelength (nm): ",WA
	  WRITE(IDLOG,*,err=921) "Wave number (cm-1): ",NU
	  WRITE(IDLOG,*,err=921) 
     &      'Spectral resolution of CKD coefficients (cm-1): ',NUSTEP
	  WRITE(IDLOG,*,err=921) 'Index CKD : ',LAMB1
	  WRITE(IDLOG,*,err=921) 'Type of absorption profile: ',
     &   ABSPROFIL
         WRITE(IDLOG,*,err=921) ' '
	 
	  IF (ABSPROFIL.NE.7) THEN
	 
	      WRITE(IDLOG,*,err=921) 
     &       'Calculation of the absorption profile'
	      WRITE(IDLOG,*,err=921) 
     &       '-------------------------------------'
	      WRITE(IDLOG,*,err=921) ' '
	      WRITE(IDLOG,*,err=921) ' Species considered'
             WRITE(IDLOG,16,err=921)'H2O:  ',IABS(1)
	      WRITE(IDLOG,16,err=921)'CO2:  ',IABS(2)
	      WRITE(IDLOG,16,err=921)'O3:   ',IABS(3)
	      WRITE(IDLOG,16,err=921)'N2O:  ',IABS(4)
	      WRITE(IDLOG,16,err=921)'CO:   ',IABS(5)
	      WRITE(IDLOG,16,err=921)'CH4:  ',IABS(6)
	      WRITE(IDLOG,16,err=921)'O2:   ',IABS(7)
	      WRITE(IDLOG,16,err=921)'NO2:  ',IABS(8) 
	      WRITE(IDLOG,*,err=921) ' '

	      IF (PSURF.NE.CTE_NOT_DEFINED_VALUE_DBLE) THEN
               WRITE(IDLOG,*,err=921) 'User Psurf (mbar)', PSURF
	        WRITE(IDLOG,*,err=921) ' '
             ENDIF

             IF (H2O.NE.CTE_NOT_DEFINED_VALUE_DBLE) THEN
               WRITE(IDLOG,*,err=921)  
     $         "Correction of the profile with the user's H2O amount:"
               WRITE(IDLOG,*,err=921) 'User H2O content (g/cm2)', H2O
             ENDIF

             !Computation of the total amount of water vapour on the atmospheric column
             QH2O=0.		!New integrated content recomputed
             DO J=1,NLEVEL
                QH2O=QH2O+RO(1,J)	! Integrated content of the profile in particles/cm2
             ENDDO
             QH2O=QH2O/6.022E+23*PDSMOL(1)	!Conversion to g/cm2
             WRITE(IDLOG,17,err=921)   
     $       '==> Used column amount of Water vapor (g/cm2):',QH2O 
	      WRITE(IDLOG,*,err=921) ' '

             IF (O3.NE.CTE_NOT_DEFINED_VALUE_DBLE) THEN
               WRITE(IDLOG,*,err=921)  
     $         "Correction of the profile with the user's O3 amount:"
               WRITE(IDLOG,*,err=921)   
     $         'User O3 content (Dobsons):', O3*1000.
               WRITE(IDLOG,*,err=921)  
     $         'Default profile Ozone content (Dobsons):',QO3*1000.
             ENDIF

             !Computation of the total amount of ozone on the atmospheric column
             QO3=0.
             DO J=1,NLEVEL
                QO3=QO3+RO(3,J)	
             ENDDO			! 
             QO3=QO3/6.022E+23*PDSMOL(3)	!Conversion particles/cm2 ==> g/cm2
             QO3=QO3*466.23			!Conversion g/cm2==> cm.atm
             WRITE(IDLOG,17,err=921)   
     $       '==> Used column amount of Ozone (Dobsons):    ',QO3*1000. 
	      WRITE(IDLOG,*,err=921) ' '

             IF (CO2.NE.CTE_NOT_DEFINED_VALUE_DBLE) THEN
               WRITE(IDLOG,*,err=921)  
     $         "Correction of the profile with the user's CO2 ratio:"
               WRITE(IDLOG,*,err=921)   
     $         'User CO2 surface concentration (ppmv):', CO2
               WRITE(IDLOG,*,err=921)  
     $         'Default profile CO2 surface concentration (ppmv):',
     $         CO2_DEFAULT_PPMV_SURF
	        WRITE(IDLOG,*,err=921) ' '
             ENDIF

             IF (CH4.NE.CTE_NOT_DEFINED_VALUE_DBLE) THEN
               WRITE(IDLOG,*,err=921)  
     $         "Correction of the profile with the user's CH4 ratio:"
               WRITE(IDLOG,*,err=921)   
     $         'User CH4 surface concentration (ppmv):', CH4
               WRITE(IDLOG,*,err=921)  
     $         'Default profile CH4 surface concentration (ppmv):',
     $         CH4_DEFAULT_PPMV_SURF
	        WRITE(IDLOG,*,err=921) ' '
             ENDIF

	      WRITE(IDLOG,*,err=921) ' '
	      WRITE(IDLOG,*,err=921) ' Gas profile'
	      WRITE(IDLOG,*,err=921)'LEVEL ALT(km) P(hPa) T(K), + 
     $       GASES (ppmv), except column 12: air density'
             DO j=1,CTE_ABS_NBLEV
                WRITE(IDLOG,19,err=921)j,(USERPROFIL(j,k),k=1,13)
             ENDDO

	      WRITE(IDLOG,*,err=921) ' '

	  ENDIF ! end of test on absprofil

      ENDIF ! end of trace
      

C* Nominal end
C-------------------
      GOTO 9999
      
C* Error case and return of status 1 to the shell
C----------------------------------------------
 925  WRITE(6,*) '  SOS_PREPA_ABSPROFILE : ERROR_925: '
      WRITE(6,*) '  => Error while getting SOS_ABS_ROOT variable'
      GOTO 9998

 926  WRITE(6,*) '  SOS_PREPA_ABSPROFILE : ERROR_926: '
      WRITE(6,*) '  => Length of $SOS_ABS_ROOT/fic is too long' 
      WRITE(6,*) '  => CTE_LENDIR has to be adjusted.' 
      GOTO 9998

 927  WRITE(6,*) '  SOS_PREPA_ABSPROFILE : ERROR_927: '
      WRITE(6,*) '  => Error while reading the file fic/SO2-NO2'
      GOTO 9998

  900 WRITE(6,*) '  ERROR in subroutine READ_CKD_COEFF'
      GOTO 9998

  905 WRITE(6,*) '  The simulation wavelength is not included ' //
     &           'in the spectral range of CKD data'
      WRITE(6,*) '  --> NU (cm-1) = ', NU
      WRITE(6,*) '  --> Covered spectral range (cm-1) : ', 
     &           CTE_CKD_NUMIN, CTE_CKD_NUMAX
      GOTO 9998
            
  912 WRITE(6,*) '  ERROR while opening logfile by SOS_PREPAPROFIL'
      GOTO 9998
      
  921 WRITE(6,*) '  ERROR while writing logfile by SOS_PREPAPROFIL'
      GOTO 9998
      
 1012 WRITE(6,*)'  ERROR while opening user profile by SOS_PREPAPROFIL'
      GOTO 9998
      
 1021 WRITE(6,*)'  ERROR while reading user profile by SOS_PREPAPROFIL'
      GOTO 9998

9998  IER=-1
      GOTO 9999

C* Format
C---------   
   16 FORMAT(a6,1x,i2)
   17 FORMAT(a46,1x,f9.3)   
   19 FORMAT(i2.2,1x,f5.1,1x,f10.5,1x,f5.1,1x,11(e11.5,1x))   

 9999 RETURN
	 
      END
      
      
      
  
