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
C* FILE: SOS_PROFIL.F
C* PROJECT: Successive orders of scattering
C* ROLE: Definition of the atmospheric profile.
C*
C* AUTHOR:
C*   Initial code: Laboratoire d'Optique atmosphérique (LOA, Villeneuve d'Ascq).
C*   Code rework and development: CS GROUP France, B.Lafrance, X.Lenot
C* DATE: 30/04/02
C* 
C* MOD:VERSION:1.0: Upgrade of the LOA OS code
C*             2.0: Correction of the application of an underlying molecular layer
C*                  for a 2-layer profile (iprofil=2 and Zmin=0): NBSC_C3.
C*             2.1: 07/07/2008 
C*                  Change of the reading format of the parameters HR and HA
C*                  to correct a decoding error for integer values.
C*		      HR and HA changed from F9.5 to F9
C*
C* MOD:VERSION:3.0: 19/11/2009 
C*       - Change of the argument passing:
C*           - passing by pairs (keyword, value).
C*           - adaptation of the error case management.
C*       - Error messages converted to English.
C*
C* MOD:VERSION:6.0: 06/06/2013
C*   - Absorption taken into account (case 3: molecular and aerosol scattering)
C*
C*   - Change of CTE_DISC to take the absorption optical thickness into account
C*
C*   - Change of the computation of the layer optical thicknesses (true computation N - (N-1))
C*     instead of the computation by derivatives
C*
C*   - The input parameters are given as a list of parameters
C*     and no longer as pairs with associated keyword => all
C*	parameters are therefore mandatory
C*
C*
C* MOD:VERSION:6.1: 19/02/2016
C*   -  Computation of the altitude at each level of the profile:
C*      Initialisation of the profile altitude to ZPROF(0)=300.D+00
C*      to avoid an Inf result that is a problem when reading the file.
C*
C*   -  Added test on the value of TA (AOT) with respect to the number
C*      of aerosol layers for CASE 2 (aerosol profile only) and CASE 4
C*      (homogeneous aerosol profile between two levels), to make sure that
C*      a variable AOT profile can be defined from one level to the next,
C*      given the 0.00001 precision of the optical thickness profile.
C*
C* MOD:VERSION:6.2: 29/05/2020
C*   -  Change of the format of the atmospheric profile file to improve
C*      the precision of the profile (format label 20)
C*
C* MOD:VERSION:7.0: 28/09/2020
C*   -  Adaptation to optimise the definition of the profile according to
C*      the absorption optical thickness.
C*      Use of a threshold (CTE_SEUIL_TAUABS) up to which the profile
C*      is defined by thin layers (of thickness CTE_TCOUCHE).
C*
C* MOD:VERSION:7.1: 31/08/2022
C*    - Constants renamed from SOS_* to CTE_*
C*    - Traces converted to English.
C*
C* MOD:VERSION:7.2: 05/04/2023
C*    - Use of the constant CTE_TOA_ALT to define the max altitude (TOA)
C*
C*
C* MOD:VERSION:8.0: 27/07/2023
C*    - Complete rewrite of the routine to handle the cases of small
C*      scattering optical thickness compared with the gaseous absorption.
C*    - Introduction of the constants CTE_OS_NT_MIN and CTE_TOA_FIRST_LAYER_OPT_THICKNESS
C*      (defined in SOS.h)
C*
C* MOD:VERSION:8.1: 26/09/2023
C*    - Change of the syntax of the test on (TI-TZMOY) in routine SOS_DISC for
C*      compatibility with Fortran 2018 for the gcc/10.2.0 compiler
C*
C* MOD:VERSION:8.2: 12/10/2023
C*    - Correction of an error in the assignment of the altitude values of the profile
C*      for a profile simulating a mixed layer of aerosols and molecules
C*      between two pure molecular layers (case IPROFIL =  2),
C*    - Correction of the proximity test between the altitude of the last computed level
C*      and the limit altitude (ZLIM) of the profile, to remove this level if very close to ZLIM.
C*     
C******************************************************************************

C----------------------------------------------------------------------------
C Definition of constants
C---------------------------------------------------------------------------- 
C Constants used:
C    CTE_LENFIC2 : Length of file names with directory tree.
C    CTE_OS_NT : Max number of layers of the atmospheric profile.
C    CTE_OS_NT_MIN : Min number of layers of the atmospheric profile.
C    CTE_TOA_ALT : Altitude of the TOA (in km).
C    CTE_TOA_FIRST_LAYER_OPT_THICKNESS : Optical thickness of the first atmospheric layer, below TOA.
C    CTE_DELTA_Z : Altitude step (in km) for the search of the first level
C                  of optical thickness CTE_TOA_FIRST_LAYER_OPT_THICKNESS.
C    CTE_THRESHOLD_DZ : Threshold (in km) for the comparison of altitudes when defining the profile.
C    CTE_ABS_NBLEV : Number of levels of the initial absorption profile.
C    CTE_TCOUCHE : Maximum optical thickness of a layer of the atmospheric profile.
C    CTE_THRESHOLD_TAUABS : Threshold on the absorption optical thickness to optimise
C                           the atmospheric profile.
C    CTE_DZTRANSI : Thickness of a transition layer (in km).
C    CTE_PROFIL_MIN_NBC : Minimum number of sub-layers of layers 1 and 3
C                         for an aerosol profile between two altitudes.
C----------------------------------------------------------------------------
#include "SOS.h"
 
      


C==============================================================================
C PROGRAM: SOS_PROFIL
C ==========
C      This program defines the atmospheric profile in optical thickness
C      and in proportion of scatterers (molecules, aerosols) versus
C      the level, taking the atmospheric absorption into account.
C      The profile can also be defined in three
C      layers: molecules, molecules + aerosols (homogeneous mixture) and
C      molecules. The intermediate layer is located by its minimum and
C      maximum altitudes.
C      
C      NB: The profile is not adjusted for a truncation of the phase
C           function. This is done by the OS code.
C
C
C Input data of the program
C ------------------------------
C     -->Parameters/arguments of the program
C
C
C      IPROFIL (I2) : (I) Type of the simulated aerosol profile.
C      		       1 : Aerosol profile with an exponential decrease
C				    given by the scale height
C				2 : Homogeneous aerosol layer between two altitudes
C
C      TR (DOUBLE) : (I) Rayleigh optical thickness.
C	       
C      HR (DOUBLE) : (I) Molecular scale height (km).
C                          Associated with TR >= 0.
C  
C      TA (DOUBLE) : (I) Aerosol optical thickness (not truncated).
C	       	       
C      HA (DOUBLE) : (I) Scale height of the aerosol profile (km).
C                          Associated with TA >= 0. and IPROFIL = 1
C 
C      ZMIN (DOUBLE) : (I) Minimum altitude of the aerosol layer (km).
C                          Associated with IPROFIL = 2
C 	   	 
C      ZMAX (DOUBLE) : (I) Maximum altitude of the aerosol layer (km).
C                          Associated with IPROFIL = 2
C
C      ABSPROFIL (I2) : (I) Type of atmospheric profile for absorption
C			         0:User, 1:TROPICAL, 2:MLS, 3:MLW, 4:SAS, 5:SAW,
C                             6:USST62, 7:No absorption
C
C      ALTABS (double) : (I)  Array of (CTE_ABS_NBLEV) containing the altitude levels of the profile
C
C      TABS (double) : (I) Array of (CTE_ABS_NBLEV-1) containing the absorption optical
C			      thicknesses (continuum + CKD) for each atmospheric layer
C
C	TRACE (LOGICAL) : (I) If true => trace files are written
C
C	IDLOG (I4) : (I) Identifier of the log file
C
C      FICPROFIL (CHARACTER) : complete name of the PROFIL result file
C                              (directory + file name + extension)
C
C      NT (I4) : (O) Number of layers in the atmospheric profile
C
C      IER (I4) : (O) error code =0 if no error, =-1 otherwise
C
C
C Results provided by the program
C ----------------------------------
C   PROFIL result file containing:
C	- the level index
C	- the molecules + aerosols optical thickness for each level
C      - the absorption optical thickness for each level
C      - the total optical thickness for each level
C	- the percentage of aerosols per level
C	- the percentage of molecules per level.
C
C   The number of levels in the profile
C
C Common used:
C --------------
C     None
C
C Error cases:
C ------------
C     -Reading or writing in a file
C     -Invalid parameters
C
C     Message displayed on screen, program stop and
C     return of status 1 to the shell
C
C  Subroutines used:
C --------------------------
C  This program calls the routines:
C     - SOS_DISC
C
C==============================================================================

      SUBROUTINE SOS_PROFILE(IPROFIL, TR, HR, TA, HA, 
     &		               ZMIN, ZMAX, ABSPROFIL, ALTABS, TABS,
     &		               TRACE, IDLOG, FICPROFIL, NT, IER) 

      
      IMPLICIT NONE


C* Definition of variables
C--------------------------------------------


C*   Profiles
      INTEGER*2 IPROFIL	 ! (I) Type of the aerosol profile
      DOUBLE PRECISION TTOT  ! total optical thickness.
      INTEGER*4 NT           ! Number of layers in the profile
      INTEGER*4 NT_NG        ! Number of layers in the profile without gaseous absorption
      
C*   molecules:
      DOUBLE PRECISION TR    ! (I) Rayleigh optical thickness.
      DOUBLE PRECISION HR    ! (I) Scale height of the molecular profile (km).

C*   Aerosols:
      DOUBLE PRECISION TA    ! (I) Aerosol optical thickness (not truncated).
      DOUBLE PRECISION HA    ! (I) Scale height of the aerosol profile (km).
      DOUBLE PRECISION ZMIN  ! (I) Minimum altitude of the molecules + aerosols layer.
      DOUBLE PRECISION ZMAX  ! (I) Maximum altitude of the molecules + aerosols layer.
      
      DOUBLE PRECISION H(0:CTE_OS_NT)     ! Optical thickness at each level.
      DOUBLE PRECISION PCAER(0:CTE_OS_NT) ! Percentage of aerosols per level.
      DOUBLE PRECISION PCMOL(0:CTE_OS_NT) ! Percentage of molecules per level.

      DOUBLE PRECISION H_NG(0:CTE_OS_NT)     ! Optical thickness at each level of the profile without gaseous absorption.
      DOUBLE PRECISION PCAER_NG(0:CTE_OS_NT) ! Percentage of aerosols per level of the profile without gaseous absorption.
      DOUBLE PRECISION PCMOL_NG(0:CTE_OS_NT) ! Percentage of molecules per level of the profile without gaseous absorption.
      	
C*   Gaseous absorption
      DOUBLE PRECISION TABS(1:CTE_ABS_NBLEV)     ! Absorption optical thickness
      DOUBLE PRECISION ALTABS(1:CTE_ABS_NBLEV)   ! Altitudes of the absorption profile
      DOUBLE PRECISION ZLIM    ! Minimum altitude of definition of the profile (not null if very strong absorption optical thickness).
      DOUBLE PRECISION TG_ZLIM ! Absorption optical thickness at altitude ZLIM
      DOUBLE PRECISION ALIN, BLIN ! Linear interpolation coefficients of the gaseous optical thickness on altitude
      INTEGER*2 ABSPROFIL    ! (I) Type of the absorption profile


      DOUBLE PRECISION TTOT_ZLIM ! Total optical thickness at altitude ZLIM
      
C*   Variables specific to the profile with an aerosol layer between two altitudes:
      DOUBLE PRECISION VR_C1	! Molecular optical thickness of layer 1.
      DOUBLE PRECISION VR_C2	! Molecular optical thickness of layer 2.
      DOUBLE PRECISION VR_C3	! Molecular optical thickness of layer 3.
      DOUBLE PRECISION VR_SC	! Molecular optical thickness between two levels.

      INTEGER*4 NBSC_C1		! Number of sub-layers in layer 1
      INTEGER*4 NBSC_C2		! Number of sub-layers in layer 2
      INTEGER*4 NBSC_C3 	       ! Number of sub-layers in layer 3
      INTEGER*4 NB_TR	 	! Number of transition layers


      DOUBLE PRECISION DELTA_Z	! Altitude difference between two levels (in km).
      
C*   Other variables
      DOUBLE PRECISION Z 		! Altitude of a level (in km).
             
      DOUBLE PRECISION Hmol(0:CTE_OS_NT)   ! Molecular optical thickness profile.
      DOUBLE PRECISION Haer(0:CTE_OS_NT)   ! Aerosol optical thickness profile.
      DOUBLE PRECISION Habs(0:CTE_OS_NT)   ! Optical thickness profile of the absorbing gases

      DOUBLE PRECISION Hmol_NG(0:CTE_OS_NT)   ! Molecular optical thickness profile without gaseous absorption.
      DOUBLE PRECISION Haer_NG(0:CTE_OS_NT)   ! Aerosol optical thickness profile without gaseous absorption.

      
      DOUBLE PRECISION ZPROF(0:CTE_OS_NT)  ! Altitudes of the profile levels (in km)
      DOUBLE PRECISION ZPROF_NG(0:CTE_OS_NT)   ! Altitudes of the profile levels without gaseous absorption (in km)
      DOUBLE PRECISION Zing			 ! Altitude value of the profile without gas
      INTEGER*4 ing                        ! Index of ZPROF_NG giving Zing

      DOUBLE PRECISION VA,VR,VG
            
      DOUBLE PRECISION  ZZ  ! Linear interpolation coefficient between 2 altitudes
      INTEGER*4 I,J  	! Level.
      INTEGER*4 II          ! Level index
      INTEGER*4 K  		! Elements of the absorbing gas profile
      
      INTEGER*4 IER     	! error code =0 if no error, =-1 otherwise
      
      LOGICAL TRACE        	! = true if writing in the trace file

      INTEGER*4 IDLOG	! Identifier number of the trace file

      CHARACTER*CTE_LENFIC2 FICPROFIL       ! path and name of the PROFIL result file

      DOUBLE PRECISION DTAU			!Optical thickness from TOA to a given altitude

      DOUBLE PRECISION T_FIRST_LAYER 	!Optical thickness of the first layer at TOA
      DOUBLE PRECISION T_LAYER 		!Optical thickness of the layers other than the first

      DOUBLE PRECISION T_FIRST_LAYER_NG 	!Optical thickness of the layer at TOA (without gaseous absorption)
      DOUBLE PRECISION T_LAYER_NG 		!Optical thickness of the layers other than the first (without gaseous absorption)

      CHARACTER*100 CAS_TXT_NG		!Text specifying the case kept without gas
      CHARACTER*100 CAS_TXT 		!Text specifying the case kept with gas
      INTEGER*2 LTXT				!Length of the character string

C* Initialisation of the error code
C------------------------------------------------------
      IER=0
  
      IF((IPROFIL.NE.1).AND.(IPROFIL.NE.2)) GOTO 940
               
      IF(IPROFIL.EQ.2) GOTO 1000	!The processing below of CASE IPROFIL=1 is skipped


C===============================================================================
C Step 1: Computation of the profile WITHOUT gaseous absorption
C ==> It is also computed when a profile with absorption is defined
C ==> Estimation of the profiles: ZPROF_NG, Hmol_NG, Haer_NG
C===============================================================================
     
      !-------------------------------------------------------------------------
      ! Determination of the number of layers and of the optical thicknesses:
      !  - Case without gaseous absorption: NT_NG
      !  - Estimation of the opt. thick. of the first layer at TOA: T_FIRST_LAYER_NG
      !  - Estimation of the opt. thick. of the underlying layers: T_LAYER_NG
      !-------------------------------------------------------------------------
      IF (CTE_TOA_FIRST_LAYER_OPT_THICKNESS.GT.CTE_TCOUCHE) GOTO 950

      ! Without gaseous absorption
      ! =======================
      TTOT = TR+TA
      IF ((TTOT/CTE_OS_NT_MIN).LE. 
     &     CTE_TOA_FIRST_LAYER_OPT_THICKNESS) THEN
         !Case of a very weakly scattering medium: all the layers have the same opt. thick., even the first one
         CAS_TXT_NG = 
     &   "NO GAS case 1: very low scattering optical thickness"      
         NT_NG = CTE_OS_NT_MIN
         T_LAYER_NG = TTOT / NT_NG
         T_FIRST_LAYER_NG = T_LAYER_NG
         
      ELSE
         IF ((TTOT/CTE_OS_NT_MIN).LT.CTE_TCOUCHE) THEN
            !Case of a weakly scattering medium:
            !    the first layer is kept separate with an opt. thick. of CTE_TOA_FIRST_LAYER_OPT_THICKNESS
            !    CTE_OS_NT_MIN layers of equal opt. thick. are taken below (< CTE_TCOUCHE)
            CAS_TXT_NG="NO GAS case 2: low scattering optical thickness"  
            NT_NG = CTE_OS_NT_MIN + 1
            T_FIRST_LAYER_NG = CTE_TOA_FIRST_LAYER_OPT_THICKNESS
            T_LAYER_NG = (TTOT-T_FIRST_LAYER_NG) / CTE_OS_NT_MIN
            
         ELSE
            !Case TTOT/CTE_TCOUCHE > CTE_OS_NT_MIN
            !    a separate first layer is defined with an opt. thick. of CTE_TOA_FIRST_LAYER_OPT_THICKNESS
            !    (NT_NG-1) layers of equal opt. thick. are taken below, T_LAYER_NG (close to CTE_TCOUCHE)
            CAS_TXT_NG="NO GAS case 3: more layers than CTE_OS_NT_MIN-"
            LTXT=INDEX(CAS_TXT_NG,'-')-1 
            CAS_TXT_NG = 
     &      CAS_TXT_NG(1:LTXT)//" for layers of opt. thick.-"
            LTXT=INDEX(CAS_TXT_NG,'-')-1 
            CAS_TXT_NG =  CAS_TXT_NG(1:LTXT)//" CTE_TCOUCHE"    

            T_FIRST_LAYER_NG = CTE_TOA_FIRST_LAYER_OPT_THICKNESS	!Opt. thick. of the 1st layer

            NT_NG = INT((TTOT-T_FIRST_LAYER_NG)/CTE_TCOUCHE)
            T_LAYER_NG =  (TTOT-T_FIRST_LAYER_NG) / NT_NG	!Opt. thick. of the layers below the ultra-thin layer at TOA
	
            NT_NG = NT_NG + 1	!Addition of the ultra-thin layer at TOA

         ENDIF
      ENDIF

    
      !-------------------------------------------------------------------------
      ! Computation of the profile without gaseous absorption
      ! ==> It is also computed when a profile with absorption is defined
      !     to store the profile without gas ZPROF_NG that will be used
      !     for the definition of the profile with gas
      ! ==> Outputs the number of layers NT_NG
      !     and the profiles Hmol_NG, Haer_NG, H_NG, ZPROF_NG, PCMOL_NG, PCAER_NG
      !-------------------------------------------------------------------------

      TTOT = TR+TA

      IF (TA.EQ.0.0D+00) THEN
          ! ---------------------------------------------------
          ! CASE 1: Case of a layer without aerosol scattering:
          !        No absorption ==> the profile is homogeneous.
          ! ---------------------------------------------------
  
          Hmol_NG(0) = 0.D+00
          Hmol_NG(1) = T_FIRST_LAYER_NG

          DO I=2,NT_NG
             Hmol_NG(I) = (I-1)*T_LAYER_NG + T_FIRST_LAYER_NG
          ENDDO

          DO I=0,NT_NG
             PCMOL_NG(I) = 1.000
             PCAER_NG(I) = 0.000
             Haer_NG(I) = 0.D+00
             H_NG(I) = Hmol_NG(I)
          ENDDO

          ZPROF_NG(0) = CTE_TOA_ALT
          DO I=1,NT_NG         
             ZPROF_NG(I)=HR*DLOG(TR/Hmol_NG(I))
          ENDDO
 
      ELSE
          ! ---------------------------------------------------
          ! CASE 2: Case of a layer with aerosol scattering
          !        and molecules: No absorption
          ! ---------------------------------------------------

          ZPROF_NG(0) = CTE_TOA_ALT
          Hmol_NG(0) = 0.D+00
          Haer_NG(0) = 0.D+00
          H_NG(0) = 0.D+00

          !Search for the altitude of level 1 whose optical thickness (mol+aer) will be T_FIRST_LAYER_NG
          DTAU = 0.
          Z = CTE_TOA_ALT 
          DO WHILE(DTAU.LT.T_FIRST_LAYER_NG)
               Z = Z - CTE_DELTA_Z !km
               DTAU = TR*DEXP(-Z/HR) + TA*DEXP(-Z/HA) 
	   ENDDO
          ZPROF_NG(1) = Z 
          VR = TR*DEXP(-Z/HR) 
          VA = TA*DEXP(-Z/HA)
          Hmol_NG(1) = VR
          Haer_NG(1) = VA
          H_NG(1) = DTAU
          PCMOL_NG(1) = VR/DTAU
          PCAER_NG(1) = VA/DTAU

          !Duplication of the proportions on level 0 at TOA (ultra-thin underlying layer)
          PCMOL_NG(0) = PCMOL_NG(1)
          PCAER_NG(0) = PCAER_NG(1)

          !Search for the altitude of levels 2 to NT_NG-1
          !and computation of the profile parameters for each level
          DO I=2,NT_NG-1
               CALL SOS_DISC(T_LAYER_NG,TA,HA,TR,HR,TABS,ALTABS,
     &                       H_NG(I-1),ZPROF_NG(1),0.D+00,0.D+00,Z)	
               ZPROF_NG(I) = Z  
               VR = TR*DEXP(-Z/HR) 
               VA = TA*DEXP(-Z/HA)
               Hmol_NG(I) = VR
               Haer_NG(I) = VA
               H_NG(I) =  VR + VA

               VR = VR - Hmol_NG(I-1)
               VA = VA - Haer_NG(I-1)
               PCMOL_NG(I) = VR / (VR+VA)
               PCAER_NG(I) = VA / (VR+VA)
	   ENDDO 
 
          !Management of the surface level (at NT_NG)
          ZPROF_NG(NT_NG) = 0.D+00 
          Hmol_NG(NT_NG) = TR
          Haer_NG(NT_NG) = TA
          H_NG(NT_NG) =  TR +TA

          VR = TR - Hmol_NG(NT_NG-1)
          VA = TA - Haer_NG(NT_NG-1)
          PCMOL_NG(NT_NG) = VR / (VR+VA)
          PCAER_NG(NT_NG) = VA / (VR+VA)

      ENDIF !End of test TA = 0.


      IF ((ABSPROFIL.EQ.7).OR.(TABS(CTE_ABS_NBLEV).EQ.0.D+00)) THEN
	  !Without absorption, the results of the profile without gas are stored
         !as final results

        NT = NT_NG
              
         DO I=0,NT
            ZPROF(I) = ZPROF_NG(I)   
	     Hmol(I) = Hmol_NG(I)
	     Haer(I) = Haer_NG(I)
	     Habs(I) = 0.D+00
	     H(I) = Hmol(I) + Haer(I)    
            PCAER(I) = PCAER_NG(I)
            PCMOL(I) = PCMOL_NG(I)
         ENDDO

      ELSE
C===============================================================================
C Step 2: Computation of the profile WITH gaseous absorption
C ==> If gaseous absorption is present
C ==> Estimation of the profiles: ZPROF, Hmol, Haer, PCMOL, PCAER
C===============================================================================
       


         IF (CTE_TCOUCHE.GT.CTE_THRESHOLD_TAUABS) GOTO 960
    
         !-------------------------------------------------------------------------
         ! Determination of the optical thicknesses with gaseous absorption:
         !  - Estimation of the opt. thick. of the first layer at TOA: T_FIRST_LAYER
         !  - Estimation of the opt. thick. of the underlying layers: T_LAYER
         !-------------------------------------------------------------------------
  
         ! Total optical thickness along the path
         TTOT = TR + TA + TABS(CTE_ABS_NBLEV)

         IF (TABS(CTE_ABS_NBLEV).GT.CTE_THRESHOLD_TAUABS) THEN   
            !Case of a strong absorption opt. thick. > threshold
            !The absorption profile is taken into account down to the altitude ZLIM to be estimated
            CAS_TXT =  
     &      "WITH GAS case 0 : Strong absorption ==> Limit layer used"
 
            !Determination of the altitude at which the gaseous absorption optical thickness
            !exceeds the threshold
            I=1
            DO WHILE (TABS(I).LT.CTE_THRESHOLD_TAUABS) 
               I=I+1
            ENDDO

            ALIN = (TABS(I) - TABS(I-1)) / (ALTABS(I) - ALTABS(I-1))
            BLIN = TABS(I) - ALIN*ALTABS(I)
          	
            TG_ZLIM = CTE_THRESHOLD_TAUABS
            ZLIM = ( TG_ZLIM - BLIN ) / ALIN


            T_FIRST_LAYER = CTE_TOA_FIRST_LAYER_OPT_THICKNESS

            TTOT_ZLIM = TA*DEXP(-ZLIM/HA)+TR*DEXP(-ZLIM/HR)+TG_ZLIM 
            T_LAYER  = (TTOT_ZLIM-T_FIRST_LAYER) / (CTE_OS_NT-NT_NG-2) !Max nb of layers - those of the profile without gas - the TOA one - the one below the opt. thick. limit
            T_LAYER  = MAX(T_LAYER, CTE_TCOUCHE)
 

         ELSE
            !Case of an absorption opt. thick. below the threshold
            !The absorption profile is taken into account down to the ground

            ZLIM = 0.		
            TG_ZLIM = TABS(CTE_ABS_NBLEV)

    
            IF ((TTOT/CTE_OS_NT_MIN).LE. 
     &           CTE_TOA_FIRST_LAYER_OPT_THICKNESS) THEN
               !Case of a medium with very small opt. thick.: all the layers have the same opt. thick., even the first one
               CAS_TXT =  "WITH GAS case 1: No limit layer -"
               LTXT=INDEX(CAS_TXT,'-')-1 
               CAS_TXT = 
     &         CAS_TXT(1:LTXT)//" + very low optical thickness"

               NT  = CTE_OS_NT_MIN
               T_LAYER  = TTOT / NT 
               T_FIRST_LAYER = T_LAYER
            ELSE
               IF ((TTOT/CTE_OS_NT_MIN).LT.CTE_TCOUCHE) THEN
                  CAS_TXT =  "WITH GAS case 2: No limit layer -"
                  LTXT=INDEX(CAS_TXT,'-')-1 
                  CAS_TXT = CAS_TXT(1:LTXT)//" + low optical thickness"

                  !Case of a medium with small opt. thick.:
                  !      the first layer is kept separate with an opt. thick. of CTE_TOA_FIRST_LAYER_OPT_THICKNESS
                  !      CTE_OS_NT_MIN layers of equal opt. thick. are taken below (< CTE_TCOUCHE)
                  NT = CTE_OS_NT_MIN + 1
                  T_FIRST_LAYER = CTE_TOA_FIRST_LAYER_OPT_THICKNESS
                  T_LAYER = (TTOT-T_FIRST_LAYER) / CTE_OS_NT_MIN

               ELSE
                  !Case TTOT/CTE_TCOUCHE > CTE_OS_NT_MIN
                  CAS_TXT= 
     &            "WITH GAS case 3: more layers than CTE_OS_NT_MIN-"
                  LTXT=INDEX(CAS_TXT,'-')-1 
                  CAS_TXT = 
     &            CAS_TXT(1:LTXT)//" for layers of opt. thick.-"
                  LTXT=INDEX(CAS_TXT,'-')-1 
                  CAS_TXT =  CAS_TXT(1:LTXT)//" CTE_TCOUCHE"

                  T_FIRST_LAYER = CTE_TOA_FIRST_LAYER_OPT_THICKNESS	!Opt. thick. of the 1st layer

                  NT = INT((TTOT-T_FIRST_LAYER)/CTE_TCOUCHE)	
                  T_LAYER =  (TTOT-T_FIRST_LAYER) / NT		!Opt. thick. of the layers below the ultra-thin layer at TOA

                  NT = NT + 1	!+1 to add the ultra-thin layer at TOA

               ENDIF
            ENDIF ! End of test on the value of TTOT

         ENDIF !End of test absorption opt. thick. > threshold
		 

    
         !-------------------------------------------------------------------------
         ! Computation of the profile WITH gaseous absorption
         !
         !    Computation of the profile altitudes ZPROF(I)
         !    Computation of the optical thicknesses Hmol(I), Haer(I) and Habs(I)
         !    and of the proportions of molecules and aerosols (PCMOL(I) and PCAER(I))
         !
         !    Computation of the number of atmospheric layers NT:
         !       Max number of layers = CTE_OS_NT
         !       Min number of layers = CTE_OS_NT_MIN
         !       Number of layers depending on a possible exceedance of the absorption optical thickness threshold
         !  ------------------------------------------

         !-------------------------------------------------------------------------
 

         ! Computation of the optical thickness and of the scattering contribution
         ! of molecules and aerosols from TOA (300 km) to the ground.
         ! ----------------------------------------------------------------
         NT = 1
         Z = CTE_TOA_ALT
         Zing = ZPROF_NG(1)
         Hmol(0) = 0.D+00
         Haer(0) = 0.D+00
         Habs(0) = 0.D+00
         H(0) = 0.D+00
	  I = 0

         TTOT_ZLIM = TA*DEXP(-ZLIM/HA)+TR*DEXP(-ZLIM/HR)+TG_ZLIM  
		
         DO WHILE ((TTOT_ZLIM-H(NT-1)).GT.T_LAYER)
 		
            I=NT
            IF (I.EQ.1) THEN
               !Search for the altitude of level 1 whose optical thickness (mol+aer+gas) will be T_FIRST_LAYER

               DTAU = 0.
               DO WHILE(DTAU.LT.T_FIRST_LAYER)

                  Z = Z - CTE_DELTA_Z

                  !Search for the index J such that ALTABS(J-1) > z >= (ALTABS(J)
	           J=2
	           DO WHILE(Z.LT.ALTABS(J))				
		       J=J+1
	           ENDDO
    
                  ! Linear interpolation of the gaseous absorption optical thickness VG between the two levels of the gas profile
                  IF (Z.LE.ALTABS(1)) THEN
                     ZZ=(Z-ALTABS(J-1))/(ALTABS(J)-ALTABS(J-1))
	              VG=(1-ZZ)*TABS(J-1)+ZZ*TABS(J)
                  ELSE
                     VG=0.
	           ENDIF
 

                  ! Molecular and aerosol optical thickness
                  VR=TR*DEXP(-Z/HR)
                  VA=TA*DEXP(-Z/HA)       

                  ! Total optical thickness
	           DTAU = VR + VA + VG
 
	        ENDDO
               ZPROF(1) = Z
               H(1) = DTAU
               ing = 1
         
            ELSE ! Determination of the altitude for the optical thickness at level I (for I > 2)

               CALL SOS_DISC(T_LAYER,TA,HA,TR,HR,TABS,ALTABS,H(I-1),
     &                       ZPROF(1),TG_ZLIM,ZLIM,Z)	      


            ENDIF !End of distinction between cases I=1 and I>1


	     IF (Z.LE.Zing) THEN
                 Z=Zing	        !If Z is below a level Zing of the profile without gas: it is forced into the profile with gas
                 ing = ing + 1      !The next (lower) level Zing is then considered
                 Zing = ZPROF_NG(ing) 
            ELSE
                IF ((Z-Zing).LE.CTE_THRESHOLD_DZ) THEN
                   ing = ing + 1           !Z very close to Zing, the next Z is not forced to the value Zing:
                   Zing = ZPROF_NG(ing)    !Zing is changed by taking the next (lower) value
                ENDIF               
            ENDIF  



            ZPROF(I)=Z

            ! Computation of the absorption optical thickness by linear interpolation
            ! Search for the altitude levels of the absorption profile around the altitude Z
	     J=2
	     DO WHILE(Z.LT.ALTABS(J))
		 J=J+1
	     ENDDO
 
            ! Linear interpolation between the two layers
            IF (Z.GT.ALTABS(1)) THEN
               ! Management of the case where Z is outside the profile (Z > 120 km)
		 !print*,"!!!Case Z=",Z," > 	ALTABS(1) =",ALTABS(1)
	        VG=TABS(J-1)
	     ELSE
	        ZZ=(Z-ALTABS(J-1))/(ALTABS(J)-ALTABS(J-1))
	        VG=(1-ZZ)*TABS(J-1)+ZZ*TABS(J)
	     ENDIF

	    
            ! Computation of the molecular and aerosol proportions.
            VR=TR*DEXP(-Z/HR)
            VA=TA*DEXP(-Z/HA)
	     Hmol(I)=VR
	     Haer(I)=VA
	     Habs(I)=VG
            H(I)=VA+VR+VG

	     VA=VA-Haer(I-1)
            VR=VR-Hmol(I-1)
	     VG=VG-Habs(I-1)
	     PCAER(I)=VA/(VA+VR+VG)
            PCMOL(I)=VR/(VA+VR+VG)

            NT = NT + 1

         ENDDO !End of loop estimating the level Z while Z > T_LAYER

         !Management of the limit depth level (ZLIM!=0. or surface)
         IF ((ZPROF(NT-1)-ZLIM).LE.CTE_THRESHOLD_DZ) THEN
            NT=NT-1	!If the difference between Z and ZLIM is small, the last computed level is ignored
         ENDIF
         ZPROF(NT) = ZLIM  
         VR = TR*DEXP(-ZLIM /HR)
         VA = TA*DEXP(-ZLIM /HA)
         VG = TG_ZLIM
         Hmol(NT) = VR
         Haer(NT) = VA
         Habs(NT) = VG
         H(NT) =  VR + VA + TG_ZLIM

         VA=VA-Haer(NT-1)
         VR=VR-Hmol(NT-1)
	  VG=VG-Habs(NT-1)
	  PCAER(NT)=VA/(VA+VR+VG)
         PCMOL(NT)=VR/(VA+VR+VG)

	  ! Definition of the profile values at TOA
         ZPROF(0)= CTE_TOA_ALT
         PCAER(0) = PCAER(1)
         PCMOL(0) = PCMOL(1)
	  Hmol(0) = 0.
	  Haer(0) = 0.
	  Habs(0) = 0.
         H(0) = 0.

	  

         IF (TABS(CTE_ABS_NBLEV).GT.CTE_THRESHOLD_TAUABS) THEN   
            !The profile is completed from the level where the gaseous optical thickness reaches the threshold
            !down to the ground by adding an extra layer.
            !This last layer (from the limit level to the ground) is given all the contents below the limit level.
            !==> assumes an attenuation such that this layer is not visible.
            NT = NT + 1
            Hmol(NT)=TR
            Haer(NT)=TA
            Habs(NT)=TABS(CTE_ABS_NBLEV)
            H(NT)=Hmol(NT)+Haer(NT)+Habs(NT)

            VR=Hmol(NT)-Hmol(NT-1)
            VA=Haer(NT)-Haer(NT-1)
            VG=Habs(NT)-Habs(NT-1)

            PCAER(NT)=VA/(VA+VR+VG)
            PCMOL(NT)=VR/(VA+VR+VG)

            ZPROF(NT)=0.D+00

         ENDIF

      ENDIF !End of test case without/with gaseous absorption


      GOTO 7777 !Case IPROFIL=2 is skipped



     
1000  CONTINUE
C======================================================
C* CASE IPROFIL=2 :
C*   Aerosol and molecular scattering.
C*   Uniform aerosol layer between two altitudes
C*   Case without gaseous absorption.
C======================================================
C*      
C*		T O A
C*	___________________________ 
C*
C*		layer 1
C*	--------------------------- Zmax+DZtransi
C*	___________________________ Zmax
C*
C*
C*		layer 2
C*
C*	___________________________ Zmin	
C*	--------------------------- Zmin-DZtransi
C*
C*		layer 3
C*	___________________________

	         
	 
C* Validity check of the arguments
C---------------------------------------
      IF ((ZMIN.LT.0.D-00).OR.(ZMAX.LE.ZMIN)) GOTO 1010


C* Number of layers of the profile
C---------------------------------------
      TTOT = TR + TA 
      NT = int(TTOT / CTE_TCOUCHE) 
      IF (NT.GT.CTE_OS_NT) NT = CTE_OS_NT
          
 
C* Molecular optical thickness of the layers
C  -----------------------------------------
      VR_C1=TR*DEXP(-(ZMAX+CTE_DZTRANSI)/HR)
      VR_C2=TR*( DEXP(-ZMIN/HR)-DEXP(-(ZMAX+CTE_DZTRANSI)/HR) )
      IF(ZMIN.EQ.0.D-00) THEN
	  VR_C3 = 0.D+00
	  NB_TR = 1
      ELSE   
	  VR_C3=TR*( 1.D+00-DEXP(-(ZMIN-CTE_DZTRANSI)/HR) )
	  NB_TR = 2     
      ENDIF 
	   
C* Number of sub-layers in each layer
C  -----------------------------------------

      NBSC_C1 = (NT-NB_TR)*VR_C1/(TR+TA) 
      NBSC_C1 = MAX( CTE_PROFIL_MIN_NBC, NBSC_C1 )
	 
      IF(ZMIN.EQ.0.D+00) THEN
	  NBSC_C3 = 0	    
      ELSE
	  NBSC_C3 = (NT-NB_TR)*VR_C3/(TR+TA)	 
	  NBSC_C3 = MAX( CTE_PROFIL_MIN_NBC, NBSC_C3 )
      ENDIF
	 	 
      NBSC_C2 = (NT-NB_TR)-NBSC_C1-NBSC_C3

      IF((TA/NBSC_C2).LT.0.00001) GOTO 1020

C* Molecular profile: layer 1
C  -----------------------------
      VR_SC=VR_C1/NBSC_C1

      DO I=1,NBSC_C1
	  Hmol(I)=Hmol(I-1)+VR_SC
	  Haer(I)=0.
      ENDDO
	       
C* Transition: ZMAX+CTE_DZTRANSI -> ZMAX
C  -------------------------------------
      I=NBSC_C1+1
      Hmol(I)=TR*DEXP(-ZMAX/HR)
      VR_SC=Hmol(I)-Hmol(I-1)
      Haer(I)=Haer(I-1)+(TA*VR_SC/VR_C2)
	      	       
C* Homogeneous mixture: aerosols + molecules: layer 2
C  ------------------------------------------------
      DELTA_Z=(ZMAX-ZMIN)/NBSC_C2
	       
      !Initialisation
      Z=ZMAX
	       
      DO I = NBSC_C1+2 , NBSC_C1+NBSC_C2+1
	         
	  Z=Z-DELTA_Z
	  Hmol(I)=TR*DEXP(-Z/HR)
	  VR_SC=Hmol(I)-Hmol(I-1)
	  Haer(I)=Haer(I-1)+(TA*VR_SC/VR_C2)
	        	  
      ENDDO

	
      IF (ZMIN.NE.0.D+00) THEN       	       
	    
C* Transition: ZMIN -> ZMIN-CTE_DZTRANSI
C  -------------------------------------
	  I=NBSC_C1+NBSC_C2+2
	  Hmol(I)=TR*DEXP(-(ZMIN-CTE_DZTRANSI)/HR)
	  Haer(I)=Haer(I-1)  

C* Molecular profile: layer 3
C  -----------------------------
	       
	  VR_SC=VR_C3/NBSC_C3
	       
	  DO I = NBSC_C1+NBSC_C2+3 , NT
            Hmol(I)=VR_SC+Hmol(I-1)
            Haer(I)=Haer(I-1)
	  ENDDO
	       
      ENDIF !(Zmin.ne.0)


	       	       	       
C* Computation of H, PCAER, PCMOL
C  -----------------------	
      ZPROF(0)=CTE_TOA_ALT
      Hmol(0)=TR*DEXP(-CTE_TOA_ALT/HR)
      Haer(0)=0.
      H(0)=Hmol(0)+Haer(0)     
      PCAER(0)=0.
      PCMOL(0)=1.
         
	
      DO I=1,NT
	       	  
         H(I)=Hmol(I)+Haer(I)
         ZPROF(I)=HR*DLOG(TR/Hmol(I))

         IF(Haer(I).EQ.Haer(I-1)) THEN
            PCAER(I)=0.
            PCMOL(I)=1
         ELSE
            PCMOL(I)=1/(1+(TA/VR_C2))
            PCAER(I)=1-PCMOL(I)
         ENDIF
		    
      ENDDO
	    
C END OF CASE IPROFIL = 2
C------------------------------------------------------ 
 7777 CONTINUE



C* Writing of the trace file
C-------------------------------------------
      
      IF (TRACE) THEN      

         WRITE(IDLOG,*,err=921) "****************"
         WRITE(IDLOG,*,err=921) "SOS_PROFIL log :"
         WRITE(IDLOG,*,err=921) ' '
  
         IF(IPROFIL.EQ.1) THEN
            
	     WRITE(IDLOG,*,err=921) "  Profile without gaz"
	     WRITE(IDLOG,*,err=921) "  -------------------"                   
	     WRITE(IDLOG,*,err=921) "    Case : ", CAS_TXT_NG   
	     WRITE(IDLOG,*,err=921) 
     &         "    Number of layers of the profile: ",NT_NG         
	     WRITE(IDLOG,*,err=921)
     &      "    First layer opt. thickness : ",T_FIRST_LAYER_NG        
	     WRITE(IDLOG,*,err=921)
     &      "    Min Opt. thickness of layers below: ",T_LAYER_NG

  	     WRITE(IDLOG,*,err=921) ' '
  	     WRITE(IDLOG,*,err=921)   
     &      "    Total optical thickness of molecules: ",Hmol_NG(NT_NG)
  	     WRITE(IDLOG,*,err=921)   
     &      "    Total optical thickness of aerosols: ",Haer_NG(NT_NG)
  	     WRITE(IDLOG,*,err=921)   
     &      "    Total optical thickness: ",H_NG(NT_NG)
  	     WRITE(IDLOG,*,err=921) ' '
	     WRITE(IDLOG,10,err=921)  
     &             'I','ALT','Hmol(I)','Haer(I)', 
     &             'Habs(I)','H(I)','PCAER(I)','PCMOL(I)'
	     DO I=0,NT_NG          
	        WRITE(IDLOG,16,err=921)  
     &               I,ZPROF_NG(I),Hmol_NG(I),Haer_NG(I),0.D+00, 
     &               H_NG(I),PCAER_NG(I),PCMOL_NG(I)
	     ENDDO

            IF ((ABSPROFIL.NE.7).AND. 
     &          (TABS(CTE_ABS_NBLEV).NE.0.D+00)) THEN

  	        WRITE(IDLOG,*,err=921) ' '
  	        WRITE(IDLOG,*,err=921) ' '
	        WRITE(IDLOG,*,err=921) 
     &         "  Profile with gaz (used for RT calculation)"
	        WRITE(IDLOG,*,err=921) 
     &         "  ------------------------------------------"                 
	        WRITE(IDLOG,*,err=921) "    Case : ", CAS_TXT         
	        WRITE(IDLOG,*,err=921) 
     &         "    Number of layers  of the profile: ",NT       
	        WRITE(IDLOG,*,err=921)
     &         "    First layer opt. thickness : ",T_FIRST_LAYER      
	        WRITE(IDLOG,*,err=921)
     &         "    Min Opt. thickness of layers below: ",T_LAYER

   	        WRITE(IDLOG,*,err=921) ' '
  	        WRITE(IDLOG,*,err=921)   
     &         "    Total optical thickness of molecules: ",Hmol(NT)
  	        WRITE(IDLOG,*,err=921)   
     &         "    Total optical thickness of aerosols: ",Haer(NT)
  	        WRITE(IDLOG,*,err=921)   
     &         "    Total optical thickness of gas absorption: ",  
     &         Habs(NT)
  	        WRITE(IDLOG,*,err=921)   
     &         "    Total optical thickness: ",H(NT)
 
   	        WRITE(IDLOG,*,err=921) ' '
  	        WRITE(IDLOG,*,err=921)   
     &         "    Min altitude (km) for profile definition: ",ZLIM
  	        WRITE(IDLOG,*,err=921)   
     &         "    (If not surface : an additional layer is added"
  	        WRITE(IDLOG,*,err=921)   
     &         "     to complete the profile)"
 
  	        WRITE(IDLOG,*,err=921) ' '
	        WRITE(IDLOG,10,err=921)  
     &             'I','ALT','Hmol(I)','Haer(I)', 
     &             'Habs(I)','H(I)','PCAER(I)','PCMOL(I)'
	        DO I=0,NT          
	           WRITE(IDLOG,16,err=921)  
     &               I,ZPROF(I),Hmol(I),Haer(I),Habs(I), 
     &               H(I),PCAER(I),PCMOL(I)
	        ENDDO
            ENDIF !End of case with gaseous absorption



         ELSE !Case IPROFIL=2

	     WRITE(IDLOG,*,err=921)     
     &      '  Profile for a variable altitude of the aerosol layer'
	     WRITE(IDLOG,*,err=921)     
     &      '  Molecular layer + aerosols : Minimal altitude  : ',ZMIN
	     WRITE(IDLOG,*,err=921)     
     &      '  Molecular layer + aerosols : Maximal altitude : ',ZMAX
            WRITE(IDLOG,*,err=921) ' '
	     WRITE(IDLOG,*,err=921) 
     &	     '  Number of transition layers      : ',NB_TR
	     WRITE(IDLOG,*,err=921) 
     &	     '  Number of sub-layers of the layer 1 : ',NBSC_C1
	     WRITE(IDLOG,*,err=921) 
     &	     '  Number of sub-layers of the layer 2 : ',NBSC_C2
	     WRITE(IDLOG,*,err=921)
     &	     '  Number of sub-layers of the layer 3 : ',NBSC_C3

   	     WRITE(IDLOG,*,err=921) ' '
	     WRITE(IDLOG,10,err=921)'I','ALT','Hmol(I)','Haer(I)',
     &      'Habs(I)','H(I)','PCAER(I)','PCMOL(I)'
	     DO I=0,NT	        
	        WRITE(IDLOG,15,err=921) I,ZPROF(I),Hmol(I),Haer(I),
     &                                  0.D+00,H(I),PCAER(I),PCMOL(I)
	     ENDDO

         ENDIF       
      
         WRITE(IDLOG,*,err=921) ' '
         WRITE(IDLOG,*,err=921) "****************"
	 	   	   	   	   	   
      ENDIF ! Trace





      
C* Recording of the atmospheric profile parameters in a file
C  ---------------------------------------------------------------------
      OPEN(2,FILE=FICPROFIL,ERR=991)
      
      DO I=0,NT
         WRITE(2,20,err=992)I,ZPROF(I),H(I),PCAER(I),PCMOL(I)
      ENDDO

C* File closing
C-------------------
      CLOSE(2)		! Closing of the PROFIL result file

C* Nominal end
C-------------------
      GOTO 9999
      
C* Error case and return of status 1 to the shell
C----------------------------------------------

  940 WRITE(6,*) '  ERROR 940 in SOS_PROFIL :'
      WRITE(6,*) '  IPROFIL must be 1 or 2' 
      WRITE(6,*) '  Current value :', IPROFIL  
      GOTO 9998

  950 WRITE(6,*) '  ERROR 950 in SOS_PROFIL :'
      WRITE(6,*) '  CTE_TOA_FIRST_LAYER_OPT_THICKNESS must be ' 
      WRITE(6,*) '  lower than CTE_TCOUCHE (see inc/SOS.h file)'
      WRITE(6,*) '     CTE_TOA_FIRST_LAYER_OPT_THICKNESS:', 
     &            CTE_TOA_FIRST_LAYER_OPT_THICKNESS    
      WRITE(6,*) '     CTE_TCOUCHE:', CTE_TCOUCHE  
      GOTO 9998

  960 WRITE(6,*) '  ERROR 960 in SOS_PROFIL :'
      WRITE(6,*) 
     &'  CTE_TCOUCHE must be lower than CTE_THRESHOLD_TAUABS'
      WRITE(6,*) '  (see inc/SOS.h file)'
      WRITE(6,*) '     CTE_TCOUCHE:', CTE_TCOUCHE  
      WRITE(6,*) 
     &'     CTE_THRESHOLD_TAUABS:', CTE_THRESHOLD_TAUABS 
      GOTO 9998

  991 WRITE(6,*) '  ERROR on SOS_PROFIL result file opening'
      GOTO 9998

  992 WRITE(6,*) '  ERROR on SOS_PROFIL result file writing'
      GOTO 9998

  912 WRITE(6,*) '  ERROR on logfile opening for SOS_PROFIL'
      GOTO 9998

  921 WRITE(6,*) '  ERROR on logfile writing for SOS_PROFIL'
      GOTO 9998

 1010 WRITE(6,*) '  ERROR on boundary altitudes of the aerosol layer :'
      WRITE(6,*) '  --> check -AP.AerLayer.Zmin and AP.AerLayer.Zmax'
      GOTO 9998

 1020 WRITE(6,*) '  AOT = ',TA, ' is lower than the minimal value',
     &           ' for the profile definition (0.00001)'
      WRITE(6,*) '  ==> Use AER.AOTref=0 if a very low value is wanted'

 9998 IER=-1
      GOTO 9999

                   
C* Format
C---------   
c   10 FORMAT(5X,A1,4X,A3,7X,A7,10X,A7,10X,A7,10X,A4,12X,2(A7,11X))
   10 FORMAT(5X,A1,4X,A3,7X,A7,3X,A7,3X,A7,3X,A4,6X,2(A8,2X))
   15 FORMAT(2X,I4,F11.5,6(F10.5))
   16 FORMAT(2X,I4,F10.5,1X,6(E15.8,2X))	
   20 FORMAT(2X,I5,F10.5,3(E15.8))
   21 FORMAT(2X,I4,6(E15.8,2X))	
   30 FORMAT(5X,A1,2X,A7,2X,A7,2X,A7,2X,A4,5X,2(A7,2X))
   35 FORMAT(2X,I5,6(F10.5))
   
 9999 RETURN
      END	!END OF PROGRAM SOS_PROFIL




C==============================================================================
C PROCEDURE: SOS_DISC
C ==========
C      This procedure estimates the altitude to be assigned to level I of the
C      atmospheric profile.
C
C      A) The optical thickness of level I is estimated from that of level I-1
C         and the optical thickness step.
C
C      B) The corresponding altitude of level I is evaluated by a
C         bisection method, computing the gaseous absorption optical thickness and
C         the molecular and aerosol ones for the tested altitudes Z.
C         The loop is left when the optical thickness at altitude Z
C         is very close to the one wanted at level I.
C
C Description of parameters
C --------------------------
C   molecules:
C      TR  (double)  (I)   Rayleigh optical thickness.
C      HR  (double)  (I)   Scale height of the molecular profile (km).
C
C   Aerosols:
C      TA  (double)  (I)   Aerosol optical thickness.
C      HA  (double)  (I)   Scale height of the aerosol profile (km).
C
C   Gas:
C      TABS(1:CTE_ABS_NBLEV)   (double)  (I)   Initial altitude profile of the absorption optical thickness
C      ALTABS(1:CTE_ABS_NBLEV) (double)  (I)   Altitudes of the initial absorption profile
C      TG_ZLIM (double)  (I) Absorption optical thickness at altitude ZLIM
C      ZLIM (double)  (I) Minimum altitude to consider (in km).
C      ZMAX_INIT  (double)  (I)  Max altitude to consider (in km).
C
C   Information on level I-1:
C      I     (I4)      (I)  Number of the processed level.
C      TIM1  (double)  (I)  Total optical thickness of level I-1.
C      DT    (double)  (I)  Optical thickness step
C
C      ZX    (double)  (O)  Estimated altitude of level I.
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
      SUBROUTINE SOS_DISC(DT,TA,HA,TR,HR,TABS,ALTABS,TIM1,
     &                    ZMAX_INIT, TG_ZLIM, ZLIM, ZX)

      IMPLICIT NONE

C* Definition of variables
C*-----------------------------------------


C*   molecules:
      DOUBLE PRECISION TR	! Rayleigh optical thickness.
      DOUBLE PRECISION HR	! Scale height of the molecular profile.

C*   Aerosols:
      DOUBLE PRECISION TA	! Aerosol optical thickness (not truncated).
      DOUBLE PRECISION HA	! Scale height of the aerosol profile.

C*   Absorbing gases
      DOUBLE PRECISION TABS(1:CTE_ABS_NBLEV)   ! Altitude profile of the absorption optical thickness
      DOUBLE PRECISION ALTABS(1:CTE_ABS_NBLEV) ! Altitudes of the absorption profile
      DOUBLE PRECISION TG      ! Absorption optical thickness
      DOUBLE PRECISION ZLIM    ! Minimum altitude for the bisection search.
      DOUBLE PRECISION TG_ZLIM ! Absorption optical thickness at altitude ZLIM

       
C*   Information on level I-1:
      DOUBLE PRECISION TIM1	! Total optical thickness of level I-1.

      DOUBLE PRECISION ZX	! Estimated altitude of level I.
      	
      DOUBLE PRECISION ZZ	! Linear interpolation coefficient between 2 altitudes
      	
      DOUBLE PRECISION DT	! Optical thickness step.
      DOUBLE PRECISION TI	! Estimated optical thickness at level I.
      DOUBLE PRECISION ZMAX	! Altitude (km).
      DOUBLE PRECISION ZMAX_INIT	! Altitude (km).
      DOUBLE PRECISION ZMIN


      DOUBLE PRECISION ZMOY	! Mean altitude.
      DOUBLE PRECISION TZMOY ! Total optical thickness for the altitude Zmoy.
      DOUBLE PRECISION XD	! Absolute difference.
      
      INTEGER*4 J		 ! Number of the processed level.


C------------------------------------------------------
C* A) Estimation of the optical thickness of level I from
C*    that of level I-1 and the step DT.
C------------------------------------------------------     
      TI=TIM1+DT

C------------------------------------------------------
C* B) Estimation of the altitude corresponding to the optical
C*    thickness TI of level I-1.
C------------------------------------------------------

C*   Initialisation of the altitudes.
C--------------------------------------
      ZMAX=ZMAX_INIT
      ZMIN=ZLIM

C*   Initialisation
C---------------------------
  706 ZMOY=(ZMAX+ZMIN)/2.
     

C*   Computation of the absorption optical thickness at altitude ZMOY
C-----------------------------------------------------------------
      IF (TG_ZLIM.GT.0.0D+00) THEN
         J=2
         DO WHILE(ZMOY.LT.ALTABS(J))
	     J=J+1
         ENDDO
	    
         ! Linear interpolation between the two layers
         IF (ZMOY.GT.ALTABS(1)) THEN !Management of the case where Z is outside the profile (Z > 120 km)
	     ZZ=0
         ELSE
	     ZZ=(ZMOY-ALTABS(J-1))/(ALTABS(J)-ALTABS(J-1))
         ENDIF
         TG=(1-ZZ)*TABS(J-1)+ZZ*TABS(J)
      
      ELSE
         TG=0.0D+00
      ENDIF	  

	
  
C*   Total optical thickness for the altitude Zmoy.
C----------------------------------------------------
      TZMOY=TA*DEXP(-ZMOY/HA)+TR*DEXP(-ZMOY/HR)+TG
 
	        
C*   Absolute difference with the estimated optical thickness of level I.
C-------------------------------------------------------------
      XD=DABS(TI-TZMOY)
  
C*   Exit test if Zmoy corresponds to level I.
C------------------------------------------------------
      IF(XD.LT..000001) GOTO 705
      IF(ZMOY.EQ..00D+00) GOTO 705
      
C*   Refinement of the bounds of the altitude of level I.
C------------------------------------------------------
      IF( (TI-TZMOY).LT.0.D+00 ) THEN
         ZMIN=ZMOY
      ELSE
         ZMAX=ZMOY
      ENDIF
      
      GOTO 706    
       
C*   Assignment of the result
C------------------------------------------------------
 705  ZX=ZMOY
 
 
      RETURN
      END	!END OF PROCEDURE SOS_DISC



