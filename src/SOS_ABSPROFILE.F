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
C* FILE: SOS_ABSPROFILE.F
C* PROJECT: Successive orders of scattering
C* ROLE: Computation of the atmospheric absorption profile for one CKD exponential
C*
C* AUTHOR:
C*   Initial code: Laboratoire d'Optique Atmosphérique (LOA, Villeneuve d'Ascq).
C*   Code rework and development: CS GROUP France, B.Lafrance, X.Lenot
C* 
C*
C* MOD:VERSION:1.0: DATE: 06/06/2013 
C*    File creation to adapt the OS to the CKD method (atmospheric absorption)
C*
C* MOD:VERSION:1.1: DATE: 24/02/2016
C*    - Array NEXP changed from Double Precision to Integer
C*    - Change of the type of variable K (exponential index number):
C*      changed from INTEGER*4 to INTEGER*2
C*      to be consistent with the type in routine coefabs.
C*
C* MOD:VERSION:2.0: DATE: 08/06/2020
C*    - Array XKI is replaced by KDIS_KI with the additional Temperature variable,
C*      as well as the CKD coefficients for the NO2 continuum
C*    - Adjustment of the input arguments (IK8 for NO2).
C*    - Removal of the computations specific to continuums (now included in CKD computations).
C*    - Removal of the output of the total transmission over the atmospheric column.
C*
C* MOD:VERSION:2.1: 31/08/2022
C*    - Constants renamed from SOS_* to CTE_*
C*    - Traces converted to English.
C*
C* MOD:VERSION:2.2: 16/06/2023
C*    - Parameters KDIS_KI_H2O, TAB_CONC_H2O and NB_CONC_H2O added as routine arguments.
C*    - The water vapour concentration is now taken into account.
C*
C******************************************************************************

C----------------------------------------------------------------------------
C Definition of constants
C---------------------------------------------------------------------------- 
C Constants used:
C    CTE_LENFIC2 :   Length of file names with directory tree.
C    CTE_ABS_NBLEV : Number of levels of the absorption profile.
c    CTE_NBABS :     Number of absorbers.
C    CTE_ABS_NBCOL : Number of columns in an absorption profile given
C		       as input.
C----------------------------------------------------------------------------
#include "SOS.h"



C==============================================================================
C PROCEDURE: SOS_ABSPROFILE.F
C ==========
C      This procedure computes the profile of absorption optical thickness
C      versus altitude (on 50 levels).
C
C      The profile is computed for a given wavenumber with the CKD method
C
C      This function is an adaptation of the TRSKD function provided by the LOA.
C      The latter computed the transmission over a given spectral
C      interval.
C
C Description of parameters
C --------------------------
C
C     ABSPROFIL (I2) : (I) Type of atmospheric profile for absorption
C			0:User, 1:TROPICAL, 2:MLS, 3:MLW, 4:SAS, 5:SAW,
C                    6:USST62, 7:No absorption
C
C     NU  (double) : (I)  Wavenumber in cm-1
C
C     LAMB1  (I4) : (I)   Index of the CKD interval
C
C     IABS  (I2) : (I) Array of CTE_NBABS elements. Each element corresponds to an
C			  absorbing gas (Absorbers: H2O CO2 O3 N2O CO CH4 O2 NO2).
C	 		  For each absorber, 1 if processed with the CKD method, 0 otherwise
C
C     USERPROFIL (double) (I)  Array of (CTE_ABS_NBLEV * CTE_ABS_NBCOL) elements
C			    containing for each absorbing gas the profile expressed in ppmv*1E+6
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
C     ALTABS (double) : (I)  Array of (CTE_ABS_NBLEV) containing the altitude levels of the profile
C
C     RO (double) (I) : Array of (CTE_NBABS * CTE_ABS_NBLEV) elements
C		          containing, for each absorbing gas of the CKD method, the profile
C			   expressed in particles/cm2 for each layer.
C			   It contains, in this order: H2O, CO2, 03, N2O, CO, CH4, O2, NO2
C                       --> Gases with fine line absorption (not only continuums)
C
C     NEXP (I4) :   (I)  Array of (CTE_NBABS * CTE_CKD_NWVL_MAX) elements containing, for each absorber
C		           and each spectral interval, the number of exponentials of the CKD method
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
C     IK1  (I4) : (I) Index of the exponential for H2O absorption
C
C     IK2  (I4) : (I) Index of the exponential for CO2 absorption
C
C     IK3  (I4) : (I) Index of the exponential for O3 absorption
C
C     IK4  (I4) : (I) Index of the exponential for N2O absorption
C
C     IK5  (I4) : (I) Index of the exponential for CO absorption
C
C     IK6  (I4) : (I) Index of the exponential for CH4 absorption
C
C     IK7  (I4) : (I) Index of the exponential for O2 absorption
C
C     IK8  (I4) : (I) Index of the exponential for NO2 absorption
C
C     TAB_PRES (double): (I) Array of the pressure values associated with the coefficients ki(P,T)
C                            of CTE_CKD_NP_MAX elements.
C     NB_PRES (integer): (I) Number of actual pressure values
C
C     TAB_TEMP (double): (I) Array of the temperature values associated with the coefficients ki(P,T)
C                            of CTE_CKD_NT_MAX elements.
C     NB_TEMP (integer): (I) Number of actual temperature values
C
C     TAB_CONC_H2O (double): (I) Array of the H2O concentration values associated with the coefficients ki(P,T,C)
C                                of water vapour, of CTE_CKD_NC_MAX elements (ppmv)
C     NB_CONC_H2O (integer): (I) Number of actual concentration values.
C
C     TAUABSTOT (double) : (O) Array of (CTE_ABS_NBLEV-1) containing the absorption optical
C			          thicknesses (continuum + CKD) for each atmospheric layer
C
C     TRACE (LOGICAL) : (I) If true => trace files are written
C
C     IDLOG (I4) : (I) Identifier of the log file
C
C     IER (I4) : (O) error code =0 if no error, =-1 otherwise
C  
C
C Common used:
C --------------
C     None
C
C Error cases:
C ------------
C     Writing of the log
C
C==============================================================================

      SUBROUTINE SOS_ABSPROFILE(ABSPROFIL, NU, LAMB1, IABS, 
     &		         USERPROFIL, ALTABS, RO, NEXP,
     &		         KDIS_KI, KDIS_KI_H2O,
     &		         IK1, IK2, IK3, IK4, IK5, IK6, IK7, IK8,
     &                 TAB_PRES, NB_PRES, TAB_TEMP, NB_TEMP,
     &                 TAB_CONC_H2O, NB_CONC_H2O,
     &			  TAUABSTOT, TRACE, IDLOG, IER) 
    
      IMPLICIT NONE

C* Definition of parameters
     
      INTEGER*4  NLEVEL    ! nb of altitude levels of the profile
      PARAMETER(NLEVEL=CTE_ABS_NBLEV)
      
      INTEGER    MAXABS   ! Number of absorbing gases (CKD method)
      PARAMETER (MAXABS=CTE_NBABS)

      INTEGER    MAXKD   ! Max number of exponentials per absorbing gas (CKD method)
      PARAMETER (MAXKD=CTE_CKD_NAI_MAX)

      INTEGER    MAXT    ! Max number of atmospheric temperature levels (CKD method)
      PARAMETER (MAXT=CTE_CKD_NT_MAX)
            
      INTEGER    MAXP    ! Max number of atmospheric pressure levels (CKD method)
      PARAMETER (MAXP=CTE_CKD_NP_MAX)
      
      INTEGER    MAXNU    ! max number of wavenumbers for the absorption computation
      PARAMETER (MAXNU=CTE_CKD_NWVL_MAX)
  
      integer ij !for test
C* Definition of variables
C--------------------------------------------

C*   Simu
      DOUBLE PRECISION NU    ! Wavenumber (cm-1)

C*   Profiles
      INTEGER*2 ABSPROFIL    ! (I) Type of the absorption profile
      CHARACTER*CTE_LENFIC2 FICABSPROFIL ! (I) Absorption profile file provided by the user
      
      DOUBLE PRECISION ALTABS(1:NLEVEL)         ! Altitudes of the absorption profile (decreasing order)
      DOUBLE PRECISION TAUABSTOT(1:NLEVEL)      ! Optical thickness profile
      DOUBLE PRECISION TRANS(1:NLEVEL)          ! Transmission profile
      
      DOUBLE PRECISION USERPROFIL(NLEVEL,CTE_ABS_NBCOL) ! absorption profile given as input (ppmv)
      DOUBLE PRECISION RO(CTE_NBABS,NLEVEL)	! Number of particles/cm2 per gas species and per layer of the atmospheric profile
      DOUBLE PRECISION TAUABS1C(NLEVEL-1)        ! Optical thickness of each layer


C*   CKD coefficients
      INTEGER*4  NEXP(CTE_NBABS,CTE_CKD_NWVL_MAX) 
                              ! (I) Number of exponentials per wavelength and species
                              
      
      DOUBLE PRECISION  KDIS_KI(CTE_CKD_NT_MAX, CTE_CKD_NP_MAX, 
     &                          CTE_CKD_NAI_MAX, CTE_NBABS,  
     &                          CTE_CKD_NWVL_MAX )
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

      INTEGER NB_TEMP, NB_PRES  ! (I) Number of temperature and pressure values in the CKD file
      DOUBLE PRECISION TAB_TEMP(CTE_CKD_NT_MAX) ! (I) Array of the temperature values read
      DOUBLE PRECISION TAB_PRES(CTE_CKD_NP_MAX) ! (I) Array of the pressure values read

      INTEGER NB_CONC_H2O  ! (I) Number of H2O concentration values in the H2O CKD file
      DOUBLE PRECISION TAB_CONC_H2O(CTE_CKD_NC_MAX) ! (I) Array of the H2O concentration values read
           
                                         
      DOUBLE PRECISION  PRS     ! Mean pressure of a layer of the profile (mbar)
      DOUBLE PRECISION  TMP     ! Mean temperature of a layer of the profile (K)
      DOUBLE PRECISION  CONC_H2O! Mean H2O concentration of a layer of the profile (in ppmv)
      DOUBLE PRECISION  XK

c*   Transmissions
      DOUBLE PRECISION TRS
      
C*   Gaseous absorption
      INTEGER*2 IABS(1:CTE_NBABS)   ! Array identifying the simulated absorption bands
                                    ! (H2O, CO2, O3, N2O, CO, CH4, O2, NO2)
                                    ! Includes absorption lines and/or continuums.

      INTEGER*2 NABS        ! Index of the absorber
      INTEGER*2 JABS		! Status of the absorber (1 if processed, 0 otherwise)
       
C*    Other variables
      INTEGER*4 I,J  					! Level.
      INTEGER*4 IK,IK1,IK2,IK3,IK4,IK5,IK6,IK7,IK8  	! Exponential index
      INTEGER*2 K	

      INTEGER*4 LAMB1

      INTEGER*4 IER     	! error code =0 if no error, =-1 otherwise

      LOGICAL TRACE        	! = true if writing in the trace file

      INTEGER*4 IDLOG	! Identifier number of the trace file
      



C* Initialisation
C------------------------------------------------------

C* Error code
      IER=0

C* Initialisation of arrays
      DO J=1,NLEVEL
         TAUABSTOT(J)=0.D+00
      ENDDO
        
C* No processing if ABSPROFIL=7
C----------------------------------

      IF (ABSPROFIL.EQ.7) RETURN
	
C* Computation of absorption optical thicknesses for an altitude profile
C* CKD method
C----------------------------------------------------------------------

     
c The transmission at the upper level is 1
      TRANS(1)=1.D+00	
      TAUABSTOT(1)=0.D+00
            	   
c Computation of the transmission + optical thickness for each layer of the profile
c There are N levels, hence N-1 layers
       !Init of the transmission
       TRS=1.0
       DO J=1,NLEVEL-1	!Loop on the layers from TOA to ground
          TAUABS1C(J)=0.		! Absorption optical thickness of layer J (profile: top layer=1, bottom layer=NLEVEL-1)
          PRS=(USERPROFIL(NLEVEL-J,2)+USERPROFIL(NLEVEL-J+1,2))/2.  !USERPROFIL : top layer=NLEVEL-1, bottom layer=NLEVEL-1
	   TMP=(USERPROFIL(NLEVEL-J,3)+USERPROFIL(NLEVEL-J+1,3))/2.
          CONC_H2O=(USERPROFIL(NLEVEL-J,4)+USERPROFIL(NLEVEL-J+1,4))/2.
          CONC_H2O=CONC_H2O*1.D-06 !Conversion to ppmv

          DO K=1,CTE_NBABS
             IF(K.EQ.1)IK=IK1
             IF(K.EQ.2)IK=IK2
             IF(K.EQ.3)IK=IK3
             IF(K.EQ.4)IK=IK4
             IF(K.EQ.5)IK=IK5
             IF(K.EQ.6)IK=IK6
             IF(K.EQ.7)IK=IK7
             IF(K.EQ.8)IK=IK8
             XK=0.
             
             IF (NEXP(K,LAMB1).GE.1) THEN       ! If nb of exponentials = 1: continuum only (no fine lines)
                
     		  CALL COEFF_ABS_CKD(K,LAMB1,IK,PRS,TMP,CONC_H2O,
     &                             TAB_PRES,NB_PRES,TAB_TEMP,NB_TEMP,
     &                             TAB_CONC_H2O, NB_CONC_H2O,
     &                             KDIS_KI,KDIS_KI_H2O,XK,IER)
		
                IF(IER.NE.0) GOTO 910
             ENDIF

             !Total absorption optical thickness (lines + continuums)
             !for atmospheric layer J and component K
             TAUABS1C(J)=TAUABS1C(J)+XK*RO(K,NLEVEL-J)

          ENDDO !End of loop on gas species for level J
  
c computation of the transmission between the top and level nLevel - J
          TRS=TRS*exp(-TAUABS1C(J))
	   TRANS(J+1)= TRS		!Transmission at bottom of layer J = level J+1


c computation of the optical thickness profile (bottom of layer)
          IF (TRS.GT.0.0D+00) THEN
	      TAUABSTOT(J+1)=-DLOG(TRS)
          ELSE
             !Undefined extreme optical thickness: set to CTE_TAUABS_MAX
	      TAUABSTOT(J+1)=CTE_TAUABS_MAX
          ENDIF
       ENDDO !End of loop on layers

       

      
C* Writing of the trace file
C-------------------------------------------
      IF (TRACE) THEN      

	     WRITE(IDLOG,*,ERR=921) ' '	 
	     WRITE(IDLOG,*,ERR=921) ' ***** Absorption profile ****'
	     WRITE(IDLOG,*,ERR=921) ' IK1 (H2O): ',IK1 
	     WRITE(IDLOG,*,ERR=921) ' IK2 (CO2): ',IK2
	     WRITE(IDLOG,*,ERR=921) ' IK3  (O3): ',IK3
	     WRITE(IDLOG,*,ERR=921) ' IK4 (N2O): ',IK4
	     WRITE(IDLOG,*,ERR=921) ' IK5  (CO): ',IK5
	     WRITE(IDLOG,*,ERR=921) ' IK6 (CH4): ',IK6
	     WRITE(IDLOG,*,ERR=921) ' IK7  (O2): ',IK7
	     WRITE(IDLOG,*,ERR=921) ' IK8 (NO2): ',IK8
	     WRITE(IDLOG,*)'LEVEL ALT(km)  OT   TRANS'
	     DO j=1,NLEVEL
               WRITE(IDLOG,18)j,ALTABS(j),TAUABSTOT(j),TRANS(j)
            ENDDO
	     WRITE(IDLOG,*,ERR=921) ' '

      ENDIF
 


C* Nominal end
C-------------------
      GOTO 9999
      
C* Error case and return of status 1 to the shell
C----------------------------------------------

  910 WRITE(6,*) '  SOS_ABSPROFILE : ERROR_910 : '
      WRITE(6,*) '  --> Error in subroutine COEFF_ABS_CKD'
      GOTO 9998
      
  921 WRITE(6,*) '  ERROR writing logfile for SOS_PREPAPROFIL'
      GOTO 9998

 9998 IER=-1
      GOTO 9999
      
C* Format
C---------   
  
!!   18 FORMAT(i2.2,1x,f7.2,2f9.5)   
   18 FORMAT(i2.2,1x,f7.2,2e15.5)   
   
 9999 RETURN
	 
      END
