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
C* FILE: SOS_PROC.F
C* PROJECT: Successive orders of scattering
C* ROLE: 
C*      This procedure checks the simulation parameters.
C*      It runs successively the different modules of the code:
C*         - Definition of the angles used: call of SOS_ANGLES
C*         - Computation of the aerosol radiative properties: call of SOS_AEROSOLS
C*         - Estimation of the AOT at the reference wavelength
C*         - Computation of the surface radiative properties: call of SOS_SURFACE
C*         - Preparation of the data of the gaseous absorption profile: call of SOS_PREPA_ABSPROFILE
C*         - Preparation of the radiative data for the OS computation: call of SOS_PREPA_OS
C*
C*         For the fine CKD mode:
C*              For each exponential of the CKD decomposition:
C*              - Generation of the gaseous absorption optical thickness profile for one exponential:
C*                call of SOS_ABSPROFILE
C*              - Generation of the optical thickness profile (scattering + absorption): call of SOS_PROFILE
C*              - Computation of the radiance field: call of SOS
C*              - Aggregation of the radiances (and transmissions) according to the CKD decomposition: call of SOS_AGGREGATE
C*
C*         For the simplified CKD mode: Computation of the absorption optical thickness by CKD
C*              - Generation of the gaseous absorption optical thickness profile: call of SOS_ABSPROFILE
C*              - Generation of the optical thickness profile (scattering + absorption): call of SOS_PROFILE
C*              - Computation of the radiance field: call of SOS
C*
C*         Recomposition of the radiance from the results of the Fourier series decomposition in azimuth:
C*         call of SOS_TRPHI_OPTION
C*
C*      Procedure that can be called from a Python program (Python binding).
C*
C* AUTHORS:
C*   CS GROUP France, X.Lenot, T.Erudel, B. Lafrance
C*
C* DATE: 23/11/2023
C*
C* MOD:VERSION:1.0: File creation
C*
C* MOD:VERSION:1.1: 05/04/2023
C*    - Introduction of a check of the sum of the CKD coefficients AIK and normalisation to obtain exactly 1.
C*   
C*  MOD:VERSION:1.2: 19/04/2023
C*    - Removal of the variable IND_ANGOUT in the call of SOS_FIN_OS
C*    - Addition of the variable IND_ANGOUT_FIN to the output parameters of SOS_PROC
C*   
C*  MOD:VERSION:1.3: 10/08/2023
C*    - Retrieval of the truncation coefficient of the aerosol phase function as output of routine SOS_AEROSOLS
C*    - Addition of the truncation coefficient of the aerosol phase function as output parameter
C*      (Python binding and routine to the Fortran main program)
C*
C*  MOD:VERSION:1.4: 15/09/2023
C*    - Correction of the index loops for the initialisation to 1 of the table KDIS_AI(I,J,K)
C*    - Removal of a comma after the command WRITE(...) in one place (log file output).
C*      No impact on the simulation results, but raised a warning at compilation.
C*    - Removal of the call of routine SOS_FIN_OS in favour of a direct call of SOS_TRPHI_OPTION
C*
C*  MOD:VERSION:1.5: 14/11/2023
C*    - Removal of the tests molecular and aerosol optical thicknesses > 0.0001 (0. instead) to check the given values HR and HA.
C*    - Correct naming of the formula of the molecular optical thickness (CNES, Perbos, 1982) instead of Hansen & Travis
C*
C*  MOD:VERSION:1.6: 09/11/2026 (T. Harmel)
C*    - Computation of the direct transmission AT_DIR of the solar beam
C*      (TOA -> surface -> output level) given to SOS_TRPHI_OPTION for the
C*      direct surface reflection terms (glitter, Fresnel, BRDF/BPDF).
C*      With -SOS.AbsModeCKD 1, AT_DIR is the CKD-weighted sum of the
C*      transmissions of each k-term, instead of the Beer law applied to
C*      the aggregated optical thickness (exact only for an air mass of 1).
C*    - Same CKD weighting for the direct transmissions used in the fluxes
C*      (TDIR_VRAI, FLUX_DIFF_DOWN, FLUX_DOWN) and in the -SOS.Trans file.
C*
C----------------------------------------------------------------------------
C Definition of constants
C---------------------------------------------------------------------------- 
C Constants used:
C    CTE_LENDIR : Length of directory names.
C     CTE_LENFIC1 : Length of file names without directory tree.
C     CTE_LENFIC2 : Length of file names with directory tree.
C    CTE_LENCOM : Max length of a command line
C    CTE_WAMIN: Minimum wavelength (in microns)
C    CTE_WAMAX: Maximum wavelength (in microns)
C    CTE_MDF : Molecular depolarisation factor
C    CTE_CKD_NWVL_MAX: Max number of wavenumbers for the absorption computation
C    CTE_CKD_NAI_MAX: Max number of exponentials per absorbing gas (CKD method)
C    CTE_CKD_NT_MAX   : Maximum number of temperatures of the CKD data
C    CTE_CKD_NP_MAX   : Maximum number of pressures of the CKD data
C    CTE_HT_STD_PSURF: Standard atmospheric pressure
C    CTE_NBABS: Number of absorbing gases
C    CTE_ABS_NBLEV : Number of layers in the profile of absorbing gases
C    CTE_ABS_NBCOL : Number of columns of the absorbing gas profile files
C    CTE_NOT_DEFINED_VALUE_INT: default value for integers
C    CTE_NOT_DEFINED_VALUE_DBLE: default value for reals
C    CTE_DEFAULT_FICGRANU: default name of the size distribution file
C    CTE_REMOVE_PREVIOUS_BIN : Option to remove or not a result file
C    CTE_OS_NBMU_MAX  : Maximum number of positive Gauss angles for the radiance computations
C    CTE_OS_NB_MAX :    Maximum value of the limit order of the expansion in Legendre polynomials
C    CTE_TOA_ALT : Altitude of the TOA level (km).
C    CTE_DEFAULT_IGMAX : Default maximum number of scattering orders
C
C----------------------------------------------------------------------------
#include "SOS.h"
#define INCTE_PI DACOS(-1.D+00)
#define INCTE_IDLOG_MAIN 77
#define INCTE_IDLOG_PROFIL 88
#define INCTE_IDLOG_OS 99     
#define INCTE_NB_Modeles 4
#define INCTE_NB_PROFIL 2


C
C Description of parameters
C ---------------------------
C
C   Inputs:
C			   
C     RESROOT (CHARACTER) : (I) Working and result output directory
C                                               (complete path)
C
C     FICMAIN_LOG (CHARACTER) : (I) Name of the trace file of the main program.
C                                   (file name + extension)
C
C     WA_SIMU (DOUBLE) : (I) Wavelength of the radiance simulation in microns.
C
C            NBMU_GAUSS_LUM (I4) : (I) Number of Gauss angles for the radiance computations.
C
C     FICANGLES_USER_LUM (CHARACTER) : (I) complete name of the "user" angles file for the radiance computations.
C				  (directory + file name + extension)
C
C 	     TETAS   (DOUBLE) : (I) Solar zenith angle in degrees (0 < tetas < 90).
C      		      (complementary angle: mus= -cos(tetas))
C
C     FICANGLES_RES_LUM (CHARACTER) : (I) Angle file for the radiance computations (in RESROOT/SOS)
C
C            NBMU_GAUSS_MIE (I4) : (I) Number of Gauss angles for the phase function computations
C
C     FICANGLES_USER_MIE (CHARACTER) : (I) Complete name of the "user" angles file for the Mie computations.
C				  (directory + file name + extension)
C
C     FICANGLES_RES_MIE : (I) Angle file for the Mie computations (in RESROOT/SOS)
C
C     FICANGLOG (I) : Complete name of the trace file (directory + file name + extension)

C     WAREF_AOT (DOUBLE) : (I) Reference wavelength for the AOT.
C
C     AOT_REF (DOUBLE) : (I) Aerosol optical thickness at the reference wavelength (not truncated).
C
C     ITRONC_AER (I4) : (I) Truncation option of the aerosol phase function
C
C     FICGRANU_LOG (CHARACTER) : (I) Name of the TRACE file of the aerosol IOP computations (in RESROOT/LOG)
C
C     FICMIE_LOG (CHARACTER) : (I) Name of the TRACE file of the Mie computations (in RESROOT/LOG)
C
C     DIR_MIE (CHARACTER) : (I) Output directory of the Mie computations
C
C     FICGRANU (CHARACTER) : (I) Name of the result file of the aerosol IOP computations (in RESROOT/SOS)
C
C     IMOD_AER (I2) : (I) Index of the mono- or multi-modal size distribution type
C			        0 : Mono-modal models.
C			        1 : WMO models.
C			        2 : Shettle & Fenn models.
C		               3 : LND bi-modal models.
C                            4 : Use of external phase functions
C
C     RN_WA, IN_WA, RN_WAREF, IN_WAREF (DOUBLE) : (I) Refractive index of the aerosols in the mono-modal case
C                    Real part and imaginary part
C                    For the radiance simulation wavelength and for the AOT reference wavelength
C
C     IGRANU (I2) : (I) Index of the mono-modal model
C			        1 : LND
C			        2 : Junge Law
C
C     LND_RADIUS_MMD_AER, LND_LNVAR_MMD_AER (DOUBLE) : (I) Modal radius and variance of the LND law of a mono-modal case.
C
C     JD_SLOPE_MMD_AER(DOUBLE) : (I) Slope of the Junge law of a mono-modal case
C
C     JD_RMIN_MMD_AER, JD_RMAX_MMD_AER (DOUBLE) : (I) Minimum and maximum radius of the Junge law.
C
C     IModele_WMO (I2) : (I) Type of the WMO model
C      		        1 : Continental WMO model.
C			        2 : Maritime WMO model.
C			        3 : Urban WMO model.
C		               4 : WMO model by user definition.
C
C     C_WMO_DL, C_WMO_WS, C_WMO_OC, C_WMO_SO (DOUBLE) : (I) User values of the user WMO mixture
C                 Proportion of Dust-Like, Water Soluble, Oceanic and Soot particles
C
C     IModele_SF (I2) : (I) Type of the Shettle & Fenn model
C      		1 : Tropospheric S&F model.
C			2 : Urban S&F model.
C			3 : Maritime S&F model.
C		       4 : Coastal S&F model.
C
C     RH (DOUBLE) : (I) Relative humidity (%) of the Shettle & Fenn model
C
C     MODE_PARAM_BILND (I2) : (I) Index of the choice of description of a mixture of 2 LND.
C      		1 : Use of user volumetric concentrations.'
C      		2 : Use of user ratio AOT_coarse / AOT_total.'
C     
C     USER_CV_COARSE, USER_CV_FINE (DOUBLE) : (I) User value of the volume concentration
C                                                 of the "LND coarse" and "LND fine"
C
C     RTAUct_WAREF (DOUBLE) : (I) Ratio of the optical thicknesses AOT_coarse / AOT_total at the reference wavelength
C
C     BMD_CM_MRWA, BMD_CM_MIWA, BMD_CM_MRWAREF, BMD_CM_MIWAREF (DOUBLE) : (I)
C                    Refractive index of the Coarse mode of a bimodal aerosol model
C                    Real part and imaginary part
C                    For the radiance simulation wavelength and for the AOT reference wavelength
C
C     BMD_CM_RMODAL, BMD_CM_VAR (DOUBLE) : (I) Modal radius and variance of the LND law
C                                              for the Coarse mode of a bimodal aerosol model.

C     BMD_FM_MRWA, BMD_FM_MIWA, BMD_FM_MRWAREF, BMD_FM_MIWAREF(DOUBLE) : (I)
C                    Refractive index of the Fine mode of a bimodal aerosol model
C                    Real part and imaginary part
C                    For the radiance simulation wavelength and for the AOT reference wavelength

C     BMD_FM_RMODAL, BMD_FM_VAR (DOUBLE) : (I) Modal radius and variance of the LND law
C                                              for the Fine mode of a bimodal aerosol model.
C
C     FICEXTDATA_AER (CHARACTER) : (I) Name of the user file giving the aerosol phase function
C                                      and the scattering and extinction coefficients.
C                 			    (complete path: directory + file name)
C
C     FICMIXTURE_AER (CHARACTER) : (I) Name of the file defining a mixture of aerosol models by the user
C                 			    (complete path: directory + file name)
C
C     FICUSER_AER (CHARACTER) : (I) Name of the aerosol file provided by the user
C                 			    (complete path: directory + file name)
C
C     FICPROFIL_LOG (CHARACTER) : (I) Name of the trace file of the profiles per CKD order (in RESROOT/LOG)
C
C     TR (DOUBLE) : (I) Molecular optical thickness at the radiance simulation wavelength
C
C     HR, HA (DOUBLE) : (I) Scale height of the molecules and of the aerosols respectively (km)
C
C     IPROFIL (I2) : (I) Type of the aerosol profile
C      		1 : Aerosol profile described by an exponential decrease using a heigth scale
C			2 : Homogeneous aerosol layer between two altitudes.
C
C     ZMIN, ZMAX (DOUBLE) : (I) Minimum and maximum altitudes of the molecules + aerosols layer (km)
C                               in the case of an aerosol mixture between two altitudes.
C
C     PSURF (DOUBLE) : (I) Surface pressure (mbar)
C
C     H2O (DOUBLE) : (I) Integrated water vapour content (g/cm2)
C
C     O3 (DOUBLE) : (I) Integrated ozone content (Dobson)
C
C     CO2 (DOUBLE) : (I) Carbon dioxide concentration at the surface (ppmv)
C
C     CH4 (DOUBLE) : (I) Methane concentration at the surface (ppmv)
C
C     ABSPROFIL (I2) : (I) Type of the absorption profile
C			      0:User
C			      1:TROPICAL
C			      2:MLS
C			      3:MLW
C			      4:SAS
C			      5:SAW
C			      6:USST62
C			      7:No absorption
C
C     FICABSPROFIL (CHARACTER) : (I) Absorption profile file provided by the user (user profile case)
C                 			    (complete path: directory + file name)
C
C     NUSTEP (DOUBLE) : (I) Wavenumber step (cm-1) in the CKD coeff file
C                           Must be 1, 5 or 10 cm-1
C
C     ISURF (I2) : (I) Surface type
C      		0 : lambertian surface for an albedo RHO.
C			1 : lambertian surface + glitter (roughness water).
C			2 : lambertian surface + Fresnel's reflexion (plane water).
C		       3 : lambertian surface + Roujean's BRDF.
C			4 : lambertian surface + Roujean's BRDF + Rondeaux's BPDF.
C			5 : lambertian surface + Roujean's BRDF + Breon's BPDF .
C			6 : lambertian surface + Roujean's BRDF + Nadal's BPDF.
C
C     DIR_SURF (CHARACTER) : (I) Storage directory of the BRDF/BPDF files of the surface models
C                 		     (complete path of the directory)
C
C     FICSURF_LOG (CHARACTER) : (I) Name of the trace file of the BRDF/BPDF computations (in RESROOT/LOG)
C
C     SURF_IND (DOUBLE) : (I) Refractive index surface / air.
C
C     WIND (DOUBLE) : (I) Wind speed at the surface (m/s)

C     K0_ROUJEAN, K1_ROUJEAN, K2_ROUJEAN (DOUBLE) : (I) Parameters of the Roujean BRDF model
C
C     ALPHA_NADAL,  BETA_NADAL (DOUBLE) : (I) Parameters of the Nadal BPDF model
C
C     COEF_C_MAIGNAN (DOUBLE) : (I) Parameter C*exp(-NDVI) of the Maignan BPDF model
C
C             RHO (DOUBLE) : (I) Surface albedo
C
C     FICSURF (CHARACTER) : (I) Name of a reflection matrix file provided by the user
C                 		    (complete path: directory + file name)
C
C     FICSOS_LOG (CHARACTER) : (I) Name of the trace file of the RTE solving by SOS_OS (in RESROOT/LOG)
C
C     FICSOS_RES_BIN (CHARACTER) : (I) Name of the result file of the successive orders computations (in RESROOT/SOS)
C                                      for a unit processing (i.e. one CKD exponential if gaseous absorption is computed)
C
C
C   Outputs:
C
C     FICTRANS  (CHARACTER) : (O) Name of the optional output file of the transmissions (in RESROOT/SOS)
C
C     FICFLUX  (CHARACTER) : (O) Name of the optional output file of the fluxes up TOA and down BOA (in RESROOT/SOS)
C
C     ZOUT (DOUBLE) : (I) Altitude of the output level of the radiance field
C			          -1 =>  Standard output: Radiance emerging at TOA + downward radiance at ground level.
C			      z (km) =>  Upward and downward radiance for the altitude z (0 <= z <= CTE_TOA_ALT).
C
C     IGMAX (I4) : (I) Maximum order of interaction
C
C             IPOLAR (I4) : (I) Option to cut the polarisation.
C                               0 : simulation without polarisation.
C			           1 : normal simulation with polarisation.
C
C     IDIRSURF (I4) : (I) Option for the direct reflection of the solar
C                         beam by the surface (sun glint, Fresnel flat sea,
C                         BRDF/BPDF), added analytically to the diffuse field:
C                           0 : diffuse field only (direct surface terms to
C                               be added separately, as in OSOAA LUTs).
C                           1 : complete field (standard).
C
C     ITRPHI (I4) : (I) Option of output type:
C                          1 : Viewing plane for one azimuth value.
C                          2 : Polar diagram versus the azimuth and viewing directions.
C
C     PHIOS (DOUBLE) : (I) Relative azimuth in degrees for the case of an output in one azimuth plane.
C
C     PAS_PHI (I4) : (I) Azimuth step in degrees (integer value) in the case of a polar diagram.
C
C     IMODE_CKD_CALCUL (I4) : (I) Option of the mode of application of the CKD decomposition
C          1 : Finest computation with one radiance computation per order of the decomposition
C              then recomposition with the CKD weighting.
C          2 : Computation for an absorption optical thickness estimated by the CKD approach
C              --> A single radiance computation
C
C     TRACE (BOOLEAN) : (I) True if a trace file of the main program is generated
C
C
C     IER (I4) : (O) Error return code (1 if error, 0 otherwise)
C
C            LUM_NBMU (I4) : (O) Number of actual angles for the radiance computations
C   
C            IND_ANGOUT_FIN (DOUBLE) : (O) Indices of the angle arrays for which there must be
C                                      an output in the result file of the user angles
C      				    (1: output required)
C
C     PHI_FIN (DOUBLE) : (O) Table of the relative azimuth values of the output field
C  
C     THETA_FIN (DOUBLE) : (O) Table of the viewing zenith angle values of the output field
C  
C     SCA_UP_FIN, SCA_DOWN_FIN (DOUBLE) : (O) Tables of the scattering angles for the output upward and downward fields
C
C     I_UP_FIN, Q_UP_FIN, U_UP_FIN (DOUBLE) : (O) Tables of the I, Q and U values for the output upward fields
C  
C     I_DOWN_FIN, Q_DOWN_FIN, U_DOWN_FINC (DOUBLE) : (O) Tables of the I, Q and U values for the output downward fields
C   
C     POL_ANG_UP_FIN, POL_ANG_DOWN_FIN (DOUBLE) : (O) Tables of the polarisation angles
C                                                     for the output upward and downward fields
C
C     POL_RATE_UP_FIN, POL_RATE_DOWN_FIN (DOUBLE) : (O) Tables of the degrees of polarisation
C                                                       for the output upward and downward fields
C
C     L_POL_UP_FIN, L_POL_DOWN_FIN (DOUBLE) : (O) Table of the polarised radiance values
C                                                 for the output upward and downward fields
C     
C     TDIR_VRAI (DOUBLE) : (O) Direct transmission for the solar incidence direction of the true atmosphere
C
C     FLUX_DIFF_DOWN (DOUBLE) : (O) Downward diffuse flux illuminating the surface, normalised by the incident solar flux at TOA
C      
C     FLUX_DOWN (DOUBLE) : (O) Total downward flux at the surface level, normalised by the incident solar flux at TOA
C     
C     EPLUS (DOUBLE) : (O) Diffuse flux emerging at TOA, normalised by the incident solar flux at TOA
C
C 
C Description of the files used
C ---------------------------------
C
C     Generation of the main trace file in text format.
C
C     File FICTRANS:
C        --> Gives the direct and diffuse transmission from TOA to the ground for the solar incidence direction.
C            Also gives the upward diffuse transmission from the surface to TOA for the different viewing angles.
C
C     File FICFLUX:
C        --> Gives the direct and diffuse downward flux at the surface level.
C            Also gives the upward diffuse flux at TOA.
C               NB: the flux values are normalised by the incident solar irradiance at TOA (Esun * cos SZA).
C            The atmospheric profile is recalled:
C                  Molecular optical thickness (MOT),
C                  Aerosol optical thickness (AOT),
C                  Gaseous absorption optical thickness (GOT).
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C
C    Various consistency checks of the parameters
C    
C    Check of the error return of a called routine:
C          SOS_ANGLES, SOS_AEROSOLS, SOS_SURFACE
C          SOS_PREPA_ABSPROFILE, SOS_PREPA_OS, SOS_ABSPROFILE, SOS_PROFILE
C          SOS, SOS_AGGREGATE and SOS_TRPHI_OPTION
C
C--------------------------------------------------------------------------------------------------------------------
      SUBROUTINE SOS_PROC(RESROOT,
     &     FICMAIN_LOG,
     &     WA_SIMU,
     &     NBMU_GAUSS_LUM,
     &     FICANGLES_USER_LUM,
     &     TETAS,
     &     FICANGLES_RES_LUM,
     &     NBMU_GAUSS_MIE,
     &     FICANGLES_USER_MIE,
     &     FICANGLES_RES_MIE,
     &     FICANGLOG,
     &     WAREF_AOT, AOT_REF,
     &     ITRONC_AER,
     &     FICGRANU_LOG, FICMIE_LOG,
     &     DIR_MIE,
     &     FICGRANU,
     &     IMOD_AER,
     &     RN_WA, IN_WA, RN_WAREF, IN_WAREF,
     &     IGRANU,
     &     LND_RADIUS_MMD_AER, LND_LNVAR_MMD_AER,
     &     JD_SLOPE_MMD_AER, JD_RMIN_MMD_AER, JD_RMAX_MMD_AER,
     &     IModele_WMO, C_WMO_DL, C_WMO_WS, C_WMO_OC, C_WMO_SO,
     &     IModele_SF, RH,
     &     MODE_PARAM_BILND, USER_CV_COARSE, USER_CV_FINE,
     &     RTAUct_WAREF,
     &     BMD_CM_MRWA, BMD_CM_MIWA, BMD_CM_MRWAREF, BMD_CM_MIWAREF,
     &     BMD_CM_RMODAL, BMD_CM_VAR,
     &     BMD_FM_MRWA, BMD_FM_MIWA, BMD_FM_MRWAREF, BMD_FM_MIWAREF,
     &     BMD_FM_RMODAL, BMD_FM_VAR,
     &     FICEXTDATA_AER,
     &     FICMIXTURE_AER,
     &     FICUSER_AER,
     &     FICPROFIL_LOG,
     &     TR, HR, HA,
     &     IPROFIL,
     &     ZMIN, ZMAX,
     &     PSURF, H2O, O3, CO2, CH4,
     &     ABSPROFIL, FICABSPROFIL,
     &     NUSTEP,
     &     ISURF,
     &     DIR_SURF,
     &     FICSURF_LOG,
     &     SURF_IND,
     &     WIND,
     &     K0_ROUJEAN, K1_ROUJEAN, K2_ROUJEAN,
     &     ALPHA_NADAL,  BETA_NADAL,
     &     COEF_C_MAIGNAN,
     &     RHO,
     &     FICSURF,
     &     FICSOS_LOG,
     &     FICSOS_RES_BIN,
     &     FICTRANS, FICFLUX,
     &     ZOUT,
     &     IGMAX,
     &     IPOLAR,
     &     IDIRSURF,
     &     ITRPHI, PHIOS, PAS_PHI,
     &     IMODE_CKD_CALCUL,
     &     IER, TRACE,
     &     LUM_NBMU,
     &     IND_ANGOUT_FIN,
     &     PHI_FIN, THETA_FIN, 
     &     SCA_UP_FIN, I_UP_FIN, Q_UP_FIN, U_UP_FIN,
     &     POL_ANG_UP_FIN, POL_RATE_UP_FIN, L_POL_UP_FIN,
     &     SCA_DOWN_FIN, I_DOWN_FIN, Q_DOWN_FIN, U_DOWN_FIN,
     &     POL_ANG_DOWN_FIN, POL_RATE_DOWN_FIN, L_POL_DOWN_FIN,
     &     TDIR_VRAI, FLUX_DIFF_DOWN, FLUX_DOWN, EPLUS,
     &     COEF_TRONCA)
      
      IMPLICIT NONE
      
C Constants
C -----------------            
      
! Identifier number of the trace file for PROC.F
      INTEGER*2 IDLOG	    
      PARAMETER(IDLOG=INCTE_IDLOG_MAIN)

! Identifier number of the trace file for PROC.F
      INTEGER*4 IDLOGPROF
      PARAMETER(IDLOGPROF=INCTE_IDLOG_PROFIL)

! Identifier number of the trace file for PROC.F
      INTEGER*4 IDLOGOS
      PARAMETER(IDLOGOS=INCTE_IDLOG_OS)
    
C AEROSOLS constants
C -------------------

!Number of possible WMO or Shettle & Fenn models.
      INTEGER CTE_NB_Modeles
      PARAMETER (CTE_NB_Modeles=INCTE_NB_Modeles)

C PROFIL constants
C -----------------

!Max value of the profile parameter
      INTEGER CTE_NB_PROFIL
      PARAMETER (CTE_NB_PROFIL=INCTE_NB_PROFIL)


     	            
C* Definition of variables
C*-----------------------------------------------------------------                                    
      
C-- Variables related to the working directory tree

C     Length of the character string RESROOT
      INTEGER*4  LDIR

C     Length of the character string DIR_SURF
      INTEGER*4  LDIRSURF

C     Length of the character string DIR_SURF/GLITTER
      INTEGER*4  LDIRGLIT

C     Length of the character string DIR_SURF/ROUJEAN
      INTEGER*4  LDIRROUJ

C     Length of the character string DIR_SURF/RH
      INTEGER*4  LDIRRH

C     Length of the character string DIR_SURF/BREON
      INTEGER*4  LDIRBREON

C     Length of the character string DIR_SURF/NADAL
      INTEGER*4  LDIRNADAL

C     Length of the character string DIR_SURF/MAIGNAN
      INTEGER*4  LDIRMAIGNAN

C     Length of the character string DIR_MIE
      INTEGER*4  LDIRMIE

C     Length of the character string DIR_SOS
      INTEGER*4  LDIRSOS

C     Length of the character string DIR_LOG
      INTEGER*4  LDIRLOG

C     Length of the character string DIR_TMP
      INTEGER*4  LDIRTMP

C     Length of the character string of a file name
      INTEGER*4 LFIC

C     Output directory
      CHARACTER*CTE_LENDIR RESROOT

C     Output directory of the surface computations
      CHARACTER*CTE_LENDIR DIR_SURF

C     Output directory of the glitter computations
      CHARACTER*CTE_LENDIR DIR_GLIT

C     Output directory of the Roujean BRDF computations
      CHARACTER*CTE_LENDIR DIR_ROUJ

C     Output directory of the Rondeaux & Herman BRDF computations
      CHARACTER*CTE_LENDIR DIR_RH

C     Output directory of the Breon BRDF computations
      CHARACTER*CTE_LENDIR DIR_BREON

C     Output directory of the Nadal BRDF computations
      CHARACTER*CTE_LENDIR DIR_NADAL

C     Output directory of the Maignan BRDF computations
      CHARACTER*CTE_LENDIR DIR_MAIGNAN

C     Output directory of the Mie computations
      CHARACTER*CTE_LENDIR DIR_MIE

C     Output directory of the successive orders
      CHARACTER*CTE_LENDIR DIR_SOS

C     Output directory of the trace files
      CHARACTER*CTE_LENDIR DIR_LOG

C     Output directory of the temporary files
      CHARACTER*CTE_LENDIR DIR_TMP


C-- Global variables

C     Wavelength of the radiance simulation (in microns)
      DOUBLE PRECISION WA_SIMU

C     Length of the character string COMMAND
      INTEGER*4 LCOM   
    
C     Name of the main trace file
      CHARACTER*CTE_LENFIC2 FICMAIN_LOG

C     Command string for system call.
      CHARACTER*CTE_LENCOM  COMMAND

C     Array indices
      INTEGER*4 I,J,K

C     Error code =0 if no error, =-1 otherwise
      INTEGER*4 IER

C     True if writing in the trace file
      LOGICAL TRACE

C     True if writing in the trace file of the profile
      LOGICAL TRACEPROF

C     True if writing in the trace file of the OS
      LOGICAL TRACEOS

C     Checks the existence of a file
      LOGICAL EX

                  
C-- Variables related to the definition of the angles

C     Solar zenith angle in degrees (0 < tetas < 90).
C     (complementary angle: mus= -cos(tetas))
      DOUBLE PRECISION TETAS

C     Number of Gauss angles used for the radiance computations
      INTEGER*4 NBMU_GAUSS_LUM

C     Number of Gauss angles used for the Mie computations
      INTEGER*4 NBMU_GAUSS_MIE  


      !Radiance computations:
C     Path and name of the "user" angles file
      CHARACTER*CTE_LENFIC2 FICANGLES_USER_LUM

C     Path and name of the file of the angles used
      CHARACTER*CTE_LENFIC2 FICANGLES_RES_LUM

      !Phase function computations:
C     Path and name of the "user" angles file
      CHARACTER*CTE_LENFIC2 FICANGLES_USER_MIE

C     Path and name of the file of the angles used
      CHARACTER*CTE_LENFIC2 FICANGLES_RES_MIE
  
C     Name of the trace file
      CHARACTER*CTE_LENFIC2 FICANGLOG
 
      
C---- Variables related to the aerosol model

C     Aerosol optical thickness at the reference wavelength (not truncated).
      DOUBLE PRECISION AOT_REF

C     Reference wavelength
      DOUBLE PRECISION WAREF_AOT


      !Mono-modal aerosol case:
      !----------------------       
C     Real part of the refractive index
      DOUBLE PRECISION RN_WA

C     Imaginary part of the refractive index
      DOUBLE PRECISION IN_WA

C     Real part of the refractive index at the reference wavelength
      DOUBLE PRECISION RN_WAREF

C     Imaginary part of the refractive index at the reference wavelength
      DOUBLE PRECISION IN_WAREF

C     Modal radius of the LND
      DOUBLE PRECISION LND_RADIUS_MMD_AER
	
C     Variance of the LND
      DOUBLE PRECISION LND_LNVAR_MMD_AER

C     Slope of the Junge law.
      DOUBLE PRECISION JD_SLOPE_MMD_AER

C     Minimum radius of the Junge law.
      DOUBLE PRECISION JD_RMIN_MMD_AER

C     Maximum radius of the Junge law.
      DOUBLE PRECISION JD_RMAX_MMD_AER
      
C     General names of the size distribution parameters for a mono-modal model
      DOUBLE PRECISION VARGRANU1_MMD_AER
      DOUBLE PRECISION VARGRANU2_MMD_AER
      DOUBLE PRECISION VARGRANU3_MMD_AER


      !Volume proportion of the WMO components
      !-----------------------------------------     
C     User value for the DL particles of the WMO model
      DOUBLE PRECISION C_WMO_DL

C     User value for the WS particles of the WMO model
      DOUBLE PRECISION C_WMO_WS

C     User value for the OC particles of the WMO model
      DOUBLE PRECISION C_WMO_OC

C     User value for the SO particles of the WMO model
      DOUBLE PRECISION C_WMO_SO

C     Percentage of relative humidity of the Shettle & Fenn models
      DOUBLE PRECISION RH

    
      !Parameters for a bimodal LND model
      !-----------------------------------------    
C     User value of the volume concentration of the "LND coarse"
      DOUBLE PRECISION USER_CV_COARSE
   
C     User value of the volume concentration of the "LND fine"
      DOUBLE PRECISION USER_CV_FINE

C     Ratio of the optical thicknesses AOT_coarse / AOT_total at the reference wavelength
      DOUBLE PRECISION RTAUct_WAREF
 
C     Value read for RN_LND(1) : coarse mode
      DOUBLE PRECISION BMD_CM_MRWA

C     Value read for IN_LND(1) : coarse mode
      DOUBLE PRECISION BMD_CM_MIWA

C     Value read for RN_LND_WAREF(1) : coarse mode
      DOUBLE PRECISION BMD_CM_MRWAREF

C     Value read for RN_LND_WAREF(2) : fine mode
      DOUBLE PRECISION BMD_CM_MIWAREF

C     Value read for RMODAL_LND(1) : coarse mode
      DOUBLE PRECISION BMD_CM_RMODAL

C     Value read for VAR_LND(1) : coarse mode
      DOUBLE PRECISION BMD_CM_VAR

C     Value read for RN_LND(2) : fine mode
      DOUBLE PRECISION BMD_FM_MRWA

C     Value read for IN_LND(2) : fine mode
      DOUBLE PRECISION BMD_FM_MIWA

C     Value read for IN_LND_WAREF(1) : coarse mode
      DOUBLE PRECISION BMD_FM_MRWAREF

C     Value read for IN_LND_WAREF(2) : fine mode
      DOUBLE PRECISION BMD_FM_MIWAREF 
   
C     Value read for RMODAL_LND(2) : fine mode
      DOUBLE PRECISION BMD_FM_RMODAL

C     Value read for VAR_LND(2) : fine mode
      DOUBLE PRECISION BMD_FM_VAR

C     Truncation option of the aerosol phase function (1 to apply it)
      INTEGER*4 ITRONC_AER

C     Truncation option (=ITRONC_AER as input of the call of SOS_AEROSOLS
C     but can switch from 1 to 0 as output if the truncation is automatically cancelled)
      INTEGER*4 ITRONC_AER_IN

C     Mono- or multi-modal size distribution
      INTEGER*2 IMOD_AER

C     Number of the size distribution model
      INTEGER*4 IGRANU

C     Type of the WMO model
      INTEGER*2 IModele_WMO

C     Type of the Shettle & Fenn model
      INTEGER*2 IModele_SF

C     Index of the choice of description of the mixture of 2 LND
      INTEGER*2 MODE_PARAM_BILND

C     Name of the TRACE file of the MIE computations
      CHARACTER*CTE_LENFIC2 FICMIE_LOG

C     Name of the TRACE file of the GRANU aerosol property computations
      CHARACTER*CTE_LENFIC2 FICGRANU_LOG

C     Name of the GRANU result file
      CHARACTER*CTE_LENFIC2 FICGRANU

C     Path and name of the user file giving the phase function
C     and the scattering and extinction coefficients
      CHARACTER*CTE_LENFIC2 FICEXTDATA_AER

C     Path and name of the aerosol file named by the user
      CHARACTER*CTE_LENFIC2 FICUSER_AER

C     Path and name of the file defining a mixture of aerosol models by the user
      CHARACTER*CTE_LENFIC2 FICMIXTURE_AER  
	
 
			
C --- Variables related to the definition of the atmospheric profile

C     Rayleigh optical thickness
      DOUBLE PRECISION TR

C     Scale height of the molecular profile (km)
      DOUBLE PRECISION HR

C     Scale height of the aerosol profile (km)
      DOUBLE PRECISION HA

C     Minimum and maximum altitudes of the molecules + aerosols layer (km)
      DOUBLE PRECISION ZMIN
      DOUBLE PRECISION ZMAX
  
C     Atmospheric pressure at the surface (mbar)
      DOUBLE PRECISION PSURF

C     Integrated water vapour content (g/cm2)
      DOUBLE PRECISION H2O  

C     Integrated ozone content (Dobson)
      DOUBLE PRECISION O3

C     Carbon dioxide concentration at the surface (ppmv)
      DOUBLE PRECISION CO2

C     Methane concentration at the surface (ppmv)
      DOUBLE PRECISION CH4

C     Direct transmission for the solar incidence direction of the equivalent atmosphere
      DOUBLE PRECISION TDIR_TRONC

C     Direct transmission for the solar incidence direction of the true atmosphere
      DOUBLE PRECISION TDIR_VRAI

C     Altitude (km)
      DOUBLE PRECISION ZALT

C     Molecular, aerosol and gaseous absorption optical thickness at ZALT
      DOUBLE PRECISION TR_ZALT
      DOUBLE PRECISION TA_ZALT
      DOUBLE PRECISION TG_ZALT

C     Type of the aerosol profile
      INTEGER*2 IPROFIL

C     Type of the absorption profile
      INTEGER*2 ABSPROFIL

C     Number of atmospheric layers
      INTEGER*4 NT

C     Absorption profile file provided by the user
      CHARACTER*CTE_LENFIC2 FICABSPROFIL

C     Name of the temporary PROFIL result file for one order of the CKD decomposition
      CHARACTER*CTE_LENFIC2 FICPROFIL_TMP

C     Path and name of the trace file of the profiles per CKD order
      CHARACTER*CTE_LENFIC2 FICPROFIL_LOG

 
C --- Variables related to the choice of the CKD coefficient data
C     Wavenumber step (cm-1) in the CKD coeff file
      DOUBLE PRECISION  NUSTEP

            
C --- Variables related to the definition of the surface model

C     Surface type
      INTEGER*2 ISURF 

C     Surface albedo
      DOUBLE PRECISION RHO

C     Wind speed (m/s)
      DOUBLE PRECISION WIND 

C     Refractive index of the air/surface interface
      DOUBLE PRECISION SURF_IND

C     Parameters of the Roujean model
      DOUBLE PRECISION K0_ROUJEAN, K1_ROUJEAN, K2_ROUJEAN

C     Parameters of the Nadal model
      DOUBLE PRECISION ALPHA_NADAL, BETA_NADAL

C     Parameter C*exp(-NDVI) of the Maignan model
      DOUBLE PRECISION COEF_C_MAIGNAN

C     Name of the SURFACE trace file
      CHARACTER*CTE_LENFIC2 FICSURF_LOG

C     Name of a file containing an existing reflection matrix
      CHARACTER*CTE_LENFIC2 FICSURF

  
C --- Variables related to the SOS computations

C     Output type: -1 =>  Standard output: Radiance emerging at TOA + downward radiance at ground level.
C                       z =>  Upward and downward radiance for the altitude z (0 <= z <= CTE_TOA_ALT).
      DOUBLE PRECISION ZOUT
 
C     Maximum order of scattering
      INTEGER*4 IGMAX

C     Option for simulation without polarisation: 0 : computations without polarisation
C                                              1 : normal computations with polarisation
      INTEGER*4 IPOLAR

C     Option for the direct reflection of the solar beam by the surface:
C          0 : diffuse field only, 1 : complete field
      INTEGER*4 IDIRSURF

C     Surface options passed to the azimuth recombination (set to 0 for
C     IDIRSURF = 0: no direct reflection term added)
      INTEGER*4 IGLI_D, IFRESNEL_D, IROUJEAN_D, IRONDEAUX_D, IBREON_D
      INTEGER*4 INADAL_D, IMAIGNAN_D

C     Option of the mode of application of the CKD decomposition
C          1 : Finest computation with one radiance computation per order of the decomposition
C              then recomposition with the CKD weighting.
C          2 : Computation for an absorption optical thickness estimated by the CKD approach
C              --> A single radiance computation
      INTEGER*4 IMODE_CKD_CALCUL

C     Index of output type: 1 =>  Viewing plane for one azimuth value.
C                                2 =>  Polar diagram versus the azimuth and viewing directions.
      INTEGER*4 ITRPHI        

C     Azimuth step in degrees (for ITRPHI=2)
      INTEGER*4 PAS_PHI

C     Relative azimuth in degrees (for ITRPHI=1)
      DOUBLE PRECISION PHIOS

C     Molecular depolarisation factor
      DOUBLE PRECISION RON

C     Path and name of the trace file of the SOS_OS computations
      CHARACTER*CTE_LENFIC2 FICSOS_LOG

C     Name of the optional output file of the transmissions
      CHARACTER*CTE_LENFIC2 FICTRANS

C     Name of the optional output file of the fluxes up TOA and down BOA
      CHARACTER*CTE_LENFIC2 FICFLUX

C     Name of the OS result file (for a unit processing)
      CHARACTER*CTE_LENFIC2 FICSOS_RES_BIN

C     Name of the temporary OS result file (aggregated computation)
      CHARACTER*CTE_LENFIC2 FICOS_AGG_TMP



C --- Internal variables

C -   For the aerosols

C     Wavelength of computation of the aerosol parameters (WA_SIMU or WAREF)
      DOUBLE PRECISION WA

C     Aerosol extinction cross section for the computed wavelength (simu or ref)
      DOUBLE PRECISION KMAT1

C     Aerosol extinction cross section for the reference wavelength of the aerosol opt. thickness
      DOUBLE PRECISION KMAT1_WAREF
	
C     Aerosol optical thickness at the simulation wavelength
      DOUBLE PRECISION TA

C     Aerosol single scattering albedo at the simulation wavelength
      DOUBLE PRECISION PIZAER

C     Truncation coefficient of the aerosol phase function
      DOUBLE PRECISION COEF_TRONCA


C -   For the surfaces

C     Result directory of the surface BRDF/BPDF computations
      CHARACTER*CTE_LENDIR  SURF_OUTDIR
      
C -   For the profiles

C     Wavenumber (cm-1)
      DOUBLE PRECISION NU

C     Absorption profile
      DOUBLE PRECISION ABS_USERPROFIL(CTE_ABS_NBLEV,CTE_ABS_NBCOL)

C     Altitudes of the levels of the absorption profile (decreasing order)
      DOUBLE PRECISION ALTABS(1:CTE_ABS_NBLEV)

C     Profile of gaseous absorption optical thickness
      DOUBLE PRECISION TAUABS(1:CTE_ABS_NBLEV)

C     Profile of vertical transmission (fictitious air mass of 1)
      DOUBLE PRECISION TRSCKD(1:CTE_ABS_NBLEV)

C     Number of particles/cm2 per gas species (H2O CO2 O3 N2O CO CH4 O2 NO2) and per layer of the atmospheric profile
      DOUBLE PRECISION RO(CTE_NBABS,CTE_ABS_NBLEV)
 
C     Array identifying the absorbing gases taken into account (H2O CO2 O3 N2O CO CH4 O2 NO2)
C     Includes the gases with continuum absorption only (NO2).
      INTEGER*2 IABS(1:CTE_NBABS)

C     Index of the CKD arrays (depends on NU)
      INTEGER*4 LAMB1

C	- CKD coefficients

C     Number of exponentials per wavelength and species
      INTEGER*4  NEXP(CTE_NBABS,CTE_CKD_NWVL_MAX)

C     Values of the coefficients ai per wavelength, species and order i of the decomposition in exponentials
      DOUBLE PRECISION  KDIS_AI(CTE_CKD_NAI_MAX, CTE_NBABS,  
     &                          CTE_CKD_NWVL_MAX) 
     
C     Values of the coefficients ki per wavelength, species, order i of the decomposition in exponentials,
C     pressure and temperature
      DOUBLE PRECISION  KDIS_KI(CTE_CKD_NT_MAX, CTE_CKD_NP_MAX, 
     &                          CTE_CKD_NAI_MAX, CTE_NBABS,  
     &                          CTE_CKD_NWVL_MAX )
 
C     Values of the coefficients ki of water vapour per wavelength, order i of the decomposition in exponentials,
C     concentration, pressure and temperature
      DOUBLE PRECISION  KDIS_KI_H2O(CTE_CKD_NT_MAX, CTE_CKD_NP_MAX, 
     &                              CTE_CKD_NC_MAX, CTE_CKD_NAI_MAX, 
     &                              CTE_CKD_NWVL_MAX )

   
C     Number of temperature and pressure values in the CKD files
      INTEGER NB_TEMP, NB_PRES

C     Array of the temperature values read
      DOUBLE PRECISION TAB_TEMP(CTE_CKD_NT_MAX)

C     Array of the pressure values read
      DOUBLE PRECISION TAB_PRES(CTE_CKD_NP_MAX)

C     Number of H2O concentration values in the H2O CKD file
      INTEGER NB_CONC_H2O 

C     Array of the H2O concentration values read in the H2O CKD file (ppmv)
      DOUBLE PRECISION TAB_CONC_H2O(CTE_CKD_NC_MAX) 

C     Exponential index
      INTEGER*4 IK1,IK2,IK3,IK4,IK5,IK6,IK7,IK8

C     Weighting coefficient of the radiances
      DOUBLE PRECISION  AIK

C     Sum of the weighting coefficients of the radiances
      DOUBLE PRECISION sumAIK
   
C -   For the OS

C     Order of the expansions of the phase functions in Legendre polynomials.
      INTEGER*4 OS_NB

C     Single scattering albedo
      DOUBLE PRECISION PIZ

C     Single scattering albedo after truncation
      DOUBLE PRECISION PIZTR

C     Truncation coefficient of the aerosol phase function
      DOUBLE PRECISION COEFTRONCA

C     Number of actual angles of the simulations
      INTEGER*4 LUM_NBMU

C     Cosines of the Gauss (and user) angles
      DOUBLE PRECISION RMU(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  

C     Weights of the Gauss (and user) angles
      DOUBLE PRECISION GA(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
				
C     Parameters alpha(k), beta(k), gamma(k) and zeta(k) of the aerosol phase matrix
      DOUBLE PRECISION AER_ALPHA(0:CTE_OS_NB_MAX)
      DOUBLE PRECISION AER_BETA(0:CTE_OS_NB_MAX)
      DOUBLE PRECISION AER_GAMMA(0:CTE_OS_NB_MAX)
      DOUBLE PRECISION AER_ZETA(0:CTE_OS_NB_MAX)

C     Diffuse flux illuminating the ground and emerging at TOA, normalised by the incident solar flux at TOA
      DOUBLE PRECISION EMOINS
      DOUBLE PRECISION EPLUS

C     Diffuse transmission for the solar incidence direction of the equivalent atmosphere
      DOUBLE PRECISION TDIFMUS

C     Diffuse transmission for the Gauss angles of the equivalent atmosphere
      DOUBLE PRECISION TDIFMUG(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  

C     Total optical thickness at ground (truncated)
      DOUBLE PRECISION TTOT_TRONC

C     Total optical thickness at ground (not truncated).
      DOUBLE PRECISION TTOT_VRAI

C     Optical thickness (truncated) at the output level of the results for the upward field
      DOUBLE PRECISION TAUOUT

C     Previous variables of fluxes and optical thicknesses computed
C     for an octet ik1*ik2*ik3*ik4*ik5*ik6*ik7*ik8 of the CKD decomposition)
      DOUBLE PRECISION EMOINS_TMP
      DOUBLE PRECISION EPLUS_TMP
      DOUBLE PRECISION TDIFMUS_TMP		
      DOUBLE PRECISION TDIFMUG_TMP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
      DOUBLE PRECISION TTOT_TRONC_TMP
      DOUBLE PRECISION TTOT_VRAI_TMP
      DOUBLE PRECISION TAUOUT_TMP

C     Direct transmission of the solar beam TOA -> surface -> output level
C     along the direction J (used for the direct surface reflection)
      DOUBLE PRECISION AT_DIR(CTE_OS_NBMU_MAX)
C     .TRUE. if AT_DIR has been accumulated over the CKD terms
      LOGICAL ATDIR_CKD
C     Direct transmission TOA -> surface for the truncated (TRONC) and the
C     true (VRAI) atmosphere: index 0 for the solar zenith angle TETAS,
C     index J for the direction RMU(J)
      DOUBLE PRECISION TDIR_TRONC_MU(0:CTE_OS_NBMU_MAX)
      DOUBLE PRECISION TDIR_VRAI_MU(0:CTE_OS_NBMU_MAX)
      DOUBLE PRECISION XMUS

C     Number of the Gauss mu for the solar zenith angle
C     complementary angle: mus = -rmu(N0)
C     N0 > 0 => to take it into account
      INTEGER*4 N0
   
C     Index specifying whether the simulations include the Fresnel reflection on a flat sea:
C         1 => for a computation of the Fresnel reflection.
C         Otherwise => not taken into account.
      INTEGER*4 IFRESNEL

C     Index specifying whether the simulations use a BRDF or BPDF reflection model:
C         1 => for BRDF or BPDF computation (N0 > 0 must be used).
C         Otherwise => no BRDF or BPDF simulation.
      INTEGER*4 IMAT_SURF

C     Index specifying whether the simulations include the Glitter:
C         1 => for Glitter computation (N0 > 0 must then be used).
C         Otherwise => no Glitter simulation.
      INTEGER*4 IGLI
          
C     Index specifying whether the simulations use the Roujean reflection model:
C         1 => for Roujean BRDF computation (N0 > 0 must be used).
C         Otherwise => no Roujean BRDF simulation.
      INTEGER*4 IROUJEAN

C     Index specifying whether the simulations use the Rondeaux - Herman reflection model:
C         1 => for RH BPDF computation (N0 > 0 must be used).
C         Otherwise => no RH BPDF simulation.
      INTEGER*4 IRONDEAUX
			 
C     Index specifying whether the simulations use the Breon reflection model:
C         1 => for Breon BPDF computation (N0 > 0 must be used).
C         Otherwise => no Breon BPDF simulation.
      INTEGER*4 IBREON

C     Index specifying whether the simulations use the Nadal reflection model:
C         1 => for Nadal BPDF computation (N0 > 0 must be used).
C         Otherwise => no Nadal BPDF simulation.
      INTEGER*4 INADAL

C     Index specifying whether the simulations use the Maignan reflection model:
C         1 => for Maignan BPDF computation (N0 > 0 must be used).
C         Otherwise => no Maignan BPDF simulation.
      INTEGER*4 IMAIGNAN

C     Indices of the angle arrays for which there must be an output in the result file of the user angles
C     ==> 1: output required
C                                      Indices starting from 0 for Python.
      INTEGER*4 IND_ANGOUT_FIN(0:CTE_OS_NBMU_MAX)  
				
C     Length of the character string FICSOS_RES_BIN
      INTEGER*4 LFICOS

C     Name of the temporary OS result file
      CHARACTER*CTE_LENFIC2 FICOS_TMP
	
C     Table of the relative azimuth values of the output field
      DOUBLE PRECISION PHI_FIN(0:360)

C     Table of the viewing zenith angle values of the output field
      DOUBLE PRECISION THETA_FIN(0:CTE_OS_NBMU_MAX)

C     Tables of the scattering angles for the output upward and downward fields
      DOUBLE PRECISION SCA_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)
      DOUBLE PRECISION SCA_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX)     

C     Tables of the I, Q and U values for the output upward and downward fields
      DOUBLE PRECISION I_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)        
      DOUBLE PRECISION Q_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)       
      DOUBLE PRECISION U_UP_FIN(0:360,0:CTE_OS_NBMU_MAX) 
      DOUBLE PRECISION I_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX)     
      DOUBLE PRECISION Q_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX)      
      DOUBLE PRECISION U_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX) 

C     Tables of the polarisation angles for the output upward and downward fields
      DOUBLE PRECISION POL_ANG_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)       
      DOUBLE PRECISION POL_ANG_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX)  

C     Tables of the degrees of polarisation for the output upward and downward fields
      DOUBLE PRECISION POL_RATE_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)  
      DOUBLE PRECISION POL_RATE_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX)  

C     Table of the polarised radiance values for the output upward and downward fields
      DOUBLE PRECISION L_POL_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)       
      DOUBLE PRECISION L_POL_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX)

C     Downward diffuse flux at the surface level
      DOUBLE PRECISION FLUX_DIFF_DOWN

C     Total downward flux at the surface level
      DOUBLE PRECISION FLUX_DOWN
  

    
C =====================================================================
C --- Definition of the I/O variables for the Python binding
C =====================================================================

Cf2py intent(in) resroot, ficmain_log,wa_simu, nbmu_gauss_lum, ficangles_user_lum, tetas
Cf2py intent(in)  ficangles_res_lum, nbmu_gauss_mie,ficangles_user_mie, ficangles_res_mie 
Cf2py intent(in) ficanglog, waref_aot, aot_ref, itronc_aer, ficgranu_log, ficmie_log, dir_mie
Cf2py intent(in)  ficgranu,imod_aer,rn_wa, in_wa, rn_waref, in_waref, igranu, lnd_radius_mmd_aer
Cf2py intent(in)  lnd_lnvar_mmd_aer, jd_slope_mmd_aer, jd_rmin_mmd_aer, jd_rmax_mmd_aer
Cf2py intent(in) imodel_wmo, c_wmo_dl, c_wmo_ws, c_wmo_oc, c_wmo_so 
Cf2py intent(in) imodele_sf, rh, mode_param_bilnd, user_cv_coarse, user_cv_fine 
Cf2py intent(in) rtauct_waref,bmd_cm_mrwa, bmd_cm_miwa, bmd_cm_mrwaref, bmd_cm_miwaref 
Cf2py intent(in) bmd_cm_rmodal, bmd_cm_var, bmd_fm_mrwa, bmd_fm_miwa, bmd_fm_mrwaref
Cf2py intent(in) bmd_fm_miwaref, bmd_fm_rmodal, bmd_fm_var, ficextdata_aer, ficmixture_aer 
Cf2py intent(in) ficuser_aer, ficprofil_log, tr, hr, ha, iprofil, zmin, zmax,psurf, h2o
Cf2py intent(in)  o3, co2, ch4, absprofil, ficabsprofil, nustep, isurf, dir_surf, ficsurf_log
Cf2py intent(in)  surf_ind, wind, k0_roujean, k1_roujean, k2_roujean, alpha_nada, beta_nada
Cf2py intent(in) coef_c_maignan, rho, ficsurf, ficsos_log, ficsos_res_bin, fic_trans 
Cf2py intent(in) fic_flux, zout, igmax, ipolar, itrphi, phios, pas_phi, imode_ckd_calcul, ier, trace

Cf2py intent(out) lum_nbmu, ind_angout_fin, phi_fin, theta_fin, sca_up_fin, i_up_fin, q_up_fin, u_up_fin, pol_ang_up_fin, pol_rate_up_fin, l_pol_up_fin
Cf2py intent(out) sca_down_fin, i_down_fin, q_down_fin, u_down_fin, pol_ang_down_fin, pol_rate_down_fin, l_pol_down_fin
Cf2py intent(out) tdir_vrai, flux_diff_down, flux_down, eplus
Cf2py intent(out) coef_tronca

                        
      
C =========================================
C --- Initialisations
C =========================================

C --- Initialisation of the processing error flag
C----------------------------------------------------------
      IER=0
 
C     Output directory
      DIR_ROUJ="NO_ROUJEAN_DIR"

C     Output directory
      SURF_OUTDIR="NO_OUTPUT_DIR"


C     Absorption profile versus altitude
      DO I=1,CTE_ABS_NBLEV
         ALTABS(I)=CTE_NOT_DEFINED_VALUE_DBLE
         TAUABS(I)=CTE_NOT_DEFINED_VALUE_DBLE
      ENDDO
 
C     Number of coefficients for each absorbing gas and wavelength
C     Set to 1 by default => case without absorption or continuum only
      DO J=1,CTE_CKD_NWVL_MAX
         DO I=1,CTE_NBABS
            NEXP(I,J)=1
         ENDDO
      ENDDO

C     Weight of each exponential for each absorbing gas and wavelength
C     Set to 1 by default => case without absorption or continuum only
      DO K=1,CTE_CKD_NWVL_MAX	
         DO J=1,CTE_NBABS
            DO I=1,CTE_CKD_NAI_MAX
               KDIS_AI(I,J,K)=1.		
            ENDDO
         ENDDO
      ENDDO

C     Default sum of the coefficients AIK
      sumAIK=1.D+00

C     Index corresponding to the wavenumber / 10
      LAMB1=1
            
C --- Initialisation of the default values of the SOS parameters

C     Molecular depolarisation factor.
      RON=CTE_MDF

C     Init of the optical thicknesses for the aggregation
      TTOT_VRAI=0.
      TTOT_TRONC=0.
      TAUOUT=0.
      
C     Init of the diffuse transmissions
      TDIFMUS=0.
      DO I=-CTE_OS_NBMU_MAX,CTE_OS_NBMU_MAX
         TDIFMUG(I)=0.
      ENDDO	  
      EMOINS=0.
      EPLUS=0.


C =========================================================
C --- Definition of the working directory tree
C =========================================================           				

C ---- Check of the definition of the working space (mandatory)
      IF (RESROOT.EQ."UNDEFINED_REPERTORY") GOTO 1000
   
      LDIR = INDEX(RESROOT,' ') - 1
      IF((LDIR.LE.0).OR.(LDIR.GT.CTE_LENDIR)) GOTO 1100

C ---- Check of the definition of the MIE result directory
C     (mandatory if AOT is not null and if no external data file is used:
C      phase functions or directly the Greek coefficients)

C     Definition of the aerosol opt. thickness at the reference wavelength: mandatory
      IF (AOT_REF.EQ.CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2301  	

      IF (AOT_REF.GT.0.D+00) THEN
		
C        Definition of the type of aerosol model: mandatory
         IF ((IMOD_AER.EQ.CTE_NOT_DEFINED_VALUE_INT)
     &        .AND.(FICUSER_AER.EQ."NO_USER_AEROSOLS" )) GOTO 2304

C        Definition of the reference wavelength for the aerosol opt. thickness: mandatory
         IF (WAREF_AOT.EQ.CTE_NOT_DEFINED_VALUE_INT) GOTO 2302  

         IF ((IMOD_AER.NE.4.)
     &        .AND.(FICUSER_AER.EQ.'NO_USER_AEROSOLS')) THEN

            IF (DIR_MIE.EQ."UNDEFINED_REPERTORY") GOTO 1001

            !Creation of the aerosol result directory tree
            LDIRMIE = INDEX(DIR_MIE,' ') - 1
            IF((LDIRMIE.LE.0).OR.(LDIRMIE.GT.CTE_LENDIR)) GOTO 1101
            COMMAND='mkdir -p '//DIR_MIE(1:LDIRMIE)
            LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') + 9
            IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1202
            CALL SYSTEM(COMMAND)
         ENDIF    
         
      ENDIF


C ---- Definition of the surface result directory tree
      IF((ISURF.NE.0).AND.(ISURF.NE.2)) THEN

         !Check of the definition of the result directory DIR_SURF (required if a BRDF or BPDF is computed)
         IF (DIR_SURF.EQ."UNDEFINED_REPERTORY")  GOTO 1002

         LDIRSURF = INDEX(DIR_SURF,' ') - 1
         IF((LDIRSURF.LE.0).OR.(LDIRSURF.GT.CTE_LENDIR)) GOTO 1102
         
         !Definition of the directories
         !--GLITTER directory
         LDIRGLIT= LDIRSURF + 8  	
         IF(LDIRGLIT.GT.CTE_LENDIR) GOTO 1103
         DIR_GLIT(1:LDIRGLIT) = DIR_SURF(1:LDIRSURF)//'/'//'GLITTER'
      
         !--ROUJEAN directory
         LDIRROUJ= LDIRSURF + 8  	
         IF(LDIRROUJ.GT.CTE_LENDIR) GOTO 1104
         DIR_ROUJ(1:LDIRROUJ) = DIR_SURF(1:LDIRSURF)//'/'//'ROUJEAN'      	
      
         !--Rondeaux & Herman directory
         LDIRRH= LDIRSURF + 3  	
         IF(LDIRRH.GT.CTE_LENDIR) GOTO 1105
         DIR_RH(1:LDIRRH) = DIR_SURF(1:LDIRSURF)//'/'//'RH'     

         !--BREON directory
         LDIRBREON= LDIRSURF + 6  	
         IF(LDIRBREON.GT.CTE_LENDIR) GOTO 1106
         DIR_BREON(1:LDIRBREON) = DIR_SURF(1:LDIRSURF)//'/'//'BREON'     
  
         !--NADAL directory
         LDIRNADAL= LDIRSURF + 6  	
         IF(LDIRNADAL.GT.CTE_LENDIR) GOTO 1107
         DIR_NADAL(1:LDIRNADAL) = DIR_SURF(1:LDIRSURF)//'/'//'NADAL'     

         !--MAIGNAN directory
         LDIRMAIGNAN= LDIRSURF + 8  	
         IF(LDIRMAIGNAN.GT.CTE_LENDIR) GOTO 11071
         DIR_MAIGNAN(1:LDIRMAIGNAN) 
     &   = DIR_SURF(1:LDIRSURF)//'/'//'MAIGNAN'

         !Creation of the surface result directory tree
         !--SURFACE directory
         COMMAND='mkdir -p '//DIR_SURF(1:LDIRSURF)
         LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') + 9
         IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1203
         CALL SYSTEM(COMMAND)

         !--GLITTER directory
         COMMAND='mkdir -p '//DIR_GLIT(1:LDIRGLIT)
         LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') + 9
         IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1204
         CALL SYSTEM(COMMAND)
     
         !--ROUJEAN directory
         COMMAND='mkdir -p '//DIR_ROUJ(1:LDIRROUJ)
         LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') + 9
         IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1205
         CALL SYSTEM(COMMAND)
     
         !--Rondeaux & Herman directory
         COMMAND='mkdir -p '//DIR_RH(1:LDIRRH)
         LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') + 9
         IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1206
         CALL SYSTEM(COMMAND)

         !--BREON directory
         COMMAND='mkdir -p '//DIR_BREON(1:LDIRBREON)
         LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') + 9
         IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1207
         CALL SYSTEM(COMMAND)
  
         !--NADAL directory
         COMMAND='mkdir -p '//DIR_NADAL(1:LDIRNADAL)
         LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') + 9
         IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1208
         CALL SYSTEM(COMMAND)

         !--MAIGNAN directory
         COMMAND='mkdir -p '//DIR_MAIGNAN(1:LDIRMAIGNAN)
         LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') + 9
         IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 12081
         CALL SYSTEM(COMMAND)

      ENDIF
		
C ---- Definition of the SOS result directory tree

      LDIRSOS= LDIR + 4  	
      IF(LDIRSOS.GT.CTE_LENDIR) GOTO 1108
      DIR_SOS(1:LDIRSOS) = RESROOT(1:LDIR)//'/'//'SOS'	

C ---- Definition of the log directory

      LDIRLOG= LDIR + 4  	
      IF(LDIRLOG.GT.CTE_LENDIR) GOTO 1109
      DIR_LOG(1:LDIRLOG) = RESROOT(1:LDIR)//'/'//'LOG'	


C ---- Definition of the temporary processing directory

      LDIRTMP= LDIR + 4  	
      IF(LDIRTMP.GT.CTE_LENDIR) GOTO 1110
      DIR_TMP(1:LDIRTMP) = RESROOT(1:LDIR)//'/'//'TMP'



C =========================================================
C --- Definition of the temporary files per CKD order
C =========================================================           				

C     Temporary profile file
      FICPROFIL_TMP = DIR_TMP(1:LDIRTMP)//'/'//"PROFIL_TMP"  

C     Temporary file of OS results (1 simulation)
      FICOS_TMP = DIR_TMP(1:LDIRTMP)//'/'//"FICOS_TMP"  

C     Temporary file of OS results (aggregated simulations)
      FICOS_AGG_TMP = DIR_TMP(1:LDIRTMP)//'/'//"FICOS_AGG_TMP.bin"


C =========================================================
C --- Creation of the working directory tree
C =========================================================           				

C ---- Working space

      COMMAND='mkdir -p '//RESROOT(1:LDIR)
      LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') + 9
      IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1201
      CALL SYSTEM(COMMAND)
	
C ---- Creation of the SOS result directory tree

      COMMAND='mkdir -p '//DIR_SOS(1:LDIRSOS)
      LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') + 9
      IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1209
      CALL SYSTEM(COMMAND)

C ---- Creation of the log directory

      COMMAND='mkdir -p '//DIR_LOG(1:LDIRLOG)
      LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') + 9
      IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1210
      CALL SYSTEM(COMMAND)


C ---- Creation of the directory of the temporary files

      COMMAND='mkdir -p '//DIR_TMP(1:LDIRTMP)
      LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') + 9
      IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1211
      CALL SYSTEM(COMMAND)

      
C =========================================================
C --- Opening of the trace file
C =========================================================        

C     ---- Definition of the main trace file
      IF (TRACE) THEN

         LFIC = INDEX(FICMAIN_LOG,' ') - 1
         IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 1301
         FICMAIN_LOG = DIR_LOG(1:LDIRLOG)//'/'//FICMAIN_LOG(1:LFIC)
       
C ---- Creation of the file
         OPEN (IDLOG,FILE=FICMAIN_LOG,ERR=1302)

         WRITE (IDLOG,*,ERR=1303) "WORKING FOLDER DEFINITION"
         WRITE (IDLOG,*,ERR=1303) "-------------------------"
         WRITE (IDLOG,*,ERR=1303) "Working folder :", RESROOT(1:LDIR)
         WRITE (IDLOG,*,ERR=1303) " "
         WRITE (IDLOG,*,ERR=1303) 
     &        "Surface result folder :", DIR_SURF(1:LDIRSURF)
         WRITE (IDLOG,*,ERR=1303) " "
         WRITE (IDLOG,*,ERR=1303) 
     &        "Mie result folder :", DIR_MIE(1:LDIRMIE)
         WRITE (IDLOG,*,ERR=1303) " "
         WRITE (IDLOG,*,ERR=1303) 
     &        "SOS result folder :", DIR_SOS(1:LDIRSOS)     
         WRITE (IDLOG,*,ERR=1303) " "
      ENDIF  ! End of Trace
  
	    
C =========================================================
C --- Existence and consistency tests of the parameters
C =========================================================     

C --------------------------------------------------------------        
C --- Check of the simulation parameters
C --------------------------------------------------------------  
C 
      IF (TRACE) THEN
         WRITE (IDLOG,*,ERR=1303) " "
         WRITE (IDLOG,*,ERR=1303) 
     &         "--> Starting the check of user parameters \n"
      ENDIF

C --- Check of the definition of the simulation wavelength
	     
      IF (WA_SIMU.EQ.CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2100           
      IF ((WA_SIMU.LT.CTE_WAMIN).OR.(WA_SIMU.GT.CTE_WAMAX))  GOTO 2101
      
      IF (TRACE) THEN

	  WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) "************** ",
     &         "Simulation parameters ***************"
         WRITE (IDLOG,*,ERR=1303) " "
          
C        Simulation wavelength
	  WRITE (IDLOG,*,ERR=1303) 
     &         "Simulation wavelength (microns): ",WA_SIMU
	  WRITE (IDLOG,*,ERR=1303) " "
	  
      ENDIF 
      
C --------------------------------------------------------------        
C --- Check of the arguments defining the angles
C --------------------------------------------------------------  
C  
C ---- Existence and validity test of the mandatory parameters

C     Solar angle: mandatory
      IF (TETAS.EQ.CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2200
      IF ((TETAS.LT.0.D-00).OR.(TETAS.GE.90.D-00)) GOTO 2201
         
	              
C ---- Display of the parameters used in the trace file
	     
      IF (TRACE) THEN

	  WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) "************** ",
     &         "Angles parameters ***************"
         WRITE (IDLOG,*,ERR=1303) " "	  
	  WRITE (IDLOG,*,ERR=1303) "***** Mandatory Parameters "
	  
C        Solar angle
	  WRITE (IDLOG,*,ERR=1303) "Solar Zenith Angle: ",TETAS
          
C        Result file for the radiances
	  WRITE (IDLOG,*,ERR=1303) "User's filename ",
     &         "for the radiance angles définition :"
         WRITE (IDLOG,*,ERR=1303) FICANGLES_RES_LUM
               
C        Result file for the phase functions
	  WRITE (IDLOG,*,ERR=1303) "User's filename ",
     &         "for the Mie angles définition :"
         WRITE (IDLOG,*,ERR=1303) FICANGLES_RES_MIE
	  	   
C        Log file for the angle computation
	  WRITE (IDLOG,*,ERR=1303) "User's log file: "
         WRITE (IDLOG,*,ERR=1303) FICANGLES_RES_MIE
               
	  WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) "***** Optional Parameters "

C        Number of Gauss angles for the radiance simulations
	  IF (NBMU_GAUSS_LUM.NE.CTE_NOT_DEFINED_VALUE_INT) THEN 
	     WRITE (IDLOG,*,ERR=1303) "Nb Gauss Angles for ",
     &            "luminance calculation (user value):",NBMU_GAUSS_LUM
	  ELSE
	     WRITE (IDLOG,*,ERR=1303) "Nb Gauss Angles for ",
     &            "luminance calculation (default value)"
	  ENDIF
	 
C        Number of Gauss angles for the Mie computations
	  IF (NBMU_GAUSS_MIE.NE.CTE_NOT_DEFINED_VALUE_INT) THEN 
	     WRITE (IDLOG,*,ERR=1303) "Nb Gauss Angles for ",
     &            "Mie calculation (user value):",NBMU_GAUSS_MIE
	  ELSE
	     WRITE (IDLOG,*,ERR=1303) "Nb Gauss Angles for ",
     &            "Mie calculation (default value)"
	  ENDIF
          
C        User angle file for the radiance computations
	  IF (FICANGLES_USER_LUM.NE.'NO_USER_ANGLES') THEN 
	     WRITE (IDLOG,*,ERR=1303) "User's file of specific angles",
     &            " for luminance calculations (user value):",
     &            FICANGLES_USER_LUM
	  ENDIF		 
		 	 
C	  User angle file for the Mie computations
	  IF (FICANGLES_USER_MIE.NE.'NO_USER_ANGLES') THEN 
	     WRITE (IDLOG,*,ERR=1303) "User's file of specific angles",
     &            " for MIE calculations (user value):",
     &            FICANGLES_USER_MIE
	  ENDIF		  	 
          
	  WRITE (IDLOG,*,ERR=1303) " "
          
      ENDIF       ! End of Trace

      
C --------------------------------------------------------------  
C --- Check of the arguments of the aerosol model
C --------------------------------------------------------------  
C
C ---- Existence and validity test of the mandatory parameters
   	
C     Case AOT_REF > 0 :
      IF (( AOT_REF.GT.0.D+00).AND.
     & ( FICUSER_AER.EQ."NO_USER_AEROSOLS" )) THEN
      
C        Check of the definition of the aerosol parameters
C        Aerosol model
	  IF ((IMOD_AER.LT.0).OR.(IMOD_AER.GT.5)) GOTO 2305

C        Option to remove the truncation
         IF ((ITRONC_AER.LT.0).OR.(ITRONC_AER.GT.1)) GOTO 23141
	 
C        ----------------------------------------------------
C        Case of a mono-modal aerosol model: IMOD_AER = 0
C        ----------------------------------------------------
	  IF (IMOD_AER.EQ.0) THEN 

C            Real and imaginary parts of the refractive index
	     IF ((RN_WA.EQ.CTE_NOT_DEFINED_VALUE_DBLE).OR.
     &          (IN_WA.EQ.CTE_NOT_DEFINED_VALUE_DBLE)) GOTO 2309

C           Validity of the imaginary part
	     IF (IN_WA.GT.0.D+00) GOTO 2310

C           Validity of the size distribution type
            IF (IGRANU.EQ.CTE_NOT_DEFINED_VALUE_INT) GOTO 2311

	     IF ((IGRANU.LT.1).OR.(IGRANU.GT.2)) GOTO 2312	 

	     ! Check of the retrieval of the parameters of the LND model
	     IF (IGRANU.EQ.1) THEN
	         IF (LND_RADIUS_MMD_AER.EQ.
     &              CTE_NOT_DEFINED_VALUE_DBLE) GOTO 23131
		  IF (LND_LNVAR_MMD_AER.EQ.
     &              CTE_NOT_DEFINED_VALUE_DBLE) GOTO 23131

                VARGRANU1_MMD_AER = LND_RADIUS_MMD_AER
	         VARGRANU2_MMD_AER = LND_LNVAR_MMD_AER
	         VARGRANU3_MMD_AER = CTE_NOT_DEFINED_VALUE_DBLE
	     ENDIF     

	     ! Check of the retrieval of the mandatory parameters of the Junge model
	     IF (IGRANU.EQ.2) THEN
	         IF (JD_SLOPE_MMD_AER.EQ.
     &              CTE_NOT_DEFINED_VALUE_DBLE) GOTO 23132
		  IF (JD_RMIN_MMD_AER.EQ.
     &              CTE_NOT_DEFINED_VALUE_DBLE) GOTO 23132

                VARGRANU1_MMD_AER = JD_RMIN_MMD_AER
	         VARGRANU2_MMD_AER = JD_SLOPE_MMD_AER
	         VARGRANU3_MMD_AER = JD_RMAX_MMD_AER
	     ENDIF	    

C	     Case where the reference wavelength differs from the simulation one
C	     ==> Check of the retrieval of the value of the refractive index
C	         for the reference wavelength.
	      IF (WA_SIMU.NE.WAREF_AOT) THEN
	         IF((RN_WAREF.EQ.CTE_NOT_DEFINED_VALUE_DBLE).OR.
     &	            (IN_WAREF.EQ.CTE_NOT_DEFINED_VALUE_DBLE)) GOTO 2314
	      ELSE
                RN_WAREF = RN_WA
                IN_WAREF = IN_WA
	     ENDIF 
	 		 	
         ENDIF
	 
C         ----------------------------------------------------
C         Case of a WMO aerosol model: IMOD_AER = 1
C         ----------------------------------------------------	 
 	  
	  IF (IMOD_AER.EQ.1) THEN 	!WMO case

	     ! Check of the validity of the parameter IModele
	     ! -----------------------------------------	
	     IF (IModele_WMO.EQ.CTE_NOT_DEFINED_VALUE_INT) GOTO 2315	

	     ! Check of the retrieval of the choice of WMO model
	     IF ((IModele_WMO.LT.1).OR.(IModele_WMO.GT.CTE_NB_Modeles))
     &         GOTO 2316	

	     ! If free user model, check of the retrieval
	     ! of the proportions of the 4 components
	     IF (IModele_WMO.EQ.4) THEN     
	         IF ((C_WMO_DL.EQ.CTE_NOT_DEFINED_VALUE_DBLE).OR.
     &		      (C_WMO_WS.EQ.CTE_NOT_DEFINED_VALUE_DBLE).OR.
     &		      (C_WMO_OC.EQ.CTE_NOT_DEFINED_VALUE_DBLE).OR.
     &		      (C_WMO_SO.EQ.CTE_NOT_DEFINED_VALUE_DBLE)) GOTO 2317
	     ENDIF !Case of a user WMO model
	  
	  ENDIF !End of case IMOD_AER = 1 (WMO)
	  
C        ----------------------------------------------------
C        Case of an S&F aerosol model: IMOD_AER = 2
C        ----------------------------------------------------
        
	  IF (IMOD_AER.EQ.2) THEN !     S&F case
         
	     ! Check of the retrieval of the choice of S&F model
	     IF (IModele_SF.EQ.CTE_NOT_DEFINED_VALUE_INT) GOTO 2318
	     
	     ! Check of the retrieval of the value of the relative humidity
	     IF (RH.EQ.CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2319

            !   Check of the validity of the parameters IMODELE and RH
            IF ((IModele_SF.LT.1).OR.(IModele_SF.GT.CTE_NB_Modeles))
     &         GOTO 2320        !Error on the model choice
            IF ((RH.LT.0.D-00).OR.(RH.GT.99.D-00)) GOTO 2321 !Error on the % humidity

	  ENDIF
         !End of case IMOD_AER = 2 (S&F)

C        ----------------------------------------------------
C        Case of a bimodal LND aerosol model: IMOD_AER = 3
C        ----------------------------------------------------	  
         !Bi-LND case
	  IF (IMOD_AER.EQ.3) THEN

             !Check of the retrieval of the parameters specific to the bimodal LND computation
             IF (MODE_PARAM_BILND.EQ.
     &           CTE_NOT_DEFINED_VALUE_INT)  GOTO 2322

             !check of the validity of the parameter MODE_PARAM_BILND
             IF ((MODE_PARAM_BILND.LT.1).OR.(MODE_PARAM_BILND.GT.2)) 
     &         GOTO 2323
  



	     
             !Retrieval of the information on the description of the mixture
	      !    Case 1: Use of the user volume concentrations
             !    --> Check that the user parameters are given for the option
             IF (MODE_PARAM_BILND.EQ.1) THEN          
	         IF (USER_CV_COARSE.EQ.
     &              CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2324
	         IF (USER_CV_FINE.EQ.
     &              CTE_NOT_DEFINED_VALUE_DBLE)   GOTO 2325  
             ENDIF
    
            !    Case 2: Use of the ratio AOT_coarse / AOT_total
            !    --> Check that all the parameters are given for the option
            IF (MODE_PARAM_BILND.EQ.2) THEN 
	         IF (RTAUct_WAREF.EQ.
     &              CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2326
            ENDIF
   
            !Retrieval of the refractive index and size distribution information of the 2 LND
            !-----------------------------------------------------------------         	 
            IF (BMD_CM_MRWA.EQ.CTE_NOT_DEFINED_VALUE_DBLE)   GOTO 2327
            IF (BMD_CM_MIWA.EQ.CTE_NOT_DEFINED_VALUE_DBLE)   GOTO 2327 
            IF (BMD_CM_RMODAL.EQ.CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2327
            IF (BMD_CM_VAR.EQ.CTE_NOT_DEFINED_VALUE_DBLE)    GOTO 2327 
            IF (BMD_FM_MRWA.EQ.CTE_NOT_DEFINED_VALUE_DBLE)   GOTO 2328
            IF (BMD_FM_MIWA.EQ.CTE_NOT_DEFINED_VALUE_DBLE)   GOTO 2328 
            IF (BMD_FM_RMODAL.EQ.CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2328
            IF (BMD_FM_VAR.EQ.CTE_NOT_DEFINED_VALUE_DBLE)    GOTO 2328            
     
            !Retrieval of the refractive index information of the 2 LND for WAREF
            !-----------------------------------------------------------------     
            IF (MODE_PARAM_BILND.EQ.2) THEN 
	     
	        IF (WA_SIMU.NE.WAREF_AOT) THEN
	           IF (BMD_CM_MRWAREF.EQ.
     &                CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2329
	           IF (BMD_CM_MIWAREF.EQ.
     &                CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2329
	           IF (BMD_FM_MRWAREF.EQ.
     &                CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2329
	           IF (BMD_FM_MIWAREF.EQ.
     &                CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2329   

	        ELSE !If waref=wa ==> The refractive index information for wa is required
		      !                and assigned to waref
		    BMD_CM_MRWAREF = BMD_CM_MRWA
		    BMD_CM_MIWAREF = BMD_CM_MIWA
		    BMD_FM_MRWAREF = BMD_FM_MRWA
		    BMD_FM_MIWAREF = BMD_FM_MIWA  

	        ENDIF ! End of test if Wa != WaRef
	        
            ENDIF !End of case of a bi-LND mixture described by an AOT ratio
	           
	  ENDIF


C         ----------------------------------------------------
C         Case of an aerosol model provided by the user: IMOD_AER = 4
C         ----------------------------------------------------	  
	  IF (IMOD_AER.EQ.4) THEN 	!Case of external functions

             ! Check of the retrieval of the name of the external data file
	     IF (FICEXTDATA_AER.EQ.
     &          "NO_USER_AEROSOLS_PHAZE_FCT") GOTO 2330

	     ! Check that the reference wavelength for the aerosols is
	     ! the same as the simulation wavelength of the radiance field
	     IF (WA_SIMU.NE.WAREF_AOT) GOTO 2331

	  ENDIF !End of case IMOD_AER = 4 (external phase functions)

	  
C         ---------------------------------------------------------------
C         Case of a mixture of aerosol models defined by the user: IMOD_AER = 5
C         ---------------------------------------------------------------	 
 	  
	  IF (IMOD_AER.EQ.5) THEN 	!Case of a mixture of models by the user
	  
	     ! Check of the retrieval of the name of the mixture definition file
	     IF (FICMIXTURE_AER.EQ."NO_USER_AEROSOLS_MIXTURE") GOTO 2340
  
	  ENDIF !End of case IMOD_AER = 5 (user mixture)

	  	  	  	  	             
      ENDIF ! END of case test on AOT > 0 and no user aerosol file



       
C     Case AOT_REF > 0 and DEF_FICUSER_AER = 1
C     => Use of the user file:
C     Error case: the result file must not be defined
C                    the simulation and reference wavelengths must be the same

      IF (( AOT_REF.GT.0.D+00)
     &      .AND.(FICUSER_AER.NE."NO_USER_AEROSOLS" )) THEN
          
	  ! Check that the reference wavelength for the aerosols is
	  ! the same as the simulation wavelength of the radiance field
         IF (WA_SIMU.NE.WAREF_AOT) GOTO 2350

	  ! Check that the result file of the IOP computations is not defined by the user
	   IF (FICGRANU.NE.CTE_DEFAULT_FICGRANU) GOTO 2351
	  
      ENDIF
      
C ---- Display of the parameters used in the trace file

      IF (TRACE) THEN

         WRITE (IDLOG,*,ERR=1303) " "
         WRITE (IDLOG,*,ERR=1303) "************** ",
     &         "Aerosols parameters ***************"
         WRITE (IDLOG,*,ERR=1303) " "
	  
         WRITE (IDLOG,*,ERR=1303) "***** Mandatory Parameters "
	  
C        Reference AOT
         WRITE (IDLOG,*,ERR=1303) 
     &         "AOT for the reference wavelength:",AOT_REF

C	  Case AOT > 0
         IF (( AOT_REF.GT.0.D+00)
     &         .AND.(FICUSER_AER.EQ."NO_USER_AEROSOLS")) THEN

C           Reference wavelength
            WRITE (IDLOG,*,ERR=1303) 
     &            "Reference Wavelength (microns): ",WAREF_AOT
      
C           Type of aerosol model
            WRITE (IDLOG,*,ERR=1303) "Aerosol Model Type: ",IMOD_AER


C	     Case of a mono-modal model
            IF (IMOD_AER.EQ.0) THEN

                WRITE (IDLOG,*,ERR=1303)
     &               "Mono-modal aerosol Model"
		   
C         	  Real part of the refractive index at WA_SIMU
                WRITE (IDLOG,*,ERR=1303)
     &               "Real part of the refractive index: ",RN_WA
     
C         	  Imaginary part of the refractive index at WA_SIMU
                WRITE (IDLOG,*,ERR=1303)
     &               "Imaginary part of the refractive index: ",IN_WA
     
C		  Size distribution type
                WRITE (IDLOG,*,ERR=1303)
     &               "Particle size distribution model: ",IGRANU

C		  Size distribution parameters
                IF (IGRANU.EQ.1) THEN
                   WRITE (IDLOG,*,ERR=1303) "LND parameters: "
                   WRITE (IDLOG,*,ERR=1303) 
     &                  "Modal radius (microns): ", VARGRANU1_MMD_AER 
                   WRITE (IDLOG,*,ERR=1303) 
     &                  "Variance (microns): ", VARGRANU2_MMD_AER 
                ELSE
                   WRITE (IDLOG,*,ERR=1303) "JUNGE's law parameters: "
                   WRITE (IDLOG,*,ERR=1303) 
     &                  "Slope: ", VARGRANU1_MMD_AER 
                   WRITE (IDLOG,*,ERR=1303) 
     &                  "Min radius (microns): ", VARGRANU2_MMD_AER 
                   WRITE (IDLOG,*,ERR=1303) 
     &                  "Max radius (microns): ", VARGRANU3_MMD_AER 
                ENDIF

	     	  IF (WA_SIMU.NE.WAREF_AOT) THEN
	     
C                  Real part of the refractive index at WA_REF
                   WRITE (IDLOG,*,ERR=1303)
     &                    "Real part of the refractive index",
     &                    " (ref wvl): ",RN_WAREF
     
C         	     Imaginary part of the refractive index at WA_SIMU
                   WRITE (IDLOG,*,ERR=1303)
     &                    "Imaginary part of the refractive index ",
     &                    "(ref wvl): ",IN_WAREF
                ENDIF
                  
            ENDIF ! End IMOD_AER.EQ.0
		
		
C	     Case of a WMO model
            IF (IMOD_AER.EQ.1) THEN

	  	  WRITE (IDLOG,*,ERR=1303) "WMO aerosol Model"
		
C         	  WMO model
                WRITE (IDLOG,*,ERR=1303)
     &                 "WMO Model type: ",IModele_WMO
                  
	     	  IF (IModele_WMO.EQ.4) THEN
	     
C         	     Volume percentage of the DL particles of the WMO model
                   WRITE (IDLOG,*,ERR=1303)
     &                    "DL particules ratio",C_WMO_DL
     
C         	     Volume percentage of the WS particles of the WMO model
                   WRITE (IDLOG,*,ERR=1303)
     &                    "WS particules ratio",C_WMO_WS

C         	     Volume percentage of the SO particles of the WMO model
                   WRITE (IDLOG,*,ERR=1303)
     &                    "SO particules ratio",C_WMO_SO

C         	     Volume percentage of the OC particles of the WMO model
                   WRITE (IDLOG,*,ERR=1303)
     &                    "OC particules ratio",C_WMO_OC

		  ENDIF
     
            ENDIF            ! End IMOD_AER.EQ.1

C	     Case of an SF model
            IF (IMOD_AER.EQ.2) THEN

	  	  WRITE (IDLOG,*,ERR=1303) "S&F aerosol Model"
					
C               S&F model
                WRITE (IDLOG,*,ERR=1303)
     &                 "S&F Model type: ",IModele_SF
          
C         	  Relative humidity
                WRITE (IDLOG,*,ERR=1303)
     &                 "Relative Humidity: ",RH

            ENDIF  ! End IMOD_AER.EQ.2
		

C	     Case of a bimodal model

            IF (IMOD_AER.EQ.3) THEN

	  	  WRITE (IDLOG,*,ERR=1303)"Bi-modal aerosol Model"
		
C         	  Bimodal computation mode
                WRITE (IDLOG,*,ERR=1303)
     &                 "Calculation mode: ",MODE_PARAM_BILND
          
		  IF (MODE_PARAM_BILND.EQ.1) THEN        

C	            volume concentration of the LND coarse
	  	     WRITE (IDLOG,*,ERR=1303)
     &                 "User's volumic concentration of LND coarse: ",
     &                    USER_CV_COARSE

C		     volume concentration of the LND fine
	  	     WRITE (IDLOG,*,ERR=1303)
     &                    "User's volumic concentration of LND fine: ",
     &                    USER_CV_FINE
 			  
            	  ENDIF

		  IF (MODE_PARAM_BILND.EQ.2) THEN        

C		     Ratio AOT coarse / AOT total
	  	     WRITE (IDLOG,*,ERR=1303)
     &                    "User's ratio (coarse AOT / total AOT) : ",
     &                    RTAUct_WAREF
 			  
            	  ENDIF

C         	  Real part of the refractive index of the Coarse mode
                WRITE (IDLOG,*,ERR=1303)
     &                 "Real part of the refractive index ",
     &                 "of Coarse model: ",BMD_CM_MRWA
     
C         	  Imaginary part of the refractive index of the Coarse mode
                WRITE (IDLOG,*,ERR=1303)
     &                 "Imaginary part of the refractive index ",
     &                 "of Coarse model: ",BMD_CM_MIWA

C         	  Modal radius of the Coarse mode
                WRITE (IDLOG,*,ERR=1303)
     &          "Coarse modal radius of Coarse model: ",BMD_CM_RMODAL
     
C         	  LND variance for the Coarse mode
                WRITE (IDLOG,*,ERR=1303)
     &          "variance LND of the Coarse Mode ",BMD_CM_VAR

C         	  Real part of the refractive index of the Fine mode
                 WRITE (IDLOG,*,ERR=1303)
     &                 "Real part of the refractive index",
     &                 " of Fine model: ",BMD_FM_MRWA
     
C         	  Imaginary part of the refractive index of the Fine mode
                WRITE (IDLOG,*,ERR=1303)
     &                 "Imaginary part of the refractive index ",
     &                 "of Fine model: ",BMD_FM_MIWA

C         	  Modal radius of the Fine mode
                WRITE (IDLOG,*,ERR=1303)
     &          "Coarse modal radius of Fine model: ",BMD_FM_RMODAL
     
C         	  LND variance for the Fine mode
                WRITE (IDLOG,*,ERR=1303)
     &          "variance LND of the Fine Mode ",BMD_FM_VAR
     
                IF (WA_SIMU.NE.WAREF_AOT) THEN

C                  Real part of the refractive index of the Coarse mode
                   WRITE (IDLOG,*,ERR=1303)
     &		     "Real part of the refractive index",
     &		     "of Coarse model (ref wvl): ",BMD_CM_MRWAREF
     
C         	     Imaginary part of the refractive index of the Coarse mode
                   WRITE (IDLOG,*,ERR=1303)
     &		     "Imaginary part of the refractive index ",
     &		     "of Coarse model (ref wvl): ",BMD_CM_MIWAREF


C         	     Real part of the refractive index of the Fine mode
                   WRITE (IDLOG,*,ERR=1303)
     &		     "Real part of the refractive index",
     &		     "of Fine model (ref wvl): ",BMD_FM_MRWAREF
     
C         	     Imaginary part of the refractive index of the Fine mode
                   WRITE (IDLOG,*,ERR=1303)
     &		     "Imaginary part of the refractive index ",
     &		     "of Fine model (ref wvl): ",BMD_FM_MIWAREF

		  ENDIF
			
	     ENDIF  ! End IMOD_AER.EQ.3
	 

C	     Case of user data
	     IF (IMOD_AER.EQ.4) THEN

	  	  WRITE (IDLOG,*,ERR=1303)
     &		  "User gives aerosol radiave properties"
					
C         	  User file
                WRITE (IDLOG,*,ERR=1303)
     &		  "User's filename: ",FICEXTDATA_AER
     
	     ENDIF    



C	     Case of a mixture model defined by the user
	     IF (IMOD_AER.EQ.5) THEN

	  	  WRITE (IDLOG,*,ERR=1303)
     &		  "User defines the mixture of aerosol models"
					
C         	  User file
                WRITE (IDLOG,*,ERR=1303)
     &		  "User's filename : ",FICMIXTURE_AER
     
	     ENDIF   

		
	  ENDIF ! end of test AOT > 0 and DEF_FICUSER_AER = 0
	 

C	  Case of an aerosol file provided by the user
	  IF (( AOT_REF.GT.0.D+00)
     &       .AND.(FICUSER_AER.NE.'NO_USER_AEROSOLS')) THEN

	     WRITE (IDLOG,*,ERR=1303) " "	  
	     WRITE (IDLOG,*,ERR=1303) 
     &      "***** Use of an aerosol User File : ",FICUSER_AER
	     WRITE (IDLOG,*,ERR=1303) "=> No Aerosol Processing"
	  
	  ELSE	  

C         Result file

	     WRITE (IDLOG,*,ERR=1303) " "	  
            WRITE (IDLOG,*,ERR=1303) "***** Optional Parameters "
	  
	     IF (FICGRANU.NE.CTE_DEFAULT_FICGRANU) THEN 
	     	  WRITE (IDLOG,*,ERR=1303) "User's filename of the ",
     &       	  "result file (user value):", FICGRANU
	     ELSE
	     	  WRITE (IDLOG,*,ERR=1303) "User's filename of the ",
     &       	  "result file (default value):", FICGRANU  
	     ENDIF		  	 
 
   	  ENDIF  ! end of test on the use of a user file

	  WRITE (IDLOG,*,ERR=1303) " "	  
	  WRITE (IDLOG,*,ERR=1303) "***** Optional Parameters "	  
	  WRITE (IDLOG,*,ERR=1303) "Phase function troncature option:",
     &       	    ITRONC_AER
      ENDIF  ! end of trace writing
          

      
C ------------------------------------------------------------------        
C --- Check of the arguments for the surface computations
C ------------------------------------------------------------------  
C  
C ---- Existence and validity test of the mandatory parameters

C     Argument RHO (surface albedo) mandatory, check of its value
      IF (RHO.EQ.CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2401        
      IF (RHO.LT.0.D-00) GOTO 2402
 
C     Argument ISURF (surface type) mandatory, check of its value
      IF (ISURF.EQ.CTE_NOT_DEFINED_VALUE_INT) GOTO 2403
     
C     Argument ISURF mandatory, check of its value
      IF((ISURF.NE.0).AND.(ISURF.NE.1).AND.(ISURF.NE.2)
     &     .AND.(ISURF.NE.3).AND.(ISURF.NE.4).AND.(ISURF.NE.5)
     &     .AND.(ISURF.NE.6).AND.(ISURF.NE.7)) GOTO 2404
            
C     Argument SURF_IND required if ISURF = 1, 2, 4, 5, 6 or 7
      IF ( (SURF_IND.EQ.CTE_NOT_DEFINED_VALUE_DBLE)
     &     .AND.((ISURF.EQ.1).OR.(ISURF.EQ.2)
     &     .OR.(ISURF.EQ.4).OR.(ISURF.EQ.5).OR.(ISURF.EQ.6)
     &     .OR.(ISURF.EQ.7)) ) GOTO 2405  

C     GLITTER case => Argument WIND mandatory
      IF (ISURF.EQ.1) THEN
         IF (WIND.EQ.CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2406
         IF (WIND.LT.0.D-00) GOTO 24061
      ENDIF	            
      
C     Case of the Roujean BRDF model => parameters K0, K1 and K2 mandatory
      IF (ISURF.GE.3) THEN
         IF ( (K0_ROUJEAN.EQ.CTE_NOT_DEFINED_VALUE_DBLE)
     &        .OR.(K1_ROUJEAN.EQ.CTE_NOT_DEFINED_VALUE_DBLE)
     &        .OR.(K2_ROUJEAN.EQ.CTE_NOT_DEFINED_VALUE_DBLE) ) THEN
             GOTO 2407
         ENDIF 	 
      ENDIF	                     

C     Case of the NADAL BPDF => parameters ALPHA_NADAL and BETA_NADAL must be defined
      IF (ISURF.EQ.6) THEN

         !Removal of the possibility to use the Nadal model (not applicable to the SOS formalism)
         IF (TRACE) THEN
            WRITE (IDLOG,*,ERR=1303) " "
            WRITE (IDLOG,*,ERR=1303) "WARNING:"
            WRITE (IDLOG,*,ERR=1303) 
     &           "The Nadal's BPDF model is not supported " 
            WRITE (IDLOG,*,ERR=1303) "==> Select another surface model"
         ENDIF 

         WRITE(6,*) "  "
         WRITE(6,*) "  "
         WRITE(6,*) "  The Nadal's BPDF model is not supported"
         WRITE(6,*) "  ==> Select another surface model"
         GOTO 9999
        
         IF ((ALPHA_NADAL.EQ.CTE_NOT_DEFINED_VALUE_DBLE)
     &   .OR.(BETA_NADAL.EQ.CTE_NOT_DEFINED_VALUE_DBLE)) GOTO 2408

C        Validity of the coefficients ALPHA_NADAL and BETA_NADAL
         IF((ALPHA_NADAL.LT.0.).OR.(ALPHA_NADAL.GE.10.)) GOTO 2409
	  IF((BETA_NADAL.LT.0.).OR.(BETA_NADAL.GE.100.)) GOTO 2410
	 
      ENDIF 


C     Case of the Maignan BPDF model => Argument COEF_C_MAIGNAN mandatory
      IF (ISURF.EQ.7) THEN
         IF (COEF_C_MAIGNAN.EQ.CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2411        
      ENDIF


C ---- Display of the parameters used in the trace file

      IF (TRACE) THEN

	  WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) "************** ",
     &       "Surface parameters ***************"
         WRITE (IDLOG,*,ERR=1303) " "

C        Surface albedo
	  WRITE (IDLOG,*,ERR=1303) "Surface Albedo: ",RHO

C        Surface model
	  WRITE (IDLOG,*,ERR=1303) "Surface Type: ",ISURF
	  
C	  Refractive index
         IF (SURF_IND.NE.CTE_NOT_DEFINED_VALUE_DBLE) THEN
	      WRITE (IDLOG,*,ERR=1303) "Refractive Index ",SURF_IND
	  ENDIF
	  
C	  GLITTER case
      	  IF (ISURF.EQ.1) THEN
	      WRITE (IDLOG,*,ERR=1303) "GLITTER case"
	      WRITE (IDLOG,*,ERR=1303) "Wind: ",WIND
      	  ENDIF	            

      	  IF (ISURF.EQ.2) THEN
	      WRITE (IDLOG,*,ERR=1303)
     &       "Fresnel's reflexion on flat water"
      	  ENDIF
  
C	  Case of the Roujean BRDF
      	  IF (ISURF.GE.3) THEN
	      WRITE (IDLOG,*,ERR=1303) "Roujean BRDF parameters"
	      WRITE (IDLOG,*,ERR=1303) "K0: ",K0_ROUJEAN           
	      WRITE (IDLOG,*,ERR=1303) "K1: ",K1_ROUJEAN 
	      WRITE (IDLOG,*,ERR=1303) "K2: ",K2_ROUJEAN 
         ENDIF	                     

      	  IF (ISURF.EQ.4) THEN
	      WRITE (IDLOG,*,ERR=1303) "Rondeaux's BPDF model"
      	  ENDIF

      	  IF (ISURF.EQ.5) THEN
	      WRITE (IDLOG,*,ERR=1303) "Breon's BPDF model"
      	  ENDIF

C     	  Case of the NADAL BPDF => parameters ALPHA_NADAL and BETA_NADAL
      	  IF (ISURF.EQ.6) THEN
	      WRITE (IDLOG,*,ERR=1303) "Nadal's BPDF model"
             WRITE (IDLOG,*,ERR=1303) "Parameters:"		
             WRITE (IDLOG,*,ERR=1303) "  ALPHA_NADAL: ",ALPHA_NADAL
             WRITE (IDLOG,*,ERR=1303) "  BETA_NADAL: ",BETA_NADAL	 
      	  ENDIF 	 

      	  IF (ISURF.EQ.7) THEN
	      WRITE (IDLOG,*,ERR=1303) "Maignan's BPDF model"
             WRITE (IDLOG,*,ERR=1303) "  C coef: ",COEF_C_MAIGNAN 
      	  ENDIF
	  
      ENDIF  ! End of trace

C ------------------------------------------------------------------        
C --- Check of the arguments for the profile computations
C ------------------------------------------------------------------  
C  
C ---- Existence and validity test of the mandatory parameters

C     Rayleigh optical thickness
      IF ((TR.NE.CTE_NOT_DEFINED_VALUE_DBLE)
     &     .AND.(TR.LT.0.D-00)) GOTO 2502
      
C     Rayleigh scale height

      IF ((TR.NE.CTE_NOT_DEFINED_VALUE_DBLE)
     &     .AND.(TR.GT.0.)) THEN
         IF (HR.EQ.CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2503  
         IF (HR.LE.0.D-00) GOTO 2504
      ENDIF
   
C     Type of profile for the aerosol scattering
      IF (IPROFIL.EQ.CTE_NOT_DEFINED_VALUE_INT) GOTO 2505   
      IF ((IPROFIL.LT.1).OR.(IPROFIL.GT.CTE_NB_PROFIL)) GOTO 2506

C     Aerosol scale height if IPROFIL=1
      IF (IPROFIL.EQ.1) THEN
         IF ((HA.EQ.CTE_NOT_DEFINED_VALUE_DBLE)
     &        .AND.(AOT_REF.GT.0.)) GOTO 2507
         IF ((HA.LE.0.D-00).AND.(AOT_REF.GE.0.))GOTO 2508	  
      ENDIF
       
C     Min and max altitudes of the profile if IPROFIL=2
      IF (IPROFIL.EQ.2) THEN
         IF ((ZMIN.EQ.CTE_NOT_DEFINED_VALUE_DBLE)
     &      .OR.(ZMAX.EQ.CTE_NOT_DEFINED_VALUE_DBLE)) GOTO 2509
      ENDIF

C     Type of gaseous absorption profile
            
      IF (ABSPROFIL.EQ.CTE_NOT_DEFINED_VALUE_INT) GOTO 2510      
      IF ((ABSPROFIL.LT.0).OR.(ABSPROFIL.GT.7)) GOTO 2511

C     User file of the absorption profile if ABSPROFIL = 0
	  
      IF (ABSPROFIL.EQ.0) THEN
      	  IF (FICABSPROFIL.EQ."NO_USER_ABS_PROFILE_FILE") GOTO 2512
      ENDIF
     
C     No gaseous absorption if IPROFIL=2
      IF ((IPROFIL.EQ.2).AND.(ABSPROFIL.NE.7)) GOTO 2513
      
C     Spectral resolution to be given if ABSPROFIL != 7
C     Radiance computation mode to be given if ABSPROFIL != 7
      IF (ABSPROFIL.NE.7) THEN  ! Case with gaseous absorption

      	  IF (NUSTEP.EQ.CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2514
		
         IF ((NUSTEP.NE.1).AND.(NUSTEP.NE.5)
     &      .AND.(NUSTEP.NE.10)) GOTO 25141

         IF (IMODE_CKD_CALCUL.EQ.CTE_NOT_DEFINED_VALUE_INT) GOTO 2515
                
      ELSE
          IMODE_CKD_CALCUL = 2	!Default mode for the case without gaseous absorption
      ENDIF


      
C ---- Display of the parameters used in the trace file
      
      IF (TRACE) THEN

	  WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) "************** ",
     &       "Profil parameters ***************"
         WRITE (IDLOG,*,ERR=1303) " "

C        Rayleigh optical thickness
	  WRITE (IDLOG,*,ERR=1303) "Molecular Optical Thickness: ",TR

         IF (TR.EQ.CTE_NOT_DEFINED_VALUE_DBLE) THEN
             WRITE (IDLOG,*,ERR=1303)
     &       "--> Calculated by CNES formulation (Perbos, 1982)"

             WRITE (IDLOG,*,ERR=1303) "Atmospheric pressure: ",PSURF
             IF(PSURF.EQ.CTE_HT_STD_PSURF) THEN
		  WRITE (IDLOG,*,ERR=1303) "--> Default value "
             ELSE
		  WRITE (IDLOG,*,ERR=1303) "--> User value "
             ENDIF
  	  ENDIF


C        Rayleigh scale height
	  IF (HR.NE.CTE_NOT_DEFINED_VALUE_DBLE) THEN
	      WRITE (IDLOG,*,ERR=1303) "Molecular Height scale ",HR
	  ENDIF
	  
C	  Profile type
      	  WRITE (IDLOG,*,ERR=1303) "Aerosol profile type: ",IPROFIL           
  
C        Aerosol scale height
	  IF (HR.NE.CTE_NOT_DEFINED_VALUE_DBLE) THEN
	      WRITE (IDLOG,*,ERR=1303) "Aerosol Height Scale ",HA
	  ENDIF                   

C        Min and max altitudes of the profile if IPROFIL=2
      	  IF (IPROFIL.EQ.2) THEN
	      WRITE (IDLOG,*,ERR=1303) "Profil min altitude ",ZMIN
	      WRITE (IDLOG,*,ERR=1303) "Profil max altitude ",ZMAX
      	  ENDIF	  

C     	  Type of gaseous absorption profile
      	  WRITE (IDLOG,*,ERR=1303) "Absorption Profil type ",ABSPROFIL
          
C        User file of the absorption profile if ABSPROFIL = 0
          IF (ABSPROFIL.EQ.0) THEN
      	      WRITE (IDLOG,*,ERR=1303) "Absorption user Profile file ",
     &	      FICABSPROFIL
          ENDIF
	  
      ENDIF  ! End of trace
     
   	
C -------------------------------------------------------------------  
C --- Check of the parameters specific to the Successive Orders computations
C     and to the production of the result files
C -------------------------------------------------------------------  
C  
C ---- Existence and validity test of the mandatory parameters

C     Check of the maximum order of interaction
      IF (IGMAX.LT.1) GOTO 2604

C     Viewing option: mandatory
      IF (ITRPHI.EQ.CTE_NOT_DEFINED_VALUE_INT ) GOTO 2605

C     Viewing option: mandatory
      IF ((ITRPHI.LT.1).OR.(ITRPHI.GT.2)) GOTO 2606

      !If ITRPHI=1 is chosen: the relative azimuth must be defined
      IF (ITRPHI.EQ.1) THEN 
         IF (PHIOS.EQ.CTE_NOT_DEFINED_VALUE_DBLE) GOTO 2607
      ENDIF
          	
      !If ITRPHI=2 is chosen: the azimuth step must be defined
      IF (ITRPHI.EQ.2) THEN 
         IF (PAS_PHI.EQ.CTE_NOT_DEFINED_VALUE_INT) GOTO 2608
         IF (PAS_PHI.LE.0) GOTO 2609
      ENDIF
      
      
C ---- Existence and validity test of the optional parameters

C     Option to cut the polarisation
      IF ((IPOLAR.LT.0).OR.(IPOLAR.GT.1)) GOTO 2610    

C     Option for the direct reflection of the solar beam by the surface
      IF ((IDIRSURF.LT.0).OR.(IDIRSURF.GT.1)) GOTO 2612
     
C     Output level of the radiance simulations: optional
C     If not defined: standard output:
C        Field emerging at TOA, downward field at GROUND: ZOUT = -1
      IF (((ZOUT.LT.0.D+00).AND.(ZOUT.NE.-1.D+00))
     &     .OR.(ZOUT.GT.CTE_TOA_ALT)) THEN
         GOTO 2611
      ENDIF      
      
      IF (TRACE) THEN

	  WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) 
     &   "************** SOS parameters ***************"
         WRITE (IDLOG,*,ERR=1303) " "

	  WRITE (IDLOG,*,ERR=1303) "***** Mandatory Parameters "
	  
C        Molecular depolarisation factor
	  WRITE (IDLOG,*,ERR=1303) "Molecular depolarization factor: ",
     &	  RON

C        Maximum order of interaction:
	  WRITE (IDLOG,*,ERR=1303) "Maximum Number of interactions: ",
     &	  IGMAX

C	  Viewing option
	  WRITE (IDLOG,*,ERR=1303) "Viewing mode: ",ITRPHI

C	  Case ITRPHI=1 => Display of the relative azimuth
	  IF (ITRPHI.EQ.1) THEN
	  	WRITE (IDLOG,*,ERR=1303) "Relative azimut: ",PHIOS
	  ENDIF	  

C	  Case ITRPHI=2 => Display of the sampling step in phi
	  IF (ITRPHI.EQ.2) THEN
	      WRITE (IDLOG,*,ERR=1303) "Azimut sampling step: ",PAS_PHI
	  ENDIF		  

	  WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) "***** Optional Parameters "	  

C        Option to cut the polarisation: optional
	  IF (IPOLAR.EQ.0) THEN
	     WRITE (IDLOG,*,ERR=1303) 
     &	     "Simulation without polarization"
	  ELSE
	     WRITE (IDLOG,*,ERR=1303) 
     &	     "Simulation including polarization"
	  ENDIF	  

C        Output level of the radiances: optional
         IF (ZOUT.NE.-1.D+00) THEN
	     WRITE (IDLOG,*,ERR=1303) "Output Altitude (km): ",
     &	     ZOUT
	  ENDIF	  

      ENDIF  ! End of trace
	 

C =========================================================
C --- Management of the directory tree of the output files
C =========================================================     
 
      
C --- Trace files

      !Angles trace file
      IF (FICANGLOG.NE."NO_LOG_FILE") THEN
          LFIC = INDEX(FICANGLOG,' ') - 1
          IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3001
          FICANGLOG = DIR_LOG(1:LDIRLOG)//'/'//FICANGLOG(1:LFIC)
      ENDIF 
       
      !Trace file of the aerosol Mie computations
      IF (FICMIE_LOG.NE."NO_LOG_FILE") THEN
          LFIC = INDEX(FICMIE_LOG,' ') - 1
          IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3002
          FICMIE_LOG = DIR_LOG(1:LDIRLOG)//'/'//FICMIE_LOG(1:LFIC)
      ENDIF 

      !Trace file of the aerosol phase computations
      IF (FICGRANU_LOG.NE."NO_LOG_FILE") THEN
          LFIC = INDEX(FICGRANU_LOG,' ') - 1
          IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3003
          FICGRANU_LOG = DIR_LOG(1:LDIRLOG)//'/'//FICGRANU_LOG(1:LFIC)
      ENDIF 
       
      !Trace file of the profile computations
      IF (FICPROFIL_LOG.NE."NO_LOG_FILE") THEN
          LFIC = INDEX(FICPROFIL_LOG,' ') - 1
          IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3004
          FICPROFIL_LOG = DIR_LOG(1:LDIRLOG)//'/'
     &	  		 //FICPROFIL_LOG(1:LFIC)
          TRACEPROF=.TRUE.
      ELSE
          TRACEPROF=.FALSE.	  
      ENDIF     
       
      !Trace file of the surface matrix computations
      IF (FICSURF_LOG.NE."NO_LOG_FILE") THEN
          LFIC = INDEX(FICSURF_LOG,' ') - 1
          IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3005
          FICSURF_LOG = DIR_LOG(1:LDIRLOG)//'/'//FICSURF_LOG(1:LFIC)
      ENDIF     
       
      !Trace file of the SOS computations
      IF (FICSOS_LOG.NE."NO_LOG_FILE") THEN
          LFIC = INDEX(FICSOS_LOG,' ') - 1
          IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3006
          FICSOS_LOG = DIR_LOG(1:LDIRLOG)//'/'//FICSOS_LOG(1:LFIC)
	   TRACEOS=.TRUE.
      ELSE
          TRACEOS=.FALSE.	   
      ENDIF           

      IF (TRACE) THEN

	  WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) 
     &       "************** Trace filenames ***************"
         WRITE (IDLOG,*,ERR=1303) " "
	  
C        Angles trace file
         WRITE (IDLOG,*,ERR=1303) "Angle trace filename: ",FICANGLOG
	
C        Trace file of the aerosol Mie computations
         WRITE (IDLOG,*,ERR=1303) "Mie trace filename: ",FICMIE_LOG
	
C        Trace file of the aerosol phase computations
         WRITE (IDLOG,*,ERR=1303) "Aerosol trace filename: ",
     &	  FICGRANU_LOG

C        Trace file of the profile computations
         WRITE (IDLOG,*,ERR=1303) "Profil trace filename: ",
     &	  FICPROFIL_LOG
	  
C        Trace file of the surface matrix computations
         WRITE (IDLOG,*,ERR=1303) "Surface trace filename: ",
     &	  FICSURF_LOG

C        Trace file of the SOS computations
         WRITE (IDLOG,*,ERR=1303) "SOS trace filename: ",FICSOS_LOG
              	  
      ENDIF	! end of trace
                                     
C ---- Output files
       
      !Angle file for the radiances
      LFIC = INDEX(FICANGLES_RES_LUM,' ') - 1
      IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3007
      FICANGLES_RES_LUM = DIR_SOS(1:LDIRSOS)//'/'
     &                    //FICANGLES_RES_LUM(1:LFIC)
         
      !Angle file for the Mie computations
      LFIC = INDEX(FICANGLES_RES_MIE,' ') - 1
      IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3008
      FICANGLES_RES_MIE = DIR_SOS(1:LDIRSOS)//'/'
     &                    //FICANGLES_RES_MIE(1:LFIC)

        
      !Aerosol phase matrix file
      LFIC = INDEX(FICGRANU,' ') - 1
      IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3009
      FICGRANU = DIR_SOS(1:LDIRSOS)//'/'//FICGRANU(1:LFIC)
             
                
      !Binary file of the radiance in Fourier series
      IF (FICSOS_RES_BIN.NE.'NO_OUTPUT') THEN
         LFIC = INDEX(FICSOS_RES_BIN,' ') - 1
         IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3011
         FICSOS_RES_BIN = DIR_SOS(1:LDIRSOS)//'/'
     &        //FICSOS_RES_BIN(1:LFIC)
      ENDIF

      !Test of the file existence
      LFICOS=INDEX(FICSOS_RES_BIN,' ')-1
      INQUIRE(FILE=FICSOS_RES_BIN(1:LFICOS),EXIST=EX)
      IF (EX.EQV..TRUE.) THEN
       
c         WRITE(6,*) 
c     &       'File -SOS.ResBin exists from previous calculations:',
c     &       FICSOS_RES_BIN(1:LFICOS)
          
	  IF (CTE_REMOVE_PREVIOUS_BIN_FILE.EQ.0) THEN
	     GOTO 30111 
	  ELSE
c	     WRITE(6,*)'  Cte CTE_REMOVE_PREVIOUS_BIN_FILE is set to 1:'
c	     WRITE(6,*)'  The file is deleted before new computations'
            OPEN(UNIT=17,FILE=FICSOS_RES_BIN,FORM='UNFORMATTED') 
            CLOSE(17,STATUS='DELETE')
	  ENDIF   
	  
      ENDIF    
       
        
      !Transmission file
      LFIC = INDEX(FICTRANS,' ') - 1
      IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3012
      IF (FICTRANS.NE.'NO_OUTPUT') THEN
          FICTRANS = DIR_SOS(1:LDIRSOS)//'/'//FICTRANS(1:LFIC)
      ENDIF
 
      !Flux file
      LFIC = INDEX(FICFLUX,' ') - 1
      IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 30121
      IF (FICFLUX.NE.'NO_OUTPUT') THEN
          FICFLUX = DIR_SOS(1:LDIRSOS)//'/'//FICFLUX(1:LFIC)
      ENDIF

      IF (TRACE) THEN

	  WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) 
     &       "************** Output filenames ***************"
         WRITE (IDLOG,*,ERR=1303) " "
	  

C        Angle file for the radiances
         WRITE (IDLOG,*,ERR=1303) "Radiance Angles: ",FICANGLES_RES_LUM
       	
C        Angle file for the Mie computations
         WRITE (IDLOG,*,ERR=1303) "Mie Angles: ",FICANGLES_RES_MIE
	
C        Aerosol phase matrix file
         WRITE (IDLOG,*,ERR=1303) "Aerosol phase matrix: ",
     &	  FICGRANU         
	  
C        Binary file of the radiance in Fourier series
         WRITE (IDLOG,*,ERR=1303) "SOS binary output: ",FICSOS_RES_BIN

C        Transmission file
         WRITE (IDLOG,*,ERR=1303) "Transmittances: ",FICTRANS

C        Flux file
         WRITE (IDLOG,*,ERR=1303) "Flux: ",FICFLUX
	      	        	  
      ENDIF	! end of trace

C =========================================================
C --- Preparation of the angles for the processing
C =========================================================   

      WRITE(6,*) "DEFINING ANGLES FOR COMPUTATIONS"

      IF (FICANGLES_USER_LUM.EQ.'NO_USER_ANGLES') THEN
          LFIC=14
      ELSE
          LFIC=INDEX(FICANGLES_USER_LUM,' ')
      ENDIF
      FICANGLES_USER_LUM=FICANGLES_USER_LUM(1:LFIC)

      IF (FICANGLES_USER_MIE.EQ.'NO_USER_ANGLES') THEN
          LFIC=14
      ELSE
          LFIC=INDEX(FICANGLES_USER_MIE,' ')
      ENDIF
      FICANGLES_USER_MIE=FICANGLES_USER_MIE(1:LFIC)

      IF (TRACE) THEN
         WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) "ANGLES CALCULATION"
	  WRITE (IDLOG,*,ERR=1303) "------------------"
         WRITE (IDLOG,*,ERR=1303) 
     &	  "--> SOS_ANGLES is called with the parameters :"
	  WRITE (IDLOG,*,ERR=1303) "NBMU_GAUSS_LUM ",NBMU_GAUSS_LUM
	  WRITE (IDLOG,*,ERR=1303) "TETAS ",TETAS
	  WRITE (IDLOG,*,ERR=1303) "FICANGLES_USER_LUM ",
     &	  FICANGLES_USER_LUM
	  WRITE (IDLOG,*,ERR=1303) "NBMU_GAUSS_MIE ",NBMU_GAUSS_MIE
	  WRITE (IDLOG,*,ERR=1303) "FICANGLES_USER_MIE ",
     &	  FICANGLES_USER_MIE
	  WRITE (IDLOG,*,ERR=1303) "FICANGLOG ",FICANGLOG
	  WRITE (IDLOG,*,ERR=1303) "FICANGLES_RES_LUM ",
     &	  FICANGLES_RES_LUM
	  WRITE (IDLOG,*,ERR=1303) "FICANGLES_RES_MIE ",
     &	  FICANGLES_RES_MIE	  
       ENDIF

       CALL SOS_ANGLES(NBMU_GAUSS_LUM, TETAS, FICANGLES_USER_LUM,                 
     &      NBMU_GAUSS_MIE, FICANGLES_USER_MIE, 
     &      FICANGLOG,
     &      FICANGLES_RES_LUM, FICANGLES_RES_MIE, IER)
     
      IF (IER.NE.0) THEN ! error case
      		
         IF (TRACE) THEN
	    WRITE (IDLOG,*,ERR=1303) 
     &	    "   ==> Error in subroutine SOS_ANGLES"
 	  ENDIF
	  GOTO 4000   
	  
      ENDIF
     
C =========================================================
C --- Preparation of the aerosol phase matrix
C =========================================================    
      

      WRITE(6,*) "AEROSOLS RADIATIVE PROPERTIES COMPUTATION"
 
      IF (TRACE) THEN
         WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) 
     &	  "-----------------------------------------"
	  WRITE (IDLOG,*,ERR=1303)
     &	  "AEROSOLS RADIATIVE PROPERTIES COMPUTATION"
	  WRITE (IDLOG,*,ERR=1303) 
     &	  "-----------------------------------------"	  
      ENDIF !End of Trace


C --- If the optical thickness is null or if the reference wavelength
C --- is the same as the simulation one, then the routine
C --- SOS_AEROSOLS is run only once. Otherwise, it must be run twice
C --- and the optical thickness at the simulation wavelength must be determined.

C     For the first run of the routine, TA
C     (simulation optical thickness) is set to AOT_REF
      TA = AOT_REF

C     Determination of the wavelength for the run of the routine
      		            
      IF (AOT_REF.EQ.0.D+00) THEN  ! Case AOT_REF null
      
         ! In this case, the computations are made directly at the simulation wavelength
	  WA = WA_SIMU 	
	  
      ELSE ! Case AOT_REF not null

	  ! The processing is made at the reference wavelength
         ! A second processing will then be made at the simulation
	  ! wavelength if it differs from the reference one
	  WA = WAREF_AOT 
			  
      ENDIF

C --- Management of the trace (display of the computation parameters)
      
      IF (TRACE) THEN
	 
         IF (( FICUSER_AER.EQ.'NO_USER_AEROSOLS' )
     &        .OR.( AOT_REF.EQ.0.D+00)) THEN
	 
           WRITE (IDLOG,*,ERR=1303) " "
	    WRITE (IDLOG,*,ERR=1303) 
     &	    "=> First Calculation of aerosol parameters"
           WRITE (IDLOG,*,ERR=1303) 
     &	    "--> SOS_AEROSOLS is called with the parameters :"
	    WRITE (IDLOG,807,ERR=1303) FICANGLES_RES_MIE
	    WRITE (IDLOG,841,ERR=1303) WA
	    WRITE (IDLOG,842,ERR=1303) TA
	    WRITE (IDLOG,843,ERR=1303) WAREF_AOT
	    WRITE (IDLOG,844,ERR=1303) AOT_REF
	    WRITE (IDLOG,845,ERR=1303) ITRONC_AER
	    WRITE (IDLOG,846,ERR=1303) IMOD_AER
	    WRITE (IDLOG,*,ERR=1303) "      ------------------"
	    WRITE (IDLOG,*,ERR=1303) "      mono-modal model :"
	    WRITE (IDLOG,*,ERR=1303) "      ------------------"
 	    WRITE (IDLOG,847,ERR=1303) RN_WAREF, IN_WAREF 
	    WRITE (IDLOG,849,ERR=1303) IGRANU
	    WRITE (IDLOG,850,ERR=1303) VARGRANU1_MMD_AER, 
     &	     			          VARGRANU2_MMD_AER, 
     &	     			          VARGRANU3_MMD_AER
	    WRITE (IDLOG,*,ERR=1303) "      -----------"
	    WRITE (IDLOG,*,ERR=1303) "      WMO model :"
	    WRITE (IDLOG,*,ERR=1303) "      -----------"
           WRITE (IDLOG,851,ERR=1303) IModele_WMO
	    WRITE (IDLOG,852,ERR=1303) C_WMO_DL,C_WMO_WS,C_WMO_OC,
     &	    			          C_WMO_SO
	    WRITE (IDLOG,*,ERR=1303) "      ----------------------"
	    WRITE (IDLOG,*,ERR=1303) "      Shettle & Fenn model :"
	    WRITE (IDLOG,*,ERR=1303) "      ----------------------"
	    WRITE (IDLOG,853,ERR=1303) IModele_SF
	    WRITE (IDLOG,854,ERR=1303) RH
	    WRITE (IDLOG,*,ERR=1303) "      -----------------------"
	    WRITE (IDLOG,*,ERR=1303) "      Bi-modal distribution :"
	    WRITE (IDLOG,*,ERR=1303) "      -----------------------"
	    WRITE (IDLOG,855,ERR=1303) MODE_PARAM_BILND
	    WRITE (IDLOG,856,ERR=1303) USER_CV_COARSE, USER_CV_FINE
	    WRITE (IDLOG,857,ERR=1303) RTAUct_WAREF
	    WRITE (IDLOG,*,ERR=1303) "      Coarse mode :"
	    WRITE (IDLOG,858,ERR=1303) BMD_CM_MRWAREF, BMD_CM_MIWAREF
	    WRITE (IDLOG,859,ERR=1303) BMD_CM_MRWAREF, BMD_CM_MIWAREF
	    WRITE (IDLOG,860,ERR=1303) BMD_CM_RMODAL, BMD_CM_VAR
	    WRITE (IDLOG,*,ERR=1303) "      Fine mode :"
	    WRITE (IDLOG,858,ERR=1303) BMD_FM_MRWAREF, BMD_FM_MIWAREF
	    WRITE (IDLOG,859,ERR=1303) BMD_FM_MRWAREF, BMD_FM_MIWAREF
	    WRITE (IDLOG,860,ERR=1303) BMD_FM_RMODAL, BMD_FM_VAR
	    WRITE (IDLOG,*,ERR=1303) "      ---------------"
 	    WRITE (IDLOG,*,ERR=1303) "      External data :"
	    WRITE (IDLOG,*,ERR=1303) "      ---------------"
 	    WRITE (IDLOG,861,ERR=1303) FICEXTDATA_AER  
	    WRITE (IDLOG,*,ERR=1303) "      --------------------"
 	    WRITE (IDLOG,*,ERR=1303) "      Mixture definition :"
	    WRITE (IDLOG,*,ERR=1303) "      --------------------"
 	    WRITE (IDLOG,8611,ERR=1303) FICMIXTURE_AER  
	    WRITE (IDLOG,*,ERR=1303) "      ------------------"
	    WRITE (IDLOG,*,ERR=1303) "      MIE computations :"
	    WRITE (IDLOG,*,ERR=1303) "      ------------------"
	    WRITE (IDLOG,864,ERR=1303) DIR_MIE 
	    WRITE (IDLOG,865,ERR=1303) FICMIE_LOG 
	    WRITE (IDLOG,*,ERR=1303) "      --------------"
	    WRITE (IDLOG,*,ERR=1303) "      Output files :"
	    WRITE (IDLOG,*,ERR=1303) "      --------------"
 	    WRITE (IDLOG,862,ERR=1303) FICGRANU_LOG 
	    WRITE (IDLOG,863,ERR=1303) FICGRANU
	    
	  ELSE
	 
	    WRITE (IDLOG,*,ERR=1303) " "
           WRITE (IDLOG,*,ERR=1303) 
     &	    "--> SOS_AEROSOLS is not launched"
     	    IF (FICUSER_AER.NE.'NO_USER_AEROSOLS') THEN
               WRITE (IDLOG,*,ERR=1303) "UserFile : ",FICUSER_AER
 	    ENDIF
	  ENDIF ! end of test on user file or AOT=0
    
      ENDIF !End of Trace
    

C     If no user file => the processing is run
 
      IF (( FICUSER_AER.EQ.'NO_USER_AEROSOLS' )
     &     .OR.( AOT_REF.EQ.0.D+00)) THEN

C --- Run of the aerosol routine

         ITRONC_AER_IN = ITRONC_AER !If ITRONC_AER=1 but the code considers
				 	 !the truncation useless, only the parameter
					 !ITRONC_AER_IN will switch from 1 to 0, without
					 !changing the value of ITRONC_AER

         CALL SOS_AEROSOLS(FICANGLES_RES_MIE, 
     &        WA, TA, WAREF_AOT, AOT_REF,	  
     &        ITRONC_AER_IN, IMOD_AER,
     &        RN_WAREF, IN_WAREF, IGRANU, VARGRANU1_MMD_AER,
     &        VARGRANU2_MMD_AER, VARGRANU3_MMD_AER,
     &        IModele_WMO, 
     &        C_WMO_DL, C_WMO_WS, C_WMO_OC, C_WMO_SO,
     &        IModele_SF, RH, 
     &        MODE_PARAM_BILND, 
     &        USER_CV_COARSE, USER_CV_FINE, RTAUct_WAREF,
     &        BMD_CM_MRWAREF, BMD_CM_MIWAREF, 
     &        BMD_CM_MRWAREF, BMD_CM_MIWAREF,
     &        BMD_CM_RMODAL, BMD_CM_VAR,
     &        BMD_FM_MRWAREF, BMD_FM_MIWAREF, 
     &        BMD_FM_MRWAREF, BMD_FM_MIWAREF,
     &        BMD_FM_RMODAL, BMD_FM_VAR,	 
     &        FICEXTDATA_AER, FICMIXTURE_AER,
     &        DIR_TMP, DIR_MIE, FICMIE_LOG,
     &        FICGRANU, FICGRANU_LOG, KMAT1,
     &        PIZAER, COEF_TRONCA, IER)
           			  
C --- Management of the error case
				  
         IF (IER.NE.0) THEN ! error case
	  			
             IF (TRACE) THEN
	         WRITE (IDLOG,*,ERR=1303) 
     &	         "   ==> Error in subroutine SOS_AEROSOLS"
 	      ENDIF
	 
	      GOTO 4100   
	 
         ENDIF 	  

      ELSE

C    With a user file => no processing

         !Aerosol phase matrix file
         LFIC = INDEX(FICUSER_AER,' ') - 1
         FICGRANU = FICUSER_AER(1:LFIC)
      	  
      ENDIF
      
      
C --- Case WAREF != WA_SIMU   (rem: not the case if a user file is used)
C --- Run of the aerosol routine for WA_SIMU

      IF ((AOT_REF.NE.0.D+00).AND.(WA_SIMU.NE.WAREF_AOT)) THEN

C ---    Backup of the coefficient KMAT computed previously
	  KMAT1_WAREF=KMAT1
	  
C ---    Simulation optical thickness arbitrarily set to 0.1
C ---    to allow the aerosol computations
	  TA=0.1	        

C ---    Management of the trace (display of the computation parameters)
      
         IF (TRACE) THEN
	 
           WRITE (IDLOG,*,ERR=1303) " "
	    WRITE (IDLOG,*,ERR=1303) 
     &	     "=> Second Calculation of aerosol parameters"
	    WRITE (IDLOG,*,ERR=1303) "Case WA = WA_SIMU"	         
           WRITE (IDLOG,*,ERR=1303) 
     &	     "--> SOS_AEROSOLS is called with the parameters :"
	    WRITE (IDLOG,807,ERR=1303) FICANGLES_RES_MIE
	    WRITE (IDLOG,841,ERR=1303) WA_SIMU
	    WRITE (IDLOG,842,ERR=1303) TA
	    WRITE (IDLOG,843,ERR=1303) WAREF_AOT
	    WRITE (IDLOG,844,ERR=1303) AOT_REF
	    WRITE (IDLOG,845,ERR=1303) ITRONC_AER
	    WRITE (IDLOG,846,ERR=1303) IMOD_AER
	    WRITE (IDLOG,*,ERR=1303) "      ------------------"
	    WRITE (IDLOG,*,ERR=1303) "      mono-modal model :"
	    WRITE (IDLOG,*,ERR=1303) "      ------------------"
 	    WRITE (IDLOG,847,ERR=1303) RN_WA, IN_WA
	    WRITE (IDLOG,849,ERR=1303) IGRANU
	    WRITE (IDLOG,850,ERR=1303) VARGRANU1_MMD_AER, 
     &	     			          VARGRANU2_MMD_AER, 
     &	     			          VARGRANU3_MMD_AER
	    WRITE (IDLOG,*,ERR=1303) "      -----------"
	    WRITE (IDLOG,*,ERR=1303) "      WMO model :"
	    WRITE (IDLOG,*,ERR=1303) "      -----------"
           WRITE (IDLOG,851,ERR=1303) IModele_WMO
	    WRITE (IDLOG,852,ERR=1303) C_WMO_DL,C_WMO_WS,C_WMO_OC,
     &	    			      C_WMO_SO
	    WRITE (IDLOG,*,ERR=1303) "      ----------------------"
	    WRITE (IDLOG,*,ERR=1303) "      Shettle & Fenn model :"
	    WRITE (IDLOG,*,ERR=1303) "      ----------------------"
	    WRITE (IDLOG,853,ERR=1303) IModele_SF
	    WRITE (IDLOG,854,ERR=1303) RH
	    WRITE (IDLOG,*,ERR=1303) "      -----------------------"
	    WRITE (IDLOG,*,ERR=1303) "      Bi-modal distribution :"
	    WRITE (IDLOG,*,ERR=1303) "      -----------------------"
	    WRITE (IDLOG,855,ERR=1303) MODE_PARAM_BILND
	    WRITE (IDLOG,856,ERR=1303) USER_CV_COARSE, USER_CV_FINE
	    WRITE (IDLOG,857,ERR=1303) RTAUct_WAREF
	    WRITE (IDLOG,*,ERR=1303) "      Coarse mode :"
	    WRITE (IDLOG,858,ERR=1303) BMD_CM_MRWA, BMD_CM_MIWA
	    WRITE (IDLOG,859,ERR=1303) BMD_CM_MRWAREF, BMD_CM_MIWAREF
	    WRITE (IDLOG,860,ERR=1303) BMD_CM_RMODAL, BMD_CM_VAR
	    WRITE (IDLOG,*,ERR=1303) "      Fine mode :"
	    WRITE (IDLOG,858,ERR=1303) BMD_FM_MRWA, BMD_FM_MIWA
	    WRITE (IDLOG,859,ERR=1303) BMD_FM_MRWAREF, BMD_FM_MIWAREF
	    WRITE (IDLOG,860,ERR=1303) BMD_FM_RMODAL, BMD_FM_VAR
	    WRITE (IDLOG,*,ERR=1303) "      ---------------"
 	    WRITE (IDLOG,*,ERR=1303) "      External data :"
	    WRITE (IDLOG,*,ERR=1303) "      ---------------"
 	    WRITE (IDLOG,861,ERR=1303) FICEXTDATA_AER  
	    WRITE (IDLOG,*,ERR=1303) "      --------------------"
 	    WRITE (IDLOG,*,ERR=1303) "      Mixture definition :"
	    WRITE (IDLOG,*,ERR=1303) "      --------------------"
 	    WRITE (IDLOG,8611,ERR=1303) FICMIXTURE_AER  
	    WRITE (IDLOG,*,ERR=1303) "      ------------------"
	    WRITE (IDLOG,*,ERR=1303) "      MIE computations :"
	    WRITE (IDLOG,*,ERR=1303) "      ------------------"
	    WRITE (IDLOG,864,ERR=1303) DIR_MIE 
	    WRITE (IDLOG,865,ERR=1303) FICMIE_LOG 
	    WRITE (IDLOG,*,ERR=1303) "      --------------"
	    WRITE (IDLOG,*,ERR=1303) "      Output files :"
	    WRITE (IDLOG,*,ERR=1303) "      --------------"
 	    WRITE (IDLOG,862,ERR=1303) FICGRANU_LOG 
	    WRITE (IDLOG,863,ERR=1303) FICGRANU
	    
         ENDIF !End of Trace
            
C ---     Run of the aerosol routine

         ITRONC_AER_IN = ITRONC_AER   !If ITRONC_AER=1 but the code considers
					   !the truncation useless, only the parameter
					   !ITRONC_AER_IN will switch from 1 to 0, without
					   !changing the value of ITRONC_AER

         CALL SOS_AEROSOLS(FICANGLES_RES_MIE, 
     &        WA_SIMU, TA, WAREF_AOT, AOT_REF,	  
     &        ITRONC_AER_IN, IMOD_AER,
     &        RN_WA, IN_WA, IGRANU, VARGRANU1_MMD_AER,
     &        VARGRANU2_MMD_AER, VARGRANU3_MMD_AER,
     &        IModele_WMO, 
     &        C_WMO_DL, C_WMO_WS, C_WMO_OC, C_WMO_SO,
     &        IModele_SF, RH, 
     &        MODE_PARAM_BILND, 
     &        USER_CV_COARSE, USER_CV_FINE, RTAUct_WAREF,
     &        BMD_CM_MRWA, BMD_CM_MIWA, 
     &        BMD_CM_MRWAREF, BMD_CM_MIWAREF,
     &        BMD_CM_RMODAL, BMD_CM_VAR,
     &        BMD_FM_MRWA, BMD_FM_MIWA, 
     &        BMD_FM_MRWAREF, BMD_FM_MIWAREF,
     &        BMD_FM_RMODAL, BMD_FM_VAR,	 
     &        FICEXTDATA_AER, FICMIXTURE_AER,
     &        DIR_TMP, DIR_MIE, FICMIE_LOG,
     &        FICGRANU, FICGRANU_LOG, KMAT1,
     &        PIZAER, COEF_TRONCA, IER)

C ---    Management of the error case
				  
         IF (IER.NE.0) THEN ! error case
	  			
             IF (TRACE) THEN
	          WRITE (IDLOG,*,ERR=1303) 
     &	          "   ==> Error in subroutine SOS_AEROSOLS"
 	      ENDIF
	 
	      GOTO 4100   
	 
         ENDIF 	 

C ---    Computation of the optical thickness at WA_SIMU
	  TA = (KMAT1 / KMAT1_WAREF) * AOT_REF
         IF (TA.LT.0.0) GOTO 2516

C ---    Management of the trace (parameters of the computation of the AOT at WA_SIMU)

	  IF (TRACE) THEN
	      WRITE (IDLOG,*,ERR=1303) " "
	      WRITE (IDLOG,*,ERR=1303) 
     &	      "-----------------------------------------"
	      WRITE (IDLOG,*,ERR=1303)
     &	      "AEROSOLS RADIATIVE PROPERTIES COMPUTATION"
	      WRITE (IDLOG,*,ERR=1303) 
     &	      "-----------------------------------------"	
	      WRITE (IDLOG,*,ERR=1303) " "
	      WRITE (IDLOG,*,ERR=1303) "=> AOT estimation at WA : "
	      WRITE (IDLOG,*,ERR=1303) "    WAREF   : ", WAREF_AOT
	      WRITE (IDLOG,*,ERR=1303) "    AOT_REF : ", AOT_REF
	      WRITE (IDLOG,*,ERR=1303) "    WA      : ", WA_SIMU         
	      WRITE (IDLOG,*,ERR=1303) 
     &	      "    Sect eff ext (WAREF) : ", KMAT1_WAREF
	      WRITE (IDLOG,*,ERR=1303) 
     &	      "    Sect eff ext (WA)    : ", KMAT1
	      WRITE (IDLOG,*,ERR=1303) 
     &	      "    ===> AOT estimated at WA : ", TA
 	  ENDIF !End of Trace

                 			  	    
      ENDIF ! end of case WAREF != WA_SIMU
      	      				  
C =========================================================
C --- Preparation of the surface properties
C =========================================================     
      				

C --- Determination of the output directory
      IF (ISURF.EQ.1) THEN
      	   ! Glitter case
          SURF_OUTDIR(1:LDIRGLIT)=DIR_GLIT(1:LDIRGLIT)
      ENDIF 
      IF (ISURF.EQ.3) THEN
      	   ! Roujean case
          SURF_OUTDIR(1:LDIRROUJ)=DIR_ROUJ(1:LDIRROUJ)
      ENDIF       
      IF (ISURF.EQ.4) THEN
      	   ! Rondeaux Herman case
          SURF_OUTDIR(1:LDIRRH)=DIR_RH(1:LDIRRH)
      ENDIF       
      IF (ISURF.EQ.5) THEN
      	   ! Breon case
          SURF_OUTDIR(1:LDIRBREON)=DIR_BREON(1:LDIRBREON)
      ENDIF       
      IF (ISURF.EQ.6) THEN
      	   ! Nadal case
          SURF_OUTDIR(1:LDIRNADAL)=DIR_NADAL(1:LDIRNADAL)
      ENDIF       
      IF (ISURF.EQ.7) THEN
      	   ! Maignan case
          SURF_OUTDIR(1:LDIRMAIGNAN)=DIR_MAIGNAN(1:LDIRMAIGNAN)
      ENDIF             

C --- Management of the trace
      
      IF (TRACE) THEN
         WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) "SURFACE CALCULATION"
	  WRITE (IDLOG,*,ERR=1303) "------------------"
	  
	  !Lambertian case
	  IF (ISURF.EQ.0) THEN	  
	  
	     WRITE (IDLOG,*,ERR=1303) "Lambertian case "
	     WRITE (IDLOG,*,ERR=1303) 
     &	     "--> SOS_SURFACE is not launched"     
     
         ELSE
	  
	     !Fresnel case
	     IF (ISURF.EQ.2) THEN
	     
	        WRITE (IDLOG,*,ERR=1303) "Lambertian case + Fresnel"
	        WRITE (IDLOG,*,ERR=1303) 
     &	        "--> SOS_SURFACE is not launched"
     
     	     ELSE  

		! User file case
                IF (FICSURF.NE."DEFAULT") THEN
		
                  WRITE (IDLOG,*,ERR=1303) "User Surface File Case"
                  WRITE (IDLOG,*,ERR=1303) "FICSURF ",FICSURF
	           WRITE (IDLOG,*,ERR=1303) 
     &                  "--> SOS_SURFACE is not launched"
                   
                ELSE            ! Other cases, run of the routine
                   
          	    WRITE (IDLOG,*,ERR=1303) 
     &	  	    "--> SOS_SURFACE is called with the parameters :"
	  	    WRITE (IDLOG,*,ERR=1303) "ISURF ",ISURF
	  	    WRITE (IDLOG,*,ERR=1303) "SURF_IND ",SURF_IND
	  	    WRITE (IDLOG,*,ERR=1303) "WIND ",WIND
	  	    WRITE (IDLOG,*,ERR=1303) "K0 ",K0_ROUJEAN
	  	    WRITE (IDLOG,*,ERR=1303) "K1 ",K1_ROUJEAN 
	  	    WRITE (IDLOG,*,ERR=1303) "K2 ",K2_ROUJEAN 
	  	    WRITE (IDLOG,*,ERR=1303) "ALPHA_NADAL ",ALPHA_NADAL
	  	    WRITE (IDLOG,*,ERR=1303) "BETA_NADAL ",BETA_NADAL
	  	    WRITE (IDLOG,*,ERR=1303) "COEF_C_MAIGNAN ",
     &	  	    COEF_C_MAIGNAN
	  	    WRITE (IDLOG,*,ERR=1303) "FICANGLES_RES_LUM ",
     &	  	    FICANGLES_RES_LUM
	  	    WRITE (IDLOG,*,ERR=1303) "FICSURF_LOG ",FICSURF_LOG
	  	    WRITE (IDLOG,*,ERR=1303) "SURF_OUTDIR ",SURF_OUTDIR
	  	    WRITE (IDLOG,*,ERR=1303) "DIR_ROUJ ",DIR_ROUJ
	     
	     	  ENDIF          ! End of user case
                 
            ENDIF             ! End of Fresnel test
	     
         ENDIF                !  End of Lambertian test
	  
      ENDIF                   ! End of trace
      
C --- Run of the tool simulating the surface properties
C --- This tool is not run if ISURF=0 or 2
C --- Or if a surface file FICSURF is defined by the user

      IF ((ISURF.NE.0).AND.(ISURF.NE.2)
     &	    .AND.(FICSURF.EQ."DEFAULT")) THEN
 
         WRITE(6,*) 'SURFACE BRDF/BPDF MATRICES COMPUTATION'   

         IF (TRACE) THEN
            WRITE (IDLOG,*,ERR=1303) " "
	     WRITE (IDLOG,*,ERR=1303) " "
	     WRITE (IDLOG,*,ERR=1303) 
     &	     "-----------------------------------------------"
	     WRITE (IDLOG,*,ERR=1303)
     &	     "SURFACE BRDF/BPDF MATRICES COMPUTATION"
	     WRITE (IDLOG,*,ERR=1303) 
     &	     "-----------------------------------------------"
            WRITE (IDLOG,*,ERR=1303) " "
         ENDIF
        
         CALL SOS_SURFACE(ISURF, FICANGLES_RES_LUM, SURF_IND, WIND,  
     &        K0_ROUJEAN, K1_ROUJEAN , K2_ROUJEAN,  
     &        ALPHA_NADAL, BETA_NADAL, COEF_C_MAIGNAN,  
     &        FICSURF_LOG,                
     &        SURF_OUTDIR, DIR_ROUJ,  DIR_TMP, IER)
     
      ENDIF

      
      
C --- Management of the error case
      
      IF (IER.NE.0) THEN ! error case
      		
         IF (TRACE) THEN
	     WRITE (IDLOG,*,ERR=1303) 
     &	     "   ==> Error in subroutine SOS_SURFACE"
 	  ENDIF
	  GOTO 4200   
	  
      ENDIF
      
C =============================================================
C --- Preparation of the atmospheric profile + radiance computation
C =============================================================     

      WRITE(6,*) 'PROFILE AND RADIANCE COMPUTATION'   
C --- Management of the trace of the main program
      
      IF (TRACE) THEN
      
C	  For the atmospheric profile
          
         WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) "PROFILE CALCULATION"
	  WRITE (IDLOG,*,ERR=1303) "-------------------"	  
         WRITE (IDLOG,*,ERR=1303) 
     &	  "--> SOS_PROFILE is called with the parameters :"
	  WRITE (IDLOG,*,ERR=1303) "WA_SIMU ",WA_SIMU
	  WRITE (IDLOG,*,ERR=1303) "IPROFIL ",IPROFIL
	  WRITE (IDLOG,*,ERR=1303) "TR ",TR
	  WRITE (IDLOG,*,ERR=1303) "HR ",HR
	  WRITE (IDLOG,*,ERR=1303) "TA ",TA
	  WRITE (IDLOG,*,ERR=1303) "HA ",HA
	  WRITE (IDLOG,*,ERR=1303) "PSURF ",PSURF
	  WRITE (IDLOG,*,ERR=1303) "ZMIN ",ZMIN
	  WRITE (IDLOG,*,ERR=1303) "ZMAX ",ZMAX
	  WRITE (IDLOG,*,ERR=1303) "H2O ",H2O
	  WRITE (IDLOG,*,ERR=1303) "O3 ",O3
	  WRITE (IDLOG,*,ERR=1303) "CO2 ",CO2
	  WRITE (IDLOG,*,ERR=1303) "CH4 ",CH4
	  WRITE (IDLOG,*,ERR=1303) "ABSPROFIL ",ABSPROFIL
	  WRITE (IDLOG,*,ERR=1303) "FICABSPROFIL ",FICABSPROFIL
	  WRITE (IDLOG,*,ERR=1303) "FICPROFIL_LOG ",FICPROFIL_LOG

C	  For the radiance computation
      
         WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) 
     &	  "--------------------"
	  WRITE (IDLOG,*,ERR=1303)
     &	  "RADIANCE COMPUTATION"
	  WRITE (IDLOG,*,ERR=1303) 
     &	  "--------------------"    
     
         IF (IGMAX.NE.CTE_DEFAULT_IGMAX)	THEN
	     WRITE (IDLOG,*,ERR=1303) " "
	     WRITE (IDLOG,*,ERR=1303)  
     &	        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            WRITE (IDLOG,*,ERR=1303)  
     &	        "WARNING : The maximal order of interactions has been"
	     WRITE (IDLOG,*,ERR=1303)  
     &	        "          changed by the user  ==> IGMAX =",IGMAX
	     WRITE (IDLOG,*,ERR=1303)  
     &	        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	     WRITE (IDLOG,*,ERR=1303) " "
         ENDIF 

	  WRITE (IDLOG,*,ERR=1303) " "
	  WRITE (IDLOG,*,ERR=1303) 
     &	  "--> SOS_OS is called with the parameters :"
	  WRITE (IDLOG,*,ERR=1303) "FICANGLES_RES_LUM ",
     &	  FICANGLES_RES_LUM
	  WRITE (IDLOG,*,ERR=1303) "FICGRANU ",FICGRANU
	  WRITE (IDLOG,*,ERR=1303) "ZOUT ",ZOUT
	  WRITE (IDLOG,*,ERR=1303) "IGMAX ",IGMAX
	  WRITE (IDLOG,*,ERR=1303) "IPOLAR ",IPOLAR
	  WRITE (IDLOG,*,ERR=1303) "IDIRSURF ",IDIRSURF
	  WRITE (IDLOG,*,ERR=1303) "ITRPHI ",ITRPHI
	  WRITE (IDLOG,*,ERR=1303) "PHIOS ",PHIOS
	  WRITE (IDLOG,*,ERR=1303) "PAS_PHI ",PAS_PHI	  
	  WRITE (IDLOG,*,ERR=1303) "RON ",RON
	  WRITE (IDLOG,*,ERR=1303) "RHO ",RHO
	  WRITE (IDLOG,*,ERR=1303) "ISURF ",ISURF
	  WRITE (IDLOG,*,ERR=1303) "SURF_IND ",SURF_IND
	  WRITE (IDLOG,*,ERR=1303) "WIND ",WIND
	  WRITE (IDLOG,*,ERR=1303) "K0 ",K0_ROUJEAN
	  WRITE (IDLOG,*,ERR=1303) "K1 ",K1_ROUJEAN 
	  WRITE (IDLOG,*,ERR=1303) "K2 ",K2_ROUJEAN 
	  WRITE (IDLOG,*,ERR=1303) "ALPHA_NADAL ",ALPHA_NADAL
	  WRITE (IDLOG,*,ERR=1303) "BETA_NADAL ",BETA_NADAL
	  WRITE (IDLOG,*,ERR=1303) "SURF_OUTDIR ",SURF_OUTDIR
	  WRITE (IDLOG,*,ERR=1303) "FICSURF ",FICSURF
	  WRITE (IDLOG,*,ERR=1303) "FICSOS_RES_BIN ",FICSOS_RES_BIN
	  WRITE (IDLOG,*,ERR=1303) "FICTRANS ",FICTRANS
	  WRITE (IDLOG,*,ERR=1303) "FICFLUX ",FICFLUX
	  WRITE (IDLOG,*,ERR=1303) "FICSOS_LOG ",FICSOS_LOG	  
	     
      ENDIF   ! End of trace

C ---- Opening of the trace files of the profile and OS subfunctions

      IF (TRACEPROF) THEN
          OPEN (IDLOGPROF,FILE=FICPROFIL_LOG,ERR=3104)
      ENDIF
       
      IF (TRACEOS) THEN
          OPEN (IDLOGOS,FILE=FICSOS_LOG,ERR=3106)
      ENDIF
      
C --- Computation of the molecular optical thickness

C     If the optical thickness is not defined, it is computed
C     from the atmospheric pressure with the
C     CNES formula (Perbos, 1982)

      IF (TR.EQ.CTE_NOT_DEFINED_VALUE_DBLE) THEN
      
         TR = (PSURF/CTE_HT_STD_PSURF) * 1.d-4 *       
     &     (84.35/WA_SIMU**4 - 1.225/WA_SIMU**5 + 1.4/WA_SIMU**6)
     	  WRITE(6,*)       
     &   '==> Rayleigh optical thickness from CNES formulation ', TR

         IF (TRACE) THEN            
            WRITE (IDLOG,*,ERR=1303) " "
	     WRITE (IDLOG,*,ERR=1303)
     &      "Molecular Optical Thickness calculated by"
	     WRITE (IDLOG,*,ERR=1303)
     &      "the CNES formulation (Perbos, 1982), TR= ",TR
            WRITE (IDLOG,*,ERR=1303) 
     &      "for the wavelength (microns)=",WA_SIMU
            WRITE (IDLOG,*,ERR=1303) 
     &      "and surface pressure (mbar)=",PSURF	      
            WRITE (IDLOG,*,ERR=1303) " "
	  ENDIF 

      ENDIF

       
C --- Run of the routine preparing the profile for the absorption
C --- Run only if the absorption is taken into account IPROFIL=1 and ABSPROFIL!=7
      
      IF ((IPROFIL.EQ.1).AND.(ABSPROFIL.NE.7)) THEN
     	 
         CALL SOS_PREPA_ABSPROFILE(WA_SIMU, NUSTEP, PSURF, H2O, O3, 
     &		           CO2, CH4, ABSPROFIL, FICABSPROFIL, TRACEPROF,
     &		           IDLOGPROF, NU, LAMB1, IABS,
     &		           ABS_USERPROFIL, ALTABS, RO,
     &			    NEXP, KDIS_AI, KDIS_KI, KDIS_KI_H2O,
     &                   TAB_PRES, NB_PRES, TAB_TEMP, NB_TEMP,
     &                   TAB_CONC_H2O,NB_CONC_H2O, IER)
  
 
C ---    Management of the error case
      
         IF (IER.NE.0) THEN ! error case
      		
            IF (TRACE) THEN
	        WRITE (IDLOG,*,ERR=1303) 
     &	        "   ==> Error in subroutine SOS_PREPA_ABSPROFILE"
 	     ENDIF
	     GOTO 4300   
	  
         ENDIF
 
C ---    Computation of the sum of the weighting coefficients AIK of the radiances per exponential
         sumAIK=0.D+00
         DO IK1=1,NEXP(1,LAMB1)
         DO IK2=1,NEXP(2,LAMB1)
         DO IK3=1,NEXP(3,LAMB1)
         DO IK4=1,NEXP(4,LAMB1)
         DO IK5=1,NEXP(5,LAMB1)
         DO IK6=1,NEXP(6,LAMB1)
         DO IK7=1,NEXP(7,LAMB1)
         DO IK8=1,NEXP(8,LAMB1)
           !Computation of the product of the distribution functions (weight given to the radiances)
            AIK=KDIS_AI(IK1,1,LAMB1)*KDIS_AI(IK2,2,LAMB1)*
     &          KDIS_AI(IK3,3,LAMB1)*KDIS_AI(IK4,4,LAMB1)*
     &          KDIS_AI(IK5,5,LAMB1)*KDIS_AI(IK6,6,LAMB1)*
     &	         KDIS_AI(IK7,7,LAMB1)*KDIS_AI(IK8,8,LAMB1)

            sumAIK=sumAIK+AIK
         ENDDO  !IK1
         ENDDO  !IK2
         ENDDO  !IK3
         ENDDO  !IK4
         ENDDO  !IK5
         ENDDO  !IK6
         ENDDO  !IK7
         ENDDO  !IK8

         IF (TRACE) THEN            
            WRITE (IDLOG,*,ERR=1303) " "
	     WRITE (IDLOG,*,ERR=1303)
     &      "Summation of AIK coefficients:", sumAIK
            WRITE (IDLOG,*,ERR=1303) " "
	  ENDIF 

         !Check of the amplitude of the difference with 1.
         IF (ABS(sumAIK-1.D+00).GE.1.D-06) GOTO 3500

     
      ENDIF ! end of test call of SOS_PREPA_ABSPROFILE


C --- Run of the preparation routine for the OS computation
C --- This routine reads the aerosol and angle files
       
      CALL SOS_PREPA_OS(FICANGLES_RES_LUM, FICGRANU, 
     &	     ISURF, SURF_IND,WIND, 
     &	     K0_ROUJEAN, K1_ROUJEAN , K2_ROUJEAN ,  
     &	     ALPHA_NADAL, BETA_NADAL, COEF_C_MAIGNAN,
     &	     SURF_OUTDIR, FICSURF, TRACEOS, IDLOGOS,
     &	     LUM_NBMU, N0, PIZ, PIZTR, COEFTRONCA, RMU, GA, OS_NB, 
     &	     IMAT_SURF, IGLI, IFRESNEL, IROUJEAN, IRONDEAUX, IBREON, 
     &	     INADAL, IMAIGNAN,
     &	     AER_ALPHA, AER_BETA, AER_GAMMA, AER_ZETA, 
     &	     IND_ANGOUT_FIN, IER)   
	
C --- Management of the error case
      
      IF (IER.NE.0) THEN ! error case
      		
         IF (TRACE) THEN
	      WRITE (IDLOG,*,ERR=1303) 
     &	      "   ==> Error in subroutine SOS_PREPA_OS"
         ENDIF
	  GOTO 4400   
	  
      ENDIF      
   
    
C* ----------------------------------------------------------------------------------------- 
C* Computation with the basic CKD approach: 1 computation per term of the decomposition in exponentials
C* -----------------------------------------------------------------------------------------    
C --- The profile and the radiances are computed for all the CKD coefficients of each species.
C --- There are 8 absorbers => the radiance is computed for each product of exponentials over the 8 absorbers
C --- When the absorption is not taken into account, the number of exponentials for each component is 1

      DO J=1,CTE_OS_NBMU_MAX
         AT_DIR(J)=0.D+00
      ENDDO
      DO J=0,CTE_OS_NBMU_MAX
         TDIR_TRONC_MU(J)=0.D+00
         TDIR_VRAI_MU(J)=0.D+00
      ENDDO
      ATDIR_CKD=.FALSE.
      XMUS=DCOS(INCTE_PI*TETAS/180.D+00)

      IF (IMODE_CKD_CALCUL.EQ.1) THEN  ! CASE RADIANCE COMPUTATION FOR EACH CKD OPT. THICKNESS

         ATDIR_CKD=.TRUE.

         WRITE(6,*)  
     &	  "==> Radiance computation for each CKD term" //
     &   " of the exponential expansion"
         DO IK1=1,NEXP(1,LAMB1)
         DO IK2=1,NEXP(2,LAMB1)
         DO IK3=1,NEXP(3,LAMB1)
         DO IK4=1,NEXP(4,LAMB1)
         DO IK5=1,NEXP(5,LAMB1)
         DO IK6=1,NEXP(6,LAMB1)
         DO IK7=1,NEXP(7,LAMB1)
         DO IK8=1,NEXP(8,LAMB1)


C ---       Display of the information on screen
            write(6,*) '     IK1 (H2O) ',IK1,'/',NEXP(1,LAMB1)
            write(6,*) '     IK2 (CO2) ',IK2,'/',NEXP(2,LAMB1)
            write(6,*) '     IK3 (O3)  ',IK3,'/',NEXP(3,LAMB1)
            write(6,*) '     IK4 (N2O) ',IK4,'/',NEXP(4,LAMB1)
            write(6,*) '     IK5 (CO)  ',IK5,'/',NEXP(5,LAMB1)
            write(6,*) '     IK6 (CH4) ',IK6,'/',NEXP(6,LAMB1)
            write(6,*) '     IK7 (O2)  ',IK7,'/',NEXP(7,LAMB1)
            write(6,*) '     IK8 (NO2) ',IK8,'/',NEXP(8,LAMB1)
            write(6,*) '     --------------------------------'
	   
C ---       Computation of the product of the distribution functions (weight given to the radiances)
            AIK=KDIS_AI(IK1,1,LAMB1)*KDIS_AI(IK2,2,LAMB1)*
     &          KDIS_AI(IK3,3,LAMB1)*KDIS_AI(IK4,4,LAMB1)*
     &          KDIS_AI(IK5,5,LAMB1)*KDIS_AI(IK6,6,LAMB1)*
     &	         KDIS_AI(IK7,7,LAMB1)*KDIS_AI(IK8,8,LAMB1)

C ---       Renormalisation (management of the cases where the sum differs slightly from 1.)
	     AIK=AIK/sumAIK

            !write(6,*) 'AIK ',AIK

	   
C ---       Computation of the profile of absorption optical thicknesses versus altitude

            CALL SOS_ABSPROFILE(ABSPROFIL, NU, LAMB1, IABS, 
     &		               ABS_USERPROFIL, ALTABS, RO, NEXP,
     &		               KDIS_KI, KDIS_KI_H2O, 
     &		               IK1,IK2,IK3,IK4,IK5,IK6,IK7,IK8,
     &                       TAB_PRES, NB_PRES, TAB_TEMP, NB_TEMP,
     &                       TAB_CONC_H2O, NB_CONC_H2O,
     &			        TAUABS, TRACEPROF, IDLOGPROF,IER) 
                  


C ---       Management of the error case
            IF (IER.NE.0) THEN ! error case
      		
               IF (TRACE) THEN
                   WRITE (IDLOG,*,ERR=1303) 
     &	            "   ==> Error in subroutine SOS_ABSPROFILE"
 	        ENDIF
	        GOTO 4500 
	  
            ENDIF

       				
            ! Run of the tool computing the atmospheric profile in optical thicknesses

            CALL SOS_PROFILE(IPROFIL, TR, HR, TA, HA, ZMIN,  
     &		               ZMAX, ABSPROFIL, ALTABS, TAUABS, 
     &			        TRACEPROF, IDLOGPROF, FICPROFIL_TMP,
     &			        NT, IER)


            ! Management of the error case
            IF (IER.NE.0) THEN ! error case
      		
                IF (TRACE) THEN
	             WRITE (IDLOG,*,ERR=1303) 
     &	             "   ==> Error in subroutine SOS_PROFILE"
                ENDIF
	         GOTO 4600   
	  
            ENDIF

 
            !Run of the radiance computation tool

            CALL SOS(FICOS_TMP, FICTRANS, FICPROFIL_TMP,
     &	             NT, ZOUT, IGMAX, IPOLAR, RON,         
     &	             SURF_IND, RHO, IMAT_SURF, IFRESNEL, FICSURF,
     &		     N0, PIZ, PIZTR, COEFTRONCA, 
     &		     RMU, GA, TETAS, OS_NB, LUM_NBMU,
     &		     AER_ALPHA, AER_BETA, AER_GAMMA, AER_ZETA,
     &		     TTOT_TRONC_TMP, TTOT_VRAI_TMP, TAUOUT_TMP, 
     &		     TDIFMUS_TMP, TDIFMUG_TMP, EMOINS_TMP, EPLUS_TMP,
     &               TRACEOS, IDLOGOS, IER)      
      
            ! Management of the error case
            IF (IER.NE.0) THEN ! error case
      		
               IF (TRACE) THEN
	            WRITE (IDLOG,*,ERR=1303) 
     &	            "   ==> Error in subroutine SOS"
               ENDIF
	        GOTO 4700   
	  
            ENDIF

            ! CKD-weighted direct transmission TOA -> surface -> output level
            IF (N0.GT.0) THEN
               DO J=1,LUM_NBMU
                  AT_DIR(J) = AT_DIR(J) + AIK *
     &                 DEXP(-TTOT_TRONC_TMP/RMU(N0)) *
     &                 DEXP(-(TTOT_TRONC_TMP-TAUOUT_TMP)/RMU(J))
               ENDDO
            ENDIF
            ! CKD-weighted direct transmission TOA -> surface
            TDIR_TRONC_MU(0) = TDIR_TRONC_MU(0) + 
     &                         AIK*DEXP(-TTOT_TRONC_TMP/XMUS)
            TDIR_VRAI_MU(0) = TDIR_VRAI_MU(0) + 
     &                         AIK*DEXP(-TTOT_VRAI_TMP/XMUS)
            DO J=1,LUM_NBMU
               TDIR_TRONC_MU(J) = TDIR_TRONC_MU(J) + 
     &                            AIK*DEXP(-TTOT_TRONC_TMP/RMU(J))
               TDIR_VRAI_MU(J) = TDIR_VRAI_MU(J) + 
     &                            AIK*DEXP(-TTOT_VRAI_TMP/RMU(J))
            ENDDO

            ! Reading of the results for the octet ik1*ik2*ik3*ik4*ik5*ik6*ik7*ik8
            ! Weighting and aggregation of the elements of the Stokes vector
            ! Weighting and aggregation of the transmissions
            ! Weighting and aggregation of the optical thicknesses
	   
            CALL SOS_AGGREGATE(LUM_NBMU, AIK, FICOS_TMP, 
     &           TTOT_TRONC_TMP, TTOT_VRAI_TMP, 
     &           TAUOUT_TMP, TDIFMUS_TMP, TDIFMUG_TMP,
     &           EMOINS_TMP, EPLUS_TMP,
     &           FICOS_AGG_TMP,
     &           FICSOS_RES_BIN,
     &           TTOT_TRONC, TTOT_VRAI,           
     &           TAUOUT, TDIFMUS, TDIFMUG, 
     &           EMOINS, EPLUS, IER)


            ! Management of the error case
            IF (IER.NE.0) THEN ! error case
      		
               IF (TRACE) THEN
	            WRITE (IDLOG,*,ERR=1303) 
     &	            "   ==> Error in subroutine SOS_AGGREGATE"
               ENDIF
	        GOTO 4800   
	  
            ENDIF

C ---    End of the loops on the CKD coeff of each absorbing element
         ENDDO  !IK1
         ENDDO  !IK2
         ENDDO  !IK3
         ENDDO  !IK4
         ENDDO  !IK5
         ENDDO  !IK6
         ENDDO  !IK7
         ENDDO  !IK8

         ! Deletion of the temporary file
	  OPEN(UNIT=3,FILE=FICOS_TMP,FORM='UNFORMATTED') 
         CLOSE(3,STATUS='DELETE')    

      ENDIF ! End of CASE RADIANCE COMPUTATION FOR EACH CKD OPT. THICKNESS


C* ----------------------------------------------------------------------------------------- 
C* Computation with the CKD approach applied to the gaseous transmission computation only:
C* ==> Simplified, approximate approach, with 1 single absorption optical thickness computation
C* -----------------------------------------------------------------------------------------    
C --- The profile and the radiances are computed for the absorption optical thickness estimated globally.

      IF (IMODE_CKD_CALCUL.EQ.2) THEN  ! CASE RADIANCE COMPUTATION FOR AN OPT. THICKNESS ESTIMATED BY CKD

         WRITE(6,*)  
     &	  "==> Radiance computation for the absorption optical " //
     &   "thickness calculated by the CKD expansion approach"
         
         DO I=1,CTE_ABS_NBLEV
            TRSCKD(I)=0.D+00
         ENDDO

         DO IK1=1,NEXP(1,LAMB1)
         DO IK2=1,NEXP(2,LAMB1)
         DO IK3=1,NEXP(3,LAMB1)
         DO IK4=1,NEXP(4,LAMB1)
         DO IK5=1,NEXP(5,LAMB1)
         DO IK6=1,NEXP(6,LAMB1)
         DO IK7=1,NEXP(7,LAMB1)
         DO IK8=1,NEXP(8,LAMB1)

	   
C ---       Computation of the product of the distribution functions (weight given to the radiances)
            AIK=KDIS_AI(IK1,1,LAMB1)*KDIS_AI(IK2,2,LAMB1)*
     &          KDIS_AI(IK3,3,LAMB1)*KDIS_AI(IK4,4,LAMB1)*
     &          KDIS_AI(IK5,5,LAMB1)*KDIS_AI(IK6,6,LAMB1)*
     &	         KDIS_AI(IK7,7,LAMB1)*KDIS_AI(IK8,8,LAMB1)

C ---       Renormalisation (management of the cases where the sum differs slightly from 1.)
	     AIK=AIK/sumAIK
	   
C ---       Computation of the profile of absorption optical thicknesses versus altitude
            CALL SOS_ABSPROFILE(ABSPROFIL, NU, LAMB1, IABS, 
     &		               ABS_USERPROFIL, ALTABS, RO, NEXP,
     &		               KDIS_KI, KDIS_KI_H2O,
     &		               IK1,IK2,IK3,IK4,IK5,IK6,IK7,IK8,
     &                       TAB_PRES, NB_PRES, TAB_TEMP, NB_TEMP,
     &                       TAB_CONC_H2O, NB_CONC_H2O,
     &			        TAUABS, TRACEPROF, IDLOGPROF,IER) 

C ---       Management of the error case
            IF (IER.NE.0) THEN ! error case
      		
               IF (TRACE) THEN
                   WRITE (IDLOG,*,ERR=1303) 
     &	            "   ==> Error in subroutine SOS_ABSPROFILE"
 	        ENDIF
	        GOTO 4500
	  
            ENDIF

            !Computation of the transmission per level I, for a fictitious air mass of 1 with the CKD method
            DO I=1,CTE_ABS_NBLEV
               TRSCKD(I)=TRSCKD(I)+AIK*DEXP(-TAUABS(I))
            ENDDO              
	  
C ---    End of the loops on the CKD coeff of each absorbing element
         ENDDO  !IK1
         ENDDO  !IK2
         ENDDO  !IK3
         ENDDO  !IK4
         ENDDO  !IK5
         ENDDO  !IK6
         ENDDO  !IK7
         ENDDO  !IK8

         ! Estimation of the absorption optical thickness profile
         DO I=1,CTE_ABS_NBLEV
             TAUABS(I)=-DLOG(TRSCKD(I))
             if (TAUABS(I).LT.0.D+00) TAUABS(I)=0.D+00	!TEST MODIF: to manage rounding errors giving a transmission very slightly above 1
         ENDDO              

          ! Run of the tool computing the atmospheric profile in optical thicknesses

         CALL SOS_PROFILE(IPROFIL, TR, HR, TA, HA, ZMIN,  
     &		            ZMAX, ABSPROFIL, ALTABS, TAUABS, 
     &			     TRACEPROF, IDLOGPROF, FICPROFIL_TMP,
     &			     NT, IER)

         ! Management of the error case
         IF (IER.NE.0) THEN ! error case
      		
              IF (TRACE) THEN
	           WRITE (IDLOG,*,ERR=1303) 
     &	           "   ==> Error in subroutine SOS_PROFILE"
              ENDIF
	       GOTO 4600   
	  
         ENDIF

 
         !Run of the radiance computation tool

         CALL SOS(FICSOS_RES_BIN, FICTRANS, FICPROFIL_TMP,
     &         NT, ZOUT, IGMAX, IPOLAR, RON,         
     &         SURF_IND, RHO, IMAT_SURF, IFRESNEL, FICSURF,
     &         N0, PIZ, PIZTR, COEFTRONCA, 
     &         RMU, GA, TETAS, OS_NB, LUM_NBMU,
     &         AER_ALPHA, AER_BETA, AER_GAMMA, AER_ZETA,
     &         TTOT_TRONC, TTOT_VRAI, TAUOUT, 
     &         TDIFMUS, TDIFMUG, EMOINS, EPLUS,
     &         TRACEOS, IDLOGOS, IER)      
      
         ! Management of the error case
         IF (IER.NE.0) THEN ! error case
      		
             IF (TRACE) THEN
	          WRITE (IDLOG,*,ERR=1303) 
     &	          "   ==> Error in subroutine SOS"
             ENDIF
	      GOTO 4700   
	  
         ENDIF

         ! Deletion of the temporary profile file
         OPEN(UNIT=3,FILE=FICPROFIL_TMP,FORM='UNFORMATTED') 
         CLOSE(3,STATUS='DELETE') 

      ENDIF ! End of CASE RADIANCE COMPUTATION FOR A SINGLE CKD OPT. THICKNESS



C* ----------------------------------------------------------------------------------------- 
C* Recomposition of the Fourier series for a given viewing azimuth
C* ==> Generation of the result files
C* -----------------------------------------------------------------------------------------    
C --- Creation of the output files (I,Q,U)
      WRITE(6,*) "FINALYSING COMPUTATIONS AND WRITING OUTPUT FILES"


      IF (TRACEOS) THEN
       WRITE(IDLOGOS,*) 
       WRITE(IDLOGOS,*)"***********************************************"
       WRITE(IDLOGOS,*)"Recombination of the radiation field by azimuth"
       WRITE(IDLOGOS,*)"***********************************************"
       WRITE(IDLOGOS,*) 
       IF (ITRPHI.EQ.1) THEN
	  WRITE(IDLOGOS,*) 
     &       "Recombination for the azimuth (degrees) : ",PHIOS
       ELSE
	  WRITE(IDLOGOS,*)
     &      "Recombination for azimuths from 0 to 360 " //
     &      "with a step (degrees): ",PAS_PHI
       ENDIF   
       WRITE(IDLOGOS,*)
      ENDIF    !TRACE
    

C --- Direct transmissions of the solar beam (Beer law) if not
C --- accumulated over the CKD terms
      IF (.NOT.ATDIR_CKD) THEN
         IF (N0.GT.0) THEN
            DO J=1,LUM_NBMU
               AT_DIR(J) = DEXP(-TTOT_TRONC/RMU(N0))*
     &                     DEXP(-(TTOT_TRONC-TAUOUT)/RMU(J))
            ENDDO
         ENDIF
         TDIR_TRONC_MU(0) = DEXP(-TTOT_TRONC/XMUS)
         TDIR_VRAI_MU(0) = DEXP(-TTOT_VRAI/XMUS)
         DO J=1,LUM_NBMU
            TDIR_TRONC_MU(J) = DEXP(-TTOT_TRONC/RMU(J))
            TDIR_VRAI_MU(J) = DEXP(-TTOT_VRAI/RMU(J))
         ENDDO
      ENDIF

C     Direct reflection of the solar beam by the surface (analytic terms
C     added to the Fourier series recombination): only for IDIRSURF = 1
      IGLI_D = IGLI * IDIRSURF
      IFRESNEL_D = IFRESNEL * IDIRSURF
      IROUJEAN_D = IROUJEAN * IDIRSURF
      IRONDEAUX_D = IRONDEAUX * IDIRSURF
      IBREON_D = IBREON * IDIRSURF
      INADAL_D = INADAL * IDIRSURF
      IMAIGNAN_D = IMAIGNAN * IDIRSURF

      CALL SOS_TRPHI_OPTION(LUM_NBMU, RMU, GA,
     &     FICSOS_RES_BIN, TTOT_TRONC, TAUOUT, AT_DIR, ZOUT,
     &     IGLI_D, N0,WIND, SURF_IND, IFRESNEL_D,
     &     IROUJEAN_D, K0_ROUJEAN, K1_ROUJEAN, K2_ROUJEAN,
     &     IRONDEAUX_D, IBREON_D,
     &     INADAL_D, ALPHA_NADAL, BETA_NADAL,
     &     IMAIGNAN_D, COEF_C_MAIGNAN,
     &     ITRPHI ,PHIOS, PAS_PHI, IPOLAR,
     &     PHI_FIN, THETA_FIN, 
     &     SCA_UP_FIN, I_UP_FIN, Q_UP_FIN, U_UP_FIN,
     &     POL_ANG_UP_FIN, POL_RATE_UP_FIN, L_POL_UP_FIN,
     &     SCA_DOWN_FIN, I_DOWN_FIN, Q_DOWN_FIN, U_DOWN_FIN,
     &     POL_ANG_DOWN_FIN, POL_RATE_DOWN_FIN, L_POL_DOWN_FIN,     
     &     IER)
           
C --- Management of the error case
      
      IF (IER.NE.0) THEN ! error case
      		
         IF (TRACEOS) THEN
	         WRITE (IDLOGOS,*,ERR=1303) 
     &	         "   ==> Error in subroutine SOS_TRPHI_OPTION"
         ENDIF
         GOTO 4900   
	  
      ENDIF

C* ----------------------------------------------------------------------------------------- 
C* Saving of the transmissions in a file
C* -----------------------------------------------------------------------------------------    
      IF (FICTRANS.NE.'NO_OUTPUT') THEN
      
         OPEN(UNIT=20,FILE=FICTRANS,FORM='FORMATTED',ERR=997)

C     --> Direct transmission for the true (not truncated) optical thickness
  
	  TDIR_TRONC = TDIR_TRONC_MU(0)
	  TDIR_VRAI = TDIR_VRAI_MU(0)
	  
	  WRITE(20,1005,err=998) TETAS 
	  WRITE(20,1006,err=998) TDIR_VRAI
	  WRITE(20,*,err=998) ' '

C     --> Diffuse transmission TOA --> Ground
	
	  WRITE(20,*,err=998)'Diffuse transmittance : TOA -> surface'
	  WRITE(20,1010,err=998)TETAS,(TDIFMUS + TDIR_TRONC - TDIR_VRAI)
	  WRITE(20,*,err=998) ' '
	  
C     --> Diffuse transmission Ground --> TOA:
C         Loop on the values IMUS: N0 = J	for computation TOA --> Ground (reciprocity)
         WRITE(20,*,err=998)'Diffuse transmittance : surface -> TOA'
	  DO J=1,LUM_NBMU
	  
	     TDIR_TRONC = TDIR_TRONC_MU(J)
	     TDIR_VRAI  = TDIR_VRAI_MU(J)
	  	      
	     WRITE(20,2010,err=998) DACOS(RMU(J))*180.D+00/INCTE_PI, 
     &	                           (TDIFMUG(J) + TDIR_TRONC - TDIR_VRAI)
         ENDDO	  

	  CLOSE(20)
  
      ENDIF !End of test whether the transmission file is to be created

C* ----------------------------------------------------------------------------------------- 
C* Saving of the TOA and BOA fluxes in a file
C* -----------------------------------------------------------------------------------------  

!     The quantities TDIR_VRAI, FLUX_DIFF_DOWN and FLUX_DOWN are retrieved to provide them
!     as outputs of SOS_PROC (even if they are not stored in an optional file).
!     --> Useful for the call of the routine SOS_PROC through the Python binding:

C     Direct transmission for the truncated optical thickness
      TDIR_TRONC = TDIR_TRONC_MU(0)
C     Direct transmission for the true (not truncated) optical thickness
      TDIR_VRAI = TDIR_VRAI_MU(0)
C     Downward diffuse flux
      FLUX_DIFF_DOWN = EMOINS + TDIR_TRONC - TDIR_VRAI
C     Total downward flux
      FLUX_DOWN = EMOINS + TDIR_TRONC

      IF (FICFLUX.NE.'NO_OUTPUT') THEN
      
          OPEN(UNIT=30,FILE=FICFLUX,FORM='FORMATTED',ERR=987)
	  
	   WRITE(30,1005,err=988)TETAS  	                       
	   WRITE(30,*,err=988) ' '

C     --> Downward diffuse flux at BOA: TOA --> Ground
	   WRITE(30,*,err=988)
     &          'Downward fluxes at BOA (normalized by TOA solar flux)'	  
	   WRITE(30,2016,err=988) TDIR_VRAI
	   WRITE(30,2017,err=988) FLUX_DIFF_DOWN
	   WRITE(30,2018,err=988) FLUX_DOWN
	   WRITE(30,*,err=988) ' '
           
C     --> Upward diffuse flux at TOA: Ground --> TOA:
          WRITE(30,*,err=988)
     &    'Upward diffuse flux at TOA (normalized by TOA solar flux):',
     &    EPLUS

          WRITE(30,*,err=988)
          WRITE(30,*,err=988)
          WRITE(30,*,err=988) 'According to the following profile'
          WRITE(30,*,err=988) 'Z(km)    MOT     AOT     GOT     TOTAL'
          DO I = CTE_ABS_NBLEV, 1 , -1
             ZALT = ABS_USERPROFIL(I,1)
             TR_ZALT = TR*DEXP(-ZALT/HR)
             TA_ZALT = TA*DEXP(-ZALT/HA)
             TG_ZALT = TAUABS(CTE_ABS_NBLEV+1-I)
             WRITE(30,2020,err=988) ZALT, TR_ZALT , TA_ZALT, TG_ZALT,
     &                              (TR_ZALT + TA_ZALT + TG_ZALT)
          ENDDO

          CLOSE(30)
	  
      ENDIF !End of test whether the flux file is to be created
     
C* Nominal end
C-------------------
      IF (TRACEPROF) CLOSE(IDLOGPROF)
      IF (TRACEOS)   CLOSE(IDLOGOS)
      IF (TRACE) THEN
         WRITE (IDLOG,*,ERR=1303) " "
         WRITE (IDLOG,*,ERR=1303) 
     &   "--> Nominal end of SOS_PROC processes"
         WRITE (IDLOG,*,ERR=1303) "JOB_STATUS=OK"
         CLOSE(IDLOG)
      ENDIF
      GOTO 9999
      

C =========================================================
C --- Error case and return of status 1 to the shell
C =========================================================  

C --- Definition of the output directory tree
 1000 WRITE(6,*) '  SOS_PROC : ERROR_1000 on parameters'
      WRITE(6,*) '     RESROOT parameter in SOS_PROC'
      WRITE(6,*) '     (-SOS_Main.ResRoot keyword parameter)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==>  Working folder of SOS computations'
      GOTO 9998 

 1001 WRITE(6,*) '  SOS_PROC : ERROR_1001 on parameters'
      WRITE(6,*) '     DIR_MIE parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.DirMie keyword parameter)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==>  Storage directory for MIE files'
      GOTO 9998 

 1002 WRITE(6,*) '  SOS_PROC : ERROR_1002 on parameters'
      WRITE(6,*) '     DIR_SURF parameter in SOS_PROC'
      WRITE(6,*) '     (-SURF.Dir keyword parameter)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==>  Storage directory for SURFACE files'
      GOTO 9998 
      
 1100 WRITE(6,*) '  SOS_PROC : ERROR_1100 : '
      WRITE(6,*) 
     &'  => Working folder name : maybe too long or not defined'
      WRITE(6,*) '     RESROOT parameter in SOS_PROC'
      WRITE(6,*) '     (-SOS_Main.ResRoot keyword parameter)'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998                  	

 1101 WRITE(6,*) '  SOS_PROC : ERROR_1101 : '
      WRITE(6,*) '  => Name of storage directory for MIE files :'
      WRITE(6,*) '     maybe too long or not defined'
      WRITE(6,*) '        DIR_MIE parameter in SOS_PROC'
      WRITE(6,*) '        (-AER.DirMie keyword parameter)'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998                  	

 1102 WRITE(6,*) '  SOS_PROC : ERROR_1102 : '
      WRITE(6,*) '  => Surface storage directory name : '
      WRITE(6,*) '     maybe too long or not defined'
      WRITE(6,*) '        DIR_SURF parameter in SOS_PROC'
      WRITE(6,*) '        (-SURF.Dir keyword parameter)'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998                  	

 1103 WRITE(6,*) '  SOS_PROC : ERROR_1103 : '
      WRITE(6,*) '  => Glitter storage directory name maybe too long'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998                  	

 1104 WRITE(6,*) '  SOS_PROC : ERROR_1104 : '
      WRITE(6,*) '  => Roujean storage directory name maybe too long'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998                  	

 1105 WRITE(6,*) '  SOS_PROC : ERROR_1105 : '
      WRITE(6,*) '  => Rondeaux storage directory name maybe too long'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998                  	

 1106 WRITE(6,*) '  SOS_PROC : ERROR_1106 : '
      WRITE(6,*) '  => Breon storage directory name maybe too long'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998                  	

 1107 WRITE(6,*) '  SOS_PROC : ERROR_1107 : '
      WRITE(6,*) '  => Nadal storage directory name maybe too long'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998                  	

11071 WRITE(6,*) '  SOS_PROC : ERROR_11071 : '
      WRITE(6,*) '  => Maignan storage directory name maybe too long'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998  

 1108 WRITE(6,*) '  SOS_PROC : ERROR_1108 : '
      WRITE(6,*) '  => SOS Working folder name : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998                  	

 1109 WRITE(6,*) '  SOS_PROC : ERROR_1109 : '
      WRITE(6,*) '  => Log folder name : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998                  	

 1110 WRITE(6,*) '  SOS_PROC : ERROR_1110 : '
      WRITE(6,*) '  => TMP folder name : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998   

C --- Creation of the output directory tree

 1201 WRITE(6,*) '  SOS_PROC : ERROR_1201 : '
      WRITE(6,*) 
     &'  => Working folder creation command : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998                  	

 1202 WRITE(6,*) '  SOS_PROC : ERROR_1202 : '
      WRITE(6,*) 
     &'  => Mie Working folder creation command : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998                  	

 1203 WRITE(6,*) '  SOS_PROC : ERROR_1203 : '
      WRITE(6,*) 
     &'  => Surface Working folder creation command : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998                  	

 1204 WRITE(6,*) '  SOS_PROC : ERROR_1204 : '
      WRITE(6,*) 
     &'  => Glitter Working folder creation command : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998                  	

 1205 WRITE(6,*) '  SOS_PROC : ERROR_1205 : '
      WRITE(6,*) 
     &'  => Roujean Working folder creation command : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998                  	

 1206 WRITE(6,*) '  SOS_PROC : ERROR_1206 : '
      WRITE(6,*) 
     &'  => Rondeaux Working folder creation command : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998                  	

 1207 WRITE(6,*) '  SOS_PROC : ERROR_1207 : '
      WRITE(6,*) 
     &'  => Breon Working folder creation command : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998                  	

 1208 WRITE(6,*) '  SOS_PROC : ERROR_1208 : '
      WRITE(6,*) 
     &'  => Nadal Working folder creation command : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998                  	

12081 WRITE(6,*) '  SOS_PROC : ERROR_12081 : '
      WRITE(6,*) 
     &'  => Maignan Working folder creation command : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998  

 1209 WRITE(6,*) '  SOS_PROC : ERROR_1209 : '
      WRITE(6,*)  
     &'  => SOS Directory results creation command : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998                  	

 1210 WRITE(6,*) '  SOS_PROC : ERROR_1210 : '
      WRITE(6,*) 
     &'  => Log folder creation command : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998                  	

 1211 WRITE(6,*) '  SOS_PROC : ERROR_1211 : '
      WRITE(6,*) 
     &'  => TMP folder creation command : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998   

C --- log management
                                           
 1301 WRITE(6,*) '  SOS_PROC : ERROR_1301 : '
      WRITE(6,*) '  => Log filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	

 1302 WRITE(6,*) '  SOS_PROC : ERROR_1302 on Log file opening'
      GOTO 9998

 1303 WRITE(6,*) '  SOS_PROC : ERROR_1303 on Log file writing'
      GOTO 9998
            
C --- Check of the simulation parameters

 2100 WRITE(6,*) '  SOS_PROC : ERROR_2100 on parameters'
      WRITE(6,*) '     WA_SIMU parameter is required in SOS_PROC'
      WRITE(6,*) '     (-SOS_Main.Wa keyword parameter)'
      WRITE(6,*) '  ==>  Wavelength for radiance simulation (microns)'
      GOTO 9998 

 2101 WRITE(6,*) '  SOS_PROC : ERROR_2101 on parameter value'
      WRITE(6,*) '     WA_SIMU in SOS_PROC (-SOS_Main.Wa keyword)'
      WRITE(6,*) '     is lower than the minimal available value'
      WRITE(6,*) '     or higher than the maximal available value'
      WRITE(6,*) '         WA_SIMU = ', WA_SIMU
      WRITE(6,*) '         CTE_WAMIN = ', CTE_WAMIN
      WRITE(6,*) '         CTE_WAMAX = ', CTE_WAMAX

C --- Check of the parameters for the angle computations

 2200 WRITE(6,*) '  SOS_PROC : ERROR_2200 on parameters'
      WRITE(6,*) '     TETAS parameter is required in SOS_PROC'
      WRITE(6,*) '     (-ANG.Thetas keyword parameter)'
      WRITE(6,*) '  ==> Solar zenithal angle (degrees)'
      GOTO 9998 

 2201 WRITE(6,*) '  SOS_PROC : ERROR_2201 on parameter value'
      WRITE(6,*) '     TETAS in SOS_PROC (-ANG.Thetas keyword)'
      WRITE(6,*) '  --> Solar Zenith Angle out of available range'
      GOTO 9998  	

C --- Check of the parameters of the aerosol computations

 2301 WRITE(6,*) '  SOS_PROC : ERROR_2301 on parameters'
      WRITE(6,*) '     AOT_REF parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.AOTref keyword parameter)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) 
     &'  ==> Aerosol optical thickness at reference wavelength '
      GOTO 9998   

 2302 WRITE(6,*) '  SOS_PROC : ERROR_2302 on parameters'
      WRITE(6,*) '     WAREF_AOT parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.Waref keyword parameter)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*)  
     &'  ==> Reference wavelength for the aerosol optical thickness'
      GOTO 9998   

 2304 WRITE(6,*) '  SOS_PROC : ERROR_2304 on parameters'
      WRITE(6,*) '     IMOD_AER parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.Model keyword parameter - Integer)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Type of aerosol model : '
      WRITE(6,*) '         0 : Mono-modal '
      WRITE(6,*) '         1 : WMO multi-modal '   
      WRITE(6,*) '         2 : Shettle & Fenn bi-modal '    
      WRITE(6,*) '         3 : Log-Normal bi-modal '
      WRITE(6,*) '         4 : Phaze functions from an external source'
      WRITE(6,*) '         5 : Multi-modes model defined by user'
      GOTO 9998   

 2305 WRITE(6,*) '  SOS_PROC : ERROR_2305 on parameter value'
      WRITE(6,*) '     IMOD_AER parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.Model keyword parameter - Integer)'
      WRITE(6,*) '  --> out of available range.'
      WRITE(6,*) '  ==> Type of aerosol model : '
      WRITE(6,*) '         0 : Mono-modal '
      WRITE(6,*) '         1 : WMO multi-modal '   
      WRITE(6,*) '         2 : Shettle & Fenn bi-modal '    
      WRITE(6,*) '         3 : Log-Normal bi-modal '
      WRITE(6,*) '         4 : Phaze functions from an external source'
      WRITE(6,*) '         5 : Mixture defined by user file'    
      GOTO 9998 
   
 2309 WRITE(6,*) '  SOS_PROC : ERROR_2309 on parameters'
      WRITE(6,*) '     RN_WA and IN_WA parameters in SOS_PROC'
      WRITE(6,*) '     keyword parameter:'
      WRITE(6,*) '      -AER.MMD.MRwa Value (format : F5.3)'
      WRITE(6,*) '      -AER.MMD.MIwa Value (format : F8.5)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) 
     &'  ==> Aerosol refractive index for mono-modal size distribution'
      WRITE(6,*) '      (real and imaginary parts)'
      WRITE(6,*) '      at the wavelength of radiance calculation.'
      GOTO 9998   

 2310 WRITE(6,*) '  SOS_PROC : ERROR_2310 on parameter value'
      WRITE(6,*) '     IN_WA parameter in SOS_PROC'
      WRITE(6,*) '     (AER.MMD.MIwa keyword parameter)'
      WRITE(6,*) '  ==> Imagery part of Aerosol refractive index '
      WRITE(6,*) '      must be a negative value'
      GOTO 9998   

 2311 WRITE(6,*) '  SOS_PROC : ERROR_2311 on parameters'
      WRITE(6,*) '     IGRANU parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.MMD.SDtype keyword parameter - Integer)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Model index for mono-modal size distribution :'
      WRITE(6,*) '         1 : LND '
      WRITE(6,*) '         2 : Junge law '
      GOTO 9998    

 2312 WRITE(6,*) '  SOS_PROC : ERROR_2312 on parameter value'
      WRITE(6,*) '     IGRANU parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.MMD.SDtype keyword parameter - Integer)'
      WRITE(6,*) '  ==> Model index for mono-modal size distribution :'
      WRITE(6,*) '         1 : LND '
      WRITE(6,*) '         2 : Junge law '
      GOTO 9998      

 2314 WRITE(6,*) '  SOS_PROC : ERROR_2314 on parameters'
      WRITE(6,*) '     RN_WAREF and IN_WAREF parameters in SOS_PROC'
      WRITE(6,*) '     keyword parameter:'
      WRITE(6,*) '      -AER.MMD.MRwaref Value (format : F5.3)'
      WRITE(6,*) '      -AER.MMD.MIwaref Value (format : F8.5)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Aerosol refractive index '
      WRITE(6,*) '      for mono-modal size distribution'
      WRITE(6,*) '      (real and imaginary parts)'
      WRITE(6,*) '      at the wavelength of reference.'
      GOTO 9998  

23131 WRITE(6,*) '  SOS_PROC : ERROR_23131 on parameters'
      WRITE(6,*) 
     &'     LND_RADIUS_MMD_AER and LND_LNVAR_MMD_AER in SOS_PROC'
      WRITE(6,*) '     (-AER.MMD.LNDradius keyword parameter)'
      WRITE(6,*) '     (-AER.MMD.LNDvar keyword parameter)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  Case IGRANU = 1 (AER.MMD.SDtype = 1)'
      WRITE(6,*) '  ==> Parameters of LND mono-modal aerosol size '
      WRITE(6,*) '      distribution :'
      WRITE(6,*) '         - modal radius (microns) '
      WRITE(6,*) '         - and standard deviation'
      GOTO 9998  	     

23132 WRITE(6,*) '  SOS_PROC : ERROR_23132 on parameters'
      WRITE(6,*) 
     &'     JD_SLOPE_MMD_AER and JD_RMIN_MMD_AER in SOS_PROC'
      WRITE(6,*) '     -AER.MMD.JD.slope keyword parameter (Float)'
      WRITE(6,*) '     -AER.MMD.JD.rmin  keyword parameter (Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  Case IGRANU = 2 (AER.MMD.SDtype = 2)'
      WRITE(6,*) '  ==> Parameters of Junge law mono-modal'
      WRITE(6,*) '      aerosol size distribution :'
      WRITE(6,*) '         - slope'
      WRITE(6,*) '         - minimal radius (microns)'
      GOTO 9998 

23141 WRITE(6,*) '  SOS_PROC : ERROR_23141 on parameters'
      WRITE(6,*) '     ITRONC_AER parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.Tronca keyword parameter)'
      WRITE(6,*) '  --> Defined value of truncation option (0 or 1)'
      GOTO 9998     

 2315 WRITE(6,*) '  SOS_PROC : ERROR_2315 on parameters'
      WRITE(6,*) '     IModele_WMO parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.WMO.Model keyword parameter - Integer)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Type of WMO model :'
      WRITE(6,*) '         1 : Continental WMO model. '
      WRITE(6,*) '         2 : Maritime WMO model. '
      WRITE(6,*) '         3 : Urban WMO model. '
      WRITE(6,*) '         4 : WMO model by user definition. '
      GOTO 9998    
         
 2316 WRITE(6,*) '  SOS_PROC : ERROR_2316 on parameter value'
      WRITE(6,*) '     IModele_WMO parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.WMO.Model keyword parameter - Integer)'
      WRITE(6,*) '  --> out of available range'
      WRITE(6,*) '  ==> Type of WMO model :'
      WRITE(6,*) '         1 : Continental WMO model. '
      WRITE(6,*) '         2 : Maritime WMO model. '
      WRITE(6,*) '         3 : Urban WMO model. '
      WRITE(6,*) '         4 : WMO model by user definition. '
      GOTO 9998      		  

 2317 WRITE(6,*) '  SOS_PROC : ERROR_2317 on parameters'
      WRITE(6,*) '     C_WMO_DL, C_WMO_WS,C_WMO_OC '
      WRITE(6,*) '     and C_WMO_SO parameters in SOS_PROC'
      WRITE(6,*) '     (-AER.WMO.DL keyword parameter - float)'
      WRITE(6,*) '     (-AER.WMO.WS keyword parameter - float)'
      WRITE(6,*) '     (-AER.WMO.OC keyword parameter - float)'
      WRITE(6,*) '     (-AER.WMO.SO keyword parameter - float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Volumetric concentration (%) for '
      WRITE(6,*) '         WMO Dust-Like particles'
      WRITE(6,*) '         WMO Water Soluble particles'	    
      WRITE(6,*) '         WMO OCeanic particles'
      WRITE(6,*) '         WMO SOot particles'	     
      GOTO 9998 

 2318 WRITE(6,*) '  SOS_PROC : ERROR_2318 on parameters'
      WRITE(6,*) '     IModele_SF parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.SF.Model keyword parameter - Integer)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Type of Shettle & Fenn model :'		  
      WRITE(6,*) '         1 : Tropospheric S&F model. '
      WRITE(6,*) '         2 : Urban S&F model. '
      WRITE(6,*) '         3 : Maritime S&F model. '
      WRITE(6,*) '         4 : Coastal S&F model. '
      GOTO 9998    		 

 2319 WRITE(6,*) '  SOS_PROC : ERROR_2319 on parameters'
      WRITE(6,*) '     RH parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.SF.RH keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> % of relative humidity'
      GOTO 9998 

 2320 WRITE(6,*) '  SOS_PROC : ERROR_2320 on parameter value'
      WRITE(6,*) '     IModele_SF parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.SF.Model keyword parameter - Integer)'
      WRITE(6,*) '  --> out of available range'
      WRITE(6,*) '  ==> Type of Shettle & Fenn model :'		  
      WRITE(6,*) '         1 : Tropospheric S&F model. '
      WRITE(6,*) '         2 : Urban S&F model. '
      WRITE(6,*) '         3 : Maritime S&F model. '
      WRITE(6,*) '         4 : Coastal S&F model. '
      GOTO 9998    		 

 2321 WRITE(6,*) '  SOS_PROC : ERROR_2321 on parameters'
      WRITE(6,*) '  --> -AER.SF.RH Value    (format : F5.2): '
      WRITE(6,*) '  ==> Relative humidity: between 0 and 99%'
      GOTO 9998   		    
     
 2322 WRITE(6,*) '  SOS_PROC : ERROR_2322 on parameters'
      WRITE(6,*) '     MODE_PARAM_BILND parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.BMD.VCdef keyword parameter - Integer)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Type of selected method for bi-modal LND'
      WRITE(6,*) '      volumetric concentration definition :'
      WRITE(6,*) '       1 : Use of user volumetric concentrations.'
      WRITE(6,*) '       2 : Use of user ratio AOT_coarse / AOT_total.'
      GOTO 9998   		

 2323 WRITE(6,*) '  SOS_PROC : ERROR_2323 on parameter value'
      WRITE(6,*) '     MODE_PARAM_BILND parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.BMD.VCdef keyword parameter - Integer)'
      WRITE(6,*) '  --> out of available range'
      WRITE(6,*) '  ==> Type of selected method for bi-modal LND'
      WRITE(6,*) '      volumetric concentration definition :'
      WRITE(6,*) '       1 : Use of user volumetric concentrations.'
      WRITE(6,*) '       2 : Use of user ratio AOT_coarse / AOT_total.'
      GOTO 9998   		    

 2324 WRITE(6,*) '  SOS_PROC : ERROR_2324 on parameters'
      WRITE(6,*) '     USER_CV_COARSE parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.BMD.CoarseVC keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Volumetric concentration for bimodal LND'      
      GOTO 9998   		  

 2325 WRITE(6,*) '  SOS_PROC : ERROR_2325 on parameters'
      WRITE(6,*) '     USER_CV_FINE parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.BMD.FineVC keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Volumetric concentration for bimodal LND'      
      GOTO 9998   		  

 2326 WRITE(6,*) '  SOS_PROC : ERROR_2326 on parameters'
      WRITE(6,*) '     RTAUct_WAREF parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.BMD.RAOT keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> User volumetric concentration for bimodal LND'  
      WRITE(6,*) '      defined by AOT ratio coarse / total'  
      GOTO 9998  

 2327 WRITE(6,*) '  SOS_PROC : ERROR_2327 on parameters'
      WRITE(6,*) '     BMD_CM_MRWA, BMD_CM_MIWA, BMD_CM_RMODAL '
      WRITE(6,*) '     and BMD_CM_VAR in SOS_PROC'
      WRITE(6,*) '     (-AER.BMD.CM.MRwa keyword parameter -Float)'
      WRITE(6,*) '     (-AER.BMD.CM.MIwa keyword parameter -Float)'
      WRITE(6,*) '     (-AER.BMD.CM.SDradius keyword parameter -Float)'
      WRITE(6,*) '     (-AER.BMD.CM.SDvar keyword parameter -Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Parameters of the bimode LND coarse mode '
      WRITE(6,*) '      for wavelength -SOS_Main.Wa (WA_SIMU)'
      GOTO 9998  

 2328 WRITE(6,*) '  SOS_PROC : ERROR_2328 on parameters'
      WRITE(6,*) '     BMD_FM_MRWA, BMD_FM_MIWA, BMD_FM_RMODAL '
      WRITE(6,*) '     and BMD_FM_VAR in SOS_PROC'
      WRITE(6,*) '     (-AER.BMD.FM.MRwa keyword parameter -Float)'
      WRITE(6,*) '     (-AER.BMD.FM.MIwa keyword parameter -Float)'
      WRITE(6,*) '     (-AER.BMD.FM.SDradius keyword parameter -Float)'
      WRITE(6,*) '     (-AER.BMD.FM.SDvar keyword parameter -Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Parameters of the bimode LND fine mode '
      WRITE(6,*) '      for wavelength -SOS_Main.Wa (WA_SIMU)'
      GOTO 9998  

 2329 WRITE(6,*) '  SOS_PROC : ERROR_2329 on parameters'
      WRITE(6,*) '     BMD_CM_MRWAREF, BMD_CM_MIWAREF,'
      WRITE(6,*) '     BMD_FM_MRWAREF, BMD_FM_MIWAREF in SOS_PROC'
      WRITE(6,*) '     (-AER.BMD.CM.MRwaref keyword parameter - Float)'
      WRITE(6,*) '     (-AER.BMD.CM.MIwaref keyword parameter - Float)'
      WRITE(6,*) '     (-AER.BMD.FM.MRwaref keyword parameter - Float)'
      WRITE(6,*) '     (-AER.BMD.FM.MIwaref keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Refractive index at the reference wavelength' 
      WRITE(6,*) '      for the bimode LND coarse and fine modes'
      GOTO 9998     

 2330 WRITE(6,*) '  SOS_PROC : ERROR_2330 on parameters'
      WRITE(6,*) '     FICEXTDATA_AER parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.ExtData keyword parameter)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> User phase function filename'  
      GOTO 9998  

 2331 WRITE(6,*) '  SOS_PROC : ERROR_2331 on parameter value'
      WRITE(6,*)  
     & '  --> Simulation and reference wavelengths must be equal'
      WRITE(6,*) '    WA_SIMU and WAREF_AOT parameters in SOS_PROC'
      WRITE(6,*) '    (-SOS_Main.Wa and -AER.Waref keyword parameters)'
      WRITE(6,*) '  ==> Case of an external User phase function'  
      GOTO 9998  

 2340 WRITE(6,*) '  SOS_PROC : ERROR_2340 on parameters'
      WRITE(6,*) '     FICMIXTURE_AER parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.DefMixture keyword parameter)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) 
     & '  ==> User definition of the mixture of aerosol models'
      GOTO 9998  

 2350 WRITE(6,*) '  SOS_PROC : ERROR_2350 on parameter value'
      WRITE(6,*)  
     & '  --> Simulation and reference wavelengths must be equal'
      WRITE(6,*) '    WA_SIMU and WAREF_AOT parameters in SOS_PROC'
      WRITE(6,*) '    (-SOS_Main.Wa and -AER.Waref keyword parameters)'
      WRITE(6,*) '  ==> Case of an aerosol User File'  
      GOTO 9998  

 2351 WRITE(6,*) '  SOS_PROC : ERROR_2351 on parameter value'
      WRITE(6,*) '     FIC_GRANU parameter in SOS_PROC'
      WRITE(6,*) '     (-AER.ResFile keyword parameter)'
      WRITE(6,*)  
     &'     is defined for the case of using an aerosol User File'
      WRITE(6,*) '     This is inconsistent.'
      WRITE(6,*)  
     &'  --> Parameter -AER.ResFile (FIC_GRANU) must not be defined !'
      GOTO 9998   


C --- Check of the surface parameters

 2401 WRITE(6,*) '  SOS_PROC : ERROR_2401 on parameters'
      WRITE(6,*) '     RHO parameter in SOS_PROC'
      WRITE(6,*) '     (-SURF.Alb keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==>  Surface albedo'
      GOTO 9998 

 2402 WRITE(6,*) '  SOS_PROC : ERROR_2402 on parameter value'
      WRITE(6,*) '     RHO parameter in SOS_PROC'
      WRITE(6,*) '     (-SURF.Alb keyword parameter - Float)'
      WRITE(6,*) '  --> must be positive '
      WRITE(6,*) '  ==>  Surface albedo'
      GOTO 9998 
 
 2403 WRITE(6,*) '  SOS_PROC : ERROR_2403 on parameters'
      WRITE(6,*) '     ISURF parameter in SOS_PROC'
      WRITE(6,*) '     (-SURF.Type keyword parameter - Integer)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==>  Type of surface model'
      WRITE(6,*) '         0 : lambertian reflectance '
      WRITE(6,*) '         1 : lambertian reflectance + glitter'
      WRITE(6,*) '         2 : lambertian reflectance + fresnel'
      WRITE(6,*) '         3 : lambertian reflectance + Roujean BRDF'
      WRITE(6,*) '         4 : lambertian reflectance + Roujean BRDF',
     &		 ' + BPDF Rondeaux'
      WRITE(6,*) '         5 : lambertian reflectance + Roujean BRDF',
     &		 ' + BPDF Breon'  
      WRITE(6,*) '         6 : lambertian reflectance + Roujean BRDF',
     &		 ' + BPDF Nadal'     
      WRITE(6,*) '         7 : lambertian reflectance + Roujean BRDF',
     &		 ' + BPDF Maignan'             
      GOTO 9998 

 2404 WRITE(6,*) '  SOS_PROC : ERROR_2404 on parameter value'
      WRITE(6,*) '     ISURF parameter in SOS_PROC'
      WRITE(6,*) '     (-SURF.Type keyword parameter - Integer)'
      WRITE(6,*) '  --> out of available range'
      WRITE(6,*) '  ==>  Type of surface model'
      WRITE(6,*) '         0 : lambertian reflectance '
      WRITE(6,*) '         1 : lambertian reflectance + glitter'
      WRITE(6,*) '         2 : lambertian reflectance + fresnel'
      WRITE(6,*) '         3 : lambertian reflectance + Roujean BRDF'
      WRITE(6,*) '         4 : lambertian reflectance + Roujean BRDF',
     &		 ' + BPDF Rondeaux'
      WRITE(6,*) '         5 : lambertian reflectance + Roujean BRDF',
     &		 ' + BPDF Breon'  
      WRITE(6,*) '         6 : lambertian reflectance + Roujean BRDF',
     &		 ' + BPDF Nadal'     
      WRITE(6,*) '         7 : lambertian reflectance + Roujean BRDF',
     &		 ' + BPDF Maignan'             
      GOTO 9998 

 2405 WRITE(6,*) '  SOS_PROC : ERROR_2405 on parameter'
      WRITE(6,*) '     SURF_IND parameter in SOS_PROC'
      WRITE(6,*) '     (-SURF.Ind keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '         Required if surface type=1 2 4 5 6 or 7'
      GOTO 9998 

 2406 WRITE(6,*) '  SOS_PROC : ERROR_2406 on parameter value'
      WRITE(6,*) '     WIND parameter in SOS_PROC'
      WRITE(6,*) '     (-SURF.Glitter.Wind keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==>  Wind speed in m/s (required if type=1)'
      GOTO 9998 

24061 WRITE(6,*) '  SOS_PROC : ERROR_24061 on parameter value'
      WRITE(6,*) '     WIND parameter in SOS_PROC'
      WRITE(6,*) '     (-SURF.Glitter.Wind keyword parameter - Float)'
      WRITE(6,*) '  ==>  Wind speed must be positive (in m/s)'
      GOTO 9998 

 2407 WRITE(6,*) '  SOS_PROC : ERROR_2407 on parameters'
      WRITE(6,*) '     K0_ROUJEAN, K1_ROUJEAN,K2_ROUJEAN in SOS_PROC'
      WRITE(6,*) '     (-SURF.Roujean.K0 keyword parameter - Float)'
      WRITE(6,*) '     (-SURF.Roujean.K1 keyword parameter - Float)'
      WRITE(6,*) '     (-SURF.Roujean.K2 keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'  
      WRITE(6,*) '  ==> Roujean BRDF parameters (required if type>=3)'
      GOTO 9998 

 2408 WRITE(6,*) '  SOS_PROC : ERROR_2408 on parameters'
      WRITE(6,*) '     ALPHA_NADAL and BETA_NADAL in SOS_PROC'
      WRITE(6,*) '     (-SURF.Nadal.Alpha keyword parameter - Float)'
      WRITE(6,*) '     (-SURF.Nadal.Beta keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Nadal BPDF parameters (required if type=6)'
      GOTO 9998 
 
 2409 WRITE(6,*) '  SOS_PROC : ERROR_2409 on parameter value'
      WRITE(6,*) '     ALPHA_NADAL parameter in SOS_PROC'
      WRITE(6,*) '     (-SURF.Nadal.Alpha keyword parameter - Float)'
      WRITE(6,*) 
     &'  --> parameter value  must be included between 0 and 10' 
      GOTO 9998 

 2410 WRITE(6,*) '  SOS_PROC : ERROR_2410 on parameter value'
      WRITE(6,*) '     BETA_NADAL parameter in SOS_PROC'
      WRITE(6,*) '     (-SURF.Nadal.Beta keyword parameter - Float)'
      WRITE(6,*) 
     &'  --> parameter value  must be included between 0 and 100' 
      GOTO 9998 

 2411 WRITE(6,*) '  SOS_PROC : ERROR_2411 on parameters'
      WRITE(6,*) '     COEF_C_MAIGNAN parameter in SOS_PROC'
      WRITE(6,*) '     (-SURF.Maignan.C  keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) 
     &"  ==>  C*exp(-NDVI) coefficient of the Maignan's BPDF model"
      GOTO 9998 

C --- Check of the profile parameters

 2502 WRITE(6,*) '  SOS_PROC : ERROR_2502 on parameter value'
      WRITE(6,*) '     TR parameter in SOS_PROC'
      WRITE(6,*) '     (-AP.MOT  keyword parameter - Float)'
      WRITE(6,*) '  --> - must be positive ' 
      WRITE(6,*) '  ==> Molecular Optical Thickness' 
      GOTO 9998 

 2503 WRITE(6,*) '  SOS_PROC : ERROR_2503 on parameters'
      WRITE(6,*) '     HR parameter in SOS_PROC'
      WRITE(6,*) '     (-AP.HR  keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Molecular height scale' 
      WRITE(6,*) '      Case of non null AP.MOT value (TR in SOS_PROC)'       
      GOTO 9998    

 2504 WRITE(6,*) '  SOS_PROC : ERROR_2504 on parameter value'
      WRITE(6,*) '     HR parameter in SOS_PROC'
      WRITE(6,*) '     (-AP.HR  keyword parameter - Float)'
      WRITE(6,*) '  --> must be positive ' 
      WRITE(6,*) '  ==> Molecular height scale' 
      WRITE(6,*) '      Case of non null AP.MOT value (TR in SOS_PROC)' 
      GOTO 9998    
 
 2505 WRITE(6,*) '  SOS_PROC : ERROR_2505 on parameters'
      WRITE(6,*) '     IPROFIL parameter in SOS_PROC'
      WRITE(6,*) 
     &'     (-AP.AerProfile.Type  keyword parameter - Integer)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Aerosol atmospheric profil type' 
      GOTO 9998    
 
 2506 WRITE(6,*) '  SOS_PROC : ERROR_2506 on parameter value'
      WRITE(6,*) '     IPROFIL parameter in SOS_PROC'
      WRITE(6,*) 
     &'     (-AP.AerProfile.Type  keyword parameter - Integer)'
      WRITE(6,*) '  --> Min value: 1'    
      WRITE(6,*) '  --> Max value: CTE_NB_PROFIL (constant)'          
      WRITE(6,*) '  CTE_NB_PROFIL: ',CTE_NB_PROFIL
      WRITE(6,*) '  -AP.AerProfile.Type: ',IPROFIL   
      WRITE(6,*) '  ==> Aerosol atmospheric profil type' 
      GOTO 9998    

 2507 WRITE(6,*) '  SOS_PROC : ERROR_2507 on parameters'
      WRITE(6,*) '     HA parameter in SOS_PROC'
      WRITE(6,*) '     (-AP.AerHS.HA  keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Height scale of the aerosol profile'
      WRITE(6,*) 
     &'  ==> Case of non null -AER.AOTref value (AOT_REF in SOS_PROC)'
      GOTO 9998    

 2508 WRITE(6,*) '  SOS_PROC : ERROR_2508 on parameter value'
      WRITE(6,*) '     HA parameter in SOS_PROC'
      WRITE(6,*) '     (-AP.AerHS.HA  keyword parameter - Float)'
      WRITE(6,*) '  --> must be positive '
      WRITE(6,*) '  ==> Height scale of the aerosol profile'
      WRITE(6,*) 
     &'  ==> Case of non null -AER.AOTref value (AOT_REF in SOS_PROC)'
      GOTO 9998    
 
 2509 WRITE(6,*) '  SOS_PROC : ERROR_2509 on parameters'
      WRITE(6,*) '     ZMIN and ZMAX parameters in SOS_PROC'
      WRITE(6,*) '     (-AP.AerLayer.Zmin keyword parameter - Float)'
      WRITE(6,*) '     (-AP.AerLayer.Zmax keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> altitudes of aerosol layer (in km)'
      WRITE(6,*) 
     &'  ==> Case of AP.AerProfile.Type = 2 (IPROFIL in SOS_PROC)'
      GOTO 9998     

 2510 WRITE(6,*) '  SOS_PROC : ERROR_2510 on parameters'
      WRITE(6,*) '     ABSPROFIL parameter in SOS_PROC'
      WRITE(6,*) 
     &'     (-AP.AbsProfile.Type  keyword parameter - Integer)'
      WRITE(6,*) '     --> must be defined'
      WRITE(6,*) '  ==> Absorption profil type'
      WRITE(6,*) '         0 : User profile '
      WRITE(6,*) '         1 : Tropical profile'
      WRITE(6,*) '         2 : MidLat Summer profile'
      WRITE(6,*) '         3 : Midlat Winter profile'
      WRITE(6,*) '         4 : SAS profile'
      WRITE(6,*) '         5 : SAW profile'
      WRITE(6,*) '         6 : USST62'
      WRITE(6,*) '         7 : No absorption'             
      GOTO 9998   
  
 2511 WRITE(6,*) '  SOS_PROC : ERROR_2511 on parameter value'
      WRITE(6,*) '     ABSPROFIL parameter in SOS_PROC'
      WRITE(6,*) 
     &'     (-AP.AbsProfile.Type  keyword parameter - Integer)'
      WRITE(6,*) '  ==> Absorption profil type'
      WRITE(6,*) '         0 : User profile '
      WRITE(6,*) '         1 : Tropical profile'
      WRITE(6,*) '         2 : MidLat Summer profile'
      WRITE(6,*) '         3 : Midlat Winter profile'
      WRITE(6,*) '         4 : SAS profile'
      WRITE(6,*) '         5 : SAW profile'
      WRITE(6,*) '         6 : USST62'
      WRITE(6,*) '         7 : No absorption'             
      GOTO 9998    
    
 2512 WRITE(6,*) '  SOS_PROC : ERROR_2512 on parameters'
      WRITE(6,*) '     FICABSPROFIL parameter in SOS_PROC'
      WRITE(6,*) '     (-AP.AbsProfile.UserFile keyword parameter)'
      WRITE(6,*) '     --> must be defined'    
      WRITE(6,*) '  ==> User absorption profile filename'
      WRITE(6,*) 
     &'  ==> Case of AP.AbsProfile.Type = 0 (ABSPROFIL in SOS_PROC)'
      GOTO 9998     

 2513 WRITE(6,*) '  SOS_PROC : ERROR_2513 on parameters'
      WRITE(6,*) '     ABSPROFIL parameter in SOS_PROC'
      WRITE(6,*) 
     &'     (-AP.AbsProfile.Type  keyword parameter - Integer)'
      WRITE(6,*) 
     &'  ==> Absorption profil type must be set to 7 (no absorption)'
      WRITE(6,*) 
     &'  ==> Case of -AP.AerProfile.Type = 2 (IPROFIl in SOS_PROC)'
      WRITE(6,*) '  ==> Aerosol layer between two altitudes'
      GOTO 9998     

 2514 WRITE(6,*) '  SOS_PROC : ERROR_2514 on parameters'
      WRITE(6,*) '     NUSTEP parameter in SOS_PROC'
      WRITE(6,*) '     (-AP.SpectralResol keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined'   
      WRITE(6,*) '  ==> Spectral resolution of CKD coefficients (cm-1)'
      GOTO 9998  

25141 WRITE(6,*) '  SOS_PROC : ERROR_25141 on parameter value'
      WRITE(6,*) '     NUSTEP parameter in SOS_PROC'
      WRITE(6,*) '     (-AP.SpectralResol keyword parameter - Float)'
      WRITE(6,*) '  ==> Spectral resolution of CKD coefficients'
      WRITE(6,*) '  --> Available values are 1, 5 or 10 cm-1'
      GOTO 9998

 2515 WRITE(6,*) '  SOS_PROC : ERROR_2515 on parameters'
      WRITE(6,*) '     IMODE_CKD_CALCUL parameter in SOS_PROC'
      WRITE(6,*) '     (-SOS.AbsModeCKD  keyword parameter - Integer)'
      WRITE(6,*) '     --> must be defined'   
      WRITE(6,*) '      1 : Fine calculation of the impact of gaz'
      WRITE(6,*) '          absorption on radiance by the CKD method.'
      WRITE(6,*) '      2 : Simplified calculation of the impact of'
      WRITE(6,*) '           gaz absorption radiance.'
      GOTO 9998  
      
 2516 WRITE(6,*) '  SOS_PROC : ERROR_2516 on AOT extrapolated to the',
     &           ' simulation wavelength'
      WRITE(6,*) '  --> AOT = ',TA, ' < 0'
      WRITE(6,*) 
     &'  ==> Check the value for -AER.AOTref (AOT_REF in SOS_PROC)'
      WRITE(6,*) 
     &'  ==> Check the parameters defining the aerosol model'
      GOTO 9998           

C --- Check of the SOS parameters
 2604 WRITE(6,*) '  SOS_PROC : ERROR_2604 on parameter value'
      WRITE(6,*) '     IGMAX parameter in SOS_PROC'
      WRITE(6,*) '     (-SOS.IGmax keyword parameter - Integer)'
      WRITE(6,*) '  --> must be higher than 1' 
      GOTO 9998 

 2605 WRITE(6,*) '  SOS_PROC : ERROR_2605 on parameters'
      WRITE(6,*) '     ITRPHI parameter in SOS_PROC'
      WRITE(6,*) '     (-SOS.View keyword parameter - Integer)'
      WRITE(6,*) '     --> must be defined' 
      WRITE(6,*) '  ==> Azimut mode  '
      GOTO 9998   

 2606 WRITE(6,*) '  SOS_PROC : ERROR_2606 on parameter value'
      WRITE(6,*) '     ITRPHI parameter in SOS_PROC'
      WRITE(6,*) '     (-SOS.View keyword parameter - Integer)'
      WRITE(6,*) '  --> must be equal to 1 or 2' 
      WRITE(6,*) '     1: Azimut plan mode' 
      WRITE(6,*) '     2: Polar diagram azimut mode'       
      GOTO 9998 
 
 2607 WRITE(6,*) '  SOS_PROC : ERROR_2607 on parameters'
      WRITE(6,*) '     PHIOS parameter in SOS_PROC'
      WRITE(6,*) '     (-SOS.View.Phi keyword parameter - Float)'
      WRITE(6,*) '     --> must be defined' 
      WRITE(6,*) 
     &'  ==> Relative azimut required when SOS.View=1 ',
     &' (ITRPHI=1 in SOS_PROC)'  
      GOTO 9998     
 
 2608 WRITE(6,*) '  SOS_PROC : ERROR_2608 on parameters'
      WRITE(6,*) '     PAS_PHI parameter in SOS_PROC'
      WRITE(6,*) '     (-SOS.View.Dphi keyword parameter - Integer)'
      WRITE(6,*) '     --> must be defined' 
      WRITE(6,*) '  ==> Azimut sampling step required when SOS.View=2'
      WRITE(6,*) '      (ITRPHI=2 in SOS_PROC)'    
      GOTO 9998   

 2609 WRITE(6,*) '  SOS_PROC : ERROR_2609 on parameter value'
      WRITE(6,*) '     PAS_PHI parameter in SOS_PROC'
      WRITE(6,*) '     (-SOS.View.Dphi keyword parameter - Integer)'
      WRITE(6,*) '  --> must be a positive value'
      WRITE(6,*) '  ==> Case SOS.View = 2 (ITRPHI in SOS_PROC)'      
      WRITE(6,*) '      Polar diagram azimut mode'      
      GOTO 9998     
   
 2610 WRITE(6,*) '  SOS_PROC : ERROR_2610 on parameter value'
      WRITE(6,*) '     IPOLAR parameter in SOS_PROC'
      WRITE(6,*) '     (-SOS.Ipolar keyword parameter - Integer)'
      WRITE(6,*) '  -->  must be equal to 0 or 1 ' 
      GOTO 9998 

 2612 WRITE(6,*) '  SOS_PROC : ERROR_2612 on parameter value'
      WRITE(6,*) '     IDIRSURF parameter in SOS_PROC'
      WRITE(6,*) '     (-SOS.DirectSurf keyword parameter - Integer)'
      WRITE(6,*) '  -->  must be equal to 0 or 1 ' 
      GOTO 9998 

 2611 WRITE(6,*) '  SOS_PROC : ERROR_2611 on parameter value'
      WRITE(6,*) '     ZOUT parameter in SOS_PROC'
      WRITE(6,*) '     (-SOS.OutputAlt keyword parameter - Float)'
      WRITE(6,*) '  -->  must be equal to:' 
      WRITE(6,*) 
     &'       -1 : standard output (Up TOA, Down BOA)'
      WRITE(6,*) 
     &'       Zout : output level altitude (in km) (Up & Down at Zout)'
      WRITE(6,*) 
     &'          for 0 < Zout < CTE_TOA_ALT (constant): TOA altitude' 
      WRITE(6,*) 
     &'          with CTE_TOA_ALT =',CTE_TOA_ALT            
      GOTO 9998 
   
C --- Definition of the trace files

 3001 WRITE(6,*) '  SOS_PROC : ERROR_3001 : '
      WRITE(6,*) '  => Angles Log filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	
   
 3002 WRITE(6,*) '  SOS_PROC : ERROR_3002 : '
      WRITE(6,*) '  => Mie Log filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	

 3003 WRITE(6,*) '  SOS_PROC : ERROR_3003 : '
      WRITE(6,*) '  => Aerosol Log filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	

 3004 WRITE(6,*) '  SOS_PROC : ERROR_3004 : '
      WRITE(6,*) '  => Profil Log filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	
   
 3005 WRITE(6,*) '  SOS_PROC : ERROR_3005 : '
      WRITE(6,*) '  => Surface Log filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	

 3006 WRITE(6,*) '  SOS_PROC : ERROR_3006 : '
      WRITE(6,*) '  => SOS Log filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	

C --- Definition of the output files

 3007 WRITE(6,*) '  SOS_PROC : ERROR_3007 : '
      WRITE(6,*) '  => Radiance Angles filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	

 3008 WRITE(6,*) '  SOS_PROC : ERROR_3008 : '
      WRITE(6,*) '  => Mie Angles filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	

 3009 WRITE(6,*) '  SOS_PROC : ERROR_3009 : '
      WRITE(6,*) '  => Aerosol Phaze function filename: maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	          	

 3011 WRITE(6,*) '  SOS_PROC : ERROR_3011 : '
      WRITE(6,*) '  => Radiance binary filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	

30111 WRITE(6,*) '  SOS_PROC : ERROR_30111 : '
      WRITE(6,*) '  => Cte CTE_REMOVE_PREVIOUS_BIN_FILE is set to 0:'
      WRITE(6,*) '    --> Remove this file' 
      WRITE(6,*) '    --> or change the CTE_RESULT directory' 
      WRITE(6,*) ''   
      GOTO 9998
        
 3012 WRITE(6,*) '  SOS_PROC : ERROR_3012 : '
      WRITE(6,*) '  => Transmittance filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998
      
30121 WRITE(6,*) '  SOS_PROC : ERROR_30121 : '
      WRITE(6,*) '  => Flux filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998   

C --- Opening of the trace files

 3104 WRITE(6,*) '  SOS_PROC : ERROR_3104 : '
      WRITE(6,*) '  => Error opening Profile log file: ',FICPROFIL_LOG
      GOTO 9998

 3106 WRITE(6,*) '  SOS_PROC : ERROR_3106 : '
      WRITE(6,*) '  => Error opening SOS log file: ',FICSOS_LOG
      GOTO 9998          
    
C --- Run of the routines
 3500 WRITE(6,*) '  SOS_PROC : ERROR_3500 : '
      WRITE(6,*) '  => The summation of AIK coefs is not 1 : ',sumAIK
      WRITE(6,*) '  => Check the CKD coefficients'
      GOTO 9998          
  
C --- Run of the routines

 4000 WRITE(6,*) '  SOS_PROC : ERROR_4000 on subroutine SOS_ANGLES'
      GOTO 9998 

 4100 WRITE(6,*) '  SOS_PROC : ERROR_4100 on subroutine SOS_AEROSOLS'
      GOTO 9998 

 4200 WRITE(6,*) '  SOS_PROC : ERROR_4200 on subroutine SOS_SURFACE'
      GOTO 9998 

 4300 WRITE(6,*) '  SOS_PROC : ERROR_4300 on subroutine 
     &            SOS_PREPA_ABSPROFILE'
      GOTO 9998 

 4400 WRITE(6,*) '  SOS_PROC : ERROR_4400 on subroutine SOS_PREPA_OS'
      GOTO 9998 

 4500 WRITE(6,*) '  SOS_PROC : ERROR_4500 on subroutine SOS_ABSPROFILE'
      GOTO 9998 

 4600 WRITE(6,*) '  SOS_PROC : ERROR_4600 on subroutine SOS_PROFIL'
      GOTO 9998 

 4700 WRITE(6,*) '  SOS_PROC : ERROR_4700 on subroutine SOS'
      GOTO 9998 

 4800 WRITE(6,*) '  SOS_PROC : ERROR_4800 on subroutine SOS_AGGREGATE'
      GOTO 9998 

 4900 WRITE(6,*) 
     &  '  SOS_PROC : ERROR_4900 on subroutine SOS_TRPHI_OPTION'
      GOTO 9998 

C -- Opening / writing of the Transmission and Flux files
  997 WRITE(6,*) '  SOS_PROC : ERROR_997 : '
      WRITE(6,*) '  => ERROR on transmission file opening'
      GOTO 9998         

  998 WRITE(6,*) '  SOS_PROC : ERROR_998 : '
      WRITE(6,*) '  => ERROR on transmission file writing '
      GOTO 9998         

  987 WRITE(6,*) '  SOS_PROC : ERROR_987 : '
      WRITE(6,*) '  => ERROR on flux file opening'
      GOTO 9998         

  988 WRITE(6,*) '  SOS_PROC : ERROR_988 : '
      WRITE(6,*) '  => ERROR on flux file writing'
      GOTO 9998         

      CALL FLUSH(6)             !Forces the flush of the buffer to the screen

       
 9998 IER=1
 
      
 9999 RETURN



C =========================================================
C --- Formats
C =========================================================            
C


C --- Formats for the trace file: inputs of SOS_AEROSOLS
  807 FORMAT(28h       FICANGLES_RES_LUM :  ,A)
  841 FORMAT(45h       Simulation wavelength (Wa)          : ,F12.5)
  842 FORMAT(45h       AOT                                 : ,F12.5) 
  843 FORMAT(45h       Reference wavelength (Waref)        : ,F12.5) 
  844 FORMAT(45h       Reference AOT                       : ,F12.5) 
  845 FORMAT(45h       Truncature option                   : ,I4) 
  846 FORMAT(45h       Type of aerosols model              : ,I4)
  847 FORMAT(45h       MMD : Refractive index for Wa       : ,
     &      (F12.3,2X,F12.5))
c$$$  848 FORMAT(45h       MMD : Refractive index for Waref    : ,
c$$$     &      (F12.3,2X,F12.5))  
  849 FORMAT(45h       Type of size distribution           : ,I4)
  850 FORMAT(45h       Size distribution parameters        : ,
     &      (F12.5,2X,F12.5,2X,F12.5))    
  851 FORMAT(45h       Type of WMO model                   : ,I4)  
  852 FORMAT(45h       WMO component conc. (DL, WS, OC, SO): ,
     &      4(F12.5,2X))    
  853 FORMAT(45h       Type of Shettle & Fenn model        : ,I4)    
  854 FORMAT(45h       S&F model : relative humidity       : ,F12.5)
  855 FORMAT(45h       BMD : Choice of mixture description : ,I4) 
  856 FORMAT(45h       BMD : user conc. Coarse / Fine      : ,
     &      (F12.5,2X,F12.5))
  857 FORMAT(45h       BMD : AOT coarse / AOTtot for Waref : ,F12.5)   
  858 FORMAT(45h         Refractive index for Wa           : ,
     &      (F12.3,2X,F12.5))      
  859 FORMAT(45h         Refractive index for Waref        : ,
     &      (F12.3,2X,F12.5))      
  860 FORMAT(45h         Size distr parameters             : ,
     &      (F12.5,2X,F12.5))  
  861 FORMAT(45h       User's file defining aerosols IOP   : ,A)  
 8611 FORMAT(45h       User's file for aerosols mixture    : ,A)  
  862 FORMAT(45h       Log file                            : ,A)  
  863 FORMAT(45h       Res file                            : ,A) 
  864 FORMAT(45h       Storage directory for MIE files     : ,A)  
  865 FORMAT(45h       Mie computations log file           : ,A)    

C --- Formats for the Transmission and Flux file:
1005  FORMAT(20hSolar Zenith Angle :,1x,f7.3)
1006  FORMAT(36hDirect transmission TOA -> surface :,1x,f8.4)
1010  FORMAT(12h    thetas =,1x,f6.3,3x,12htd(thetas) =,1x,f7.4)  
2010  FORMAT(12h    thetav =,1x,f6.3,3x,12htd(thetav) =,1x,f7.4)
2016  FORMAT(34h   - Downward direct flux at BOA :,1x,f9.5)
2017  FORMAT(34h   - Downward diffuse flux at BOA:,1x,f9.5)
2018  FORMAT(34h   ==> Downward total flux at BOA:,1x,f9.5)
2020  FORMAT(f7.2,2x,f7.4,1x,f7.4,1x,f7.4,1x,f7.4) 

     
      END     !END OF PROGRAM SOS_PROC
