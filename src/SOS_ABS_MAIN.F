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
C* FILE: SOS_ABS_MAIN.F
C* PROJECT: Successive orders of scattering
C* ROLE: Launch program of the OS code
C*           - manages the inputs / outputs and the processing options,
C*           - runs successively the different modules of the OS
C*           - the radiances are computed for each product of the
C*             decomposition of the K-Correlated method
C*
C*
C* AUTHOR:
C*   Initial code: Laboratoire d'Optique Atmospherique (LOA, Villeneuve d'Ascq).
C*   Code rework and development: CS GROUP France, B.Lafrance, X.Lenot
C*
C*
C* MOD:VERSION:1.0: DATE: 06/06/2013
C*    Creation of the function (replaces the former script main_SOS.ksh)
C*
C* MOD:VERSION:1.1: DATE: 22/02/2016
C*   - Added test on the existence of the file FICSOS_RES_BIN (-SOS.ResBin)
C*     to avoid, when a processing is restarted, continuing to fill
C*     this file.
C*
C*    - Removal of the formats when reading the input parameters in favour
C*      of READ(PROG_ARG,*,ERR=NNN)
C*      ==> Allows a compilation with gfortran and avoids wrong conversions
C*          if the input value does not strictly respect the expected format
C*          (for instance, an integer N given instead of a real N.
C*           leading to a wrongly interpreted value).
C*
C*    - Minor adjustments to strictly respect the 72 columns.
C*    - Minor syntax corrections of error messages
C*      (AP.HA corrected to AP.AerHS.HA)
C*    - No aerosol LOG file by default.
C*
C*    - Array NEXP changed from Double precision to Integer
C*
C*    - Change of a check of TA (AOT) after spectral extrapolation
C*
C*    - Management of the case where the pressure is not given as input (assignment
C*      of the default value)
C*      
C*    - Correction of the passing of the refractive indices to routine SOS_AEROSOLS
C*      for the reference wavelength (first call of SOS_AEROSOLS).
C*
C*    - Correction of the management of the passing of the truncation option index
C*      to SOS_AEROSOLS to avoid switching the option to zero if the truncation is
C*      cancelled by the code.
C*
C*    - Adjustment to correctly locate the files FICUP_ANGLES_USER
C*      and FICDOWN_ANGLES_USER in the directory tree of the result files.
C*
C*    - Change of the management of the variances for the LND size distributions:
C*      the use of SIG = log10(variance) is abandoned in favour of a
C*      classical formulation of the LND model versus the variance.
C*		N(r) = exp(-0.5 *(log10(r/rm)/sig)**2) / r
C*	becomes
C*		N(r) = exp(-0.5 *(log(r/rm)/sig)**2) / r
C*      ==> Impacts the values of the input parameters of SOS_AEROSOLS and thus
C*          given on the command line to SOS_ABS_MAIN.
C*
C*
C* MOD:VERSION:1.2: DATE: 10/12/2018
C*   - Fix to no longer compute systematically the transmission file
C*     even when the option is not requested (generation of an
C*     unwanted NO_OUTPUT file).
C*
C*   - Removal of the reading of the parameter AP.ResFile defining the name of the
C*     generated profile file, since the profile is specific to each term of the
C*     CKD decomposition. Storing it is misleading about its content.
C*     Use of a temporary profile file
C*     (default name, file deleted at the end of the processing).
C*
C*   - By default the file FICABSPROFIL is named "NO_USER_ABS_PROFILE_FILE"
C*
C*
C* MOD:VERSION:1.3: DATE: 01/08/2019
C*   - Fix to avoid a conflict, between the main program
C*     SOS_ABS_MAIN and the routine SOS_TRPHI_OPTION, on the declaration
C*     of the name of a radiance output file for the user angles
C*     when the option is not requested (i.e. when there is no definition
C*     of -SOS.ResFileUp.UserAng or -SOS.ResFileDown.UserAng but
C*     -ANG.Rad.UserAngFile is defined).
C*     ==> The initialisation of the names of the files FICUP_ANGLES_USER
C*         and FICDOWN_ANGLES_USER changes from "NO_USER_ANGLES" to "NO_OUTPUT".
C*     ==> Adjustment of the test completing the name of these files.
C*
C*
C* MOD:VERSION:2.0: DATE: 20/04/2021
C*   - Adjustment for the use of new CKD data
C*     ==> Removal of the distinction between absorptions by fine lines or continuums.
C*     ==> Removal of the tables AA and BB
C*     ==> Adjustment of the tables XKI and AL
C*     ==> Addition of constants defining the dimensions of the tables of CKD coefs
C*     ==> Addition of the new input parameter NUSTEP giving the requested spectral step
C*         of the CKD coefficient files to read.
C*     ==> Abandon of the SO2 (absorption beyond 3.7 microns)
C*         and N2 (absorption for wavelengths below 0.33 microns) continuums
C*     ==> Adjustments to process the NO2 continuum from CKD coefficients
C*     ==> Removal of the test skipping the simulation if the transmission over the
C*         atmospheric column is too low, to keep the scattering in
C*         the upper atmospheric layers.
C*   - The keyword AP.H2o becomes AP.H2O
C*   - Introduction of the possibility to use the Maignan BPDF model.
C*   - Deactivation of the option to use the Nadal BPDF model.
C*   - Added test warning if a read parameter is not one of the
C*     possible parameters.
C*   - Addition of a directory for the temporary processing files (DIR_TMP)
C*     to allow several simultaneous runs from the same directory
C*     ==> temporary file FICPROFIL_TMP: profile file for a unit processing of the routine SOS.
C*     ==> temporary file FICOS_TMP: result of a unit processing of the routine SOS.
C*     ==> temporary file FICOS_AGG_TMP: aggregation of the results of the unit processings.
C*   - Clarification of the use of the mono-modal models
C*     ==> Removal of the possibility to define the Mie file to use (the code defines them by itself).
C*         The parameter -AER.MMD.Mie.Filename is no longer proposed.
C*     ==> Removal of the parameters -SDparam1 and -SDparam2 of the mono-modal models,
C*         in favour of new, more explicit parameters:
C*         -AER.MMD.LNDradius :  Modal radius of the LND.
C*         -AER.MMD.LNDvar : Variance of the LND for a mono-modal distribution.
C*         -AER.MMD.JD.slope : Slope of the Junge's law.
C*         -AER.MMD.JD.rmin  : Minimal radius of the Junge's law.
C*         -AER.MMD.JD.rmax  : Maximal radius of the Junge's law (optional)
C*     ==> Addition of a default value of the maximum radius of a Junge model
C*         (SOS_DEFAULT_AER_JUNGE_RMAX)
C*     ==> Adaptation of the call of routine SOS_AEROSOLS to pass the argument associated with -AER.MMD.JD.rmax
C*     ==> Removal of the parameter -AER.MMD.Mie.AlphaMax
C*         and of the argument ALPHAF passed to routine SOS_AEROSOLS.
C*         The maximum value of the size parameter for the Mie computations is defined by routine SOS_AEROSOLS.
C*    - Addition of a processing option IMOD=5 to use the definition
C*      of a mixture of modes given in a file by the user.
C*      ==> Introduction of the parameter FICMIXTURE_AER, which gives the location of the file.
C*      ==> Adaptation of the call of routine SOS_AEROSOLS.
C*    - Addition of parameters to define the surface concentration of carbon dioxide (AP.CO2)
C*      and methane (AP.CH4).
C*
C* MOD:VERSION:2.1: 31/08/2022
C*    - Constants renamed from SOS_* to CTE_*
C*    - Removal of the arguments TRACEOS and IDLOGOS in the call of SOS_AGGREGATE
C*    - Traces converted to English.
C*
C* MOD:VERSION:2.2: 24/11/2022
C*    - Rework of the program SOS_ABS_MAIN to introduce the routine SOS_PROC, allowing a call of the code from Python.
C*
C* MOD:VERSION:2.3: 05/04/2023
C*    - Use of the constant CTE_TOA_ALT to define the max altitude (TOA)
C*    
C*  MOD:VERSION:2.4: 19/04/2023
C*    - Addition of the variable IND_ANGOUT_FIN to the output parameters of SOS_PROC
C*
C*  MOD:VERSION:2.5: 10/08/2023
C*    - Addition of the truncation coefficient of the aerosol phase function as output parameter
C*      of the routine SOS_PROC
C*
C******************************************************************************

C----------------------------------------------------------------------------
C Definition of constants
C---------------------------------------------------------------------------- 
C Constants used:
C    CTE_LENDIR : Length of directory names.
C     CTE_LENFIC1 : Length of file names without directory tree.
C     CTE_LENFIC2 : Length of file names with directory tree.
C    CTE_LENKEYWORD : Maximum size of the keywords for argument passing.
C    CTE_HT_STD_PSURF: Standard atmospheric pressure
C    CTE_NOT_DEFINED_VALUE_INT: default value for integers
C    CTE_NOT_DEFINED_VALUE_DBLE: default value for reals
C    CTE_DEFAULT_FICANGLES_RES_LUM: default name of the angle file for the radiance computation
C    CTE_DEFAULT_FICANGLES_RES_MIE: default name of the angle file for the aerosol computation
C    CTE_DEFAULT_FICGRANU: default name of the size distribution file
C    CTE_DEFAULT_RESBIN: default name of the binary result file
C    CTE_DEFAULT_RESUP: default name of the upward radiance file
C    CTE_DEFAULT_RESDOWN: default name of the downward radiance file
C    CTE_OS_NBMU_MAX  : Maximum number of positive Gauss angles for the radiance computations
C    CTE_DEFAULT_AER_JUNGE_RMAX : Maximum radius (in microns) of a Junge model size distribution
C    CTE_DEFAULT_IGMAX : Default maximum number of scattering orders

C----------------------------------------------------------------------------
#include "SOS.h"
#define INCTE_PI DACOS(-1.D+00)
#define INCTE_IDLOG_MAIN 77

C==============================================================================
C PROGRAM: SOS_ABS_MAIN
C ==========  
C
C Main program of the OS code with absorption
C
C
C Input data of the program
C ------------------------------
C    -->Parameters/arguments of the program
C
C        The arguments are passed as a list of pairs (-Keyword Value)
C
C        Definition of the keywords, mandatory / optional status, format of the associated value:
C
C 	 ------- Arguments specific to the main program
C 	 ------------------------------------------------
C
C         -SOS_Main.ResRoot ==> Mandatory
C             RESROOT   (CHARACTER*CTE_LENDIR) : working directory of the processing
C                                               (complete path)
C 	
C         -SOS_Main.Log  ==> Optional:
C            FICMAIN_LOG (CHARACTER*CTE_LENFIC2) : Name of the trace file of the main program.
C                                   (file name + extension)
C						         ==> No trace if not defined.
C
C         -SOS_Main.Wa ==> Mandatory
C	     PROG_ARG  (CHARACTER*10)  : Wavelength of the radiance simulation
C                                       in microns.
C                                       (decoding format: free --> Float)
C
C
C
C 	 ------- Arguments specific to the definition of the angles
C 	 ----------------------------------------------------	
C
C    Keyword  : -ANG.Rad.NbGauss  
C    ---------------------------   
C       Status : Optional (if not defined, the default value CTE_DEFAULT_NBMU_LUM is used : see SOS.h )       
C       Value :  Number of Gauss angles to be used for radiance computations 
C       Value format : No format constraint --> read as Integer
C
C    Keyword  : -ANG.Rad.UserAngFile  
C    -------------------------------   
C       Status : Optional      
C       Value :  Filename of the complementary list of user's angles to complete the 
C                ANG.Rad.NbGauss angles (complete path). 
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C
C    Keyword  : -ANG.Thetas  
C    --------------------------   
C       Status : Required     
C       Solar zenithal angle (degrees) 
C       Value format : No format constraint --> read as Float
C
C    Keyword  : -ANG.Rad.ResFile  
C    ---------------------------   
C       Status : Optional (if not defined, the default value CTE_DEFAULT_FICANGLES_RES_LUM is used : see SOS.h )      
C       Value :  Filename of list of angles used to BRDF/BPDF and radiance computations. 
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C
C    Keyword  : -ANG.Aer.NbGauss  
C    ---------------------------   
C       Status : Optional (if not defined, the default value CTE_DEFAULT_NBMU_MIE is used : see SOS.h )       
C       Value :  Number of Gauss angles to be used for phaze matrix computations 
C       Value format : No format constraint --> read as Integer
C
C    Keyword  : -ANG.Aer.UserAngFile 
C    -------------------------------   
C       Status : Optional      
C       Value :  Filename of the complementary list of user's angles to complete the 
C                ANG.Aer.NbGauss angles (complete path). 
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C
C    Keyword  : -ANG.Aer.ResFile  
C    ----------------------------   
C       Status : Optional (if not defined, the default value CTE_DEFAULT_FICANGLES_RES_MIE is used : see SOS.h )       
C       Value :  Filename of list of angles used to the matrix phase function computations. 
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C
C    Keyword  : -ANG.Log
C    -------------------   
C       Status : Optional
C       Value :  log filename for ANGLES computations. 
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C        ==> Already existing file will be overwritted
C
C
C 	 ------- Arguments specific to the definition of the surface model
C 	 ---------------------------------------------------------------	
C
C    Keyword  : -SURF.Dir
C    -------------------   
C       Status : Conditional (required if a surface BRDF/BPDF is simulated)
C       Value :  Storage directory for SURFACE files producted by SOS_SURFACE computations (complete path). 
C       Value format : Encoding as string (CTE_LENDIR characters max : see SOS.h).
C
C    Keyword  : -SURF.Type  
C    ---------------------   
C       Status : Required      
C       Value :  Kind of surface 
C      		0 : lambertian surface for an albedo RHO.
C			1 : lambertian surface + glitter (roughness water).
C			2 : lambertian surface + Fresnel's reflexion (plane water).
C		       3 : lambertian surface + Roujean's BRDF.
C			4 : lambertian surface + Roujean's BRDF + Rondeaux's BPDF.
C			5 : lambertian surface + Roujean's BRDF + Breon's BPDF .
C			6 : lambertian surface + Roujean's BRDF + Nadal's BPDF.
C       Value format : No format constraint --> read as Integer.
C
C    Keyword  : -SURF.Log  
C    --------------------   
C       Status : Optional    
C       Value :  log filename for SOS_SURFACE computations.
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C        ==> Already existing file will be overwritted
C
C    Keyword  : -SURF.Ind  
C    --------------------   
C       Status : Associated to  -SURF.Type 1 2 4 5 6 or 7     
C       Value :  Surface / atmosphere refractive index. 
C       Value format : No format constraint --> read as Float
C
C    Keyword  : -SURF.Glitter.Wind  
C    -----------------------------   
C       Status : Associated to  -SURF.Type 1    
C       Value :  Wind velocity (m/s). 
C       Value format : No format constraint --> read as Float
C
C    Keyword  : -SURF.Roujean.K0  
C    ---------------------------   
C       Status : Associated to  -SURF.Type 3 4 5 6 or 7      
C       Value :  Parameter K0 of Roujean's BRDF model. 
C       Value format : No format constraint --> read as Float
C
C    Keyword  : -SURF.Roujean.K1  
C    --------------------------- 
C       Status : Associated to  -SURF.Type 3 4 5 6 or 7      
C       Value :  Parameter K1 of Roujean's BRDF model. 
C       Value format : No format constraint --> read as Float
C
C    Keyword  : -SURF.Roujean.K2  
C    ---------------------------   
C       Status : Associated to  -SURF.Type 3 4 5 6 or 7      
C       Value :  Parameter K2 of Roujean's BRDF model. 
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -SURF.Nadal.Alpha 
C    ----------------------------   
C       Status : Associated to  -SURF.Type 6     
C       Value :  Parameter Alpha of Nadal's BPDF model. 
C       Value format : No format constraint --> read as Float
C
C    Keyword  : -SURF.Nadal.Beta  
C    ---------------------------   
C       Status : Associated to  -SURF.Type 6      
C       Value :  Parameter Beta of Nadal's BPDF model. 
C       Value format : No format constraint --> read as Float
C
C    Keyword  : -SURF.Maignan.C   
C    ---------------------------   
C       Status : Associated to  -SURF.Type 7      
C       Value :  C*exp(-NDVI) parameter of the Maignan's surface BPDF model. 
C       Value format : No format constraint --> read as Float
C
C    Keyword  : -SURF.Alb 
C    --------------------   
C       Status : Required      
C       Value :  Surface albedo for the wavelength of radiance simulation
C       Value format : No format constraint --> read as Float
C
C    Keyword  : -SURF.File  
C    ---------------------   
C       Status : Optional. Associated to  -SURF.Type 1 3 4 5 6 or 7     
C       Value :  User filename for the surface reflexion matrix (complete path).
C                --> Glitter or others BRDF/BPDF. 
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C
C 	 ------- Arguments specific to the definition of the aerosol models
C 	 ----------------------------------------------------------------	
C
C    Keyword  : -AER.Waref  
C    ---------------------   
C       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 2 (for common case with AER.AOT > 0)    
C       Value :  Reference wavelength for the aerosol optical thickness (microns).
C       Value format : No format constraint --> read as Float
C	  
C    Keyword  : -AER.AOTref  
C    ----------------------   
C       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 2 (for common case with AER.AOT > 0)    
C       Value :  Aerosol optical thickness for the reference wavelength.
C       Value format : No format constraint --> read as Float
C     		       
C    Keyword  : -AER.Tronca 
C    ----------------------   
C       Status : Default value = 1 (for common case with AER.AOT > 0)     
C       Value :  Option for aerosol phaze function troncation 
C                (1 by default to apply a truncation, 0 to skip the truncation)
C       Value format : No format constraint --> read as Integer.
C
C    Keyword  : -AER.Log
C    -------------------   
C       Status : Optional (for common case with AER.AOT > 0)     
C       Value :  log filename for SOS_AEROSOLS computations (complete path).
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C        ==> Already existing file will be overwritted
C
C    Keyword  : -AER.MieLog  
C    ----------------------   
C       Status : Optional  
C       Value :  log filename for MIE computations (complete path).
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C        ==> Already existing file will be overwritted
C
C    Keyword  : -AER.DirMie
C    ----------------------   
C       Status : Required      
C       Value :  Storage directory for MIE files producted by SOS_AEROSOLS computations (complete path). 
C       Value format : Encoding as string (CTE_LENDIR characters max : see SOS.h).
C
C    Keyword  : -AER.ResFile
C    -----------------------   
C       Status : Optional (if not defined, neither AER.UserFile, 
C                          the default value CTE_DEFAULT_FICGRANU is used : see SOS.h )   
C       Value :  Filename of the result SOS_AEROSOLS computations. 
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C
C    Keyword  : -AER.UserFile
C    -----------------------   
C       Status : Optional (required if AER.ResFile is not defined)      
C       Value :  Filename of the user's aerosol filename (complete path). 
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C
C    Keyword  : -AER.Model 
C    ----------------------   
C       Status : Required (for common case with AER.AOT > 0)     
C       Value :  Type of size distribution (mono or multi-mode)
C      		        0 : Mono-modal models.
C			        1 : WMO models.
C			        2 : Shettle & Fenn models.
C		               3 : LND bi-modal models.
C                            4 : Use of external phase functions
C       Value format : No format constraint --> read as Integer.
C
C					
C    Keyword  : -AER.MMD.MRwa  
C    ------------------------ 
C       Status : Associated to -AER.Model 0 (for common case with AER.AOT > 0)     
C       Value :  Real part of the refractive index at the wavelength of aerosol properties calculation.
C       Value format : No format constraint --> read as Float
C	  
C    Keyword  : -AER.MMD.MIwa  
C    ------------------------ 
C       Status : Associated to -AER.Model 0 (for common case with AER.AOT > 0)     
C       Value :  Imaginary part of the refractive index for the wavelength of aerosol properties calculation.
C       Value format : No format constraint --> read as Float
C
C    Keyword  : -AER.MMD.MRwaref  
C    ------------------------ 
C       Status : Associated to -AER.Model 0 (for common case with AER.AOT > 0)     
C       Value :  Real part of the refractive index at the reference wavelength of aerosol properties calculation.
C       Value format : No format constraint --> read as Float
C	  
C    Keyword  : -AER.MMD.MIwaref  
C    ------------------------ 
C       Status : Associated to -AER.Model 0 (for common case with AER.AOT > 0)     
C       Value :  Imaginary part of the refractive index for the reference wavelength of aerosol properties calculation.
C       Value format : No format constraint --> read as Float
C	  
C    Keyword  : -AER.MMD.SDtype 
C    --------------------------   
C       Status : Associated to -AER.Model 0 (for common case with AER.AOT > 0)  
C       Value :  Type of mono-modal size distribution 
C			1 : LND.
C			2 : Junge's law.
C       Value format : No format constraint --> read as Integer.			
C	    
C    Keyword  : -AER.MMD.LNDradius  
C    ----------------------------- 
C       Status : Associated to -AER.Model 0 and -AER.MMD.SDtype 1 (for common case with AER.AOT > 0)     
C       Value :  Modal radius of the LND (microns).
C			-AER.MMD.SDtype 2 : minimal radius of the Junge's law (microns).
C       Value format : No format constraint --> read as Float	  
C	    
C    Keyword  : -AER.MMD.LNDvar  
C    -------------------------- 
C       Status : Associated to -AER.Model 0 and -AER.MMD.SDtype 1 (for common case with AER.AOT > 0)     
C       Value :  Variance of the LND (microns) for a mono-modal distribution.
C       Value format : No format constraint --> read as Float	  
C	  
C    Keyword  : -AER.MMD.JD.slope 
C    ---------------------------- 
C       Status : Associated to -AER.Model 0 and -AER.MMD.SDtype 2 (for common case with AER.AOT > 0)     
C       Value :  Slope of the Junge's law.
C       Value format : No format constraint --> read as Float	  
C	  
C    Keyword  : -AER.MMD.JD.rmin 
C    ---------------------------- 
C       Status : Associated to -AER.Model 0 and -AER.MMD.SDtype 2 (for common case with AER.AOT > 0)     
C       Value :  Minimal radius of the Junge's law (microns).
C       Value format : No format constraint --> read as Float	  
C	  
C    Keyword  : -AER.MMD.JD.rmax 
C    ---------------------------- 
C       Status : Associated to -AER.Model 0 and -AER.MMD.SDtype 2 (for common case with AER.AOT > 0)     
C       Value :  Maximal radius of the Junge's law (microns). 
C                If not defined : the default value CTE_DEFAULT_AER_JUNGE_RMAX is used
C       Value format : No format constraint --> read as Float	  
C  
C    Keyword  : -AER.WMO.Model  
C    ------------------------- 
C       Status : Associated to -AER.Model 1 (for common case with AER.AOT > 0)     
C       Value :  Type of WMO model.
C      		        1 : Continental WMO model.
C			        2 : Maritime WMO model.
C			        3 : Urban WMO model.
C		               4 : WMO model by user definition.
C       Value format : No format constraint --> read as Integer.	  
C
C    Keyword  : -AER.WMO.DL  
C    ---------------------- 
C       Status : Associated to -AER.Model 1 and -AER.WMO.Model 4 (for common case with AER.AOT > 0)     
C       Value :  Volumetric concentration (%) for WMO Dust-Like particles.
C       Value format : No format constraint --> read as Float
C
C    Keyword  : -AER.WMO.WS  
C    ---------------------- 
C       Status : Associated to -AER.Model 1 and -AER.WMO.Model 4 (for common case with AER.AOT > 0)     
C       Value :  Volumetric concentration (%) for WMO Water Soluble particles.
C       Value format : Encoding as real F9.5 (10 characters max).
C
C    Keyword  : -AER.WMO.OC  
C    ---------------------- 
C       Status : Associated to -AER.Model 1 and -AER.WMO.Model 4 (for common case with AER.AOT > 0)     
C       Value :  Volumetric concentration (%) for WMO OCeanic particles.
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AER.WMO.SO  
C    ---------------------- 
C       Status : Associated to -AER.Model 1 and -AER.WMO.Model 4 (for common case with AER.AOT > 0)     
C       Value :  Volumetric concentration (%) for WMO SOot particles.
C       Value format : No format constraint --> read as Float.	  
C	  
C    Keyword  : -AER.SF.Model  
C    ------------------------ 
C       Status : Associated to -AER.Model 2 (for common case with AER.AOT > 0)     
C       Value :  Type of Shettle & Fenn model.
C      		1 : Tropospheric S&F model.
C			2 : Urban S&F model.
C			3 : Maritime S&F model.
C		       4 : Coastal S&F model.
C       Value format : No format constraint --> read as Integer.	  
C
C    Keyword  : -AER.SF.HR  
C    --------------------- 
C       Status : Associated to -AER.Model 2 (for common case with AER.AOT > 0)     
C       Value :  Relative humidity (%) for Shettle & Fenn model.
C       Value format : No format constraint --> read as Float.
C	  
C    Keyword  : -AER.BMD.VCdef 
C    ---------------------   
C       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)  
C       Value :  Type of selected method for bi-modal LND volumetric concentration definition
C			1 : Use of user volumetric concentrations.
C			2 : Use of user ratio AOT_coarse / AOT_total.
C       Value format : No format constraint --> read as Integer.			
C
C    Keyword  : -AER.BMD.CoarseVC  
C    ---------------------------- 
C       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 1 (for common case with AER.AOT > 0)     
C       Value :  User volumetric concentration for the "LND coarse mode".
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AER.BMD.FineVC  
C    -------------------------- 
C       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 1 (for common case with AER.AOT > 0)     
C       Value :  User volumetric concentration for the "LND fine mode".
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AER.BMD.RAOT  
C    ------------------------ 
C       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 2 (for common case with AER.AOT > 0)     
C       Value :  User value of the ration AOT_coarse / AOT_total for the aerosol reference wavelength.
C       Value format : No format constraint --> read as Float.
C	  	  
C    Keyword  : -AER.BMD.CM.MRwa  
C    ---------------------------   
C       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)    
C       Value :  Real part of the refractive index for the "LND coarse mode" 
C                at the wavelength of aerosol properties calculation.
C       Value format : No format constraint --> read as Float.
C	  
C    Keyword  : -AER.BMD.CM.MIwa  
C    ---------------------------   
C       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)    
C       Value :  Imaginary part of the refractive index for the "LND coarse mode" 
C                at the wavelength of aerosol properties calculation.
C       Value format : No format constraint --> read as Float.
C	  
C    Keyword  : -AER.BMD.CM.MRwaref  
C    ------------------------------   
C       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 2 (for common case with AER.AOT > 0)    
C       Value :  Real part of the refractive index for the "LND coarse mode" 
C                at the aerosol reference wavelength.
C       Value format : No format constraint --> read as Float.
C	  
C    Keyword  : -AER.BMD.CM.MIwaref  
C    -----------------------------   
C       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 2 (for common case with AER.AOT > 0)    
C       Value :  Imaginary part of the refractive index for the "LND coarse mode" 
C                at the aerosol reference wavelength.
C       Value format : No format constraint --> read as Float.
C	  
C    Keyword  : -AER.BMD.CM.SDradius  
C    -------------------------------   
C       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)    
C       Value :  Modal radius of the "LND coarse mode" (microns).
C       Value format : No format constraint --> read as Float.
C	  
C    Keyword  : -AER.BMD.CM.SDvar  
C    ----------------------------   
C       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)    
C       Value :  Log10 of the "LND coarse mode" variance.
C       Value format : No format constraint --> read as Float.
C	  
C    Keyword  : -AER.BMD.FM.MRwa  
C    ---------------------------   
C       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)    
C       Value :  Real part of the refractive index for the "LND fine mode" 
C                at the wavelength of aerosol properties calculation.
C       Value format : No format constraint --> read as Float.
C	  
C    Keyword  : -AER.BMD.FM.MIwa  
C    ---------------------------   
C       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)    
C       Value :  Imaginary part of the refractive index for the "LND fine mode" 
C                at the wavelength of aerosol properties calculation.
C       Value format : No format constraint --> read as Float.
C	  
C    Keyword  : -AER.BMD.FM.MRwaref  
C    ------------------------------   
C       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 2 (for common case with AER.AOT > 0)    
C       Value :  Real part of the refractive index for the "LND fine mode" 
C                at the aerosol reference wavelength.
C       Value format : No format constraint --> read as Float.
C	  
C    Keyword  : -AER.BMD.FM.MIwaref  
C    -----------------------------   
C       Status : Associated to -AER.Model 3 and -AER.BMD.VCdef 2 (for common case with AER.AOT > 0)    
C       Value :  Imaginary part of the refractive index for the "LND fine mode" 
C                at the aerosol reference wavelength.
C       Value format : No format constraint --> read as Float.
C	  
C    Keyword  : -AER.BMD.FM.SDradius  
C    -------------------------------   
C       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)    
C       Value :  Modal radius of the "LND fine mode" (microns).
C       Value format : No format constraint --> read as Float.
C	  
C    Keyword  : -AER.BMD.FM.SDvar  
C    ----------------------------   
C       Status : Associated to -AER.Model 3 (for common case with AER.AOT > 0)    
C       Value :  Log10 of the "LND fine mode" variance.
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AER.ExtData
C    -----------------------   
C       Status : Associated to -AER.Model 4 (for common case with AER.AOT > 0)       
C       Value :  Filename (complete path) of user's external phaze functions 
C                and radiative parameters (extinction and scattering coefficients)
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C       Special requirement : the reference aerosol wavelength (-AER.Waref) and 
C                             the radiance simulation wavelength (-SOS.Wa) have to be equal.	 
C
C    Keyword  : -AER.DefMixture
C    --------------------------   
C       Status : Associated to -AER.Model 5 (for common case with AER.AOT > 0)       
C       Value :  Filename (complete path) of user's file defining the mixture of aerosol models.
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C    
C
C 	 ------- Arguments specific to the definition of the profiles
C 	 ------------------------------------------------------------------	
C
C    Keyword  : -AP.Log
C    ------------------   
C       Status : Optional
C       Value :  log filename for SOS_PROFIL computations.
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C        ==> Already existing file will be overwritted
C
C    Keyword  : -AP.MOT  
C    ------------------   
C       Status : Optional      
C       Value :  Molecular optical thickness for the wavelength of radiance simulation
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AP.HR  
C    ------------------   
C       Status : Associated to -AP.MOT >= 0.    
C       Value :  Molecular optical thickness for the wavelength of radiance simulation
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AP.AerProfile.Type  
C    -------------------   
C       Status : Required      
C       Value :  Type of aerosol profile
C      		     1 : Aerosol profile described by an exponential decrease using a heigth scale
C		     2 : Homogeneous aerosol layer between two altitudes.
C       Value format : No format constraint --> read as Integer.
C
C    Keyword  : -AP.AerHS.HA  
C    -----------------------   
C       Status : Associated to -AP.Type 1 and -AP.AOT >= 0.   
C       Value :  Aerosol heigth scale (km).
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AP.AerLayer.Zmin  
C    ----------------------------   
C       Status : Associated to -AP.Type 2      
C       Value :  Minimal altitude of the aerosol layer (km).
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AP.AerLayer.Zmax  
C    ----------------------------   
C       Status : Associated to -AP.Type 2      
C       Value :  Maximal altitude of the aerosol layer (km).
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AP.Psurf  
C    ----------------------------   
C       Status : Optional (default is McClatchey Value)      
C       Value :  Surface Pressure (mbar).
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AP.H2O  
C    ----------------------------   
C       Status : Optional (default is McClatchey Value)	 
C       Value :  Water vapor content (g/cm2).
C       Value format : No format constraint --> read as Float.
C  
C    Keyword  : -AP.O3  
C    ----------------------------   
C       Status : Optional (default is McClatchey Value)	 
C       Value :  Ozone content (Dobson).
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AP.CO2 
C    ----------------------------   
C       Status : Optional (default is McClatchey Value)	 
C       Value :  Carbon dioxyde concentration at surface level (ppmv).
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AP.CH4  
C    ----------------------------   
C       Status : Optional (default is McClatchey Value)	 
C       Value :  Methane concentration at surface level (ppmv).
C       Value format : No format constraint --> read as Float.
C
C    Keyword  : -AP.AbsProfile.Type  
C    ----------------------------   
C       Status : Required	 
C       Value :  Type of asborption atmospheric profile
C			 0:User
C			 1:TROPICAL
C			 2:MLS
C			 3:MLW
C			 4:SAS
C			 5:SAW
C			 6:USST62
C			 7:No absorption
C       Value format : No format constraint --> read as Integer.
C 
C
C    Keyword  : -AP.SpectralResol  
C    ----------------------------   
C       Status : Conditional (required in case of simulation of absorptions, i.e. if AP.AbsProfile.Type != 7)	 
C       Value :  Spectral resolution in cm-1 of CKD datafiles
C       Value format : No format constraint --> read as Float.
C 

C    Keyword  : -AP.AbsProfile.UserFile  
C    ----------------------------   
C       Status : Required in association to  AP.AbsProfile.Type = 0	 
C       Value :  Filename of user's absorption profile file (complete path)
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h)	 
C
C
C 	 ------- Arguments specific to the Successive Orders computations
C 	 ------- and to the production of the result files
C 	 ---------------------------------------------------------	
C
C    Keyword  : -SOS.ResBin  
C    ----------------------   
C       Status : Optional (if not defined, the default value CTE_DEFAULT_RESBIN is used : see SOS.h)        
C       Value :  Filename of the binary file resulting from SOS computations . 
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C
C    Keyword  : -SOS.ResFileUp  
C    -------------------------   
C       Status : Optional (if not defined, the default value CTE_DEFAULT_RESUP is used : see SOS.h)     
C       Value :  Filename of the ascii file resulting from SOS computations  
C                --> upward radiance.
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C
C    Keyword  : -SOS.ResFileDown  
C    ---------------------------   
C       Status : Optional (if not defined, the default value CTE_DEFAULT_RESDOWN is used : see SOS.h)      
C       Value :  Filename of the ascii file resulting from SOS computations  
C                --> downward radiance.
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C
C    Keyword  : -SOS.ResFileUp.UserAng 
C    ---------------------------------  
C       Status : Optional  (requires the use of -SOS.UserAng.File)       
C       Value :  Filename of the ascii file resulting from SOS computations (complete path). 
C                giving radiance values only for user's angles 
C                --> upward radiance.
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C        ==> Already existing file will be overwritted     
C
C    Keyword  : -SOS.ResFileDown.UserAng 
C    -----------------------------------  
C       Status : Optional  (requires the use of -SOS.UserAng.File)       
C       Value :  Filename of the ascii file resulting from SOS computations (complete path). 
C                giving radiance values only for user's angles 
C                --> downward radiance.
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C        ==> Already existing file will be overwritted     
C
C    Keyword  : -SOS.Log
C    -------------------   
C       Status : Optional
C       Value :  log filename for SOS computations.
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C        ==> Already existing file will be overwritted
C
C    Keyword  : -SOS.AbsModeCKD  
C    --------------------------   
C       Status : Required      
C       Value :  Number defining the CKD mode calculation
C                1 : Computation of the radiance from CKD expansion
C                2 : Computation of the absorption optical thickness from CKD expansion 
C       Value format : No format constraint --> read as Integer
C
C    Keyword  : -SOS.Trans
C    ---------------------  
C       Status : Optional 
C       Value :  Filename of the file given information about the simulated 
C                atmospheric transmissions (will be located in the repertory Working_dir/SOS).
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C        ==> Already existing file will be overwritted
C
C    Keyword  : -SOS.Flux
C    ---------------------  
C       Status : Optional 
C       Value :  Filename of the file given information about the fluxes at BOA and TOA
C                (will be located in the repertory Working_dir/SOS).
C       Value format : Encoding as string (CTE_LENFIC2 characters max : see SOS.h).
C        ==> Already existing file will be overwritted
C
C    Keyword  : -SOS.OutputAlt  
C    ---------------------------   
C       Status : Optionnal   
C       Value :  Altitude (in km) for radiance results.
C                  If not defined : -1 is the default value attributed 
C                  -1 : standard results 
C			--> Upward radiance at TOA
C                       --> Downward radiance at surface level
C                  z  : results for the altitude z of the atmospheric profile
C			(0 <= z <= CTE_TOA_ALT)
C			--> Upward and downward radiance are given for the
C			    altitude z of the profile.
C       Value format : No format constraint --> read as Float.
C			
C    Keyword  : -SOS.Ipolar 
C    ----------------------   
C       Status : Optionnal   
C       Value :  Profile level for radiance results.
C                  If not defined : 1 is the default value attributed
C                  0 : radiance calculation without taking account of polarization. 
C                  1 : radiance calculation with polarization processes.
C       Value format : No format constraint --> read as Integer.
C
C    Keyword  : -SOS.IGmax  
C    ---------------------   
C       Status : Optionnal      
C       Value :  Maximal order of interaction (scattering & surface reflexion). 
C       Value format : No format constraint --> read as Integer.   	 	 	
C
C    Keyword  : -SOS.View  
C    --------------------   
C       Status : Required      
C       Value :  Index for output type
C                  1 : radiance for viewing directions in a constante azimutal plane. 
C                  2 : radiance for viewing directions covering half the space, with a constant step on azimut.  
C       Value format : No format constraint --> read as Integer.   	 	 	
C		 	 	 
C    Keyword  : -SOS.View.Phi  
C    ------------------------   
C       Status : Associated to -SOS.View 1      
C       Value :  Relative azimuth angle (degrees).
C       Value format : No format constraint --> read as Float.  	 	 	
C
C    Keyword  : -SOS.View.Dphi  
C    -------------------------   
C       Status : Associated to -SOS.View 2      
C       Value :  Constant step on azimut (degrees). Integer value.  
C       Value format : No format constraint --> read as Integer.   
C  					
C       
C Results provided by the program
C ----------------------------------
C   Outputs of SOS_PROC:
C      - Normalised radiance fields, by default emerging at TOA and downward at the surface,
C        or upward and downward for the selected level of the profile.
C      - Polarisation data:
C        Parameters Q and U of the Stokes vector, polarisation angle, polarised intensity, degree of polarisation.
C      - Direct and diffuse transmissions, downward and upward.
C      - Direct and diffuse downward flux at the surface level, upward diffuse flux at TOA.
C      - List of indices of the user angles to filter the radiance data.
C      - Truncation coefficient of the aerosol phase function.
C
C   -->Content of the result file for the upward field
C       Ascii file, record format:
C          2(2X,F7.2),2X,3(E13.6,2X),2(F7.2,2X),E13.6
C       or 3(2X,F7.2),1X,3(E13.6,2X),2(1X,F7.2),E13.6
C
C       2 cases: ITRPHI = 1,
C               the viewing plane is obtained for PHIos.
C               Header of 31 lines.
C               Number of records:
C                   one for each value from -NBMU to NBMU (useful angles)
C               Each record contains 8 values:
C                   viewing angle (in degrees), scattering angle (in degrees),
C                   Stokes parameters I,Q,U (in sr-1),
C                   polarisation angle (in degrees), degree of polarisation (in %),
C                   polarisation intensity (in sr-1).
C
C               ITRPHI = 2, 
C               a polar diagram is obtained versus the viewing
C               angle (teta) and the azimuth angle (phi)
C               Header of 26 lines.
C               Number of records:
C                   one for each value from 1 to NBMU (useful angles)
C               Each record contains 9 values:
C                   coordinates phi, teta (in degrees), scattering angle (in degrees),
C                   Stokes parameters I,Q,U (in sr-1),
C                   polarisation angle (in degrees), degree of polarisation (in %),
C                   polarisation intensity (in sr-1).
C
C   -->Content of the result file for the downward field
C      Ascii file, record format (F7.2,2X,3(F13.6,2X))
C      or (2(F7.2,2X),3(F13.6,2X))
C      Same description as for the result file of the upward field
C


C
C   Trace files (optional)
C
C
C Environment variable
C ------------------------
C   SOS_ABS_ROOT : Root of the SOS-ABS code directory
C
C
C Common used:
C --------------
C None
C
C
C Error cases:
C ------------        
C     -Incorrect decoding of parameters
C     -Error in a subroutine
C     -Error when opening the trace file
C     -Invalid parameters
C
C     Message displayed on screen, program stop and return of status 1 to the shell
C   
C
C  Subroutines used:
C --------------------------
C  This program calls the routine:
C     - SOS_PROC
C     - SOS_OUTPUT_HEADER
C     - SOS_OUTPUT_HEADER_POLAR_DIAG
C
C============================================================================


      PROGRAM SOS_ABS_MAIN
      
      IMPLICIT NONE
      

C Constants
C -----------------            
      DOUBLE PRECISION CONVDEGRAD         !Conversion coefficient from degrees to radians.
      
      INTEGER*2 IDLOG			! Identifier number of the trace file
      PARAMETER(IDLOG=INCTE_IDLOG_MAIN)
 



     	            
C* Definition of variables
C*-----------------------------------------------------------------
    
C-- Variables related to the input parameters
      INTEGER*2 NBARGS		 ! Number of input arguments
      INTEGER*4 IARGC     	        ! Type of the function IARGC
      INTEGER*4 NPARAM 		 ! Number of the input parameter
                                    ! (type INTEGER*4 mandatory for f77)
	
      CHARACTER*10 PROG_ARG,DARG,IARG   ! Program arguments
					
      CHARACTER*CTE_LENKEYWORD KEYWORD  ! Keyword of the argument passing.
      INTEGER*2 KEYWORD_VALID	 ! Validity flag of the keyword read
                                    
      
C-- Variables related to the working directory tree

      INTEGER*4 LDIR     	       ! length of the character string RESROOT
      INTEGER*4 LDIRSOS     	! length of the character string DIR_SOS
      INTEGER*4 LFIC    	       ! length of the character string of a file name

      CHARACTER*CTE_LENDIR RESROOT      ! Output directory
      CHARACTER*CTE_LENDIR DIR_SURF     ! Output directory of the surface computations
      CHARACTER*CTE_LENDIR DIR_MIE      ! Output directory of the Mie computations
      CHARACTER*CTE_LENDIR DIR_SOS      ! Output directory of the successive orders
      CHARACTER*CTE_LENDIR DIR_LOG      ! Output directory of the trace files

C-- Global variables

      DOUBLE PRECISION WA_SIMU	! Wavelength of the radiance simulation (in microns)
      
      CHARACTER*CTE_LENFIC2 FICMAIN_LOG ! Name of the main trace file

      INTEGER*4 I,J,K,JJ,IP ! Array indices

      INTEGER*4 IER        	! error code =0 if no error, =-1 otherwise

      LOGICAL TRACE     	! = true if writing in the trace file

                  
C---- Variables related to the definition of the angles
      		
      DOUBLE PRECISION TETAS	! Solar zenith angle in degrees (0 < tetas < 90).
      				       ! (complementary angle: mus= -cos(tetas))
 
      INTEGER*4 NBMU_GAUSS_LUM     ! Number of Gauss angles used for the radiance computations
      INTEGER*4 NBMU_GAUSS_MIE     ! Number of Gauss angles used for the Mie computations

      !Radiance computations:
      CHARACTER*CTE_LENFIC2 FICANGLES_USER_LUM  ! path and name of the "user" angles file
      CHARACTER*CTE_LENFIC2 FICANGLES_RES_LUM   ! path and name of the file of the angles used

      !Phase function computations:
      CHARACTER*CTE_LENFIC2 FICANGLES_USER_MIE  ! path and name of the "user" angles file
      CHARACTER*CTE_LENFIC2 FICANGLES_RES_MIE   ! path and name of the file of the angles used
  
      CHARACTER*CTE_LENFIC2 FICANGLOG           ! name of the trace file


C---- Variables related to the aerosol model
 
      DOUBLE PRECISION AOT_REF	! Aerosol optical thickness at the reference wavelength (not truncated).
      DOUBLE PRECISION WAREF_AOT	! Reference wavelength
    
   	!Mono-modal aerosol case:
      DOUBLE PRECISION RN_WA		! real part of the refractive index at the simulation wavelength
      DOUBLE PRECISION IN_WA		! imaginary part of the refractive index at the simulation wavelength
      DOUBLE PRECISION RN_WAREF	! real part of the refractive index at the reference wavelength
      DOUBLE PRECISION IN_WAREF	! imaginary part of the refractive index at the reference wavelength
      DOUBLE PRECISION LND_RADIUS_MMD_AER ! Modal radius of the LND.
      DOUBLE PRECISION LND_LNVAR_MMD_AER  ! Variance of the LND.
      DOUBLE PRECISION JD_SLOPE_MMD_AER   ! Slope of the Junge law.
      DOUBLE PRECISION JD_RMIN_MMD_AER    ! Minimum radius of the Junge law.
      DOUBLE PRECISION JD_RMAX_MMD_AER    ! Maximum radius of the Junge law.
    	!Volume proportion of the WMO components
      DOUBLE PRECISION C_WMO_DL	    	! User value for the DL particles of the WMO model
      DOUBLE PRECISION C_WMO_WS	    	! User value for the WS particles of the WMO model
      DOUBLE PRECISION C_WMO_OC	    	! User value for the OC particles of the WMO model
      DOUBLE PRECISION C_WMO_SO	    	! User value for the SO particles of the WMO model
      DOUBLE PRECISION RH                 ! % of relative humidity of the Shettle & Fenn models

	!Parameters for a bimodal LND model
      DOUBLE PRECISION USER_CV_COARSE	!User value of the volume concentration of the "LND coarse"
      DOUBLE PRECISION USER_CV_FINE	!User value of the volume concentration of the "LND fine"
      DOUBLE PRECISION RTAUct_WAREF       ! Ratio of the optical thicknesses AOT_coarse / AOT_total
      					       ! at the reference wavelength.
      DOUBLE PRECISION BMD_CM_MRWA	       ! Value read for RN_LND(1) : coarse mode
      DOUBLE PRECISION BMD_CM_MIWA	       ! Value read for IN_LND(1) : coarse mode
      DOUBLE PRECISION BMD_CM_MRWAREF	! Value read for RN_LND_WAREF(1) : coarse mode
      DOUBLE PRECISION BMD_CM_MIWAREF	! Value read for RN_LND_WAREF(2) : fine mode
      DOUBLE PRECISION BMD_CM_RMODAL	! Value read for RMODAL_LND(1) : coarse mode
      DOUBLE PRECISION BMD_CM_VAR	       ! Value read for VAR_LND(1) : coarse mode
      DOUBLE PRECISION BMD_FM_MRWA	       ! Value read for RN_LND(2) : fine mode
      DOUBLE PRECISION BMD_FM_MIWA	       ! Value read for IN_LND(2) : fine mode .
      DOUBLE PRECISION BMD_FM_MRWAREF	! Value read for IN_LND_WAREF(1) : coarse mode
      DOUBLE PRECISION BMD_FM_MIWAREF	! Value read for IN_LND_WAREF(2) : fine mode
      DOUBLE PRECISION BMD_FM_RMODAL	! Value read for RMODAL_LND(2) : fine mode
      DOUBLE PRECISION BMD_FM_VAR	       ! Value read for VAR_LND(2) : fine mode
 			     
      DOUBLE PRECISION COEF_TRONCA	! Truncation coefficient of the aerosol phase function
                                  
      INTEGER*4 ITRONC_AER   	! Truncation option of the aerosol phase function
					! (1 to apply it).
      INTEGER*2 IMOD_AER   		! Mono- or multi-modal size distribution.
      INTEGER*4 IGRANU	   	! Number of the size distribution model
      INTEGER*2 IModele_WMO 	! Type of the WMO model.
      INTEGER*2 IModele_SF  	! Type of the Shettle & Fenn model.
      INTEGER*2 MODE_PARAM_BILND   ! Index of the choice of description of the mixture of 2 LND.
                      
      CHARACTER*CTE_LENFIC2 FICMIE_LOG      ! name of the MIE TRACE file
      CHARACTER*CTE_LENFIC2 FICGRANU_LOG    ! name of the GRANU TRACE file
      CHARACTER*CTE_LENFIC2 FICGRANU        ! name of the GRANU result file
      CHARACTER*CTE_LENFIC2 FICEXTDATA_AER  ! path and name of the user file giving the
      						  ! phase function and the scattering
						  ! and extinction coefficients.
      CHARACTER*CTE_LENFIC2 FICUSER_AER     ! path and name of the aerosol file named by the user
      CHARACTER*CTE_LENFIC2 FICMIXTURE_AER  ! path and name of the file defining a mixture of aerosol models by the user
				
C --- Variables related to the definition of the atmospheric profile
  
      DOUBLE PRECISION TR     ! Rayleigh optical thickness.
      DOUBLE PRECISION HR     ! Scale height of the molecular profile (km)
      DOUBLE PRECISION HA     ! Scale height of the aerosol profile (km)
      DOUBLE PRECISION ZMIN   ! Minimum altitude of the molecules + aerosols layer.
      DOUBLE PRECISION ZMAX   ! Maximum altitude of the molecules + aerosols layer.
      DOUBLE PRECISION PSURF  ! Atmospheric pressure at the surface (mbar).
      DOUBLE PRECISION H2O    ! Integrated water vapour content (g/cm2).
      DOUBLE PRECISION O3     ! Integrated ozone content (Dobson).
      DOUBLE PRECISION CO2    ! Carbon dioxide concentration at the surface (ppmv).
      DOUBLE PRECISION CH4    ! Methane concentration at the surface (ppmv).

      DOUBLE PRECISION TDIR_VRAI   !Direct transmission for the solar incidence direction
      				   	!of the true atmosphere

      INTEGER*2 IPROFIL	  ! Type of the aerosol profile
      INTEGER*2 ABSPROFIL     ! Type of the absorption profile
        
      CHARACTER*CTE_LENFIC2 FICABSPROFIL    ! Absorption profile file provided by the user
      CHARACTER*CTE_LENFIC2 FICPROFIL_LOG   ! Path and name of the trace file of the profiles per CKD order
      
 
C --- Variables related to the choice of the CKD coefficient data
        
      DOUBLE PRECISION  NUSTEP  ! Wavenumber step (cm-1) in the CKD coeff file
                
C --- Variables related to the definition of the surface model
  
      INTEGER*2 ISURF	       ! Surface type

      DOUBLE PRECISION RHO	       ! Surface albedo.
      DOUBLE PRECISION WIND	       ! Wind speed (m/s).
      DOUBLE PRECISION SURF_IND	! Refractive index of the interface.
      DOUBLE PRECISION K0_ROUJEAN, K1_ROUJEAN, K2_ROUJEAN 	! Parameters of the Roujean model
      DOUBLE PRECISION ALPHA_NADAL	! Parameter alpha of the Nadal model.
      DOUBLE PRECISION BETA_NADAL	! Parameter beta of the Nadal model.
      DOUBLE PRECISION COEF_C_MAIGNAN  ! Parameter C*exp(-NDVI) of the Maignan model.
       
      CHARACTER*CTE_LENFIC2 FICSURF_LOG   ! name of the SURFACE trace file.
      CHARACTER*CTE_LENFIC2 FICSURF       ! name of a file containing an existing reflection matrix
      

C --- Variables related to the SOS computations

      DOUBLE PRECISION ZOUT	! Output type:
				!      -1 =>  Standard output: Radiance emerging at TOA
				!           + downward radiance at ground level.
				!       z =>  Upward and downward radiance
				!             for altitude z (0 <= z <= CTE_TOA_ALT).
			
      INTEGER*4 IGMAX	! Maximum order of scattering.
      
      INTEGER*4 IPOLAR	! Option for simulation without polarisation.
      			       ! 0 : computations without polarisation.
			       ! 1 : normal computations with polarisation.

      INTEGER*4 IMODE_CKD_CALCUL ! Option of the mode of application of the CKD decomposition
                                 ! 1 : Finest computation with one radiance computation per order of the decomposition
                                 !     then recomposition with the CKD weighting.
                                 ! 2 : Computation for an absorption optical thickness estimated by the CKD approach
                                 !     --> A single radiance computation

      INTEGER*4 ITRPHI	! Index of output type.
			       !   1 =>  Viewing plane for one azimuth value.
			       !   2 =>  Polar diagram versus the azimuth
			       !         and viewing directions.
			
      INTEGER*4 PAS_PHI	! Azimuth step in degrees (for ITRPHI=2)
      			
      DOUBLE PRECISION PHIOS ! Relative azimuth in degrees (for ITRPHI=1).
       
      CHARACTER*CTE_LENFIC2 FICSOS_LOG ! path and name of the trace file
      CHARACTER*CTE_LENFIC2 FICUP	    ! Name of the result file for the upward field.
      CHARACTER*CTE_LENFIC2 FICDOWN    ! Name of the result file for the downward field.
      CHARACTER*CTE_LENFIC2 FICUP_ANGLES_USER ! Name of the result file for the upward
      				                  ! field (for the user angles).
      CHARACTER*CTE_LENFIC2 FICDOWN_ANGLES_USER	! Name of the result file for the downward
      				                     ! field (for the user angles).
      CHARACTER*CTE_LENFIC2 FICTRANS    ! Name of the optional output file of the transmissions.
      CHARACTER*CTE_LENFIC2 FICFLUX     ! Name of the optional output file of the fluxes up TOA and down BOA.

      CHARACTER*CTE_LENFIC2 FICSOS_RES_BIN  ! Name of the OS result file (for a unit processing)

   
C --- Internal variables
      
C -   For the OS

      INTEGER*4 LUM_NBMU     ! Number of actual angles of the simulations

      INTEGER*4 IND_ANGOUT_FIN(0:CTE_OS_NBMU_MAX) 
                        ! Indices of the angle arrays for which there must be
      			   ! an output in the result file of the user angles
			   ! (1: output required)
                        				 

      DOUBLE PRECISION PHI_FIN(0:360)		   !Table of the relative azimuths of the output field
      DOUBLE PRECISION THETA_FIN(0:CTE_OS_NBMU_MAX) !Table of the viewing zenith angles of the output field
            
      DOUBLE PRECISION SCA_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)   !Table of the scattering angles for the output upward field
      DOUBLE PRECISION SCA_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX) !Table of the scattering angles for the output downward field

      DOUBLE PRECISION I_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)     !Table of the I values for the output upward field
      DOUBLE PRECISION Q_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)     !Table of the Q values for the output upward field
      DOUBLE PRECISION U_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)     !Table of the U values for the output upward field
      DOUBLE PRECISION I_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX)   !Table of the I values for the output downward field
      DOUBLE PRECISION Q_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX)   !Table of the Q values for the output downward field
      DOUBLE PRECISION U_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX)   !Table of the U values for the output downward field
      
      DOUBLE PRECISION POL_ANG_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)   !Table of the polarisation angles for the output upward field
      DOUBLE PRECISION POL_ANG_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX) !Table of the polarisation angles for the output downward field
     
      DOUBLE PRECISION POL_RATE_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)   !Table of the degrees of polarisation for the output upward field
      DOUBLE PRECISION POL_RATE_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX) !Table of the degrees of polarisation for the output downward field
  
      DOUBLE PRECISION L_POL_UP_FIN(0:360,0:CTE_OS_NBMU_MAX)   !Table of the polarised radiance values for the output upward field
      DOUBLE PRECISION L_POL_DOWN_FIN(0:360,0:CTE_OS_NBMU_MAX) !Table of the polarised radiance values for the output downward field
         
      DOUBLE PRECISION FLUX_DIFF_DOWN	!Downward diffuse flux at the surface level
      DOUBLE PRECISION FLUX_DOWN		!Total downward flux at the surface level

      DOUBLE PRECISION EPLUS    ! Diffuse flux emerging at TOA, normalised by the incident solar flux at TOA.
    
C       
      INTEGER*4 IPHI        ! Loop index on the relative azimuth.
   
      DOUBLE PRECISION XI   !  Value of the Stokes parameter I for an angle
      DOUBLE PRECISION XQ   ! Value of the Stokes parameter Q for an angle
      DOUBLE PRECISION XU   ! Value of the Stokes parameter U for an angle

      DOUBLE PRECISION XAN  ! Polarisation angle (degrees)
      DOUBLE PRECISION LPOL ! Polarised intensity normalised to the solar irradiance at TOA
      DOUBLE PRECISION TPOL ! Degree of polarisation in %
      DOUBLE PRECISION PHI1, PHI2 ! Relative azimuth (in radians).
   
   
C =========================================
C --- Initialisations 
C =========================================

C --- Initialisation of the processing error flag
C----------------------------------------------------------
      IER=0
	 
C --- Initialisation of the optional parameters with their default value
C-----------------------------------------------------------------------------

      FICMAIN_LOG="NO_LOG_FILE"

C --- Initialisation of the default values of the angle parameters

C     Number of Gauss angles for the radiance computations: not defined
      NBMU_GAUSS_LUM=CTE_NOT_DEFINED_VALUE_INT
      
C     Number of Gauss angles for the Mie computations: not defined
      NBMU_GAUSS_MIE=CTE_NOT_DEFINED_VALUE_INT
      
C     No "user" angle file for the radiance computations: not defined
      FICANGLES_USER_LUM="NO_USER_ANGLES"
      
C     No "user" angle file for the Mie computations: not defined
      FICANGLES_USER_MIE="NO_USER_ANGLES"
 
C     No trace file of the angle definition
      FICANGLOG = "NO_LOG_FILE"
       
C     Result files
      FICANGLES_RES_LUM = CTE_DEFAULT_FICANGLES_RES_LUM
      FICANGLES_RES_MIE = CTE_DEFAULT_FICANGLES_RES_MIE

C --- Initialisation of the mandatory parameters
      WA_SIMU=CTE_NOT_DEFINED_VALUE_DBLE   
      TETAS=CTE_NOT_DEFINED_VALUE_DBLE   
      AOT_REF=CTE_NOT_DEFINED_VALUE_DBLE

      RESROOT="UNDEFINED_REPERTORY"
      DIR_SURF="UNDEFINED_REPERTORY"
      DIR_MIE="UNDEFINED_REPERTORY"
      DIR_SOS="UNDEFINED_REPERTORY"
      DIR_LOG="UNDEFINED_REPERTORY"

      IMOD_AER=CTE_NOT_DEFINED_VALUE_INT
      WAREF_AOT=CTE_NOT_DEFINED_VALUE_DBLE
      RHO=CTE_NOT_DEFINED_VALUE_DBLE
      ISURF=CTE_NOT_DEFINED_VALUE_INT
      IPROFIL=CTE_NOT_DEFINED_VALUE_INT



C --- Initialisation of the default values of the aerosol parameters

c     Parameters for the mono-modal model
      RN_WA=CTE_NOT_DEFINED_VALUE_DBLE
      IN_WA=CTE_NOT_DEFINED_VALUE_DBLE
      RN_WAREF=CTE_NOT_DEFINED_VALUE_DBLE
      IN_WAREF=CTE_NOT_DEFINED_VALUE_DBLE 
      IGRANU=CTE_NOT_DEFINED_VALUE_INT
      LND_RADIUS_MMD_AER=CTE_NOT_DEFINED_VALUE_DBLE	
      LND_LNVAR_MMD_AER=CTE_NOT_DEFINED_VALUE_DBLE
      JD_SLOPE_MMD_AER=CTE_NOT_DEFINED_VALUE_DBLE
      JD_RMIN_MMD_AER=CTE_NOT_DEFINED_VALUE_DBLE
      JD_RMAX_MMD_AER=CTE_DEFAULT_AER_JUNGE_RMAX		!Default Rmax

C     WMO model parameters
      C_WMO_DL=CTE_NOT_DEFINED_VALUE_DBLE	
      C_WMO_WS=CTE_NOT_DEFINED_VALUE_DBLE	
      C_WMO_OC=CTE_NOT_DEFINED_VALUE_DBLE     
      C_WMO_SO=CTE_NOT_DEFINED_VALUE_DBLE	
      IModele_WMO=CTE_NOT_DEFINED_VALUE_INT

C     S&F parameters
      IModele_SF=CTE_NOT_DEFINED_VALUE_INT          
      RH=CTE_NOT_DEFINED_VALUE_DBLE

C     Bimodal LND parameters
      MODE_PARAM_BILND = CTE_NOT_DEFINED_VALUE_INT 
      USER_CV_COARSE = CTE_NOT_DEFINED_VALUE_DBLE 
      USER_CV_FINE 	= CTE_NOT_DEFINED_VALUE_DBLE  
      RTAUct_WAREF 	= CTE_NOT_DEFINED_VALUE_DBLE  
      BMD_CM_MRWA 	= CTE_NOT_DEFINED_VALUE_DBLE 
      BMD_CM_MIWA 	= CTE_NOT_DEFINED_VALUE_DBLE 
      BMD_CM_MRWAREF = CTE_NOT_DEFINED_VALUE_DBLE  
      BMD_CM_MIWAREF = CTE_NOT_DEFINED_VALUE_DBLE 
      BMD_CM_RMODAL 	= CTE_NOT_DEFINED_VALUE_DBLE  
      BMD_CM_VAR 	= CTE_NOT_DEFINED_VALUE_DBLE 
      BMD_FM_MRWA	= CTE_NOT_DEFINED_VALUE_DBLE  
      BMD_FM_MIWA	= CTE_NOT_DEFINED_VALUE_DBLE
      BMD_FM_MRWAREF = CTE_NOT_DEFINED_VALUE_DBLE 
      BMD_FM_MIWAREF = CTE_NOT_DEFINED_VALUE_DBLE
      BMD_FM_RMODAL	= CTE_NOT_DEFINED_VALUE_DBLE 
      BMD_FM_VAR 	= CTE_NOT_DEFINED_VALUE_DBLE 	 

C     User data parameters
      FICEXTDATA_AER = "NO_USER_AEROSOLS_PHAZE_FCT" 

C     Parameters of the user mixture of models
      FICMIXTURE_AER = "NO_USER_AEROSOLS_MIXTURE" 
      
C     Name of the user aerosol file
      FICUSER_AER="NO_USER_AEROSOLS"	 
              
C     Default name of the file of the aerosol radiative properties
      FICGRANU = CTE_DEFAULT_FICGRANU
     
C     No trace file of the Mie computations for the aerosols
      FICMIE_LOG = "NO_LOG_FILE"

C     No trace file of the definition of the aerosol phase matrices
      FICGRANU_LOG = "NO_LOG_FILE" 

C     Truncation of the forward peak of the aerosol phase function
      ITRONC_AER = 1
 
C --- Initialisation of the default values of the surface parameters
      
C     No trace file of the definition of the surface matrices
      FICSURF_LOG = "NO_LOG_FILE"

C	  Refractive index
      SURF_IND=CTE_NOT_DEFINED_VALUE_DBLE
      
C     Glitter parameters
      WIND=CTE_NOT_DEFINED_VALUE_DBLE
      
C     Roujean parameters
      K0_ROUJEAN=CTE_NOT_DEFINED_VALUE_DBLE
      K1_ROUJEAN=CTE_NOT_DEFINED_VALUE_DBLE
      K2_ROUJEAN=CTE_NOT_DEFINED_VALUE_DBLE 
      
C     Nadal parameters
      ALPHA_NADAL=CTE_NOT_DEFINED_VALUE_DBLE
      BETA_NADAL=CTE_NOT_DEFINED_VALUE_DBLE   

C     Maignan parameters
      COEF_C_MAIGNAN=CTE_NOT_DEFINED_VALUE_DBLE     

C     BRDF/BPDF result file
      FICSURF="DEFAULT"
      		
C --- Initialisation of the default values of the profile parameters
     
C     Molecular optical thickness
      TR=CTE_NOT_DEFINED_VALUE_DBLE
           
C     Scale heights
      HA=CTE_NOT_DEFINED_VALUE_DBLE
      HR=CTE_NOT_DEFINED_VALUE_DBLE
      
C     Pressure (default value)
      PSURF=CTE_HT_STD_PSURF
 
C     Min and max altitudes of the aerosol profile
      ZMIN=CTE_NOT_DEFINED_VALUE_DBLE
      ZMAX=CTE_NOT_DEFINED_VALUE_DBLE     
 
C     Absorption profile
      ABSPROFIL=CTE_NOT_DEFINED_VALUE_INT

C     Contents of absorbing gases
      H2O=CTE_NOT_DEFINED_VALUE_DBLE
      O3=CTE_NOT_DEFINED_VALUE_DBLE
      CO2=CTE_NOT_DEFINED_VALUE_DBLE
      CH4=CTE_NOT_DEFINED_VALUE_DBLE

C     File of the absorbing gas profile
      FICABSPROFIL = "NO_USER_ABS_PROFILE_FILE"             

C     CKD data
      NUSTEP=CTE_NOT_DEFINED_VALUE_DBLE
      IMODE_CKD_CALCUL=CTE_NOT_DEFINED_VALUE_INT
  
C     No trace file
      FICPROFIL_LOG = "NO_LOG_FILE"

 
            
C --- Initialisation of the default values of the SOS parameters
      
C     No trace file
      FICSOS_LOG = "NO_LOG_FILE"

C     Option to cut the polarisation: optional (computation with polarisation by default)
      IPOLAR=1   

C     Maximum number of scatterings
      IGMAX=CTE_DEFAULT_IGMAX

      
C     Output level of the radiance simulations: optional
      ZOUT=-1.0

C     Output files on user angles
      FICUP_ANGLES_USER="NO_OUTPUT"
      FICDOWN_ANGLES_USER="NO_OUTPUT"            

C     Output file of the transmissions
      FICTRANS="NO_OUTPUT"

C     Output file of the fluxes
      FICFLUX="NO_OUTPUT"

C     Result files
      FICSOS_RES_BIN=CTE_DEFAULT_RESBIN
      FICUP=CTE_DEFAULT_RESUP
      FICDOWN=CTE_DEFAULT_RESDOWN
   
C     Definition of the azimuths
      ITRPHI=CTE_NOT_DEFINED_VALUE_INT  
      PHIOS=CTE_NOT_DEFINED_VALUE_DBLE
      PAS_PHI=CTE_NOT_DEFINED_VALUE_INT    
  
             		       
C =========================================================
C --- Reading of the input parameters/arguments and decoding
C =========================================================           				
      
C --- Number of input arguments
      NBARGS = IARGC()
      
C --- Reading loop of the pairs of arguments
      DO NPARAM=1,NBARGS, 2
         CALL GETARG(NPARAM,KEYWORD)
         KEYWORD_VALID=0
		
C ------- Elementary validity check of the keyword
         IF (KEYWORD(1:1).NE.'-') GOTO 100

C ------- Arguments specific to the main program

	 IF (KEYWORD.EQ."-SOS_Main.ResRoot") THEN
	     CALL GETARG (NPARAM+1,RESROOT)		!Name of the working directory.
            KEYWORD_VALID=1
	 ENDIF	

	 IF (KEYWORD.EQ."-SOS_Main.Log") THEN
	     CALL GETARG (NPARAM+1,FICMAIN_LOG)		!Name of the global trace file.
            KEYWORD_VALID=1
	 ENDIF		 
		
	 IF (KEYWORD.EQ."-SOS_Main.Wa") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)     !Wavelength of the radiance simulation
	     READ(PROG_ARG,*,ERR=101) WA_SIMU    
            KEYWORD_VALID=1
	 ENDIF	
	 
C ------- Arguments specific to the definition of the angles
	
	 IF (KEYWORD.EQ."-ANG.Rad.NbGauss") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)	!Nb of Gauss angles for the radiance simulations
	     READ(PROG_ARG,*,ERR=201) NBMU_GAUSS_LUM
            KEYWORD_VALID=1
	 ENDIF	 
	 	 
	 IF (KEYWORD.EQ."-ANG.Rad.UserAngFile") THEN
	     CALL GETARG (NPARAM+1,FICANGLES_USER_LUM)	!Name of the "user" angles file (complete path)
            KEYWORD_VALID=1			              !for the radiance computations
	 ENDIF	  	 
	 		 
	 IF (KEYWORD.EQ."-ANG.Thetas") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Solar zenith angle (degrees)
	     READ(PROG_ARG,*,ERR=202) TETAS   
            KEYWORD_VALID=1
	 ENDIF	

	 IF (KEYWORD.EQ."-ANG.Rad.ResFile") THEN
	     CALL GETARG (NPARAM+1,FICANGLES_RES_LUM)	!Name of the result file (complete path)
            KEYWORD_VALID=1	     				!for the radiance computations
	 ENDIF	  	
				
	 IF (KEYWORD.EQ."-ANG.Aer.NbGauss") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)	!Nb of Gauss angles for the Mie computations
	     READ(PROG_ARG,*,ERR=203) NBMU_GAUSS_MIE
            KEYWORD_VALID=1
	 ENDIF	 
	 	 
	 IF (KEYWORD.EQ."-ANG.Aer.UserAngFile") THEN
	     CALL GETARG (NPARAM+1,FICANGLES_USER_MIE)	!Name of the "user" angles file (complete path)
            KEYWORD_VALID=1		        	       !for the Mie computations
	 ENDIF	  	 
	 	 	 
	 IF (KEYWORD.EQ."-ANG.Aer.ResFile") THEN
	     CALL GETARG (NPARAM+1,FICANGLES_RES_MIE)	!Name of the result file (complete path)
            KEYWORD_VALID=1	     				!for the Mie computations
	 ENDIF	  	
	 	 
	 IF (KEYWORD.EQ."-ANG.Log") THEN
	     CALL GETARG (NPARAM+1,FICANGLOG)		!Name of the trace file.
            KEYWORD_VALID=1
	 ENDIF		

C ------- Arguments specific to the definition of the aerosol models

	 IF (KEYWORD.EQ."-AER.Waref") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)     	!Reference wavelength
	     READ(PROG_ARG,*,ERR=301) WAREF_AOT 	!for the aerosol opt. thickness.
            KEYWORD_VALID=1
	 ENDIF

	 IF (KEYWORD.EQ."-AER.AOTref") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)   !Total AOT at the wavelength WAREF
	     READ(PROG_ARG,*,ERR=302) AOT_REF
            KEYWORD_VALID=1
	 ENDIF 	 

	 IF (KEYWORD.EQ."-AER.Tronca") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)       	!Truncation index
	     READ(PROG_ARG,*,ERR=303) ITRONC_AER
            KEYWORD_VALID=1
	 ENDIF	 
	 
	 IF (KEYWORD.EQ."-AER.Log") THEN
	     CALL GETARG (NPARAM+1,FICGRANU_LOG) !Trace file of the Aerosols program
            KEYWORD_VALID=1
	 ENDIF	 	  	 

	 IF (KEYWORD.EQ."-AER.MieLog") THEN
	     CALL GETARG (NPARAM+1,FICMIE_LOG)	!Trace file of the Mie computations
            KEYWORD_VALID=1
	 ENDIF
         
         IF (KEYWORD.EQ."-AER.DirMie") THEN
	     CALL GETARG (NPARAM+1,DIR_MIE)	!Storage directory of the Mie files
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AER.ResFile") THEN
	     CALL GETARG (NPARAM+1,FICGRANU)	!Result file
            KEYWORD_VALID=1
	 ENDIF	 	  	
	 	 	  
	 IF (KEYWORD.EQ."-AER.Model") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)	!Type of aerosol model
	     READ(PROG_ARG,*,ERR=304) IMOD_AER
            KEYWORD_VALID=1
	 ENDIF	 				   	 	  	  

	 IF (KEYWORD.EQ."-AER.MMD.MRwa") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!real part of the refractive index for a mono-modal model
	     READ(PROG_ARG,*,ERR=306) RN_WA
            KEYWORD_VALID=1
	 ENDIF	 

	 IF (KEYWORD.EQ."-AER.MMD.MIwa") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!imaginary part of the refractive index for a mono-modal model
	     READ(PROG_ARG,*,ERR=307) IN_WA
            KEYWORD_VALID=1
	 ENDIF

	 IF (KEYWORD.EQ."-AER.MMD.MRwaref") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Real part of the refractive index for a mono-modal
	     READ(PROG_ARG,*,ERR=308) RN_WAREF		!aerosol model at the reference wavelength
            KEYWORD_VALID=1                      	
	 ENDIF	 					

	 IF (KEYWORD.EQ."-AER.MMD.MIwaref") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Imaginary part of the refractive index for a mono-modal
	     READ(PROG_ARG,*,ERR=309) IN_WAREF	 	!aerosol model at the reference wavelength
            KEYWORD_VALID=1                      	
	 ENDIF	 					
	 
	 IF (KEYWORD.EQ."-AER.MMD.SDtype") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!index of the size distribution type for a mono-modal model
	     READ(PROG_ARG,*,ERR=310) IGRANU
            KEYWORD_VALID=1
	 ENDIF

	 IF (KEYWORD.EQ."-AER.MMD.LNDradius") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Modal radius (microns) of the LND size distribution for a mono-modal aerosol
	     READ(PROG_ARG,*,ERR=311) LND_RADIUS_MMD_AER	
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AER.MMD.LNDvar") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Variance of the LND size distribution for a mono-modal aerosol
	     READ(PROG_ARG,*,ERR=312) LND_LNVAR_MMD_AER	
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AER.MMD.JD.slope") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Slope of the Junge law for a mono-modal aerosol model
	     READ(PROG_ARG,*,ERR=3121) JD_SLOPE_MMD_AER	
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AER.MMD.JD.rmin") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Minimum radius of the Junge law (microns) for a mono-modal aerosol model
	     READ(PROG_ARG,*,ERR=3122) JD_RMIN_MMD_AER	
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AER.MMD.JD.rmax") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Maximum radius of the Junge law (microns) for a mono-modal aerosol model
	     READ(PROG_ARG,*,ERR=3123) JD_RMAX_MMD_AER	
            KEYWORD_VALID=1
	 ENDIF							
	 
	 IF (KEYWORD.EQ."-AER.WMO.Model") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)		!Index of the choice of WMO model
	     READ(PROG_ARG,*,ERR=313) IModele_WMO 
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AER.WMO.DL") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)	!Volume percentage of the DL particles of the WMO model
	     READ(PROG_ARG,*,ERR=314) C_WMO_DL
            KEYWORD_VALID=1
	 ENDIF

	 IF (KEYWORD.EQ."-AER.WMO.WS") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)    	!Volume percentage of the WS particles of the WMO model
	     READ(PROG_ARG,*,ERR=314) C_WMO_WS
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AER.WMO.OC") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)   	!Volume percentage of the OC particles
	     READ(PROG_ARG,*,ERR=314) C_WMO_OC
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AER.WMO.SO") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)   	!Volume percentage of the SO particles
	     READ(PROG_ARG,*,ERR=314) C_WMO_SO
            KEYWORD_VALID=1
	 ENDIF	
	 
	 IF (KEYWORD.EQ."-AER.SF.Model") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)	!Index of the choice of S&F model
	     READ(PROG_ARG,*,ERR=315) IModele_SF  
            KEYWORD_VALID=1
	 ENDIF 	
	 
	 IF (KEYWORD.EQ."-AER.SF.RH") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)	!Relative humidity for the S&F model
	     READ(PROG_ARG,*,ERR=316) RH 
            KEYWORD_VALID=1
	 ENDIF 	 
	 	 
	 IF (KEYWORD.EQ."-AER.BMD.VCdef") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)   !Index of the type of description of the bimodal LND mixture
	     READ(PROG_ARG,*,ERR=317) MODE_PARAM_BILND 
            KEYWORD_VALID=1
	 ENDIF 	 	 

	 IF (KEYWORD.EQ."-AER.BMD.CoarseVC") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)   !Volume concentration of the "coarse" mode
	     READ(PROG_ARG,*,ERR=318) USER_CV_COARSE
            KEYWORD_VALID=1
	 ENDIF 	 

	 IF (KEYWORD.EQ."-AER.BMD.FineVC") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)   !Volume concentration of the "fine" mode
	     READ(PROG_ARG,*,ERR=319) USER_CV_FINE 
            KEYWORD_VALID=1
	 ENDIF 	 	 

	 IF (KEYWORD.EQ."-AER.BMD.RAOT") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)	 !AOT_coarse / AOT_tot at the wavelength WAREF
	     READ(PROG_ARG,*,ERR=320) RTAUct_WAREF 
            KEYWORD_VALID=1
	 ENDIF 	 
	 

	 IF (KEYWORD.EQ."-AER.BMD.CM.MRwa") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=321) BMD_CM_MRWA
            KEYWORD_VALID=1
	 ENDIF 
	 
	 IF (KEYWORD.EQ."-AER.BMD.CM.MIwa") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=322) BMD_CM_MIWA
            KEYWORD_VALID=1
	 ENDIF	 

	 IF (KEYWORD.EQ."-AER.BMD.CM.MRwaref") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=323) BMD_CM_MRWAREF
            KEYWORD_VALID=1
	 ENDIF 
	 
	 IF (KEYWORD.EQ."-AER.BMD.CM.MIwaref") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=324) BMD_CM_MIWAREF
            KEYWORD_VALID=1
	 ENDIF	 

	 IF (KEYWORD.EQ."-AER.BMD.CM.SDradius") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=325) BMD_CM_RMODAL
            KEYWORD_VALID=1
	 ENDIF 

	 IF (KEYWORD.EQ."-AER.BMD.CM.SDvar") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=326) BMD_CM_VAR
            KEYWORD_VALID=1
	 ENDIF 
	 	 
	 IF (KEYWORD.EQ."-AER.BMD.FM.MRwa") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=327) BMD_FM_MRWA
            KEYWORD_VALID=1
	 ENDIF 
	 
	 IF (KEYWORD.EQ."-AER.BMD.FM.MIwa") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=328) BMD_FM_MIWA	
            KEYWORD_VALID=1
	 ENDIF	 

	 IF (KEYWORD.EQ."-AER.BMD.FM.MRwaref") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=329) BMD_FM_MRWAREF
            KEYWORD_VALID=1
	 ENDIF 
	 
	 IF (KEYWORD.EQ."-AER.BMD.FM.MIwaref") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=330) BMD_FM_MIWAREF
            KEYWORD_VALID=1
	 ENDIF	 

	 IF (KEYWORD.EQ."-AER.BMD.FM.SDradius") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=331) BMD_FM_RMODAL  
            KEYWORD_VALID=1
	 ENDIF 

	 IF (KEYWORD.EQ."-AER.BMD.FM.SDvar") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) 
	     READ(PROG_ARG,*,ERR=332) BMD_FM_VAR	 
            KEYWORD_VALID=1
	 ENDIF 	
	 
	 IF (KEYWORD.EQ."-AER.ExtData") THEN
	     CALL GETARG (NPARAM+1,FICEXTDATA_AER) 
            KEYWORD_VALID=1
	 ENDIF 	 	 	 	 	 	 		  	 	 	 	 	 	 	 		  
	 
	 IF (KEYWORD.EQ."-AER.DefMixture") THEN
	     CALL GETARG (NPARAM+1,FICMIXTURE_AER) 
            KEYWORD_VALID=1
	 ENDIF 	
	 	 	 
	 IF (KEYWORD.EQ."-AER.UserFile") THEN
	     CALL GETARG (NPARAM+1,FICUSER_AER) 
            KEYWORD_VALID=1
	 ENDIF
	 
C ------- Arguments specific to the definition of the atmospheric profile
	 IF (KEYWORD.EQ."-AP.Log") THEN
	     CALL GETARG (NPARAM+1,FICPROFIL_LOG) !Trace file of the program
            KEYWORD_VALID=1
	 ENDIF	 	  	 
	 	 
	 IF (KEYWORD.EQ."-AP.MOT") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Molecular optical thickness
	     READ(PROG_ARG,*,ERR=401) TR
            KEYWORD_VALID=1
	 ENDIF	
	 	 
	 IF (KEYWORD.EQ."-AP.HR") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Molecular scale height (km)
	     READ(PROG_ARG,*,ERR=402) HR
            KEYWORD_VALID=1
	 ENDIF	

	 IF (KEYWORD.EQ."-AP.AerHS.HA") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Aerosol scale height (km)
	     READ(PROG_ARG,*,ERR=403) HA
            KEYWORD_VALID=1
	 ENDIF	

	 IF (KEYWORD.EQ."-AP.AerProfile.Type") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Type of the simulated aerosol profile
	     READ(PROG_ARG,*,ERR=404) IPROFIL 
            KEYWORD_VALID=1
	 ENDIF	
	 	 	 	 
	 IF (KEYWORD.EQ."-AP.AerLayer.Zmin") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Minimum altitude of the aerosol layer (km)
	     READ(PROG_ARG,*,ERR=405) ZMIN
            KEYWORD_VALID=1
	 ENDIF		 
	 	 
	 IF (KEYWORD.EQ."-AP.AerLayer.Zmax") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Maximum altitude of the aerosol layer (km)
	     READ(PROG_ARG,*,ERR=406) ZMAX
            KEYWORD_VALID=1
	 ENDIF	
	 
	 IF (KEYWORD.EQ."-AP.Psurf") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Surface pressure
	     READ(PROG_ARG,*,ERR=407) PSURF
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AP.H2O") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Water vapour content
	     READ(PROG_ARG,*,ERR=408) H2O
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AP.O3") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Ozone content
	     READ(PROG_ARG,*,ERR=409) O3
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AP.CO2") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Carbon dioxide concentration at the surface
	     READ(PROG_ARG,*,ERR=4091) CO2
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AP.CH4") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Methane concentration at the surface
	     READ(PROG_ARG,*,ERR=4092) CH4
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AP.AbsProfile.Type") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Profile type for the absorption
	     READ(PROG_ARG,*,ERR=410) ABSPROFIL 
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AP.AbsProfile.UserFile") THEN
	     CALL GETARG (NPARAM+1,FICABSPROFIL)	!User absorption profile file
            KEYWORD_VALID=1
	 ENDIF
	 
	 IF (KEYWORD.EQ."-AP.SpectralResol") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG) !Spectral resolution of the CKD files to use (cm-1)
	     READ(PROG_ARG,*,ERR=411) NUSTEP 
            KEYWORD_VALID=1
	 ENDIF
     
     
C ------- Arguments specific to the surface model

	 IF (KEYWORD.EQ."-SURF.Type") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
	     READ(PROG_ARG,*,ERR=501) ISURF
            KEYWORD_VALID=1
	 ENDIF

	 IF (KEYWORD.EQ."-SURF.Dir") THEN
	     CALL GETARG (NPARAM+1,DIR_SURF)
            KEYWORD_VALID=1
	 ENDIF	  	

	 IF (KEYWORD.EQ."-SURF.Log") THEN
	     CALL GETARG (NPARAM+1,FICSURF_LOG)
            KEYWORD_VALID=1
	 ENDIF	  	 

	 IF (KEYWORD.EQ."-SURF.Ind") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
	     READ(PROG_ARG,*,ERR=502) SURF_IND
            KEYWORD_VALID=1
	 ENDIF
	  
	 IF (KEYWORD.EQ."-SURF.Glitter.Wind") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
            READ(PROG_ARG,*,ERR=503) WIND    
            KEYWORD_VALID=1
	 ENDIF

	 IF (KEYWORD.EQ."-SURF.Roujean.K0") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
            READ(PROG_ARG,*,ERR=504) K0_ROUJEAN    
            KEYWORD_VALID=1
	 ENDIF
	  
	 IF (KEYWORD.EQ."-SURF.Roujean.K1") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
            READ(PROG_ARG,*,ERR=505) K1_ROUJEAN    
            KEYWORD_VALID=1
	 ENDIF

	 IF (KEYWORD.EQ."-SURF.Roujean.K2") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
             READ(PROG_ARG,*,ERR=506) K2_ROUJEAN    
            KEYWORD_VALID=1
	 ENDIF
	  	 
	 IF (KEYWORD.EQ."-SURF.Nadal.Alpha") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
            READ(PROG_ARG,*,ERR=507) ALPHA_NADAL    
            KEYWORD_VALID=1
	 ENDIF	 

        IF (KEYWORD.EQ."-SURF.Nadal.Beta") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
            READ(PROG_ARG,*,ERR=508) BETA_NADAL	
            KEYWORD_VALID=1
        ENDIF	  

        IF (KEYWORD.EQ."-SURF.Maignan.C") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)
            READ(PROG_ARG,*,ERR=5081) COEF_C_MAIGNAN	
            KEYWORD_VALID=1
        ENDIF	

	 IF (KEYWORD.EQ."-SURF.Alb") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)	!Surface albedo
	     READ(PROG_ARG,*,ERR=509) RHO
            KEYWORD_VALID=1
	 ENDIF 	 	 

	 IF (KEYWORD.EQ."-SURF.File") THEN
	     CALL GETARG (NPARAM+1,FICSURF)	!Name of the reflection matrix file
            KEYWORD_VALID=1
	 ENDIF
	 
C ------- Arguments specific to the Successive Orders computations and to the production
C         of the result files

	 IF (KEYWORD.EQ."-SOS.Log") THEN
	     CALL GETARG (NPARAM+1,FICSOS_LOG)	!Name of the trace file.
            KEYWORD_VALID=1
	 ENDIF	 

	 IF (KEYWORD.EQ."-SOS.ResBin") THEN
	     CALL GETARG (NPARAM+1,FICSOS_RES_BIN) !Name of the binary OS file.
            KEYWORD_VALID=1
	 ENDIF

	 IF (KEYWORD.EQ."-SOS.ResFileUp") THEN
	     CALL GETARG (NPARAM+1,FICUP)	!Name of the output file for the UP field.
            KEYWORD_VALID=1
	 ENDIF

	 IF (KEYWORD.EQ."-SOS.ResFileDown") THEN
	     CALL GETARG (NPARAM+1,FICDOWN) !Name of the output file for the DOWN field.
            KEYWORD_VALID=1
	 ENDIF

	 IF (KEYWORD.EQ."-SOS.ResFileUp.UserAng") THEN
            CALL GETARG (NPARAM+1,FICUP_ANGLES_USER) !Name of the output file for the UP field
            KEYWORD_VALID=1     !at the "user" angles.
	 ENDIF

	 IF (KEYWORD.EQ."-SOS.ResFileDown.UserAng") THEN
	     CALL GETARG (NPARAM+1,FICDOWN_ANGLES_USER)	!Name of the output file for the DOWN field.
            KEYWORD_VALID=1                             !at the "user" angles.
	 ENDIF

	 IF (KEYWORD.EQ."-SOS.Trans") THEN
	     CALL GETARG (NPARAM+1,FICTRANS)  !Name of the output file of the transmission values.
            KEYWORD_VALID=1
	 ENDIF

	 IF (KEYWORD.EQ."-SOS.Flux") THEN
	     CALL GETARG (NPARAM+1,FICFLUX)  !Name of the output file of the flux values
            KEYWORD_VALID=1
	 ENDIF
          	 	 		 
         IF (KEYWORD.EQ."-SOS.OutputAlt") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)  !Output altitude
            READ(PROG_ARG,*,ERR=601) ZOUT  
            KEYWORD_VALID=1
         ENDIF	 

         IF (KEYWORD.EQ."-SOS.IGmax") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)  !Maximum order of interaction
            READ(PROG_ARG,*,ERR=602) IGMAX
            KEYWORD_VALID=1
         ENDIF	 

         IF (KEYWORD.EQ."-SOS.Ipolar") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)  !Option to cut the polarisation
            READ(PROG_ARG,*,ERR=603) IPOLAR
            KEYWORD_VALID=1
         ENDIF	 

         IF (KEYWORD.EQ."-SOS.View") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)  !Viewing option
            READ(PROG_ARG,*,ERR=604) ITRPHI
            KEYWORD_VALID=1
         ENDIF	 

         IF (KEYWORD.EQ."-SOS.View.Phi") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)  !Relative azimuth angle
            READ(PROG_ARG,*,ERR=605) PHIOS 
            KEYWORD_VALID=1
         ENDIF	 

         IF (KEYWORD.EQ."-SOS.View.Dphi") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)  !Step in azimuth angles
            READ(PROG_ARG,*,ERR=606) PAS_PHI     
            KEYWORD_VALID=1
         ENDIF	 		 

         IF (KEYWORD.EQ."-SOS.AbsModeCKD") THEN
	     CALL GETARG (NPARAM+1,PROG_ARG)  !Option to cut the polarisation
            READ(PROG_ARG,*,ERR=608) IMODE_CKD_CALCUL
            KEYWORD_VALID=1
         ENDIF
	 

         IF (KEYWORD_VALID.EQ.0) THEN
	     WRITE(6,*) 
	     WRITE(6,*)	     
     &      "**************************************************"
	     WRITE(6,*) "WARNING:"
	     WRITE(6,*) "The defined parameter:",KEYWORD
	     WRITE(6,*) "is not part of the SOS-ABS parameters."
	     WRITE(6,*) 
     &      "==> Its definition does not act on the simulation."
	     WRITE(6,*) 
     &      "==> Be sure not to confuse the keyword parameters!"
	     WRITE(6,*)
     &      "**************************************************"
	     WRITE(6,*) 
         ENDIF
	 		
      ENDDO	!End of the reading loop of the arguments


C =========================================================
C --- Definition of the working directory tree
C =========================================================           				

C ---- Check of the definition of the working space (mandatory)
      IF (RESROOT.EQ."UNDEFINED_REPERTORY") GOTO 1000
      
      LDIR = INDEX(RESROOT,' ') - 1
      IF((LDIR.LE.0).OR.(LDIR.GT.CTE_LENDIR)) GOTO 1100
		
C ---- Definition of the SOS result directory tree
      LDIRSOS= LDIR + 4  	
      IF(LDIRSOS.GT.CTE_LENDIR) GOTO 1108
      DIR_SOS(1:LDIRSOS) = RESROOT(1:LDIR)//'/'//'SOS'	


C =========================================================
C --- Definition of the main trace file
C =========================================================        

       !Main trace file
       IF (FICMAIN_LOG.EQ."NO_LOG_FILE") THEN
	   TRACE=.FALSE.
       ELSE
          TRACE=.TRUE.	
       ENDIF ! End of Trace
 
 
      
C -------------------------------------------------------------------  
C --- Check of the parameters specific to the Successive Orders computations
C     and to the production of the result files
C -------------------------------------------------------------------  
C  
C ---- Existence and validity test of the mandatory parameters
      
C     Viewing option: mandatory
      IF (ITRPHI.EQ.CTE_NOT_DEFINED_VALUE_INT) THEN
          GOTO 2605
      ELSE
	  
          !If ITRPHI=1 is chosen: the relative azimuth must be defined
          IF (ITRPHI.EQ.1) THEN
	      IF (PHIOS.EQ.CTE_NOT_DEFINED_VALUE_DBLE) THEN
	         GOTO 2607
	      ELSE
	        !In the case of a definition of the observation azimuth plane,
	        !PAS_PHI does not need to be defined since it is not used by the SOS code,
	        !but the argument is expected as input of SOS_TRPHI_OPTION
	        ! ==> any (integer) value is defined
	        PAS_PHI=-9999
	     ENDIF
	  ENDIF
	  
	  !If ITRPHI=2 is chosen: the azimuth step must be defined
	  IF (ITRPHI.EQ.2) THEN
	     IF (PAS_PHI.EQ.CTE_NOT_DEFINED_VALUE_INT) THEN
	        GOTO 2608
	     ELSE
	        !In the case of a definition of the observation field with an azimuth step,
	        !PHIOS does not need to be defined since it is not used by the SOS code,
	        !but the argument is expected as input of SOS_TRPHI_OPTION
	        ! ==> any value is defined
	        PHIOS=-9999.0
	     ENDIF
	  ENDIF
      ENDIF ! ITRPHI
 

C =========================================================
C --- Management of the directory tree of the output files
C =========================================================     
                              
C ---- Output files
         	     	
      !Upward radiance file
      LFIC = INDEX(FICUP,' ') - 1
      IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3013
      FICUP = DIR_SOS(1:LDIRSOS)//'/'//FICUP(1:LFIC)

      IF(FICUP_ANGLES_USER.NE."NO_OUTPUT") THEN
	  LFIC = INDEX(FICUP_ANGLES_USER,' ') - 1
         IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3013
         FICUP_ANGLES_USER 
     &    = DIR_SOS(1:LDIRSOS)//'/'//FICUP_ANGLES_USER(1:LFIC)
      ENDIF

      !Downward radiance file
      LFIC = INDEX(FICDOWN,' ') - 1
      IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3014
      FICDOWN = DIR_SOS(1:LDIRSOS)//'/'//FICDOWN(1:LFIC)

      IF(FICDOWN_ANGLES_USER.NE."NO_OUTPUT") THEN
	  LFIC = INDEX(FICDOWN_ANGLES_USER,' ') - 1
         IF((LFIC.LE.0).OR.(LFIC.GT.CTE_LENFIC1)) GOTO 3013
         FICDOWN_ANGLES_USER 
     &    = DIR_SOS(1:LDIRSOS)//'/'//FICDOWN_ANGLES_USER(1:LFIC)
      ENDIF


C     Call of the subroutine SOS_PROC
      CALL SOS_PROC(RESROOT, FICMAIN_LOG, WA_SIMU, NBMU_GAUSS_LUM,
     &     FICANGLES_USER_LUM, TETAS, FICANGLES_RES_LUM,
     &     NBMU_GAUSS_MIE, FICANGLES_USER_MIE, FICANGLES_RES_MIE,
     &     FICANGLOG, WAREF_AOT, AOT_REF, ITRONC_AER, FICGRANU_LOG,
     &     FICMIE_LOG, DIR_MIE, FICGRANU, IMOD_AER, 
     &     RN_WA, IN_WA, RN_WAREF, IN_WAREF, 
     &     IGRANU, LND_RADIUS_MMD_AER, LND_LNVAR_MMD_AER,
     &     JD_SLOPE_MMD_AER, JD_RMIN_MMD_AER, JD_RMAX_MMD_AER,
     &     IModele_WMO, C_WMO_DL, C_WMO_WS, C_WMO_OC, C_WMO_SO,
     &     IModele_SF, RH, MODE_PARAM_BILND, USER_CV_COARSE,
     &     USER_CV_FINE, RTAUct_WAREF, BMD_CM_MRWA, BMD_CM_MIWA,
     &     BMD_CM_MRWAREF, BMD_CM_MIWAREF, BMD_CM_RMODAL, BMD_CM_VAR,
     &     BMD_FM_MRWA, BMD_FM_MIWA, BMD_FM_MRWAREF, BMD_FM_MIWAREF,
     &     BMD_FM_RMODAL, BMD_FM_VAR, FICEXTDATA_AER, FICMIXTURE_AER,
     &     FICUSER_AER, FICPROFIL_LOG, TR, HR, HA, IPROFIL, ZMIN, ZMAX,
     &     PSURF, H2O, O3, CO2, CH4, ABSPROFIL, FICABSPROFIL, NUSTEP,
     &     ISURF, DIR_SURF, FICSURF_LOG, SURF_IND, WIND, K0_ROUJEAN,
     &     K1_ROUJEAN, K2_ROUJEAN, ALPHA_NADAL, BETA_NADAL,
     &     COEF_C_MAIGNAN, RHO, FICSURF, FICSOS_LOG, FICSOS_RES_BIN,
     &     FICTRANS, FICFLUX, ZOUT, IGMAX, IPOLAR, ITRPHI, PHIOS,
     &     PAS_PHI, IMODE_CKD_CALCUL,
     &     IER, TRACE,
     &     LUM_NBMU,IND_ANGOUT_FIN, PHI_FIN, THETA_FIN, 
     &     SCA_UP_FIN, I_UP_FIN, Q_UP_FIN, U_UP_FIN,
     &     POL_ANG_UP_FIN, POL_RATE_UP_FIN, L_POL_UP_FIN,
     &     SCA_DOWN_FIN, I_DOWN_FIN, Q_DOWN_FIN, U_DOWN_FIN,
     &     POL_ANG_DOWN_FIN, POL_RATE_DOWN_FIN, L_POL_DOWN_FIN,
     &     TDIR_VRAI, FLUX_DIFF_DOWN, FLUX_DOWN, EPLUS,
     &     COEF_TRONCA)


C --- Management of the error case
      
      IF (IER.NE.0) THEN        ! error case
         WRITE (6,*)
     &        "   ==> Error in subroutine SOS_PROC"
         GOTO 9998
      ENDIF

C --- Trace file
      IF (TRACE) THEN
         OPEN (IDLOG,FILE=FICMAIN_LOG,STATUS='OLD',
     &         ACCESS='APPEND',ERR=1302)
         WRITE (IDLOG,*,ERR=1303) " "

    	  WRITE (IDLOG,*,ERR=1303) "FICUP_ANGLES_USER ",
     &	  FICUP_ANGLES_USER
    	  WRITE (IDLOG,*,ERR=1303) "FICDOWN_ANGLES_USER ",
     &	  FICDOWN_ANGLES_USER
 	  WRITE (IDLOG,*,ERR=1303) "FICUP ",FICUP
    	  WRITE (IDLOG,*,ERR=1303) "FICDOWN ",FICDOWN
      ENDIF

C* Opening of the files
C------------------------------------------------------
C     Result files for the upward field and the downward field
C     If the file already exists, it is overwritten by the new file
      OPEN(UNIT=14,FILE=FICUP,FORM='FORMATTED',ERR=991)
      
      OPEN(UNIT=15,FILE=FICDOWN,FORM='FORMATTED',ERR=991)

C     Result files for the upward field and the downward field
C     If the file already exists, it is overwritten by the new file
C     --> Case of the optional files containing the results for the
C         user angles only
      IF (FICUP_ANGLES_USER.NE.'NO_OUTPUT') THEN
          OPEN(UNIT=24,FILE=FICUP_ANGLES_USER,FORM='FORMATTED',ERR=991)
      ENDIF
      	  
      IF (FICDOWN_ANGLES_USER.NE.'NO_OUTPUT') THEN
          OPEN(UNIT=25,FILE=FICDOWN_ANGLES_USER,FORM='FORMATTED',
     &	       ERR=991)
      ENDIF

    
C* Output of one azimuth plane
C------------------------------------------------------
      IF(ITRPHI.EQ.1) THEN

         PHI1 = PHIOS
         PHI2 = PHIOS+180.D+00

         !Writing of the header information of the file
         !------------------------------------------------
         IF (ZOUT.EQ.-1) THEN   

            CALL SOS_OUTPUT_HEADER(14, 1, PHI1, PHI2, CTE_TOA_ALT) 
            CALL SOS_OUTPUT_HEADER(15, 2, PHI1, PHI2, 0.D+00)

            IF (FICUP_ANGLES_USER.NE.'NO_OUTPUT') THEN
               CALL SOS_OUTPUT_HEADER(24, 1, PHI1, PHI2, CTE_TOA_ALT) 
            ENDIF

            IF (FICDOWN_ANGLES_USER.NE.'NO_OUTPUT') THEN
               CALL SOS_OUTPUT_HEADER(25, 2, PHI1, PHI2, 0.D+00) 
            ENDIF
         
         ELSE   
            CALL SOS_OUTPUT_HEADER(14, 1, PHI1, PHI2, ZOUT)
            CALL SOS_OUTPUT_HEADER(15, 2, PHI1, PHI2, ZOUT)

            IF (FICUP_ANGLES_USER.NE.'NO_OUTPUT') THEN
               CALL SOS_OUTPUT_HEADER(24, 1, PHI1, PHI2, ZOUT) 
            ENDIF
            
            IF (FICDOWN_ANGLES_USER.NE.'NO_OUTPUT') THEN
               CALL SOS_OUTPUT_HEADER(25, 2, PHI1, PHI2, ZOUT) 
            ENDIF

         ENDIF

         !Half plane PHIos + 180 degrees: negative viewing angles by convention
         !----------------------------------------------------------------------

         DO J=LUM_NBMU,1,-1
	 
            JJ=J-1
	    	    
            !Recording for the upward field
            !------------------------------------
            XI = I_UP_FIN(0, JJ)
            XQ = Q_UP_FIN(0, JJ)
            XU = U_UP_FIN(0, JJ)
            XAN = POL_ANG_UP_FIN(0, JJ)
            TPOL = POL_RATE_UP_FIN(0, JJ)
            LPOL = L_POL_UP_FIN(0, JJ)
      
            !Nominal case: upward field for all the angles to output
            WRITE(14,55,err=992) 
     &         -THETA_FIN(JJ),SCA_UP_FIN(0, JJ),XI,XQ,XU,XAN,TPOL,LPOL

            !Optional case: upward field for the user angles only
            IF (FICUP_ANGLES_USER.NE.'NO_OUTPUT') THEN
               IF (IND_ANGOUT_FIN(JJ).EQ.1) THEN
                  WRITE(24,55,err=992) 
     &                -THETA_FIN(JJ),SCA_UP_FIN(0,JJ),XI,XQ,XU,XAN,
     &                 TPOL,LPOL
               ENDIF
            ENDIF  

           
            !Recording for the downward field
            !---------------------------------------
            XI = I_DOWN_FIN(0, JJ)
            XQ = Q_DOWN_FIN(0, JJ)
            XU = U_DOWN_FIN(0, JJ)
            XAN = POL_ANG_DOWN_FIN(0, JJ)
            TPOL = POL_RATE_DOWN_FIN(0, JJ)
            LPOL = L_POL_DOWN_FIN(0, JJ)

            !Nominal case: downward field for all the angles
            WRITE(15,55,err=992)  
     &           -THETA_FIN(JJ),SCA_DOWN_FIN(0, JJ),XI,XQ,XU,XAN,
     &           TPOL,LPOL   	     

            !Optional case: downward field for the user angles only
            IF (FICDOWN_ANGLES_USER.NE.'NO_OUTPUT') THEN
               IF (IND_ANGOUT_FIN(JJ).EQ.1) THEN
                  WRITE(25,55,err=992) 
     &                 -THETA_FIN(J-1),SCA_DOWN_FIN(0, JJ),XI,XQ,XU,
     &                 XAN,TPOL,LPOL
               ENDIF
            ENDIF 	    
         ENDDO

         DO J=1,LUM_NBMU
	 
            JJ = J -1

            !Recording for the upward field
            !------------------------------------
            XI = I_UP_FIN(1, JJ)
            XQ = Q_UP_FIN(1, JJ)
            XU = U_UP_FIN(1, JJ)
            XAN = POL_ANG_UP_FIN(1, JJ)
            TPOL = POL_RATE_UP_FIN(1, JJ)
            LPOL = L_POL_UP_FIN(1, JJ)

            !Nominal case: upward field for all the angles
            WRITE(14,55,err=992) 
     &            THETA_FIN(JJ),SCA_UP_FIN(1, JJ),XI,XQ,XU,XAN,
     &            TPOL,LPOL
     
            !Optional case: upward field for the user angles only
            IF (FICUP_ANGLES_USER.NE.'NO_OUTPUT') THEN
                IF (IND_ANGOUT_FIN(JJ).EQ.1) THEN
	            WRITE(24,55,err=992) 
     &                  THETA_FIN(JJ),SCA_UP_FIN(1, JJ),XI,XQ,XU,
     &                  XAN,TPOL,LPOL
                ENDIF
	     ENDIF  

          
            !Recording for the downward field
            !---------------------------------------
            XI = I_DOWN_FIN(1, JJ)
            XQ = Q_DOWN_FIN(1, JJ)
            XU = U_DOWN_FIN(1, JJ)
            XAN = POL_ANG_DOWN_FIN(1, JJ)
            TPOL = POL_RATE_DOWN_FIN(1, JJ)
            LPOL = L_POL_DOWN_FIN(1, JJ)
	    
            !Nominal case: downward field for all the angles to output
	     WRITE(15,55,err=992)  
     &            THETA_FIN(JJ),SCA_DOWN_FIN(1, JJ),XI,XQ,XU,XAN,
     &            TPOL,LPOL   	     

            !Optional case: downward field for the user angles only
            IF (FICDOWN_ANGLES_USER.NE.'NO_OUTPUT') THEN
                IF (IND_ANGOUT_FIN(JJ).EQ.1) THEN
	            WRITE(25,55,err=992) 
     &                  THETA_FIN(JJ),SCA_DOWN_FIN(1, JJ),XI,XQ,XU,XAN,
     &                  TPOL,LPOL
                ENDIF
	     ENDIF 	    

         ENDDO	    
	 
	    
      ENDIF  !ITRPHI=1
      


C* Output of a polar diagram
C---------------------------------
      IF(ITRPHI.EQ.2) THEN        
	 
         !Writing of the header information of the file
         !------------------------------------------------
         IF (ZOUT.EQ.-1) THEN   

            CALL SOS_OUTPUT_HEADER_POLAR_DIAG(14, 1, CTE_TOA_ALT) 
            CALL SOS_OUTPUT_HEADER_POLAR_DIAG(15, 2, 0.D+00)

            IF (FICUP_ANGLES_USER.NE.'NO_OUTPUT') THEN
               CALL SOS_OUTPUT_HEADER_POLAR_DIAG(24, 1, CTE_TOA_ALT) 
            ENDIF

            IF (FICDOWN_ANGLES_USER.NE.'NO_OUTPUT') THEN
               CALL SOS_OUTPUT_HEADER_POLAR_DIAG(25, 2, 0.D+00) 
            ENDIF
         
         ELSE   
            CALL SOS_OUTPUT_HEADER_POLAR_DIAG(14, 1, ZOUT)
            CALL SOS_OUTPUT_HEADER_POLAR_DIAG(15, 2, ZOUT)

            IF (FICUP_ANGLES_USER.NE.'NO_OUTPUT') THEN
               CALL SOS_OUTPUT_HEADER_POLAR_DIAG(24, 1,  ZOUT) 
            ENDIF

            IF (FICDOWN_ANGLES_USER.NE.'NO_OUTPUT') THEN
               CALL SOS_OUTPUT_HEADER_POLAR_DIAG(25, 2, ZOUT) 
            ENDIF

         ENDIF

         IP=0
         DO IPHI=0,360,PAS_PHI

            DO J=1,LUM_NBMU

               JJ=J-1
	       
               !Recording for the upward field
               !------------------------------------
               XI = I_UP_FIN(IP, JJ)
               XQ = Q_UP_FIN(IP, JJ)
               XU = U_UP_FIN(IP, JJ)
               XAN = POL_ANG_UP_FIN(IP, JJ)
               TPOL = POL_RATE_UP_FIN(IP, JJ)
               LPOL = L_POL_UP_FIN(IP, JJ)

               !Nominal case: upward field for all the angles
               WRITE(14,56,err=992) 
     &              PHI_FIN(IP),THETA_FIN(JJ),SCA_UP_FIN(IPHI, JJ),
     &              XI,XQ,XU,XAN,TPOL, LPOL
     
               !Optional case: upward field for the user angles only
               IF (FICUP_ANGLES_USER.NE.'NO_OUTPUT') THEN
                  IF (IND_ANGOUT_FIN(JJ).EQ.1) THEN
                     WRITE(24,56,err=992) 
     &                    PHI_FIN(IP),THETA_FIN(JJ),
     &                    SCA_UP_FIN(IP, JJ),XI,XQ,XU, XAN,TPOL,LPOL
                     
                  ENDIF
               ENDIF  

               !Nominal case: downward field for all the angles
               XI = I_DOWN_FIN(IP, JJ)
               XQ = Q_DOWN_FIN(IP, JJ)
               XU = U_DOWN_FIN(IP, JJ)
               XAN = POL_ANG_DOWN_FIN(IP, JJ)
               TPOL = POL_RATE_DOWN_FIN(IP, JJ)
               LPOL = L_POL_DOWN_FIN(IP, JJ) 
               WRITE(15,56,err=992) 
     &              PHI_FIN(IP),THETA_FIN(JJ),SCA_DOWN_FIN(IP, JJ),
     &              XI,XQ,XU,XAN,TPOL,LPOL

	 
               !Optional case: downward field for the user angles only
                IF (FICDOWN_ANGLES_USER.NE.'NO_OUTPUT') THEN
                  IF (IND_ANGOUT_FIN(JJ).EQ.1) THEN
                      WRITE(25,56,err=992) 
     &                    PHI_FIN(IP),THETA_FIN(JJ),
     &                    SCA_DOWN_FIN(IP, JJ),XI,XQ,XU,XAN,TPOL,LPOL
                    ENDIF
	        ENDIF  

      
            ENDDO !MU
	     
	     IP=IP+1
	     
	  ENDDO !PHI

      ENDIF   !ITRPHI=2


      
C* File closing
C-------------------
      CLOSE(14)
      CLOSE(15)
      CLOSE(21)

      IF (FICUP_ANGLES_USER.NE.'NO_OUTPUT') CLOSE(24)
      IF (FICDOWN_ANGLES_USER.NE.'NO_OUTPUT') CLOSE(25)	 
      
C* Nominal end
C-------------------
      IF (TRACE) THEN
         WRITE (IDLOG,*,ERR=1303) " "
         WRITE (IDLOG,*,ERR=1303) 
     &   "--> Nominal end of SOS_ABS_MAIN processes"
         WRITE (IDLOG,*,ERR=1303) "JOB_STATUS=OK"
         CLOSE(IDLOG)
      ENDIF
      
      GOTO 9999
      

C =========================================================
C --- Error case and return of status 1 to the shell
C =========================================================  

C --- Format of the inputs

 100  WRITE(6,*) '  SOS_MAIN : ERROR_100 on argument loading'
      WRITE(6,*) '  --> Check and correct the keyword : ', KEYWORD
      GOTO 9998  

C --- General parameters

 101  WRITE(6,*) '  SOS_MAIN : ERROR_101 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SOS_Main.Wa value : Wavelength for '
      WRITE(6,*) '       radiance simulation (microns) '
      GOTO 9998

C --- Parameters for the angle computations

 201  WRITE(6,*) '  SOS_MAIN : ERROR_201 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -ANG.Rad.NbGauss value : Number of Gauss'
      WRITE(6,*) '       angles for radiance computations (I4).' 
      GOTO 9998    

 202  WRITE(6,*) '  SOS_MAIN : ERROR_202 on argument loading'
      WRITE(6,*) '  --> Check required formats for -ANG.Thetas value :'
      WRITE(6,*) '      Solar zenithal angle (degrees).'
      GOTO 9998  

 203  WRITE(6,*) '  SOS_MAIN : ERROR_203 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -ANG.Aer.NbGauss value : Number of Gauss'
      WRITE(6,*) '       angles for Mie computations (I4).' 
      
C --- Aerosol parameters

 301  WRITE(6,*) '  SOS_MAIN : ERROR_301 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.Waref : Reference wavelength for the '
      WRITE(6,*) '       aerosol optical thickness .' 
      GOTO 9998       

 302  WRITE(6,*) '  SOS_MAIN : ERROR_302 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.AOTref : Aerosol optical thickness at '
      WRITE(6,*) '       reference wavelength (F9.5).' 
      GOTO 9998   
     
 303  WRITE(6,*) '  SOS_MAIN : ERROR_303 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.Tronca : Option of aerosols phase '
      WRITE(6,*) '       function truncature (I4).' 
      GOTO 9998   

 304  WRITE(6,*) '  SOS_MAIN : ERROR_304 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.Model : Type of aerosols model (I2).' 
      GOTO 9998                 

 306  WRITE(6,*) '  SOS_MAIN : ERROR_306 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.MMD.MRwa : Mono-modal aerosol model -'
      WRITE(6,*) '       Real part of the aerosol refractive index '
      WRITE(6,*) '       for the simulation wavelength (F5).' 
      GOTO 9998            

 307  WRITE(6,*) '  SOS_MAIN : ERROR_307 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.MMD.MIwa : Mono-modal aerosol model -'
      WRITE(6,*) '      Imaginary part of the aerosol refractive index'
      WRITE(6,*) '      for the simulation wavelength (F8).' 
      GOTO 9998             

 308  WRITE(6,*) '  SOS_MAIN : ERROR_308 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.MMD.MRwaref : Mono-modal aerosol model -'
      WRITE(6,*) '       Real part of the aerosol refractive index '
      WRITE(6,*) '       for the aerosol reference wavelength (F5).' 
      GOTO 9998            

 309  WRITE(6,*) '  SOS_MAIN : ERROR_309 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.MMD.MIwaref : Mono-modal aerosol model -'
      WRITE(6,*) '      Imaginary part of the aerosol refractive index'
      WRITE(6,*) '      for the aerosol reference wavelength (F8).' 
      GOTO 9998             

 310  WRITE(6,*) '  SOS_MAIN : ERROR_310 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.MMD.SDtype (I4) : Mono-modal aerosol model'
      WRITE(6,*) '       Model index : 1 LND '
      WRITE(6,*) '                     2 Junge Law '
      GOTO 9998            

 311  WRITE(6,*) '  SOS_MAIN : ERROR_311 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.MMD.LNDradius : Mono-modal aerosol model -'
      WRITE(6,*) '       Modal radius (microns) of the LND (Float)'
      GOTO 9998           
            
 312  WRITE(6,*) '  SOS_MAIN : ERROR_312 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.MMD.LNDvar : Mono-modal aerosol model -'
      WRITE(6,*) '       Variance of the LND (Float)'
      GOTO 9998           

 3121 WRITE(6,*) '  SOS_MAIN : ERROR_3121 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.MMD.JD.slope : Aerosol Junge model - '
      WRITE(6,*) '       Slope (Float)'
      GOTO 9998  

 3122 WRITE(6,*) '  SOS_MAIN : ERROR_3122 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.MMD.JD.rmin : Aerosol Junge model - '
      WRITE(6,*) '       Minimal radius in microns (Float)'
      GOTO 9998  

 3123 WRITE(6,*) '  SOS_MAIN : ERROR_3123 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.MMD.JD.rmax : Aerosol Junge model - '
      WRITE(6,*) '       Maximal radius in microns (Float)'
      GOTO 9998  
	              
 313  WRITE(6,*) '  SOS_MAIN : ERROR_313 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.WMO.Model : Type of WMO Model (I2)'
      GOTO 9998       	

 314  WRITE(6,*) '  SOS_MAIN : ERROR_314 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      Volumetric concentration of WMO user model :'
      WRITE(6,*) '      -AER.WMO.DL :  Dust-Like component (F9.5)'
      WRITE(6,*) '      -AER.WMO.WS :  Water Soluble component (F9.5)'
      WRITE(6,*) '      -AER.WMO.OC :  Oceanic component (F9.5)'
      WRITE(6,*) '      -AER.WMO.SO :  Soot component (F9.5)'
      GOTO 9998      	 

 315  WRITE(6,*) '  SOS_MAIN : ERROR_315 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.SF.Model : Type of Shettle & Fenn Model'
      WRITE(6,*) '       (I2)'
      GOTO 9998       	

 316  WRITE(6,*) '  SOS_MAIN : ERROR_316 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.SF.RH : % of humidity relative (F5.2)'
      GOTO 9998       	

 317  WRITE(6,*) '  SOS_MAIN : ERROR_317 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.VCdef : Type of selected method for'
      WRITE(6,*) '       bi-modal LND volumetric concentration'
      WRITE(6,*) '       definition (I2)'
      GOTO 9998       	
	
 318  WRITE(6,*) '  SOS_MAIN : ERROR_318 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.CoarseVC : User volumetric '
      WRITE(6,*) '       concentration for "LND coarse mode" (F8.5)'
      GOTO 9998       	 
	
 319  WRITE(6,*) '  SOS_MAIN : ERROR_319 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.FineVC : User volumetric '
      WRITE(6,*) '       concentration for "LND fine mode" (F8.5)'
      GOTO 9998       

 320  WRITE(6,*) '  SOS_MAIN : ERROR_320 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.RAOT : User value of the ration '
      WRITE(6,*) '       AOT_coarse / AOT_total '
      WRITE(6,*) '       for the aerosol reference wavelength (F9.5)'
      GOTO 9998       
      
 321  WRITE(6,*) '  SOS_MAIN : ERROR_321 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.CM.MRwa : '
      WRITE(6,*) '       Real part of the refractive '
      WRITE(6,*) '       index for the "LND coarse mode" '
      WRITE(6,*) '       at the wavelength of aerosol properties '
      WRITE(6,*) '       calculation (F5.3)'
      GOTO 9998       	 		

 322  WRITE(6,*) '  SOS_MAIN : ERROR_322 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.CM.MIwa : '
      WRITE(6,*) '       Imaginary part of the refractive '
      WRITE(6,*) '       index for the "LND coarse mode" '
      WRITE(6,*) '       at the wavelength of aerosol properties '
      WRITE(6,*) '       calculation (F8.5)'
      GOTO 9998       
	 
 323  WRITE(6,*) '  SOS_MAIN : ERROR_323 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.CM.MRwaref : '
      WRITE(6,*) '       Real part of the refractive '
      WRITE(6,*) '       index for the "LND coarse mode" '
      WRITE(6,*) '       at the aerosol reference wavelength'
      WRITE(6,*) '       calculation (F5.3)'
      GOTO 9998       	 		

 324  WRITE(6,*) '  SOS_MAIN : ERROR_324 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.CM.MIwaref : '
      WRITE(6,*) '       Imaginary part of the refractive '
      WRITE(6,*) '       index for the "LND coarse mode" '
      WRITE(6,*) '       at the aerosol reference wavelength'
      WRITE(6,*) '       calculation (F8.5)'
      GOTO 9998  

 325  WRITE(6,*) '  SOS_MAIN : ERROR_325 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.CM.SDradius : '
      WRITE(6,*) '       Modal radius of the "LND coarse mode"'
      WRITE(6,*) '       in microns (F9.5)'
      GOTO 9998       	 		
           
 326  WRITE(6,*) '  SOS_MAIN : ERROR_326 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.CM.SDvar : '
      WRITE(6,*) '       "LND coarse mode" variance (F9.5)'
      GOTO 9998       	 			

 327  WRITE(6,*) '  SOS_MAIN : ERROR_327 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.FM.MRwa : '
      WRITE(6,*) '       Real part of the refractive '
      WRITE(6,*) '       index for the "LND fine mode" '
      WRITE(6,*) '       at the wavelength of aerosol properties'
      WRITE(6,*) '       calculation (F5.3)'
      GOTO 9998       	 		

 328  WRITE(6,*) '  SOS_MAIN : ERROR_328 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.FM.MIwa : '
      WRITE(6,*) '       Imaginary part of the refractive '
      WRITE(6,*) '       index for the "LND fine mode" '
      WRITE(6,*) '       at the wavelength of aerosol properties'
      WRITE(6,*) '       calculation (F8.5)'
      GOTO 9998       
	 
 329  WRITE(6,*) '  SOS_MAIN : ERROR_329 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.FM.MRwaref : '
      WRITE(6,*) '       Real part of the refractive '
      WRITE(6,*) '       index for the "LND fine mode" '
      WRITE(6,*) '       at the aerosol reference wavelength'
      WRITE(6,*) '       calculation (F5.3)'
      GOTO 9998       	 		

 330  WRITE(6,*) '  SOS_MAIN : ERROR_330 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.FM.MIwaref : '
      WRITE(6,*) '       Imaginary part of the refractive '
      WRITE(6,*) '       index for the "LND fine mode" '
      WRITE(6,*) '       at the aerosol reference wavelength'
      WRITE(6,*) '       calculation (F8.5)'
      GOTO 9998  

 331  WRITE(6,*) '  SOS_MAIN : ERROR_331 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.FM.SDradius : '
      WRITE(6,*) '       Modal radius of the "LND fine mode"'
      WRITE(6,*) '       in microns (F9.5)'
      GOTO 9998       	 		
           
 332  WRITE(6,*) '  SOS_MAIN : ERROR_332 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AER.BMD.FM.SDvar : '
      WRITE(6,*) '       "LND fine mode" variance (F9.5)'
      GOTO 9998       	 				

C --- Profile parameters

 401  WRITE(6,*) '  SOS_MAIN : ERROR_401 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.MOT : Molecular atmos optical thickness '
      WRITE(6,*) '    (for the radiance simulation wavelength) (F9.5).'
      GOTO 9998   

 402  WRITE(6,*) '  SOS_MAIN : ERROR_402 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.HR : Height scale of the molecular '
      WRITE(6,*) '        atmospheric profile (F9).' 
      GOTO 9998   

 403  WRITE(6,*) '  SOS_MAIN : ERROR_403 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.AerHS.HA : Height scale of the aerosol '
      WRITE(6,*) '       atmospheric profile (F9).' 
      GOTO 9998   

 404  WRITE(6,*) '  SOS_MAIN : ERROR_404 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.AerProfile.Type : type of '
      WRITE(6,*) '        atmospheric profile (I2).' 
      GOTO 9998   

 405  WRITE(6,*) '  SOS_MAIN : ERROR_405 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.AerLayer.Zmin : Minimal altitude (km)'
      WRITE(6,*) '        of aerosol layer (F9.5).' 
      GOTO 9998   

 406  WRITE(6,*) '  SOS_MAIN : ERROR_406 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.AerLayer.Zmax : Maximal altitude (km)'
      WRITE(6,*) '       of aerosol layer (F9.5).' 
      GOTO 9998   

 407  WRITE(6,*) '  SOS_MAIN : ERROR_407 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.Psurf : Atmospheric pressure at sea '
      WRITE(6,*) '       level (mbar) (F9.5).' 
      GOTO 9998   

 408  WRITE(6,*) '  SOS_MAIN : ERROR_408 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.H2O : Water vapor content (g/cm2) (F9) '
      GOTO 9998   

 409  WRITE(6,*) '  SOS_MAIN : ERROR_409 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.O3 : Ozone content (Dobson) (F9) '
      GOTO 9998   

 4091 WRITE(6,*) '  SOS_MAIN : ERROR_4091 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.CO2 : Carbon dioxyde concentration'
      WRITE(6,*) '                at surface level (ppmv) (F9)'
      GOTO 9998   

 4092 WRITE(6,*) '  SOS_MAIN : ERROR_4092 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.CH4 : Methane concentration'
      WRITE(6,*) '                at surface level (ppmv) (F9)'
      GOTO 9998   

 410  WRITE(6,*) '  SOS_MAIN : ERROR_410 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.AbsProfile.Type : ',
     &		 'Absorption profile type (I2)'
      GOTO 9998   

 411  WRITE(6,*) '  SOS_MAIN : ERROR_411 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -AP.SpectralResol : ',
     &		 'Spectral resolution in cm-1 of CKD datafiles (Float)'
      GOTO 9998   

C --- Surface parameters

 501  WRITE(6,*) '  SOS_MAIN : ERROR_501 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SURF.Type : Type of surface model (I2) '
      GOTO 9998   

 502  WRITE(6,*) '  SOS_MAIN : ERROR_502 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SURF.Ind : Refraction indice (F5.3) '
      GOTO 9998   

 503  WRITE(6,*) '  SOS_MAIN : ERROR_503 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SURF.Glitter.Wind : Wind speed (m/s) (F4) '
      GOTO 9998   

 504  WRITE(6,*) '  SOS_MAIN : ERROR_504 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SURF.Roujean.K0 : Roujean K0 parameter (F7.3)'
      GOTO 9998   

 505  WRITE(6,*) '  SOS_MAIN : ERROR_505 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SURF.Roujean.K1 : Roujean K1 parameter (F7.3)'
      GOTO 9998   

 506  WRITE(6,*) '  SOS_MAIN : ERROR_506 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SURF.Roujean.K2 : Roujean K2 parameter (F7.3)'
      GOTO 9998   

 507  WRITE(6,*) '  SOS_MAIN : ERROR_507 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SURF.Nadal.Alpha : Nadal Alpha parameter (F6)'
      GOTO 9998   
 
 508  WRITE(6,*) '  SOS_MAIN : ERROR_508 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SURF.Nadal.Beta : Nadal Beta parameter (F6) '
      GOTO 9998   
 
 5081 WRITE(6,*) '  SOS_MAIN : ERROR_5081 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SURF.Maignan.C : C*exp(-NDVI) parameter (F6) '
      WRITE(6,*) "      of the Maignan's surface BPDF model "
      GOTO 9998 

 509  WRITE(6,*) '  SOS_MAIN : ERROR_509 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SURF.Alb : Surface Albedo'
      GOTO 9998   

C --- SOS parameters

 601  WRITE(6,*) '  SOS_MAIN : ERROR_601 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SOS.OutputAlt : Output Altitude in km (F9) '
      GOTO 9998   

 602  WRITE(6,*) '  SOS_MAIN : ERROR_602 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SOS.IGmax : maximal order of interaction (I4)'
      GOTO 9998	 

 603  WRITE(6,*) '  SOS_MAIN : ERROR_603 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '   -SOS.Ipolar : Polarization activation option (I2)'
      GOTO 9998	 

 604  WRITE(6,*) '  SOS_MAIN : ERROR_604 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SOS.View : Output Type (I4)'
      GOTO 9998	 

 605  WRITE(6,*) '  SOS_MAIN : ERROR_605 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SOS.View.Phi : Relative azimut angle (F9)'
      GOTO 9998	 

 606  WRITE(6,*) '  SOS_MAIN : ERROR_606 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SOS.View.Dphi : Azimut sampling step (I4)'
      GOTO 9998	 

 608  WRITE(6,*) '  SOS_MAIN : ERROR_608 on argument loading'
      WRITE(6,*) '  --> Check required formats for '
      WRITE(6,*) '      -SOS.AbsModeCKD : Selected ',
     & 		 ' mode for application of CKD method factor (I4)'
      GOTO 9998	 

C --- Existence of the directory tree parameters

 1000 WRITE(6,*) '  SOS_MAIN : ERROR_1000 on parameters'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SOS_Main.ResRoot Value ' 
      WRITE(6,*) '  ==>  Working folder of SOS computations'
      GOTO 9998 

C --- Definition of the output directory tree

 1100 WRITE(6,*) '  SOS_MAIN : ERROR_1100 : '
      WRITE(6,*) '  => Working folder name : maybe too long'
      WRITE(6,*) '     -SOS_Main.ResRoot parameter'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998                  	             	             	                 	              	

 1108 WRITE(6,*) '  SOS_MAIN : ERROR_1108 : '
      WRITE(6,*) '  => SOS Working folder name : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998                  	               	  

 1109 WRITE(6,*) '  SOS_MAIN : ERROR_1109 : '
      WRITE(6,*) '  => Main log folder name : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENDIR has to be adjusted.' 
      GOTO 9998     

C --- log management
                                           
 1301 WRITE(6,*) '  SOS_MAIN : ERROR_1301 : '
      WRITE(6,*) '  => Log filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	

 1302 WRITE(6,*) '  SOS_MAIN : ERROR_1302 on Log file opening'
      GOTO 9998

 1303 WRITE(6,*) '  SOS_MAIN : ERROR_1303 on Log file writing'
      GOTO 9998
            

C --- Check of the profile parameters

 2515 WRITE(6,*) '  SOS_MAIN : ERROR_2515 on parameters'
      WRITE(6,*) '  --> -SOS.AbsModeCKD value must be defined'
      WRITE(6,*) '      1 : Fine calculation of the impact of gaz'
      WRITE(6,*) '          absorption on radiance by the CKD method.'
      WRITE(6,*) '      2 : Simplified calculation of the impact of'
      WRITE(6,*) '           gaz absorption radiance.'
      GOTO 9998  
      

C --- Check of the SOS parameters

 2605 WRITE(6,*) '  SOS_MAIN : ERROR_2605 on parameters'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SOS.View  (I4) '
      WRITE(6,*) '  ==> Azimut mode  '
      GOTO 9998     

 2607 WRITE(6,*) '  SOS_MAIN : ERROR_2607 on parameters'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SOS.View.Phi  (F9) '
      WRITE(6,*) '  ==> Relative azimut (required when SOS.View=1) '      
      GOTO 9998     
 
 2608 WRITE(6,*) '  SOS_MAIN : ERROR_2608 on parameters'
      WRITE(6,*) '  --> Required (-Keyword Value) definition : ',
     &           '     -SOS.View.Dphi  (Integer) '
      WRITE(6,*) '  ==> Azimut sampling step'
      WRITE(6,*) '      (required when SOS.View=2)'    
      GOTO 9998     

C --- Definition of the output files

 3013 WRITE(6,*) '  SOS_MAIN : ERROR_3013 : '
      WRITE(6,*) '  => Upward radiance filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	

 3014 WRITE(6,*) '  SOS_MAIN : ERROR_3014 : '
      WRITE(6,*) '  => Downward radiance filename : maybe too long'
      WRITE(6,*) '  => Check if CTE_LENFIC1 has to be adjusted.' 
      GOTO 9998                  	
     
      
C* Error case: return of status 1 to the calling program
C------------------------------------------------------
  991 WRITE(6,*) "SOS_TRPHI : Error while opening a file"
      GOTO 9998

  992 WRITE(6,*) 
     &  "SOS_TRPHI : Error while reading or writing on a file"
      GOTO 9998

                    
 9998 IF (TRACE) THEN
         !OPEN (IDLOG,FILE=FICMAIN_LOG,STATUS='OLD',ERR=1302)
         WRITE (IDLOG,*,ERR=1303) " "
         WRITE (IDLOG,*,ERR=1303) "--> Processing error"
	  WRITE (IDLOG,*,ERR=1303) "JOB_STATUS=ERROR"
	  CLOSE(IDLOG)
      ENDIF 

      CALL FLUSH(6)        !Forces the flush of the buffer to the screen

      CALL EXIT(1)    
      
 9999 STOP



C =========================================================
C --- Formats
C =========================================================            
C
     
C --- Formats for the files SOS_Up and SOS_Down
55    FORMAT(2(2X,F7.2),2X,3(E13.6,2X),2(F7.2,2X),E13.6)
56    FORMAT(3(2X,F7.2),1X,3(E13.6,2X),2(1X,F7.2),E13.6)

     
      END   	!END OF PROGRAM SOS_MAIN
