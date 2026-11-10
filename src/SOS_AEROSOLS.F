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
C* FILE: SOS_AEROSOLS.F
C* PROJECT: Successive orders of scattering
C* ROLE: This module manages the computation of the radiative properties of
C*       different aerosol size distributions.
C*
C* AUTHOR:
C*   Initial code: Laboratoire d'Optique Atmospherique (LOA, Villeneuve d'Ascq).
C*   Code rework and development: CS GROUP France, B.Lafrance, X.Lenot
C*
C* DATE: 30/04/2002
C*
C* MOD:VERSION:1.1: Upgrade of the LOA OS code
C* MOD:VERSION:1.2: Update of the comments.
C*
C* MOD:VERSION:2.0: 24/08/2008 
C*       - Introduction of a bi-LND modelling.
C*       - Conversion of the outputs KMAT1 and KMAT2 of routine SOS_GRANU
C*         into extinction and scattering cross sections (micron^2)
C*         instead of extinction and scattering coefficients (/mic)
C*         that were not normalised to one particle.
C*       - Update of the comments.
C*
C*
C* MOD:VERSION:2.1: 08/09/2008 
C*       - Correction of indentation overflows for f77 compilation.
C*       - Removal of the loop exit when the end of the MIE file is reached.
C*         Indeed, if this happens, it is an error: the MIE file is
C*         incomplete (probably lack of disk space when it was generated).
C*         This case is then correctly handled as an error.
C*       - Correction of the computation of the proportions of the components of the
C*         Dubovik bimodal models: case of the computation from the optical thickness ratios
C*         of the two modes (fine and coarse) and from the total optical thickness for a reference
C*         wavelength ==> use of the refractive index of the
C*         particles for the reference wavelength (only).
C*
C* MOD:VERSION:3.0: 19/11/2009 
C*       - Change of the argument passing:
C*       - passing by pairs (keyword, value).
C*       - adaptation of the error case management.
C*       - Error messages converted to English.
C*       - Introduction of the parameter CTE_MIE_ALPHAMIN giving the minimum value
C*         of the size parameter ALPHAO of the Mie computations.
C*
C* MOD:VERSION:4.0: 10/01/2010 
C*   - Removal of the constant CTE_OS_FICANGLE (from SOS.h)
C*     in favour of passing the name of the angle file (produced by SOS_ANGLES) as argument.
C*
C*   - Change of the size of the radiance arrays:
C*     ==> Changed from CTE_MIE_NBMU to CTE_MIE_NBMU_MAX.
C*     ==> Loop on the angles up to MIE_NBMU, equal to the actual number of angles used.
C*                                   
C*   - Change of the size of the phase function arrays:
C*     ==> Changed from CTE_MIE_NBMU to CTE_MIE_NBMU_MAX.
C*     ==> Loop on the angles up to MIE_NBMU, equal to the actual number of angles used.
C*
C*   - Change of the size of the arrays of the decomposition in Legendre polynomials:
C*     ==> Changed from CTE_OS_NB to CTE_OS_NB_MAX.
C*     ==> Loop on the orders up to OS_NB (read in the angle file)
C*
C*   - CTE_MIE : Added input arguments: NBMU_MIE (nb of useful angles) and XMU (array of angles)
C*
C*   - SOS_NOMFIC_MIE : Added input arguments: NBMU_GAUSS (nb of Gauss angles)
C*                                                     and FICANGLES_USER (user angle file)
C*
C*   - SOS_GRANU:
C*         * Added input argument: NBMU_MIE (nb of useful angles)
C*         * Change of the reading line of the arrays
C*           of phase functions per value of the size parameter Alpha
C*           ==> limited to the useful bounds -MIE_NBMU:MIE_NBMU
C*
C*   - SOS_DECOMPO_LEGENDRE:
C*         * Added input arguments: NBMU_MIE (nb of useful angles) and OS_NB
C*         * Loops limited to NBMU_MIE instead of CTE_MIE_NBMU_MAX
C*  
C*   - Removal of the use of the environment variable SOS_RACINE_MU
C*     ==> becomes SOS_ABS_FIC since it is only used for the WMO and S&F files.
C*
C*   - Addition of a processing option IMOD=4 to ingest an external phase matrix
C*     --> External phase functions and scattering and extinction coefficients
C*     ==> Introduction of the routines:
C*            	- SOS_INTERPO_SPLINT
C*     		- SOS_SPLINE
C*    		- SOS_SPLINT
C*     ==> Addition of the constant (SOS.h): CTE_MAXNB_ANG_EXT
C*
C* MOD:VERSION:5.0: 06/06/2013
C*     
C*   -  The input parameters are given as a list of parameters
C*      and no longer as pairs with associated keyword => all
C*	 parameters are therefore mandatory
C*
C*   -  Change of the value of Pi, declared as a constant in the file
C*      (INCTE_PI replaces CTE_PI)
C* 
C* MOD:VERSION:6.0: 02/03/2016
C*
C*    - Minor adjustments to strictly respect the 72 columns
C*      (including when decoding the constants): required for gfortran compilation
C*
C*    - Change of the way the formatted files are read
C*      (Angles, external data file for the aerosols)
C*      to be compatible with the gfortran compiler.
C*
C*    - Change of the format of the result GRANU file to allow its reading
C*      by the code SOS_PREPA_OS.F compiled with gfortran and to clarify its content
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
C* MOD:VERSION:6.1: 12/06/2020
C* 
C*    - Addition of the parameter DIR_TMP defining the temporary directory of the
C*      Mie computations.
C*
C*    - Definition of a temporary Mie file FICMIE_TMP generated
C*      under DIR_TMP before being moved to the storage directory
C*      of the Mie files (DIRS) with the formalised file name (FICMIE).
C*
C*    - Addition of the argument VARGRANU3_MMD giving the max radius of a
C*      size distribution following the Junge law.
C*      VARGRANU1 is renamed VARGRANU1_MMD and VARGRANU2 VARGRANU2_MMD
C*      ==> Adaptation of the call of CTE_GRANU
C*
C*    - The maximum size parameter of the Mie computations (ALPHAF)
C*      for a mono-modal model is no longer given but computed by the routine.
C*
C*    - Addition of a processing option IMOD=5 to use the definition
C*      of a mixture of modes given in a file by the user.
C*      ==> Introduction of the parameter FICMIXTURE_AER, which gives the location of the file.
C*      ==> Implementation of the processing block for IMOD=5.
C*
C*
C* MOD:VERSION:6.2: 31/08/2022
C*    - Removal of the use of the environment variable SOS_ABS_FIC.
C*    - Constants renamed from SOS_* to CTE_*
C*    - Traces converted to English.
C*    - Introduction of the constants CTE_AER_MU1_TRONCA and CTE_AER_MU2_TRONCA
C*      to define the angles used to estimate the truncation of the forward
C*      scattering peak of the aerosol phase function.
C*    - The information about the units for extinction and scattering coefficients,
C*      and for the number of particles, are clarified in the log file.
C*
C* MOD:VERSION:6.3: 08/08/2023
C*    - Addition of a test on the sum of the optical thickness rates defining an aerosol mixture
C*       ==> Normalisation applied if the difference to 1 is below the threshold CTE_GAP_TOLER_SUM_RATES
C*    - Clarification of the name of the variable giving the truncation coefficient (COEF_TRONCA)
C*      and this variable becomes an output parameter of routine SOS_AEROSOLS.
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
C    CTE_LENCOM: length of the system command strings
C    INCTE_PI : Value of PI.
C    CTE_MIE_NBMU_MAX : Maximum number of positive angles.
C    CTE_MIE_ALPHAMIN : Default minimum value of the size parameter
C                       Alpha of the Mie computations.
C    CTE_ALPHAMAX_WMO_DL : Limit value of the size parameter for the computation
C                          of the files of the WMO Dust-Like component.
C    CTE_ALPHAMAX_WMO_WS : Same for the Water Soluble component.
C    CTE_ALPHAMAX_WMO_OC : Same for the Oceanic component.
C    CTE_ALPHAMAX_WMO_SO : Same for the Soot component.
C    CTE_ALPHAMAX_SF_SR : Same for the Small Rural component of Shettle & Fenn.
C    CTE_ALPHAMAX_SF_SU : Same for the Small Urban component of Shettle & Fenn.
C    CTE_COEF_NRMAX : Constant for the determination of the max size parameter.
C    CTE_OS_NB_MAX :  Maximum value of the limit order of the expansion in
C                     Legendre polynomials.
C    CTE_PH_SEUIL_TRONCA : Threshold for the truncation.
C    CTE_AER_DATAWMO : Name of the file containing the WMO parameters (rm,sig,...).
C    CTE_AER_DATASF  : Name of the file containing the S&F modal radius and log variance.
C    CTE_AER_XX_SF   : Name of the files containing the S&F refractive indices
C		       for each component XX = SR, LR, SU, LU, OM.
C    CTE_MAXNB_ANG_EXT : Maximum size of the arrays of external phase functions
C    CTE_AER_MU1_TRONCA : cos of the angle Theta1 for the estimation of the segment of
C                         truncation of the aerosol phase function.
C    CTE_AER_MU2_TRONCA : cos of the angle Theta2 for the estimation of the segment of
C                         truncation of the aerosol phase function.
C
C
C Specific constants:
C    LENLINE  : Maximum length of a data line read in a file
C    CTE_NOT_DEFINED_VALUE_DBLE: default value for reals
C----------------------------------------------------------------------------
#include "SOS.h"
#define INCTE_PI DACOS(-1.D+00)
#define LENLINE 1000

C==============================================================================
C PROGRAM: SOS_AEROSOLS
C ==========  
C  Simulated models:
C  ----------------
C      This program computes the radiative properties of a mono-modal,
C      multi-modal (WMO) or bi-modal (S&F) aerosol size distribution.
C
C      The basic components of the WMO models are:
C      dust-like, oceanic, water-soluble and soot.
C      Four typical WMO models can be simulated: continental, maritime,
C      urban and a composition defined by the user.
C
C      The basic components of the Shettle & Fenn models are:
C      small rural, large rural, small urban, large urban and oceanic.
C      The typical S&F models that can be simulated are: tropospheric, urban,
C      maritime and coastal.
C
C      If imodalite = 0 => a single aerosol size distribution is used
C			   (Log Normal or Junge law).
C
C      If imodalite = 1 => case of multi-modal size distributions, use
C                          of the four standard WMO models.
C
C      If imodalite = 2 => case of bi-modal size distributions, use
C                          of the standard S&F models.
C
C      If imodalite = 3 => case of bi-modal LND size distributions.
C
C      If imodalite = 4 => case of the use of an external phase matrix
C                         (phase functions in intensity and polarisation).
C
C      If imodalite = 5 => case of the use of a user file
C                          defining a mixture of aerosol models.
C
C
C  Mie computations and application:
C  -----------------------------
C      The MIE computations give the radiative properties of a component.
C      They are made by the procedure "SOS_MIE".
C      The MIE files are stored in a specific directory.
C
C      The step on the size parameter alpha is:
C         if        alpha  <= 0.1  ===> step = 0.0001
C         if 0.1  < alpha  <= 1.0  ===> step = 0.001
C         if 1.0  < alpha  <= 10   ===> step = 0.01
C         if 10   < alpha  <= 30   ===> step = 0.05
C         if 30   < alpha  <= 100  ===> step = 0.1
C         if 100  < alpha          ===> step = 1.0
C
C      From the Mie computations, the phase function
C      P(mu) and the polarised phase functions Q(mu) and U(mu) are determined for a
C      given wavelength and a given particle size distribution
C      for each model.
C
C
C      A truncation is applied to the phase function P(mu) if the option
C      itronc = 1 is selected.
C
C      The phase functions are then expanded in the basis of the
C      Legendre polynomials. This gives:
C         - For P(mu) => the coefficients beta(k) from 0 to OS_NB,
C         - For Q(mu) => the coefficients gamma(k) from 0 to OS_NB,
C         - For U(mu) => the coefficients delta(k) from 0 to OS_NB.
C
C      The coefficients alp(k) and zeta(k) are computed from beta(k)
C      and delta(k). They are used for the expression of the phase matrix
C      decomposed in Fourier series (in the OS code).
C
C 
C  Remark:
C  --------    
C      If the aerosol optical thickness is null, an Aerosols file
C      is created with values set to zero.
C
C
C Input data of the program
C ------------------------------
C    -->Parameters/arguments of the program
C
C             FICANGLES (CHARACTER) : (I) complete name of the file of the angles used
C                                         for the radiance and BRDF computations
C				  (directory + file name + extension)
C
C      	TA   (DOUBLE)  : (I) Aerosol optical thickness (not truncated).
C                             (decoding format: real F9.5)
C      		       
C	       WA   (DOUBLE)  : (I) Simulation wavelength in microns.
C                             (decoding format: real F9.5)
C
C             WAREF (DOUBLE) : (I) reference wavelength WAREF.
C                             (decoding format: real F9.5)
C	 
C      	TAWAREF (DOUBLE) : (I) total AOT at the wavelength WAREF.
C                             (decoding format: real F9.5)
C     		       
C	       ITRONC  (I4)   : (I) Truncation option (1 to apply it).
C			       (decoding format: integer I4)
C
C	       IMOD  (I4)  : (I) Mono-modal or multi-modal size distribution.
C			       (decoding format: integer I2)
C      		        	0 : Mono-modal models
C				       1 : WMO models
C				       2 : Shettle & Fenn models
C		        	       3 : Bimodal LND models
C                       	       4 : Use of external phase functions
C                                  5 : Use of a file defining a mixture of modes
C
C	       RN   (DOUBLE)  : (I) ==> Associated with IMOD = 0
C				   Real part of the refractive index for a mono-modal model
C                                  at the radiance simulation wavelength (WA)
C                                 (decoding format: real F5.3)
C
C             IN   (DOUBLE)  : (I) ==> Associated with IMOD = 0
C				   Imaginary part of the refractive index for a mono-modal model
C                                  at the radiance simulation wavelength (WA)
C                                (decoding format: real F8.5)
C	  
C             IGRANU  (I4)   : (I) ==> Associated with IMOD = 0
C			      	   Index of the size distribution type for a mono-modal model
C                          		1 : LND (parameters: modal radius in microns and variance).
C                          		2 : Junge law (parameters: minimum radius
C                              		in microns and power applied to the radius).
C			       (decoding format: integer I4)
C
C             VARGRANU1_MMD  (DOUBLE) (I) :  1st size distribution parameter for a mono-modal model
C             VARGRANU2_MMD  (DOUBLE) (I) :  2nd size distribution parameter for a mono-modal model
C             VARGRANU3_MMD  (DOUBLE) (I) :  3rd size distribution parameter for a mono-modal model
C	 		   If LND       : modal radius and standard deviation
C	  		   If Junge law : min radius, slope and max radius
C
C             IMODELE_WMO   (I4)  : (I) ==> Associated with IMOD = 1
C				    Index of the choice of WMO model
C			       (decoding format: integer I2)
C      		        		1 : Continental WMO model
C					2 : Maritime WMO model
C					3 : Urban WMO model
C					4 : WMO model "user definition"
C
C             C_WMO_DL (DOUBLE)  : (I) ==> Associated with IMOD = 1
C				       Volume percentage C of the DL particles (WMO).
C                             (decoding format: real F9.5)
C
C             C_WMO_WS (DOUBLE)  : (I) ==> Associated with IMOD = 1
C				       Volume percentage C of the WS particles (WMO).
C                             (decoding format: real F9.5)
C
C             C_WMO_OC (DOUBLE)  : (I) ==> Associated with IMOD = 1
C				       Volume percentage C of the OC particles (WMO).
C                             (decoding format: real F9.5)
C
C             C_WMO_SO (DOUBLE)  : (I) ==> Associated with IMOD = 1
C				       Volume percentage C of the SO particles (WMO).
C                             (decoding format: real F9.5)
C
C             IMODELE_SF (DOUBLE) : (I) ==> Associated with IMOD = 2
C					Index of the choice of Shettle & Fenn model
C			       (decoding format: integer I2)
C      		        		1 : Tropospheric SF model
C					       2 : Urban SF model
C					       3 : Maritime SF model
C					       4 : Coastal SF model
C
C             RH  (DOUBLE) : (I) ==> Associated with IMOD = 2
C				 Relative humidity for the Shettle & Fenn model
C				(decoding format: real F5.2)
C
C             MODE_PARAM_BILND (I4) : (I) ==> Associated with IMOD = 3
C					  Index of the type of description of the bimodal LND mixture.
C			       (decoding format: integer I2)
C      		        		 1 : Use of the user volume concentrations
C					        2 : Use of the ratio AOT_coarse / AOT_total
C
C             USER_CV_COARSE (DOUBLE) : (I) ==> Associated with IMOD = 3 and MODE_PARAM_BILND = 1
C      	     					User volume concentration for the "LND coarse" mode.
C                                		       (decoding format: real F8.5)
C	 
C             USER_CV_FINE (DOUBLE) : (I) ==> Associated with IMOD = 3 and MODE_PARAM_BILND = 1
C      	     			      	  User volume concentration for the "LND fine" mode.
C                                	         (decoding format: real F8.5)
C	 
C             RTAUct_WAREF (DOUBLE) : (I) ==> Associated with IMOD = 3 and MODE_PARAM_BILND = 2
C      	     			     	  ratio AOT_coarse / AOT_tot at the wavelength WAREF.
C                             (decoding format: real F9.5)
C
C             BMD_CM_MRwa  (DOUBLE) : (I) ==> Associated with IMOD = 3
C      	     				  Real part of the refractive index
C                                   	  of the particles of the "LND coarse" mode for the simulation WA
C                               	         (decoding format: real F5.3)
C
C             BMD_CM_MIwa  (DOUBLE) : (I) ==> Associated with IMOD = 3
C      	     				  Imaginary part of the refractive index
C                                   	  of the particles of the "LND coarse" mode for the simulation WA
C                               	         (decoding format: real F8.5)
C
C             BMD_CM_MRwaref  (DOUBLE) : (I) ==> Associated with IMOD = 3
C      	     				  Real part of the refractive index
C                                   	     of the particles of the "LND coarse" mode for WAREF
C                               	            (decoding format: real F5.3)
C
C             BMD_CM_MIwaref (DOUBLE) : (I) ==> Associated with IMOD = 3
C      	     				  Imaginary part of the refractive index
C                                   	    of the particles of the "LND coarse" mode for WAREF
C                               	           (decoding format: real F8.5)
C
C             BMD_CM_RMODAL (DOUBLE) : (I)  ==> Associated with IMOD =  3
C      	     				    Modal radius of the "LND coarse"
C                             (decoding format: real F9.5)
C
C             BMD_CM_VAR (DOUBLE) : (I) ==> Associated with IMOD = 3
C      	     				    Variance of the "LND coarse"
C                             (decoding format: real F9.5)
C
C             BMD_FM_MRwa  (DOUBLE) : (I) ==> Associated with IMOD = 3
C      	     				  Real part of the refractive index
C                                   	  of the particles of the "LND fine" mode for the simulation WA
C                               	  (decoding format: real F5.3)
C
C             BMD_FM_MIwa  (DOUBLE) : (I) ==> Associated with IMOD = 3
C      	     				  Imaginary part of the refractive index
C                                   	  of the particles of the "LND fine" mode for the simulation WA
C                               	         (decoding format: real F8.5)
C
C             BMD_FM_MRwaref  (DOUBLE) : (I) ==> Associated with IMOD = 3
C      	     				  Real part of the refractive index
C                                   	     of the particles of the "LND fine" mode for WAREF
C                               	            (decoding format: real F5.3)
C
C             BMD_FM_MIwaref (DOUBLE) : (I) ==> Associated with IMOD = 3
C      	     				  Imaginary part of the refractive index
C                                   	    of the particles of the "LND fine" mode for WAREF
C                               	           (decoding format: real F8.5)
C
C             BMD_FM_RMODAL (DOUBLE) : (I)  ==> Associated with IMOD =  3
C      	     				    Modal radius of the "LND fine"
C                             (decoding format: real F9.5)
C
C             BMD_FM_VAR (DOUBLE) : (I) ==> Associated with IMOD = 3
C      	     				    Variance of the "LND fine"
C                             (decoding format: real F9.5)
C
C             FICEXTDATA (CHARACTER) : (I) ==> Associated with IMOD = 4
C                     			   Name of the external data file of phase functions
C                     			   and radiative parameters (extinction and scattering coefficients)
C				  (directory + file name + extension)
C
C             FICMIXTURE_AER (CHARACTER) : (I) ==> Associated with IMOD = 5
C                     			   Name of the file defining a mixture of aerosol modes
C				  (directory + file name + extension)
C
C             DIR_TMP  (CHARACTER)  : (I) temporary directory where the files are generated
C
C             DIRS   (CHARACTER)  : (I) storage directory of the MIE files
C
C             FICLOGMIE (CHARACTER) : (I) MIE trace file
C				              (0 or complete path to the trace file)
C				  (directory + file name + extension)
C
C             FICGRANU (CHARACTER) : (I) complete name of the GRANU result file
C				  (directory + file name + extension)
C 
C             FICLOGGRANU (CHARACTER) : (I) Trace file of the Aerosols program
C				              (0 or complete path to the trace file)
C				  (directory + file name + extension)
C
C      	KMAT1 (DOUBLE) : (O) Extinction cross section (micron^2). Not truncated.
C
C      	PIZ (DOUBLE) : (O) Single scattering albedo.
C
C      	COEF_TRONCA (DOUBLE) : Truncation coefficient of the phase function (O)
C
C     IER (I4) : (O) error code =0 if no error, =-1 otherwise
C 
C                          
C
C    -->Files
C       The MIE file
C       The Granu trace file
C       The Mie trace file
C       The file of angles and internal parameters for the expansion limits.
C       The file of external phase functions (optional input)
C       The file defining a mixture of aerosol modes (optional input)
C
C Results provided by the program
C ----------------------------------
C   GRANU result file containing:
C	- extinction and scattering coefficient, Beta(0)/3
C	- truncation coefficient
C	- truncated single scattering albedo
C	- the parameters of the decomposition of the phase function in
C                     Legendre polynomials.
C
C Environment variable
C ------------------------
C   SOS_ABS_ROOT : defines the access path to the main directory of the SOS_ABS code
C   This environment variable is used to locate the files under $SOS_ABS_ROOT/fic
C
C
C Description of the files used
C ---------------------------------
C   -->Content of the MIE file
C       Unformatted binary file, created by SOS_MIE
C       1 header record with
C                Refractive index (real part)
C                Refractive index (imaginary part)
C                Limit size parameter
C                Number of positive Gauss angles
C
C      1 record for each value of the size parameter with:
C        - Size parameter (real),
C        - Extinction efficiency factor (real),
C        - Scattering efficiency factor (real),
C        - Asymmetry factor (double precision),
C        - Array of the phase function Mie(mu) (real),
C        - Array of the polarised phase function Qmie(mu) (real),
C        - Array of the polarised phase function Umie(mu) (real).
C
C    --> The arrays have a useful size limited to MIE_NBMU:
C        they are created and read with this limitation.
C
C
C    -->Content of the angle file:
C       Formatted ascii file, created by SOS_ANGLES.
C       It contains the information on the angles to be used
C       and on the adapted expansion orders:
C		NB_TOTAL_ANGLES :    total nb of angles to be used
C               NB_GAUSS_ANGLES :  nb of Gauss angles
C               ANGLES_USERFILE :  Name of the user angles file
C               INTERNAL_OS_NB :   Order of the expansions of the phase functions
C		           in Legendre polynomials to be used.
C               INDEX COS_ANGLE WEIGHT : lines giving the angle number, its cosine and its weight
C                                        with format I4,1X,2D21.14
C 
C
C    -->Content of the file of external phase functions:
C       Formatted ascii file, created by the user.
C       The format of this file is inspired by that of the DLS tool of Oleg Dubovik (LOA)
C       giving the primary scattering properties of homogeneous
C       randomly oriented spheroid particles.
C       The first lines had to be adapted.
C 
C          1st line : EXTINCTION_COEF : Value
C              Values : Extinction cross section (micron²)
C              Format  : Comment + separator « : » + value (unformatted real)
C                
C          2nd line : SCATTERING_COEF : Value
C              Values : Scattering cross section (micron²)
C              Format  : Comment + separator « : » + value (unformatted real)
C
C	  3rd line : NB_LINES : nbAng
C              Value : Number of angles describing the phase functions.
C              Format : Comment + separator « : » + value (unformatted integer)
C              
C          4th line : Comment line
C          
C          Lines 5 to 5+nbAng-1
C               Values : ANGLE  F11  -F12/F11  F22/F1   F33/F11
C                           ANGLE 	: angle in degrees.
C                           F11 	: phase function in intensity
C 			   -F12/F11 	: ratio of the functions -F12 and F11
C 			   F22/F11 	: ratio of the functions F22 and F11
C 			   F33/F11 	: ratio of the functions F33 and F11
C               ----------------------   ===========================================
C               NB SIGN CONVENTION: F12 is negative for Rayleigh scattering.
C                                        = OS convention
C               ----------------------   ===========================================
C 		 Format  : Real values, unformatted
C 
C          The following lines are not read.
C
C 
C
C    -->Content of the file defining a mixture of aerosol models:
C       Formatted ascii file, created by the user.
C       
C          1st line : NUMBER OF AEROSOLS MODES: Value
C              Values : Number of blocks of mode definition
C              Format : Comment + separator « : » + value (unformatted integer)
C               
C          Per data block of a mode:
C                 Block line 1 : SIZE DISTRIBUTION  : Value
C                 Values : Size distribution type (LND or JUNGE)
C                 Format  : Comment + separator « : » + value (character string)
C
C                 Case of an LND size distribution:
C	              Block line 2 : MODAL RADIUS  : Value
C                    Value : Modal radius of the LND  (in microns)
C              Format  : Comment + separator « : » + value (unformatted real)
C              
C	              Block line 3 : STANDARD DEVIATION  : Value
C                     Value : Variance of the LND (in microns)
C              Format  : Comment + separator « : » + value (unformatted real)
C  
C                 Case of a JUNGE size distribution:
C	               Block line 2 : SLOPE : Value
C                     Value : Slope of the Junge distribution
C              Format  : Comment + separator « : » + value (unformatted real)
C              
C	               Block line 3 : MIN RADIUS : Value
C                     Value : Minimum radius of the Junge distribution (in microns)
C              Format  : Comment + separator « : » + value (unformatted real)
C              
C	               Block line 4 : MAX RADIUS : Value
C                     Value : Maximum radius of the Junge distribution (in microns)
C              Format  : Comment + separator « : » + value (unformatted real)
C    
C                 Block line + 1 : REFRACTIVE INDEX at WA_SIMU - Real part  : Value
C                 Values : Real part of the refractive index at the radiance simulation wavelength
C                 Format  : Comment + separator « : » + value (character string)
C    
C                 Block line + 2 : REFRACTIVE INDEX at WA_SIMU - Imaginary part  : Value
C                 Values : Imaginary part (<0) of the refractive index at the radiance simulation wavelength
C                 Format  : Comment + separator « : » + value (character string)
C    
C                 Block line + 3 : REFRACTIVE INDEX at WA_REF - Real part  : Value
C                 Values : Real part of the refractive index at the AOT reference wavelength
C                 Format  : Comment + separator « : » + value (character string)
C    
C                 Block line + 4 : REFRACTIVE INDEX at WA_REF - Imaginary part  : Value
C                 Values : Imaginary part (<0) of the refractive index at the AOT reference wavelength
C                 Format  : Comment + separator « : » + value (character string)
C    
C                 Block line + 5 : PROPORTION TO THE TOTAL AOT at WA_REF  : Value
C                 Values : Proportion of the mode (between 0 and 1) in the contribution to the AOT at the reference wavelength
C                 Format  : Comment + separator « : » + value (character string)
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
C     Message displayed on screen, program stop and
C     return of status 1 to the shell
C   
C
C  Subroutines used:
C --------------------------
C  This program calls the routines:
C     - SOS_NOM_FICMIE
C     - SOS_MIE
C     - SOS_INIT_PARAMWMO
C     - SOS_INIT_PARAMSF
C     - SOS_INTERPOL
C     - SOS_GRANU
C     - SOS_DECOMPO_LEGENDRE
C     - SOS_INTERPO_SPLINT
C     - SOS_SPLINE
C     - SOS_SPLINT
C
C============================================================================
      SUBROUTINE SOS_AEROSOLS(FICANGLES, 
     &	    	        WA, TA, WAREF, TA_WAREF,	  
     &	  	        ITRONC, IMOD,
     &	  	        RN, IN, IGRANU, 
     &	  	        VARGRANU1_MMD, VARGRANU2_MMD, VARGRANU3_MMD,
     &                IMODELE_WMO, 
     &	  	        C_WMO_DL, C_WMO_WS, C_WMO_OC, C_WMO_SO,
     &	  	        IMODELE_SF, RH, 
     &	  	        MODE_PARAM_BILND, 
     &	  	        USER_CV_COARSE, USER_CV_FINE, RTAUct_WAREF,
     &	  	        BMD_CM_MRWA, BMD_CM_MIWA, 
     &	  	        BMD_CM_MRWAREF, BMD_CM_MIWAREF,
     &	  	        BMD_CM_RMODAL, BMD_CM_VAR,
     &	  	        BMD_FM_MRWA, BMD_FM_MIWA, 
     &	  	        BMD_FM_MRWAREF, BMD_FM_MIWAREF,
     &	  	        BMD_FM_RMODAL, BMD_FM_VAR,	 
     &	  	        FICEXTDATA, FICMIXTURE_AER,
     &			 DIR_TMP, DIRS, FICLOGMIE,
     &	  	        FICGRANU, FICLOGGRANU, KMAT1,
     &			 PIZ, COEF_TRONCA, IER)
 
      
      IMPLICIT NONE

C AEROSOLS constants
C -------------------

      INTEGER CTE_MODALITE	!Max value of the modality parameter
      PARAMETER (CTE_MODALITE=4)

      INTEGER CTE_NB_MODELE	!Number of possible WMO or Shettle & Fenn models.
      PARAMETER (CTE_NB_MODELE=4)

      INTEGER CTE_NB_WMO	    !Number of WMO aerosol types used:
      PARAMETER (CTE_NB_WMO=4)  !DL,WS,OC,SO.
      
      INTEGER CTE_NB_SF	    !Number of aerosol types of the Shettle & Fenn model
      PARAMETER (CTE_NB_SF=5)   !used: SR,LR,SU,LU,OM.
          
      INTEGER CTE_NB_AERO	    !Initialisation with the largest number of aerosol
      PARAMETER (CTE_NB_AERO=5) !types (Shettle & Fenn case).
     

C* Definition of variables
C*-----------------------------------------------------------------


      DOUBLE PRECISION XMU(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)  
                                                  !Cosines of the Gauss (and user) angles.
      DOUBLE PRECISION XHR(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)  
                                                  !Weights of the Gauss angles.

      DOUBLE PRECISION TA	! Aerosol optical thickness (not truncated).

      DOUBLE PRECISION WA	! Selected wavelength in microns.
            						  						  
      DOUBLE PRECISION C(CTE_NB_WMO)	 !Volume proportion of the WMO components
      DOUBLE PRECISION C_WMO_DL	    	 ! --> User value for the DL particles
      DOUBLE PRECISION C_WMO_WS	    	 ! --> User value for the WS particles
      DOUBLE PRECISION C_WMO_OC	    	 ! --> User value for the OC particles
      DOUBLE PRECISION C_WMO_SO	    	 ! --> User value for the SO particles
      DOUBLE PRECISION V(CTE_NB_WMO)	 !Volume concentration of the WMO components
      DOUBLE PRECISION N(CTE_NB_WMO)	 ! C/V
      DOUBLE PRECISION NTOT		        ! Sum of the weighting coefficients
                                           ! of the aerosol components.
      DOUBLE PRECISION ALPHAF_PARTICULE(4) ! Alpha Max value for each
      					        ! WMO component (DL, WS, OC, SO)
					   
      DOUBLE PRECISION Ni(CTE_NB_SF)	 !particle density (Shettle & Fenn)
      DOUBLE PRECISION RH                  ! % of relative humidity of the Shettle & Fenn models



 	!   Mono-modal model
	!------------------- 
 	! Refractive index for the mono-modal model
      DOUBLE PRECISION RN	! real part of the refractive index
      DOUBLE PRECISION IN	! imaginary part of the refractive index
 	! Size distribution parameters for the mono-modal model
      DOUBLE PRECISION VARGRANU1_MMD    ! Size distribution parameter 1.
      DOUBLE PRECISION VARGRANU2_MMD    ! Size distribution parameter 2.
      DOUBLE PRECISION VARGRANU3_MMD    ! Size distribution parameter 3.



      DOUBLE PRECISION ALPHAO   ! minimum alpha value of the Mie computations
      DOUBLE PRECISION ALPHAF   ! maximum alpha value of the Mie computations

      DOUBLE PRECISION SOMME_NR	! number of particles: integration over the radii of NR * PR.
      
      DOUBLE PRECISION MR(CTE_NB_AERO)   ! real part of the refractive index
      					      ! for each aerosol component
      DOUBLE PRECISION MI(CTE_NB_AERO)   ! imaginary part of the refractive index.
      					      ! for each aerosol component.
 



      DOUBLE PRECISION VARGRANU1I(CTE_NB_AERO)   ! Size distribution parameter 1
      						       ! for each component model.
      DOUBLE PRECISION VARGRANU2I(CTE_NB_AERO)   ! Size distribution parameter 2
      						       ! for each component model.
	      				

	!Parameters for a bimodal LND model
      INTEGER*2 MODE_PARAM_BILND   !Index of the choice of description of the mixture of 2 LND.
      DOUBLE PRECISION CVI(2)	!Volume concentration of the "LND coarse" (1) and
      				      	!of the "LND fine" (2)
      DOUBLE PRECISION USER_CV_COARSE	!User value of the volume concentration of the "LND coarse"
      DOUBLE PRECISION USER_CV_FINE	!User value of the volume concentration of the "LND fine"
      DOUBLE PRECISION RTAUct_WAREF	!Ratio of the optical thicknesses AOT_coarse / AOT_total
      DOUBLE PRECISION WAREF		!Wavelength of the ratio AOT_coarse / AOT_total
      DOUBLE PRECISION TA_WAREF		!Total aerosol optical thickness
      					       !at the reference wavelength.
					
      DOUBLE PRECISION RN_LND(2)	!Real part of the refractive index of the
      					!"LND coarse" (1) and of the "LND fine" (2) at the wavelength WA.
      DOUBLE PRECISION BMD_CM_MRWA	! Value read for RN_LND(1) : coarse mode
      DOUBLE PRECISION BMD_FM_MRWA	! Value read for RN_LND(2) : fine mode
       
      DOUBLE PRECISION IN_LND(2) 	!Imaginary part of the refractive index of the
      					!"LND coarse" (1) and of the "LND fine" (2) at the wavelength WA.
      DOUBLE PRECISION BMD_CM_MIWA	! Value read for IN_LND(1) : coarse mode
      DOUBLE PRECISION BMD_FM_MIWA	! Value read for IN_LND(2) : fine mode
      
      DOUBLE PRECISION RN_LND_WAREF(2)	!Real part of the refractive index of the
      					       !"LND coarse" (1) and of the "LND fine" (2) at the wavelength WAREF.
      DOUBLE PRECISION BMD_CM_MRWAREF	! Value read for RN_LND_WAREF(1) : coarse mode
      DOUBLE PRECISION BMD_CM_MIWAREF	! Value read for RN_LND_WAREF(2) : fine mode
      
      DOUBLE PRECISION IN_LND_WAREF(2) 	!Imaginary part of the refractive index of the
      					       !"LND coarse" (1) and of the "LND fine" (2) at the wavelength WAREF.
      DOUBLE PRECISION BMD_FM_MRWAREF	! Value read for IN_LND_WAREF(1) : coarse mode
      DOUBLE PRECISION BMD_FM_MIWAREF	! Value read for IN_LND_WAREF(2) : fine mode
      
      DOUBLE PRECISION RMODAL_LND(2)  	!Modal radius of the "LND coarse" (1) and "fine" (2).
      DOUBLE PRECISION BMD_CM_RMODAL	! Value read for RMODAL_LND(1) : coarse mode
      DOUBLE PRECISION BMD_FM_RMODAL	! Value read for RMODAL_LND(2) : fine mode
      
      DOUBLE PRECISION VAR_LND(2)	!Variance of the "LND coarse" (1) and "fine" (2).
      DOUBLE PRECISION BMD_CM_VAR	! Value read for VAR_LND(1) : coarse mode
      DOUBLE PRECISION BMD_FM_VAR	! Value read for VAR_LND(2) : fine mode
          
      
      DOUBLE PRECISION KMAT1c      ! Extinction cross section (micron^2) for the
      					! "coarse" mode at the wavelength WAREF. Not truncated.
      DOUBLE PRECISION KMAT1f	! Extinction cross section (micron^2) for the
      					! "fine" mode at the wavelength WAREF. Not truncated.
   

  
      DOUBLE PRECISION KMAT1I	! Extinction cross section (micron^2)
      				       ! for a component I. Not truncated.
      DOUBLE PRECISION KMAT2I 	! Scattering cross section (micron^2)
      				       ! for a component I. Not truncated.
      DOUBLE PRECISION TAUEXTI     ! Extinction optical thicknesses of an aerosol mode I
      DOUBLE PRECISION TAUSCAI     ! Scattering optical thicknesses of an aerosol mode I

      DOUBLE PRECISION KMAT1	! Extinction cross section (micron^2). Not truncated.
      DOUBLE PRECISION KMAT2 	! Scattering cross section (micron^2). Not truncated.
      DOUBLE PRECISION COEF_TRONCA ! Truncation coefficient.
      DOUBLE PRECISION PIZ	       ! Single scattering albedo.
      DOUBLE PRECISION PIZTR	! Single scattering albedo after truncation.
      
      
      DOUBLE PRECISION P11I(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Phase function computed from the functions of
		         !the MIE theory for the component I.
      DOUBLE PRECISION P12I(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Element P12 of the phase matrix for the component I.
      DOUBLE PRECISION P33I(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Element P33 of the phase matrix for the component I.
      
      DOUBLE PRECISION P11(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Phase function of the simulated model.
      DOUBLE PRECISION P12(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Element P12 of the phase matrix.
      DOUBLE PRECISION P22(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Element P22 of the phase matrix.
      DOUBLE PRECISION P33(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Element P33 of the phase matrix.
		       
      DOUBLE PRECISION TTT(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Copy of the phase function (P11) before truncation.
     
      DOUBLE PRECISION Z1	!Value of BETA11(0) before normalisation.
      DOUBLE PRECISION TETA !Scattering angle.
      
                  
      DOUBLE PRECISION ALP(0:CTE_OS_NB_MAX)	! Parameters alpha(k).
      DOUBLE PRECISION BETA11(0:CTE_OS_NB_MAX)	! Parameters beta11(k).
      DOUBLE PRECISION BETA22(0:CTE_OS_NB_MAX)	! Parameters beta22(k).
      DOUBLE PRECISION GAMMA12(0:CTE_OS_NB_MAX)	! Parameters gamma12(k).
      DOUBLE PRECISION ZETA(0:CTE_OS_NB_MAX)	! Parameters zeta(k).
      DOUBLE PRECISION DELTA33(0:CTE_OS_NB_MAX)	! Parameters delta33(k).
     
	    
      DOUBLE PRECISION X, RMAX   ! Parameters for the computation of AlphaMax
      				     ! according to the particle size.
               
      DOUBLE PRECISION VAL	! Variable to read unused values
      
      DOUBLE PRECISION ANG_EXT(CTE_MAXNB_ANG_EXT)     ! Table of the angles of the external phase functions
      DOUBLE PRECISION MU_EXT(CTE_MAXNB_ANG_EXT)	     ! Table of the cosines of the angles of the external phase functions
      DOUBLE PRECISION F11_EXT(CTE_MAXNB_ANG_EXT)     ! Phase function P11 of external data
      DOUBLE PRECISION F12_EXT(CTE_MAXNB_ANG_EXT)     ! Phase function P11 of external data
      DOUBLE PRECISION F22_EXT(CTE_MAXNB_ANG_EXT)     ! Phase function P11 of external data
      DOUBLE PRECISION F33_EXT(CTE_MAXNB_ANG_EXT)     ! Phase function P11 of external data
      DOUBLE PRECISION F12sF11_EXT(CTE_MAXNB_ANG_EXT) ! Phase function P12/P11 of external data
      DOUBLE PRECISION F22sF11_EXT(CTE_MAXNB_ANG_EXT) ! Phase function P22/P11 of external data
      DOUBLE PRECISION F33sF11_EXT(CTE_MAXNB_ANG_EXT) ! Phase function P33/P11 of external data
      
      INTEGER*2 IMOD	       ! Index of the aerosol modelling case.
      INTEGER*2 IMODELE_WMO ! Type of the WMO model.
      INTEGER*2 IMODELE_SF  ! Type of the Shettle & Fenn model.

      INTEGER*4 NBMU_GAUSS  ! Number of Gauss angles used
      INTEGER*4 MIE_NBMU   ! Number of actual angles of the simulations
      INTEGER*4 OS_NB   ! Order of the expansions of the phase functions
			   ! in Legendre polynomials.
      INTEGER*4 IGRANU	   ! Number of the size distribution model
      INTEGER*4 ITRONC	   ! Truncation option (1 to apply it).
      INTEGER*4 LDIRTMP        ! length of the character string DIRTMP
      INTEGER*4 LDIRS	   ! length of the character string DIRS
      INTEGER*4 LFICS	   ! length of the character string FICS
      INTEGER*4 LFICTMP	   ! length of the character string FICMIE_TMP
      INTEGER*4 LFIC           ! length of a character string
      INTEGER*4 IFIN           ! nb of significant characters of a string
      INTEGER*4 LCOM           ! length of the character string COMMAND
      INTEGER*4 I	      ! Array index
      INTEGER*4 J	      ! Array index
      INTEGER*4 K          ! Array index
      INTEGER*4 IER        ! error code =0 if no error, =-1 otherwise
      INTEGER*4 IARGC      ! Type of the function IARGC
      INTEGER*4 NPARAM     ! Number of the input parameter (type INTEGER*4 mandatory for f77)
      INTEGER*4 NB_FICEXT_ANGLES ! Number of angles of the external phase function file
  
      INTEGER*4 NB_MODES_MELANGE   ! Number of aerosol modes for a mixture defined by the user
      INTEGER*4 NMOD			! Number of an aerosol mode
    
      CHARACTER*CTE_LENDIR DIR_TMP        ! access path to the MIE files being generated
      CHARACTER*CTE_LENDIR DIRS	       ! access path to the MIE file
      CHARACTER*CTE_LENFIC1 FICS		! name of the MIE file
      CHARACTER*CTE_LENFIC2 FICMIE	       ! path and name of the MIE file
      CHARACTER*CTE_LENFIC2 FICMIE_TMP    ! path and name of the MIE file being generated
      CHARACTER*CTE_LENCOM  COMMAND       ! Command string for system call.
      CHARACTER*CTE_LENFIC2 FICLOGMIE           ! name of the MIE TRACE file
      CHARACTER*CTE_LENFIC2 FICLOGGRANU         ! name of the GRANU TRACE file
      CHARACTER*CTE_LENFIC2 FICGRANU	      ! name of the GRANU result file
      CHARACTER*CTE_LENFIC2 FICANGLES           ! path and name of the file of the actual angles
      CHARACTER*CTE_LENFIC2 FICANGLES_USER      ! path and name of the "user" angles file
      CHARACTER*CTE_LENFIC2 FICEXTDATA          ! path and name of the external data file of phase
      						      ! functions and radiative parameters (extinction and
						      ! scattering coefficients)
 
      CHARACTER*CTE_LENFIC2 FICMIXTURE_AER      ! path and name of the file defining a mixture of aerosol
                                                ! models by the user
     
      CHARACTER*CTE_LENDIR DIRFIC	             ! access path to the file of the parameters of the WMO and S&F models
      CHARACTER*CTE_LENFIC2 FICDATAWMO          ! path and name of the file Data_WMO
      
      CHARACTER*15 PARTICULE(CTE_NB_AERO)       ! name of the different WMO particles
      					             ! or of the S&F models

      CHARACTER*LENLINE LIGNE_TEXTE             ! Character line for decoding parameters
      CHARACTER*1000 CVAL	                    ! Character string to read unused values

      CHARACTER*10 CMODEL       ! Character string giving the name of an aerosol mode (LND or JUNGE)
 
      INTEGER*4 MEL_MODEL_TYPE(CTE_MAX_NB_MODE_MIXTURE)           ! Array listing the model types of a mixture (coef IGRANU)
      DOUBLE PRECISION MEL_MODEL_PARAM(CTE_MAX_NB_MODE_MIXTURE,3) ! Array listing the parameters of the models of a mixture
      DOUBLE PRECISION MEL_MODEL_INDEX_WAREF(CTE_MAX_NB_MODE_MIXTURE,2)  ! Array listing the refractive indices of the models of a mixture at WAREF
      DOUBLE PRECISION MEL_MODEL_INDEX_WASIMU(CTE_MAX_NB_MODE_MIXTURE,2) ! Array listing the refractive indices of the models of a mixture at WASIMU
      DOUBLE PRECISION MEL_MODEL_AOT(CTE_MAX_NB_MODE_MIXTURE)     ! Array listing the AOT per model of a mixture
      DOUBLE PRECISION AOT_TAUX					   ! Proportion of the total AOT of the model of a mixture
                                                                  ! at the reference wavelength.
      DOUBLE PRECISION SOM_TAUX_AOT				   ! Sum of the rates of the different modes
      DOUBLE PRECISION COEF_ALPHA(CTE_MAX_NB_MODE_MIXTURE)        ! Array listing the weighting coefficients alpha per model
      DOUBLE PRECISION SOM_COEF_ALPHA                             ! Sum of the coefficients before normalisation
      DOUBLE PRECISION RN_WAREF_LU, IN_WAREF_LU			   ! Refractive index (real part and imaginary part)
                                                                  ! read in the user file defining a mixture
                                                                  ! at the reference wavelength.
      DOUBLE PRECISION RN_WASIMU_LU, IN_WASIMU_LU		   ! Same at the radiance simulation wavelength.

      LOGICAL TRACE     ! = true if writing in the trace file
      LOGICAL EX	   ! Checks the existence of a file


C*-------------
C* Processing
C*-------------    
 
 
      
C* Initialisation
C*-----------------------------------------------------------------
      IER=0

C --- Naming of the temporary MIE file
C     for its creation in the work space before moving it
C     to the archive space.
      LDIRTMP = INDEX(DIR_TMP,'/TMP') +3
      FICMIE_TMP=DIR_TMP(1:LDIRTMP)//'/'//"MIE_FILE_BUILDING_TMP"   
      LFICTMP = INDEX(FICMIE_TMP,' ') - 1
	
C --- Assignment of the minimum value of the size parameter
C     (constant converted to double precision)
      ALPHAO = CTE_MIE_ALPHAMIN      
      
      IF ((ALPHAO.LT.1D-4).OR.(ALPHAO.GE.10D-00)) GOTO 1010
     	
C* Opening of the trace file
C------------------------------------------------------
      IF (FICLOGGRANU.EQ.'NO_LOG_FILE') THEN
	 TRACE=.FALSE.
      ELSE
        TRACE=.TRUE.
        OPEN (88,FILE=FICLOGGRANU,ERR=911)
        
        WRITE (88,*,ERR=921) ""
	 WRITE (88,*,ERR=921) ""
	 WRITE (88,*,ERR=921) "****************************"	 
	 WRITE (88,*,ERR=921) "SOS_AEROSOLS :  "
	 WRITE (88,*,ERR=921) "****************************"
	 WRITE (88,*,ERR=921) ""
	 WRITE (88,*,ERR=921) "Reference wavelength :", WAREF
	 WRITE (88,*,ERR=921) "Wavelength for radiance simulation:", WA
	 WRITE (88,*,ERR=921) ""
	 
      ENDIF
      
C* Reading of the angle and weight tables used by the simulations
C* -----------------------------------------------------------------     
      OPEN(UNIT=10,FILE=FICANGLES,STATUS='OLD',ERR=930)

      READ(10,'(a)',ERR=931) LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) MIE_NBMU
      
      READ(10,'(a)',ERR=931)  LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) NBMU_GAUSS
      
      READ(10,'(a)',ERR=931)  LIGNE_TEXTE
      !Only the part after the separator ':' is kept
      IFIN=INDEX(LIGNE_TEXTE,':')	
      LIGNE_TEXTE=LIGNE_TEXTE(IFIN+1:LENLINE)	
      !Search for the end of the useful line ' ' to reduce the string size
      IFIN=INDEX(LIGNE_TEXTE,' ')
      FICANGLES_USER=LIGNE_TEXTE(1:IFIN)
      
      READ(10,'(a)',ERR=931)  LIGNE_TEXTE
      IFIN=INDEX(LIGNE_TEXTE,':')
      READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) OS_NB
      
      READ(10,*,ERR=931) CVAL
      DO J=1,MIE_NBMU
	  READ(10,*,ERR=931) CVAL, XMU(J), XHR(J)
      ENDDO 
      
      CLOSE(10)
      
      DO J=1,MIE_NBMU
         XMU(-J) = -XMU(J)
	 XHR(-J) =  XHR(J)
      ENDDO 
                
        
      IF (TRACE) THEN
        WRITE (88,*,ERR=921) "USED ANGLES"
	 WRITE (88,*,ERR=921) "****************"
	 WRITE (88,*,ERR=921) " "
	 WRITE (88,*,ERR=921) "Total number of angles : ", MIE_NBMU
	 WRITE (88,*,ERR=921) "Number of Gauss angles : ", NBMU_GAUSS
	 WRITE (88,*,ERR=921) " "
        IF (FICANGLES_USER.NE.'NO_USER_ANGLES') THEN
	    WRITE(88,*,ERR=921)  
     &	       "A file of user's angles is used :",
     &         FICANGLES_USER
     	    WRITE(88,*,ERR=921) "Number of user's angles :",
     &         (MIE_NBMU-NBMU_GAUSS)
	 ELSE
	    WRITE(88,*,ERR=921)  "No file of user's angles"
	 ENDIF
        WRITE (88,*,ERR=921) " " 
	 WRITE (88,*,ERR=921) "Max order NB : ", OS_NB
	 WRITE (88,*,ERR=921) " " 
	 WRITE (88,*,ERR=921) " " 
	 WRITE (88,*,err=921) "List of used angles :"
	 WRITE (88,*,err=921) "  Index, mu, angle (degrees), weigth"
	 DO 1000 J=-MIE_NBMU,MIE_NBMU
	    IF (J.EQ.0) GOTO 1000
	    IF (J.LT.0) THEN 
	        WRITE(88,222,ERR=921) J, XMU(J),
     &	                     -DACOS(-XMU(J))*180./INCTE_PI, XHR(J)
           ELSE
	        WRITE(88,222,ERR=921) J, XMU(J),
     &	                      DACOS(XMU(J))*180./INCTE_PI, XHR(J)
	    ENDIF	    
1000    CONTINUE
	 WRITE (88,*,ERR=921) " " 
	 WRITE (88,*,ERR=921) " "  	
      ENDIF !Trace file


C* Retrieval of the name of the directory containing the WMO
C* and Shettle & Fenn parameter files
C----------------------------------------------------------------------------
      CALL GETENV('SOS_ABS_ROOT',DIRFIC)
      IFIN=INDEX(DIRFIC,' ')
      IFIN=IFIN-1
      IF (IFIN.LE.0) GOTO 925
      
      IF ((IFIN+4).GT.CTE_LENDIR) GOTO 926
      
      DIRFIC(IFIN+1:IFIN+4)='/fic'
      IFIN=IFIN+4
      IF (IFIN.LE.0) IFIN=CTE_LENDIR
            
C* Initialisation of the proportions of the components of the WMO and Shettle & Fenn models
C----------------------------------------------------------------------------------      
      DO I=1,CTE_NB_WMO
	  C(I)=0.
	  N(I)=0.
      ENDDO
      
       DO I=1,CTE_NB_SF
	  Ni(I)=0.
      ENDDO

      NTOT=0.
            
C* Initialisation of the arrays of the phase matrix parameters
C--------------------------------------------------------------------
      DO 101 K=0,OS_NB
         BETA11(K)=0.D+00
	  BETA22(K)=0.D+00
         GAMMA12(K)=0.D+00
         DELTA33(K)=0.D+00
         ALP(K)=0.D+00
         ZETA(K)=0.D+00
  101 CONTINUE

C* Initialisation of the phase function arrays
C--------------------------------------------------
      DO 102 J=-MIE_NBMU,MIE_NBMU
         P11(J)=0.
         P12(J)=0.
         P22(J)=0.
	  P33(J)=0.
         TTT(J)=0.
  102 CONTINUE

C* Initialisation of the extinction and scattering cross sections
C--------------------------------------------------------------------
      KMAT1=0.
      KMAT2=0.
                        
C* Computation of the parameters of the aerosol phase matrix
C* and of the single scattering albedo.
C------------------------------------------------------------
      IF (TA.EQ.0.D+00) THEN	!No computation on the size distribution
      
      	PIZ = 0.D+00
	GOTO 1111	!Recording of the GRANU result file
      
      ENDIF

      !Processing goes on for TA > 0
      !=====================================
      
      	
C* Length of the name of the directory of the Mie files
C---------------------------------------------------
      LDIRS = INDEX(DIRS,' ')
      LDIRS = LDIRS - 1
      IF(LDIRS.LE.0) LDIRS = 100
      

C ===================
C* Mono-modal computation
C ===================

      IF(IMOD.EQ.0) THEN
	 	 	
      
C* Values of RN and IN forced to the formats F5.3 and F8.5
C  -----------------------------------------
	 RN = DNINT(RN*1000.D+00) / 1000.D+00
	 IN = -DNINT(-IN*100000.D+00) / 100000.D+00 	 


        IF (IGRANU.EQ.1) THEN	! LND case
	 
	     !!!!!IGRANU_MMD_CALL = 1	! LND index for routine OSOAA_GRANU
	     
	     X = VARGRANU2_MMD   
	     RMAX = VARGRANU1_MMD * DEXP(X*X) *
     &	            DEXP(X*DSQRT(-2.D+00*LOG(CTE_COEF_NRMAX)))
	     
	 ELSE	! Junge law case for aerosols
	 
	     !!!!!IGRANU_MMD_CALL = 3	! Junge index for aerosols for OSOAA_GRANU

	     RMAX = VARGRANU3_MMD
	     
	 ENDIF ! End of case IGRANU_MMD = 2


	 ALPHAF 
     &	 =REAL(100+100*DINT(2.*INCTE_PI*RMAX / (100.*CTE_WAMIN)))

	 IF ((ALPHAO.GT.ALPHAF).OR.(ALPHAF.GE.1D5)) GOTO 1009 


        IF (TRACE) THEN
           WRITE(88,*,ERR=921) "MONO-MODAL MODEL"
	    WRITE(88,*,ERR=921) "****************"
	    WRITE(88,*,ERR=921) " "
	    WRITE(88,*,ERR=921) "Wavelength in microns : ",WA
	    WRITE(88,*,ERR=921) " "
	    WRITE(88,*,ERR=921) 
     & 	 "Refractive index (real part) : ",RN
	    WRITE(88,*,ERR=921) 
     & 	 "Refractive index (imaginary part) : ",IN
	    IF(IGRANU.EQ.1) THEN
               WRITE(88,*,ERR=921) "Size distribution : Log-normal"
	        WRITE(88,*,ERR=921) "  Modal radius (microns) : " , 
     &	        VARGRANU1_MMD 
	        WRITE(88,*,ERR=921) "  Variance (sig) : ", VARGRANU2_MMD
	        WRITE(88,*,ERR=921) "  sig / ln10  : " , 
     &         VARGRANU2_MMD  / DLOG(10.D+00)
               WRITE(88,*,ERR=921) 
     &	        "  Calculated Rmax limit (mic) : ", RMAX
	    ELSE
	       WRITE(88,*,ERR=921) "Size distribution : Junge's law"
              WRITE(88,*,ERR=921) "  Slope : " , VARGRANU2_MMD 
              WRITE(88,*,ERR=921) "  Minimal radius (microns) : ", 
     &	       VARGRANU1_MMD
              WRITE(88,*,ERR=921) "  Maximal radius (microns) : ", RMAX
	    ENDIF   

           WRITE(88,*,ERR=921) " "
           WRITE(88,*,ERR=921) "  AlphaF limit : ", ALPHAF
	    WRITE(88,*,ERR=921) " "	
	    WRITE(88,*,ERR=921) " "  		
         ENDIF !Trace file
 

     	 
C* Building of the name of the MIE file
C  ----------------------------------

	 CALL SOS_NOM_FICMIE(NBMU_GAUSS,FICANGLES_USER,
     &                      RN,IN,ALPHAO,ALPHAF,FICS) 	     

	    
	 FICMIE = DIRS(1:LDIRS)//'/'//FICS	 
	 LFICS=INDEX(FICMIE,' ')-1

	 IF (TRACE) THEN 
	    WRITE (88,*,ERR=921) "MIE file : ", 
     &	    FICMIE(1:LFICS)	
	 ENDIF	
	 
	 !Screen display
	 WRITE(6,*) "  MIE file : ", FICMIE(1:LFICS)

	 
C* If the MIE file already exists, the program is not run
C  -----------------------------------------------------------	 
        INQUIRE(FILE=FICMIE,EXIST=EX)
	 
        IF (EX) THEN  
            WRITE(6,*)'  ==> The file has already been calculated'
	  
	 ELSE   
      
C* Computation of the MIE file
C  ------------------------ 
            CALL SOS_MIE(MIE_NBMU,XMU,XHR,RN,IN,ALPHAO,ALPHAF,
     &	                 FICMIE_TMP,FICLOGMIE,IER)
            IF (IER.NE.0) GOTO 993   ! error case

            COMMAND='mv '//FICMIE_TMP(1:LFICTMP)//' '//FICMIE(1:LFICS)
            LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   !The search starts after the last blank at the beginning of the string
            LCOM = LCOM + (4+LFICTMP)          
            IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020

            INQUIRE(FILE=FICMIE,EXIST=EX)
            IF (EX) THEN !If the MIE file exists, the temporary file (generated by another process) is deleted
                OPEN(UNIT=3,FILE=FICMIE_TMP,FORM='UNFORMATTED')
                CLOSE(3,STATUS='DELETE') 
            ELSE !Otherwise it is moved to the storage space of the MIE files
               CALL SYSTEM(COMMAND)  
            ENDIF

        ENDIF
	 
C* Computation of the extinction and scattering cross sections
C  and of the phase functions I, Q, U
C  ----------------------------------------------------------          
        CALL SOS_GRANU(FICMIE,IGRANU,VARGRANU1_MMD,
     &                 VARGRANU2_MMD,VARGRANU3_MMD,
     &                 WA,MIE_NBMU,XMU,TRACE,
     &                 KMAT1,KMAT2,SOMME_NR,
     &                 P11,P12,P33,IER)
        IF (IER.NE.0) GOTO 995   ! error case

C* Case of spherical particles: P22(J) = P11(J)
C  -----------------------------------------------------
        DO 1021 J=-MIE_NBMU,MIE_NBMU
           P22(J)=P11(J)
 1021   CONTINUE	 
  
C* Computation of the parameters alpha, beta, gamma, delta, zeta
C  with or without truncation
C  -----------------------------------------------------
	 CALL SOS_DECOMPO_LEGENDRE(ITRONC,TRACE,MIE_NBMU,XMU,XHR,OS_NB,
     &				      P11,TTT,P12,P22,P33,COEF_TRONCA,Z1,
     &				      ALP,BETA11,BETA22,
     &                            GAMMA12,DELTA33,ZETA,IER)  
        IF (IER.NE.0) GOTO 996   ! error case

      ENDIF	!End of test on modality IMOD=0 (Mono-modal)



C ====================================================
C* Multi-modal computation with the WMO components
C ====================================================
	    
      IF(IMOD.EQ.1) THEN

	 
C* Initialisation
C  --------------
        ! Name of the different WMO particles
	 PARTICULE(1)='DUST LIKE'
        PARTICULE(2)='WATER SOLUBLE'
        PARTICULE(3)='OCEANIC'
        PARTICULE(4)='SOOT'

        ALPHAF_PARTICULE(1) = CTE_ALPHAMAX_WMO_DL
	 ALPHAF_PARTICULE(2) = CTE_ALPHAMAX_WMO_WS
	 ALPHAF_PARTICULE(3) = CTE_ALPHAMAX_WMO_OC
	 ALPHAF_PARTICULE(4) = CTE_ALPHAMAX_WMO_SO
	 
	 IGRANU =   1	! Log Normal distribution by default
	    
C* Initialisation of the WMO parameters
C  ---------------------------------------
	 FICDATAWMO=DIRFIC(1:IFIN)//'/'//CTE_AER_DATAWMO

	 CALL SOS_INIT_PARAMWMO(FICDATAWMO,WA,VARGRANU1I,
     &                          VARGRANU2I,MR,MI,V,IER)
        IF (IER.NE.0) GOTO 994   ! error case

C* Definition of the volume proportion of each component
C  for the selected model
C  -----------------------------------------------------------	 
        IF(IMODELE_WMO.EQ.1) THEN	!Continental model
	    C(1) = 0.70
	    C(2) = 0.29
	    C(4) = 0.01      
	    
        ELSE IF(IMODELE_WMO.EQ.2) THEN 	!Maritime model
	    C(2) = 0.05
	    C(3) = 0.95
	 
        ELSE IF(IMODELE_WMO.EQ.3) THEN	!Urban model
	    C(1) = 0.17
	    C(2) = 0.61
	    C(4) = 0.22
	       
        ELSE IF(IMODELE_WMO.EQ.4) THEN	!Model defined by the user
	 
C*          Check of the retrieval of the volume proportions
C*          defined by the user
C           ------------------------------------------------------
	    
	    C(1) = C_WMO_DL    !Volume percentage of the DL particles
	    C(2) = C_WMO_WS    !Volume percentage of the WS particles
	    C(3) = C_WMO_OC    !Volume percentage of the OC particles
	    C(4) = C_WMO_SO    !Volume percentage of the SO particles

        ENDIF	!End of the definition of the models
	    
        DO 103 I=1,CTE_NB_WMO	    
            N(I) = C(I)/V(I)
	     NTOT = NTOT + N(I)
  103   CONTINUE


        IF (TRACE) THEN
          WRITE (88,*,ERR=921) "WMO MODEL"
	   WRITE (88,*,ERR=921) "*********"
	   WRITE (88,*,ERR=921) " "
	   IF(IMODELE_WMO.EQ.1) 
     &	      WRITE(88,*,ERR=921) "CONTINENTAL Model"
	   IF(IMODELE_WMO.EQ.2) WRITE(88,*,ERR=921) "MARITIME Model"
	   IF(IMODELE_WMO.EQ.3) WRITE(88,*,ERR=921) "URBAN Model"
	   WRITE (88,*,ERR=921) " "
	   WRITE (88,*,ERR=921) "  Coefficients alphaK of the mixture:"
	   DO I=1,CTE_NB_WMO
	      WRITE (88,*,ERR=921) "    ",PARTICULE(I)," : ", N(I)/NTOT
          ENDDO
	   WRITE (88,*,ERR=921) " "
	   IF(IMODELE_WMO.EQ.4) THEN
	     WRITE (88,*,ERR=921) "USER Model "
	     DO I=1,4
	        WRITE (88,*,ERR=921) "    ",PARTICULE(I)," : ",
     &		                     100*C(I)," % volume"
            ENDDO
	   ENDIF 
	   WRITE (88,*,err=921) "Wavelength in microns : ",WA
	   WRITE (88,*,ERR=921) " " 
	   WRITE (88,*,ERR=921) " "  		
        ENDIF !Trace file
      
      

C* For each WMO aerosol type
C  -------------------------------------	       
        DO 105 I=1,CTE_NB_WMO
	    
	    IF(C(I).EQ.0.) GOTO 105
		
C* Building of the name of the MIE file
C  ----------------------------------
	    ALPHAF = ALPHAF_PARTICULE(I)

           CALL SOS_NOM_FICMIE(NBMU_GAUSS,FICANGLES_USER,
     &	                        MR(I),MI(I),ALPHAO,ALPHAF,FICS)       
	    FICMIE = DIRS(1:LDIRS)//'/'//FICS
           LFICS=INDEX(FICMIE,' ')-1
	    
	    IF (TRACE) THEN 
	       WRITE (88,*,ERR=921) "Fichier MIE : ", 
     &	       FICMIE(1:LFICS)	
	    ENDIF	
	 
	    !Screen display
	    WRITE(6,*) "  MIE file : ", FICMIE(1:LFICS)
	    
C* If the MIE file already exists, the program is not run
C  -----------------------------------------------------------		              
           INQUIRE(FILE=FICMIE,EXIST=EX)
	    
	    IF (EX) THEN  
               WRITE(6,*)'  ==> The file has already been calculated'
	  
	    ELSE   
      
C* Computation of the MIE file
C  ------------------------ 
              CALL SOS_MIE(MIE_NBMU,XMU,XHR,MR(I),MI(I),ALPHAO,ALPHAF,
     &                      FICMIE_TMP,FICLOGMIE,IER)
    	       IF (IER.NE.0) GOTO 993   ! error case
               
              COMMAND=
     &        'mv '//FICMIE_TMP(1:LFICTMP)//' '//FICMIE(1:LFICS)
              LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   !The search starts after the last blank at the beginning of the string
              LCOM = LCOM + (4+LFICTMP)
              IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020
              
              INQUIRE(FILE=FICMIE,EXIST=EX)
              IF (EX) THEN !If the MIE file exists, the temporary file (generated by another process) is deleted
                 OPEN(UNIT=3,FILE=FICMIE_TMP,FORM='UNFORMATTED')
                 CLOSE(3,STATUS='DELETE') 
              ELSE !Otherwise it is moved to the storage space of the MIE files
                 CALL SYSTEM(COMMAND)  
              ENDIF
   
           ENDIF
	 
C* Computation of the extinction and scattering cross sections
C  and of the phase functions I, Q, U of the component I
C  ---------------------------------------------------------- 
	    IF (TRACE) THEN
	       WRITE(88,*,err=921) PARTICULE(I)
	       WRITE(88,*,err=921)'============='
	       WRITE(88,*,err=921) ' '
	    ENDIF
	       
           CALL SOS_GRANU(FICMIE,IGRANU,VARGRANU1I(I),
     &                    VARGRANU2I(I),CTE_NOT_DEFINED_VALUE_DBLE,
     &                    WA,MIE_NBMU,XMU,TRACE,
     &			     KMAT1I,KMAT2I,SOMME_NR,P11I,P12I,P33I,IER)
           IF (IER.NE.0) GOTO 995   ! error case
	    

C* Computation of the extinction and scattering cross sections of the mixture
C  ---------------------------------------------------------------------	     
	    KMAT1 = KMAT1 + (N(I)/NTOT)*KMAT1I
	    KMAT2 = KMAT2 + (N(I)/NTOT)*KMAT2I
	    
	    DO 104 J=-MIE_NBMU,MIE_NBMU	     
	       P11(J) = P11(J) + (N(I)/NTOT)*P11I(J)*KMAT2I
	       P12(J) = P12(J) + (N(I)/NTOT)*P12I(J)*KMAT2I
	       P33(J) = P33(J) + (N(I)/NTOT)*P33I(J)*KMAT2I
  104	    CONTINUE  !End of loop on the Gauss angles
    
  105	 CONTINUE	 !End of loop on the WMO components


C* Normalisation of the phase functions of the mixture
C  -----------------------------------------------	    
	 DO 106 J=-MIE_NBMU,MIE_NBMU   
	    P11(J) = P11(J) / KMAT2
	    P12(J) = P12(J) / KMAT2
	    P33(J) = P33(J) / KMAT2
  106	 CONTINUE

C* Case of spherical particles: P22(J) = P11(J)
C  -----------------------------------------------------
        DO 1022 J=-MIE_NBMU,MIE_NBMU
           P22(J)=P11(J)
 1022   CONTINUE	 
 
C* Computation of the parameters alpha, beta, gamma, delta, zeta
C  with or without truncation
C  -----------------------------------------------------  
	 CALL SOS_DECOMPO_LEGENDRE(ITRONC,TRACE,MIE_NBMU,XMU,XHR,OS_NB,
     &				      P11,TTT,P12,P22,P33,COEF_TRONCA,Z1,
     &				      ALP,BETA11,BETA22,
     &                            GAMMA12,DELTA33,ZETA,IER)  
        IF (IER.NE.0) GOTO 996   ! error case


      ENDIF	!End of test on modality IMOD=1 (WMO)
  
       
C =====================================================================
C* Bi-modal computation with the components of the Shettle & Fenn (S&F) model
C =====================================================================
      
      IF(IMOD.EQ.2) THEN
      	 
C* Initialisation
C  --------------
         ! Name of the different particles of the S&F model
	  PARTICULE(1)='SMALL_RURAL'
         PARTICULE(2)='LARGE RURAL'
         PARTICULE(3)='SMALL_URBAN'
         PARTICULE(4)='LARGE_URBAN'
	  PARTICULE(5)='OCEANIC_MODEL'

	  IGRANU =   1	! Log Normal distribution by default
	    
C* Initialisation of the parameters of the S&F model
C  -------------------------------------------
	  CALL SOS_INIT_PARAMSF(DIRFIC,WA,RH,VARGRANU1I,
     &                         VARGRANU2I,MR,MI,IER)
         IF (IER.NE.0) GOTO 997   ! error case

C* Definition of the density of each component for the
C  selected model
C  -----------------------------------------------------------	 
         IF(IMODELE_SF.EQ.1) THEN		!Tropospheric model
	    Ni(1) = 1.0     
	    
         ELSE IF(IMODELE_SF.EQ.2) THEN 	!Urban model
	    Ni(3) = 0.999875
	    Ni(4) = 0.000125
	 
         ELSE IF(IMODELE_SF.EQ.3) THEN	!Maritime model
	    Ni(1) = 0.99
	    Ni(5) = 0.01
	       
         ELSE IF(IMODELE_SF.EQ.4) THEN	!Coastal model
	    Ni(1) = 0.995
	    Ni(5) = 0.005
	       
         ENDIF	!End of the definition of the S&F models


      IF (TRACE) THEN
        WRITE (88,*,ERR=921) "SHETTLE & FENN MODEL"
	 WRITE (88,*,ERR=921) "********************"
	 WRITE (88,*,ERR=921) " "
	 IF(IMODELE_SF.EQ.1) WRITE(88,*,ERR=921) "TROPOSPHERIC Model"
	 IF(IMODELE_SF.EQ.2) WRITE(88,*,ERR=921) "URBAN Model"
	 IF(IMODELE_SF.EQ.3) WRITE(88,*,ERR=921) "MARITIME Model"
	 IF(IMODELE_SF.EQ.4) WRITE(88,*,ERR=921) "COASTAL Model"
	 WRITE (88,*,ERR=921) " "
	 WRITE (88,*,ERR=921) "  Coefficients alphaK of the mixture:"
	 DO I=1,CTE_NB_SF
	    WRITE (88,*,ERR=921) "    ",PARTICULE(I)," : ", Ni(I)
        ENDDO
        WRITE (88,*,ERR=921) " "
	 WRITE (88,*,err=921) "Percentage of relative humidity : ",RH
	 WRITE (88,*,err=921) "Wavelength in microns : ",WA 
	 WRITE (88,*,ERR=921) " "   
	 WRITE (88,*,ERR=921) " "		
      ENDIF !Trace file

      
C* For each aerosol type of the S&F model
C  ----------------------------------------     
	 DO 109 I=1,CTE_NB_SF
	    
	    IF(Ni(I).EQ.0.) GOTO 109
		
C* Building of the name of the MIE file
C  ----------------------------------	   
           IF (I.EQ.1) ALPHAF = CTE_ALPHAMAX_SF_SR
	    IF (I.EQ.3) ALPHAF = CTE_ALPHAMAX_SF_SU
	    IF ((I.NE.1).AND.(I.NE.3)) THEN
	        X = VARGRANU2I(I)  
	        RMAX = VARGRANU1I(I) * DEXP(X*X) *
     &	               DEXP(X*DSQRT(-2.D+00*LOG(CTE_COEF_NRMAX)))
     	
	        ALPHAF = REAL( 100 +
     &		100*DINT(2.*INCTE_PI*RMAX / (100.*WA)) )
	    ENDIF

	    IF ((ALPHAO.GT.ALPHAF).OR.(ALPHAF.GE.1D5)) GOTO 1009 

           CALL SOS_NOM_FICMIE(NBMU_GAUSS,FICANGLES_USER,
     &	                        MR(I),MI(I),ALPHAO,ALPHAF,FICS)       
	    FICMIE = DIRS(1:LDIRS)//'/'//FICS
           LFICS=INDEX(FICMIE,' ')-1
	    
	    IF (TRACE) THEN 
	       WRITE (88,*,ERR=921) "MIE file : ", 
     &	       FICMIE(1:LFICS)	
	    ENDIF	
	 
	    !Screen display
	    WRITE(6,*) "  MIE file : ", FICMIE(1:LFICS)
	    
C* If the MIE file already exists, the program is not run
C  -----------------------------------------------------------           	              
           INQUIRE(FILE=FICMIE,EXIST=EX)
	    
	    IF (EX) THEN  
               WRITE(6,*)'  ==> The file has already been calculated'
	  
	    ELSE   
      
C* Computation of the MIE file
C  ------------------------ 
              CALL SOS_MIE(MIE_NBMU,XMU,XHR,
     &	                    MR(I),MI(I),ALPHAO,ALPHAF,
     &                     FICMIE_TMP,FICLOGMIE,IER)
     	       IF (IER.NE.0) GOTO 993   ! error case
       
              COMMAND=
     &        'mv '//FICMIE_TMP(1:LFICTMP)//' '//FICMIE(1:LFICS)
              LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   !The search starts after the last blank at the beginning of the string
              LCOM = LCOM + (4+LFICTMP)
              IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020

              INQUIRE(FILE=FICMIE,EXIST=EX)
              IF (EX) THEN !If the MIE file exists, the temporary file (generated by another process) is deleted
                 OPEN(UNIT=3,FILE=FICMIE_TMP,FORM='UNFORMATTED')
                 CLOSE(3,STATUS='DELETE') 
              ELSE !Otherwise it is moved to the storage space of the MIE files
                 CALL SYSTEM(COMMAND)  
              ENDIF

            ENDIF
	 
C* Computation of the extinction and scattering cross sections
C  and of the phase functions I, Q, U of the component I
C  ---------------------------------------------------------- 
	    IF (TRACE) THEN
	       WRITE(88,*,err=921) PARTICULE(I)
	       WRITE(88,*,err=921)'============='
	       WRITE(88,*,err=921) ' '
	       IF ((I.NE.1).AND.(I.NE.3)) THEN
	         WRITE(88,*,err=921) 'RMAX (microns) :',RMAX
		  WRITE(88,*,err=921) 'Wavelength (microns) : ',WA
		  WRITE(88,*,err=921) 'Max size parameter : ',ALPHAF
	       ENDIF
	    ENDIF
	       
           CALL SOS_GRANU(FICMIE,IGRANU,VARGRANU1I(I),
     &                    VARGRANU2I(I),CTE_NOT_DEFINED_VALUE_DBLE,
     &                    WA,MIE_NBMU,XMU,TRACE,
     &			     KMAT1I,KMAT2I,SOMME_NR,P11I,P12I,P33I,IER)
           IF (IER.NE.0) GOTO 995   ! error case
	     
C* Computation of the extinction and scattering cross sections of the mixture
C  ---------------------------------------------------------------------
	    KMAT1 = KMAT1 + Ni(I)*KMAT1I
	    KMAT2 = KMAT2 + Ni(I)*KMAT2I
	    
	    DO 108 J=-MIE_NBMU,MIE_NBMU	     
	       P11(J) = P11(J) + Ni(I)*P11I(J)*KMAT2I
	       P12(J) = P12(J) + Ni(I)*P12I(J)*KMAT2I
	       P33(J) = P33(J) + Ni(I)*P33I(J)*KMAT2I
  108	    CONTINUE  !End of loop on the Gauss angles
    
  109	 CONTINUE	 !End of loop on the components of the Shettle and Fenn model
	    

C* Normalisation of the phase functions of the mixture
C  -----------------------------------------------
	 DO 110 J=-MIE_NBMU,MIE_NBMU   
	    P11(J) = P11(J) / KMAT2
	    P12(J) = P12(J) / KMAT2
	    P33(J) = P33(J) / KMAT2
  110	 CONTINUE

C* Case of spherical particles: P22(J) = P11(J)
C  -----------------------------------------------------
        DO 1023 J=-MIE_NBMU,MIE_NBMU
           P22(J)=P11(J)
 1023   CONTINUE	 
 
C* Computation of the parameters alpha, beta, gamma, delta, zeta
C  with or without truncation
C  -----------------------------------------------------  
	 CALL SOS_DECOMPO_LEGENDRE(ITRONC,TRACE,MIE_NBMU,XMU,XHR,OS_NB,
     &				      P11,TTT,P12,P22,P33,COEF_TRONCA,Z1,
     &				      ALP,BETA11,BETA22,
     &                            GAMMA12,DELTA33,ZETA,IER)  
        IF (IER.NE.0) GOTO 996   ! error case

      	 
      ENDIF	!End of test on modality IMOD=2 (Shettle&Fenn)



C =====================================================================
C* Computation of the radiative properties of a bimodal LND size distribution
C =====================================================================
      
      IF(IMOD.EQ.3) THEN
	
      
C*       Retrieval of the information on the description of the mixture
C        ----------------------------------------------------------      
C        Use of the user volume concentrations
C         --> Check that the user parameters are given for the option
         IF (MODE_PARAM_BILND.EQ.1) THEN 
	         
	     CVI(1) = USER_CV_COARSE !"coarse" mode
	     CVI(2) = USER_CV_FINE   !"fine" mode
	  
         ENDIF
       	    
C*       Retrieval of the refractive index and size distribution information of the 2 LND
C        -----------------------------------------------------------------         	 
      
         RN_LND(1) = BMD_CM_MRWA  !"coarse" mode
         RN_LND(2) = BMD_FM_MRWA  !"fine" mode
      
         IN_LND(1) = BMD_CM_MIWA  !"coarse" mode
         IN_LND(2) = BMD_FM_MIWA  !"fine" mode
      
         RMODAL_LND(1) = BMD_CM_RMODAL  !"coarse" mode
         RMODAL_LND(2) = BMD_FM_RMODAL  !"fine" mode
      
         VAR_LND(1) = BMD_CM_VAR  !"coarse" mode
         VAR_LND(2) = BMD_FM_VAR  !"fine" mode
      
      
         IF (MODE_PARAM_BILND.EQ.2) THEN 
	  
	    RN_LND_WAREF(1) = BMD_CM_MRWAREF !"coarse" mode
	    RN_LND_WAREF(2) = BMD_FM_MRWAREF !"fine" mode
	  
	    IN_LND_WAREF(1) = BMD_CM_MIWAREF !"coarse" mode
	    IN_LND_WAREF(2) = BMD_FM_MIWAREF !"fine" mode
	    
         ENDIF
 
         

         DO 111 I=1,2
	 
C* Check of the sign of the argument IN
C  -----------------------------------------
	    IF (IN_LND(I).GT.0.D+00) GOTO 10100	 

C* Values of RN and IN forced to the formats F5.3 and F8.5
C  -----------------------------------------
	    RN_LND(I) = DNINT(RN_LND(I)*1000.D+00) / 1000.D+00
	    IN_LND(I) = -DNINT(-IN_LND(I)*100000.D+00) / 100000.D+00 
	    
111      CONTINUE	 
	 

         IF (TRACE) THEN
           WRITE (88,*,ERR=921) "LOG-NORMAL BIMODAL MODEL"
	    WRITE (88,*,ERR=921) "************************"
	    WRITE (88,*,ERR=921) " "
	    WRITE (88,*,ERR=921) "Wavelength in microns : ", WA
	    WRITE (88,*,ERR=921) " "
	    WRITE(88,*,ERR=921) "Size distribution : log-normal coarse"
	    WRITE(88,*,ERR=921) "                    -----------------"
	    WRITE(88,*,ERR=921) 
     &     "   Refractive index (real part) : ",RN_LND(1)
	    WRITE(88,*,ERR=921) 
     &     "   Refractive index (imaginary part) : ",IN_LND(1)
	    WRITE(88,*,ERR=921) 
	    WRITE(88,*,ERR=921) "   Modal radius (microns) : " , 
     &	    RMODAL_LND(1) 
	    WRITE(88,*,ERR=921) "   Standard deviation : " , VAR_LND(1)
	    WRITE(88,*,ERR=921) " "
	    WRITE(88,*,ERR=921) "Size distribution : log-normal fine"
	    WRITE(88,*,ERR=921) "                    ---------------"
	    WRITE(88,*,ERR=921) 
     &     "   Refractive index (real part) : ",RN_LND(2)
	    WRITE(88,*,ERR=921) 
     &     "   Refractive index (imaginary part) : ",IN_LND(2)
	    WRITE(88,*,ERR=921) 
	    WRITE(88,*,ERR=921) "   Modal radius (microns) : " , 
     &	    RMODAL_LND(2) 
	    WRITE(88,*,ERR=921) "   Standard deviation : " , VAR_LND(2) 
	    WRITE(88,*,ERR=921) " "
	 
	    IF (MODE_PARAM_BILND.EQ.1) THEN 
	       WRITE(88,*,ERR=921)  
     &         "Proportion of components defined by", 
     &	       " volume concentrations :"		    
	       WRITE(88,*,ERR=921) 
     & 	"   Volume concentration for LND coarse mode : ",CVI(1)
	       WRITE(88,*,ERR=921) 
     & 	"   Volume concentration for LND fine mode   : ",CVI(2)
	       WRITE(88,*,ERR=921) " "
	    ENDIF 
	  
	    IF (MODE_PARAM_BILND.EQ.2) THEN 
	       WRITE(88,*,ERR=921) 
     &         "Proportion of components defined by the", 
     &	       "by the ratio AOT_coarse / AOT_totale :"
     	       WRITE(88,*,ERR=921) 
     &	       "   Reference wavelength (microns) : ",WAREF
     	       WRITE(88,*,ERR=921) 
     &	       "   --> Total aerosols extinction optical thickness : ",
     &         TA_WAREF	    
	       WRITE(88,*,ERR=921) 
     &         "   --> ratio AOT coarse / totale : ",RTAUct_WAREF
              WRITE(88,*,ERR=921) " "
	       WRITE(88,*,ERR=921) 
     &	       "   --> Refractive index of the coarse mode "
	       WRITE(88,*,ERR=921) 
     & 	       "          - Real part : ",RN_LND_WAREF(1)
	       WRITE(88,*,ERR=921) 
     & 	       "          - Imaginary part : ",IN_LND_WAREF(1)
	       WRITE(88,*,ERR=921) 
     &	       "   --> Refractive index of the fine mode "
	       WRITE(88,*,ERR=921) 
     & 	       "          - Real part : ",RN_LND_WAREF(2)
	       WRITE(88,*,ERR=921) 
     & 	       "          - Imaginary part : ",IN_LND_WAREF(2)
	       WRITE(88,*,ERR=921) " "
	    ENDIF  	        	    
	 
         ENDIF !Trace file
      	    
	 
 
C* Initialisation
C  --------------
 	  IGRANU =   1	! Log Normal distribution by default


C* Computation of the coefficients CVI for an expression of the
C  mixture as AOT_coarse / AOT_total
C  --------------------------------------------------
         IF (MODE_PARAM_BILND.EQ.2) THEN
	    !Use of the ratio AOT_coarse / AOT_total
 
	    IF (TRACE) THEN
	       WRITE (88,*,ERR=921) " "
	       WRITE (88,*,ERR=921) 
     &	       "Calculation of the mixture coefficients"
	       WRITE (88,*,ERR=921) 
     &	       "======================================="
               WRITE (88,*,ERR=921) "  "
	    ENDIF
	    
C* For each of the two modes (1 : coarse, 2 : fine)
C  -------------------------------------------------  
           DO 114 I=1,2  
		
C* Building of the name of the MIE file (for WAREF)
C  -----------------------------------------------	               
	       X = VAR_LND(I)  
	       RMAX = RMODAL_LND(I) * DEXP(X*X) *
     &	              DEXP(X*DSQRT(-2.D+00*LOG(CTE_COEF_NRMAX)))

              ALPHAF = REAL( 100 +
     &	       	 100*DINT(2.*INCTE_PI*RMAX / (100.*WAREF)) )  

	       IF ((ALPHAO.GT.ALPHAF).OR.(ALPHAF.GE.1D5)) GOTO 1009 	      
        
              CALL SOS_NOM_FICMIE(NBMU_GAUSS,FICANGLES_USER,
     &	                           RN_LND_WAREF(I),IN_LND_WAREF(I),
     &	                           ALPHAO,ALPHAF,FICS)       
	       FICMIE = DIRS(1:LDIRS)//'/'//FICS
              LFICS=INDEX(FICMIE,' ')-1
	       
	       IF (TRACE) THEN 
	          WRITE (88,*,ERR=921) "MIE file: ", 
     &	          FICMIE(1:LFICS)	
	       ENDIF	
	 
	       !Screen display
	       WRITE(6,*) "  MIE file : ", FICMIE(1:LFICS)
	    
C* If the MIE file already exists, the program is not run
C  -----------------------------------------------------------           	              
              INQUIRE(FILE=FICMIE,EXIST=EX)
	    
	       IF (EX) THEN  
                  WRITE(6,*)'  ==> The file has already been calculated'
	  
	       ELSE   
      
C* Computation of the MIE file (for WAREF)
C  ------------------------------------ 
                 CALL SOS_MIE(MIE_NBMU,XMU,XHR,
     &		                RN_LND_WAREF(I),IN_LND_WAREF(I),
     &                        ALPHAO,ALPHAF,
     &                        FICMIE_TMP,FICLOGMIE,IER)
     	          IF (IER.NE.0) GOTO 993   ! error case

                 COMMAND=
     &           'mv '//FICMIE_TMP(1:LFICTMP)//' '//FICMIE(1:LFICS)
                 LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   !The search starts after the last blank at the beginning of the string
                 LCOM = LCOM + (4+LFICTMP)
                 IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020

                 INQUIRE(FILE=FICMIE,EXIST=EX)
                 IF (EX) THEN !If the MIE file exists, the temporary file (generated by another process) is deleted
                    OPEN(UNIT=3,FILE=FICMIE_TMP,FORM='UNFORMATTED')
                    CLOSE(3,STATUS='DELETE') 
                 ELSE !Otherwise it is moved to the storage space of the MIE files
                    CALL SYSTEM(COMMAND)  
                 ENDIF
         
              ENDIF
	 
C* Computation of the extinction and scattering cross sections
C  and of the phase functions I, Q, U of the component I for WAREF
C  ------------------------------------------------------------- 
            
	       IF (TRACE) THEN
	          IF (I.EQ.1) THEN
	             WRITE(88,*,err=921) "   LND coarse mode"
		      WRITE(88,*,err=921) '   ==============='
	          ELSE !I=2 : fine mode
	             WRITE(88,*,err=921) "   LND fine mode" 
		      WRITE(88,*,err=921) '   =============' 
	          ENDIF	  
	          WRITE(88,*,err=921) ' '
	          WRITE(88,*,err=921) '      RMAX (microns) :',RMAX
	          WRITE(88,*,err=921) 
     &	          '      Wavelength (microns) : ',WAREF
	          WRITE(88,*,err=921) 
     &	          '      Max size parameter : ',ALPHAF	  
	          WRITE(88,*,err=921) ' '     
	       ENDIF
	       
              CALL SOS_GRANU(FICMIE,IGRANU,RMODAL_LND(I),
     &                    VAR_LND(I),CTE_NOT_DEFINED_VALUE_DBLE,
     &                    WAREF,MIE_NBMU,XMU,TRACE,
     &			     KMAT1I,KMAT2I,SOMME_NR,P11I,P12I,P33I,IER)
              IF (IER.NE.0) GOTO 995   ! error case

    	     
C* Storage of the extinction coefficients of each mode
C  ---------------------------------------------------------------
	       IF (I.EQ.1)  KMAT1c = KMAT1I	!coarse mode
	       IF (I.EQ.2)  KMAT1f = KMAT1I	!fine mode
      	    	        
  114	    CONTINUE	 !End of loop on the coarse and fine modes of the bimodal LND


C* Computation of the coefficients CVI of each mode
C  ---------------------------------------------------------------
           CVI(1)=   (RTAUct_WAREF * TA_WAREF) / KMAT1c
	    CVI(2)=   ( (1.-RTAUct_WAREF)*TA_WAREF) / KMAT1f

	    IF (TRACE) THEN
	       WRITE(88,*,err=921) 
     &	       '      Estimate : CV coarse mode (non-normalized):',
     &        CVI(1)
	       WRITE(88,*,err=921) 
     &	       '      Estimate : CV fine mode (non-normalized):',CVI(2)
	       WRITE(88,*,err=921) ' '  
	       WRITE(88,*,err=921) ' '   
	    ENDIF
			    
	  ENDIF !End of test on the use of the ratio AOT_coarse / AOT_tot

		
C* Normalisation to 1 of the sum of the coefficients CVI of each mode
C  -----------------------------------------------------------------
         NTOT = CVI(1) + CVI(2)
	  CVI(1) = CVI(1) / NTOT
         CVI(2) = CVI(2) / NTOT
	
         IF (TRACE) THEN
	      WRITE(88,*,err=921)
     &	      '      CV coarse mode (normalized):', CVI(1)
	      WRITE(88,*,err=921)
     &	      '      CV fine mode (normalized):', CVI(2)	  
	      WRITE(88,*,err=921) ' '  
	      WRITE(88,*,err=921) ' '   
	     WRITE (88,*,ERR=921) 
     &	     "Calculation of the radiative properties of the mixture"
	     WRITE (88,*,ERR=921) 
     &	     "======================================================"
            WRITE (88,*,ERR=921) "  "
	  ENDIF
	    
C* For each of the two modes (1 : coarse, 2 : fine)
C  ------------------------------------------------- 
	        
	  DO 112 I=1,2
	    	    
	    IF(CVI(I).EQ.0.) GOTO 112
		
C* Building of the name of the MIE file (for the simulation wavelength)
C  ----------------------------------	               
	    X = VAR_LND(I)   
	    RMAX = RMODAL_LND(I) * DEXP(X*X) *
     &	           DEXP(X*DSQRT(-2.D+00*LOG(CTE_COEF_NRMAX)))
     	
	    ALPHAF = REAL( 100 +
     &	      100*DINT(2.*INCTE_PI*RMAX / (100.*WA)) )

	    IF ((ALPHAO.GT.ALPHAF).OR.(ALPHAF.GE.1D5)) GOTO 1009 

           CALL SOS_NOM_FICMIE(NBMU_GAUSS,FICANGLES_USER,
     & 	                 RN_LND(I),IN_LND(I),ALPHAO,ALPHAF,FICS)     
	    FICMIE = DIRS(1:LDIRS)//'/'//FICS
           LFICS=INDEX(FICMIE,' ')-1
	    
	    IF (TRACE) THEN 
	        WRITE (88,*,ERR=921) "MIE file: ", 
     &	        FICMIE(1:LFICS)	
	    ENDIF	
	 
	    !Screen display
	    WRITE(6,*) "  MIE file : ", FICMIE(1:LFICS)
	    
C* If the MIE file already exists, the program is not run
C  -----------------------------------------------------------           	              
           INQUIRE(FILE=FICMIE,EXIST=EX)
	    
	    IF (EX) THEN  
               WRITE(6,*)'  ==> The file has already been calculated'
	  
	    ELSE   
      
C* Computation of the MIE file (for the simulation wavelength)
C  ------------------------ 
              CALL SOS_MIE(MIE_NBMU,XMU,XHR,
     &                     RN_LND(I),IN_LND(I),ALPHAO,ALPHAF,
     &                     FICMIE_TMP,FICLOGMIE,IER)
     	       IF (IER.NE.0) GOTO 993   ! error case
 
              COMMAND=
     &        'mv '//FICMIE_TMP(1:LFICTMP)//' '//FICMIE(1:LFICS)
              LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   !The search starts after the last blank at the beginning of the string
              LCOM = LCOM + (4+LFICTMP)
              IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020

              INQUIRE(FILE=FICMIE,EXIST=EX)
              IF (EX) THEN !If the MIE file exists, the temporary file (generated by another process) is deleted
                 OPEN(UNIT=3,FILE=FICMIE_TMP,FORM='UNFORMATTED')
                 CLOSE(3,STATUS='DELETE') 
              ELSE !Otherwise it is moved to the storage space of the MIE files
                 CALL SYSTEM(COMMAND)  
              ENDIF
        
           ENDIF
	 
C* Computation of the extinction and scattering cross sections
C  and of the phase functions I, Q, U of the component I
C   (for the simulation wavelength)
C  ---------------------------------------------------------- 
	    IF (TRACE) THEN
	       IF (I.EQ.1) THEN
	          WRITE(88,*,err=921) "   1) LND coarse mode"
		   WRITE(88,*,err=921) '   =================='
	       ELSE !I=2 : fine mode
	          WRITE(88,*,err=921) "   2) LND fine mode" 
		   WRITE(88,*,err=921) '   ================' 
	       ENDIF	  
	       WRITE(88,*,err=921) ' '
	       WRITE(88,*,err=921) '      RMAX (microns) :',RMAX
	       WRITE(88,*,err=921) 
     &	       '      Wavelength (microns) : ',WA
	       WRITE(88,*,err=921) 
     &	       '      Max size parameter : ',ALPHAF	  
	       WRITE(88,*,err=921) ' '     
	    ENDIF
	       
           CALL SOS_GRANU(FICMIE,IGRANU,RMODAL_LND(I),
     &                    VAR_LND(I),CTE_NOT_DEFINED_VALUE_DBLE,
     &                    WA,MIE_NBMU,XMU,TRACE,
     &			     KMAT1I,KMAT2I,SOMME_NR,P11I,P12I,P33I,IER)
           IF (IER.NE.0) GOTO 995   ! error case

		 
C* Computation of the extinction and scattering coefficients of the mixture
C  ---------------------------------------------------------------
	    KMAT1 = KMAT1 + CVI(I)*KMAT1I
	    KMAT2 = KMAT2 + CVI(I)*KMAT2I
	    
	    DO 113 J=-MIE_NBMU,MIE_NBMU	     
	       P11(J) = P11(J) + CVI(I)*P11I(J)*KMAT2I
	       P12(J) = P12(J) + CVI(I)*P12I(J)*KMAT2I
	       P33(J) = P33(J) + CVI(I)*P33I(J)*KMAT2I	       
  113	    CONTINUE  !End of loop on the Gauss angles
    
  112	  CONTINUE	 !End of loop on the coarse and fine modes of the bimodal LND
	    
	 		
		
C* Normalisation of the phase functions of the mixture
C  -----------------------------------------------
	DO 120 J=-MIE_NBMU,MIE_NBMU   
	    P11(J) = P11(J) / KMAT2
	    P12(J) = P12(J) / KMAT2
	    P33(J) = P33(J) / KMAT2
  120	CONTINUE

C* Case of spherical particles: P22(J) = P11(J)
C  -----------------------------------------------------
       DO 1024 J=-MIE_NBMU,MIE_NBMU
          P22(J)=P11(J)
 1024  CONTINUE	 
 
C* Computation of the parameters alpha, beta, gamma, delta, zeta
C  with or without truncation
C  -----------------------------------------------------  
	CALL SOS_DECOMPO_LEGENDRE(ITRONC,TRACE,MIE_NBMU,XMU,XHR,OS_NB,
     &				     P11,TTT,P12,P22,P33,COEF_TRONCA,Z1,
     &				     ALP,BETA11,BETA22,
     &                           GAMMA12,DELTA33,ZETA,IER)     
       IF (IER.NE.0) GOTO 996   ! error case

          	 
      ENDIF	!End of test on modality IMOD=3 (bimodal LND)




C ==============================================================================
C* External data:
C* --> External phase functions and scattering and extinction coefficients
C*
C*     NB: Developed on the basis of the outputs of the Oleg Dubovik tool (LOA)
C*          which gives the primary scattering properties of homogeneous
C*          randomly oriented spheroid particles.
C*          A parameter giving the number of angles describing the functions
C*          had to be added, though.
C*
C*          Be careful that a ratio -F12 / F11 is read.
C*          Importance of the sign: F12 < 0 for Rayleigh scattering.

C ==============================================================================

      IF(IMOD.EQ.4) THEN
	 
    	
C* Reading of the external data file
C  -----------------------------------------------------------------
         OPEN(UNIT=20,FILE=FICEXTDATA,STATUS='OLD',ERR=940)

	  READ(20,'(a)',ERR=941) LIGNE_TEXTE
         IFIN=INDEX(LIGNE_TEXTE,':')
         READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) KMAT1
	 
         READ(20,'(a)',ERR=941) LIGNE_TEXTE
         IFIN=INDEX(LIGNE_TEXTE,':')
         READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) KMAT2
	 
	  READ(20,'(a)',ERR=941) LIGNE_TEXTE
         IFIN=INDEX(LIGNE_TEXTE,':')
         READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) NB_FICEXT_ANGLES



	  !Check of the number of angles with respect to the size of the arrays
	  IF (NB_FICEXT_ANGLES.GT.CTE_MAXNB_ANG_EXT) GOTO 950
	
	  READ(20,*,ERR=941) 

	  DO J=1,NB_FICEXT_ANGLES
            READ(20,*,ERR=941,END=942) 
     &	                 ANG_EXT(J), F11_EXT(J), F12sF11_EXT(J),
     &	                 F22sF11_EXT(J), F33sF11_EXT(J)
         ENDDO 
	
	  CLOSE(20)
	 
	  IF (TRACE) THEN
            WRITE (88,*,ERR=921) "EXTERNAL DATA (USER'S DATA)"
	     WRITE (88,*,ERR=921) "***************************"
	     WRITE (88,*,ERR=921) " "
	     WRITE (88,*,ERR=921)  "Filename : ",FICEXTDATA
	     WRITE (88,*,ERR=921) " "  		 
            WRITE (88,*,ERR=921) "Extinction coef : ", KMAT1
	     WRITE (88,*,ERR=921) "Scattering coef  : ", KMAT2
	     WRITE (88,*,ERR=921) "Nb of angles : ", NB_FICEXT_ANGLES
	     WRITE (88,*,ERR=921) 
     &	     "  ANGLE     F11         -F12/F11      F22/F11      F33/F11"
	     DO J=1,NB_FICEXT_ANGLES
               WRITE (88,402,ERR=921) ANG_EXT(J), F11_EXT(J), 
     &	       			   F12sF11_EXT(J), F22sF11_EXT(J), 
     &	                               F33sF11_EXT(J)
            ENDDO 
	     WRITE (88,*,ERR=921) " "
         ENDIF !Trace file
	 

C* Conversion to the functions F12, F22 and F33 for the mu of the angles
C* describing the external functions
C  -----------------------------------------------------------------	 
	  DO J=1,NB_FICEXT_ANGLES
	     MU_EXT(J)  = DCOS(ANG_EXT(J)*INCTE_PI/180.D+00)
            F12_EXT(J) = -F12sF11_EXT(J)*F11_EXT(J)
	     F22_EXT(J) = F22sF11_EXT(J)*F11_EXT(J)
	     F33_EXT(J) = F33sF11_EXT(J)*F11_EXT(J)
         ENDDO 
	 
C* Conversion to the functions F11, F12, F22 and F33 for the Gauss mu
C  -----------------------------------------------------------------	

C*      P11 : The routine is called by passing XMU and P11
C*            which point to the first useful element
C*            (since the covered domain is larger than the useful one
C*             and the array starts at a negative index)
C*      ---	
	  CALL SOS_INTERPO_SPLINT(NB_FICEXT_ANGLES, MU_EXT, F11_EXT,
     &	                          2*MIE_NBMU+1, XMU(-MIE_NBMU), 
     &                           P11(-MIE_NBMU), IER)	
	  IF (IER.NE.0) GOTO 945   ! error case
	
C*      P12 :
C*      ---	
	  CALL SOS_INTERPO_SPLINT(NB_FICEXT_ANGLES, MU_EXT, F12_EXT,
     &	                          2*MIE_NBMU+1, XMU(-MIE_NBMU), 
     &                           P12(-MIE_NBMU), IER)	
	  IF (IER.NE.0) GOTO 945   ! error case
	
C*      P22 :
C*      ---	
	  CALL SOS_INTERPO_SPLINT(NB_FICEXT_ANGLES, MU_EXT, F22_EXT,
     &	                          2*MIE_NBMU+1, XMU(-MIE_NBMU), 
     &                           P22(-MIE_NBMU), IER)	
	  IF (IER.NE.0) GOTO 945   ! error case
		
C*      P33 :
C*      ---	
	  CALL SOS_INTERPO_SPLINT(NB_FICEXT_ANGLES, MU_EXT, F33_EXT,
     &	                          2*MIE_NBMU+1, XMU(-MIE_NBMU), 
     &                           P33(-MIE_NBMU), IER)		
	  IF (IER.NE.0) GOTO 945   ! error case
	
	!C* Add a recomposition module for the TRACE file
	!C* File headers
	!ADJUST DEFINITION OF THE TABLE LIMIT F11, etc. hard-coded ==> in SOS.h
	!Careful with the sign: -F12sF11
	 
	  IF (TRACE) THEN
	    WRITE (88,*,ERR=921) " Functions F11, F12, F22 and F33 for"
	    WRITE (88,*,ERR=921) " the original angles of the datafile"
	    WRITE (88,*,ERR=921) " ==================================="
	    WRITE (88,*,ERR=921) 
     &	    "  ANGLE     F11          F12          F22          F33"
	    DO J=1,NB_FICEXT_ANGLES
              WRITE (88,402,ERR=921) ANG_EXT(J), F11_EXT(J), 
     &	       			  F12_EXT(J), F22_EXT(J), F33_EXT(J)
           ENDDO 
	    WRITE (88,*,ERR=921) " " 
	    WRITE (88,*,ERR=921) " Functions  F11, F12, F22 and F33"
	    WRITE (88,*,ERR=921) " for the simulation angles "
	    WRITE (88,*,ERR=921) " ================================" 
	    DO J=-MIE_NBMU, MIE_NBMU
	       TETA=DACOS(XMU(J))*180./INCTE_PI
              WRITE (88,402,ERR=921) TETA,P11(J),P12(J),P22(J),P33(J)
           ENDDO 	  		
         ENDIF !Trace file
	 
	 	
	
	
	 
C* Computation of the parameters alpha, beta, gamma, delta, zeta
C  with or without truncation
C  -----------------------------------------------------
	  CALL SOS_DECOMPO_LEGENDRE(ITRONC,TRACE,MIE_NBMU,XMU,XHR,OS_NB,
     &				       P11,TTT,P12,P22,P33,COEF_TRONCA,Z1,
     &				       ALP,BETA11,BETA22,
     &                             GAMMA12,DELTA33,ZETA,IER)  
         IF (IER.NE.0) GOTO 996   ! error case

      ENDIF	!End of test on modality IMOD=4 (external data)
      
      
      


C ==============================================================================
C* Definition of a mixture of models by the user
C ==============================================================================

      IF(IMOD.EQ.5) THEN

         SOM_TAUX_AOT = 0.	 

    	  LFIC=INDEX(FICMIXTURE_AER,' ')-1

C* Reading of the file defining the mixture of aerosol modes
C  -----------------------------------------------------------------
         OPEN(UNIT=20,FILE=FICMIXTURE_AER,STATUS='OLD',ERR=960)
	 
	  READ(20,'(a)',ERR=961) LIGNE_TEXTE
         IFIN=INDEX(LIGNE_TEXTE,':')
         READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) NB_MODES_MELANGE
	
         DO 122 NMOD=1,NB_MODES_MELANGE 

            !Reading of the model type and check against the available models
	     READ(20,'(a)',ERR=961) LIGNE_TEXTE
            IFIN=INDEX(LIGNE_TEXTE,':')
            READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) CMODEL
            IFIN=INDEX(CMODEL,' ')-1
            CMODEL=CMODEL(1:IFIN)
 
            IF (CMODEL.EQ."LND") THEN
               MEL_MODEL_TYPE(NMOD)=1
            ELSE IF (CMODEL.EQ."JUNGE") THEN 
               MEL_MODEL_TYPE(NMOD)=2
            ELSE  !Neither LND nor Junge
               GOTO 962
            ENDIF

            !Reading of the model parameters
            IF (CMODEL.EQ."LND") THEN

               READ(20,'(a)',ERR=961) LIGNE_TEXTE
               IFIN=INDEX(LIGNE_TEXTE,':')
               READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) 
     &         MEL_MODEL_PARAM(NMOD,1)			!LND: MODAL RADIUS

               READ(20,'(a)',ERR=961) LIGNE_TEXTE
               IFIN=INDEX(LIGNE_TEXTE,':')
               READ(LIGNE_TEXTE(IFIN+1:LENLINE),*)
     &         MEL_MODEL_PARAM(NMOD,2)			!LND: STANDARD DEVIATION

            ELSE 
               IF (CMODEL.EQ."JUNGE") THEN

                  READ(20,'(a)',ERR=961) LIGNE_TEXTE
                  IFIN=INDEX(LIGNE_TEXTE,':')
                  READ(LIGNE_TEXTE(IFIN+1:LENLINE),*)
     &            MEL_MODEL_PARAM(NMOD,2)			!JUNGE: SLOPE

                  READ(20,'(a)',ERR=961) LIGNE_TEXTE
                  IFIN=INDEX(LIGNE_TEXTE,':')
                  READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) 
     &            MEL_MODEL_PARAM(NMOD,1)		       !JUNGE: MIN RADIUS

                  READ(20,'(a)',ERR=961) LIGNE_TEXTE
                  IFIN=INDEX(LIGNE_TEXTE,':')
                  READ(LIGNE_TEXTE(IFIN+1:LENLINE),*)
     &            MEL_MODEL_PARAM(NMOD,3)			!JUNGE: MAX RADIUS

               ENDIF
            ENDIF

            !Reading of the refractive indices
            !at the wavelength of the radiance simulations
	     READ(20,'(a)',ERR=961) LIGNE_TEXTE
            IFIN=INDEX(LIGNE_TEXTE,':')
	     READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) RN_WASIMU_LU
	     READ(20,'(a)',ERR=961) LIGNE_TEXTE
            IFIN=INDEX(LIGNE_TEXTE,':')
	     READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) IN_WASIMU_LU

            !Reading of the refractive indices
            !at the reference wavelength of the AOT
	     READ(20,'(a)',ERR=961) LIGNE_TEXTE
            IFIN=INDEX(LIGNE_TEXTE,':')
	     READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) RN_WAREF_LU
	     READ(20,'(a)',ERR=961) LIGNE_TEXTE
            IFIN=INDEX(LIGNE_TEXTE,':')
	     READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) IN_WAREF_LU

            MEL_MODEL_INDEX_WAREF(NMOD,1) = RN_WAREF_LU
            MEL_MODEL_INDEX_WAREF(NMOD,2) = IN_WAREF_LU
            MEL_MODEL_INDEX_WASIMU(NMOD,1) = RN_WASIMU_LU
            MEL_MODEL_INDEX_WASIMU(NMOD,2) = IN_WASIMU_LU
 
            !Reading of the AOT at WAREF
	     READ(20,'(a)',ERR=961) LIGNE_TEXTE
            IFIN=INDEX(LIGNE_TEXTE,':')
	     READ(LIGNE_TEXTE(IFIN+1:LENLINE),*) AOT_TAUX

            MEL_MODEL_AOT(NMOD) = TA_WAREF*AOT_TAUX

            SOM_TAUX_AOT = SOM_TAUX_AOT + AOT_TAUX
 122     CONTINUE !End of loop on the models of the user file
	
	  CLOSE(20)

	 
	  IF (TRACE) THEN
            WRITE (88,*,ERR=921) "AEROSOLS MIXTURE"
	     WRITE (88,*,ERR=921) "****************"
	     WRITE (88,*,ERR=921) " "
	     WRITE (88,*,ERR=921)  "Filename : ",FICMIXTURE_AER(1:LFIC)
	     WRITE (88,*,ERR=921) " "  		 
            WRITE (88,*,ERR=921) "Number of modes :", NB_MODES_MELANGE
            DO NMOD=1,NB_MODES_MELANGE
	        WRITE (88,*,ERR=921) ""
	        WRITE (88,*,ERR=921) " - Mode: ",NMOD
               IF (MEL_MODEL_TYPE(NMOD).EQ.1) THEN
	           WRITE (88,*,ERR=921) 
     &            "   - Type of size distribution : LND"
	           WRITE(88,*,ERR=921) 
     &            "        Modal radius (microns) : " , 
     &	           MEL_MODEL_PARAM(NMOD,1)
	           WRITE(88,*,ERR=921)  
     &            "        Variance (sig) : ", 
     &	           MEL_MODEL_PARAM(NMOD,2)
               ELSE
	           WRITE (88,*,ERR=921) 
     &            "   - Type of size distribution : JUNGE"
	           WRITE(88,*,ERR=921)  
     &            "        Slope : " , MEL_MODEL_PARAM(NMOD,2)
	           WRITE(88,*,ERR=921)  
     &            "        Minimal radius (microns) : ",
     &	           MEL_MODEL_PARAM(NMOD,1)
	           WRITE(88,*,ERR=921)  
     &            "        Maximal radius (microns) : ", 
     &	           MEL_MODEL_PARAM(NMOD,3)
               ENDIF
	        WRITE (88,*,ERR=921) 
     &         "   - Wavelength WA for aerosols simu (microns) : ", WA
	        WRITE (88,*,ERR=921) 
     &         "   - Refractive index (mr+j*mi) at WA: ",
     &         MEL_MODEL_INDEX_WASIMU(NMOD,1), 
     &         MEL_MODEL_INDEX_WASIMU(NMOD,2)
	        WRITE (88,*,ERR=921) 
     &         "   - Reference wavelength for AOT (mic.) ==> WAREF: ",
     &         WAREF
	        WRITE (88,*,ERR=921) 
     &         "   - Refractive index (mr+j*mi) at WAREF: ",
     &         MEL_MODEL_INDEX_WAREF(NMOD,1), 
     &         MEL_MODEL_INDEX_WAREF(NMOD,2)

	        WRITE (88,*,ERR=921) 
     &         "   - AOT ratio at WAREF: ",MEL_MODEL_AOT(NMOD)/TA_WAREF
	        WRITE (88,*,ERR=921) 
     &         "   - AOT at WAREF: ",MEL_MODEL_AOT(NMOD)
            ENDDO
	     WRITE (88,*,ERR=921) " "
            WRITE (88,*,ERR=921)  
     &      "Sommation of AOT ratios:",SOM_TAUX_AOT
	     WRITE (88,*,ERR=921) " "
         ENDIF !Trace file


         !Test that the sum of the AOT proportion rates over the different modes is 1
         IF (DABS(SOM_TAUX_AOT-1.D+00).GT.CTE_GAP_TOLER_SUM_RATES) THEN
            IF (TRACE) THEN
               WRITE (88,*,ERR=921)   
     &         "The sum of AOT rates is not equal to 1"
               WRITE (88,*,ERR=921)  
     &         "Difference to 1 exceeding the threshold ", 
     &         "CTE_GAP_TOLER_SUM_RATES:", CTE_GAP_TOLER_SUM_RATES
	        WRITE (88,*,ERR=921) "Processing stopped"
            ENDIF 
            GOTO 963	!Difference to 1 too large ==> exit with error
         ELSE
            IF (SOM_TAUX_AOT.NE.1.D+00) THEN
                 DO NMOD=1,NB_MODES_MELANGE 
                    MEL_MODEL_AOT(NMOD) =  
     &              MEL_MODEL_AOT(NMOD) / SOM_TAUX_AOT
                 ENDDO
                 IF (TRACE) THEN
                    WRITE (88,*,ERR=921)  
     &              "==> Normalisation of AOT ratios ",
     &              "so that their sum is strictly equal to 1"
                    DO NMOD=1,NB_MODES_MELANGE	      
	                WRITE (88,*,ERR=921) "    - Mode: ",NMOD
	                WRITE (88,*,ERR=921)                
     &                 "      Normalized AOT ratio at WAREF: ",                
     &                 MEL_MODEL_AOT(NMOD)/TA_WAREF
	                WRITE (88,*,ERR=921) 
     &                 "      ==> AOT at WAREF: ",MEL_MODEL_AOT(NMOD)
                    ENDDO
	             WRITE (88,*,ERR=921) " "
                ENDIF
            ENDIF  
         ENDIF

         


C* Computation loops on the modes:
C  --------------------------------  
           !----------------
           !COMPUTATIONS AT WAREF: Estimation of the weighting coefficients Alpha
           !----------------
           SOM_COEF_ALPHA = 0.
		
           DO 124 NMOD=1,NB_MODES_MELANGE
		
              IF (MEL_MODEL_TYPE(NMOD).EQ.1) THEN   !LND size distribution case
	          X = MEL_MODEL_PARAM(NMOD,2)
	          RMAX = MEL_MODEL_PARAM(NMOD,1) * DEXP(X*X) *
     &	                 DEXP(X*DSQRT(-2.D+00*LOG(CTE_COEF_NRMAX)))

	       ELSE					  !JUNGE size distribution case
	          RMAX = MEL_MODEL_PARAM(NMOD,3)
	       ENDIF

              ALPHAF = REAL( 100 +
     &	       	100*DINT(2.*INCTE_PI*RMAX / (100.*CTE_WAMIN)) )

	       IF ((ALPHAO.GT.ALPHAF).OR.(ALPHAF.GE.1D5)) GOTO 1009   
 
              IF (TRACE) THEN
	          WRITE(88,*,err=921) ' '
	          WRITE(88,*,err=921) '      RMAX (microns) :',RMAX
	          WRITE(88,*,err=921) 
     &	          '      Maximal size parameter : ',ALPHAF	  
	          WRITE(88,*,err=921) ' '     
              ENDIF

              !Building of the name of the MIE file (for refractive indices of the Mode
              !                                    at the reference wavelength
              !                                    and limit ALPHAF estimated at WAMIN)
              CALL SOS_NOM_FICMIE(NBMU_GAUSS,FICANGLES_USER,
     &	                           MEL_MODEL_INDEX_WAREF(NMOD,1),
     &	                           MEL_MODEL_INDEX_WAREF(NMOD,2),
     &	                           ALPHAO,ALPHAF,FICS)       
	       FICMIE = DIRS(1:LDIRS)//'/'//FICS
              LFICS=INDEX(FICMIE,' ')-1
	       
	       IF (TRACE) THEN 
	           WRITE (88,*,ERR=921) "Calculations for WA = ", WAREF
	           WRITE (88,*,ERR=921) "MIE file: ", FICMIE(1:LFICS)
	       ENDIF	
	 
	       !Screen display
	       WRITE(6,*) "  MIE file : ", FICMIE(1:LFICS)
	    
              !If the MIE file already exists, the program is not run
              !-----------------------------------------------------------     
              INQUIRE(FILE=FICMIE,EXIST=EX)
	    
	       IF (EX) THEN  
                  WRITE(6,*)
     &	           '  ==> The file has already been calculated'
	  
	       ELSE   
      
                 !Computation of the MIE file (for limit ALPHAF)
                 !--------------------------------------------- 
                 CALL SOS_MIE(MIE_NBMU,XMU,XHR,
     &	                       MEL_MODEL_INDEX_WAREF(NMOD,1),
     &	                       MEL_MODEL_INDEX_WAREF(NMOD,2),
     &                        ALPHAO,ALPHAF,
     &                        FICMIE_TMP,FICLOGMIE,IER)
     	          IF (IER.NE.0) GOTO 993   ! error case

                 COMMAND=
     &           'mv '//FICMIE_TMP(1:LFICTMP)//' '//FICMIE(1:LFICS)
                 LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   
                 LCOM = LCOM + (4+LFICTMP)
                 IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020

                 INQUIRE(FILE=FICMIE,EXIST=EX)
                 IF (EX) THEN !If the MIE file exists, the temporary file (generated by another process) is deleted
                    OPEN(UNIT=3,FILE=FICMIE_TMP,FORM='UNFORMATTED')
                    CLOSE(3,STATUS='DELETE') 
                 ELSE !Otherwise it is moved to the storage space of the MIE files
                    CALL SYSTEM(COMMAND)  
                 ENDIF
         
              ENDIF
	 
              ! Computation of the extinction and scattering cross sections
              ! and of the phase functions I, Q, U of the component I for WAREF
              ! ------------------------------------------------------------- 
 	       
              !Call of SOS_GRANU for the model IGRANU=MEL_MODEL_TYPE(NMOD)
              !and the parameters VARGRANU1 = MEL_MODEL_PARAM(NMOD,1)
              !and VARGRANU2 = MEL_MODEL_PARAM(NMOD,2)
              CALL SOS_GRANU(FICMIE,MEL_MODEL_TYPE(NMOD),
     &                       MEL_MODEL_PARAM(NMOD,1),
     &                       MEL_MODEL_PARAM(NMOD,2),
     &                       MEL_MODEL_PARAM(NMOD,3),
     &                       WAREF,MIE_NBMU,XMU,TRACE,
     &			        KMAT1I,KMAT2I,SOMME_NR,P11I,P12I,P33I,IER)

              IF (IER.NE.0) GOTO 995   ! error case

    	     
              ! Computation of the coefficients Alpha before normalisation
              ! The computation must be made at WAREF since the reference AOT is at this wavelength
              ! ---------------------------------------------------------------------	
              TAUEXTI = MEL_MODEL_AOT(NMOD)
              COEF_ALPHA(NMOD) = TAUEXTI/KMAT1I
              SOM_COEF_ALPHA = SOM_COEF_ALPHA + COEF_ALPHA(NMOD)
   
  124	    CONTINUE	 !End of loop on the modes

           ! Finalisation of the computation of the weighting coefficients of the modes
           ! ----------------------------------------------------------------
           DO NMOD=1,NB_MODES_MELANGE
              COEF_ALPHA(NMOD) = COEF_ALPHA(NMOD) / SOM_COEF_ALPHA
           ENDDO

            
           IF (TRACE) THEN
	       WRITE(88,*,err=921) ' '
	       WRITE(88,*,err=921) '      Coef alpha :'
              DO NMOD=1,NB_MODES_MELANGE
	          WRITE(88,*,err=921) 
     &	          '      Mode, coef alpha : ',NMOD, COEF_ALPHA(NMOD)
	       ENDDO  
	       WRITE(88,*,err=921) ' '     
           ENDIF



           !----------------
           !COMPUTATIONS AT WASIMU: Application of the mixture of modes
           !----------------
           KMAT1 = 0.
           KMAT2 = 0.

           DO 125 NMOD=1,NB_MODES_MELANGE
		
              IF (MEL_MODEL_TYPE(NMOD).EQ.1) THEN   !LND size distribution case
	          X = MEL_MODEL_PARAM(NMOD,2)
	          RMAX = MEL_MODEL_PARAM(NMOD,1) * DEXP(X*X) *
     &	                 DEXP(X*DSQRT(-2.D+00*LOG(CTE_COEF_NRMAX)))

	       ELSE					  !JUNGE size distribution case
	          RMAX = MEL_MODEL_PARAM(NMOD,3)
	       ENDIF

              ALPHAF = REAL( 100 +
     &	       	100*DINT(2.*INCTE_PI*RMAX / (100.*CTE_WAMIN)) )

	       IF ((ALPHAO.GT.ALPHAF).OR.(ALPHAF.GE.1D5)) GOTO 1009   
 
              IF (TRACE) THEN
	          WRITE(88,*,err=921) ' '
	          WRITE(88,*,err=921) '      RMAX (microns) :',RMAX
	          WRITE(88,*,err=921) 
     &	          '      Maximal size parameter : ',ALPHAF	  
	          WRITE(88,*,err=921) ' '     
              ENDIF

              IF (WA.EQ.WAREF) THEN	!If WA == WAREF, the routine is called for WAREF
					!==> The Mie computations must be made for the indices at WAREF
                 RN = MEL_MODEL_INDEX_WAREF(NMOD,1) 
                 IN = MEL_MODEL_INDEX_WAREF(NMOD,2) 
              ELSE !WASIMU		!If WA != WAREF, the routine is called for WA
					!==> The Mie computations must be made for the indices at WA
                 RN = MEL_MODEL_INDEX_WASIMU(NMOD,1) 
                 IN = MEL_MODEL_INDEX_WASIMU(NMOD,2) 
              ENDIF

              !Building of the name of the MIE file (for refractive indices of the Mode
              !                                    at the simulation wavelength
              !                                    and limit ALPHAF estimated at WAMIN)
              CALL SOS_NOM_FICMIE(NBMU_GAUSS,FICANGLES_USER,
     &	                           RN,IN,ALPHAO,ALPHAF,FICS)       
	       FICMIE = DIRS(1:LDIRS)//'/'//FICS
              LFICS=INDEX(FICMIE,' ')-1
	       
	       IF (TRACE) THEN 
	           WRITE (88,*,ERR=921) "Calculations for WA = ", WA
	           WRITE (88,*,ERR=921) "MIE file: ", FICMIE(1:LFICS)
	       ENDIF	
	 
	       !Screen display
	       WRITE(6,*) "  MIE file : ", FICMIE(1:LFICS)
	    
              !If the MIE file already exists, the program is not run
              !-----------------------------------------------------------     
              INQUIRE(FILE=FICMIE,EXIST=EX)
	    
	       IF (EX) THEN  
                  WRITE(6,*)
     &	           '  ==> The file has already been calculated'
	  
	       ELSE   
      
                 !Computation of the MIE file (for limit ALPHAF)
                 !--------------------------------------------- 
                 CALL SOS_MIE(MIE_NBMU,XMU,XHR,
     &	                       RN,IN,ALPHAO,ALPHAF,
     &                        FICMIE_TMP,FICLOGMIE,IER)
     	          IF (IER.NE.0) GOTO 993   ! error case

                 COMMAND=
     &           'mv '//FICMIE_TMP(1:LFICTMP)//' '//FICMIE(1:LFICS)
                 LCOM = INDEX(COMMAND((5+LFICTMP):CTE_LENCOM),' ')   
                 LCOM = LCOM + (4+LFICTMP)
                 IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1020

                 INQUIRE(FILE=FICMIE,EXIST=EX)
                 IF (EX) THEN !If the MIE file exists, the temporary file (generated by another process) is deleted
                    OPEN(UNIT=3,FILE=FICMIE_TMP,FORM='UNFORMATTED')
                    CLOSE(3,STATUS='DELETE') 
                 ELSE !Otherwise it is moved to the storage space of the MIE files
                    CALL SYSTEM(COMMAND)  
                 ENDIF
         
              ENDIF
	 
              ! Computation of the extinction and scattering cross sections
              ! and of the phase functions I, Q, U of the component I for WAREF
              ! ------------------------------------------------------------- 
            	       
              !Call of SOS_GRANU for the model IGRANU=MEL_MODEL_TYPE(NMOD)
              !and the parameters VARGRANU1 = MEL_MODEL_PARAM(NMOD,1)
              !and VARGRANU2 = MEL_MODEL_PARAM(NMOD,2)
              CALL SOS_GRANU(FICMIE,MEL_MODEL_TYPE(NMOD),
     &                       MEL_MODEL_PARAM(NMOD,1),
     &                       MEL_MODEL_PARAM(NMOD,2),
     &                       MEL_MODEL_PARAM(NMOD,3),
     &                       WA,MIE_NBMU,XMU,TRACE,
     &			        KMAT1I,KMAT2I,SOMME_NR,P11I,P12I,P33I,IER)

              IF (IER.NE.0) GOTO 995   ! error case

    	     
              ! Computation of the terms of the extinction and scattering cross sections of the mixture
              ! Computation of the functions of the phase matrix of the mixture
              ! ---------------------------------------------------------------------	
              KMAT1 = KMAT1 + COEF_ALPHA(NMOD)*KMAT1I
              KMAT2 = KMAT2 + COEF_ALPHA(NMOD)*KMAT2I
	    
	       DO 134 J=-MIE_NBMU,MIE_NBMU	     
	          P11(J) = P11(J) + P11I(J)*COEF_ALPHA(NMOD)*KMAT2I
	          P12(J) = P12(J) + P12I(J)*COEF_ALPHA(NMOD)*KMAT2I
	          P33(J) = P33(J) + P33I(J)*COEF_ALPHA(NMOD)*KMAT2I
  134	       CONTINUE  !End of loop on the Gauss angles
             
  125	    CONTINUE	 !End of loop on the modes

           ! Finalisation of the computation of the functions of the phase matrix of the mixture
	    ! ----------------------------------------------------------------------	     
	    DO 135 J=-MIE_NBMU,MIE_NBMU	     
	       P11(J) = P11(J) / KMAT2
	       P12(J) = P12(J) / KMAT2
	       P33(J) = P33(J) / KMAT2
  135	    CONTINUE  !End of loop on the Gauss angles

C* Case of spherical particles: P22(J) = P11(J)
C  -----------------------------------------------------
       DO 1034 J=-MIE_NBMU,MIE_NBMU
          P22(J)=P11(J)
 1034  CONTINUE
	 
C* Computation of the parameters alpha, beta, gamma, delta, zeta
C  with or without truncation
C  -----------------------------------------------------
	  CALL SOS_DECOMPO_LEGENDRE(ITRONC,TRACE,MIE_NBMU,XMU,XHR,OS_NB,
     &				       P11,TTT,P12,P22,P33,COEF_TRONCA,Z1,
     &				       ALP,BETA11,BETA22,
     &                             GAMMA12,DELTA33,ZETA,IER)  
         IF (IER.NE.0) GOTO 996   ! error case

      ENDIF	!End of test on modality IMOD=5 (user mixture)
      
      
            
      
      
C ======================================
C* Computations common to all the models
C ======================================
	 
C* Computation of the single scattering albedo
C-------------------------------------------
      PIZ=KMAT2/KMAT1

C* Truncation of the aerosol single scattering albedo
C---------------------------------------------------------
      PIZTR = PIZ*(1.-COEF_TRONCA/2.)/(1.-PIZ*COEF_TRONCA/2.)      
    
C* Writing of the trace file
C-------------------------------------------
      IF (TRACE) THEN
        WRITE(88,*,err=921) ' '
	 WRITE(88,*,err=921) ' '
	 WRITE(88,*,err=921) '=============================='
        WRITE(88,*,err=921) '* ROUTINE AEROSOLS : RESULTS *'
	 WRITE(88,*,err=921) '=============================='
        WRITE(88,*,err=921) ' '
	 WRITE(88,*,err=921) ' '
	
	 IF(ITRONC.EQ.1) THEN
	      WRITE(88,*,err=921) '  Calculation with truncation'
	 ELSE
	      WRITE(88,*,err=921) '  Calculation without truncation'	   
	 ENDIF
	 
        WRITE(88,*,err=921) ' '
            	 
	 IF((COEF_TRONCA.NE.0.D+00).OR.(IMOD.NE.0)) THEN
	    WRITE(88,*,err=921) '  Truncation coefficient :',COEF_TRONCA
	    WRITE(88,*,err=921) 
     &	    '  Extinction cross-section (mic^2) - no trunca. applied:'
     &          ,KMAT1
	    WRITE(88,*,err=921) 
     &	    '  Scattering cross-section (mic^2) - no trunca. applied: '
     &          ,KMAT2
	    IF (COEF_TRONCA.NE.0.D+00) THEN
	       WRITE(88,*,err=921) 
     &		 '  Single scattering albedo '//
     &         ' - adjusted to the truncation: ',PIZTR 
	       WRITE(88,*,err=921) 
     &		 '  Single scattering albedo - no trunca. applied:',
     &             PIZTR / (1.-COEF_TRONCA/2.+PIZTR*COEF_TRONCA/2.)      
	       WRITE(88,*,err=921)
     &		 '  Beta11(0) non-normalized '//
     &         ' - adjusted to the truncation: ',Z1  
              WRITE(88,*,err=921) 
     &	        '  Asymmetry factor '//
     &         ' - adjusted to the truncation: ',BETA11(1)/3.
              WRITE(88,*,err=921) 
     &	        '  Asymmetry factor - no trunca. applied: ',
     &          COEF_TRONCA/2.+(1.-COEF_TRONCA/2.)*BETA11(1)/3.     
	    ELSE
	       WRITE(88,*,err=921) 
     & 	 ' Single scattering albedo - no trunca. applied: ',PIZTR
	       WRITE(88,*,err=921) 
     & 	'  Beta11(0) non-normalized - no trunca. applied: ',Z1
              WRITE(88,*,err=921) 
     &	       '  Asymmetry factor - no trunca. applied: ',BETA11(1)/3.
     	    ENDIF !End of test on truncation
           WRITE(88,*,err=921) ' '
	    WRITE(88,*,err=921) '  1st column : Scattering angle '
	    WRITE(88,*,err=921) 
     &	       '  2nd column : Truncated phase function P11' 	 
	    WRITE(88,*,err=921)
     &        '  3rd column : Original phase function P11 '
	    WRITE(88,*,err=921) 
     &        '  4th column : Truncated polarized phase function '
           WRITE(88,*,err=921) '		     P12(mu) truncated'
	    WRITE(88,*,err=921)
     &        '  5th column : Original polarized phase function '
           WRITE(88,*,err=921) '		     P12(mu) '
	    WRITE(88,*,err=921)
     &        '  6th column : Truncated polarized phase function '
           WRITE(88,*,err=921) '                 P33(mu) truncated'
	    WRITE(88,*,err=921) 
     &        '  7th column : Original polarized phase function '
           WRITE(88,*,err=921) '                 P33(mu) '
	    WRITE(88,*,err=921) ' '	    
           DO 44 J=-MIE_NBMU,MIE_NBMU
	       TETA=DACOS(XMU(J))*180./INCTE_PI
	       WRITE(88,10,err=921) TETA,P11(J),TTT(J),
     &                              P12(J)*P11(J)/TTT(J),P12(J),
     &                              P33(J)*P11(J)/TTT(J),P33(J)
   44      CONTINUE
        END IF
    
        WRITE(88,*,err=921) ' '  
	 WRITE(88,35,err=921)'I','ALP(I)','BETA11(I)','GAMMA12(I)',
     &                       'DELTA33(I)','ZETA(I)'  
        DO 45 K=0,OS_NB
	    WRITE(88,25,err=921) K,ALP(K),BETA11(K),GAMMA12(K),
     &	                         DELTA33(K),ZETA(K)
   45   CONTINUE
   	 WRITE(88,*,err=921) ' '
	 	   	   	   	   	   
      END IF

C* Recording of the size distribution parameters in the Result file
C  -----------------------------------------------------------------------
 1111 OPEN(1,FILE=FICGRANU,ERR=913)

      ! Extinction cross section
      WRITE(1,40,ERR=923) KMAT1

      ! Scattering cross section
      WRITE(1,41,ERR=923) KMAT2
     
      ! Asymmetry factor (not truncated)
      WRITE(1,42,ERR=923) 
     &	 ( COEF_TRONCA/2.+(1.-COEF_TRONCA/2.)*BETA11(1)/3. )

      ! Truncation coefficient
      WRITE(1,46,ERR=923) COEF_TRONCA
      
      ! Single scattering albedo (possibly truncated)
      WRITE(1,47,ERR=923) PIZTR

      ! Information on the content of the next columns
      WRITE(1,39,ERR=923) 
      WRITE(1,48,ERR=923) OS_NB
      WRITE(1,49,ERR=923)

      ! Parameters of the phase matrix (possibly truncated)
      DO K=0,OS_NB
         WRITE(1,50,ERR=923) ALP(K),BETA11(K),GAMMA12(K),ZETA(K)
      ENDDO

 

C* File closing
C-------------------
      CLOSE(99)	! Closing of the MIE trace file
      CLOSE(88)      ! Closing of the GRANU trace file
      CLOSE(1)	! Closing of the GRANU result file
           
C* Nominal end
C-------------------
      GOTO 9999
            
C* Error case and return of status 1 to the shell
C----------------------------------------------
 925  WRITE(6,*) '  SOS_AEROSOLS : ERROR_925: '
      WRITE(6,*) '  => Error while getting SOS_ABS_ROOT variable'
      GOTO 9998

 926  WRITE(6,*) '  SOS_AEROSOLS : ERROR_926: '
      WRITE(6,*) '  => Length of $SOS_ABS_ROOT/fic is too long' 
      WRITE(6,*) '  => CTE_LENDIR has to be adjusted.' 
      GOTO 9998

  911 WRITE(6,*) '  ERROR on logfile opening for GRANU'
      GOTO 9998

  930 WRITE(6,*) '  ERROR on angles file opening for SOS_AEROSOLS'
      GOTO 9998

  931 WRITE(6,*) '  ERROR on angles file reading for SOS_AEROSOLS'
      GOTO 9998

  940 WRITE(6,*) '  ERROR on external aerosols IOP file opening'
      WRITE(6,*) '  for SOS_AEROSOLS'
      GOTO 9998

  950 WRITE(6,*) '  ERROR : The number of angles for the external phase'
      WRITE(6,*) '  functions description is higher than the defined '
      WRITE(6,*) '  size of internal tables'
      WRITE(6,*) '  --> Check the value of CTE_MAXNB_ANG_EXT in SOS.h'
      WRITE(6,*) '      compared to the value associated to NbAng '
      WRITE(6,*) '      in the external phase functions file'
      GOTO 9998

  941 WRITE(6,*) '  ERROR on external phaze functions file reading'
      WRITE(6,*) '  for SOS_AEROSOLS'
      GOTO 9998

  942 WRITE(6,*) '  ERROR on external phaze functions file reading'
      WRITE(6,*) '  for SOS_AEROSOLS : End of file encountered'
      WRITE(6,*) '  --> Check if the value associated to NbAng is '
      WRITE(6,*) '      correct, in the external phase functions file'
      GOTO 9998

  945 WRITE(6,*) '  ERROR while interpolating external phaze functions'
      WRITE(6,*) '  for SOS_AEROSOLS'
      WRITE(6,*) '  Check that the abscissas of the phaze functions'
      WRITE(6,*) '  are correct : do not use twice same values'
      GOTO 9998

  960 WRITE(6,*) '  ERROR while opening the user file for the aerosols'
      WRITE(6,*) '  mixture definition (SOS_AEROSOLS) :'
      WRITE(6,*) FICMIXTURE_AER(1:LFIC)
      GOTO 9998

  961 WRITE(6,*) '  ERROR while reading the user file for the aerosols'
      WRITE(6,*) '  mixture definition (SOS_AEROSOLS) :'
      WRITE(6,*) FICMIXTURE_AER(1:LFIC)
      GOTO 9998

  962 WRITE(6,*) '  Non expected model name in the file for the'
      WRITE(6,*) '  mixture definition (SOS_AEROSOLS) :'
      WRITE(6,*) '  --> Wrong name : ', CMODEL
      WRITE(6,*) '  in file : ', FICMIXTURE_AER(1:LFIC)
      WRITE(6,*) '  Available models are : LND, JUNGE'
      GOTO 9998

  963 WRITE(6,*) "  The summation of AOT rates is not equal to 1"
      WRITE(6,*) "  in the user's file:", FICMIXTURE_AER(1:LFIC)
      WRITE(6,*) "  The summation is:", SOM_TAUX_AOT
      WRITE(6,*) "  --> Check the AOT rate by mode"
      WRITE(6,*) "  --> Or adjust the value of the threshold"
      WRITE(6,*) "      CTE_GAP_TOLER_SUM_RATES in SOS.h"
      GOTO 9998

  913 WRITE(6,*) '  ERROR on GRANU result file opening'
      GOTO 9998

  921 WRITE(6,*) '  ERROR on logfile writing for GRANU'
      GOTO 9998

  922 WRITE(6,*) '  ERROR on Gauss angles file reading'
      GOTO 9998

  923 WRITE(6,*) '  ERROR on GRANU result file writing'
      GOTO 9998

  993 WRITE(6,*) '  ERROR on subroutine SOS_MIE'
      GOTO 9998

  994 WRITE(6,*) '  ERROR on subroutine SOS_INIT_PARAMWMO'
      GOTO 9998

  995 WRITE(6,*) '  ERROR on subroutine SOS_GRANU'
      GOTO 9998

  996 WRITE(6,*) '  ERROR on subroutine SOS_DECOMPO_LEGENDRE'
      GOTO 9998

  997 WRITE(6,*) '  ERROR on subroutine SOS_INIT_PARAMSF'
      GOTO 9998

 1006 WRITE(6,*) '  ERROR on IMOD value : out of available range.'
      GOTO 9998

 1007 WRITE(6,*) '  ERROR on TA value : wrong negative value.'
      GOTO 9998

 1008 WRITE(6,*) '  ERROR on WA value : '
      WRITE(6,*) '  --> wavelenth has to be defined in ]0.2, 4.0[ mic.'
      GOTO 9998

 1009 WRITE(6,*) '  ERROR : CTE_MIE_ALPHAMIN > ALPHAF'
      WRITE(6,*) '  or ALPHAF > 1D5 (max value) '
      WRITE(6,*) '  ALPHAF = ', ALPHAF
      WRITE(6,*) 
     &'  --> check the definition of the size distribution model'
      WRITE(6,*) '  --> check value in SOS.h and compile again'
      GOTO 9998

 1010 WRITE(6,*) '  ERROR : CTE_MIE_ALPHAMIN F6.4'
      WRITE(6,*) '          ALPHAF F9.2'
      WRITE(6,*) '  --> check SOS.h and argument value & compile again'
      GOTO 9998  

 1020 WRITE(6,*) '  ERROR : mv command'
      WRITE(6,*) '  => The command to move the temporary MIE file', 
     &           '  to the final depository is maybe too long'
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      GOTO 9998
  
10100 WRITE(6,*) '  ERROR : --> imaginary part of refractive'
      WRITE(6,*) '              indexes have to be negative'

 9998 IER=-1
      GOTO 9999          
      
       
C* Format
C---------   
   10 FORMAT(1X,F6.2,6(2X,E12.4))
   25 FORMAT(1X,I3,5(2X,E13.5))
   35 FORMAT(3X,A1,5X,A6,9X,A9,6X,A10,5X,A10,5X,A7)
  
   
   39 FORMAT(33h---------------------------------)
   40 FORMAT(38hEXTINCTION CROSS SECTION (mic^2)     :,E13.5)
   41 FORMAT(38hSCATTERING CROSS SECTION (mic^2)     :,E13.5) 
   42 FORMAT(38hASYMMETRY FACTOR (no truncation)     :,E13.5)   
   46 FORMAT(38hTRUNCATION COEFFICIENT               :,F9.5)  
   47 FORMAT(38hSINGLE SCATTERING ALBEDO (truncation):,F9.5) 
   48 FORMAT(36hPHASE MATRIX COEFFICIENTS FOR K=0 TO,I4)
   49 FORMAT(55hALPHA(K)        BETA11(K)       GAMMA12(K)      ZETA(K))  
   50 FORMAT(E15.8,3(1X,E15.8))  

  402 FORMAT(F7.2,4E14.5)
     
  222 FORMAT(I4,1X,D21.14,2X,F7.3,2X,D21.14,2X)
  
 9999 RETURN
      END   	!END OF PROGRAM SOS_AEROSOLS

 
 
 
C==============================================================================
C PROCEDURE: SOS_NOM_FICMIE
C ==========
C      This procedure generates the name of the Mie file from:
C      	- the refractive index: real part RN
C	- the refractive index: imaginary part IN
C	- the minimum size parameter ALPHAO
C	- the maximum size parameter ALPHAF
C      - the number of Gauss angles used NBMU_GAUSS
C      - the name of the user angle file FICANGLES_USER
C          (without its directory tree)
C
C      The following formats are imposed:
C         RN with format F5.3 : precision to 1.E-O3.
C         IN with format F8.5 : precision to 1.E-O5 for a negative value.
C         ALPHAO with format F6.4 : 1.E-04 <= ALPHA0 < 10
C         ALPHAF with format F8.2 : 1.E-02 <= ALPHA0 < 1.E+05
C         NBMU_GAUSS with format I4 : 1 < NBMU_GAUSS < 9999
C
C      The name of the Mie file will be:
C
C      MIE-RN-IN-ALPHAO-ALPHAF-MUnb-NomFicUser with
C      - RN with format F5.3 (e.g. 1.400)
C      - IN (opposite of the true negative value) with format F7.5 (e.g. 0.00800)
C      - ALPHAO with format F6.4 (e.g. 0.0001)
C      - ALPHAF with format F8.2 (e.g. 00500.00)
C      - MUnb the number of Gauss angles (e.g. MU40)
C      - NomFicUser the name of the user file (without directory tree) (e.g. MyAngles)
C
C      RN is converted into a character string CRN(1:4)
C      For this, Format F5.3 * 1000 => A4
C           The first character CRN(1) is the units digit of RN,
C           The second character CRN(2) is the first decimal of RN,
C           etc.
C
C
C Description of input / output parameters
C ------------------------------------------
C	NBMU_GAUSS (I4)	: (I) Number of Gauss angles used
C       FICANGLES_USER (CHARACTER*CTE_LENFIC2)	: (I) Name of the "user"
C                                                 angles file (complete path)
C	RN (double)	: (I) real part of the refractive index
C	IN (double)	: (I) imaginary part of the refractive index
C       ALPHAO (double) : (I) minimum size parameter
C       ALPHAF (double) : (I) maximum size parameter
C	FICS (CHARACTER*CTE_LENFIC1) : (O) name of the Mie file (without directory tree)
C
C Constant used:
C ------------------
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
      SUBROUTINE SOS_NOM_FICMIE(NBMU_GAUSS,FICANGLES_USER,
     &                          RN,IN,ALPHAO,ALPHAF,FICS) 

      IMPLICIT NONE
         
      DOUBLE PRECISION RN	! real part of the refractive index.
      DOUBLE PRECISION IN	! imaginary part of the refractive index.
      DOUBLE PRECISION ALPHAO   ! value of alpha from which
      				! the Mie computation starts.
      DOUBLE PRECISION ALPHAF   ! value of alpha at which the computation stops.

      CHARACTER*CTE_LENFIC2 FICANGLES_USER  ! path and name of the "user" angles file
      CHARACTER*CTE_LENFIC1 FICS	    ! name of the MIE file
          
      CHARACTER*4 CRN		! Character strings
      CHARACTER*6 CIN 
      CHARACTER*5 CALPHAO   
      CHARACTER*7 CALPHAF
      CHARACTER*4 CNBMU
      
      INTEGER*2 I
      INTEGER*4 IDEB, IDEB_CH, IFIN, IFIN_CH  ! start and end positions of a character string
      INTEGER*4 NBMU_GAUSS  ! Number of Gauss angles used
      			      	 
	    
C --------- Conversion of the input parameters into character strings
      WRITE(CRN,'(I4)') int(RN*1000)  		
      WRITE(CIN,'(I6)') int(-IN*100000)		!Positive value
      
      WRITE(CALPHAO,'(I5)') int(ALPHAO*10000)
      WRITE(CALPHAF,'(I7)') int(ALPHAF*100)
      
      WRITE(CNBMU,'(I4)') NBMU_GAUSS

C --------- Blanks are replaced by 0
      DO I=1,4
         IF( CRN(I:I).EQ.' ') CRN(I:I) = '0'
      ENDDO
      
      DO I=1,6
         IF( CIN(I:I).EQ.' ') CIN(I:I) = '0'
      ENDDO
      
      DO I=1,5
         IF( CALPHAO(I:I).EQ.' ') CALPHAO(I:I) = '0'
      ENDDO
      
      DO I=1,7
         IF( CALPHAF(I:I).EQ.' ') CALPHAF(I:I) = '0'
      ENDDO
      	 
	 
      FICS = 'MIE'
      
C -------- Building of the character string for the index RN
      FICS(4:4) = CRN(1:1)     ! Units digit
      FICS(5:5) = '.'	       ! Decimal point
      FICS(6:8) = CRN(2:4)     ! First three decimals
      FICS(9:9) = '-'
      
C -------- Building of the character string for the index IN
      FICS(10:10) = CIN(1:1)     ! Units digit
      FICS(11:11) = '.'	         ! Decimal point
      FICS(12:16) = CIN(2:6)     ! First five decimals
      FICS(17:17) = '-'	
      
C -------- Building of the character string for AlphaO
      FICS(18:18) = CALPHAO(1:1)     ! Units digit
      FICS(19:19) = '.'	             ! Decimal point
      FICS(20:23) = CALPHAO(2:5)     ! First four decimals
      FICS(24:24) = '-'	

C -------- Building of the character string for AlphaF
      FICS(25:29) = CALPHAF(1:5)     ! Units digits
      FICS(30:30) = '.'	             ! Decimal point
      FICS(31:32) = CALPHAF(6:7)     ! First two decimals
      FICS(33:33) = '-'	            

C -------- Building of the character string for MUgauss
C --------Structured to have NBMU_GAUSS on the left (without "blanks")
      
      FICS(34:35) = 'MU'
      
      IF (NBMU_GAUSS.LT.10) THEN
          FICS(36:36) = CNBMU(4:4)	!MUn
      ELSE 
          IF (NBMU_GAUSS.LT.100) THEN
              FICS(36:37) = CNBMU(3:4)  !MUnn
	  ELSE 
	     IF (NBMU_GAUSS.LT.1000) THEN
                 FICS(36:38) = CNBMU(2:4) !MUnnn
	     ELSE !NBMU_GAUSS < 10000
	     	 FICS(36:39) = CNBMU(1:4) !MUnnnn
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

C       Addition of the name of the "user" angle file to the name of the MIE file
         IFIN=INDEX(FICS,' ')  
	  FICS(IFIN:IFIN) = '-'
	  
	  IFIN=IFIN+1
	  FICS(IFIN:IFIN+IFIN_CH-IDEB_CH+1) 
     &	    = FICANGLES_USER(IDEB_CH:IFIN_CH)
	  
      ENDIF 
     	
	
      END !SOS_NOM_FICMIE
 
 
 
 
C==============================================================================
C PROCEDURE: SOS_INIT_PARAMWMO
C ==========
C      This procedure initialises all the parameters of each
C      WMO component for the simulation wavelength:
C      	- the modal radius
C	- the variance
C	- volume concentration
C	- real part of the refractive index
C	- imaginary part
C
C      Linear interpolation of the indices for the input wavelength.
C
C      Forces the result indices to the format F5.3 for the real part
C      and F8.5 for the imaginary part (negative value).
C
C Description of input / output parameters
C ------------------------------------------
C	FICDATAWMO (CHARACTER*CTE_LENFIC2)	: (I) path and name of the file Data_WMO
C	WA (double)	: (I) selected wavelength
C       VARGRANU1I(CTE_NB_WMO) (double) : (O) modal radius of each component
C       VARGRANU2I(CTE_NB_WMO) (double) : (O) decimal log of the variance of
C					  each component
C	MR(CTE_NB_WMO) (double)	: (O) real part of the refractive index
C       MI(CTE_NB_WMO) (double) : (O) imaginary part of the refractive index
C	V(CTE_NB_WMO) (double)	: (O) volume concentration micron3/cm3
C	IER (I4)	: (O) error code = 0 if no error, =1 otherwise
C
C 
C Description of the files used
C ---------------------------------
C   -->Content of the WMO data file
C      Ascii file, record format:
C      1st line : 4(1X,F9.5)
C                   value of the modal radius (mic)
C		    for the components DL,WS,OC,SO
C      2nd line : 4(1X,F9.5)
C                   decimal log of the variance
C		    for the components DL,WS,OC,SO
C      3rd line : 2(1X,F9.5,1X,E12.5)
C                   volume concentration (mic**3/cm**3)
C		    for the components DL,WS,OC,SO
C      4th line to the end : 9(1X,F9.5)
C                            wavelength, real and imaginary part
C			     of the refractive index of each component.
C
C
C Constants used:
C --------------------
C     CTE_LENFIC2 : Length of file names with directory tree.
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------    
C     -Opening of a file
C       Error if the WMO data file does not exist
C     -Reading in a file
C
C     Message displayed on screen, the routine stops its computations and
C     returns status -1 to the calling program
C
C
C  Subroutines used:
C --------------------------
C  This program calls the routine: SOS_INTERPOL
C
C
C  MOD:VERSION 2.0: 02/03/2016
C      Conversion of the variance of the WMO data file (in decimal log)
C      into a variance consistent with the classical definition of the LND: VAR_WMO * LOG(10)
C
C==============================================================================
      SUBROUTINE SOS_INIT_PARAMWMO(FICDATAWMO,WA,VARGRANU1I,
     &				   VARGRANU2I,MR,MI,V,IER)

      IMPLICIT NONE
      
C AEROSOLS constants
C -------------------

      INTEGER CTE_NB_WMO	!Number of WMO aerosol types used: DL,WS,OC,SO
      PARAMETER (CTE_NB_WMO=4)

      INTEGER CTE_NB_AERO	!Dimension of the result arrays
      PARAMETER (CTE_NB_AERO=5)	!(4 WMO components, 5 Shettle & Fenn components)
      
            
C* Definition of variables
C*-----------------------------------------

      CHARACTER*CTE_LENFIC2 FICDATAWMO	! path and name of the file Data_WMO
      
      DOUBLE PRECISION WA	! Selected wavelength in microns.
      DOUBLE PRECISION WA1	! Wavelength in microns read in the
      DOUBLE PRECISION WA2	! file Data_WMO.
      
      DOUBLE PRECISION VARGRANU1I(CTE_NB_AERO)   ! Size distribution parameter 1
      						       ! for each component model.
      DOUBLE PRECISION VARGRANU2I(CTE_NB_AERO)   ! Size distribution parameter 2
      						       ! for each component model.
      
      DOUBLE PRECISION MR(CTE_NB_AERO)   ! real part of the refractive index
      					      ! for each aerosol model.
      DOUBLE PRECISION MI(CTE_NB_AERO)   ! imaginary part of the refractive index.
      					      ! for each aerosol model.
      DOUBLE PRECISION MR1(CTE_NB_WMO)  !
      DOUBLE PRECISION MI1(CTE_NB_WMO)  ! used to store
      DOUBLE PRECISION MR2(CTE_NB_WMO)  ! intermediate values
      DOUBLE PRECISION MI2(CTE_NB_WMO)  !
      
      DOUBLE PRECISION V(CTE_NB_WMO)	!Volume concentration of the components
      
      INTEGER*4 I	   !Index.
      INTEGER*4 IER        ! error code =0 if no error, =-1 otherwise
      
C* Initialisation
C  --------------
       IER=0	 
       
       DO I=1,CTE_NB_WMO
          VARGRANU1I(I)=0.
          VARGRANU2I(I)=0.
	   MR(I)=0.
    	   MI(I)=0.
	   V(I)=0.
       ENDDO

C* Opening of the file containing the WMO data
C  ---------------------------------------------------	 
       OPEN(8,FILE=FICDATAWMO,STATUS='OLD',ERR=991)

C* Reading of the modal radius
C  ---------------------- 	 
       READ(8,333,err=992) VARGRANU1I(1),VARGRANU1I(2),
     &                     VARGRANU1I(3),VARGRANU1I(4)
     	  
C* Reading of the decimal log of the variance
C  ------------------------------------- 	 
       READ(8,333,err=992) VARGRANU2I(1),VARGRANU2I(2),
     &                     VARGRANU2I(3),VARGRANU2I(4)      

C* Conversion to a variance for LND consistent with the classical definition
C  ------------------------------------- 	 
       DO I=1,4
          VARGRANU2I(I) = VARGRANU2I(I) * DLOG(10.D+00)
       ENDDO

C* Reading of the volume concentration
C  ------------------------------------- 	 
       READ(8,444,err=992) V(1),V(2),V(3),V(4) 

C* Determination of the refractive indices by linear interpolation
C  ------------------------------------------------------------------

C* Initialisation with the first values of the file
C  ---------------------------------------------------- 
       READ(8,555,END=999,err=992)WA1,MR1(1),MI1(1),MR1(2),MI1(2),
     &                            MR1(3),MI1(3),MR1(4),MI1(4)

C* For each line of the file, search for the wavelength
C  just above the selected WA
C  --------------------------------------------------------------
  882  READ(8,555,END=999,err=992)WA2,MR2(1),MI2(1),MR2(2),MI2(2),
     &                            MR2(3),MI2(3),MR2(4),MI2(4)
     
       ! WA1<=WA<WA2, interpolation of the indices
       IF((WA1.LE.WA).AND.(WA.LE.WA2)) THEN
       
	  DO I=1,CTE_NB_WMO	  
	     CALL SOS_INTERPOL(MR1(I),MR2(I),WA1,WA2,WA,MR(I))
	     CALL SOS_INTERPOL(MI1(I),MI2(I),WA1,WA2,WA,MI(I))
	     
C******** Conversion to the formats F5.3 and F8.3
	     MR(I) =  DNINT( MR(I)*1000.D+00) / 1000.D+00
	     MI(I) = -DNINT(-MI(I)*100000.D+00) / 100000.D+00	    
	  ENDDO
	  GOTO 881
	       
       ELSE
       	 
          WA1=WA2	! Re-initialisation: line1 receives line2
	  DO I=1,CTE_NB_WMO
	     MR1(I)=MR2(I)
	     MI1(I)=MI2(I)
     	  ENDDO
	  GOTO 882	! next line
	       
       ENDIF

  999  CONTINUE
  	    
  881  CONTINUE
       CLOSE(8)

C* Nominal end
C-------------------
       GOTO 9999
              
C* Error case and return of status -1 to the calling program
C----------------------------------------------------------
  991  WRITE(6,*) 'Error while opening the WMO datafile'
       WRITE(6,*) 'FICDATAWMO :',FICDATAWMO
       IER=-1
       GOTO 9999
  992  WRITE(6,*) 'Error while reading the WMO datafile'
       IER=-1
       GOTO 9999

C* Format
C---------   
  333  FORMAT(4(1X,F9.5))
  444  FORMAT(1X,F9.5,1X,E12.5,1X,F9.5,1X,E12.5)
  555  FORMAT(9(1X,F9.5))

       
 9999  CONTINUE
 
       RETURN
       END	!END OF PROCEDURE SOS_INIT_PARAMWMO


C==============================================================================
C PROCEDURE: SOS_INIT_PARAMSF
C ==========
C      This procedure initialises all the parameters of each
C      component of the Shettle and Fenn models for the wavelength
C      and the relative humidity selected:
C      	- the modal radius
C	- the variance
C	- real part of the refractive index
C	- imaginary part
C
C Description of parameters
C --------------------------
C	DIRFIC (CHARACTER*CTE_LENDIR)  : (I) access path to the data file
C	WA (double)	      : (I) selected wavelength
C	RH (double)	      : (I) relative humidity
C       VARGRANU1I(CTE_NB_AERO) (double) : (O) modal radius of each component
C       VARGRANU2I(CTE_NB_AERO) (double) : (O) decimal log of the variance of
C					  each component
C	MR(CTE_NB_AERO) (double) : (O) real part of the refractive index
C       MI(CTE_NB_AERO) (double) : (O) imaginary part of the refractive index
C	IER (I4)	: (O) error code = 0 if no error, =1 otherwise
C
C 
C Description of the files used
C ---------------------------------
C   -->Content of the file of the size distribution parameters of the
C      Shettle and Fenn components:
C      Ascii file, record format:
C      1st line : 5(1X,F9.5)
C                   decimal log of the variance
C		    for the components SR,LR,SU,LU and OM
C      2nd line to the end : 1X,F5.2,5(1X,F9.5)
C			     % of relative humidity, modal radius
C			     of each component.
C
C   -->Content of the files IRefrac_XX (with XX = SR, LR ...)
C      Ascii file, record format 17(1X,F9.5)
C      for each line: wavelength, real part and imaginary part
C			  of the refractive index for 8 relative
C			  humidities: 0%, 50%, 70%, 80%, 90%, 95%, 98%, 99%.
C
C
C Constants used:
C --------------------
C    CTE_LENDIR : Length of directory names.
C     CTE_LENFIC2 : Length of file names with directory tree.
C    CTE_AER_DATASF  : Name of the file containing the S&F modal radius and log variance.
C    CTE_AER_XX_SF   : Name of the files containing the S&F refractive indices
C		       for each component XX = SR, LR, SU, LU, OM.
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------    
C     -Opening of a file
C       Error if one of the data files does not exist
C     -Reading in a file
C
C     Message displayed on screen, the routine stops its computations and
C     returns status -1 to the calling program
C
C  Subroutines used:
C --------------------------
C  This program calls the routine: SOS_INTERPOL
C
C
C  MOD:VERSION 2.0: 02/03/2016
C      Conversion of the variance of the S&F data file (in decimal log)
C      into a variance consistent with the classical definition of the LND: VAR_WMO * LOG(10)
C
C==============================================================================
      SUBROUTINE SOS_INIT_PARAMSF(DIRFIC,WA,RH,VARGRANU1I,
     &                         	  VARGRANU2I,MR,MI,IER)

      IMPLICIT NONE
      
C AEROSOLS constants
C -------------------

      INTEGER CTE_NB_SF	        !Number of components of the Shettle & Fenn models
      PARAMETER (CTE_NB_SF=5)

      INTEGER CTE_NB_RH		!Number of relative humidity values
      PARAMETER (CTE_NB_RH=8)
 
      INTEGER CTE_NB_AERO	!Dimension of the result arrays
      PARAMETER (CTE_NB_AERO=5)	!(4 WMO components, 5 Shettle & Fenn components)
      
            
C* Definition of variables
C*-----------------------------------------

      CHARACTER*CTE_LENDIR DIRFIC	! Access path to the S&F files
      CHARACTER*CTE_LENFIC2 FICDATASF	! Name of the file Data_SF
      CHARACTER*CTE_LENFIC2 FIC_SR_SF 	! File IRefrac_SR
      CHARACTER*CTE_LENFIC2 FIC_LR_SF 	! File IRefrac_LR
      CHARACTER*CTE_LENFIC2 FIC_SU_SF 	! File IRefrac_SU
      CHARACTER*CTE_LENFIC2 FIC_LU_SF 	! File IRefrac_LU
      CHARACTER*CTE_LENFIC2 FIC_OM_SF 	! File IRefrac_OM
      
      DOUBLE PRECISION WA	! Selected wavelength in microns.
      DOUBLE PRECISION WA1	! Wavelength in microns read in the files
      DOUBLE PRECISION WA2	
      
      DOUBLE PRECISION RH	! Relative humidity (%)
      DOUBLE PRECISION RH1	! Relative humidity read in the files
      DOUBLE PRECISION RH2	!
      
      DOUBLE PRECISION VARGRANU1I(CTE_NB_AERO)  ! Size distribution parameter 1
      						      ! for each component model.
      DOUBLE PRECISION VARGRANU2I(CTE_NB_AERO)  ! Size distribution parameter 2
      						      ! for each component model.
      
      DOUBLE PRECISION MR(CTE_NB_AERO)  ! real part of the refractive index
      					     ! for each aerosol model.
      DOUBLE PRECISION MI(CTE_NB_AERO)  ! imaginary part of the refractive index.
      					     ! for each aerosol model.
      
      DOUBLE PRECISION MR_CPT1		!
      DOUBLE PRECISION MR_CPT2		! Intermediate values
      DOUBLE PRECISION MI_CPT1		!
      DOUBLE PRECISION MI_CPT2		!
      
      DOUBLE PRECISION RM1(CTE_NB_SF)	!
      DOUBLE PRECISION RM2(CTE_NB_SF)	! used to store
      DOUBLE PRECISION MR1(CTE_NB_RH)   ! intermediate values
      DOUBLE PRECISION MR2(CTE_NB_RH)   ! 
      DOUBLE PRECISION MI1(CTE_NB_RH)   ! 
      DOUBLE PRECISION MI2(CTE_NB_RH)   !
	
      INTEGER*4 I,H	   ! Index
      INTEGER*4 CPT	   ! Counter
      INTEGER*4 IFIN	   ! End flag of a character string
      
      INTEGER*4 IER        ! error code =0 if no error, =-1 otherwise
      
C* Initialisation
C  --------------
       IER=0
       
       DO I=1,CTE_NB_SF
          VARGRANU1I(I)=0.
          VARGRANU2I(I)=0.	  
	   MR(I)=0.
    	   MI(I)=0. 
       ENDDO

       
       ! Initialisation of the names of the different files containing
       ! the data of the S&F model
       
       IFIN=INDEX(DIRFIC,' ')
       IFIN=IFIN-1
       IF (IFIN.LE.0) IFIN=CTE_LENDIR
      
       FICDATASF = DIRFIC(1:IFIN)//'/'//CTE_AER_DATASF
       FIC_SR_SF = DIRFIC(1:IFIN)//'/'//CTE_AER_SR_SF
       FIC_LR_SF = DIRFIC(1:IFIN)//'/'//CTE_AER_LR_SF
       FIC_SU_SF = DIRFIC(1:IFIN)//'/'//CTE_AER_SU_SF
       FIC_LU_SF = DIRFIC(1:IFIN)//'/'//CTE_AER_LU_SF
       FIC_OM_SF = DIRFIC(1:IFIN)//'/'//CTE_AER_OM_SF

C* Opening of the file containing the size distribution parameters
C* of the S&F components
C  --------------------------------------------------------------
       OPEN(9,FILE=FICDATASF,STATUS='OLD',ERR=991)
   	  
C* Reading of the decimal log of the variance
C  ------------------------------------- 	 
       READ(9,222,err=993) VARGRANU2I(1),VARGRANU2I(2),VARGRANU2I(3),
     &                     VARGRANU2I(4),VARGRANU2I(5) 

C* Conversion to a variance for LND consistent with the classical definition
C  ------------------------------------- 	 
       DO I=1,5
          VARGRANU2I(I) = VARGRANU2I(I) * DLOG(10.D+00)
       ENDDO
                    
C* Reading of the modal radius versus the relative humidity
C  --------------------------------------------------------- 	 
       READ(9,333,err=993) RH1,RM1(1),RM1(2),
     &                      RM1(3),RM1(4),RM1(5)
        
	CPT=1
	
	! Case RH = 0%
	IF(RH1.EQ.RH) THEN
	
	   DO I=1,CTE_NB_SF
	      VARGRANU1I(I)=RM1(I)
	   ENDDO
	   GOTO 884
	   
	ENDIF

C* Search for the RH values around the selected value
C  ------------------------------------------------------------
  883   READ(9,333,err=993) RH2,RM2(1),RM2(2),
     &                      RM2(3),RM2(4),RM2(5)
     
        CPT=CPT+1	! Counter giving the position of
			! the value RH2 just above RH
	
	! Case RH1<RH<=RH2, interpolation of the refractive indices
	IF((RH1.LT.RH).AND.(RH.LE.RH2)) THEN
	
	   DO I=1,CTE_NB_SF	  
	      CALL SOS_INTERPOL(RM1(I),RM2(I),RH1,RH2,RH,VARGRANU1I(I))
     	   ENDDO
	   GOTO 884
	       
        ELSE	! Re-initialisation of the values
       
          RH1=RH2
	  DO I=1,CTE_NB_SF
	     RM1(I)=RM2(I)
     	  ENDDO
	  GOTO 883	! next line in the file
	  
       ENDIF
       
  884  CONTINUE


C* Opening of the files containing the refractive index for
C  each component of the S&F models
C  -----------------------------------------------------------
       OPEN(11,FILE=FIC_SR_SF,STATUS='OLD',ERR=992)
       OPEN(12,FILE=FIC_LR_SF,STATUS='OLD',ERR=992)
       OPEN(13,FILE=FIC_SU_SF,STATUS='OLD',ERR=992)
       OPEN(14,FILE=FIC_LU_SF,STATUS='OLD',ERR=992)
       OPEN(15,FILE=FIC_OM_SF,STATUS='OLD',ERR=992)
  
       DO I=1,CTE_NB_SF

C* Determination of the refractive indices by linear interpolation
C  ------------------------------------------------------------------

C* Initialisation with the first values of the file
C  ----------------------------------------------------
          READ((I+10),555,END=999,err=993)WA1,MR1(1),MI1(1),MR1(2),
     &                             MI1(2),MR1(3),MI1(3),MR1(4),MI1(4),
     &				   MR1(5),MI1(5),MR1(6),MI1(6),MR1(7),
     &				   MI1(7),MR1(8),MI1(8)

C* For each line of the file, search for the wavelength
C  just above the selected one
C  --------------------------------------------------------------
  885     READ((I+10),555,END=999,err=993)WA2,MR2(1),MI2(1),MR2(2),
     &                             MI2(2),MR2(3),MI2(3),MR2(4),MI2(4),
     &				   MR2(5),MI2(5),MR2(6),MI2(6),MR2(7),
     &				   MI2(7),MR2(8),MI2(8)
     
          ! Case WA1<=WA<=WA2, interpolation of the data
	  IF((WA1.LE.WA).AND.(WA.LE.WA2)) THEN
	  
	     ! Case RH = 0%, interpolation in WA but not in RH
	     IF(CPT.EQ.1) THEN
	        CALL SOS_INTERPOL(MR1(CPT),MR2(CPT),WA1,WA2,WA,MR(I))
		 CALL SOS_INTERPOL(MI1(CPT),MI2(CPT),WA1,WA2,WA,MI(I))
	     
	     ELSE
	      ! Interpolations for the real part of the refractive index
		! Interpolation in WA for RH1(CPT-1) and RH2(CPT) around RH
		 CALL SOS_INTERPOL(MR1(CPT-1),MR2(CPT-1),
     &		                  WA1,WA2,WA,MR_CPT1)
	        CALL SOS_INTERPOL(MR1(CPT),MR2(CPT),WA1,WA2,WA,MR_CPT2)
		! Interpolation in RH
	        CALL SOS_INTERPOL(MR_CPT1,MR_CPT2,RH1,RH2,RH,MR(I))
	      
	      ! Interpolations for the imaginary part
	        CALL SOS_INTERPOL(MI1(CPT-1),MI2(CPT-1),
     &				  WA1,WA2,WA,MI_CPT1)
	        CALL SOS_INTERPOL(MI1(CPT),MI2(CPT),WA1,WA2,WA,MI_CPT2)
	        CALL SOS_INTERPOL(MI_CPT1,MI_CPT2,RH1,RH2,RH,MI(I))
	     ENDIF
		
          ELSE
       
            WA1=WA2	! Re-initialisation
	     DO H=1,CTE_NB_RH
	        MR1(H)=MR2(H)
	        MI1(H)=MI2(H)
     	     ENDDO
	     GOTO 885	! Next line of the file
	       
          ENDIF
	  
C* Conversion to the formats F5.3 and F8.3
	  MR(I) =  DNINT( MR(I)*1000.D+00) / 1000.D+00
	  MI(I) = -DNINT(-MI(I)*100000.D+00) / 100000.D+00
	  	    
       ENDDO
       
  999  CONTINUE
  	    
       CLOSE(9)
       CLOSE(11)
       CLOSE(12)
       CLOSE(13)
       CLOSE(14)
       CLOSE(15)

C* Nominal end
C-------------------
       GOTO 9999
              
C* Error case and return of status -1 to the calling program
C----------------------------------------------------------
  991  WRITE(6,*) 
     &     'Error while opening the Shettle&Fenn component  datafile'
       IER=-1
       GOTO 9999
  992  WRITE(6,*) 
     &	  'Error while opening the S&F  refractive indexes datafile'
       IER=-1
       GOTO 9999     
  993  WRITE(6,*) 
     &    'Error while reading the Shettle&Fenn component  datafile'
       IER=-1
       GOTO 9999
  994  WRITE(6,*) 
     &	  'Error while reading the S&F  refractive indexes datafile'
       IER=-1
       GOTO 9999
       
C* Format
C---------   
  222  FORMAT(5(1X,F9.5))
  333  FORMAT(1X,F5.2,5(1X,F9.5))
  555  FORMAT(17(1X,F9.5))

       
 9999  CONTINUE
 
       RETURN
       END	!END OF PROCEDURE SOS_INIT_PARAMSF




C==============================================================================
C PROCEDURE: SOS_INTERPOL
C ==========
C      This procedure linearly interpolates in X a function f
C      whose values (X1,Y1) and (X2,Y2) are known.
C
C Description of parameters
C --------------------------
C	X1 (double) : (I)
C	X2 (double) : (I)
C	Y1 (double) : (I) Value f(X1)
C	Y2 (double) : (I) Value f(X2)
C	X (double)  : (I) Abscissa at which the value f(X) is interpolated
C   Y (double)  : (O) Interpolated value f(X)
C 
C
C==============================================================================
      SUBROUTINE SOS_INTERPOL(Y1,Y2,X1,X2,X,Y)

      IMPLICIT NONE

C* Definition of variables
C*-----------------------------------------

      DOUBLE PRECISION Y	
      DOUBLE PRECISION Y1
      DOUBLE PRECISION Y2
      DOUBLE PRECISION X
      DOUBLE PRECISION X1
      DOUBLE PRECISION X2      
      
C* Computation of the interpolated parameter
C  -----------------------------
      Y = ((Y2-Y1)/(X2-X1))*(X-X2) + Y2

      END	!END OF PROCEDURE SOS_INTERPOL


            
C==============================================================================
C PROCEDURE: SOS_DECOMPO_LEGENDRE
C ==========
C       This procedure computes the parameters alpha, beta, gamma delta and zeta
C	of the decomposition of the phase function in Legendre polynomials.
C	Possible application of the truncation.
C
C       The phase functions are stored in arrays defined up to
C       CTE_MIE_NBMU_MAX but are used up to MIE_NBMU (actual angles).
C
C Description of parameters
C --------------------------
C	ITRONC (I4) : (I/O) Truncation option (1 to apply it).
C                           Set to 0 as output if the truncation is too small.
C	TRACE (logical) : (I)  =true if writing in the trace file
C      MIE_NBMU (I4) : (I) Number of actually useful (positive) angles
C	 XMU(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX) (double) : (I) Cosines of the Gauss angles.
C	XHR(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX) (double) : (I) Weights of the Gauss angles.
C	P11(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX) (double)  : (I) Phase function computed
C					         from the functions of the MIE theory.
C	P12(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX) (double)  : (I) Element P12 of the phase matrix.
C	P33(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	(double) : (I) Element P33 of the phase matrix.
C	TTT(-CTE_MIE_NBMU_MAX:CTE_MIE_NBM_MAXU) (double) : (O) Copy of the phase function (P11)
C					 	               before truncation.
C	COEF_TRONCA (double)	     : (O) Truncation coefficient.
C	Z1 (double)		     : (O) Value of BETA11(0) before normalisation.
C	ALP(0:CTE_OS_NB_MAX) (double)     : (O) Parameters alpha(k).
C	BETA11(0:CTE_OS_NB_MAX) (double)  : (O) Parameters beta11(k).
C	BETA22(0:CTE_OS_NB_MAX) (double)  : (O) Parameters beta22(k).
C	GAMMA12(0:CTE_OS_NB_MAX) (double) : (O) Parameters gamma12(k).
C	ZETA(0:CTE_OS_NB_MAX) (double)    : (O) Parameters zeta(k).
C	DELTA33(0:CTE_OS_NB_MAX) (double)  : (O) Parameters delta33(k).
C	IER (I4)		     : (O) Error code =0 if no error, =1 otherwise.
C
C 
C Description of the files used
C ---------------------------------
C None
C
C Constants used:
C --------------------
C    CTE_MIE_NBMU_MAX : Maximum number of positive angles.
C    CTE_OS_NB_MAX :  Maximum value of the limit order of the expansion in
C                     Legendre polynomials.
C    SOS_PH_SEUIL_TRONCA : Threshold for the truncation.
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C     -Writing in a file
C
C     Message displayed on screen, the routine stops its computations and
C     returns status -1 to the calling program
C
C==============================================================================
      SUBROUTINE SOS_DECOMPO_LEGENDRE(ITRONC,TRACE,MIE_NBMU,XMU,XHR,
     &				          OS_NB,P11,TTT,P12,P22,P33,
     &				          COEF_TRONCA,Z1,
     &				          ALP,BETA11,BETA22,
     &                                GAMMA12,DELTA33,ZETA,IER)

      IMPLICIT NONE

C* Definition of variables
C*-----------------------------------------
      DOUBLE PRECISION XMU(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)  
                                                  !Cosines of the Gauss angles.
      DOUBLE PRECISION XHR(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)  
                                                  !Weights of the Gauss angles.
      
      DOUBLE PRECISION P11(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Phase function computed from the functions of
		       !the MIE theory.
      DOUBLE PRECISION P12(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Element P12 of the phase matrix.
      DOUBLE PRECISION P22(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Element P22 of the phase matrix.
      DOUBLE PRECISION P33(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Element P33 of the phase matrix.

      DOUBLE PRECISION TTT(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Copy of the phase function (P11) before truncation.
		       
		
      DOUBLE PRECISION PL(-1:CTE_OS_NB_MAX+1) 
                       !Legendre polynomials of order L computed for one angle.
      DOUBLE PRECISION POL(0:CTE_OS_NB_MAX+1) 
                       !Legendre function of order S=2 and variable order L
		       !(computed for one angle).
		       			
      DOUBLE PRECISION D	!Variable for the computation of POL.
      DOUBLE PRECISION E	!Variable for the computation of POL.
      		       
      DOUBLE PRECISION X
      DOUBLE PRECISION XRMU	!Value used from the array of Gauss angles.
				
      DOUBLE PRECISION XX	
      DOUBLE PRECISION XXX

      DOUBLE PRECISION COEF_TRONCA	!Truncation coefficient.
      DOUBLE PRECISION AA	!Slope of the truncation line.
      DOUBLE PRECISION X1	
      DOUBLE PRECISION X2


      DOUBLE PRECISION CO1	!Variable for the computation of ALP(k) and ZETA(k).
      DOUBLE PRECISION CO2	!Variable for the computation of ALP(k) and ZETA(k).
      DOUBLE PRECISION CO3	!Variable for the computation of ALP(k) and ZETA(k).

      DOUBLE PRECISION SOM1  	!Variable for the computation of ALP(k) and ZETA(k).
      DOUBLE PRECISION SOM2	!Variable for the computation of ALP(k) and ZETA(k).
      DOUBLE PRECISION SOM3	!Variable for the computation of ALP(k) and ZETA(k).
      DOUBLE PRECISION SOM4	!Variable for the computation of ALP(k) and ZETA(k).
      
      DOUBLE PRECISION Z1	!Value of BETA11(0) before normalisation.

      DOUBLE PRECISION ALP(0:CTE_OS_NB_MAX)	! Parameters alpha(k).
      DOUBLE PRECISION BETA11(0:CTE_OS_NB_MAX)	! Parameters beta11(k).
      DOUBLE PRECISION BETA22(0:CTE_OS_NB_MAX)	! Parameters beta22(k).
      DOUBLE PRECISION GAMMA12(0:CTE_OS_NB_MAX)	! Parameters gamma12(k).
      DOUBLE PRECISION ZETA(0:CTE_OS_NB_MAX)	! Parameters zeta(k).
      DOUBLE PRECISION DELTA33(0:CTE_OS_NB_MAX)	! Parameters delta33(k).
                        
	!   Bounding angles for the truncation (for conversion to double precision)
	!-----------------------------------
      DOUBLE PRECISION MU1_TRONCA ! cosines of the angles Theta1 and Theta2 for the definition
      DOUBLE PRECISION MU2_TRONCA ! of the segment truncating the phase function

      INTEGER*4 MIE_NBMU   ! Number of actually useful (positive) angles
      INTEGER*4 ITRONC	   ! Truncation option (1 to apply it).
      INTEGER*4 OS_NB 	   ! Order of the expansions of the phase functions
			   ! in Legendre polynomials.
      INTEGER*4 K          ! Array index.
      INTEGER*4 KK	   ! Index of the array of Gauss angles for the mu
      			   ! closest to 0.94 from above.
      INTEGER*4 NN
      INTEGER*4 MM				
      INTEGER*4 I	   ! Index.
      INTEGER*4 J	   ! Index.
      INTEGER*4 IER	   ! Error code =0 if no error, =1 otherwise.

      LOGICAL TESS	   !Variable for the computation of the truncation.
      LOGICAL TRACE	   !=true if writing in the trace file
    



C* Copy of the phase function in intensity before the truncation.
C-----------------------------------------------------------------
      DO 55 J=-MIE_NBMU,MIE_NBMU
         TTT(J)=P11(J)
  55  CONTINUE 
      
C==============
C* TRUNCATION
C============== 
      IF(ITRONC.EQ.0) GO TO 59

      MU1_TRONCA = CTE_AER_MU1_TRONCA	
      MU2_TRONCA = CTE_AER_MU2_TRONCA


C* Truncation of the phase function in intensity for the angles
C* below 19.95 degrees (mu = 0.94).
C------------------------------------------------------------------------------

C* Determination of the Gauss angles just below 36.87 degrees (mu=0.8)
C* and 19.95 degrees (mu = 0.94)
C* Determination of the Gauss angles just below 36.87 degrees (mu=0.8 : MU1_TRONCA)
C* and 19.95 degrees (mu = 0.94 : MU2_TRONCA), for instance, typically for aerosols.
C--------------------------------------------------------------------
      TESS=.TRUE.
      DO 56 J=1,MIE_NBMU
         IF (.NOT.(TESS.AND.(XMU(J).GT.MU1_TRONCA)
     &	     .AND.(XHR(J).NE.0.))) GO TO 56
         K=J-1
         TESS=.FALSE.
  56  CONTINUE

      TESS=.TRUE.
      DO 57 J=1,MIE_NBMU
         IF (.NOT.(TESS.AND.(XMU(J).GT.MU2_TRONCA)
     &	     .AND.(XHR(J).NE.0.)))  GO TO 57
         KK=J-1
         TESS=.FALSE.
  57  CONTINUE

      IF (TRACE) THEN
         WRITE(88,*,err=921) 
     &	    'Truncation applied for a segment defined between :'
         WRITE(88,*,err=921) 'Mu1 :',MU1_TRONCA, 
     &	       ' --> Theta1 = ', DACOS(MU1_TRONCA)*180./INCTE_PI
         WRITE(88,*,err=921) 'Mu1_effective :',XMU(K), 
     &	       ' --> Theta1_eff = ', DACOS(XMU(K))*180./INCTE_PI
	  WRITE(88,*,err=921) 'Mu2 :',MU2_TRONCA, 
     &	       ' --> Theta2 = ', DACOS(MU2_TRONCA)*180./INCTE_PI
         WRITE(88,*,err=921) 'Mu2_effective :',XMU(KK), 
     &	       ' --> Theta2_eff = ', DACOS(XMU(KK))*180./INCTE_PI
	 WRITE(88,*,err=921) ""
      ENDIF

C* Determination of the slope of the line "log10(p(cos(teta))) = AA * teta + BB"
C* going through the values of log10 (p) for the Gauss mu above.
C--------------------------------------------------------------------------
      AA=
     &(DLOG10(P11(KK))-DLOG10(P11(K)))/(DACOS(XMU(KK))-DACOS(XMU(K)))

C* Application of the truncation
C-------------------------------
      X1=DLOG10(P11(KK))
      X2=DACOS(XMU(KK))
      
      DO 58 J=KK+1,MIE_NBMU
         IF(DABS(XMU(J)-1.).LE.1D-08) COEF_TRONCA=X1-AA*X2
         COEF_TRONCA=X1+AA*(DACOS(XMU(J))-X2)
         P11(J)=10**COEF_TRONCA
  58  CONTINUE

  59  CONTINUE


C* Computation of the parameters BETA11(k) of the decomposition of the
C* (possibly truncated) phase function in Legendre polynomials.
C------------------------------------------------------------------
      DO 61 J=-MIE_NBMU,MIE_NBMU
         IF(J.EQ.0) GO TO 61
         X=P11(J)*XHR(J)
         XRMU=XMU(J)

C* Computation of the Legendre polynomials of order K, for the fixed value of mu.
C------------------------------------------------------------------------
         DO 60 K=0,OS_NB
            PL(-1)=0.
            PL(0)=1.
            PL(K+1)=((2*K+1.)*XRMU*PL(K)-K*PL(K-1))/(K+1.)
            BETA11(K)=BETA11(K)+X*PL(K)
  60     CONTINUE
  61  CONTINUE

      DO 62 K=0,OS_NB
         BETA11(K)=(2*K+1)*BETA11(K)*.5
  62  CONTINUE



C* Computation of the truncation coefficient
C--------------------------------------
      IF (ITRONC.EQ.1) THEN
         COEF_TRONCA=2*(1-BETA11(0))
      ELSE
         COEF_TRONCA=0.0
      ENDIF

C* If the truncation coefficient is below the threshold Seuil_Tronca
C* the truncation is cancelled.
C----------------------------------------------------------------------	
      IF((ITRONC.EQ.1).AND.(COEF_TRONCA.LT.CTE_PH_SEUIL_TRONCA)) THEN
      
         IF (TRACE) THEN
            WRITE(88,*,err=921) 'Truncation coefficient = ',
     &      COEF_TRONCA,' Value too weak : TRUNCATION IS CANCELED'
   
	 ENDIF
	 
C* Copy back of the phase function in intensity before the truncation
C* and re-initialisation of the parameters of the decomposition of the matrix
C* computations
C-------------------------------------------------------------------
      	 DO 66 J=-MIE_NBMU,MIE_NBMU
            P11(J)=TTT(J)
  66  	 CONTINUE
   
         ITRONC=0
	 DO 67 K=0,OS_NB
             BETA11(K)=0.D+00
	     BETA22(K)=0.D+00
             GAMMA12(K)=0.D+00
             DELTA33(K)=0.D+00
             ALP(K)=0.D+00
             ZETA(K)=0.D+00
 67      CONTINUE
  
	 GOTO 59 
      ENDIF    !End of test whether the truncation is too small.
  

C* Computation of the parameters GAMMA12(k), BETA22(K) and DELTA33(k) of the decomposition of the
C* polarised phase functions P12, P22 and P33 in Legendre polynomials or functions.
C* Application of the possible truncation through the ratio of the truncated
C* phase function in intensity to the initial one: P11(J)/TTT(J).
C-----------------------------------------------------------------------------
      DO 11 J=-MIE_NBMU,MIE_NBMU       
	 IF(J.EQ.0)  GO TO 11
       
         XRMU=XMU(J)
	 
         PL(-1)=0.
         PL(0)=1.
	  POL(0)=0.
         POL(1)=0.
         POL(2)=3.*(1.-XRMU**2)/2./DSQRT(6.0D+00)
	 
	  XXX=XHR(J)*P12(J)*P11(J)/TTT(J)
         DO 12 K=2,OS_NB
            D=(2.*K+1.)/DSQRT(1.D+00*(K+3.)*(K-1.))
            E=DSQRT(1.D+00*(K+2.)*(K-2.))/(2.*K+1.)
            POL(K+1)=D*(XRMU*POL(K)-E*POL(K-1))
            GAMMA12(K)=GAMMA12(K)+XXX*POL(K)
  12     CONTINUE
    	
	 
         X=XHR(J)*P22(J)*(P11(J)/TTT(J))	!Parentheses needed to ensure strict non-regression
         XX=XHR(J)*P33(J)*P11(J)/TTT(J)
         DO 14 K=0,OS_NB
            PL(K+1)=((2.*K+1.)*XRMU*PL(K)-K*PL(K-1))/(K+1.)
	     BETA22(K)=BETA22(K)+X*PL(K)
            DELTA33(K)=DELTA33(K)+XX*PL(K)
  14     CONTINUE
  
  11  CONTINUE
   
      DO 15 K=0,OS_NB
         BETA22(K)  = BETA22(K) *(2.*K+1.)*.5
         DELTA33(K) = DELTA33(K)*(2.*K+1.)*.5
         GAMMA12(K) = GAMMA12(K)*(2.*K+1.)*.5
   15 CONTINUE


C* Computation of the parameters ALP(k) and ZETA(k) from BETA22(k) and DELTA33(k).
C----------------------------------------------------------------------------------
      DO 16 I=2,OS_NB
         CO1=4*(2*I+1.)/I/(I-1.)/(I+1.)/(I+2.)
         CO2=I*(I-1.)/((I+1.)*(I+2.))
         CO3=CO2*DELTA33(I)
         CO2=CO2*BETA22(I)
         NN=I*.5
         MM=(I-1)*.5
         SOM1=0.
         SOM2=0.
         SOM3=0.
         SOM4=0.
         DO 17 J=1,NN
            X2=(I-1.)*(I-1.)-3.*(2*J-1.)*(I-J)
            SOM1=SOM1+X2*BETA22(I-2*J)
            SOM2=SOM2+X2*DELTA33(I-2*J)
  17     CONTINUE
  
         DO 18 J=0,MM
            X2=(I-1.)*(I-1.)-3.*J*(2*I-2*J-1.)
            SOM3=SOM3+X2*BETA22(I-2*J-1)
            SOM4=SOM4+X2*DELTA33(I-2*J-1)
  18     CONTINUE
  
         ZETA(I)=CO3-CO1*(SOM2-SOM3)
         ALP(I)=CO2-CO1*(SOM1-SOM4)
  16  CONTINUE

C* Normalisation of the parameters ALP(K), BETA(K), GAMMA(K), DELTA(K) and ZETA(K)
C* of the phase matrix.
C----------------------------------------------------------------------------
      Z1=BETA11(0)
      DO 20 K=0,OS_NB
         ALP(K)=ALP(K)/Z1
         BETA11(K)=BETA11(K)/Z1
	  BETA22(K)=BETA22(K)/Z1
         GAMMA12(K)=GAMMA12(K)/Z1
         DELTA33(K)=DELTA33(K)/Z1
         ZETA(K)=ZETA(K)/Z1
  20  CONTINUE

C* Nominal end
C-------------------
       GOTO 9999

C* Error case and return of status -1 to the calling program
C----------------------------------------------------------
  921  WRITE(6,*) 'Error while writing GRANU logfile'
       IER=-1
       GOTO 9999

       
 9999  CONTINUE
 
       RETURN
       END	!END OF PROCEDURE SOS_DECOMPO_LEGENDRE
       
       
       
       

C*==============================================================================
C* PROCEDURE: SOS_GRANU
C* ==========
C      This procedure computes the phase function P(mu) and the polarised
C      phase functions Q(mu) and U(mu) from the MIE computations made
C      by the program "SOS_MIE", for a given wavelength and for
C      a given particle size distribution (given by igranu).
C
C      If igranu = 1 => Log-Normal distribution:
C                N(r) = exp(-0.5 *(ln(r/rm)/sig)**2) / r
C		 varGranu1 = rm : modal radius in microns,
C		 varGranu2 = sig : variance.
C
C      If igranu = 2 => Junge law: N(r) = Nr0 = R0**(-xnu) for r < R0
C				           N(r) = r**(-xnu) for r > R0
C				           N(r) = 0 for r > Rmax
C		 varGranu1 = R0 : minimum radius in microns,
C		 varGranu2 = xnu : power applied to the radius (slope of the Junge law).
C		 varGranu3 = Rmax : maximum radius in microns,

C
C       The phase functions are stored in arrays defined up to
C       CTE_MIE_NBMU_MAX but are used up to MIE_NBMU (actual angles).
C
C
C Description of parameters
C ----------------------------
C       FICMIE(CHARACTER*CTE_LENFIC2)  (I) : Name of the MIE file
C                                 (directory + file name + extension).
C       IGRANU  (I4)         (I) : Number of the size distribution model to use.
C       VARGRANU1 (double)   (I) : Size distribution parameter 1.
C	 VARGRANU2 (double)   (I) : Size distribution parameter 2
C	 VARGRANU3 (double)   (I) : Size distribution parameter 3
C	 WA  (double)         (I) : Wavelength in microns.
C       MIE_NBMU (I4) 	 (I) : Number of actually useful (positive) angles
C	 XMU(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX) (double) : (I) Cosines of the Gauss angles.
C	 TRACE (logical)      (I) : = true if writing in the trace file
C                                 (logical unit number of the trace file =88).
C	 KMAT1 (double)       (O) : Extinction cross section (micron^2).
C	 KMAT2 (double)       (O) : Scattering cross section (micron^2).
C	 SOMME_NR (double)    (O) : Number of particles integrated over the size distribution.
C       P11 (double)         (O) : Phase function computed from the
C                                  functions of the MIE theory.
C       P12 (double)         (O) : Element P12 of the phase matrix.
C	 P33 (double)         (O) : Element P33 of the phase matrix.
C       IER  (I4)               (O) : Error code =0 if no error, =-1 otherwise.
C
C
C Other result provided by the procedure
C ----------------------------------------
C      Information in the trace file (optional)
C   (logical unit number of the trace file =88,
C    file created/opened by the calling program)
C
C
C
C Description of the files used
C ---------------------------------
C   -->Content of the MIE file
C       Unformatted binary file, created by SOS_MIE
C       1 header record with
C                Refractive index (real part)
C                Refractive index (imaginary part)
C                Limit size parameter
C                Number of positive Gauss angles
C
C      1 record for each value of the size parameter with:
C        - Size parameter (real),
C        - Extinction efficiency factor (real),
C        - Scattering efficiency factor (real),
C        - Asymmetry factor (double precision),
C        - Array of the phase function Imie(mu) (real),
C        - Array of the polarised phase function Qmie(mu) (real),
C        - Array of the polarised phase function Umie(mu) (real).
C
C    --> The arrays have a useful size limited to MIE_NBMU:
C        they are read with this limitation.
C
C
C Constants used:
C----------------------
C     CTE_LENFIC2 : Length of file names with directory tree.
C    INCTE_PI : Value of PI
C    CTE_MIE_NBMU_MAX : Maximum number of positive angles.
C    CTE_OS_NB_MAX :  Maximum value of the limit order of the expansion in
C                     Legendre polynomials.
C
C Common used:
C --------------
C None
C
C
C Error cases:
C ------------
C     -Inconsistency of the number of Gauss angles
C     -Opening of a file
C       Error if the MIE file or the Gauss angle file does not exist
C     -Reading/writing in a file
C
C     Message displayed on screen, the routine stops its computations and
C     returns status -1 to the calling program
C        
C   
C
C
C* MOD:VERSION 2.0: 24/06/2008 
C	   - The coefficients KMAT1 and KMAT2 are normalised to one particle
C            ==> The routine outputs extinction and scattering cross sections
C                (in mic^2) instead of extinction and scattering coefficients (/mic).
C          - Introduction of the computation of the effective radius and effective variance.
C
C MOD:VERSION 3.0: 10/01/2010
C          - One dimension of the arrays changed from CTE_MIE_NBMU to CTE_MIE_NBMU_MAX
C          - Loop on the phase functions up to MIE_NBMU
C          - Change of the reading line of the arrays of phase functions
C            per value of the size parameter Alpha
C            ==> limited to the useful bounds -MIE_NBMU:MIE_NBMU of the actual record
C                of the MIE files
C
C MOD:VERSION 4.0: 03/02/2016
C          - The formulation N(r) = exp(-0.5 *(log10(r/rm)/sig)**2) / r
C            for sig = log10(variance) is abandoned in favour
C            of a classical formulation of the LND model:
C			    N(r) = exp(-0.5 *(log(r/rm)/sig)**2) / r	
C	     for sig = variance.
C
C MOD:VERSION 4.1: 08/10/2020
C          - Addition of the argument VARGRANU3 to define the radius Rmax of a
C            Junge size distribution
C          
C==============================================================================


      SUBROUTINE SOS_GRANU(FICMIE,IGRANU,VARGRANU1,VARGRANU2,
     &                     VARGRANU3,WA,MIE_NBMU,XMU,TRACE,
     &                     KMAT1,KMAT2,SOMME_NR,P11,P12,P33,IER)

      IMPLICIT NONE


C* Definition of variables
C----------------------------------------------------------------------------
      
      CHARACTER*CTE_LENFIC2 FICMIE	! Path and name of the MIE file
      
      DOUBLE PRECISION VARGRANU1   	! Size distribution parameter 1.
      DOUBLE PRECISION VARGRANU2   	! Size distribution parameter 2.
      DOUBLE PRECISION VARGRANU3   	! Size distribution parameter 3.
      DOUBLE PRECISION RM		!Modal radius of an LND size distribution.
      DOUBLE PRECISION SIG		!Variance of an LND size distribution
      DOUBLE PRECISION B		!Intermediate variable.
      DOUBLE PRECISION R0		!Min threshold radius of a Junge law size distribution
      DOUBLE PRECISION RMAX		!Max threshold radius of a Junge law size distribution
      DOUBLE PRECISION XNU		!Power of a Junge law size distribution.
      DOUBLE PRECISION NR0		!Number of particles for a radius < r0
                                	!of a Junge law.
					
      DOUBLE PRECISION REFF		!Effective radius of the size distribution (microns).
      DOUBLE PRECISION VEFF		!Effective variance of the size distribution.
      
      DOUBLE PRECISION WA		!Wavelength in microns.

      DOUBLE PRECISION XMU(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)
                                        !Cosines of the Gauss angles.
      DOUBLE PRECISION TETA     	!Scattering angle.
      
      DOUBLE PRECISION KMAT1		! Extinction cross section (micron^2).
      DOUBLE PRECISION KMAT2 		! Scattering cross section (micron^2).
      
      DOUBLE PRECISION P11(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Phase function computed from the functions of
		       !the MIE theory.
      DOUBLE PRECISION P12(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Element P12 of the phase matrix.
      DOUBLE PRECISION P33(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)	
                       !Element P33 of the phase matrix.
		      
      DOUBLE PRECISION XG	!Asymmetry factor.


      DOUBLE PRECISION ALPHAF	!Maximum value of the size parameter contained
      				!in the file of the MIE computations.

      DOUBLE PRECISION R	!Particle radius for the size parameter
      				!ALPHA and the processing wavelength WA.
      DOUBLE PRECISION PR	!Step between two successive values of the particle
      				!radii (in microns).
      DOUBLE PRECISION NR	!NR*PR = Number of particles with a radius between
      				!R and R + PR (per micron**3).

      DOUBLE PRECISION SOMME_NR	!Integration over the radii of NR * PR
      				!==> nb of particles per unit volume.
      DOUBLE PRECISION SOMME_R2	!Integration over the radii of R^2 * NR * PR.
      DOUBLE PRECISION SOMME_R3	!Integration over the radii of R^3 * NR * PR.
      DOUBLE PRECISION SOMME_R4	!Integration over the radii of R^4 * NR * PR.
      DOUBLE PRECISION R2NRPR	!Value of R*R * NR * PR.
      
      DOUBLE PRECISION RN	!real part of the refractive index.
      DOUBLE PRECISION RI	!imaginary part of the refractive index.
	
      DOUBLE PRECISION X1		

      
      REAL*4 ALPHA		!Value of the size parameter.
      REAL*4 PAS		!Step between two values of ALPHA.
      REAL*4 QEXT		!Extinction efficiency factor.
      REAL*4 QSCA		!Scattering efficiency factor.
      REAL*4 IMIE(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)         
                                !Phase function: Imie(mu) for ALPHA.
      REAL*4 QMIE(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)         
                                !Polarised phase function: Qmie(mu).
      REAL*4 UMIE(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)         
                                !Polarised phase function: Umie(mu).

      INTEGER*4 MIE_NBMU  	! Number of actually useful (positive) angles
      INTEGER*4 IGRANU  	! Number of the size distribution model.
      INTEGER*4 J	   	! Index.
      INTEGER*4 NBMUMIE	   	! Number of positive GAUSS angles used for
      			   	! the MIE computations (record of the MIE file).
      INTEGER*4 IER        	! Error code =0 if no error, =-1 otherwise

      LOGICAL TRACE        	! = true if writing in the trace file.
      

C* Initialisation of the error code
C------------------------------------------------------
      IER=0

C* Initialisation of the phase function arrays
C------------------------------------------------------
      DO 102 J=-MIE_NBMU,MIE_NBMU
         P11(J)=0.
         P12(J)=0.
         P33(J)=0.
  102 CONTINUE

C* Initialisation of the extinction and scattering coefficients
C-------------------------------------------------------------
      KMAT1=0.
      KMAT2=0.

C* Initialisation of the integration sums over the size distribution
C-------------------------------------------------------------
      SOMME_NR=0.
      SOMME_R2=0.
      SOMME_R3=0.
      SOMME_R4=0.
                  
C* Opening of the files
C------------------------------------------------------
C     MIE file
      OPEN (1,FILE=FICMIE,FORM='UNFORMATTED',
     &      STATUS='OLD',ERR=991)
 
C* Reading of the header of the MIE file
C-------------------------------------------      
      READ(1,ERR=992) RN,RI,ALPHAF,NBMUMIE
        
C* Parameter check
C--------------------
      IF (NBMUMIE.NE.MIE_NBMU) GOTO 996

C* Reading of the functions computed by the MIE theory
C* for each value of the size parameter ALPHA
C*
C*      Reading of the arrays of phase functions
C*      that were written with bounds
C*      limited to the useful dimension -MIE_NBMU:MIE_NBMU
C*       
C  (loop on label 123,
C   stop condition: end of the MIE file)
C------------------------------------------------------
      PAS = 0.0001

 123     READ(1,err=9992) ALPHA,QEXT,QSCA,XG,
     &                    (IMIE(J),J=-MIE_NBMU,MIE_NBMU),
     &                    (QMIE(J),J=-MIE_NBMU,MIE_NBMU),
     &                    (UMIE(J),J=-MIE_NBMU,MIE_NBMU)
		
C* Computation of the particle radius associated with ALPHA for the wavelength WA
C--------------------------------------------------------------------------
         R=ALPHA*WA/2.D+00/INCTE_PI
	 
C* Test on the limit of the MIE file
C--------------------------------------
         IF (ALPHA.GE.(ALPHAF-PAS)) GOTO 40
		
C* Estimation of the step for the next value of ALPHA
C------------------------------------------------------
         IF (ALPHA.GT.0.10)  PAS=0.001
         IF (ALPHA.GT.1.00 ) PAS=0.01
         IF (ALPHA.GT.10.)   PAS=0.05
         IF (ALPHA.GT.30.)   PAS= 0.10
         IF (ALPHA.GT.100.)  PAS=1.00


C* Expression of the different possible size distributions:
C------------------------------------------------------

C*  Log-normal distribution:
C***************************
C*  RM : Modal radius (microns).
C*  SIG : Variance
C------------------------------
		
         IF (IGRANU.eq.1) THEN
            RM = VARGRANU1
            SIG= VARGRANU2

            B = dlog(R/RM)/SIG
            NR = dexp(-B*B/2.)/(R*SIG*sqrt(2*INCTE_PI))
         ENDIF

C*  Junge law: NR = plateau Rmin**(-XNU) for r <= Rmin
C****************
C*  R0 : Minimum radius (microns).
C---------------------------
         IF (IGRANU.eq.2) THEN
            R0  = VARGRANU1
            XNU = VARGRANU2
            RMAX = VARGRANU3

            NR0= (R0)**(-XNU)

	     IF (R.GT.RMAX) GOTO 40 !Exit of the loop on the particle size
            IF (R.LE.R0) THEN
	       NR = NR0
	     ELSE
	       NR = R**(-XNU)   
	     ENDIF            
         ENDIF


C* Computation of the extinction (kMAT1) and scattering (KMAT2) coefficients.
C* Computation of the phase functions in intensity (P11) and polarisation (P12 and P33).
C-------------------------------------------------------------------------------
         PR=WA*PAS/2./INCTE_PI
         X1=NR*PR*INCTE_PI*R**2
         KMAT1=KMAT1+X1*QEXT
	  KMAT2=KMAT2+QSCA*X1
	
         X1=QSCA*X1
         DO 1 J=-MIE_NBMU,MIE_NBMU
            P11(J)=P11(J)+IMIE(J)*X1
            P12(J)=P12(J)+QMIE(J)*X1
            P33(J)=P33(J)+UMIE(J)*X1
   1     CONTINUE
   
C* Computation of the total number of particles
C-------------------------------------------------------------------------------
      SOMME_NR = SOMME_NR + NR*PR
 
C* Computation of the terms for the effective radius and mean volume
C-------------------------------------------------------------------------------
      R2NRPR = R*R*NR*PR
      SOMME_R2 = SOMME_R2 + R2NRPR
      SOMME_R3 = SOMME_R3 + R2NRPR*R
      SOMME_R4 = SOMME_R4 + R2NRPR*R*R
        
C* Next value of ALPHA
C-----------------------------------------
      GOTO 123 	
   40 CONTINUE

C* Computation of the effective radius and effective variance
C-------------------------------------------------------------------------------      
      IF (TRACE) THEN
         REFF = SOMME_R3 / SOMME_R2
	  VEFF =  1.D+00 +  
     &	         (SOMME_R4 - 2.*REFF*SOMME_R3) / (REFF*REFF * SOMME_R2)
      END IF  

C* Computation of the volume of a mean particle of the size distribution
C-------------------------------------------------------------------------------    
c      VMOY = (4.D+00*CTE_PI/3.D+00)*SOMME_R3 / SOMME_NR
      
C* Normalisation of the phase functions.
C* Copy of the phase function in intensity.
C------------------------------------------------------
      DO 55 J=-MIE_NBMU,MIE_NBMU
         P11(J)=P11(J)/KMAT2
         P12(J)=P12(J)/KMAT2
         P33(J)=P33(J)/KMAT2
  55  CONTINUE
  

C* Normalisation to one particle of the scattering
C* and extinction coefficients.
C* ==> Conversion to scattering and extinction
C*     cross sections
C------------------------------------------------------  
       KMAT1= KMAT1 / SOMME_NR
       KMAT2= KMAT2 / SOMME_NR

       
C* Writing of the trace file
C-------------------------------------------
      IF (TRACE) THEN
        WRITE(88,*,err=921) ' '
        WRITE(88,*,err=921) '* SUBROUTINE GRANU'
	 WRITE(88,*,err=921) '=================='
        WRITE(88,*,err=921) ' '
        WRITE(88,*,err=921) 'SIMULATION PARAMETERS'
	 WRITE(88,*,err=921) '---------------------'
	 WRITE(88,*,err=921) '  MIE file : ',FICMIE
	 WRITE(88,*,err=921) '  Wavelength in microns : ',WA
	 IF (IGRANU.EQ.1) THEN
	      WRITE(88,*,err=921)
     &	           '  Model 1 : Log-Normale Distribution'
	      WRITE(88,*,err=921) '     rm in microns : ',VARGRANU1
	      WRITE(88,*,err=921) '     sig : ',VARGRANU2
	      WRITE(88,*,err=921) '     sig/ln10 : ',
     &	                                VARGRANU2/DLOG(10.D+00)
	 ELSE
	      WRITE(88,*,err=921) "  Model 2 : Junge's Law "
	      WRITE(88,*,err=921) '     R0 in microns : ',VARGRANU1
	      WRITE(88,*,err=921) '     Slope  : ',VARGRANU2	      	   
	 ENDIF
        WRITE(88,*,err=921) ' '
        WRITE(88,*,err=921) 'RESULTS'
	 WRITE(88,*,err=921) '-------'
	 WRITE(88,*,err=921)
     &	    ' Extinction coefficient - no truncation applied'//
     &	    ' (undefined unit: inverse dimension to a length) :',
     &	    KMAT1*SOMME_NR
	 WRITE(88,*,err=921) 
     &	    ' Scattering coefficient - no truncation applied'//
     &	    ' (undefined unit: inverse dimension to a length) :',
     &     KMAT2*SOMME_NR
        WRITE(88,*,err=921) 
     &	    ' Number of particles (undefined unit: dimension'//
     &	    ' proportional to a length^3) :', SOMME_NR
        WRITE(88,*,err=921)
     &	    ' Extinction cross-section (mic^2) - no trunca. applied : '
     &     ,KMAT1
	 WRITE(88,*,err=921) 
     &	    ' Scattering cross-section (mic^2) - no trunca. applied : ' 
     &     ,KMAT2
	 WRITE(88,*,err=921) 
     &	    ' Effective radius of the size distribution (microns) : ',
     &	    REFF
        WRITE(88,*,err=921) 
     &	    ' Effective variance of the size distribution : ',VEFF
        WRITE(88,*,err=921) ' '
	 WRITE(88,*,err=921) '  1st column : Scattering angle'
	 WRITE(88,*,err=921) '  2nd column : Original phase function'
        WRITE(88,*,err=921) '                 P11(mu) '
	 WRITE(88,*,err=921) 
     &       '  3rd column : Original polarized phase function'
        WRITE(88,*,err=921) '                 P12(mu) '
	 WRITE(88,*,err=921) 
     &       '  4th column : Original polarized phase function '
        WRITE(88,*,err=921) '                 P33(mu) '
	 WRITE(88,*,err=921) ' '
        DO 44 J=-MIE_NBMU,MIE_NBMU
	    TETA=DACOS(XMU(J))*180./INCTE_PI
	    WRITE(88,10,err=921) TETA,P11(J),P12(J),P33(J)

   44   CONTINUE
   	 WRITE(88,*,err=921) ' '

      END IF  

 


C* File closing
C-------------------
      CLOSE(1)
      
C* Nominal end
C-------------------
      GOTO 9999

      
C* Error case and return of status -1 to the calling program
C----------------------------------------------------------
  991 WRITE(6,*) 'Erreur a l''ouverture du fichier MIE'
      IER=-1
      GOTO 9999
  992 WRITE(6,*) 'Erreur de lecture de la premiere ligne du fichier MIE'
      WRITE(6,*) '  --> Fichier de MIE probablement incomplet'
      WRITE(6,*) '      ==> Controlez l''espace disque'
      IER=-1
      GOTO 9999
 9992 WRITE(6,*) 'Erreur de lecture sur le fichier MIE'
      WRITE(6,*) '  --> Probable tentative de lecture apres fin fichier'
      WRITE(6,*) '  --> Fichier de MIE probablement incomplet'
      WRITE(6,*) '      ==> Controlez l''espace disque'
      IER=-1
      GOTO 9999    
  921 WRITE(6,*) 'Erreur d''écriture sur le fichier trace GRANU'
      IER=-1
      GOTO 9999      
  996 WRITE(6,*) 'Incoherence sur le nombre d''angles de Gauss'
      write(6,*) 'NBMUMIE=',NBMUMIE,' MIE_NBMU=',MIE_NBMU
      IER=-1  
      GOTO 9999

C* Format
C---------   
   10 FORMAT(1X,F6.2,6(2X,E12.4))
   
         
 9999 CONTINUE
 
  
      RETURN
      END    !END OF PROCEDURE SOS_GRANU
      
      
      
      
      


C*==============================================================================
C* PROCEDURE: SOS_INTERPO_SPLINT
C* ==========
C      This procedure interpolates a function defined by (XIN, YIN) on NBIN nodes
C      at the abscissas XOUT for NBOUT nodes (XOUT, YOUT).
C
C      The approach relies on the functions SPLINE and SPLINT
C         - The nodes (XIN, YIN) are first sorted by increasing abscissa XIN.
C         - The function SPLINE computes the second derivatives of YIN_TRIE at XIN_TRIE
C           from the values YIN_TRIE of the function at the abscissas XIN_TRIE and the
C           values of the first derivative of the function at the first and last nodes
C           (min and max of XIN_TRIE).
C	  - The function SPLINT computes the values YOUT at the abscissas XOUT from
C           the values of the function and of its second derivative previously computed
C           at the abscissas XIN_TRIE.
C
C Description of parameters
C ----------------------------
C       NBIN (I4)		(I) : Number of elements of the arrays XIN and YIN
C       XIN(1:NBIN)   (double)	(I) : Array of the abscissas of the function for which the ordinate is known.
C       YIN(1:NBIN)   (double)	(I) : Array of the known values of the function (at the abscissas XIN).
C       NBOUT (I4)		(I) : Number of elements of the arrays XOUT and YOUT
C       XOUT(1:NBOUT) (double)	(I) : Array of the requested abscissas of the function.
C       YOUT(1:NBOUT) (double)	(O) : Array of the estimated values of the function at the abscissas XOUT.
C       IER  (I4)               (O) : Error code =0 if no error, =-1 otherwise.
C
C Constants used:
C----------------------
C   None
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C     - Error return of the routine SOS_SPLINT.
C
C     Message displayed on screen, the routine stops its computations and
C     returns status -1 to the calling program
C        
C  Subroutines used:
C --------------------------
C  This program calls the routines: SOS_SPLINE and SOS_SPLINT
C
C==============================================================================

      SUBROUTINE SOS_INTERPO_SPLINT(NBIN,XIN,YIN,NBOUT,XOUT,YOUT,IER)

      IMPLICIT NONE


C* Definition of variables
C----------------------------------------------------------------------------
      INTEGER*4 NBIN  	! Number of nodes (XIN,YIN) of the function Y
      INTEGER*4 NBOUT  	! Number of nodes (XOUT,YOUT) of the function Y
      INTEGER*4 J	! Index.
      INTEGER*4 K	! Index.
      INTEGER*4 IER     ! Error code =0 if no error, =-1 otherwise
      
      DOUBLE PRECISION XIN(NBIN)   !Values of the abscissas of the nodes (XIN,YIN)
      DOUBLE PRECISION YIN(NBIN)   !Values of the ordinates of the nodes (XIN,YIN)
      DOUBLE PRECISION XOUT(NBOUT)  !Values of the abscissas of the nodes (XOUT,YOUT)
      DOUBLE PRECISION YOUT(NBOUT)  !Requested values of the ordinates of the nodes (XOUT,YOUT)

      DOUBLE PRECISION XIN_TRIE(NBIN)   !Increasing values of the abscissas of the nodes (XIN,YIN)
      DOUBLE PRECISION YIN_TRIE(NBIN)   !Values of the ordinates of the nodes (XIN,YIN) for XIN sorted in increasing order
      DOUBLE PRECISION VAL		!Temporary value
      
      DOUBLE PRECISION DY1	!First derivative of YIN_TRIE at the abscissa XIN_TRIE(1)
      DOUBLE PRECISION DYN	!First derivative of YIN_TRIE at the abscissa XIN_TRIE(NBIN)
      DOUBLE PRECISION DSEC_YIN_TRIE(NBIN)  !Values of the second derivatives of the function YIN_TRIE
      					    !at the nodes XIN_TRIE
    

C* Initialisation of the error code
C------------------------------------------------------
      IER=0

      
C* Sorting of the tables (XIN, YIN) by increasing XIN
C------------------------------------------------------
      DO J=1,NBIN
         XIN_TRIE(J)=XIN(J)
	 YIN_TRIE(J)=YIN(J)
      ENDDO
      
      DO J=1,NBIN	     
	 DO K=J+1,NBIN
	    IF (XIN_TRIE(J).GT.XIN_TRIE(K)) THEN		    		    
	       !Swap of the values XIN(J) and XIN(K)
	        VAL=XIN_TRIE(J)
	        XIN_TRIE(J)=XIN_TRIE(K)
	        XIN_TRIE(K)=VAL
	       !Swap of the values YIN(J) and YIN(K)
		 VAL=YIN_TRIE(J)
		 YIN_TRIE(J)=YIN_TRIE(K)
		 YIN_TRIE(K)=VAL   	    	    
	    ENDIF
	 ENDDO 
      ENDDO 
   
      
C* Computation of the second derivatives of YIN_TRIE at XIN_TRIE
C------------------------------------------------------
      DO J=1,NBIN
         DSEC_YIN_TRIE(J)=0.D+00
      ENDDO

C     First derivatives of YIN_TRIE at XIN_TRIE(1) and XIN_TRIE(NBIN)
      DY1 = (YIN_TRIE(2) - YIN_TRIE(1)) / (XIN_TRIE(2) - XIN_TRIE(1))
      DYN = (YIN_TRIE(NBIN) - YIN_TRIE(NBIN-1)) / 
     &      (XIN_TRIE(NBIN) - XIN_TRIE(NBIN-1))
      
      CALL SOS_SPLINE(NBIN,XIN_TRIE,YIN_TRIE,DY1,DYN,DSEC_YIN_TRIE)

C* Application of the SPLINT interpolation at the abscissas XOUT
C------------------------------------------------------------       
       DO J=1,NBOUT
	 
	 CALL SOS_SPLINT(NBIN,XIN_TRIE,YIN_TRIE,DSEC_YIN_TRIE,
     &	                 XOUT(J),YOUT(J),IER)
     	 
         IF (IER.NE.0) GOTO 901
      ENDDO
            
      
C* Nominal end
C-------------------
      GOTO 9999

      
C* Error case and return of status -1 to the calling program
C----------------------------------------------------------
  901 WRITE(6,*) '  ERROR for SOS_INTERPO_SPLINT'
      GOTO 9999  
     
 9999 CONTINUE
  
      RETURN
      END    !END OF PROCEDURE SOS_INTERPO_SPLINT
      



C*==============================================================================
C* PROCEDURE: SOS_SPLINE
C* ==========
C      This procedure computes the values of the second derivatives of a function defined
C      by (X, Y) on N nodes (assumed sorted by increasing X)
C      and by its first derivatives at the extreme nodes (X(1),DY1) and (X(N),DYN).
C
C      The routine returns the values of the second derivatives at each point X.
C
C Description of parameters
C ----------------------------
C       N  (I4)		(I) : Number of elements of the arrays X, Y and DSEC_Y
C       X(1:N) (double)	(I) : Array of the abscissas of the function for which the ordinate is known.
C       Y(1:N) (double)	(I) : Array of the known values of the function (at the abscissas X).
C       DY1    (double)	(I) : First derivative of the function at X(1).
C       DYN    (double)	(I) : First derivative of the function at X(N).
C       DSEC_Y(1:N) (double) (O) : Array of the values of the second derivatives of the function at the abscissas X.
C
C Constants used:
C----------------------
C   None
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

      SUBROUTINE SOS_SPLINE(N,X,Y,DY1,DYN,DSEC_Y)

      IMPLICIT NONE

      INTEGER*4 N		 !Nb of elements of the input (X,Y) and output (DSEC_Y) arrays
      DOUBLE PRECISION X(N)	 !Array of the abscissas of the function
      DOUBLE PRECISION Y(N)	 !Array of the values of the function at the abscissas X
      DOUBLE PRECISION DY1   !First derivative of the function at X(1)
      DOUBLE PRECISION DYN	 !First derivative of the function at X(N)
      DOUBLE PRECISION DSEC_Y(N) !Array of the values of the second derivatives of the function at the abscissas X
      
      DOUBLE PRECISION U(N)

      DOUBLE PRECISION SIG
      DOUBLE PRECISION P
      DOUBLE PRECISION QN
      DOUBLE PRECISION UN
      
      INTEGER*4 K ! Index.
      
C* Definition of variables
C----------------------------------------------------------------------------


      IF (DY1.GT..99E30) THEN
        DSEC_Y(1)=0.
        U(1)=0.
      ELSE
        DSEC_Y(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-DY1)
      ENDIF
      DO 11 K=2,N-1
        SIG=(X(K)-X(K-1))/(X(K+1)-X(K-1))
        P=SIG*DSEC_Y(K-1)+2.
        DSEC_Y(K)=(SIG-1.)/P
        U(K)=(6.*((Y(K+1)-Y(K))/(X(K+1)-X(K))-(Y(K)-Y(K-1))
     *      /(X(K)-X(K-1)))/(X(K+1)-X(K-1))-SIG*U(K-1))/P
11    CONTINUE
      IF (DYN.GT..99E30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(DYN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      DSEC_Y(N)=(UN-QN*U(N-1))/(QN*DSEC_Y(N-1)+1.)
      DO 12 K=N-1,1,-1
        DSEC_Y(K)=DSEC_Y(K)*DSEC_Y(K+1)+U(K)
12    CONTINUE
      RETURN
      END   !END OF PROCEDURE SOS_SPLINE





C*==============================================================================
C* PROCEDURE: SOS_SPLINT
C* ==========
C
C      This procedure interpolates a function at a point XVAL from
C      three vectors of the same dimension X(1:N), Y(1:N) and DSEC_Y(1:N)
C      giving respectively the abscissas, the ordinates and the
C      second derivatives of the function.
C
C Description of parameters
C ----------------------------
C       N  (I4)		(I) : Number of elements of the arrays X, Y and DSEC_Y
C       X(1:N) (double)	(I) : Array of the abscissas of the function for which the ordinate is known.
C       Y(1:N) (double)	(I) : Array of the known values of the function (at the abscissas X).
C       DSEC_Y(1:N) (double) (I) : Array of the values of the second derivatives of the function at the abscissas X.
C       XVAL (double)	(I) : Value of the abscissa at which the function must be interpolated.
C       YVAL (double)	(O) : Value of the function interpolated at XVAL.
C       IER  (I4)               (O) : Error code =0 if no error, =-1 otherwise.
C
C Constants used:
C----------------------
C   None
C
C Common used:
C --------------
C None
C
C Error cases:
C ------------
C  Identical abscissas of the function for the indices bounding the
C  requested abscissa.
C
C==============================================================================

      SUBROUTINE SOS_SPLINT(N,X,Y,DSEC_Y,XVAL,YVAL,IER)

      IMPLICIT NONE

      INTEGER*4 N	!Nb of elements of the input arrays (X,Y and DSEC_Y)
      INTEGER*4 KLO	!Lower value of index
      INTEGER*4 KHI	!Higher value of index   
      INTEGER*4 K	! Index
      INTEGER*4 IER     ! Error code =0 if no error, =-1 otherwise
      
      DOUBLE PRECISION X(N)	 !Array of the abscissas of the function
      DOUBLE PRECISION Y(N)	 !Array of the values of the function at the abscissas X
      DOUBLE PRECISION DSEC_Y(N) !Array of the values of the second derivatives of the function at the abscissas X

      DOUBLE PRECISION XVAL	!Value of the abscissa at which the function must be interpolated.
      DOUBLE PRECISION YVAL	!Value of the function interpolated at XVAL.
      
      DOUBLE PRECISION A,B
      DOUBLE PRECISION H


C* Initialisation of the error code
C------------------------------------------------------
      IER=0
      
C* Definition of variables
C----------------------------------------------------------------------------
      
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(X(K).GT.XVAL)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=X(KHI)-X(KLO)
      IF (H.NE.0.) THEN
         A=(X(KHI)-XVAL)/H
         B=(XVAL-X(KLO))/H
         YVAL=A*Y(KLO)+B*Y(KHI)+
     *      ((A**3-A)*DSEC_Y(KLO)+(B**3-B)*DSEC_Y(KHI))*(H**2)/6.
      ELSE
         IER=-1
      ENDIF
      
      IF (IER.EQ.-1) GOTO 901
      
C* Nominal end
C-------------------
      GOTO 9999
      
C* Error case and return of status -1 to the calling program
C----------------------------------------------------------
  901 WRITE(6,*) '  ERROR for SPLINT interpolation'
      WRITE(6,*) '  --> Bad X table input'
      GOTO 9999  
     
 9999 CONTINUE      
      RETURN
      END    !END OF PROCEDURE SOS_SPLINT
