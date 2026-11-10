C******************************************************************************
C* FILE: SOS.h
C* PROJECT: Successive orders of scattering
C* ROLE: Definition of the constants of the programs of the OS code.
C*
C* AUTHOR: B.Lafrance ( CS )
C* DATE: 30/01/2002
C*
C* MOD:VERSION:3.0
C*
C* MOD:VERSION:3.1 22/02/2016
C*    - Introduction of the constant CTE_SOLAR_DISC_SOLID_ANGLE
C*    - Introduction of the constant CTE_REMOVE_PREVIOUS_BIN_FILE
C*
C* MOD:VERSION:3.2 10/11/2020
C*    - Introduction of the constant CTE_SEUIL_NUM
C*    - Adjustment of the values of the dimensions of the CKD data:
C*      CTE_CKD_NWVL_MAX, CTE_CKD_NAI_MAX, CTE_CKD_NT_MAX, CTE_CKD_NP_MAX
C*    - Introduction of the constant CTE_HT_STD_PSURF
C*    - Introduction of the constant CTE_MAX_NB_MODE_MIXTURE
C*
C* MOD:VERSION:3.3 31/08/2022
C*    - Introduction of the constant CTE_THRESHOLD_Q_U_NULL
C*    - Introduction of the constants

C* MOD:VERSION:3.4 15/06/2023
C*    - Adjustment of the sizing constants of the CKD files:
C*      CTE_CKD_NUMAX (from 25000 to 27500 cmm1) and CTE_CKD_NB_NU_PER_FILE (from 225 to 50)
C*      CTE_CKD_NWVL_MAX (from 225 to 50), CTE_CKD_NAI_MAX (from 20 to 5)
C*      CTE_CKD_NT_MAX  (from 8 to 9) and CTE_CKD_NP_MAX  (from 25 to 31)
C*    - Addition of the parameter CTE_CKD_NC_MAX giving the max number of water vapour concentrations

C* MOD:VERSION:3.5 27/07/2023
C*    - Addition of constants for the definition of the profile:
C*      CTE_TOA_FIRST_LAYER_OPT_THICKNESS : Optical thickness at the first level of the profile
C*      CTE_OS_NT_MIN : Minimum number of layers
C*      CTE_DELTA_Z : Altitude step for the search of the altitude of the first level
C*    - The parameter CTE_SEUIL_TAUABS is renamed CTE_THRESHOLD_TAUABS
C*
C*    - Addition of the parameter CTE_GAP_TOLER_SUM_RATES defining the maximum tolerated error on the difference to 1
C*      of the sum of the optical thickness rates of an aerosol mixture
C*
C* MOD:VERSION:3.6 18/09/2023
C*    - Conversion of CTE_TOA_ALT from 210. to 120.D+00 for compilation with gcc/10.2.0
C*
C#######################################################################
C    1. Common constants
C#######################################################################

C Length of the character strings
C-----------------------------------
C Directories
#define CTE_LENDIR 350

C Name of the files without directory tree
#define CTE_LENFIC1 150

C Complete name of the files (with directory tree)
#define CTE_LENFIC2 500

C Maximum size of the keywords for argument passing
#define CTE_LENKEYWORD 30

C Length of command system chains
#define CTE_LENCOM 400


C Validity domain (in microns)
C---------------------------------
#define CTE_WAMIN 0.364
#define CTE_WAMAX 4.0

C ID for values not defined by the user
C--------------------------------------
C --> Integer type :
#define CTE_NOT_DEFINED_VALUE_INT -999
C --> Double type :
#define CTE_NOT_DEFINED_VALUE_DBLE -999.D+00

C Default files
C--------------------
#define CTE_DEFAULT_FICANGLES_RES_LUM  "SOS_UsedAngles.txt"
#define CTE_DEFAULT_FICANGLES_RES_MIE  "Aer_UsedAngles.txt"
#define CTE_DEFAULT_FICGRANU  "Aerosols.txt"
#define CTE_DEFAULT_RESBIN   "SOS_Result.bin"
#define CTE_DEFAULT_RESUP  "SOS_Up.txt"
#define CTE_DEFAULT_RESDOWN "SOS_Down.txt"
      
C#######################################################################
C    2. Constants specific to the routines of the program SOS_AEROSOLS
C#######################################################################

C Dimension of the arrays --> constant specific to the program SOS_MIE
C Must be larger than 2*AlphaMax+20 (AlphaMax: max size parameter).
C----------------------------------------------------------------
#define CTE_MIE_DIM 10000

C Maximum size of the arrays of external phase functions
C --> constant specific to the program SOS_AEROSOLS
C----------------------------------------------------------------
#define CTE_MAXNB_ANG_EXT 200


C Default value for the maximum radius of aerosols particles
C for a Junge model of size distribution
C    --> Specific constant to the file SOS_AEROSOLS.F
C----------------------------------------------------------------
#define CTE_DEFAULT_AER_JUNGE_RMAX 50.


C Minimum value of the size parameter Alpha of the MIE computations.
C (must be between 0.0001 and 10.0
C  NB: necessarily in double precision)
C    --> constant specific to the program SOS_AEROSOLS
C----------------------------------------------------------------
#define CTE_MIE_ALPHAMIN 0.0001D+00

C Limit value of the size parameter for the computation of the MIE
C files of the components of the WMO and Shettle & Fenn models
C    --> constant specific to the program SOS_AEROSOLS
C----------------------------------------------------------------
#define CTE_ALPHAMAX_WMO_DL   4000.
#define CTE_ALPHAMAX_WMO_WS   50.
#define CTE_ALPHAMAX_WMO_OC   800.
#define CTE_ALPHAMAX_WMO_SO   10.
#define CTE_ALPHAMAX_SF_SR    70.
#define CTE_ALPHAMAX_SF_SU    90.

C Minimum value of the size distribution ratio n(r) / Nmax
C for the determination of the limit size parameter
C needed by the Mie computations.
C    --> constant specific to the program SOS_AEROSOLS
C----------------------------------------------------------------
#define CTE_COEF_NRMAX 0.0001


C Name of the file containing the modal radius, the log of the variance,
C the volume concentration and the values of the refractive indices
C (real part and imaginary part) versus the
C wavelength.
C    --> constant specific to the program SOS_AEROSOLS
C---------------------------------------------------------------
#define CTE_AER_DATAWMO		"Data_WMO_cor_2015_12_16"

C Name of the file containing the log of the variance and the values
C of the modal radius versus the relative humidity.
C    --> constant specific to the program SOS_AEROSOLS
C---------------------------------------------------------------
#define CTE_AER_DATASF		"Data_SF_cor_2015_12_16"

C Name of the files containing the values of the refractive indices
C (real part and imaginary part) versus the
C wavelength and the relative humidity.
C    --> constants specific to the program SOS_AEROSOLS
C---------------------------------------------------------------
#define CTE_AER_SR_SF		"IRefrac_SR_cor_2015_12_16"
#define CTE_AER_LR_SF		"IRefrac_LR"
#define CTE_AER_SU_SF		"IRefrac_SU_cor_2015_12_16"
#define CTE_AER_LU_SF		"IRefrac_LU_cor_2015_12_16"
#define CTE_AER_OM_SF		"IRefrac_OM_cor_2015_12_16"

C Angles used to define the angular range for the phaze function 
C linearization while applying a truncation.
C    --> Specific constant to the files SOS_AEROSOLS.F
C---------------------------------------------------------------
#define CTE_AER_MU1_TRONCA 0.8
#define CTE_AER_MU2_TRONCA 0.94

C Threshold for the truncation
C    --> constants specific to the program SOS_AEROSOLS
C---------------------------------------------------------------
#define CTE_PH_SEUIL_TRONCA 0.1


C Maximum number of modes for an aerosol mixture
C    --> constants specific to the program SOS_AEROSOLS
C---------------------------------------------------------------
#define CTE_MAX_NB_MODE_MIXTURE 20

C Tolerance threshold on the difference to 1 of the sum of the
C optical thickness rates of an aerosol mixture
C    --> constants specific to the program SOS_AEROSOLS
C---------------------------------------------------------------
#define CTE_GAP_TOLER_SUM_RATES 0.000001

      
C#######################################################################
C    3. Constants specific to the definition of the atmospheric profile
C#######################################################################

C Standard atmospheric pressure
#define CTE_HT_STD_PSURF 1013.0

C Max altitude (in km) of the atmospheric profile (TOA level)
C       --> constant specific to the programs SOS, SOS_PROFIL and SOS_OS
C-------------------------------------------------------------------
#define CTE_TOA_ALT 120.D+00

C Maximum number of layers of the atmospheric profile.
C       --> constant specific to the programs SOS, SOS_PROFIL and SOS_OS
C-------------------------------------------------------------------
#define CTE_OS_NT 600

C Minimum optical thickness of a layer of the atmospheric profile
C (except the first level).
C       --> constant specific to the programs SOS, SOS_PROFIL and SOS_OS
C-------------------------------------------------------------------
#define CTE_TCOUCHE 0.005

C Optical thickness of the first layer of the atmospheric profile.
C       --> constant specific to the program SOS_PROFIL
C-------------------------------------------------------------------
#define CTE_TOA_FIRST_LAYER_OPT_THICKNESS 0.0002

C Step (in km) for the search of the altitude of the first level
C       --> constant specific to the program SOS_PROFIL
C-------------------------------------------------------------------
#define CTE_DELTA_Z 0.05

C Threshold (in km) for the comparison of the altitudes when defining
C the profile
C       --> constant specific to the program SOS_PROFIL
C-------------------------------------------------------------------
#define CTE_THRESHOLD_DZ 0.001

C Minimum number of layers of the atmospheric profile.
C       --> constant specific to the programs SOS, SOS_PROFIL and SOS_OS
C-------------------------------------------------------------------
#define CTE_OS_NT_MIN 100

C Minimum number of molecular sub-layers for a profile with
C variable positioning of the aerosols between two altitudes.
C       --> constant specific to the program SOS_PROFIL
C-------------------------------------------------------------------
#define CTE_PROFIL_MIN_NBC 3

C Definition of the thickness of the transition layer (km)
C    --> constant specific to the program SOS_PROFIL
C-------------------------------------------------------------------
#define CTE_DZTRANSI	0.010

C Number of absorbers (there are 8 species to take into account,
C including the absorption lines and/or the continuums):
C --> H2O CO2 O3 N2O CO CH4 O2 NO2
C----------------------------------------------------------------------------
#define CTE_NBABS  8

C Number of levels of the absorbing gas profiles
C-------------------------------------------------------
#define CTE_ABS_NBLEV  50

C Number of columns of the absorbing gas profile files
C----------------------------------------------------------------
#define CTE_ABS_NBCOL  13

C Definition of the absorption bands to consider (1 to take into account, 0 otherwise)
C Absorbers: H2O CO2 O3 N2O CO CH4 O2 NO2
C The NO2 absorptions are only defined by continuums.
C The other gases include lines and continuums.
C-------------------------------------------------------------------------------
#define CTE_ABS_H2O  1
#define CTE_ABS_CO2  1
#define CTE_ABS_O3   1
#define CTE_ABS_N2O  1
#define CTE_ABS_CO   1
#define CTE_ABS_CH4  1
#define CTE_ABS_O2   1
#define CTE_CONT_NO2 1

C Maximum dimensions of the CKD coefficient files
C-------------------------------------------------------------------------------
C CTE_CKD_NWVL_MAX : Maximum number of wavelengths in the CKD coeff tables
C CTE_CKD_NAI_MAX  : Maximum number of exponentials of the expansions
C CTE_CKD_NT_MAX   : Maximum number of temperatures
C CTE_CKD_NP_MAX   : Maximum number of pressures
C CTE_CKD_NC_MAX   : Maximum number of H2O concentrations (specific to water vapour)
C-------------------------------------------------------------------------------
#define CTE_CKD_NWVL_MAX  50
#define CTE_CKD_NAI_MAX  5
#define CTE_CKD_NT_MAX  9
#define CTE_CKD_NP_MAX  31
#define CTE_CKD_NC_MAX  12

C Definition of the spectral splitting of the CKD coefficient files
C-------------------------------------------------------------------------------
C CTE_CKD_NUMAX   : Maximum wavenumber covered by the CKD coeff files (in cm-1)
C CTE_CKD_NUMIN   : Minimum wavenumber covered by the CKD coeff files (in cm-1)
C CTE_CKD_NB_NU_PER_FILE : Number of wavelengths per CKD file (must be the same for all the spectral resolutions)
C-----------------------------------------------
#define CTE_CKD_NUMAX  27500
#define CTE_CKD_NUMIN   2500
#define CTE_CKD_NB_NU_PER_FILE  50


C Maximum absorption optical thickness (applies to extreme cases if the transmission is null)
C-------------------------------------------------------------------------------
#define CTE_TAUABS_MAX   999.

C Threshold on the absorption optical thickness to optimise the atmospheric profile
C ------------------------------------------------------------------------------------
#define CTE_THRESHOLD_TAUABS 1.5


C#######################################################################
C    4. Constants specific to the programs computing the surface
C       reflection matrices
C#######################################################################

C Factor for the comparison of GMIN and GMAX
C    --> constant specific to the program SOS_GLITTER
C--------------------------------------------------------------
#define CTE_PH_TEST 10000


C Number of azimuth angles (2**CTE_PH_NQ)
C    --> constant specific to the programs
C        SOS_GLITTER, SOS_ROUJEAN and SOS_NADAL
C--------------------------------------------------------------
#define CTE_PH_NU 1024


C Exponent of the power of 2 of CTE_PH_NU (NU=2**NQ)
C    --> constant specific to the program SOS_GLITTER
C--------------------------------------------------------------
#define CTE_PH_NQ 10


C Threshold to stop the Fourier series decomposition
C of the function (F21 Nadal) / (F21 Fresnel)
C    --> constant specific to the program SOS_NADAL
C--------------------------------------------------------------
#define CTE_SEUIL_SF_NADAL 0.001


C Threshold to stop the Fourier series decomposition
C of the Roujean BRDF
C    --> constant specific to the program SOS_ROUJEAN
C--------------------------------------------------------------
#define CTE_SEUIL_SF_ROUJEAN 0.001


C Limit solar zenith angle for the computation of the Roujean BRDF
C (in degrees)
C       --> constant specific to the programs SOS_ROUJEAN
C           and SOS_NOM_FIC_SURFACE
C--------------------------------------------------------------
#define CTE_TETAS_LIM_ROUJEAN 60


C Limit viewing zenith angle for the computation of the Roujean BRDF
C (in degrees)
C       --> constant specific to the programs SOS_ROUJEAN
C           and SOS_NOM_FIC_SURFACE
C--------------------------------------------------------------
#define CTE_TETAV_LIM_ROUJEAN 60


C Threshold for the tests of numerical rounding errors
C    --> constant specific to the program SOS_ROUJEAN
C--------------------------------------------------------------
#define CTE_SEUIL_NUM 1.D-10




C#######################################################################
C    5. Constants specific to the routines of the program SOS.
C#######################################################################

C Factor of molecular depolarization in the air
C    --> constant specific to the program SOS_OS
C--------------------------------------------------------------
#define CTE_MDF 0.0279

C Minimum order of the Fourier series expansion
C    --> constant specific to the program SOS_OS
C------------------------------------------------------------
#define CTE_OS_IBOR 0

C Default maximum number of scattering orders
C    --> constant specific to the program SOS_OS
C------------------------------------------------------------
#define CTE_DEFAULT_IGMAX 100


C Convergence threshold of the geometric series
C    --> constant specific to the program SOS_OS
C------------------------------------------------------------
#define CTE_PH_SEUIL_CV_SG  0.00001

C Stop thresholds of the multiple scattering
C    --> constants specific to the program SOS_OS
C------------------------------------------------------------
#define CTE_PH_SEUIL_SUMDIF  0.00001
#define CTE_PH_SEUIL_VALDIF  1.0D-50

C Threshold to stop the Fourier series decomposition
C    --> constant specific to the program SOS_OS
C------------------------------------------------------------
#define CTE_PH_SEUIL_SF  0.00001



C Threshold for the computation of the rotation angles
C    --> constant specific to the file SOS_TRPHI (routine SOS_ANGLE)
C--------------------------------------------------------------
#define CTE_SEUIL_Z 0.0001

C Threshold for the computation of the elements M21 and M31 of the
C reflection matrix expressed in a frame tied to the meridian plane
C    --> constant specific to the file SOS_TRPHI (routine SOS_MATRIC)
C--------------------------------------------------------------
#define CTE_SEUIL_X 0.00001

C Threshold value under which Q or U is fixed to be null
C    --> Specific constant to the file SOS_TRPHI.F
C--------------------------------------------------------------
#define CTE_THRESHOLD_Q_U_NULL 1.D-15



C Solid angle of the solar disc (sr) for
C the mean Earth - Sun distance
C    --> constant specific to the file SOS_TRPHI
C--------------------------------------------------------------
#define CTE_SOLAR_DISC_SOLID_ANGLE 6.8D-05

C Option to remove or not the file -SOS.ResBin if it already exists
C from a previous simulation.
C    0 : no removal --> The code stops with a stop message.
C    1 : the code removes the file with a warning message.
#define CTE_REMOVE_PREVIOUS_BIN_FILE 1


C#######################################################################
C    6. Constants specific to the definition of the angles and orders
C       of the expansions
C#######################################################################

C ===================================================
C 6.1 Constants for the sizing of the arrays
C     of angles and limit orders of the expansions
C     in Fourier series and Legendre polynomials.
C ===================================================

C    Maximum number of (positive) angles for the sizing
C    of the phase function arrays.
C    Advice: * Must be at least the default value CTE_DEFAULT_NBMU_MIE
C                + the max number of user angles CTE_NBMAX_USER_ANGLES
C              * To be adjusted if the common use is to have more Gauss
C                angles than the default value.
C              * But do not define an excessively large value, so as not to
C                increase the processing time.
C    --> constant specific to the programs SOS_AEROSOLS, SOS_MIE
C        and SOS_ANGLES
C--------------------------------------------------------------
#define CTE_MIE_NBMU_MAX 100		

C    Maximum number of (positive) angles for the sizing of the
C    arrays of radiance fields and of surface reflection
C    Advice: * Must be at least the default value CTE_DEFAULT_NBMU_LUM
C                + the max number of user angles CTE_NBMAX_USER_ANGLES
C              * To be adjusted if the common use is to have more Gauss
C                angles than the default value.
C              * But do not define an excessively large value, so as not to
C                increase the processing time.
C    --> constant specific to the programs:
C        SOS_SURFACE, SOS_GLITTER, SOS_ROUJEAN, SOS_RONDEAUX_BREON
C        SOS_NADAL, SOS, SOS_OS, SOS_TRPHI and SOS_ANGLES
C-----------------------------------------------------------------------
#define CTE_OS_NBMU_MAX 80    

C    Maximum order of the Fourier series expansion
C    and limit order of the expansion in Legendre polynomials
C    for the sizing of the arrays
C    Advice: Take CTE_OS_NB_MAX >= 2*CTE_MIE_NBMU_MAX
C       --> constant specific to the programs
C           SOS, SOS_OS, SOS_AEROSOLS, SOS_SURFACE, SOS_ROUJEAN and SOS_NADAL
C---------------------------------------------------------------------------
#define CTE_OS_NB_MAX 200

C    Maximum order of the expansion in Legendre polynomials
C    of the elements of the Fresnel matrix
C    for the sizing of the arrays
C    Advice: Take CTE_OS_NS_MAX = 2*CTE_OS_NBMU_MAX
C       --> constant specific to the program SOS_SURFACE
C---------------------------------------------------------------------------
#define CTE_OS_NS_MAX 136

C    Maximum order of the Fourier series expansion
C    of the G function for the sizing of the arrays
C    Advice: Take CTE_OS_NM_MAX = CTE_OS_NB_MAX + CTE_OS_NS_MAX
C       --> constant specific to the programs
C           SOS_SURFACE, SOS_GLITTER and SOS_RONDEAUX_BREON
C---------------------------------------------------------------------------
#define CTE_OS_NM_MAX 336


C ===================================================
C 6.2 Number of angles to use and maximum orders
C     of the expansions when the angles are not defined
C     by the user.
C ===================================================

C    Default number of positive Gauss angles for the aerosol computations
C    --> constant specific to the program SOS_ANGLES
C--------------------------------------------------------------
#define CTE_DEFAULT_NBMU_MIE 40

C    Default number of positive Gauss angles for the computations of
C    radiance fields
C    --> constant specific to the program SOS_ANGLES
C--------------------------------------------------------------
#define CTE_DEFAULT_NBMU_LUM 24

C    Default value of the limit order of the expansion in Legendre
C    polynomials and Fourier series of the radiance field:
C    Usually: CTE_DEFAULT_OS_NB = 2*CTE_DEFAULT_NBMU_MIE
C       --> constant specific to the program SOS_ANGLES
C---------------------------------------------------------------------------
#define CTE_DEFAULT_OS_NB 80

C    Default value of the limit order of the expansion in Legendre
C    polynomials of the elements of the Fresnel matrix
C    Usually: CTE_DEFAULT_OS_NS = 2*CTE_DEFAULT_NBMU_LUM
C       --> constant specific to the program SOS_ANGLES
C---------------------------------------------------------------------------
#define CTE_DEFAULT_OS_NS 48

C    Default value of the limit order of the expansion in Fourier series
C    of the G function.
C    It is required that CTE_DEFAULT_OS_NM >= CTE_DEFAULT_OS_NB + CTE_DEFAULT_OS_NS
C       --> constant specific to the program SOS_ANGLES
C---------------------------------------------------------------------------
#define CTE_DEFAULT_OS_NM 128


C ===================================================
C 6.3 Limitation of the number of angles for the computation
C     of the files defining the angles to use
C ===================================================

C    Maximum number of "user" angles
C    complementing the Gauss angles
C       --> constant specific to the program SOS_ANGLES
C--------------------------------------------------------------
#define CTE_NBMAX_USER_ANGLES 20

C    Maximum number of positive angles for the computation of the angles
C    used by the codes (Gauss angles + user angles
C    + 1 = solar zenith angle)
C    Advice: corresponds to the max of CTE_MIE_NBMU_MAX and CTE_OS_NBMU_MAX
C       --> constant specific to the program SOS_ANGLES
C--------------------------------------------------------------
#define CTE_NBANGLES_MAX 100

C    Threshold on the difference between cosine(X) and cosine(TETAS) to
C    consider that TETAS = X
C       --> constant specific to the program SOS_ANGLES
C--------------------------------------------------------------
#define CTE_SEUIL_ECART_MUS 0.00001
