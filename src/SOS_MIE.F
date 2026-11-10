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
C* FILE: SOS_MIE.F
C* PROJECT: Successive orders of scattering
C* ROLE: This module manages the procedures for the computations of the
C*       MIE theory.
C*
C* AUTHOR:
C*   Initial code: Laboratoire d'Optique Atmospherique (LOA, Villeneuve d'Ascq).
C*   Code rework and development: CS GROUP France, B.Lafrance
C* DATE: 30/04/02
C*
C* MOD:VERSION:1.0 : Upgrade of the LOA OS code
C*             1.1 : Adjustment of the comments
C*
C* MOD:VERSION:2.0: 10/01/2010 
C*   - Angle management:
C*         * Removal of the reading of a Gauss angle file
C*         * Input arguments passed:
C*                 MIE_NBMU : number of useful angles
C*                 RMU      : array of the cosines of the angles
C*                 CHR      : array of the weights of the angles
C*
C*   - Removal of the constant SOS_LENDIR
C*
C*   - Change of the size of the phase function arrays:
C*     ==> Changed from CTE_MIE_NBMU to CTE_MIE_NBMU_MAX.
C*     ==> Loop on the angles up to MIE_NBMU, equal to the actual number of angles used.
C*
C*   - SOS_FPHASE_MIE:
C*         * Added input argument: MIE_NBMU (nb of useful angles)
C*         * Change of the record line of the arrays
C*           of phase functions per value of the size parameter Alpha
C*           ==> limited to the useful bounds -MIE_NBMU:MIE_NBMU
C* 
C* MOD:VERSION:2.1: 24/06/2010 
C*     FORMAT labels ,X, changed to ,1X, to avoid compilation errors.
C*
C* MOD:VERSION:2.2: 06/06/2013
C*
C*   -  Change of the value of Pi, declared as a constant in the file
C*      (INCTE_PI replaces CTE_PI)
C*  
C* MOD:VERSION:2.3: 24/02/2016 
C*    - Minor adjustments to strictly respect the 72 columns
C*  
C* MOD:VERSION:2.4: 30/01/2019
C*    - Management of the no-trace option with the name "NO_LOG_FILE"
C*
C* MOD:VERSION:2.5: 12/06/2020
C*    - Removal of the status "NEW" when opening the MIE file
C*      ==> Allows a restart that overwrites the existing file, if any.
C*
C* MOD:VERSION:2.6: 31/08/2022
C*    - Constants renamed from SOS_* to CTE_*
C*    - Traces converted to English.
C*
C* MOD:VERSION:2.7: 26/09/2023
C*    - Change of the syntax of a zeroing loop in routine SOS_XERO for
C*      compatibility with Fortran 2018 for the gcc/10.2.0 compiler
C*         
C*                  
C******************************************************************************

C----------------------------------------------------------------------------
C Definition of constants
C---------------------------------------------------------------------------- 
C Constants used:
C    CTE_LENFIC2 : Length of file names with directory tree.
C    INCTE_PI : Value of PI
C    CTE_MIE_NBMU_MAX : Maximum number of positive Gauss angles
C    CTE_MIE_DIM : Dimension of the arrays
C                              
C----------------------------------------------------------------------------
#include "SOS.h"
#define INCTE_PI DACOS(-1.D+00)

C==============================================================================
C PROCEDURE: SOS_MIE
C ==========      
C      This procedure applies the MIE theory to spherical particles.
C
C      The chemical nature of the particles is characterised by their refractive
C      index relative to the surrounding medium: m.
C      m = mr + i * mi, with "mr" the real part and "mi" the imaginary part.
C      The imaginary part is given negative. It is responsible for the
C      absorption effects.
C
C      The Mie computations are made versus the size parameter:
C      alpha = 2*pi * r / wa, where r is the particle radius and wa is the
C      wavelength of the radiation.
C
C      The step on alpha is:
C         if        alpha  <= 0.1  ===> step = 0.0001
C         if 0.1  < alpha  <= 1.0  ===> step = 0.001
C         if 1.0  < alpha  <= 10   ===> step = 0.01
C         if 10   < alpha  <= 30   ===> step = 0.05
C         if 30   < alpha  <= 100  ===> step = 0.1
C         if 100  < alpha          ===> step = 1.0
C
C
C      This program computes, for the index m:
C      
C         - The extinction efficiency factor: Qext(alpha)
C           from alpha0 to alphaF
C         - The scattering efficiency factor:  Qsca(alpha)
C           from alpha0 to alphaF
C         - The asymmetry factor: G
C         - The phase functions: Imie(alpha,mu),
C                                     Qmie(alpha,mu),
C				          Umie(alpha,mu)
C 	     from alpha0 to alphaF and for mu, cosine of the scattering angle, of the
C           Gauss angles between 0 and pi.
C
C Description of parameters
C ----------------------------
C      MIE_NBMU (I4) : (I) Number of actually useful (positive) angles
C	RMU(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX) (double) : (I) Cosines of the angles.
C	CHR(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX) (double) : (I) Weights of the angles.
C	RN (double) (I) : real part of the refractive index.
C	IN (double) (I) : imaginary part of the refractive index.
C	ALPHAO (double) (I) : value of alpha from which
C			        the computation starts.
C	ALPHAF (double) (I) : value of alpha at which the computation stops.
C	FICMIE (CHARACTER*CTE_LENFIC2) (I) : path and name of the MIE file.
C	FICLOG (CHARACTER*CTE_LENFIC2) (I) : name of the MIE TRACE file.
C	IER (I4) (O) : error code = 0 if no error, =1 otherwise
C
C
C Result provided by the procedure
C --------------------------------
C   MIE file
C   Trace file (optional)
C
C
C Environment variable
C ------------------------
C   None
C 
C Description of the files used
C ---------------------------------
C   -->Content of the MIE file (created)
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
C
C Common used:
C --------------
C  None
C
C
C Error cases:
C ------------
C     -Incorrect decoding of parameters
C     -Incorrect program parameters ALPHA0 >ALPHAF
C     -Opening of a file
C      Error if the MIE file already exists
C     -Reading/writing in a file
C
C     Message displayed on screen, program stop and
C     return of status 1 to the shell
C        
C   
C
C  Subroutines used:
C --------------------------
C  This program calls the routines:
C      - SOS_XERO
C      - SOS_FPHASE_MIE
C
C==============================================================================

      SUBROUTINE SOS_MIE(MIE_NBMU,RMU,CHR,RN,IN,ALPHAO,ALPHAF,
     &                   FICMIE,FICLOG,IER)


      IMPLICIT NONE


C* Definition of variables
C----------------------------------------------------------------------------        
      DOUBLE PRECISION RN	! real part of the refractive index.
      DOUBLE PRECISION IN	! imaginary part of the refractive index.
      DOUBLE PRECISION ALPHAO   ! value of alpha from which
      				    ! the computation starts.
      DOUBLE PRECISION ALPHAF   ! value of alpha at which the computation stops.
                                          
      DOUBLE PRECISION RMU(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)   !Cosines of the angles.
      DOUBLE PRECISION CHR(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)   !Weights of the angles.

      DOUBLE PRECISION CNA(-1:CTE_MIE_DIM)   !Real function Cn(alpha).
      DOUBLE PRECISION SNA(-1:CTE_MIE_DIM)   !Real function Sn(alpha).
           
      DOUBLE PRECISION RGNA(-1:CTE_MIE_DIM)  !Real part of the function Gn(alpha).
      DOUBLE PRECISION IGNA(-1:CTE_MIE_DIM)  !Imaginary part of Gn(alpha).
      
      DOUBLE PRECISION RDNA(0:CTE_MIE_DIM)   !Real function Dn(alpha).
      DOUBLE PRECISION RDNB(0:CTE_MIE_DIM)   !Real part of the function Dn(alpha*m).
      DOUBLE PRECISION IDNB(0:CTE_MIE_DIM)   !Imaginary part of Dn(alpha*m).
      
      DOUBLE PRECISION RA(0:CTE_MIE_DIM)     !Real part of the function An(alpha).
      DOUBLE PRECISION IA(0:CTE_MIE_DIM)     !Imaginary part of An(alpha).
      DOUBLE PRECISION RB(0:CTE_MIE_DIM)     !Real part of the function Bn(alpha).
      DOUBLE PRECISION IB(0:CTE_MIE_DIM)     !Imaginary part of Bn(alpha).
      
      

      
      DOUBLE PRECISION QEXT		!Extinction efficiency factor.
      DOUBLE PRECISION QSCA		!Scattering efficiency factor.
      DOUBLE PRECISION G		!Asymmetry factor.
      
      DOUBLE PRECISION ALPHA	!Current value of the size parameter.
      DOUBLE PRECISION PAS		!Step between two values of ALPHA.
      
      DOUBLE PRECISION RBETA	!Real part of m*alpha.
      DOUBLE PRECISION IBETA       !Imaginary part of m*alpha.
      

      DOUBLE PRECISION X
      DOUBLE PRECISION Y
      DOUBLE PRECISION Z
      DOUBLE PRECISION W
      
      DOUBLE PRECISION X1
      DOUBLE PRECISION X2
      DOUBLE PRECISION X3
      DOUBLE PRECISION X4
      DOUBLE PRECISION X5
      DOUBLE PRECISION X6
      DOUBLE PRECISION X7
           
      DOUBLE PRECISION Y1
      DOUBLE PRECISION Y2
      DOUBLE PRECISION Y3
      DOUBLE PRECISION Y4
      DOUBLE PRECISION Y5
      DOUBLE PRECISION Y6
      DOUBLE PRECISION Y7
      DOUBLE PRECISION Y8
      
      DOUBLE PRECISION Z3
      DOUBLE PRECISION Z4
      DOUBLE PRECISION Z5
      DOUBLE PRECISION Z6
      DOUBLE PRECISION Z7
      DOUBLE PRECISION Z8
      DOUBLE PRECISION Z9
      
      DOUBLE PRECISION T
      DOUBLE PRECISION Q
      
      DOUBLE PRECISION W6		!Value of 2 / ALPHA**2
      DOUBLE PRECISION XX
      DOUBLE PRECISION YY
      DOUBLE PRECISION ZZ
      DOUBLE PRECISION TT
      
      DOUBLE PRECISION A2

      INTEGER*4 MIE_NBMU   		! Number of actual angles of the simulations
      INTEGER*4 UN			!Value 1 or -1
      INTEGER*4 TEST
      
      INTEGER*4 N1			! Maximum order of computation of the series Dn and Sn
      INTEGER*4 N2			! Maximum order of computation of the series Gn and Cn
      INTEGER*4 N2P1			
      INTEGER*4 M
       
      INTEGER*4 J				
      INTEGER*4 I
      INTEGER*4 N
      INTEGER*4 DIM                   ! dimension
      INTEGER*4 IFIN                  ! nb of significant characters of a string
      
      LOGICAL TRACE                   ! = true if writing in the trace file
       
      CHARACTER*CTE_LENFIC2 FICMIE    ! access path to the MIE file, name of the MIE file
      CHARACTER*CTE_LENFIC2 FICLOG    ! name of the MIE TRACE file
      
      INTEGER*4 IER       	      ! error code =0 if no error, =-1 otherwise
      
                       
C* Initialisation of the error code
C------------------------------------------------------
      IER=0     

C* Compatibility test between alphaF and the maximum order
C* of the expansions for the computation of the
C* Ricatti-Bessel functions (and of their derivatives)
C------------------------------------------------------
      IF (DINT(ALPHAF+ALPHAF+20).GT.CTE_MIE_DIM) GOTO 997



C* Reading of the input parameters/arguments and decoding
C------------------------------------------------------      
      IF (FICLOG.EQ.'NO_LOG_FILE') THEN
	TRACE=.FALSE.
      ELSE
        TRACE=.TRUE.
      END IF     
      
C* Opening of the files
C------------------------------------------------------
      OPEN (1,FILE=FICMIE,FORM='UNFORMATTED',ERR=911)            
            
C     Trace file
      IF (TRACE) THEN
         OPEN (99,FILE=FICLOG,ERR=913)
	  WRITE(99,*,err=923) ' '
         WRITE(99,*,err=923) 'Calculation of MIE file'   
      END IF
 
               
C* Initialisation of the arrays to 0
C------------------------------------------------------  
      DIM=CTE_MIE_DIM 
      CALL SOS_XERO(CNA,-1,DIM)
      CALL SOS_XERO(SNA,-1,DIM)
      CALL SOS_XERO(RGNA,-1,DIM)  
      CALL SOS_XERO(IGNA,-1,DIM)
      CALL SOS_XERO(RDNA,0,DIM) 
      CALL SOS_XERO(RDNB,0,DIM)           
      CALL SOS_XERO(IDNB,0,DIM)    
      
                                                       
C* Trace file
C------------------------------------------------------  
      IF (TRACE) THEN
        WRITE(99,*,err=923)'SIMULATION PARAMETERS'
        WRITE(99,*,err=923)
     *      ' Refractive index (real part): ',RN
        WRITE(99,*,err=923)
     *      ' Refractive index (imaginary part): ',IN
        WRITE(99,*,err=923)' Limit values of the size parameter:',
     *           ALPHAO,' and ',ALPHAF
        WRITE(99,*,err=923)' '
        WRITE(99,*,err=923)' List of used angles  : '
	WRITE(99,*,err=923) "   Index, mu, angle (degrees), weight"
	DO J=-MIE_NBMU,MIE_NBMU
	    WRITE(99,222,ERR=923) J, RMU(J),
     &	                 DACOS(RMU(J))*180./INCTE_PI, CHR(J)
	ENDDO
	WRITE(99,*,err=923)' '
	WRITE(99,*,err=923)'RESULTS'
	WRITE(99,*,err=923)' Size parameter alpha'
	WRITE(99,*,err=923)
     *      ' Mie extinction efficiency factor Qext'
	WRITE(99,*,err=923)
     *      ' Mie scattering efficiency factor Qsca'
	WRITE(99,*,err=923)' Anisotropy factor G'
	WRITE(99,*,err=923)' '
        WRITE(99,10,err=923)'  alpha','Qext','Qsca','G'
      END IF

C* Writing of the header of the output file
C------------------------------------------------------                         
      WRITE(1,err=921)RN,IN,ALPHAF,MIE_NBMU
      
  
  
C* Loop on the values of the size parameter ALPHA from
C* ALPHA0 to ALPHAF with variable step. (label 7777)
C------------------------------------------------------          
      ALPHA=ALPHAO
           
 7777 CONTINUE   
 
C* Estimation of the step for the next value of ALPHA
C------------------------------------------------------
         PAS=0.0001
         IF(ALPHA.gt.0.1)  PAS =0.001
         IF(ALPHA.gt.1.00) PAS =0.01
         IF(ALPHA.gt.10.)  PAS =0.05
         IF(ALPHA.gt.30.)  PAS =0.10
         IF(ALPHA.gt.100.) PAS =1.00
 
 
      
C* Computation of the derived Ricatti-Bessel functions:
C* Gn(alpha), Dn(alpha), Dn(alpha*m) - real and imaginary parts
C*
C* Computation of the real functions Sn(alpha) and Cn(alpha)
C------------------------------------------------------    

C*  Limited expansion orders: N1 for Dn and Sn
C*                            N2 for Gn and Cn
                 
         N1   = DINT(ALPHA+ALPHA+20)
	  IF(N1.GT.CTE_MIE_DIM) GOTO 997  
         N2   = DINT(ALPHA+ALPHA+5)  
         N2P1 = N2+1  
	 
	 

C*  Recurrence on Gn(alpha): real part RGNA and imaginary part IGNA,
C*  and on Cn(alpha) from n=1 to N2.
C------------------------------------------------------  
      
C*   Initialisation of the series Cn (real CNA)
C*   and Gn (real part RGNA and imaginary part IGNA)

         CNA(-1) = -DSIN(ALPHA)   
	  CNA(0)  =  DCOS(ALPHA)
	     
         RGNA(-1) =  0.D+00
         RGNA(0)  =  0.D+00         
         IGNA(-1) =  0.D+00
         IGNA(0)  = -1.D+00


C*  Recurrence loop
    
         DO 25 I=1,N2
            X=RGNA(I-1)
            Z=I/ALPHA
            Y=IGNA(I-1)
            W=((Z-X)*(Z-X)+(Y*Y)) 
            RGNA(I)=(Z-X)/W-Z 
            IGNA(I)=Y/W

            CNA(I)=(2*I-1.D+00)*CNA(I-1)/ALPHA-CNA(I-2)	
	    
C*    Check:
C*      If Cn(alpha) does not diverge => the recurrence goes on.
C*      Otherwise, the limit orders N2 and N1 are redefined, then the
C*      recurrence loop is left.
            IF(CNA(I).LT.1.D+304) GO TO 25	    
              N2=I 
              N2P1=I+1
              N1=I+15
            GOTO 100
	    
   25    CONTINUE
  100    CONTINUE



C*  Recurrence on Dn(alpha)   : real function RDNA
C*                 Dn(alpha*m) : real part RDNB and imaginary part IDNB
C*                 Sn(alpha)   : real function SNA
C*  Loop from N1 to 0
C------------------------------------------------------         

         RBETA=RN*ALPHA       
         IBETA=IN*ALPHA    
         X1=RBETA*RBETA+IBETA*IBETA 
         X2=RBETA/X1
         X3=IBETA/X1
	 
C*   Initialisation at N1
         RDNA(N1)=0.D+00
         RDNB(N1)=0.D+00
         IDNB(N1)=0.D+00
      
         SNA(N1)=0.D+00 
         SNA(N1-1)=1.D+00

C*   Recurrence loop from N1-1 to 0
         DO 30 I=N1-1,0,-1
	    
            X=RDNB(I+1) 
            Y=IDNB(I+1)
            Z=X+(I+1.D+00)*X2
            W=Y-(I+1.D+00)*X3 
            X4=Z*Z+W*W 
            RDNB(I)=(I+1.D+00)*X2-Z/X4  
            IDNB(I)=-(I+1.D+00)*X3+W/X4
	    
            Z=(I+1.D+00)/ALPHA
            X=RDNA(I+1)  
            RDNA(I)=Z-1.D+00/(X+Z)
	    
            SNA(I-1)=(2.D+00*I+1.D+00)*SNA(I)/ALPHA-SNA(I+1)

C*   If the values of Sn become too large, all the computed terms
C*   are normalised by Sn-1.
            IF (SNA(I-1).GT.1D+304) THEN    	     	
                TEST=I-1 
                X=SNA(TEST)
                DO 35 J=TEST,N2
                   SNA(J)=SNA(J)/X 
   35           CONTINUE
            ENDIF
	    		  
   30    CONTINUE 
 
 
C*   Normalisation of Sn(alpha) by sin(alpha)/S0(alpha)
C*   to satisfy the condition n=0: S0(alpha) = sin(alpha)
         Q=SNA(0)/DSIN(ALPHA)    
                                                           
         DO 45 I=0,N2       
             SNA(I)=SNA(I)/Q 
   45    CONTINUE
   
         
          


C* Computation of the functions An(alpha) and Bn(alpha):
C* real (RA, RB) and imaginary (IA, IB) parts
C------------------------------------------------------    
         UN=1 
	                                        
         DO 50 I=1,N2     
            X1=SNA(I)
            X2=CNA(I)
            X3=RDNB(I)
            X4=IDNB(I)                
            X5=RDNA(I)
            X6=RGNA(I)         
            X7=IGNA(I)           
            Y1=X3-RN*X5          
            Y2=X4-IN*X5         
            Y3=X3-RN*X6+IN*X7         
            Y4=X4-RN*X7-IN*X6       
            Y5=RN*X3-IN*X4-X5        
            Y6=IN*X3+RN*X4        
            Y7=RN*X3-IN*X4-X6       
            Y8=IN*X3+RN*X4-X7       
            Z4=Y2*Y3-Y1*Y4       
            Z3=Y1*Y3+Y2*Y4        
            Z5=X1*X1+X2*X2       
            Z6=Y3*Y3+Y4*Y4
            Z7=Y5*Y7+Y6*Y8       
            Z8=Y6*Y7-Y5*Y8         
            Z9=Y7*Y7+Y8*Y8          
            Q=(I+I+1.D+00)/I/(I+1.D+00)*UN

C*   Test on X2 (function Cn(alpha)): tends to diverge for large n.
C*   Z5 tends to infinity as X2**2. The value of X2 / Z5 then tends to 0.
C*   This test avoids a numerical problem of Inf / Inf ratio.
	    IF (X2.GT.1D+300) THEN
		Y1=0.D+00
		Y2=0.D+00
		Y3=0.D+00
		Y4=0.D+00		
	    ELSE
                Y1=X1*(X1*Z3+X2*Z4)/Z5/Z6          
                Y2=X1*(X1*Z4-X2*Z3)/Z5/Z6         
                Y3=X1*(X1*Z7+X2*Z8)/Z5/Z9          
                Y4=X1*(X1*Z8-X2*Z7)/Z5/Z9  
	    ENDIF
            
            RA(I)=Y2*Q               
            IB(I)=Y3*Q           
            Q=-Q              
            RB(I)=Y4*Q         
            IA(I)=Y1*Q        
            UN=-UN    
   50    CONTINUE
    	 
                                                                 
         RA(0)=0.D+00
         IA(0)=0.D+00
         RB(0)=0.D+00
         IB(0)=0.D+00
         RA(N2P1)=0.D+00 
         IA(N2P1)=0.D+00
         RB(N2P1)=0.D+00
         IB(N2P1)=0.D+00


C* Computation of the scattering (Qsca)
C* and extinction (Qext) efficiency factors.
C* Computation of the asymmetry factor (G).
C------------------------------------------------------
         QEXT=0.D+00
         QSCA=0.D+00
         G=0.D+00
         J=-1
         X=RA(1)      
         Y=IA(1)       
         Z=RB(1)        
         T=IB(1)
                                                                         
         DO 55 N=1,N2         
            M=N+1
            XX=RA(M) 
            YY=IA(M)
            ZZ=RB(M) 
            TT=IB(M) 
            A2=(N+1.D+00)
            QEXT=QEXT+N*A2*J*(Y-T) 
            QSCA=QSCA+N*N*A2*A2/(N+A2)*(X*X+Y*Y+Z*Z+T*T)    
            J=-J  
            G=G-A2*N/(A2+N)*(N*(A2+1.)*(A2+1.)/(2.*N+3.)*
     1        (Y*YY+X*XX+T*TT+Z*ZZ)+Y*T+X*Z)
            X=XX
            Y=YY
            Z=ZZ
            T=TT
   55    CONTINUE 
                                                             
         W6=2.D+00/ALPHA/ALPHA  
         QEXT=W6*QEXT         
         QSCA=W6*QSCA             
         G=4.D+00*G/QSCA/ALPHA/ALPHA

C* Call of the routine computing the Mie phase functions
C* and storing the results for ALPHA.
C------------------------------------------------------
         CALL SOS_FPHASE_MIE(MIE_NBMU,RMU,RA,IA,RB,IB,ALPHA,
     *                      QEXT,QSCA,G,N2,TRACE,IER) 
         IF (IER.NE.0) GOTO 996   ! error case
      
C* Next value of the size parameter ALPHA (stop if the
C* limit ALPHAF is exceeded).
C------------------------------------------------------
         ALPHA=ALPHA+PAS 
      IF (ALPHA.LE.ALPHAF) GOTO 7777 
                                                   
 
C* File closing
C-------------------
      CLOSE(1)
      CLOSE(4)

      
C* Nominal end
C-------------------
      GOTO 9999
      
C* Error case and return of status 1 to the shell
C----------------------------------------------
  911 WRITE(6,*) 'Erreur a l''ouverture du fichier MIE'
      IER=-1  
      GOTO 9999
  913 WRITE(6,*) 'Erreur a l''ouverture du fichier trace MIE'
      IER=-1  
      GOTO 9999        
  921 WRITE(6,*) 'Erreur d''Ã©criture sur le fichier MIE'
      IER=-1  
      GOTO 9999
  922 WRITE(6,*) 'Erreur de lecture sur le fichier angles de Gauss'
      IER=-1  
      GOTO 9999
  923 WRITE(6,*) 'Erreur d''Ã©criture sur le fichier trace MIE'
      IER=-1  
      GOTO 9999    
  996 WRITE(6,*) 'Erreur dans la routine SOS_PHASE_MIE'
      IER=-1
      GOTO 9999
  997 WRITE(6,*) 'Valeur AlphaMax trop grande devant CTE_MIE_DIM'
      IER=-1
      GOTO 9999      
  
C* Format
C---------   
   10 FORMAT(A7,17X,A4,16X,A4,15X,A1)          
  222 FORMAT(I4,1X,D21.14,2X,F7.3,2X,D21.14,2X)   
                                                   
 9999 CONTINUE
 
      RETURN
      END    !END OF PROCEDURE SOS_MIE
      

 
                                                                      

C==============================================================================
C PROCEDURE: SOS_XERO
C ==========      
C  This procedure sets to zero the elements of an array IT(N:M)
C
C
C Description of parameters
C --------------------------
C    N,M (I4) (I)  Dimensions of the array IT
C    IT(N:M) (double) (I/O)  Array to initialise
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
      SUBROUTINE SOS_XERO(IT,N,M) 
   
      IMPLICIT NONE
      
C* Definition of variables
C---------------------------
      INTEGER*4 N,M			! Dimensions of the array
      DOUBLE PRECISION IT(N:M)		! Array
               
      INTEGER*4 I			!Array index
      
C*   Initialisation of the array
C-----------------------------           
      DO I=N,M  
         IT(I)=0.D+00
      ENDDO
      
      RETURN 
      END  	!END OF PROCEDURE SOS_XERO
      
 
                         
                                                        
C==============================================================================
C PROCEDURE: SOS_FPHASE_MIE
C ==========      
C      This procedure computes the phase functions Imie(mu), Qmie(mu),
C      Umie(mu) for a given value of the size parameter ALPHA and
C      for the Gauss mu of the array RMU.
C
C      The results are recorded (unformatted):
C        - Size parameter (real),
C        - Extinction efficiency factor (real),
C        - Scattering efficiency factor (real),
C        - Asymmetry factor (double precision),
C        - Array of the phase function Imie(mu) (real),
C        - Array of the polarised phase function Qmie(mu) (real),
C        - Array of the polarised phase function Umie(mu) (real).
C
C    --> The arrays have a useful size limited to MIE_NBMU:
C        they are written with this limitation.
C
C     
C Description of parameters
C --------------------------
C    MIE_NBMU (I4)   (I)   Number of actually useful (positive) angles
C    RMU(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX) (double) (I) Cosines of the Gauss angles
C    RA(0:CTE_MIE_DIM)  (double) (I)   Real part of the function An(alpha).
C    IA(0:CTE_MIE_DIM)  (double) (I)   Imaginary part of An(alpha).
C    RB(0:CTE_MIE_DIM)  (double) (I)   Real part of the function Bn(alpha).
C    IB(0:CTE_MIE_DIM)  (double) (I)   Imaginary part of Bn(alpha).
C    ALPHA (double)  (I)   Value of the size parameter.
C    KMA1  (double)  (I)   Extinction efficiency factor.
C    KMA2  (double)  (I)   Scattering efficiency factor.
C    G     (double)  (I)   Asymmetry factor.
C    N2    (I4)      (I)   Maximum order of summation on N.
C    TRACE (logical) (I)   true if writing in the trace file
C    IER (I4) (O) : error code = 0 if no error, =1 otherwise
C
C
C Result stored in the MIE file:
C ------------------------------------
C    AY1   (R4)	        Value of the size parameter ALPHA.
C    AY2   (R4)	        Extinction efficiency factor.
C    AY3   (R4)	        Scattering efficiency factor.
C    AY4   (double)     Asymmetry factor.
C    IMIE(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)    (R4)	      
C                 phase function: Imie(mu) for ALPHA.
C    QMIE(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)    (R4)	     
C                  polarised phase function: Qmie(mu).
C    UMIE(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)    (R4)     
C                  polarised phase function: Umie(mu).
C      
C
C Common used:
C --------------
C     None
C
C Error cases:
C ------------
C     Message displayed on screen, program stop and
C     return of status 1 to the shell
C
C==============================================================================
      SUBROUTINE SOS_FPHASE_MIE(MIE_NBMU,RMU,RA,IA,RB,IB,ALPHA,
     *                          KMA1,KMA2,G,N2,TRACE,IER)
 

      IMPLICIT NONE
      
C* Definition of variables
C----------------------------------------------------------------------------
      DOUBLE PRECISION RMU(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX) 
                                        !Cosines of the Gauss angles
      
      DOUBLE PRECISION ALPHA       !Value of the size parameter.
      DOUBLE PRECISION KMA1	   !Extinction efficiency factor.
      DOUBLE PRECISION KMA2        !Scattering efficiency factor.
      DOUBLE PRECISION G 	   !Asymmetry factor.


      DOUBLE PRECISION RA(0:CTE_MIE_DIM)     !Real part of the function An(alpha).
      DOUBLE PRECISION IA(0:CTE_MIE_DIM)     !Imaginary part of An(alpha).
      DOUBLE PRECISION RB(0:CTE_MIE_DIM)     !Real part of the function Bn(alpha).
      DOUBLE PRECISION IB(0:CTE_MIE_DIM)     !Imaginary part of Bn(alpha).
      DOUBLE PRECISION AR	!Value of the array RA for a given n.
      DOUBLE PRECISION AI       !Value of the array IA for a given n.
      DOUBLE PRECISION BR       !Value of the array RB for a given n.
      DOUBLE PRECISION BI       !Value of the array IB for a given n.
            
      DOUBLE PRECISION COEF
      DOUBLE PRECISION X        !Value mu = -cos(tetad)
      				     !with tetad the scattering angle.
				           
      DOUBLE PRECISION PIM      !Function PI(mu) of order n-1.
      DOUBLE PRECISION PIV      !Function PI(mu) of order n.
      DOUBLE PRECISION PIP     	!Function PI(mu) of order n+1.
      DOUBLE PRECISION TAU      !Function TAU(mu) of order n.
      
      DOUBLE PRECISION RES1     !Real part of the function S1(alpha,m,mu).
      DOUBLE PRECISION IMS1     !Imaginary part of S1(alpha,m,mu).
      DOUBLE PRECISION RES2     !Real part of the function S2(alpha,m,mu).
      DOUBLE PRECISION IMS2     !Imaginary part of S2(alpha,m,mu).
      
      DOUBLE PRECISION Y1
      DOUBLE PRECISION Y2
      DOUBLE PRECISION Y3
      DOUBLE PRECISION Y4
      DOUBLE PRECISION AY4	!Asymmetry factor.
      
      REAL*4 AY1		!Value of the size parameter ALPHA.
      REAL*4 AY2		!Extinction efficiency factor.
      REAL*4 AY3		!Scattering efficiency factor.
      REAL*4 IMIE(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)     
                             !phase function: Imie(mu) for ALPHA.
      REAL*4 QMIE(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)     
                             !polarised phase function: Qmie(mu).
      REAL*4 UMIE(-CTE_MIE_NBMU_MAX:CTE_MIE_NBMU_MAX)
                             !polarised phase function: Umie(mu).

      LOGICAL TRACE          ! = true if writing in the trace file

      INTEGER*4 MIE_NBMU     ! Number of actual angles of the simulations
      INTEGER*4 N2      !Maximum order of summation on N
      INTEGER*4 J
      INTEGER*4 N
      
      INTEGER*4 IER        ! error code =0 if no error, =-1 otherwise


      AY1=ALPHA           
      AY2=KMA1      
      AY3=KMA2    
      AY4=G
      COEF=2.D+00/KMA2/ALPHA**2
      
      DO 1 J=-MIE_NBMU,MIE_NBMU              
         X=-RMU(J)            
         PIM=0.D+00 
         PIV=1.D+00  
         TAU=X          
         RES1=0.D+00
         RES2=0.D+00
         IMS1=0.D+00
         IMS2=0.D+00
	 
         DO 2 N=1,N2           
            AI=IA(N)
            BI=IB(N)
            AR=RA(N)
            BR=RB(N)
            RES1=RES1-AI*PIV-BI*TAU  
            RES2=RES2+AI*TAU+BI*PIV  
            IMS1=IMS1+AR*PIV+BR*TAU 
            IMS2=IMS2-AR*TAU-BR*PIV  
            PIP=((2.D+00*N+1.D+00)*X*PIV-(N+1.D+00)*PIM)/N 
            PIM=PIV
            PIV=PIP  
            TAU=(N+1.D+00)*X*PIV-(N+2.D+00)*PIM 
    2    CONTINUE 
                                                                     
         Y1=RES1*RES1+IMS1*IMS1  
         Y2=RES2*RES2+IMS2*IMS2 
         Y3=2.D+00*RES2*RES1   
         Y4=2.D+00*IMS2*IMS1
         IMIE(J)=COEF*(Y1+Y2)      
         QMIE(J)=COEF*(Y2-Y1)
         UMIE(J)=COEF*(Y3+Y4)
    1 CONTINUE    
    
C* Recording of the results for ALPHA
C*    Recording of the phase function arrays
C*    with record bounds limited to the useful
C*    dimension -CTE_MIE_NBMU:CTE_MIE_NBMU
C------------------------------------------------------
      WRITE(1,err=921) AY1,AY2,AY3,AY4,
     &                (IMIE(J),J=-MIE_NBMU,MIE_NBMU),
     &                (QMIE(J),J=-MIE_NBMU,MIE_NBMU),
     &                (UMIE(J),J=-MIE_NBMU,MIE_NBMU)
      
C* Trace file: alpha,Qext,Qsca,G
C------------------------------------------------------  
      IF (TRACE) THEN
	   WRITE(99,20,err=923) AY1,AY2,AY3,AY4
      END IF


C* Nominal end
C-------------------
      GOTO 9999
      
C* Error case and return of status -1 to the calling program
C----------------------------------------------------------
  921 WRITE(6,*) 'Erreur d''Ã©criture sur le fichier MIE'
      IER=-1 
      GOTO 9999
  923 WRITE(6,*) 'Erreur d''Ã©criture sur le fichier trace MIE'
      IER=-1 
      GOTO 9999

C* Format 
C----------        
20    FORMAT(1X,E15.7,5X,E15.7,5X,E15.7,5X,E15.7)  

 9999 CONTINUE
 
      RETURN                                                                    
      END         !END OF PROCEDURE SOS_FPHASE_MIE

