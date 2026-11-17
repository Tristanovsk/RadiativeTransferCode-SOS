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
C* FILE: SOS_AGGREGATE.F
C* PROJECT: Successive orders of scattering
C* ROLE: Weighting and aggregation of the outputs of SOS.F for each exponential of the CKD method
C*      - Aggregation of the parameters I Q U for each order of the Fourier series decomposition
C*	 - Aggregation of the atmospheric transmissions
C*	 - Aggregation of the optical thicknesses
C*
C* AUTHOR:
C*   CS GROUP France, X.Lenot, B.Lafrance
C* DATE: 06/06/2013
C*
C* MOD:VERSION:6.0: File creation to adapt the OS to the CKD method (atmospheric absorption)
C*
C* MOD:VERSION:6.1: 01/09/2020
C*     The name of the temporary file FICOS_AGG_TMP is passed as argument to allow
C*     several simultaneous runs from the same directory.
C*
C* MOD:VERSION:6.2: 31/08/2022
C*    - Constants renamed from SOS_* to CTE_*
C*    - Removal of the arguments TRACE and IDLOG (not used)
C*
C******************************************************************************

C----------------------------------------------------------------------------
C Definition of constants
C---------------------------------------------------------------------------- 
C Constants used:
C    CTE_LENFIC2 : Length of file names with directory tree.
C    CTE_OS_NBMU_MAX : Maximum number of positive angles for solving the RTE
C----------------------------------------------------------------------------
#include "SOS.h"
#define INCTE_PI DACOS(-1.D+00)

C==============================================================================
C PROGRAM: SOS_AGGREGATE
C ==========      
C      Program that weights the parameters I, Q, U, the transmissions and optical thicknesses
C      computed for one exponential of the CKD method, and aggregates them with the same parameters
C      computed previously.
C
C Input data of the program
C ------------------------------
C    -->Parameters/arguments of the program
C
C             NBMU (I4) : (I) Number of actual angles for radiance computations
C
C             AIK (I4) : (I) Weighting coefficient of the processed CKD exponential
C
C             FICOS_TMP (CHARACTER) : (I) Complete name of the OS result file obtained previously
C                                     (directory + file name + extension)
C
C             TTOT_TRONC_TMP (DOUBLE) : (I) Total optical thickness at ground
C				            (obtained previously) (truncated)
C
C             TTOT_VRAI_TMP (DOUBLE) : (I) Total optical thickness at ground
C				           (obtained previously) (not truncated)
C
C             TAUOUT_TMP (DOUBLE) : (I) Optical thickness at the output level for the upward field
C			               (obtained previously) (truncated)
C
C             TDIFMUS_TMP (DOUBLE) : (I) Diffuse transmission for the solar incidence direction
C      				      of the equivalent atmosphere (obtained previously)
C
C             TDIFMUG_TMP (DOUBLE) : (I) 1D array of 2*CTE_OS_NBMU_MAX + 1 elements
C			     	             containing the diffuse transmission of the equivalent atmosphere
C				             for the Gauss angles and the angles provided by the user
C   				         (obtained previously)
C
C             EMOINS_TMP (DOUBLE) : (I) Diffuse flux illuminating the ground, normalised by the incident solar flux at TOA.
C
C             EPLUS_TMP (DOUBLE) : (I)  Diffuse flux emerging at TOA, normalised by the incident solar flux at TOA.
C        
C             FICOS (CHARACTER) : (I) Complete name of the OS result file (aggregated results)
C                                     (directory + file name + extension)
C
C             TTOT_TRONC (DOUBLE) : (O) Total optical thickness at ground
C				           (aggregated results) (truncated)
C
C             TTOT_VRAI (DOUBLE) : (O) Total optical thickness at ground
C				          (aggregated results) (not truncated)
C
C             TAUOUT (DOUBLE) : (O) Optical thickness at the output level for the upward field
C			               (aggregated results) (truncated)
C
C             TDIFMUS (DOUBLE) : (O) Diffuse transmission for the solar incidence direction
C      				  of the equivalent atmosphere (aggregated results)
C
C             TDIFMUG (DOUBLE) : (O) 1D array of 2*CTE_OS_NBMU_MAX + 1 elements
C			     	         containing the diffuse transmission of the equivalent atmosphere
C				         for the Gauss angles and the angles provided by the user
C   				         (aggregated results)
C
C             EMOINS (DOUBLE) : (O) Diffuse flux illuminating the ground, normalised by the incident solar flux at TOA
C   				        (aggregated results)
C
C             EPLUS (DOUBLE) : (O) Diffuse flux emerging at TOA, normalised by the incident solar flux at TOA
C   				       (aggregated results)
C    
C      	IER (I4) : (O) error code =0 if no error, =-1 otherwise
C
C    -->Files
C   OS result file (Fourier series decomposition)
C
C
C Results provided by the program
C ----------------------------------
C   OS result file (Fourier series decomposition)
C   Diffuse transmissions
C   Optical thicknesses (ground and output level)
C
C   If one or several files with the same names already exist, they are overwritten
C
C
C 
C Description of the files used
C ---------------------------------
C
C   -->Content of the OS result file
C       Unformatted binary file, created by SOS_OS
C       It contains the Stokes parameters Q, U and I for an order S
C       of the Fourier series expansion, versus the angle mu
C       (double precision values).
C
C   -->Content of the temporary file FICOS_TMP
C       Same as the OS result file.
C       Input file with the parameters computed previously
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
C  None
C
C============================================================================

      SUBROUTINE SOS_AGGREGATE(NBMU, AIK, FICOS_TMP, TTOT_TRONC_TMP,
     &	 	    	      TTOT_VRAI_TMP, TAUOUT_TMP, TDIFMUS_TMP,
     &		    	      TDIFMUG_TMP, EMOINS_TMP, EPLUS_TMP,
     &                     FICOS_AGG_TMP, FICOS,
     &		    	      TTOT_TRONC, TTOT_VRAI, TAUOUT,            
     &	 	    	      TDIFMUS, TDIFMUG, EMOINS, EPLUS,
     &                     IER)
     

C* Definition of variables
C*------------------------------------------------------
				
     
      DOUBLE PRECISION TDIFMUS  !Diffuse transmission for the solar incidence direction
      				    !of the equivalent atmosphere
				
      DOUBLE PRECISION TDIFMUG(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
      				!Diffuse transmission for the Gauss angles of the equivalent atmosphere

      DOUBLE PRECISION EMOINS   ! Diffuse flux illuminating the ground, normalised by
      				    ! the incident solar flux at TOA.
                                
      DOUBLE PRECISION EPLUS    ! Diffuse flux emerging at TOA, normalised by
      				    ! the incident solar flux at TOA.
                                
      DOUBLE PRECISION TDIFMUS_TMP  !Diffuse transmission for the solar incidence direction
      				        !of the equivalent atmosphere (computed for a septet ik1*ik2*ik3*ik4*ik5*ik6*ik7)
				
      DOUBLE PRECISION TDIFMUG_TMP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
      				!Diffuse transmission for the Gauss angles of the equivalent atmosphere
				!  (computed for an octet ik1*ik2*ik3*ik4*ik5*ik6*ik7*ik8 of the CKD decomposition)

      DOUBLE PRECISION EMOINS_TMP   ! Diffuse flux illuminating the ground, normalised by
      				        ! the incident solar flux at TOA.
                                    ! (computed for an octet ik1*ik2*ik3*ik4*ik5*ik6*ik7*ik8)

      DOUBLE PRECISION EPLUS_TMP    ! Diffuse flux emerging at TOA, normalised by
      				        ! the incident solar flux at TOA.
                                    ! (computed for an octet ik1*ik2*ik3*ik4*ik5*ik6*ik7*ik8)

      DOUBLE PRECISION TTOT_TRONC ! Total optical thickness at ground (truncated).
      DOUBLE PRECISION TTOT_VRAI  ! Total optical thickness at ground (not truncated).

      DOUBLE PRECISION TTOT_TRONC_TMP ! Total optical thickness at ground (truncated).
                                      ! (computed for an octet ik1*ik2*ik3*ik4*ik5*ik6*ik7*ik8)

      DOUBLE PRECISION TTOT_VRAI_TMP  ! Total optical thickness at ground (not truncated).
					   ! (computed for an octet ik1*ik2*ik3*ik4*ik5*ik6*ik7*ik8)

      DOUBLE PRECISION TAUOUT	! Optical thickness at the
      					! output level of the results for the upward field (truncated).
      DOUBLE PRECISION TAUOUT_TMP 	! Optical thickness at the
      				  	! output level of the results for the upward field (truncated).
				 	! (computed for an octet ik1*ik2*ik3*ik4*ik5*ik6*ik7*ik8)


      ! Stokes parameters resulting from the aggregation
      DOUBLE PRECISION I3_RES(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Stokes parameter I for an order S
      				! of the Fourier series expansion,
				! versus the angle mu.
      DOUBLE PRECISION Q3_RES(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Stokes parameter Q.
      DOUBLE PRECISION U3_RES(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Stokes parameter U.

      ! Stokes parameters of the current iteration
      DOUBLE PRECISION I3_TMP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Stokes parameter I for an order S
      				! of the Fourier series expansion,
				! versus the angle mu.
      DOUBLE PRECISION Q3_TMP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Stokes parameter Q.
      DOUBLE PRECISION U3_TMP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Stokes parameter U.

      ! Stokes parameters from the previous iterations
      DOUBLE PRECISION I3_RES_TMP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Stokes parameter I for an order S
      				! of the Fourier series expansion,
				! versus the angle mu.
      DOUBLE PRECISION Q3_RES_TMP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Stokes parameter Q.
      DOUBLE PRECISION U3_RES_TMP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Stokes parameter U.

      DOUBLE PRECISION AIK  !Weighting coefficient

      DOUBLE PRECISION TRANS  !Atmospheric transmission

      INTEGER*4 NBMU    ! number of useful angles

      CHARACTER*CTE_LENFIC2 FICOS_TMP      ! Name of the temporary OS result file (computation for one iteration)

      CHARACTER*CTE_LENFIC2 FICOS_AGG_TMP   ! Name of the temporary OS result file (aggregated computation)

      CHARACTER*CTE_LENFIC2 FICOS       ! Name of the final OS result file (including all aggregations)

      INTEGER*4 IER        !error code =0 if no error, =-1 otherwise

      LOGICAL EX	! Checks the existence of a file

      INTEGER*4 I,J   	! Array index

      INTEGER*4 IOS16,IOS17   	! File statuses

      INTEGER*4 LFICOS ! length of the character string FICOS
      INTEGER*4 LFICOSTMP ! length of the character string FICOS_TMP
      INTEGER*4 LFICOSAGGTMP ! length of the character string FICOS_AGG_TMP

      INTEGER*4 LCOM     	        ! length of the character string COMMAND
      
      CHARACTER*CTE_LENCOM  COMMAND     ! Command string for system call.

 	  
      	        		 

C* Initialisation of values
C ---------------------------

C*    Init of the file status values
      IOS16=0
      IOS17=0


C*    Length of character strings
      LFICOS=INDEX(FICOS,' ')-1 
      LFICOSTMP=INDEX(FICOS_TMP,' ')-1 
      LFICOSAGGTMP=INDEX(FICOS_AGG_TMP,' ')-1 

C* Check of the existence of the reflection matrix file
C  ----------------------------------------------------------	
      INQUIRE(FILE=FICOS(1:LFICOS),EXIST=EX)

	  
C--------------------------------------------------------------
C
C       ******  AGGREGATION OF RESULTS   ******
C
C--------------------------------------------------------------
      
     	 	
 
C************************************************************************ 
C* Stokes parameters
C* The IQU values of the Fourier series decompositions are read
C* in the temporary file, weighted by the weighting coefficient AIK
C* and added to the IQU values stored in the result file
C************************************************************************ 
 
C* Opening of the files to read
C------------------------------------------------     
C     Opening of the temporary OS file
      OPEN(UNIT=16,FILE=FICOS_TMP,FORM='UNFORMATTED',
     &     STATUS='OLD',ERR=1101) 

      IF (EX.EQV..FALSE.) THEN  

C* Case of a first iteration: the result file does not exist yet
C* The value IOSTAT17 of the output file is set to a value
C* different from 0 => the file will not be read
	
	  IOS17=-1

      ELSE

C* Opening of the OS result file
          OPEN(UNIT=17,FILE=FICOS,FORM='UNFORMATTED',
     &         STATUS='OLD',ERR=1201) 
      
      ENDIF 
      

C* Opening of the output file
C------------------------------------------------     
C     Opening of the temporary output file
      OPEN(UNIT=18,FILE=FICOS_AGG_TMP,FORM='UNFORMATTED',
     &     STATUS='NEW',ERR=1301) 
     
C* Aggregation of the IQU parameters for each order of the Fourier series
C-------------------------------------------------------------------------
    
C*    Loop on order IS
      IS=0

C*    Initialisation of the values read in the files
      DO J=-NBMU,NBMU 
      
	     Q3_TMP(J)=0.
	     U3_TMP(J)=0.
	     I3_TMP(J)=0.
	     Q3_RES_TMP(J)=0.
	     U3_RES_TMP(J)=0.
	     I3_RES_TMP(J)=0.
	  
      ENDDO
      
C*    Reading of the TMP file that was written with
C*    array bounds limited to the useful
C*    dimension -NBMU:NBMU
   10 IF (IOS16.EQ.0) THEN
   
   	    READ(16,IOSTAT=IOS16,err=1102) (Q3_TMP(J),J=-NBMU,NBMU),
     & 		                           (U3_TMP(J),J=-NBMU,NBMU),
     &		                           (I3_TMP(J),J=-NBMU,NBMU)
 
      ENDIF
      
C*    Reading of the Result file that was written with
C*    array bounds limited to the useful
C*    dimension -NBMU:NBMU
      IF (IOS17.EQ.0) THEN
      
      	    READ(17,IOSTAT=IOS17,err=1202) (Q3_RES_TMP(J),J=-NBMU,NBMU),
     & 		                     (U3_RES_TMP(J),J=-NBMU,NBMU),
     &		                            (I3_RES_TMP(J),J=-NBMU,NBMU)

      ELSE
             ! Case where both files are finished
	     IF (IOS16.NE.0) GOTO 9

      ENDIF

C*    Both files are finished (end of file reached by the reads above):
C*    no record to write. Without this test, each aggregation after the
C*    first one appended an empty order to the result file, until the
C*    number of orders exceeded CTE_OS_NS_MAX (ERROR_4900 in SOS_TRPHI_MULTI
C*    when many CKD terms are combined).
      IF ((IOS16.NE.0).AND.(IOS17.NE.0)) GOTO 9
      
C*    Computation of the aggregation
      	
      DO J=-NBMU,NBMU 
      
            ! The weighted value of the new computation is added
	     ! to the result of the previous iteration
	     Q3_RES(J)=Q3_RES_TMP(J)+AIK*Q3_TMP(J)
	     U3_RES(J)=U3_RES_TMP(J)+AIK*U3_TMP(J)
	     I3_RES(J)=I3_RES_TMP(J)+AIK*I3_TMP(J)
	  
	     !Re-initialisation of the arrays in case the file is finished
	     Q3_TMP(J)=0.
	     U3_TMP(J)=0.
	     I3_TMP(J)=0.
	     Q3_RES_TMP(J)=0.
	     U3_RES_TMP(J)=0.
	     I3_RES_TMP(J)=0.
	  
      ENDDO

C*    Writing of the results

      WRITE(18,err=1302) (Q3_RES(J),J=-NBMU,NBMU),
     & 		    (U3_RES(J),J=-NBMU,NBMU),
     &		           (I3_RES(J),J=-NBMU,NBMU)
        
      IS=IS+1
      GOTO 10
      
 9    CLOSE(16) 
      CLOSE(17)
      CLOSE(18)
      
C --- Deletion of the OS result file and replacement by the new one

C*    Deletion of the old file

      OPEN(UNIT=17,FILE=FICOS,FORM='UNFORMATTED') 
      CLOSE(17,STATUS='DELETE')    

C*    Renaming of the result file

      COMMAND='mv '//FICOS_AGG_TMP(1:LFICOSAGGTMP)//
     &	      ' '//FICOS(1:LFICOS)
      LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') - 1
      IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1001
      CALL SYSTEM(COMMAND)
      

C*    Test of the file existence
      INQUIRE(FILE=FICOS(1:LFICOS),EXIST=EX)
      IF (EX.EQV..FALSE.) GOTO 1002
      
C************************************************************************ 
C* Atmospheric transmission
C* Simple aggregation of content with weighting of the transmission
C************************************************************************ 
      TDIFMUS=TDIFMUS+AIK*TDIFMUS_TMP

      EMOINS=EMOINS+AIK*EMOINS_TMP
      EPLUS=EPLUS+AIK*EPLUS_TMP
      
      DO J=-NBMU,NBMU 
      	TDIFMUG(J)=TDIFMUG(J)+AIK*TDIFMUG_TMP(J)
      ENDDO
      
C***************************************************************************************** 
C* Optical thicknesses
C* Conversion to transmissions (with mus=1), aggregation, then back to optical thickness
C***************************************************************************************** 
  
C*    Truncated total optical thickness
      IF (TTOT_TRONC.NE.0) THEN
      	   TRANS = AIK*DEXP(-TTOT_TRONC_TMP)+DEXP(-TTOT_TRONC)
      ELSE
          TRANS = AIK*DEXP(-TTOT_TRONC_TMP)
      ENDIF	
      TTOT_TRONC = -DLOG(TRANS)
      
C*    True total optical thickness
      IF (TTOT_VRAI.NE.0) THEN
      	   TRANS = AIK*DEXP(-TTOT_VRAI_TMP)+DEXP(-TTOT_VRAI)
      ELSE
          TRANS = AIK*DEXP(-TTOT_VRAI_TMP)
      ENDIF	
      TTOT_VRAI = -DLOG(TRANS)

C*    Optical thickness at the output level for the upward field
      IF (TAUOUT.NE.0) THEN
      	   TRANS = AIK*DEXP(-TAUOUT_TMP)+DEXP(-TAUOUT)
      ELSE
          TRANS = AIK*DEXP(-TAUOUT_TMP)
      ENDIF	
      TAUOUT = -DLOG(TRANS)

C* Nominal end
C-------------------
      GOTO 9999
      
C* Error case
C-------------------

 1001 WRITE(6,*) '  SOS_AGGREGATE : ERROR_1001 : '
      WRITE(6,*) '  => File copy command maybe too long:' 
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      IER=-1      
      GOTO 9999                  	

 1002 WRITE(6,*) '  SOS_AGGREGATE : ERROR_1002 : '
      WRITE(6,*) '  => Ouput OS file has note been created: ',FICOS 
      IER=-1      
      GOTO 9999                  	

 1101 WRITE(6,*) '  SOS_AGGREGATE : ERROR_1101 : '
      WRITE(6,*) '  => Error opening file: ',FICOS_TMP 
      IER=-1      
      GOTO 9999                  	

 1102 WRITE(6,*) '  SOS_AGGREGATE : ERROR_1102 : '
      WRITE(6,*) '  => Error reading file: ',FICOS_TMP 
      IER=-1      
      GOTO 9999                  	

 1201 WRITE(6,*) '  SOS_AGGREGATE : ERROR_1201 : '
      WRITE(6,*) '  => Error opening file: ',FICOS 
      IER=-1      
      GOTO 9999                  	

 1202 WRITE(6,*) '  SOS_AGGREGATE : ERROR_1202 : '
      WRITE(6,*) '  => Error reading file: ',FICOS 
      IER=-1      
      GOTO 9999                  	

 1301 WRITE(6,*) '  SOS_AGGREGATE : ERROR_1301 : '
      WRITE(6,*) '  => Error opening file: ',FICOS_AGG_TMP 
      IER=-1      
      GOTO 9999                  	


 1302 WRITE(6,*) '  SOS_AGGREGATE : ERROR_1302 : '
      WRITE(6,*) '  => Error writing file: ',FICOS_AGG_TMP 
      IER=-1      
      GOTO 9999                  	

C* Exit of the function

 9999 RETURN
 
      END	!END OF PROGRAM SOS_AGGREGATE
