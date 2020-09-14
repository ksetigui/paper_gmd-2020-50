! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Main Program File
! 
! Generated by KPP-2.1 symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : saprc99_mosaic_4bin_vbs2_Main.f90
! Time                 : Wed Jun 14 17:35:30 2017
! Working directory    : /data/ksetigui/setigui/WRFnew/chem/KPP/mechanisms/saprc99_mosaic_4bin_vbs2
! Equation file        : saprc99_mosaic_4bin_vbs2.kpp
! Output root filename : saprc99_mosaic_4bin_vbs2
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! MAIN - Main program - driver routine
!   Arguments :
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM saprc99_mosaic_4bin_vbs2_Driver

  USE saprc99_mosaic_4bin_vbs2_Model
  USE saprc99_mosaic_4bin_vbs2_Initialize, ONLY: Initialize

      REAL(kind=dp) :: T, DVAL(NSPEC)
      REAL(kind=dp) :: RSTATE(20)
      INTEGER :: i
  
!~~~> TIME VARIABLES 

      STEPMIN = 0.0d0
      STEPMAX = 0.0d0

      DO i=1,NVAR
        RTOL(i) = 1.0d-4
        ATOL(i) = 1.0d-3
      END DO
     
      CALL Initialize()
      CALL InitSaveData()

!~~~> Time loop
      T = TSTART
kron: DO WHILE (T < TEND)

        CALL GetMass( C, DVAL )
        WRITE(6,991) (T-TSTART)/(TEND-TSTART)*100, T,       &
                   ( TRIM(SPC_NAMES(MONITOR(i))),           &
                     C(MONITOR(i))/CFACTOR, i=1,NMONITOR ), &
                   ( TRIM(SMASS(i)), DVAL(i)/CFACTOR, i=1,NMASS )
        TIME = T
        CALL SaveData()
        CALL Update_SUN() 
        CALL Update_RCONST()

        CALL INTEGRATE( TIN = T, TOUT = T+DT, RSTATUS_U = RSTATE, &
        ICNTRL_U = (/ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) )
        T = RSTATE(1)

      END DO kron
!~~~> End Time loop

      CALL GetMass( C, DVAL )
      WRITE(6,991) (T-TSTART)/(TEND-TSTART)*100, T,     &
               ( TRIM(SPC_NAMES(MONITOR(i))),           &
                 C(MONITOR(i))/CFACTOR, i=1,NMONITOR ), &
               ( TRIM(SMASS(i)), DVAL(i)/CFACTOR, i=1,NMASS )
      TIME = T
      CALL SaveData()
      CALL CloseSaveData()

990   FORMAT('Done[%]. Time ',20(4X,A12))
991   FORMAT(F6.1,'%. T=',E9.3,2X,200(A,'=',E11.4,'; '))

END PROGRAM saprc99_mosaic_4bin_vbs2_Driver

! End of MAIN function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


