! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! The ODE Function of Chemical Model File
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
! File                 : racm_Function.f90
! Time                 : Wed Jun 14 17:35:28 2017
! Working directory    : /data/ksetigui/setigui/WRFnew/chem/KPP/mechanisms/racm
! Equation file        : racm.kpp
! Output root filename : racm
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE racm_Function

  USE racm_Parameters
  IMPLICIT NONE

! A - Rate for each equation
  REAL(kind=dp) :: A(NREACT)

CONTAINS



END MODULE racm_Function

