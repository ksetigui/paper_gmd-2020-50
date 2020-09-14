! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Parameter Module File
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
! File                 : nmhc9_Parameters.f90
! Time                 : Mon Feb 27 11:08:45 2017
! Working directory    : /data/ksetigui/setigui/WRFV3_171015/chem/KPP/mechanisms/nmhc9
! Equation file        : nmhc9.kpp
! Output root filename : nmhc9
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE nmhc9_Parameters

  USE nmhc9_Precision
  PUBLIC
  SAVE


! NSPEC - Number of chemical species
  INTEGER, PARAMETER :: NSPEC = 59 
! NVAR - Number of Variable species
  INTEGER, PARAMETER :: NVAR = 57 
! NVARACT - Number of Active species
  INTEGER, PARAMETER :: NVARACT = 57 
! NFIX - Number of Fixed species
  INTEGER, PARAMETER :: NFIX = 2 
! NREACT - Number of reactions
  INTEGER, PARAMETER :: NREACT = 145 
! NVARST - Starting of variables in conc. vect.
  INTEGER, PARAMETER :: NVARST = 1 
! NFIXST - Starting of fixed in conc. vect.
  INTEGER, PARAMETER :: NFIXST = 58 
! NONZERO - Number of nonzero entries in Jacobian
  INTEGER, PARAMETER :: NONZERO = 527 
! LU_NONZERO - Number of nonzero entries in LU factoriz. of Jacobian
  INTEGER, PARAMETER :: LU_NONZERO = 580 
! CNVAR - (NVAR+1) Number of elements in compressed row format
  INTEGER, PARAMETER :: CNVAR = 58 
! NLOOKAT - Number of species to look at
  INTEGER, PARAMETER :: NLOOKAT = 0 
! NMONITOR - Number of species to monitor
  INTEGER, PARAMETER :: NMONITOR = 0 
! NMASS - Number of atoms to check mass balance
  INTEGER, PARAMETER :: NMASS = 1 
! PI - Value of pi
  REAL(kind=dp), PARAMETER :: PI = 3.14159265358979 

! Index declaration for variable species in C and VAR
!   VAR(ind_spc) = C(ind_spc)

  INTEGER, PARAMETER :: ind_C2H6 = 1 
  INTEGER, PARAMETER :: ind_C4H10 = 2 
  INTEGER, PARAMETER :: ind_MeCOCO = 3 
  INTEGER, PARAMETER :: ind_C3H8 = 4 
  INTEGER, PARAMETER :: ind_N2O5 = 5 
  INTEGER, PARAMETER :: ind_MeO2NO2 = 6 
  INTEGER, PARAMETER :: ind_NACA = 7 
  INTEGER, PARAMETER :: ind_H2O2 = 8 
  INTEGER, PARAMETER :: ind_PAA = 9 
  INTEGER, PARAMETER :: ind_CH3COOH = 10 
  INTEGER, PARAMETER :: ind_ISOOH = 11 
  INTEGER, PARAMETER :: ind_MeOOH = 12 
  INTEGER, PARAMETER :: ind_PAN = 13 
  INTEGER, PARAMETER :: ind_MPAN = 14 
  INTEGER, PARAMETER :: ind_HCOOH = 15 
  INTEGER, PARAMETER :: ind_EtOOH = 16 
  INTEGER, PARAMETER :: ind_PrOOH = 17 
  INTEGER, PARAMETER :: ind_PrONO2 = 18 
  INTEGER, PARAMETER :: ind_C2H4 = 19 
  INTEGER, PARAMETER :: ind_MEKOOH = 20 
  INTEGER, PARAMETER :: ind_O1D = 21 
  INTEGER, PARAMETER :: ind_HNO4 = 22 
  INTEGER, PARAMETER :: ind_ACETP = 23 
  INTEGER, PARAMETER :: ind_CH4 = 24 
  INTEGER, PARAMETER :: ind_HNO3 = 25 
  INTEGER, PARAMETER :: ind_C3H6OOH = 26 
  INTEGER, PARAMETER :: ind_C4H9OOH = 27 
  INTEGER, PARAMETER :: ind_ACET = 28 
  INTEGER, PARAMETER :: ind_MVKOOH = 29 
  INTEGER, PARAMETER :: ind_MEK = 30 
  INTEGER, PARAMETER :: ind_CO = 31 
  INTEGER, PARAMETER :: ind_MeOH = 32 
  INTEGER, PARAMETER :: ind_ISON = 33 
  INTEGER, PARAMETER :: ind_PrO2 = 34 
  INTEGER, PARAMETER :: ind_MGLO = 35 
  INTEGER, PARAMETER :: ind_ACETO2 = 36 
  INTEGER, PARAMETER :: ind_C3H6 = 37 
  INTEGER, PARAMETER :: ind_ISOP = 38 
  INTEGER, PARAMETER :: ind_ACETOL = 39 
  INTEGER, PARAMETER :: ind_C3H6O2 = 40 
  INTEGER, PARAMETER :: ind_MVK = 41 
  INTEGER, PARAMETER :: ind_MEKO2 = 42 
  INTEGER, PARAMETER :: ind_C4H9O2 = 43 
  INTEGER, PARAMETER :: ind_ONIT = 44 
  INTEGER, PARAMETER :: ind_HCHO = 45 
  INTEGER, PARAMETER :: ind_ALD = 46 
  INTEGER, PARAMETER :: ind_EtO2 = 47 
  INTEGER, PARAMETER :: ind_ISO2 = 48 
  INTEGER, PARAMETER :: ind_MVKO2 = 49 
  INTEGER, PARAMETER :: ind_MeO2 = 50 
  INTEGER, PARAMETER :: ind_PA = 51 
  INTEGER, PARAMETER :: ind_NO3 = 52 
  INTEGER, PARAMETER :: ind_O3 = 53 
  INTEGER, PARAMETER :: ind_NO2 = 54 
  INTEGER, PARAMETER :: ind_OH = 55 
  INTEGER, PARAMETER :: ind_NO = 56 
  INTEGER, PARAMETER :: ind_HO2 = 57 

! Index declaration for fixed species in C
!   C(ind_spc)

  INTEGER, PARAMETER :: ind_H2O = 58 
  INTEGER, PARAMETER :: ind_M = 59 

! Index declaration for fixed species in FIX
!    FIX(indf_spc) = C(ind_spc) = C(NVAR+indf_spc)

  INTEGER, PARAMETER :: indf_H2O = 1 
  INTEGER, PARAMETER :: indf_M = 2 

END MODULE nmhc9_Parameters

