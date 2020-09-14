! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Utility Data Module File
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
! File                 : radm2sorg_Monitor.f90
! Time                 : Mon Feb 27 11:08:59 2017
! Working directory    : /data/ksetigui/setigui/WRFV3_171015/chem/KPP/mechanisms/radm2sorg
! Equation file        : radm2sorg.kpp
! Output root filename : radm2sorg
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE radm2sorg_Monitor


  CHARACTER(LEN=12), PARAMETER, DIMENSION(61) :: SPC_NAMES = (/ &
     'SULF        ','ORA1        ','ORA2        ', &
     'CO2         ','SO2         ','ETH         ', &
     'O1D         ','HC5         ','HC8         ', &
     'TOL         ','XYL         ','TPAN        ', &
     'HONO        ','H2O2        ','N2O5        ', &
     'HC3         ','CH4         ','PAA         ', &
     'O3P         ','HNO4        ','OP1         ', &
     'CSL         ','PAN         ','OL2         ', &
     'HNO3        ','CO          ','ISO         ', &
     'OLT         ','OLI         ','DCB         ', &
     'GLY         ','XNO2        ','KET         ', &
     'MGLY        ','TOLP        ','XYLP        ', &
     'OLTP        ','OLN         ','XO2         ', &
     'OL2P        ','HC5P        ','OP2         ', &
     'HCHO        ','HC8P        ','TCO3        ', &
     'O3          ','ONIT        ','ALD         ', &
     'OLIP        ','KETP        ','HO2         ', &
     'MO2         ','OH          ','NO3         ', &
     'ACO3        ','HC3P        ','ETHP        ', &
     'NO          ','NO2         ','H2O         ', &
     'M           ' /)

  INTEGER, DIMENSION(1) :: LOOKAT
  INTEGER, DIMENSION(1) :: MONITOR
  CHARACTER(LEN=12), DIMENSION(1) :: SMASS
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_0 = (/ &
     '        NO2 --> O3P + NO                                                                            ', &
     '         O3 --> O1D                                                                                 ', &
     '         O3 --> O3P                                                                                 ', &
     '       HONO --> OH + NO                                                                             ', &
     '       HNO3 --> OH + NO2                                                                            ', &
     '       HNO4 --> 0.65 HO2 + 0.35 OH + 0.35 NO3 + 0.65 NO2                                            ', &
     '        NO3 --> NO                                                                                  ', &
     '        NO3 --> O3P + NO2                                                                           ', &
     '       H2O2 --> 2 OH                                                                                ', &
     '       HCHO --> CO                                                                                  ', &
     '       HCHO --> CO + 2 HO2                                                                          ', &
     '        ALD --> CO + HO2 + MO2                                                                      ', &
     '        OP1 --> HCHO + HO2 + OH                                                                     ', &
     '        OP2 --> ALD + HO2 + OH                                                                      ', &
     '        PAA --> MO2 + OH                                                                            ', &
     '        KET --> ACO3 + ETHP                                                                         ', &
     '        GLY --> 1.87 CO + 0.13 HCHO                                                                 ', &
     '        GLY --> 1.55 CO + 0.45 HCHO + 0.8 HO2                                                       ', &
     '       MGLY --> CO + HO2 + ACO3                                                                     ', &
     '        DCB --> TCO3 + HO2                                                                          ', &
     '       ONIT --> 0.8 KET + 0.2 ALD + HO2 + NO2                                                       ', &
     '    O3P + M --> O3                                                                                  ', &
     '  O3P + NO2 --> NO                                                                                  ', &
     '    O1D + M --> O3P                                                                                 ', &
     '  O1D + H2O --> 2 OH                                                                                ', &
     '    O3 + NO --> NO2                                                                                 ', &
     '    O3 + OH --> HO2                                                                                 ', &
     '   O3 + HO2 --> OH                                                                                  ', &
     '   HO2 + NO --> OH + NO2                                                                            ', &
     '  HO2 + NO2 --> HNO4                                                                                ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_1 = (/ &
     '       HNO4 --> HO2 + NO2                                                                           ', &
     '      2 HO2 --> H2O2                                                                                ', &
     '2 HO2 + H2O --> H2O2                                                                                ', &
     '  H2O2 + OH --> HO2 + H2O                                                                           ', &
     '    OH + NO --> HONO                                                                                ', &
     '   2 NO + M --> 2 NO2                                                                               ', &
     '   O3 + NO2 --> NO3                                                                                 ', &
     '   NO3 + NO --> 2 NO2                                                                               ', &
     '  NO3 + NO2 --> NO + NO2                                                                            ', &
     '  HO2 + NO3 --> HNO3                                                                                ', &
     '  NO3 + NO2 --> N2O5                                                                                ', &
     '       N2O5 --> NO3 + NO2                                                                           ', &
     '       N2O5 --> 2 HNO3                                                                              ', &
     '   OH + NO2 --> HNO3                                                                                ', &
     '  HNO3 + OH --> NO3 + H2O                                                                           ', &
     '  HNO4 + OH --> NO2 + H2O                                                                           ', &
     '   HO2 + OH --> H2O                                                                                 ', &
     '   SO2 + OH --> SULF + HO2                                                                          ', &
     '    CO + OH --> CO2 + HO2                                                                           ', &
     '   CH4 + OH --> MO2 + H2O                                                                           ', &
     '   ETH + OH --> ETHP + H2O                                                                          ', &
     '   HC3 + OH --> 0.025 KET + 0.009 HCHO + 0.075 ALD + 0.17 HO2 + 0.83 HC3P + H2O                     ', &
     '   HC5 + OH --> 0.25 XO2 + HC5P + H2O                                                               ', &
     '   HC8 + OH --> 0.75 XO2 + HC8P + H2O                                                               ', &
     '   OL2 + OH --> OL2P                                                                                ', &
     '   OLT + OH --> OLTP                                                                                ', &
     '   OLI + OH --> OLIP                                                                                ', &
     '   TOL + OH --> 0.25 CSL + 0.75 TOLP + 0.25 HO2                                                     ', &
     '   XYL + OH --> 0.17 CSL + 0.83 XYLP + 0.17 HO2                                                     ', &
     '   CSL + OH --> 0.9 XO2 + 0.9 TCO3 + 0.1 HO2 - 0.9 OH                                               ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_2 = (/ &
     '  HCHO + OH --> CO + HO2 + H2O                                                                      ', &
     '   ALD + OH --> ACO3 + H2O                                                                          ', &
     '   KET + OH --> KETP + H2O                                                                          ', &
     '   GLY + OH --> 2 CO + HO2 + H2O                                                                    ', &
     '  MGLY + OH --> CO + ACO3 + H2O                                                                     ', &
     '   DCB + OH --> TCO3 + H2O                                                                          ', &
     '   OP1 + OH --> 0.5 HCHO + 0.5 MO2 + 0.5 OH                                                         ', &
     '   OP2 + OH --> 0.5 ALD + 0.5 OH + 0.5 HC3P                                                         ', &
     '   PAA + OH --> ACO3 + H2O                                                                          ', &
     '   PAN + OH --> XO2 + HCHO + NO3                                                                    ', &
     '  ONIT + OH --> HC3P + NO2                                                                          ', &
     '   ISO + OH --> OLTP                                                                                ', &
     ' ACO3 + NO2 --> PAN                                                                                 ', &
     '        PAN --> ACO3 + NO2                                                                          ', &
     ' TCO3 + NO2 --> TPAN                                                                                ', &
     '       TPAN --> TCO3 + NO2                                                                          ', &
     '   MO2 + NO --> HCHO + HO2 + NO2                                                                    ', &
     '  HC3P + NO --> 0.25 KET + 0.09 HCHO + 0.036 ONIT + 0.75 ALD + 0.964 HO2 + 0.964 NO2                ', &
     '  HC5P + NO --> 0.69 KET + 0.08 ONIT + 0.38 ALD + 0.92 HO2 + 0.92 NO2                               ', &
     '  HC8P + NO --> 1.06 KET + 0.04 HCHO + 0.24 ONIT + 0.35 ALD + 0.76 HO2 + 0.76 NO2                   ', &
     '  OL2P + NO --> 1.6 HCHO + 0.2 ALD + HO2 + NO2                                                      ', &
     '  OLTP + NO --> HCHO + ALD + HO2 + NO2                                                              ', &
     '  OLIP + NO --> 0.1 KET + 0.28 HCHO + 1.45 ALD + HO2 + NO2                                          ', &
     '  ACO3 + NO --> MO2 + NO2                                                                           ', &
     '  TCO3 + NO --> 0.95 CO + 0.89 GLY + 0.11 MGLY + 2 XO2 + 0.92 HO2 + 0.05 ACO3 + NO2                 ', &
     '  TOLP + NO --> 0.7 DCB + 0.16 GLY + 0.17 MGLY + HO2 + NO2                                          ', &
     '  XYLP + NO --> 0.806 DCB + 0.45 MGLY + HO2 + NO2                                                   ', &
     '  ETHP + NO --> ALD + HO2 + NO2                                                                     ', &
     '  KETP + NO --> MGLY + HO2 + NO2                                                                    ', &
     '   OLN + NO --> HCHO + ALD + 2 NO2                                                                  ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_3 = (/ &
     ' HCHO + NO3 --> HNO3 + CO + HO2                                                                     ', &
     '  ALD + NO3 --> HNO3 + ACO3                                                                         ', &
     '  GLY + NO3 --> HNO3 + 2 CO + HO2                                                                   ', &
     ' MGLY + NO3 --> HNO3 + CO + ACO3                                                                    ', &
     '  DCB + NO3 --> HNO3 + TCO3                                                                         ', &
     '  CSL + NO3 --> 0.5 CSL + HNO3 + XNO2                                                               ', &
     '  OL2 + NO3 --> OLN                                                                                 ', &
     '  OLT + NO3 --> OLN                                                                                 ', &
     '  OLI + NO3 --> OLN                                                                                 ', &
     '  ISO + NO3 --> OLN                                                                                 ', &
     '   OL2 + O3 --> 0.4 ORA1 + 0.42 CO + HCHO + 0.12 HO2                                                ', &
     '   OLT + O3 --> 0.2 ORA1 + 0.2 ORA2 + 0.06 CH4 + 0.33 CO + 0.53 HCHO + 0.5 ALD + 0.23 HO2 + 0.22 MO2', &
     '   OLI + O3 --> 0.06 ORA1 + 0.29 ORA2 + 0.09 CH4 + 0.23 CO + 0.1 KET + 0.18 HCHO + 0.72 ALD + 0.26 H', &
     '   ISO + O3 --> 0.2 ORA1 + 0.2 ORA2 + 0.33 CO + 0.53 HCHO + 0.5 ALD + 0.23 HO2 + 0.22 MO2 + 0.1 OH  ', &
     '  HO2 + MO2 --> OP1                                                                                 ', &
     ' HO2 + ETHP --> OP2                                                                                 ', &
     ' HO2 + HC3P --> OP2                                                                                 ', &
     ' HC5P + HO2 --> OP2                                                                                 ', &
     ' HC8P + HO2 --> OP2                                                                                 ', &
     ' OL2P + HO2 --> OP2                                                                                 ', &
     ' OLTP + HO2 --> OP2                                                                                 ', &
     ' OLIP + HO2 --> OP2                                                                                 ', &
     ' KETP + HO2 --> OP2                                                                                 ', &
     ' HO2 + ACO3 --> PAA                                                                                 ', &
     ' TOLP + HO2 --> OP2                                                                                 ', &
     ' XYLP + HO2 --> OP2                                                                                 ', &
     ' TCO3 + HO2 --> OP2                                                                                 ', &
     '  OLN + HO2 --> ONIT                                                                                ', &
     '      2 MO2 --> 1.5 HCHO + HO2                                                                      ', &
     ' MO2 + ETHP --> 0.75 HCHO + 0.75 ALD + HO2                                                          ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_4 = (/ &
     ' MO2 + HC3P --> 0.6 KET + 0.75 HCHO + 0.15 ALD + HO2                                                ', &
     ' HC5P + MO2 --> 0.75 KET + 0.77 HCHO + 0.41 ALD + HO2                                               ', &
     ' HC8P + MO2 --> 1.39 KET + 0.8 HCHO + 0.46 ALD + HO2                                                ', &
     ' OL2P + MO2 --> 1.55 HCHO + 0.35 ALD + HO2                                                          ', &
     ' OLTP + MO2 --> 1.25 HCHO + 0.75 ALD + HO2                                                          ', &
     ' OLIP + MO2 --> 0.55 KET + 0.89 HCHO + 0.725 ALD + HO2                                              ', &
     ' KETP + MO2 --> 0.75 MGLY + 0.75 HCHO + HO2                                                         ', &
     ' MO2 + ACO3 --> 0.5 ORA2 + HCHO + 0.5 HO2 + 0.5 MO2                                                 ', &
     ' TOLP + MO2 --> 0.7 DCB + 0.16 GLY + 0.17 MGLY + HCHO + 2 HO2                                       ', &
     ' XYLP + MO2 --> 0.806 DCB + 0.45 MGLY + HCHO + 2 HO2                                                ', &
     ' TCO3 + MO2 --> 0.5 ORA2 + 0.475 CO + 0.445 GLY + 0.055 MGLY + XO2 + 0.5 HCHO + 0.46 HO2 + 0.025 ACO', &
     'ACO3 + ETHP --> 0.5 ORA2 + ALD + 0.5 HO2 + 0.5 MO2                                                  ', &
     'ACO3 + HC3P --> 0.5 ORA2 + 0.8 KET + 0.2 ALD + 0.5 HO2 + 0.5 MO2                                    ', &
     'HC5P + ACO3 --> 0.5 ORA2 + 0.86 KET + 0.14 ALD + 0.5 HO2 + 0.5 MO2                                  ', &
     'HC8P + ACO3 --> 0.5 ORA2 + 0.9 KET + 0.1 ALD + 0.5 HO2 + 0.5 MO2                                    ', &
     'OL2P + ACO3 --> 0.5 ORA2 + 0.8 HCHO + 0.6 ALD + 0.5 HO2 + 0.5 MO2                                   ', &
     'OLTP + ACO3 --> 0.5 ORA2 + 0.5 HCHO + ALD + 0.5 HO2 + 0.5 MO2                                       ', &
     'OLIP + ACO3 --> 0.5 ORA2 + 0.55 KET + 0.14 HCHO + 0.725 ALD + 0.5 HO2 + 0.5 MO2                     ', &
     'KETP + ACO3 --> 0.5 ORA2 + MGLY + 0.5 HO2 + 0.5 MO2                                                 ', &
     '     2 ACO3 --> 2 MO2                                                                               ', &
     'TOLP + ACO3 --> DCB + 0.2 GLY + 0.8 MGLY + HO2 + MO2                                                ', &
     'XYLP + ACO3 --> DCB + MGLY + HO2 + MO2                                                              ', &
     'TCO3 + ACO3 --> 0.95 CO + 0.89 GLY + 0.11 MGLY + 2 XO2 + 0.92 HO2 + MO2 + 0.05 ACO3                 ', &
     '  XO2 + HO2 --> OP2                                                                                 ', &
     '  XO2 + MO2 --> HCHO + HO2                                                                          ', &
     ' XO2 + ACO3 --> MO2                                                                                 ', &
     '      2 XO2 --> H2O                                                                                 ', &
     '   XO2 + NO --> NO2                                                                                 ', &
     ' XNO2 + NO2 --> ONIT                                                                                ', &
     ' XNO2 + HO2 --> OP2                                                                                 ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(6) :: EQN_NAMES_5 = (/ &
     ' XNO2 + MO2 --> HCHO + HO2                                                                          ', &
     'XNO2 + ACO3 --> MO2                                                                                 ', &
     '     2 XNO2 --> H2O                                                                                 ', &
     '  OLN + MO2 --> 1.75 HCHO + ALD + 0.5 HO2 + NO2                                                     ', &
     ' OLN + ACO3 --> 0.5 ORA2 + HCHO + ALD + 0.5 MO2 + NO2                                               ', &
     '      2 OLN --> 2 HCHO + 2 ALD + 2 NO2                                                              ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(156) :: EQN_NAMES = (/&
    EQN_NAMES_0, EQN_NAMES_1, EQN_NAMES_2, EQN_NAMES_3, EQN_NAMES_4, &
    EQN_NAMES_5 /)

! INLINED global variables

! End INLINED global variables


END MODULE radm2sorg_Monitor
