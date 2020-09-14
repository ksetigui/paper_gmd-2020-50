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
! File                 : mozart_mosaic_4bin_vbs0_Monitor.f90
! Time                 : Mon Feb 27 11:08:41 2017
! Working directory    : /data/ksetigui/setigui/WRFV3_171015/chem/KPP/mechanisms/mozart_mosaic_4bin_vbs0
! Equation file        : mozart_mosaic_4bin_vbs0.kpp
! Output root filename : mozart_mosaic_4bin_vbs0
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE mozart_mosaic_4bin_vbs0_Monitor


  CHARACTER(LEN=12), PARAMETER, DIMENSION(90) :: SPC_NAMES_0 = (/ &
     'SO4         ','NUME        ','DEN         ', &
     'BIOG1_c     ','BIOG1_o     ','SMPA        ', &
     'VOCA        ','SMPBB       ','VOCBB       ', &
     'NH3         ','C2H6        ','C3H8        ', &
     'BIGENE      ','BIGALK      ','N2O         ', &
     'SO2         ','TOLUENE     ','CRESOL      ', &
     'H2O2        ','EO          ','N2O5        ', &
     'XOH         ','XOOH        ','DMS         ', &
     'C2H5OH      ','MACROOH     ','H2          ', &
     'C2H5OOH     ','C3H7OOH     ','ROOH        ', &
     'ENEO2       ','MEKOOH      ','CH3OOH      ', &
     'HO2NO2      ','HYDRALD     ','C2H4        ', &
     'ONIT        ','TOLOOH      ','EO2         ', &
     'CH3COOOH    ','MEK         ','CH3COOH     ', &
     'GLYOXAL     ','TERPOOH     ','HNO3        ', &
     'POOH        ','PAN         ','CH4         ', &
     'ISOPOOH     ','MPAN        ','CH3OH       ', &
     'TOLO2       ','BIGALD      ','C10H16      ', &
     'O1D_CB4     ','ALKOOH      ','MEKO2       ', &
     'CH3COCH3    ','PO2         ','O           ', &
     'ISOPNO3     ','GLYALD      ','ALKO2       ', &
     'C2H5O2      ','CO          ','C3H7O2      ', &
     'ISOP        ','ONITR       ','C3H6        ', &
     'CH3CHO      ','TERPO2      ','RO2         ', &
     'HYAC        ','MACR        ','CH3COCHO    ', &
     'XO2         ','MVK         ','MACRO2      ', &
     'CH2O        ','ISOPO2      ','NO3         ', &
     'CH3O2       ','O3          ','OH          ', &
     'MCO3        ','HO2         ','NO          ', &
     'NO2         ','CH3CO3      ','H2O         ' /)
  CHARACTER(LEN=12), PARAMETER, DIMENSION(1) :: SPC_NAMES_1 = (/ &
     'M           ' /)
  CHARACTER(LEN=12), PARAMETER, DIMENSION(91) :: SPC_NAMES = (/&
    SPC_NAMES_0, SPC_NAMES_1 /)

  INTEGER, DIMENSION(1) :: LOOKAT
  INTEGER, DIMENSION(1) :: MONITOR
  CHARACTER(LEN=12), DIMENSION(1) :: SMASS
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_0 = (/ &
     '              M --> 2 O                                                                             ', &
     '             O3 --> O1D_CB4                                                                         ', &
     '             O3 --> O                                                                               ', &
     '            N2O --> O1D_CB4                                                                         ', &
     '            NO2 --> O + NO                                                                          ', &
     '           N2O5 --> NO3 + NO2                                                                       ', &
     '           HNO3 --> OH + NO2                                                                        ', &
     '            NO3 --> 0.89 O3 + 0.11 NO + 0.89 NO2                                                    ', &
     '         HO2NO2 --> 0.33 NO3 + 0.33 OH + 0.66 HO2 + 0.66 NO2                                        ', &
     '         CH3OOH --> CH2O + OH + HO2                                                                 ', &
     '           CH2O --> CO + 2 HO2                                                                      ', &
     '           CH2O --> H2 + CO                                                                         ', &
     '           H2O2 --> 2 OH                                                                            ', &
     '         CH3CHO --> CO + CH3O2 + HO2                                                                ', &
     '           POOH --> CH3CHO + CH2O + OH + HO2                                                        ', &
     '       CH3COOOH --> CH3O2 + OH                                                                      ', &
     '            PAN --> 0.4 NO3 + 0.4 CH3O2 + 0.6 NO2 + 0.6 CH3CO3                                      ', &
     '           MPAN --> MCO3 + NO2                                                                      ', &
     '           MACR --> 0.67 CO + 0.67 CH2O + 0.33 OH + 0.33 MCO3 + 0.67 HO2 + 0.67 CH3CO3              ', &
     '            MVK --> 0.7 CO + 0.7 C3H6 + 0.3 CH3O2 + 0.3 CH3CO3                                      ', &
     '        C2H5OOH --> CH3CHO + OH + HO2                                                               ', &
     '        C3H7OOH --> 0.82 CH3COCH3 + OH + HO2                                                        ', &
     '           ROOH --> CH2O + OH + CH3CO3                                                              ', &
     '       CH3COCH3 --> CH3O2 + CH3CO3                                                                  ', &
     '       CH3COCHO --> CO + HO2 + CH3CO3                                                               ', &
     '           XOOH --> OH                                                                              ', &
     '          ONITR --> CO + CH2O + HO2 + NO2                                                           ', &
     '        ISOPOOH --> 0.288 MACR + 0.402 MVK + 0.69 CH2O + HO2                                        ', &
     '           HYAC --> CH2O + HO2 + CH3CO3                                                             ', &
     '         GLYALD --> CO + CH2O + 2 HO2                                                               ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_1 = (/ &
     '            MEK --> C2H5O2 + CH3CO3                                                                 ', &
     '         BIGALD --> 0.13 GLYOXAL + 0.45 CO + 0.18 CH3COCHO + 0.56 HO2 + 0.13 CH3CO3                 ', &
     '        GLYOXAL --> 2 CO + 2 HO2                                                                    ', &
     '         ALKOOH --> 0.8 MEK + 0.25 CH3COCH3 + 0.4 CH3CHO + 0.1 CH2O + OH + 0.9 HO2                  ', &
     '         MEKOOH --> CH3CHO + OH + CH3CO3                                                            ', &
     '         TOLOOH --> 0.45 GLYOXAL + 0.9 BIGALD + 0.45 CH3COCHO + OH                                  ', &
     '        TERPOOH --> 0.1 CH3COCH3 + MACR + MVK + OH + HO2                                            ', &
     '          O + M --> O3                                                                              ', &
     '         O + O3 --> M                                                                               ', &
     '    O1D_CB4 + M --> O                                                                               ', &
     '  O1D_CB4 + H2O --> 2 OH                                                                            ', &
     '   H2 + O1D_CB4 --> OH + HO2                                                                        ', &
     '        H2 + OH --> HO2 + H2O                                                                       ', &
     '         O + OH --> HO2                                                                             ', &
     '        O + HO2 --> OH                                                                              ', &
     '        O3 + OH --> HO2                                                                             ', &
     '       O3 + HO2 --> OH                                                                              ', &
     '    2 HO2 + H2O --> H2O2                                                                            ', &
     '      H2O2 + OH --> HO2 + H2O                                                                       ', &
     '       OH + HO2 --> H2O                                                                             ', &
     '           2 OH --> O + H2O                                                                         ', &
     '           2 OH --> H2O2                                                                            ', &
     '  N2O + O1D_CB4 --> 2 NO                                                                            ', &
     '  N2O + O1D_CB4 --> M                                                                               ', &
     '       HO2 + NO --> OH + NO2                                                                        ', &
     '        O3 + NO --> NO2                                                                             ', &
     '        O + NO2 --> NO                                                                              ', &
     '       O3 + NO2 --> NO3                                                                             ', &
     '      NO3 + HO2 --> OH + NO2                                                                        ', &
     '      NO3 + NO2 --> N2O5                                                                            ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_2 = (/ &
     '           N2O5 --> NO3 + NO2                                                                       ', &
     '       OH + NO2 --> HNO3                                                                            ', &
     '      HNO3 + OH --> NO3 + H2O                                                                       ', &
     '       NO3 + NO --> 2 NO2                                                                           ', &
     '      HO2 + NO2 --> HO2NO2                                                                          ', &
     '    HO2NO2 + OH --> NO2 + H2O                                                                       ', &
     '         HO2NO2 --> HO2 + NO2                                                                       ', &
     '       N2O5 + M --> 2 HNO3                                                                          ', &
     '            NO3 --> HNO3                                                                            ', &
     '            NO2 --> 0.5 HNO3 + 0.5 OH + 0.5 NO                                                      ', &
     '       CH4 + OH --> CH3O2 + H2O                                                                     ', &
     '  CH4 + O1D_CB4 --> 0.05 H2 + 0.25 CH2O + 0.75 CH3O2 + 0.75 OH + 0.4 HO2                            ', &
     '     CH3O2 + NO --> CH2O + HO2 + NO2                                                                ', &
     '        2 CH3O2 --> 2 CH2O + 2 HO2                                                                  ', &
     '        2 CH3O2 --> CH3OH + CH2O                                                                    ', &
     '    CH3O2 + HO2 --> CH3OOH                                                                          ', &
     '    CH3OOH + OH --> 0.3 CH2O + 0.7 CH3O2 + 0.3 OH + H2O                                             ', &
     '     CH2O + NO3 --> HNO3 + CO + HO2                                                                 ', &
     '      CH2O + OH --> CO + HO2 + H2O                                                                  ', &
     '        CO + OH --> HO2                                                                             ', &
     '      C2H4 + OH --> 0.75 EO2 + 0.5 CH2O + 0.25 HO2                                                  ', &
     '      C2H4 + O3 --> 0.25 CH3COOH + 0.5 CO + CH2O + 0.12 OH + 0.12 HO2                               ', &
     '       SO2 + OH --> SO4                                                                             ', &
     '   GLYOXAL + OH --> CO + HO2                                                                        ', &
     '       EO2 + NO --> EO + NO2                                                                        ', &
     '         EO + M --> GLYALD + HO2                                                                    ', &
     '             EO --> 2 CH2O + HO2                                                                    ', &
     '      C2H6 + OH --> C2H5O2                                                                          ', &
     '    C2H5O2 + NO --> CH3CHO + HO2 + NO2                                                              ', &
     '   C2H5O2 + HO2 --> C2H5OOH                                                                         ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_3 = (/ &
     ' C2H5O2 + CH3O2 --> 0.2 C2H5OH + 0.3 CH3OH + 0.8 CH3CHO + 0.7 CH2O + HO2                            ', &
     '   C2H5OOH + OH --> 0.5 C2H5O2 + 0.5 CH3CHO + 0.5 OH                                                ', &
     '      C3H6 + OH --> PO2                                                                             ', &
     '      C3H6 + O3 --> 0.25 CH3COOH + 0.08 CH4 + 0.56 CO + 0.5 CH3CHO + 0.54 CH2O + 0.31 CH3O2 + 0.33 O', &
     '     C3H6 + NO3 --> ONIT                                                                            ', &
     '       PO2 + NO --> CH3CHO + CH2O + HO2 + NO2                                                       ', &
     '      PO2 + HO2 --> POOH                                                                            ', &
     '      POOH + OH --> 0.5 PO2 + 0.5 HYAC + 0.5 OH                                                     ', &
     '    CH3CHO + OH --> CH3CO3                                                                          ', &
     '   CH3CHO + NO3 --> HNO3 + CH3CO3                                                                   ', &
     '    NO + CH3CO3 --> CH3O2 + NO2                                                                     ', &
     '   NO2 + CH3CO3 --> PAN                                                                             ', &
     '   HO2 + CH3CO3 --> 0.75 CH3COOOH + 0.25 CH3COOH + 0.25 O3                                          ', &
     ' CH3O2 + CH3CO3 --> 0.1 CH3COOH + CH2O + 0.9 CH3O2 + 0.9 HO2                                        ', &
     '  CH3COOOH + OH --> 0.5 CH2O + 0.5 CH3CO3                                                           ', &
     '            PAN --> NO2 + CH3CO3                                                                    ', &
     '       2 CH3CO3 --> 2 CH3O2                                                                         ', &
     '      C3H8 + OH --> C3H7O2                                                                          ', &
     '    C3H7O2 + NO --> 0.82 CH3COCH3 + 0.27 CH3CHO + HO2 + NO2                                         ', &
     '   C3H7O2 + HO2 --> C3H7OOH                                                                         ', &
     ' C3H7O2 + CH3O2 --> 0.82 CH3COCH3 + CH2O + HO2                                                      ', &
     '   C3H7OOH + OH --> C3H7O2                                                                          ', &
     '  CH3COCH3 + OH --> RO2                                                                             ', &
     '       RO2 + NO --> NUME + CH2O + NO2 + CH3CO3                                                      ', &
     '      RO2 + HO2 --> DEN + ROOH                                                                      ', &
     '    RO2 + CH3O2 --> 0.5 CH3OH + 0.2 HYAC + 0.5 CH3COCHO + 0.8 CH2O + 0.3 HO2 + 0.3 CH3CO3           ', &
     '      ROOH + OH --> RO2                                                                             ', &
     '    BIGENE + OH --> ENEO2                                                                           ', &
     '     ENEO2 + NO --> 0.5 CH3COCH3 + CH3CHO + 0.5 CH2O + HO2 + NO2                                    ', &
     '    BIGALK + OH --> ALKO2                                                                           ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_4 = (/ &
     '     ALKO2 + NO --> 0.1 ONIT + 0.75 MEK + 0.25 CH3COCH3 + 0.4 CH3CHO + 0.1 CH2O + 0.9 HO2 + 0.9 NO2 ', &
     '    ALKO2 + HO2 --> ALKOOH                                                                          ', &
     '    ALKOOH + OH --> ALKO2                                                                           ', &
     '      ONIT + OH --> CH3COCHO + NO2                                                                  ', &
     '       MEK + OH --> MEKO2                                                                           ', &
     '     MEKO2 + NO --> CH3CHO + NO2 + CH3CO3                                                           ', &
     '    MEKO2 + HO2 --> MEKOOH                                                                          ', &
     '    MEKOOH + OH --> MEKO2                                                                           ', &
     '   TOLUENE + OH --> 0.25 CRESOL + 0.7 TOLO2 + 0.25 HO2                                              ', &
     '    CRESOL + OH --> XOH                                                                             ', &
     '      XOH + NO2 --> 0.7 BIGALD + 0.7 HO2 + 0.7 NO2                                                  ', &
     '     TOLO2 + NO --> 0.45 GLYOXAL + 0.9 BIGALD + 0.45 CH3COCHO + 0.9 HO2 + 0.9 NO2                   ', &
     '    TOLO2 + HO2 --> TOLOOH                                                                          ', &
     '    TOLOOH + OH --> TOLO2                                                                           ', &
     '      ISOP + OH --> ISOPO2                                                                          ', &
     '      ISOP + O3 --> 0.2 CH3COOH + 0.3 CO + 0.07 C3H6 + 0.4 MACR + 0.2 MVK + 0.6 CH2O + 0.1 O3 + 0.27', &
     '    ISOPO2 + NO --> 0.37 HYDRALD + 0.08 ONITR + 0.23 MACR + 0.32 MVK + 0.55 CH2O + HO2 + 0.92 NO2   ', &
     '   ISOPO2 + NO3 --> 0.4 HYDRALD + 0.25 MACR + 0.35 MVK + 0.6 CH2O + HO2 + NO2                       ', &
     '   ISOPO2 + HO2 --> ISOPOOH                                                                         ', &
     '   ISOPOOH + OH --> 0.5 XO2 + 0.5 ISOPO2                                                            ', &
     ' ISOPO2 + CH3O2 --> 0.3 HYDRALD + 0.25 CH3OH + 0.19 MACR + 0.26 MVK + 1.2 CH2O + HO2                ', &
     'ISOPO2 + CH3CO3 --> 0.4 HYDRALD + 0.25 MACR + 0.35 MVK + 0.6 CH2O + CH3O2 + HO2                     ', &
     '       MVK + OH --> MACRO2                                                                          ', &
     '       MVK + O3 --> 0.05 CO + 0.04 CH3CHO + 0.95 CH3COCHO + 0.8 CH2O + 0.2 O3 + 0.08 OH + 0.06 HO2  ', &
     '      MACR + OH --> 0.5 MACRO2 + 0.5 MCO3                                                           ', &
     '      MACR + O3 --> 0.2 CO + 0.8 CH3COCHO + 0.7 CH2O + 0.2 O3 + 0.215 OH + 0.275 HO2                ', &
     '    MACRO2 + NO --> 0.53 GLYALD + 0.22 CO + 0.22 HYAC + 0.25 CH3COCHO + 0.25 CH2O + 0.47 HO2 + NO2 +', &
     '    MACRO2 + NO --> 0.8 ONITR                                                                       ', &
     '   MACRO2 + NO3 --> 0.53 GLYALD + 0.22 CO + 0.22 HYAC + 0.25 CH3COCHO + 0.25 CH2O + 0.47 HO2 + NO2 +', &
     '   MACRO2 + HO2 --> MACROOH                                                                         ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_5 = (/ &
     ' MACRO2 + CH3O2 --> 0.25 CH3OH + 0.26 GLYALD + 0.11 CO + 0.23 HYAC + 0.24 CH3COCHO + 0.88 CH2O + 0.7', &
     'MACRO2 + CH3CO3 --> 0.53 GLYALD + 0.22 CO + 0.22 HYAC + 0.25 CH3COCHO + 0.25 CH2O + CH3O2 + 0.47 HO2', &
     '   MACROOH + OH --> 0.2 MACRO2 + 0.1 OH + 0.5 MCO3 + 0.2 HO2                                        ', &
     '      MCO3 + NO --> CH2O + NO2 + CH3CO3                                                             ', &
     '     NO3 + MCO3 --> CH2O + NO2 + CH3CO3                                                             ', &
     '     MCO3 + HO2 --> 0.75 CH3COOOH + 0.25 CH3COOH + 0.25 O3                                          ', &
     '   CH3O2 + MCO3 --> 2 CH2O + HO2 + CH3CO3                                                           ', &
     '  MCO3 + CH3CO3 --> CH2O + CH3O2 + CH3CO3                                                           ', &
     '         2 MCO3 --> 2 CH2O + 2 CH3CO3                                                               ', &
     ' MCO3 + NO2 + M --> MPAN                                                                            ', &
     '       MPAN + M --> MCO3 + NO2                                                                      ', &
     '    C10H16 + OH --> TERPO2                                                                          ', &
     '    C10H16 + O3 --> MACR + MVK + 0.7 OH + HO2                                                       ', &
     '   C10H16 + NO3 --> TERPO2 + NO2                                                                    ', &
     '    TERPO2 + NO --> 0.1 CH3COCH3 + MACR + MVK + HO2 + NO2                                           ', &
     '   TERPO2 + HO2 --> TERPOOH                                                                         ', &
     '   TERPOOH + OH --> TERPO2                                                                          ', &
     '   CH3COOH + OH --> CH3O2                                                                           ', &
     '     ISOP + NO3 --> ISOPNO3                                                                         ', &
     '   ISOPNO3 + NO --> 0.794 ONITR + 0.167 MACR + 0.039 MVK + 0.072 CH2O + 0.794 HO2 + 1.206 NO2       ', &
     '  ISOPNO3 + NO3 --> 0.794 ONITR + 0.167 MACR + 0.039 MVK + 0.072 CH2O + 0.794 HO2 + 1.206 NO2       ', &
     '  ISOPNO3 + HO2 --> 0.794 ONITR + 0.167 MACR + 0.039 MVK + 0.008 CH2O + 0.794 HO2 + 0.206 NO2       ', &
     '  CH3COCHO + OH --> CO + CH3CO3                                                                     ', &
     ' CH3COCHO + NO3 --> HNO3 + CO + CH3CO3                                                              ', &
     '     ONITR + OH --> HYDRALD + HO2 + 0.4 NO2                                                         ', &
     '    ONITR + NO3 --> HYDRALD + HO2 + NO2                                                             ', &
     '   HYDRALD + OH --> XO2                                                                             ', &
     '       XO2 + NO --> 0.25 GLYALD + CO + 0.25 HYAC + 0.25 CH3COCHO + 1.5 HO2 + NO2                    ', &
     '      XO2 + NO3 --> 0.25 GLYALD + CO + 0.25 HYAC + 0.25 CH3COCHO + 1.5 HO2 + NO2                    ', &
     '      XO2 + HO2 --> XOOH                                                                            ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(20) :: EQN_NAMES_6 = (/ &
     '    XO2 + CH3O2 --> 0.3 CH3OH + 0.1 GLYALD + 0.4 CO + 0.1 HYAC + 0.1 CH3COCHO + 0.7 CH2O + HO2      ', &
     '   XO2 + CH3CO3 --> 0.25 GLYALD + CO + 0.25 HYAC + 0.25 CH3COCHO + CH3O2 + 1.5 HO2                  ', &
     '      XOOH + OH --> XO2                                                                             ', &
     '      XOOH + OH --> OH                                                                              ', &
     '     CH3OH + OH --> CH2O + HO2                                                                      ', &
     '    C2H5OH + OH --> CH3CHO + HO2                                                                    ', &
     '      MPAN + OH --> 0.5 HYAC + 0.5 CH2O + 0.5 NO3 + 0.5 HO2                                         ', &
     '       PAN + OH --> CH2O + NO3                                                                      ', &
     '      HYAC + OH --> CH3COCHO + HO2                                                                  ', &
     '    GLYALD + OH --> 0.2 GLYOXAL + 0.8 CH2O + HO2                                                    ', &
     '       DMS + OH --> SO2                                                                             ', &
     '       DMS + OH --> 0.5 SO2 + 0.5 HO2                                                               ', &
     '      DMS + NO3 --> SO2 + HNO3                                                                      ', &
     '       NH3 + OH --> M                                                                               ', &
     '            HO2 --> 0.5 H2O2                                                                        ', &
     '       2 C2H5O2 --> 0.4 C2H5OH + 1.6 CH3CHO + 1.2 HO2                                               ', &
     '      VOCA + OH --> SMPA + OH                                                                       ', &
     '     VOCBB + OH --> SMPBB + OH                                                                      ', &
     '      ISOP + OH --> BIOG1_c + ISOP + OH                                                             ', &
     '    C10H16 + OH --> BIOG1_o + C10H16 + OH                                                           ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(200) :: EQN_NAMES = (/&
    EQN_NAMES_0, EQN_NAMES_1, EQN_NAMES_2, EQN_NAMES_3, EQN_NAMES_4, &
    EQN_NAMES_5, EQN_NAMES_6 /)

! INLINED global variables

! End INLINED global variables


END MODULE mozart_mosaic_4bin_vbs0_Monitor
