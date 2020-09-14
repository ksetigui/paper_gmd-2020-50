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
! File                 : saprc99_mosaic_8bin_vbs2_aq_Monitor.f90
! Time                 : Mon Feb 27 11:09:03 2017
! Working directory    : /data/ksetigui/setigui/WRFV3_171015/chem/KPP/mechanisms/saprc99_mosaic_8bin_vbs2_aq
! Equation file        : saprc99_mosaic_8bin_vbs2_aq.kpp
! Output root filename : saprc99_mosaic_8bin_vbs2_aq
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE saprc99_mosaic_8bin_vbs2_aq_Monitor


  CHARACTER(LEN=12), PARAMETER, DIMENSION(90) :: SPC_NAMES_0 = (/ &
     'H2SO4       ','HCOOH       ','CCO_OH      ', &
     'RCO_OH      ','CO2         ','CCO_OOH     ', &
     'RCO_OOH     ','XN          ','XC          ', &
     'NUME        ','DEN         ','ANT1_c      ', &
     'ANT2_c      ','ANT3_c      ','ANT4_c      ', &
     'BIOG1_c     ','BIOG2_c     ','BIOG3_c     ', &
     'BIOG4_c     ','PSD1        ','O1D         ', &
     'CH4         ','C2H6        ','PAN         ', &
     'PAN2        ','PBZN        ','MA_PAN      ', &
     'H2O2        ','C3H8        ','BACL        ', &
     'ETOH        ','N2O5        ','SO2         ', &
     'DMS         ','HONO        ','ALK3        ', &
     'TBU_O       ','ALK5        ','ARO2        ', &
     'COOH        ','HOCOO       ','BZNO2_O     ', &
     'HNO4        ','MEOH        ','ARO1        ', &
     'DCB3        ','DCB2        ','CRES        ', &
     'ALK4        ','DCB1        ','C2H2        ', &
     'ROOH        ','NPHE        ','BALD        ', &
     'PHEN        ','MGLY        ','CO          ', &
     'HNO3        ','ETHENE      ','ACET        ', &
     'C3H6        ','BZ_O        ','OLE1        ', &
     'ISOPRENE    ','R2O2        ','METHACRO    ', &
     'TERP        ','SESQ        ','GLY         ', &
     'ISOPROD     ','OLE2        ','MVK         ', &
     'CCHO        ','HCHO        ','RNO3        ', &
     'O3P         ','RCHO        ','MEK         ', &
     'PROD2       ','O3          ','RCO_O2      ', &
     'NO2         ','HO2         ','BZCO_O2     ', &
     'C_O2        ','NO3         ','NO          ', &
     'RO2_R       ','RO2_N       ','CCO_O2      ' /)
  CHARACTER(LEN=12), PARAMETER, DIMENSION(4) :: SPC_NAMES_1 = (/ &
     'MA_RCO3     ','OH          ','H2O         ', &
     'M           ' /)
  CHARACTER(LEN=12), PARAMETER, DIMENSION(94) :: SPC_NAMES = (/&
    SPC_NAMES_0, SPC_NAMES_1 /)

  INTEGER, DIMENSION(1) :: LOOKAT
  INTEGER, DIMENSION(1) :: MONITOR
  CHARACTER(LEN=12), DIMENSION(1) :: SMASS
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_0 = (/ &
     '              NO2 --> O3P + NO                                                                      ', &
     '          O3P + M --> O3                                                                            ', &
     '         O3P + O3 --> M                                                                             ', &
     '     O3P + NO + M --> NO2                                                                           ', &
     '        O3P + NO2 --> NO                                                                            ', &
     '        O3P + NO2 --> NO3                                                                           ', &
     '          O3 + NO --> NO2                                                                           ', &
     '         O3 + NO2 --> NO3                                                                           ', &
     '         NO3 + NO --> 2 NO2                                                                         ', &
     '         2 NO + M --> 2 NO2                                                                         ', &
     '        NO2 + NO3 --> N2O5                                                                          ', &
     '             N2O5 --> NO2 + NO3                                                                     ', &
     '       N2O5 + H2O --> 2 HNO3                                                                        ', &
     '        NO2 + NO3 --> NO2 + NO                                                                      ', &
     '              NO3 --> NO                                                                            ', &
     '              NO3 --> O3P + NO2                                                                     ', &
     '               O3 --> O3P                                                                           ', &
     '               O3 --> O1D                                                                           ', &
     '        O1D + H2O --> 2 OH                                                                          ', &
     '          O1D + M --> O3P                                                                           ', &
     '          NO + OH --> PSD1 + HONO                                                                   ', &
     '             HONO --> NO + OH                                                                       ', &
     '             HONO --> NO2 + HO2                                                                     ', &
     '        HONO + OH --> PSD1 + NO2                                                                    ', &
     '         NO2 + OH --> PSD1 + HNO3                                                                   ', &
     '         NO3 + OH --> PSD1 + NO2 + HO2                                                              ', &
     '        HNO3 + OH --> PSD1 + NO3                                                                    ', &
     '             HNO3 --> NO2 + OH                                                                      ', &
     '          CO + OH --> PSD1 + HO2                                                                    ', &
     '          O3 + OH --> PSD1 + HO2                                                                    ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_1 = (/ &
     '         HO2 + NO --> NO2 + OH                                                                      ', &
     '        NO2 + HO2 --> HNO4                                                                          ', &
     '             HNO4 --> NO2 + HO2                                                                     ', &
     '             HNO4 --> 0.61 NO2 + 0.61 HO2 + 0.39 NO3 + 0.39 OH                                      ', &
     '        HNO4 + OH --> NO2                                                                           ', &
     '         O3 + HO2 --> OH                                                                            ', &
     '            2 HO2 --> H2O2                                                                          ', &
     '      2 HO2 + H2O --> H2O2                                                                          ', &
     '        HO2 + NO3 --> 0.2 HNO3 + 0.8 NO2 + 0.8 OH                                                   ', &
     '            2 NO3 --> 2 NO2                                                                         ', &
     '             H2O2 --> 2 OH                                                                          ', &
     '        H2O2 + OH --> HO2                                                                           ', &
     '         HO2 + OH --> PSD1 + H2O + M                                                                ', &
     '         SO2 + OH --> H2SO4 + PSD1 + HO2                                                            ', &
     '           OH + M --> PSD1 + HO2                                                                    ', &
     '        C_O2 + NO --> NUME + HCHO + NO2 + HO2                                                       ', &
     '       HO2 + C_O2 --> DEN + COOH                                                                    ', &
     '       C_O2 + NO3 --> NUME + HCHO + NO2 + HO2                                                       ', &
     '           2 C_O2 --> DEN + MEOH + HCHO                                                             ', &
     '           2 C_O2 --> DEN + 2 HCHO + 2 HO2                                                          ', &
     '       NO + RO2_R --> NUME + NO2 + HO2                                                              ', &
     '      HO2 + RO2_R --> DEN + ROOH                                                                    ', &
     '      NO3 + RO2_R --> NUME + NO2 + HO2                                                              ', &
     '     C_O2 + RO2_R --> DEN + 0.25 MEOH + 0.75 HCHO + HO2                                             ', &
     '          2 RO2_R --> DEN + HO2                                                                     ', &
     '        R2O2 + NO --> NUME + NO2                                                                    ', &
     '       R2O2 + HO2 --> DEN + HO2                                                                     ', &
     '       R2O2 + NO3 --> NUME + NO2                                                                    ', &
     '      R2O2 + C_O2 --> DEN + C_O2                                                                    ', &
     '     R2O2 + RO2_R --> DEN + RO2_R                                                                   ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_2 = (/ &
     '           2 R2O2 --> DEN + 2 R2O2                                                                  ', &
     '       NO + RO2_N --> NUME + RNO3                                                                   ', &
     '      HO2 + RO2_N --> DEN + ROOH                                                                    ', &
     '     C_O2 + RO2_N --> DEN + 0.25 MEOH + 0.75 HCHO + 0.5 MEK + 0.5 PROD2 + HO2                       ', &
     '      NO3 + RO2_N --> NUME + MEK + NO2 + HO2                                                        ', &
     '    RO2_R + RO2_N --> DEN + 0.5 MEK + 0.5 PROD2 + HO2                                               ', &
     '     R2O2 + RO2_N --> DEN + RO2_N                                                                   ', &
     '          2 RO2_N --> DEN + MEK + PROD2 + HO2                                                       ', &
     '     NO2 + CCO_O2 --> PAN                                                                           ', &
     '              PAN --> NO2 + CCO_O2                                                                  ', &
     '      NO + CCO_O2 --> NO2 + C_O2                                                                    ', &
     '     HO2 + CCO_O2 --> 0.25 CCO_OH + 0.75 CCO_OOH + 0.25 O3                                          ', &
     '     NO3 + CCO_O2 --> NO2 + C_O2                                                                    ', &
     '    C_O2 + CCO_O2 --> CCO_OH + HCHO                                                                 ', &
     '   RO2_R + CCO_O2 --> CCO_OH                                                                        ', &
     '    R2O2 + CCO_O2 --> CCO_O2                                                                        ', &
     '   RO2_N + CCO_O2 --> CCO_OH + PROD2                                                                ', &
     '         2 CCO_O2 --> 2 C_O2                                                                        ', &
     '     RCO_O2 + NO2 --> PAN2                                                                          ', &
     '             PAN2 --> RCO_O2 + NO2                                                                  ', &
     '      RCO_O2 + NO --> CCHO + NO2 + RO2_R                                                            ', &
     '     RCO_O2 + HO2 --> 0.25 RCO_OH + 0.75 RCO_OOH + 0.25 O3                                          ', &
     '     RCO_O2 + NO3 --> CCHO + NO2 + RO2_R                                                            ', &
     '    RCO_O2 + C_O2 --> RCO_OH + HCHO                                                                 ', &
     '   RCO_O2 + RO2_R --> RCO_OH                                                                        ', &
     '    R2O2 + RCO_O2 --> RCO_O2                                                                        ', &
     '   RCO_O2 + RO2_N --> RCO_OH + PROD2                                                                ', &
     '  RCO_O2 + CCO_O2 --> CCHO + C_O2 + RO2_R                                                           ', &
     '         2 RCO_O2 --> 2 CCHO + 2 RO2_R                                                              ', &
     '    NO2 + BZCO_O2 --> PBZN                                                                          ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_3 = (/ &
     '             PBZN --> NO2 + BZCO_O2                                                                 ', &
     '     BZCO_O2 + NO --> BZ_O + R2O2 + NO2                                                             ', &
     '    HO2 + BZCO_O2 --> 0.25 RCO_OH + 0.75 RCO_OOH + 0.25 O3                                          ', &
     '    BZCO_O2 + NO3 --> BZ_O + R2O2 + NO2                                                             ', &
     '   BZCO_O2 + C_O2 --> RCO_OH + HCHO                                                                 ', &
     '  BZCO_O2 + RO2_R --> RCO_OH                                                                        ', &
     '   R2O2 + BZCO_O2 --> BZCO_O2                                                                       ', &
     '  BZCO_O2 + RO2_N --> RCO_OH + PROD2                                                                ', &
     ' BZCO_O2 + CCO_O2 --> BZ_O + R2O2 + C_O2                                                            ', &
     ' RCO_O2 + BZCO_O2 --> BZ_O + R2O2 + CCHO + RO2_R                                                    ', &
     '        2 BZCO_O2 --> 2 BZ_O + 2 R2O2                                                               ', &
     '    NO2 + MA_RCO3 --> MA_PAN                                                                        ', &
     '           MA_PAN --> NO2 + MA_RCO3                                                                 ', &
     '     NO + MA_RCO3 --> HCHO + NO2 + CCO_O2                                                           ', &
     '    HO2 + MA_RCO3 --> 0.25 RCO_OH + 0.75 RCO_OOH + 0.25 O3                                          ', &
     '    NO3 + MA_RCO3 --> HCHO + NO2 + CCO_O2                                                           ', &
     '   C_O2 + MA_RCO3 --> RCO_OH + HCHO                                                                 ', &
     '  RO2_R + MA_RCO3 --> RCO_OH                                                                        ', &
     '   R2O2 + MA_RCO3 --> MA_RCO3                                                                       ', &
     '  RO2_N + MA_RCO3 --> 2 RCO_OH                                                                      ', &
     ' CCO_O2 + MA_RCO3 --> HCHO + C_O2 + CCO_O2                                                          ', &
     ' RCO_O2 + MA_RCO3 --> CCHO + HCHO + RO2_R + CCO_O2                                                  ', &
     'BZCO_O2 + MA_RCO3 --> BZ_O + R2O2 + HCHO + CCO_O2                                                   ', &
     '        2 MA_RCO3 --> 2 HCHO + 2 CCO_O2                                                             ', &
     '      TBU_O + NO2 --> RNO3                                                                          ', &
     '            TBU_O --> ACET + C_O2                                                                   ', &
     '       BZ_O + NO2 --> NPHE                                                                          ', &
     '       BZ_O + HO2 --> PHEN                                                                          ', &
     '             BZ_O --> PHEN                                                                          ', &
     '    BZNO2_O + NO2 --> 2 XN + 6 XC                                                                   ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_4 = (/ &
     '    BZNO2_O + HO2 --> NPHE                                                                          ', &
     '          BZNO2_O --> NPHE                                                                          ', &
     '             HCHO --> CO + 2 HO2                                                                    ', &
     '             HCHO --> CO                                                                            ', &
     '        HCHO + OH --> PSD1 + CO + HO2                                                               ', &
     '       HCHO + HO2 --> HOCOO                                                                         ', &
     '            HOCOO --> HCHO + HO2                                                                    ', &
     '       HOCOO + NO --> HCOOH + NO2 + HO2                                                             ', &
     '       HCHO + NO3 --> CO + HNO3 + HO2                                                               ', &
     '        CCHO + OH --> PSD1 + CCO_O2                                                                 ', &
     '             CCHO --> CO + HO2 + C_O2                                                               ', &
     '       CCHO + NO3 --> HNO3 + CCO_O2                                                                 ', &
     '        RCHO + OH --> PSD1 + 0.034 CO + 0.034 CCHO + 0.965 RCO_O2 + 0.034 RO2_R + 0.001 RO2_N       ', &
     '             RCHO --> CO + CCHO + HO2 + RO2_R                                                       ', &
     '       RCHO + NO3 --> HNO3 + RCO_O2                                                                 ', &
     '        ACET + OH --> PSD1 + R2O2 + HCHO + CCO_O2                                                   ', &
     '             ACET --> C_O2 + CCO_O2                                                                 ', &
     '         MEK + OH --> PSD1 + 0.616 R2O2 + 0.482 CCHO + 0.115 HCHO + 0.37 RCHO + 0.096 RCO_O2 + 0.37 ', &
     '              MEK --> CCHO + RO2_R + CCO_O2                                                         ', &
     '        MEOH + OH --> PSD1 + HCHO + HO2                                                             ', &
     '        ETOH + OH --> PSD1 + 0.96 CCHO + 0.081 HCHO + 0.95 HO2 + 0.05 RO2_R                         ', &
     '        COOH + OH --> PSD1 + 0.35 HCHO + 0.65 C_O2 + 0.35 OH                                        ', &
     '             COOH --> HCHO + HO2 + OH                                                               ', &
     '        ROOH + OH --> PSD1 + RCHO + 0.34 RO2_R + 0.66 OH                                            ', &
     '             ROOH --> RCHO + HO2 + OH                                                               ', &
     '              GLY --> 2 CO + 2 HO2                                                                  ', &
     '              GLY --> CO + HCHO                                                                     ', &
     '         GLY + OH --> PSD1 + 1.26 CO + 0.37 RCO_O2 + 0.63 HO2                                       ', &
     '        GLY + NO3 --> 1.26 CO + HNO3 + 0.37 RCO_O2 + 0.63 HO2                                       ', &
     '             MGLY --> CO + HO2 + CCO_O2                                                             ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_5 = (/ &
     '        MGLY + OH --> PSD1 + CO + CCO_O2                                                            ', &
     '       MGLY + NO3 --> CO + HNO3 + CCO_O2                                                            ', &
     '             BACL --> 2 CCO_O2                                                                      ', &
     '        PHEN + OH --> PSD1 + 0.24 BZ_O + 0.23 GLY + 0.76 RO2_R                                      ', &
     '       PHEN + NO3 --> HNO3 + BZ_O                                                                   ', &
     '        CRES + OH --> PSD1 + 0.23 MGLY + 0.24 BZ_O + 0.76 RO2_R                                     ', &
     '       CRES + NO3 --> HNO3 + BZ_O                                                                   ', &
     '       NPHE + NO3 --> BZNO2_O + HNO3                                                                ', &
     '        BALD + OH --> PSD1 + BZCO_O2                                                                ', &
     '             BALD --> 7 XC                                                                          ', &
     '       BALD + NO3 --> HNO3 + BZCO_O2                                                                ', &
     '    METHACRO + OH --> PSD1 + 0.084 MGLY + 0.416 CO + 0.084 HCHO + 0.416 MEK + 0.5 RO2_R + 0.5 MA_RCO', &
     '    METHACRO + O3 --> 0.333 HCOOH + 0.9 MGLY + 0.45 CO + 0.2 HCHO + 0.1 RCO_O2 + 0.008 HO2 + 0.1 RO2', &
     '   METHACRO + NO3 --> 0.5 CO + 0.5 HNO3 + 0.5 RO2_R + 0.5 MA_RCO3                                   ', &
     '   METHACRO + O3P --> RCHO                                                                          ', &
     '         METHACRO --> 0.67 CO + 0.67 HCHO + 0.34 HO2 + 0.33 RO2_R + 0.67 CCO_O2 + 0.33 MA_RCO3 + 0.3', &
     '         MVK + OH --> PSD1 + 0.3 MGLY + 0.675 R2O2 + 0.3 HCHO + 0.675 RCHO + 0.3 RO2_R + 0.025 RO2_N', &
     '         MVK + O3 --> 0.351 HCOOH + 0.95 MGLY + 0.475 CO + 0.1 HCHO + 0.05 RCO_O2 + 0.064 HO2 + 0.05', &
     '        MVK + O3P --> 0.45 RCHO + 0.55 MEK                                                          ', &
     '              MVK --> 0.7 CO + 0.7 PROD2 + 0.3 C_O2 + 0.3 MA_RCO3                                   ', &
     '     ISOPROD + OH --> PSD1 + 0.174 MGLY + 0.336 CO + 0.15 GLY + 0.129 CCHO + 0.055 HCHO + 0.013 RCHO', &
     '     ISOPROD + O3 --> 0.1 HCOOH + 0.372 RCO_OH + 0.742 MGLY + 0.498 CO + 0.023 GLY + 0.047 CCHO + 0.', &
     '    ISOPROD + NO3 --> 0.008 MGLY + 0.572 CO + 0.15 HNO3 + 0.227 HCHO + 0.572 RNO3 + 0.218 RCHO + 0.7', &
     '          ISOPROD --> 1.233 CO + 0.467 CCHO + 0.3 HCHO + 0.233 MEK + 0.3 RCO_O2 + 1.233 HO2 + 0.467 ', &
     '       PROD2 + OH --> PSD1 + 0.084 CCHO + 0.213 HCHO + 0.558 RCHO + 0.115 MEK + 0.329 PROD2 + 0.049 ', &
     '            PROD2 --> 0.515 R2O2 + 0.246 CCHO + 0.506 HCHO + 0.71 RCHO + 0.333 RCO_O2 + 0.96 RO2_R +', &
     '        RNO3 + OH --> PSD1 + 0.006 ACET + 0.596 R2O2 + 0.439 CCHO + 0.01 HCHO + 0.31 RNO3 + 0.213 RC', &
     '             RNO3 --> 0.02 ACET + 0.152 R2O2 + 0.431 CCHO + 0.134 HCHO + 0.147 RCHO + 0.243 MEK + 0.', &
     '        DCB1 + OH --> PSD1 + CO + RCHO + RO2_R                                                      ', &
     '        DCB1 + O3 --> 1.5 CO + GLY + 1.5 HO2 + 0.5 OH                                               ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_6 = (/ &
     '        DCB2 + OH --> PSD1 + R2O2 + RCHO + CCO_O2                                                   ', &
     '             DCB2 --> 0.5 MGLY + CO + R2O2 + 0.5 GLY + 0.5 HO2 + RO2_R + 0.5 CCO_O2                 ', &
     '        DCB3 + OH --> PSD1 + R2O2 + RCHO + CCO_O2                                                   ', &
     '             DCB3 --> 0.5 MGLY + CO + R2O2 + 0.5 GLY + 0.5 HO2 + RO2_R + 0.5 CCO_O2                 ', &
     '         CH4 + OH --> PSD1 + C_O2                                                                   ', &
     '      ETHENE + OH --> PSD1 + 0.195 CCHO + 1.61 HCHO + RO2_R                                         ', &
     '      ETHENE + O3 --> 0.37 HCOOH + 0.5 CO + HCHO + 0.12 HO2 + 0.12 OH                               ', &
     '     ETHENE + NO3 --> RCHO + RO2_R                                                                  ', &
     '     ETHENE + O3P --> 0.491 CO + 0.009 GLY + 0.25 CCHO + 0.191 HCHO + 0.5 HO2 + 0.3 C_O2 + 0.2 RO2_R', &
     '    ISOPRENE + OH --> PSD1 + 0.079 R2O2 + 0.23 METHACRO + 0.357 ISOPROD + 0.32 MVK + 0.624 HCHO + 0.', &
     '    ISOPRENE + O3 --> 0.204 HCOOH + 0.15 RCO_OH + 0.275 CO + 0.126 R2O2 + 0.39 METHACRO + 0.16 MVK +', &
     '   ISOPRENE + NO3 --> 0.187 R2O2 + 0.936 ISOPROD + 0.187 NO2 + 0.749 RO2_R + 0.064 RO2_N            ', &
     '   ISOPRENE + O3P --> 0.24 R2O2 + 0.24 HCHO + 0.75 PROD2 + 0.25 C_O2 + 0.01 RO2_N + 0.24 MA_RCO3    ', &
     '        TERP + OH --> PSD1 + 0.5 R2O2 + 0.276 HCHO + 0.474 RCHO + 0.276 PROD2 + 0.75 RO2_R + 0.25 RO', &
     '        TERP + O3 --> 0.103 HCOOH + 0.189 RCO_OH + 0.031 BACL + 0.157 CO + 0.13 ACET + 0.729 R2O2 + ', &
     '       TERP + NO3 --> 0.75 R2O2 + 0.276 RNO3 + 0.474 RCHO + 0.474 NO2 + 0.276 RO2_R + 0.25 RO2_N    ', &
     '       TERP + O3P --> 0.147 RCHO + 0.853 PROD2                                                      ', &
     '        SESQ + OH --> PSD1 + 0.5 R2O2 + 0.276 HCHO + 0.474 RCHO + 0.276 PROD2 + 0.75 RO2_R + 0.25 RO', &
     '        SESQ + O3 --> 0.103 HCOOH + 0.189 RCO_OH + 0.031 BACL + 0.157 CO + 0.13 ACET + 0.729 R2O2 + ', &
     '       SESQ + NO3 --> 0.75 R2O2 + 0.276 RNO3 + 0.474 RCHO + 0.474 NO2 + 0.276 RO2_R + 0.25 RO2_N    ', &
     '       SESQ + O3P --> 0.147 RCHO + 0.853 PROD2                                                      ', &
     '        C2H6 + OH --> PSD1 + CCHO + RO2_R                                                           ', &
     '        C3H8 + OH --> PSD1 + 0.704 ACET + 0.261 RCHO + 0.965 RO2_R + 0.035 RO2_N                    ', &
     '        C2H2 + OH --> 0.297 HCOOH + PSD1 + 0.393 CO + 0.607 GLY + 0.096 HCHO + 0.297 HO2 + 0.1 RO2_R', &
     '        ALK3 + OH --> PSD1 + 0.236 TBU_O + 0.024 ACET + 0.559 R2O2 + 0.445 CCHO + 0.026 HCHO + 0.122', &
     '        ALK4 + OH --> PSD1 + 0.002 CO + 0.452 ACET + 0.936 R2O2 + 0.455 CCHO + 0.024 HCHO + 0.244 RC', &
     '        ALK5 + OH --> PSD1 + 0.072 ACET + 0.948 R2O2 + 0.099 CCHO + 0.026 HCHO + 0.204 RCHO + 0.089 ', &
     '        ARO1 + OH --> PSD1 + 0.051 DCB3 + 0.108 DCB2 + 0.207 CRES + 0.491 DCB1 + 0.059 BALD + 0.017 ', &
     '        ARO2 + OH --> PSD1 + 0.087 BACL + 0.093 DCB3 + 0.099 DCB2 + 0.187 CRES + 0.561 DCB1 + 0.05 B', &
     '        OLE1 + OH --> PSD1 + 0.005 ACET + 0.205 R2O2 + 0.294 CCHO + 0.732 HCHO + 0.497 RCHO + 0.119 ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(30) :: EQN_NAMES_7 = (/ &
     '        OLE1 + O3 --> 0.185 HCOOH + 0.05 CCO_OH + 0.119 RCO_OH + 0.345 CO + 0.001 ACET + 0.154 CCHO ', &
     '       OLE1 + NO3 --> 0.024 ACET + 0.488 R2O2 + 0.009 CCHO + 0.511 RNO3 + 0.037 RCHO + 0.824 RO2_R +', &
     '       OLE1 + O3P --> 0.45 RCHO + 0.437 MEK + 0.113 PROD2                                           ', &
     '        OLE2 + OH --> PSD1 + 0.061 BALD + 0.127 ACET + 0.001 R2O2 + 0.025 METHACRO + 0.025 ISOPROD +', &
     '        OLE2 + O3 --> 0.073 HCOOH + 0.129 CCO_OH + 0.247 RCO_OH + 0.042 BALD + 0.265 CO + 0.045 ACET', &
     '       OLE2 + NO3 --> 0.015 BALD + 0.102 ACET + 0.711 R2O2 + 0.048 MVK + 0.507 CCHO + 0.079 HCHO + 0', &
     '       OLE2 + O3P --> 0.012 CO + 0.012 METHACRO + 0.069 RCHO + 0.659 MEK + 0.259 PROD2 + 0.013 HO2 +', &
     '        C2H2 + O3 --> 0.5 CO2 + 1.5 CO + 1.5 HO2 + 0.5 OH                                           ', &
     '        C3H6 + OH --> 0.048 XC + PSD1 + 0.984 CCHO + 0.984 HCHO + 0.984 RO2_R + 0.016 RO2_N         ', &
     '        C3H6 + O3 --> 0.185 HCOOH + 0.17 CCO_OH + 0.135 CO2 + 0.07 XC + 0.51 CO + 0.5 CCHO + 0.5 HCH', &
     '       C3H6 + NO3 --> XN + 2.693 XC + 0.949 RO2_R + 0.051 RO2_N                                     ', &
     '       C3H6 + O3P --> 0.55 XC + 0.45 RCHO + 0.55 MEK                                                ', &
     '        ARO1 + OH --> ANT1_c + ARO1 + OH                                                            ', &
     '        ARO1 + OH --> ANT2_c + ARO1 + OH                                                            ', &
     '        ARO1 + OH --> ANT3_c + ARO1 + OH                                                            ', &
     '        ARO1 + OH --> ANT4_c + ARO1 + OH                                                            ', &
     '        ARO2 + OH --> ANT1_c + ARO2 + OH                                                            ', &
     '        ARO2 + OH --> ANT2_c + ARO2 + OH                                                            ', &
     '        ARO2 + OH --> ANT3_c + ARO2 + OH                                                            ', &
     '        ARO2 + OH --> ANT4_c + ARO2 + OH                                                            ', &
     '    ISOPRENE + OH --> BIOG1_c + ISOPRENE + OH                                                       ', &
     '    ISOPRENE + OH --> BIOG2_c + ISOPRENE + OH                                                       ', &
     '    ISOPRENE + OH --> BIOG3_c + ISOPRENE + OH                                                       ', &
     '        TERP + OH --> BIOG2_c + TERP + OH                                                           ', &
     '        TERP + OH --> BIOG3_c + TERP + OH                                                           ', &
     '        TERP + OH --> BIOG4_c + TERP + OH                                                           ', &
     '        SESQ + OH --> BIOG2_c + SESQ + OH                                                           ', &
     '        SESQ + OH --> BIOG3_c + SESQ + OH                                                           ', &
     '    ISOPRENE + O3 --> BIOG3_c + ISOPRENE + O3                                                       ', &
     '        TERP + O3 --> BIOG2_c + TERP + O3                                                           ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(15) :: EQN_NAMES_8 = (/ &
     '        TERP + O3 --> BIOG3_c + TERP + O3                                                           ', &
     '        TERP + O3 --> BIOG4_c + TERP + O3                                                           ', &
     '        SESQ + O3 --> BIOG2_c + SESQ + O3                                                           ', &
     '        SESQ + O3 --> BIOG3_c + SESQ + O3                                                           ', &
     '        SESQ + O3 --> BIOG4_c + SESQ + O3                                                           ', &
     '   ISOPRENE + NO3 --> BIOG2_c + ISOPRENE + NO3                                                      ', &
     '   ISOPRENE + NO3 --> BIOG3_c + ISOPRENE + NO3                                                      ', &
     '   ISOPRENE + NO3 --> BIOG4_c + ISOPRENE + NO3                                                      ', &
     '       TERP + NO3 --> BIOG3_c + TERP + NO3                                                          ', &
     '       TERP + NO3 --> BIOG4_c + TERP + NO3                                                          ', &
     '       SESQ + NO3 --> BIOG3_c + SESQ + NO3                                                          ', &
     '       SESQ + NO3 --> BIOG4_c + SESQ + NO3                                                          ', &
     '         DMS + OH --> SO2                                                                           ', &
     '         DMS + OH --> 0.5 SO2 + 0.5 HO2                                                             ', &
     '        DMS + NO3 --> SO2 + HNO3                                                                    ' /)
  CHARACTER(LEN=100), PARAMETER, DIMENSION(255) :: EQN_NAMES = (/&
    EQN_NAMES_0, EQN_NAMES_1, EQN_NAMES_2, EQN_NAMES_3, EQN_NAMES_4, &
    EQN_NAMES_5, EQN_NAMES_6, EQN_NAMES_7, EQN_NAMES_8 /)

! INLINED global variables

! End INLINED global variables


END MODULE saprc99_mosaic_8bin_vbs2_aq_Monitor
