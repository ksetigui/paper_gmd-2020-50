MODULE nmhc9_Parameters
  USE nmhc9_Precision
  PUBLIC
  SAVE
  INTEGER, PARAMETER :: NSPEC = 59
  INTEGER, PARAMETER :: NVAR = 57
  INTEGER, PARAMETER :: NVARACT = 57
  INTEGER, PARAMETER :: NFIX = 2
  INTEGER, PARAMETER :: NREACT = 145
  INTEGER, PARAMETER :: NVARST = 1
  INTEGER, PARAMETER :: NFIXST = 58
  INTEGER, PARAMETER :: NONZERO = 527
  INTEGER, PARAMETER :: LU_NONZERO = 580
  INTEGER, PARAMETER :: CNVAR = 58
  INTEGER, PARAMETER :: NLOOKAT = 0
  INTEGER, PARAMETER :: NMONITOR = 0
  INTEGER, PARAMETER :: NMASS = 1
  REAL(kind=dp), PARAMETER :: PI = 3.14159265358979
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
  INTEGER, PARAMETER :: ind_H2O = 58
  INTEGER, PARAMETER :: ind_M = 59
  INTEGER, PARAMETER :: indf_H2O = 1
  INTEGER, PARAMETER :: indf_M = 2
END MODULE nmhc9_Parameters