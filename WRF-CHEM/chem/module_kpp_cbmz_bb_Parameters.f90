MODULE cbmz_bb_Parameters
  USE cbmz_bb_Precision
  PUBLIC
  SAVE
  INTEGER, PARAMETER :: NSPEC = 58
  INTEGER, PARAMETER :: NVAR = 56
  INTEGER, PARAMETER :: NVARACT = 54
  INTEGER, PARAMETER :: NFIX = 2
  INTEGER, PARAMETER :: NREACT = 136
  INTEGER, PARAMETER :: NVARST = 1
  INTEGER, PARAMETER :: NFIXST = 57
  INTEGER, PARAMETER :: NONZERO = 509
  INTEGER, PARAMETER :: LU_NONZERO = 576
  INTEGER, PARAMETER :: CNVAR = 57
  INTEGER, PARAMETER :: NLOOKAT = 0
  INTEGER, PARAMETER :: NMONITOR = 0
  INTEGER, PARAMETER :: NMASS = 1
  REAL(kind=dp), PARAMETER :: PI = 3.14159265358979
  INTEGER, PARAMETER :: ind_H2SO4 = 1
  INTEGER, PARAMETER :: ind_HCl = 2
  INTEGER, PARAMETER :: ind_NH3 = 3
  INTEGER, PARAMETER :: ind_HCOOH = 4
  INTEGER, PARAMETER :: ind_RCOOH = 5
  INTEGER, PARAMETER :: ind_O1D = 6
  INTEGER, PARAMETER :: ind_SO2 = 7
  INTEGER, PARAMETER :: ind_ANOL = 8
  INTEGER, PARAMETER :: ind_H2O2 = 9
  INTEGER, PARAMETER :: ind_PAN = 10
  INTEGER, PARAMETER :: ind_TOL = 11
  INTEGER, PARAMETER :: ind_N2O5 = 12
  INTEGER, PARAMETER :: ind_XYL = 13
  INTEGER, PARAMETER :: ind_CH4 = 14
  INTEGER, PARAMETER :: ind_CRO = 15
  INTEGER, PARAMETER :: ind_HNO4 = 16
  INTEGER, PARAMETER :: ind_C2H6 = 17
  INTEGER, PARAMETER :: ind_TO2 = 18
  INTEGER, PARAMETER :: ind_XPAR = 19
  INTEGER, PARAMETER :: ind_CH3OOH = 20
  INTEGER, PARAMETER :: ind_ETHOOH = 21
  INTEGER, PARAMETER :: ind_HONO = 22
  INTEGER, PARAMETER :: ind_ETH = 23
  INTEGER, PARAMETER :: ind_CH3OH = 24
  INTEGER, PARAMETER :: ind_O3P = 25
  INTEGER, PARAMETER :: ind_CRES = 26
  INTEGER, PARAMETER :: ind_HNO3 = 27
  INTEGER, PARAMETER :: ind_CO = 28
  INTEGER, PARAMETER :: ind_PAR = 29
  INTEGER, PARAMETER :: ind_OPEN = 30
  INTEGER, PARAMETER :: ind_ISOPN = 31
  INTEGER, PARAMETER :: ind_ISOPP = 32
  INTEGER, PARAMETER :: ind_ISOPO2 = 33
  INTEGER, PARAMETER :: ind_OLET = 34
  INTEGER, PARAMETER :: ind_ISOP = 35
  INTEGER, PARAMETER :: ind_HCHO = 36
  INTEGER, PARAMETER :: ind_XO2 = 37
  INTEGER, PARAMETER :: ind_AONE = 38
  INTEGER, PARAMETER :: ind_OLEI = 39
  INTEGER, PARAMETER :: ind_ETHP = 40
  INTEGER, PARAMETER :: ind_NAP = 41
  INTEGER, PARAMETER :: ind_MGLY = 42
  INTEGER, PARAMETER :: ind_ALD2 = 43
  INTEGER, PARAMETER :: ind_CH3O2 = 44
  INTEGER, PARAMETER :: ind_ISOPRD = 45
  INTEGER, PARAMETER :: ind_ANO2 = 46
  INTEGER, PARAMETER :: ind_ROOH = 47
  INTEGER, PARAMETER :: ind_RO2 = 48
  INTEGER, PARAMETER :: ind_ONIT = 49
  INTEGER, PARAMETER :: ind_OH = 50
  INTEGER, PARAMETER :: ind_O3 = 51
  INTEGER, PARAMETER :: ind_NO2 = 52
  INTEGER, PARAMETER :: ind_C2O3 = 53
  INTEGER, PARAMETER :: ind_NO3 = 54
  INTEGER, PARAMETER :: ind_HO2 = 55
  INTEGER, PARAMETER :: ind_NO = 56
  INTEGER, PARAMETER :: ind_H2O = 57
  INTEGER, PARAMETER :: ind_M = 58
  INTEGER, PARAMETER :: indf_H2O = 1
  INTEGER, PARAMETER :: indf_M = 2
END MODULE cbmz_bb_Parameters
