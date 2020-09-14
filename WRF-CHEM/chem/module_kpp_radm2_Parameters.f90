MODULE radm2_Parameters
  USE radm2_Precision
  PUBLIC
  SAVE
  INTEGER, PARAMETER :: NSPEC = 61
  INTEGER, PARAMETER :: NVAR = 59
  INTEGER, PARAMETER :: NVARACT = 55
  INTEGER, PARAMETER :: NFIX = 2
  INTEGER, PARAMETER :: NREACT = 156
  INTEGER, PARAMETER :: NVARST = 1
  INTEGER, PARAMETER :: NFIXST = 60
  INTEGER, PARAMETER :: NONZERO = 572
  INTEGER, PARAMETER :: LU_NONZERO = 659
  INTEGER, PARAMETER :: CNVAR = 60
  INTEGER, PARAMETER :: NLOOKAT = 0
  INTEGER, PARAMETER :: NMONITOR = 0
  INTEGER, PARAMETER :: NMASS = 1
  REAL(kind=dp), PARAMETER :: PI = 3.14159265358979
  INTEGER, PARAMETER :: ind_SULF = 1
  INTEGER, PARAMETER :: ind_ORA1 = 2
  INTEGER, PARAMETER :: ind_ORA2 = 3
  INTEGER, PARAMETER :: ind_CO2 = 4
  INTEGER, PARAMETER :: ind_SO2 = 5
  INTEGER, PARAMETER :: ind_ETH = 6
  INTEGER, PARAMETER :: ind_O1D = 7
  INTEGER, PARAMETER :: ind_HC5 = 8
  INTEGER, PARAMETER :: ind_HC8 = 9
  INTEGER, PARAMETER :: ind_TOL = 10
  INTEGER, PARAMETER :: ind_XYL = 11
  INTEGER, PARAMETER :: ind_TPAN = 12
  INTEGER, PARAMETER :: ind_HONO = 13
  INTEGER, PARAMETER :: ind_H2O2 = 14
  INTEGER, PARAMETER :: ind_N2O5 = 15
  INTEGER, PARAMETER :: ind_HC3 = 16
  INTEGER, PARAMETER :: ind_CH4 = 17
  INTEGER, PARAMETER :: ind_PAA = 18
  INTEGER, PARAMETER :: ind_O3P = 19
  INTEGER, PARAMETER :: ind_HNO4 = 20
  INTEGER, PARAMETER :: ind_OP1 = 21
  INTEGER, PARAMETER :: ind_CSL = 22
  INTEGER, PARAMETER :: ind_PAN = 23
  INTEGER, PARAMETER :: ind_OL2 = 24
  INTEGER, PARAMETER :: ind_HNO3 = 25
  INTEGER, PARAMETER :: ind_CO = 26
  INTEGER, PARAMETER :: ind_ISO = 27
  INTEGER, PARAMETER :: ind_OLT = 28
  INTEGER, PARAMETER :: ind_OLI = 29
  INTEGER, PARAMETER :: ind_DCB = 30
  INTEGER, PARAMETER :: ind_GLY = 31
  INTEGER, PARAMETER :: ind_XNO2 = 32
  INTEGER, PARAMETER :: ind_KET = 33
  INTEGER, PARAMETER :: ind_MGLY = 34
  INTEGER, PARAMETER :: ind_TOLP = 35
  INTEGER, PARAMETER :: ind_XYLP = 36
  INTEGER, PARAMETER :: ind_OLTP = 37
  INTEGER, PARAMETER :: ind_OLN = 38
  INTEGER, PARAMETER :: ind_XO2 = 39
  INTEGER, PARAMETER :: ind_OL2P = 40
  INTEGER, PARAMETER :: ind_HC5P = 41
  INTEGER, PARAMETER :: ind_OP2 = 42
  INTEGER, PARAMETER :: ind_HCHO = 43
  INTEGER, PARAMETER :: ind_HC8P = 44
  INTEGER, PARAMETER :: ind_TCO3 = 45
  INTEGER, PARAMETER :: ind_O3 = 46
  INTEGER, PARAMETER :: ind_ONIT = 47
  INTEGER, PARAMETER :: ind_ALD = 48
  INTEGER, PARAMETER :: ind_OLIP = 49
  INTEGER, PARAMETER :: ind_KETP = 50
  INTEGER, PARAMETER :: ind_HO2 = 51
  INTEGER, PARAMETER :: ind_MO2 = 52
  INTEGER, PARAMETER :: ind_OH = 53
  INTEGER, PARAMETER :: ind_NO3 = 54
  INTEGER, PARAMETER :: ind_ACO3 = 55
  INTEGER, PARAMETER :: ind_HC3P = 56
  INTEGER, PARAMETER :: ind_ETHP = 57
  INTEGER, PARAMETER :: ind_NO = 58
  INTEGER, PARAMETER :: ind_NO2 = 59
  INTEGER, PARAMETER :: ind_H2O = 60
  INTEGER, PARAMETER :: ind_M = 61
  INTEGER, PARAMETER :: indf_H2O = 1
  INTEGER, PARAMETER :: indf_M = 2
END MODULE radm2_Parameters
