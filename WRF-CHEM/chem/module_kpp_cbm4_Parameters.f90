MODULE cbm4_Parameters
  USE cbm4_Precision
  PUBLIC
  SAVE
  INTEGER, PARAMETER :: NSPEC = 33
  INTEGER, PARAMETER :: NVAR = 32
  INTEGER, PARAMETER :: NVARACT = 32
  INTEGER, PARAMETER :: NFIX = 1
  INTEGER, PARAMETER :: NREACT = 81
  INTEGER, PARAMETER :: NVARST = 1
  INTEGER, PARAMETER :: NFIXST = 33
  INTEGER, PARAMETER :: NONZERO = 276
  INTEGER, PARAMETER :: LU_NONZERO = 300
  INTEGER, PARAMETER :: CNVAR = 33
  INTEGER, PARAMETER :: NLOOKAT = 0
  INTEGER, PARAMETER :: NMONITOR = 0
  INTEGER, PARAMETER :: NMASS = 1
  REAL(kind=dp), PARAMETER :: PI = 3.14159265358979
  INTEGER, PARAMETER :: ind_O1D_CB4 = 1
  INTEGER, PARAMETER :: ind_H2O2 = 2
  INTEGER, PARAMETER :: ind_PAN = 3
  INTEGER, PARAMETER :: ind_CRO = 4
  INTEGER, PARAMETER :: ind_TOL = 5
  INTEGER, PARAMETER :: ind_N2O5 = 6
  INTEGER, PARAMETER :: ind_XYL = 7
  INTEGER, PARAMETER :: ind_XO2N = 8
  INTEGER, PARAMETER :: ind_HONO = 9
  INTEGER, PARAMETER :: ind_PNA = 10
  INTEGER, PARAMETER :: ind_TO2 = 11
  INTEGER, PARAMETER :: ind_HNO3 = 12
  INTEGER, PARAMETER :: ind_ROR = 13
  INTEGER, PARAMETER :: ind_CRES = 14
  INTEGER, PARAMETER :: ind_MGLY = 15
  INTEGER, PARAMETER :: ind_CO = 16
  INTEGER, PARAMETER :: ind_ETH = 17
  INTEGER, PARAMETER :: ind_XO2 = 18
  INTEGER, PARAMETER :: ind_OPEN = 19
  INTEGER, PARAMETER :: ind_PAR = 20
  INTEGER, PARAMETER :: ind_HCHO = 21
  INTEGER, PARAMETER :: ind_ISOP = 22
  INTEGER, PARAMETER :: ind_OLE = 23
  INTEGER, PARAMETER :: ind_ALD2 = 24
  INTEGER, PARAMETER :: ind_O3 = 25
  INTEGER, PARAMETER :: ind_NO2 = 26
  INTEGER, PARAMETER :: ind_HO = 27
  INTEGER, PARAMETER :: ind_HO2 = 28
  INTEGER, PARAMETER :: ind_O = 29
  INTEGER, PARAMETER :: ind_NO3 = 30
  INTEGER, PARAMETER :: ind_NO = 31
  INTEGER, PARAMETER :: ind_C2O3 = 32
  INTEGER, PARAMETER :: ind_H2O = 33
  INTEGER, PARAMETER :: indf_H2O = 1
END MODULE cbm4_Parameters
