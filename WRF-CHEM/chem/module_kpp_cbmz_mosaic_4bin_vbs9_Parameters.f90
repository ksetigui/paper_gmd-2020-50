MODULE cbmz_mosaic_4bin_vbs9_Parameters
  USE cbmz_mosaic_4bin_vbs9_Precision
  PUBLIC
  SAVE
  INTEGER, PARAMETER :: NSPEC = 138
  INTEGER, PARAMETER :: NVAR = 136
  INTEGER, PARAMETER :: NVARACT = 124
  INTEGER, PARAMETER :: NFIX = 2
  INTEGER, PARAMETER :: NREACT = 210
  INTEGER, PARAMETER :: NVARST = 1
  INTEGER, PARAMETER :: NFIXST = 137
  INTEGER, PARAMETER :: NONZERO = 957
  INTEGER, PARAMETER :: LU_NONZERO = 1029
  INTEGER, PARAMETER :: CNVAR = 137
  INTEGER, PARAMETER :: NLOOKAT = 0
  INTEGER, PARAMETER :: NMONITOR = 0
  INTEGER, PARAMETER :: NMASS = 1
  REAL(kind=dp), PARAMETER :: PI = 3.14159265358979
  INTEGER, PARAMETER :: ind_H2SO4 = 1
  INTEGER, PARAMETER :: ind_HCl = 2
  INTEGER, PARAMETER :: ind_NH3 = 3
  INTEGER, PARAMETER :: ind_PCG1_B_C = 4
  INTEGER, PARAMETER :: ind_PCG1_B_O = 5
  INTEGER, PARAMETER :: ind_PCG1_F_C = 6
  INTEGER, PARAMETER :: ind_PCG1_F_O = 7
  INTEGER, PARAMETER :: ind_HCOOH = 8
  INTEGER, PARAMETER :: ind_RCOOH = 9
  INTEGER, PARAMETER :: ind_ARO1 = 10
  INTEGER, PARAMETER :: ind_ARO2 = 11
  INTEGER, PARAMETER :: ind_ALK1 = 12
  INTEGER, PARAMETER :: ind_OLE1 = 13
  INTEGER, PARAMETER :: ind_API1 = 14
  INTEGER, PARAMETER :: ind_API2 = 15
  INTEGER, PARAMETER :: ind_LIM1 = 16
  INTEGER, PARAMETER :: ind_LIM2 = 17
  INTEGER, PARAMETER :: ind_PSD1 = 18
  INTEGER, PARAMETER :: ind_PSD2 = 19
  INTEGER, PARAMETER :: ind_PCG2_B_O = 20
  INTEGER, PARAMETER :: ind_PCG3_B_O = 21
  INTEGER, PARAMETER :: ind_PCG4_B_O = 22
  INTEGER, PARAMETER :: ind_PCG5_B_O = 23
  INTEGER, PARAMETER :: ind_PCG6_B_O = 24
  INTEGER, PARAMETER :: ind_PCG7_B_O = 25
  INTEGER, PARAMETER :: ind_PCG8_B_O = 26
  INTEGER, PARAMETER :: ind_PCG9_B_O = 27
  INTEGER, PARAMETER :: ind_PCG2_F_O = 28
  INTEGER, PARAMETER :: ind_PCG3_F_O = 29
  INTEGER, PARAMETER :: ind_PCG4_F_O = 30
  INTEGER, PARAMETER :: ind_PCG5_F_O = 31
  INTEGER, PARAMETER :: ind_PCG6_F_O = 32
  INTEGER, PARAMETER :: ind_PCG7_F_O = 33
  INTEGER, PARAMETER :: ind_PCG8_F_O = 34
  INTEGER, PARAMETER :: ind_PCG9_F_O = 35
  INTEGER, PARAMETER :: ind_O1D = 36
  INTEGER, PARAMETER :: ind_SO2 = 37
  INTEGER, PARAMETER :: ind_PCG2_B_C = 38
  INTEGER, PARAMETER :: ind_OPCG1_B_C = 39
  INTEGER, PARAMETER :: ind_PCG3_B_C = 40
  INTEGER, PARAMETER :: ind_PCG4_B_C = 41
  INTEGER, PARAMETER :: ind_PCG5_B_C = 42
  INTEGER, PARAMETER :: ind_PCG6_B_C = 43
  INTEGER, PARAMETER :: ind_PCG7_B_C = 44
  INTEGER, PARAMETER :: ind_PCG8_B_C = 45
  INTEGER, PARAMETER :: ind_PCG9_B_C = 46
  INTEGER, PARAMETER :: ind_OPCG8_B_O = 47
  INTEGER, PARAMETER :: ind_ANOL = 48
  INTEGER, PARAMETER :: ind_OPCG8_B_C = 49
  INTEGER, PARAMETER :: ind_OPCG7_B_O = 50
  INTEGER, PARAMETER :: ind_OPCG7_B_C = 51
  INTEGER, PARAMETER :: ind_OPCG6_B_O = 52
  INTEGER, PARAMETER :: ind_OPCG6_B_C = 53
  INTEGER, PARAMETER :: ind_OPCG5_B_O = 54
  INTEGER, PARAMETER :: ind_OPCG5_B_C = 55
  INTEGER, PARAMETER :: ind_OPCG4_B_O = 56
  INTEGER, PARAMETER :: ind_OPCG4_B_C = 57
  INTEGER, PARAMETER :: ind_OPCG3_B_O = 58
  INTEGER, PARAMETER :: ind_OPCG3_B_C = 59
  INTEGER, PARAMETER :: ind_OPCG2_B_C = 60
  INTEGER, PARAMETER :: ind_OPCG1_B_O = 61
  INTEGER, PARAMETER :: ind_OPCG2_B_O = 62
  INTEGER, PARAMETER :: ind_PCG2_F_C = 63
  INTEGER, PARAMETER :: ind_OPCG1_F_C = 64
  INTEGER, PARAMETER :: ind_PCG3_F_C = 65
  INTEGER, PARAMETER :: ind_PCG4_F_C = 66
  INTEGER, PARAMETER :: ind_PCG5_F_C = 67
  INTEGER, PARAMETER :: ind_PCG6_F_C = 68
  INTEGER, PARAMETER :: ind_PCG7_F_C = 69
  INTEGER, PARAMETER :: ind_PCG8_F_C = 70
  INTEGER, PARAMETER :: ind_PCG9_F_C = 71
  INTEGER, PARAMETER :: ind_OPCG8_F_O = 72
  INTEGER, PARAMETER :: ind_OPCG8_F_C = 73
  INTEGER, PARAMETER :: ind_OPCG7_F_O = 74
  INTEGER, PARAMETER :: ind_OPCG7_F_C = 75
  INTEGER, PARAMETER :: ind_OPCG6_F_O = 76
  INTEGER, PARAMETER :: ind_OPCG6_F_C = 77
  INTEGER, PARAMETER :: ind_OPCG5_F_O = 78
  INTEGER, PARAMETER :: ind_OPCG5_F_C = 79
  INTEGER, PARAMETER :: ind_OPCG4_F_O = 80
  INTEGER, PARAMETER :: ind_OPCG4_F_C = 81
  INTEGER, PARAMETER :: ind_OPCG3_F_O = 82
  INTEGER, PARAMETER :: ind_OPCG3_F_C = 83
  INTEGER, PARAMETER :: ind_OPCG2_F_C = 84
  INTEGER, PARAMETER :: ind_OPCG1_F_O = 85
  INTEGER, PARAMETER :: ind_OPCG2_F_O = 86
  INTEGER, PARAMETER :: ind_PAN = 87
  INTEGER, PARAMETER :: ind_H2O2 = 88
  INTEGER, PARAMETER :: ind_TOL = 89
  INTEGER, PARAMETER :: ind_N2O5 = 90
  INTEGER, PARAMETER :: ind_XYL = 91
  INTEGER, PARAMETER :: ind_CRO = 92
  INTEGER, PARAMETER :: ind_CH4 = 93
  INTEGER, PARAMETER :: ind_API = 94
  INTEGER, PARAMETER :: ind_LIM = 95
  INTEGER, PARAMETER :: ind_HNO4 = 96
  INTEGER, PARAMETER :: ind_TO2 = 97
  INTEGER, PARAMETER :: ind_C2H6 = 98
  INTEGER, PARAMETER :: ind_XPAR = 99
  INTEGER, PARAMETER :: ind_CH3OOH = 100
  INTEGER, PARAMETER :: ind_ETHOOH = 101
  INTEGER, PARAMETER :: ind_HONO = 102
  INTEGER, PARAMETER :: ind_ETH = 103
  INTEGER, PARAMETER :: ind_CH3OH = 104
  INTEGER, PARAMETER :: ind_CRES = 105
  INTEGER, PARAMETER :: ind_O3P = 106
  INTEGER, PARAMETER :: ind_HNO3 = 107
  INTEGER, PARAMETER :: ind_CO = 108
  INTEGER, PARAMETER :: ind_PAR = 109
  INTEGER, PARAMETER :: ind_OPEN = 110
  INTEGER, PARAMETER :: ind_ISOPN = 111
  INTEGER, PARAMETER :: ind_ISOPO2 = 112
  INTEGER, PARAMETER :: ind_ISOPP = 113
  INTEGER, PARAMETER :: ind_OLET = 114
  INTEGER, PARAMETER :: ind_ISOP = 115
  INTEGER, PARAMETER :: ind_HCHO = 116
  INTEGER, PARAMETER :: ind_XO2 = 117
  INTEGER, PARAMETER :: ind_AONE = 118
  INTEGER, PARAMETER :: ind_OLEI = 119
  INTEGER, PARAMETER :: ind_NAP = 120
  INTEGER, PARAMETER :: ind_MGLY = 121
  INTEGER, PARAMETER :: ind_ETHP = 122
  INTEGER, PARAMETER :: ind_ALD2 = 123
  INTEGER, PARAMETER :: ind_CH3O2 = 124
  INTEGER, PARAMETER :: ind_ISOPRD = 125
  INTEGER, PARAMETER :: ind_ONIT = 126
  INTEGER, PARAMETER :: ind_O3 = 127
  INTEGER, PARAMETER :: ind_RO2 = 128
  INTEGER, PARAMETER :: ind_NO2 = 129
  INTEGER, PARAMETER :: ind_C2O3 = 130
  INTEGER, PARAMETER :: ind_OH = 131
  INTEGER, PARAMETER :: ind_HO2 = 132
  INTEGER, PARAMETER :: ind_NO = 133
  INTEGER, PARAMETER :: ind_ANO2 = 134
  INTEGER, PARAMETER :: ind_ROOH = 135
  INTEGER, PARAMETER :: ind_NO3 = 136
  INTEGER, PARAMETER :: ind_H2O = 137
  INTEGER, PARAMETER :: ind_M = 138
  INTEGER, PARAMETER :: indf_H2O = 1
  INTEGER, PARAMETER :: indf_M = 2
END MODULE cbmz_mosaic_4bin_vbs9_Parameters
