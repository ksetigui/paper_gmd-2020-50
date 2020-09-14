MODULE cri_mosaic_4bin_aq_Parameters
  USE cri_mosaic_4bin_aq_Precision
  PUBLIC
  SAVE
  INTEGER, PARAMETER :: NSPEC = 235
  INTEGER, PARAMETER :: NVAR = 233
  INTEGER, PARAMETER :: NVARACT = 230
  INTEGER, PARAMETER :: NFIX = 2
  INTEGER, PARAMETER :: NREACT = 637
  INTEGER, PARAMETER :: NVARST = 1
  INTEGER, PARAMETER :: NFIXST = 234
  INTEGER, PARAMETER :: NONZERO = 1887
  INTEGER, PARAMETER :: LU_NONZERO = 2254
  INTEGER, PARAMETER :: CNVAR = 234
  INTEGER, PARAMETER :: NLOOKAT = 0
  INTEGER, PARAMETER :: NMONITOR = 0
  INTEGER, PARAMETER :: NMASS = 1
  REAL(kind=dp), PARAMETER :: PI = 3.14159265358979
  INTEGER, PARAMETER :: ind_HCl = 1
  INTEGER, PARAMETER :: ind_NH3 = 2
  INTEGER, PARAMETER :: ind_ClNO2 = 3
  INTEGER, PARAMETER :: ind_H2SO4 = 4
  INTEGER, PARAMETER :: ind_SO3 = 5
  INTEGER, PARAMETER :: ind_MSA = 6
  INTEGER, PARAMETER :: ind_DMSO2 = 7
  INTEGER, PARAMETER :: ind_O1D = 8
  INTEGER, PARAMETER :: ind_C2H6 = 9
  INTEGER, PARAMETER :: ind_HSO3 = 10
  INTEGER, PARAMETER :: ind_NC4H10 = 11
  INTEGER, PARAMETER :: ind_CH3CL = 12
  INTEGER, PARAMETER :: ind_CH2CL2 = 13
  INTEGER, PARAMETER :: ind_CHCL3 = 14
  INTEGER, PARAMETER :: ind_CH3CCL3 = 15
  INTEGER, PARAMETER :: ind_CDICLETH = 16
  INTEGER, PARAMETER :: ind_TDICLETH = 17
  INTEGER, PARAMETER :: ind_TRICLETH = 18
  INTEGER, PARAMETER :: ind_TCE = 19
  INTEGER, PARAMETER :: ind_CH4 = 20
  INTEGER, PARAMETER :: ind_TM135B = 21
  INTEGER, PARAMETER :: ind_OETHTOL = 22
  INTEGER, PARAMETER :: ind_METHTOL = 23
  INTEGER, PARAMETER :: ind_PETHTOL = 24
  INTEGER, PARAMETER :: ind_DIME35EB = 25
  INTEGER, PARAMETER :: ind_C3H8 = 26
  INTEGER, PARAMETER :: ind_OXYL = 27
  INTEGER, PARAMETER :: ind_TM123B = 28
  INTEGER, PARAMETER :: ind_TM124B = 29
  INTEGER, PARAMETER :: ind_BENZENE = 30
  INTEGER, PARAMETER :: ind_TOLUENE = 31
  INTEGER, PARAMETER :: ind_N2O5 = 32
  INTEGER, PARAMETER :: ind_MEK = 33
  INTEGER, PARAMETER :: ind_TNCARB15 = 34
  INTEGER, PARAMETER :: ind_CH3O2NO2 = 35
  INTEGER, PARAMETER :: ind_C2H2 = 36
  INTEGER, PARAMETER :: ind_CH3OH = 37
  INTEGER, PARAMETER :: ind_CARB15 = 38
  INTEGER, PARAMETER :: ind_CCARB12 = 39
  INTEGER, PARAMETER :: ind_RCOOH25 = 40
  INTEGER, PARAMETER :: ind_HONO = 41
  INTEGER, PARAMETER :: ind_DMSO = 42
  INTEGER, PARAMETER :: ind_C2H5OH = 43
  INTEGER, PARAMETER :: ind_IPROPOL = 44
  INTEGER, PARAMETER :: ind_CARB12 = 45
  INTEGER, PARAMETER :: ind_HCOOH = 46
  INTEGER, PARAMETER :: ind_AROH14 = 47
  INTEGER, PARAMETER :: ind_RAROH14 = 48
  INTEGER, PARAMETER :: ind_AROH17 = 49
  INTEGER, PARAMETER :: ind_RAROH17 = 50
  INTEGER, PARAMETER :: ind_NPROPOL = 51
  INTEGER, PARAMETER :: ind_CH3CO2H = 52
  INTEGER, PARAMETER :: ind_CH3SO3 = 53
  INTEGER, PARAMETER :: ind_RN9NO3 = 54
  INTEGER, PARAMETER :: ind_RN12NO3 = 55
  INTEGER, PARAMETER :: ind_RN10OOH = 56
  INTEGER, PARAMETER :: ind_RN19OOH = 57
  INTEGER, PARAMETER :: ind_RN8OOH = 58
  INTEGER, PARAMETER :: ind_RN11OOH = 59
  INTEGER, PARAMETER :: ind_NRN6OOH = 60
  INTEGER, PARAMETER :: ind_NRN12OOH = 61
  INTEGER, PARAMETER :: ind_C2H5CO3H = 62
  INTEGER, PARAMETER :: ind_HNO4 = 63
  INTEGER, PARAMETER :: ind_RU14NO3 = 64
  INTEGER, PARAMETER :: ind_NRU14OOH = 65
  INTEGER, PARAMETER :: ind_HOC2H4NO3 = 66
  INTEGER, PARAMETER :: ind_RTN28NO3 = 67
  INTEGER, PARAMETER :: ind_RTN28OOH = 68
  INTEGER, PARAMETER :: ind_RTN26OOH = 69
  INTEGER, PARAMETER :: ind_NRTN28OOH = 70
  INTEGER, PARAMETER :: ind_RTN25OOH = 71
  INTEGER, PARAMETER :: ind_RTN24OOH = 72
  INTEGER, PARAMETER :: ind_C2H5OOH = 73
  INTEGER, PARAMETER :: ind_IC3H7OOH = 74
  INTEGER, PARAMETER :: ind_RTX24OOH = 75
  INTEGER, PARAMETER :: ind_CH3CO3H = 76
  INTEGER, PARAMETER :: ind_ANHY = 77
  INTEGER, PARAMETER :: ind_TNCARB12 = 78
  INTEGER, PARAMETER :: ind_RN16NO3 = 79
  INTEGER, PARAMETER :: ind_RN16OOH = 80
  INTEGER, PARAMETER :: ind_RN17OOH = 81
  INTEGER, PARAMETER :: ind_RN12OOH = 82
  INTEGER, PARAMETER :: ind_NRN9OOH = 83
  INTEGER, PARAMETER :: ind_IC3H7NO3 = 84
  INTEGER, PARAMETER :: ind_HOCH2CO3H = 85
  INTEGER, PARAMETER :: ind_RU10OOH = 86
  INTEGER, PARAMETER :: ind_RU12PAN = 87
  INTEGER, PARAMETER :: ind_NRU12OOH = 88
  INTEGER, PARAMETER :: ind_RA13OOH = 89
  INTEGER, PARAMETER :: ind_CH3NO3 = 90
  INTEGER, PARAMETER :: ind_C2H5NO3 = 91
  INTEGER, PARAMETER :: ind_RA16OOH = 92
  INTEGER, PARAMETER :: ind_CH3OOH = 93
  INTEGER, PARAMETER :: ind_RN10NO3 = 94
  INTEGER, PARAMETER :: ind_RTN23OOH = 95
  INTEGER, PARAMETER :: ind_RTN10OOH = 96
  INTEGER, PARAMETER :: ind_RTN25NO3 = 97
  INTEGER, PARAMETER :: ind_RN19NO3 = 98
  INTEGER, PARAMETER :: ind_RTX28NO3 = 99
  INTEGER, PARAMETER :: ind_RTX24NO3 = 100
  INTEGER, PARAMETER :: ind_HOC2H4OOH = 101
  INTEGER, PARAMETER :: ind_RTX22NO3 = 102
  INTEGER, PARAMETER :: ind_RTX22OOH = 103
  INTEGER, PARAMETER :: ind_RN15NO3 = 104
  INTEGER, PARAMETER :: ind_RN18NO3 = 105
  INTEGER, PARAMETER :: ind_RA25OOH = 106
  INTEGER, PARAMETER :: ind_H2O2 = 107
  INTEGER, PARAMETER :: ind_TXCARB22 = 108
  INTEGER, PARAMETER :: ind_DMS = 109
  INTEGER, PARAMETER :: ind_MSIA = 110
  INTEGER, PARAMETER :: ind_CH3SCH2OO = 111
  INTEGER, PARAMETER :: ind_PPN = 112
  INTEGER, PARAMETER :: ind_PHAN = 113
  INTEGER, PARAMETER :: ind_RU14OOH = 114
  INTEGER, PARAMETER :: ind_RU12OOH = 115
  INTEGER, PARAMETER :: ind_RN14OOH = 116
  INTEGER, PARAMETER :: ind_MPAN = 117
  INTEGER, PARAMETER :: ind_RN9OOH = 118
  INTEGER, PARAMETER :: ind_RA13NO3 = 119
  INTEGER, PARAMETER :: ind_ARNOH14 = 120
  INTEGER, PARAMETER :: ind_RA16NO3 = 121
  INTEGER, PARAMETER :: ind_ARNOH17 = 122
  INTEGER, PARAMETER :: ind_RTN26PAN = 123
  INTEGER, PARAMETER :: ind_RTN14OOH = 124
  INTEGER, PARAMETER :: ind_RTX28OOH = 125
  INTEGER, PARAMETER :: ind_NRTX28OOH = 126
  INTEGER, PARAMETER :: ind_SO2 = 127
  INTEGER, PARAMETER :: ind_RA25NO3 = 128
  INTEGER, PARAMETER :: ind_PAN = 129
  INTEGER, PARAMETER :: ind_RN18OOH = 130
  INTEGER, PARAMETER :: ind_RA19OOH = 131
  INTEGER, PARAMETER :: ind_CARB17 = 132
  INTEGER, PARAMETER :: ind_TXCARB24 = 133
  INTEGER, PARAMETER :: ind_RA22OOH = 134
  INTEGER, PARAMETER :: ind_RN13NO3 = 135
  INTEGER, PARAMETER :: ind_TNCARB11 = 136
  INTEGER, PARAMETER :: ind_RN15OOH = 137
  INTEGER, PARAMETER :: ind_CARB14 = 138
  INTEGER, PARAMETER :: ind_RA19NO3 = 139
  INTEGER, PARAMETER :: ind_CARB10 = 140
  INTEGER, PARAMETER :: ind_CARB11A = 141
  INTEGER, PARAMETER :: ind_RA22NO3 = 142
  INTEGER, PARAMETER :: ind_CH3S = 143
  INTEGER, PARAMETER :: ind_NRU14O2 = 144
  INTEGER, PARAMETER :: ind_NRN12O2 = 145
  INTEGER, PARAMETER :: ind_NRTN28O2 = 146
  INTEGER, PARAMETER :: ind_CARB9 = 147
  INTEGER, PARAMETER :: ind_NRN6O2 = 148
  INTEGER, PARAMETER :: ind_RN13AO2 = 149
  INTEGER, PARAMETER :: ind_C3H6 = 150
  INTEGER, PARAMETER :: ind_CARB3 = 151
  INTEGER, PARAMETER :: ind_NUCARB12 = 152
  INTEGER, PARAMETER :: ind_UDCARB17 = 153
  INTEGER, PARAMETER :: ind_C2H4 = 154
  INTEGER, PARAMETER :: ind_UDCARB8 = 155
  INTEGER, PARAMETER :: ind_IC3H7O2 = 156
  INTEGER, PARAMETER :: ind_RA22AO2 = 157
  INTEGER, PARAMETER :: ind_RTN24O2 = 158
  INTEGER, PARAMETER :: ind_BPINENE = 159
  INTEGER, PARAMETER :: ind_TNCARB10 = 160
  INTEGER, PARAMETER :: ind_TBUT2ENE = 161
  INTEGER, PARAMETER :: ind_RTN23NO3 = 162
  INTEGER, PARAMETER :: ind_NOA = 163
  INTEGER, PARAMETER :: ind_RA13O2 = 164
  INTEGER, PARAMETER :: ind_APINENE = 165
  INTEGER, PARAMETER :: ind_TNCARB26 = 166
  INTEGER, PARAMETER :: ind_NRN9O2 = 167
  INTEGER, PARAMETER :: ind_CARB13 = 168
  INTEGER, PARAMETER :: ind_NRTX28O2 = 169
  INTEGER, PARAMETER :: ind_UDCARB14 = 170
  INTEGER, PARAMETER :: ind_CH3SO2 = 171
  INTEGER, PARAMETER :: ind_CH3SO = 172
  INTEGER, PARAMETER :: ind_O3P = 173
  INTEGER, PARAMETER :: ind_RN18O2 = 174
  INTEGER, PARAMETER :: ind_RN13OOH = 175
  INTEGER, PARAMETER :: ind_RA25O2 = 176
  INTEGER, PARAMETER :: ind_C5H8 = 177
  INTEGER, PARAMETER :: ind_RN16AO2 = 178
  INTEGER, PARAMETER :: ind_CARB7 = 179
  INTEGER, PARAMETER :: ind_RU14O2 = 180
  INTEGER, PARAMETER :: ind_CARB16 = 181
  INTEGER, PARAMETER :: ind_RA16O2 = 182
  INTEGER, PARAMETER :: ind_UDCARB11 = 183
  INTEGER, PARAMETER :: ind_RTN10O2 = 184
  INTEGER, PARAMETER :: ind_RTX28O2 = 185
  INTEGER, PARAMETER :: ind_NRU12O2 = 186
  INTEGER, PARAMETER :: ind_CO = 187
  INTEGER, PARAMETER :: ind_HOCH2CO3 = 188
  INTEGER, PARAMETER :: ind_RA19AO2 = 189
  INTEGER, PARAMETER :: ind_CARB6 = 190
  INTEGER, PARAMETER :: ind_UCARB12 = 191
  INTEGER, PARAMETER :: ind_HNO3 = 192
  INTEGER, PARAMETER :: ind_RTN28O2 = 193
  INTEGER, PARAMETER :: ind_RTN14O2 = 194
  INTEGER, PARAMETER :: ind_CH3COCH3 = 195
  INTEGER, PARAMETER :: ind_RA19CO2 = 196
  INTEGER, PARAMETER :: ind_RA22BO2 = 197
  INTEGER, PARAMETER :: ind_HOCH2CH2O2 = 198
  INTEGER, PARAMETER :: ind_HOCH2CHO = 199
  INTEGER, PARAMETER :: ind_RU10O2 = 200
  INTEGER, PARAMETER :: ind_RU12O2 = 201
  INTEGER, PARAMETER :: ind_RN18AO2 = 202
  INTEGER, PARAMETER :: ind_UCARB10 = 203
  INTEGER, PARAMETER :: ind_C2H5CO3 = 204
  INTEGER, PARAMETER :: ind_RTN23O2 = 205
  INTEGER, PARAMETER :: ind_RN19O2 = 206
  INTEGER, PARAMETER :: ind_HCHO = 207
  INTEGER, PARAMETER :: ind_RN11O2 = 208
  INTEGER, PARAMETER :: ind_RN15O2 = 209
  INTEGER, PARAMETER :: ind_CH3OO = 210
  INTEGER, PARAMETER :: ind_RTX22O2 = 211
  INTEGER, PARAMETER :: ind_RTN25O2 = 212
  INTEGER, PARAMETER :: ind_RN12O2 = 213
  INTEGER, PARAMETER :: ind_C2H5O2 = 214
  INTEGER, PARAMETER :: ind_RTX24O2 = 215
  INTEGER, PARAMETER :: ind_RN10O2 = 216
  INTEGER, PARAMETER :: ind_RN14O2 = 217
  INTEGER, PARAMETER :: ind_RN9O2 = 218
  INTEGER, PARAMETER :: ind_CH3CO3 = 219
  INTEGER, PARAMETER :: ind_RN17O2 = 220
  INTEGER, PARAMETER :: ind_RN8O2 = 221
  INTEGER, PARAMETER :: ind_RTN26O2 = 222
  INTEGER, PARAMETER :: ind_O3 = 223
  INTEGER, PARAMETER :: ind_RN16O2 = 224
  INTEGER, PARAMETER :: ind_RN15AO2 = 225
  INTEGER, PARAMETER :: ind_RN13O2 = 226
  INTEGER, PARAMETER :: ind_C2H5CHO = 227
  INTEGER, PARAMETER :: ind_HO2 = 228
  INTEGER, PARAMETER :: ind_NO = 229
  INTEGER, PARAMETER :: ind_NO3 = 230
  INTEGER, PARAMETER :: ind_OH = 231
  INTEGER, PARAMETER :: ind_NO2 = 232
  INTEGER, PARAMETER :: ind_CH3CHO = 233
  INTEGER, PARAMETER :: ind_H2O = 234
  INTEGER, PARAMETER :: ind_M = 235
  INTEGER, PARAMETER :: indf_H2O = 1
  INTEGER, PARAMETER :: indf_M = 2
END MODULE cri_mosaic_4bin_aq_Parameters
